"""Validate faster exact joins, then recover the one failed scoring export.

No change to ROOT, physical generation, inference or method selection. This
stage does not mutate the production manifest while G4 generation is running.
"""
from datetime import datetime
import json
from pathlib import Path
import shutil
import sys
import time

import numpy as np

project=Path(__file__).resolve().parents[3];sys.path.insert(0,str(project))
from pileup.analysis.alpha_rate_campaign import digest,save
from pileup.analysis.alpha_rate_data_hashed import export_window
from pileup.analysis.root_redundant_photon_fields import exact_redundant_fields
from pileup.analysis.alpha_correlated_study import EXPERIMENT,OUTPUT,RECORD


def same_arrays(a,b):
    with np.load(a) as x,np.load(b) as y:
        assert set(x.files)==set(y.files)
        for name in x.files:
            assert x[name].dtype==y[name].dtype and x[name].shape==y[name].shape
            assert x[name].tobytes()==y[name].tobytes(),(str(a),name)
        return len(x.files)


def perform(local,target,row):
    checks={}
    try:meta=export_window(local,target,expected=row)
    except Exception as error:
        if 'decompress' not in str(error).lower():raise
        checks['original_error']=repr(error)
        with exact_redundant_fields(export_window,checks):meta=export_window(local,target,expected=row)
    return meta,checks


def main():
    snapshot=json.loads((RECORD/'manifest-after-export-failure.json').read_text())
    cache=Path('/tmp/gagg-correlated-io-validation-20261021');cache.mkdir(exist_ok=True)
    validations=[]
    for rate in (100000,10000000,40000000):
        row=min((r for r in snapshot['runs'] if r['rate_hz']==rate and r['status']=='completed'),key=lambda r:r['seed'])
        root=project/row['root'];local=cache/root.name
        if not local.exists():shutil.copyfile(root,local)
        assert digest(root)==digest(local)==row['root_sha256']
        target=cache/('validated-'+root.stem);tick=time.monotonic()
        meta,checks=perform(local,target,row)
        counts={}
        for folder in ('waveforms','truth'):
            counts[folder]=same_arrays(target/folder/(root.stem+'.npz'),OUTPUT/'inputs'/folder/(root.stem+'.npz'))
        old=json.loads((OUTPUT/'inputs/metadata'/(root.stem+'.json')).read_text())
        ignored={'root','io_note','photon_join'}
        assert all(meta[k]==old[k] for k in meta.keys()&old.keys()-ignored)
        validations.append(dict(root=row['root'],arrays_bitwise_equal=counts,metadata_equal_except_io_provenance=True,
            join=meta['photon_join'],seconds=time.monotonic()-tick,branch_recovery=checks))
        print('HASHED_EXPORT_EQUIVALENT',rate,counts,round(time.monotonic()-tick,2),flush=True)
        local.unlink()
    save(RECORD/'hashed-export-validation.json',dict(rows=validations,inference_changed=False,
        original_exporter_preserved=True,implementation='Hash uniqueness plus exact full-record equality; rare collisions fall back to original structured join'))
    failed=next(r for r in snapshot['runs'] if r['status']=='interrupted')
    root=project/failed['root'];local=Path('/tmp/gagg-correlated-io-probe-20261021')/root.name
    assert local.exists() and digest(local)==digest(root)
    target=cache/'recovered-export';tick=time.monotonic()
    meta,checks=perform(local,target,failed)
    assert checks['one_to_one_complete_coverage'] and checks['recovered_global_fields']==['all_photon_arrival_x_mm']
    voltage_arrays=same_arrays(target/'waveforms'/(root.stem+'.npz'),OUTPUT/'inputs/waveforms'/(root.stem+'.npz'))
    source_paths=[project/'pileup/analysis/alpha_rate_data_hashed.py',project/'pileup/analysis/root_redundant_photon_fields.py',Path(__file__)]
    copies=EXPERIMENT/'io-source'
    for path in source_paths:
        destination=copies/path.relative_to(project);destination.parent.mkdir(parents=True,exist_ok=True)
        assert not destination.exists();shutil.copyfile(path,destination)
    save(RECORD/'recovered-export.json',dict(root=failed['root'],root_sha256=digest(root),local_export=str(target),
        checks=checks,measured_voltage_arrays_bitwise_equal=voltage_arrays,seconds=time.monotonic()-tick,
        inference_changed=False,original_root_preserved=True,not_a_resimulation=True))
    save(EXPERIMENT/'io-recovery.json',dict(created=datetime.now().astimezone().isoformat(),
        purpose='Scoring-only exact recovery of unreadable redundant photon fields; faster equivalent lookup',
        frozen_inference_changed=False,physics_changed=False,healthy_equivalence_checks=validations,
        failed_root=failed['root'],recovery=checks,
        source_sha256={str(p.relative_to(project)):digest(p) for p in source_paths},
        snapshot_sha256={str(p.relative_to(project)):digest(p) for p in copies.rglob('*.py')}))
    print('EXACT_GLOBAL_X_RECOVERED',checks,round(time.monotonic()-tick,2),flush=True)


if __name__=='__main__':main()
