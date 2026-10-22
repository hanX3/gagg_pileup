"""Publish exact checked scoring exports after all predeclared G4 runs finish.

Use ordinary local byte copies and local export staging; publish only complete
arrays/metadata. Failed prior export and original ROOTs remain preserved.
"""
from datetime import datetime
import json
from pathlib import Path
import shutil
import sys
import time

import awkward as ak
import numpy as np

project=Path(__file__).resolve().parents[3];sys.path.insert(0,str(project))
from pileup.analysis.alpha_correlated_study import EXPERIMENT,OUTPUT,RECORD
from pileup.analysis.alpha_correlated_evaluate import frozen
from pileup.analysis.alpha_rate_campaign import digest,save
from pileup.analysis.alpha_rate_data_hashed import export_window,open_root
from pileup.analysis.root_redundant_photon_fields import exact_redundant_fields


def equal_npz(a,b):
    with np.load(a) as x,np.load(b) as y:
        assert set(x.files)==set(y.files)
        for k in x.files:assert x[k].dtype==y[k].dtype and x[k].shape==y[k].shape and x[k].tobytes()==y[k].tobytes(),k


def main():
    frozen();assert (RECORD/'hashed-export-validation.json').exists()
    assert (EXPERIMENT/'io-recovery.json').exists()
    io=json.loads((EXPERIMENT/'io-recovery.json').read_text())
    for name,sha in io['source_sha256'].items():assert digest(project/name)==sha,name
    manifest=json.loads((EXPERIMENT/'manifest.json').read_text())
    assert len(manifest['runs'])==23 and all(r['returncode']==0 for r in manifest['runs'])
    assert all(r['status'] in ('completed','generated','interrupted') for r in manifest['runs'])
    recovery=json.loads((RECORD/'recovered-export.json').read_text())
    cache=Path('/tmp/gagg-correlated-final-exports-20261021');cache.mkdir(exist_ok=True)
    io_rows=[]
    for row in manifest['runs']:
        if row['status']=='completed' and 'export_metadata' in row:continue
        root=project/row['root'];sha=digest(root);tick=time.monotonic();checks={}
        if row['root']==recovery['root']:
            assert sha==recovery['root_sha256']
            stage=Path(recovery['local_export']);checks=recovery['checks']
            meta=json.loads((stage/'metadata'/(root.stem+'.json')).read_text())
        else:
            local=cache/root.name;assert not local.exists();shutil.copyfile(root,local)
            try:
                assert digest(local)==sha
                stage=cache/('export-'+root.stem)
                try:meta=export_window(local,stage,expected=row)
                except Exception as error:
                    if 'decompress' not in str(error).lower():raise
                    checks['original_error']=repr(error)
                    with exact_redundant_fields(export_window,checks):meta=export_window(local,stage,expected=row)
                with open_root(local) as f:
                    info=ak.to_list(f['RunInfo'].arrays())[0]
                    assert info['random_seed']==row['seed'] and info['alpha_p11b_rate_hz']==row['rate_hz']
                    assert not info['brems_enabled'] and not info['c12_capture_enabled']
                    assert all(f[k].num_entries==1 for k in ('RunInfo','WaveformEvent','SiPMEvent'))
            finally:
                local.unlink(missing_ok=True)
        assert meta['root_sha256']==sha
        for folder in ('waveforms','truth'):
            incoming=stage/folder/(root.stem+'.npz');destination=OUTPUT/'inputs'/folder/incoming.name
            if destination.exists():equal_npz(incoming,destination)
            else:shutil.copyfile(incoming,destination)
        meta['root']=row['root'];meta['io_note']='Verified ordinary local byte copy and local export staging; checked hash join plus exact full records; original ROOT retained'
        if checks:meta['redundant_field_recovery']=checks
        destination=OUTPUT/'inputs/metadata'/(root.stem+'.json');assert not destination.exists();save(destination,meta)
        previous=row['status'];row.update(status='completed',root_sha256=sha,root_bytes=root.stat().st_size,
            export_metadata=str(destination.relative_to(project)),exported_at=datetime.now().astimezone().isoformat(),
            scoring_io_source=str((EXPERIMENT/'io-recovery.json').relative_to(project)))
        if checks:
            row['scoring_field_recovery']=checks
            save(RECORD/(root.stem+'-redundant-fields.json'),checks)
        if previous=='interrupted':
            row['legacy_status']=previous;row['legacy_interruption_scope']='Scoring export only; G4 finished successfully and was not rerun'
        io_rows.append(dict(root=row['root'],previous_status=previous,seconds=time.monotonic()-tick,recovery=checks))
        save(EXPERIMENT/'manifest.json',manifest)
        print('EXACT_EXPORT_COMPLETED',row['label'],round(time.monotonic()-tick,1),checks.get('recovered_global_fields',[]),checks.get('recovered_primary_fields',[]),flush=True)
    manifest['exports_completed_at']=datetime.now().astimezone().isoformat()
    manifest['io_transition']=dict(inference_changed=False,physics_changed=False,original_root_files_preserved=True,
        hashed_join_validated_against_legacy=True,stage_exports_locally=True,
        source=str(Path(__file__).relative_to(project)),source_sha256=digest(Path(__file__)),rows=io_rows)
    save(EXPERIMENT/'manifest.json',manifest);save(RECORD/'final-export-transition.json',manifest['io_transition'])
    assert sum(r['status']=='completed' and 'export_metadata' in r for r in manifest['runs'])==23
    print('ALL_23_NATIVE_EXPORTS_COMPLETE',flush=True)


if __name__=='__main__':main()
