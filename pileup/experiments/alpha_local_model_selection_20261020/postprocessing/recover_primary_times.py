"""Validate exact redundant-field recovery, then recover one scoring export."""
import importlib.util
import json
from pathlib import Path
import shutil
import sys

import numpy as np

P=Path.cwd();R=P/'jupyter/run_records/alpha_local_model_selection_20261020'
E=P/'pileup/experiments/alpha_local_model_selection_20261020'
O=P/'jupyter/artifacts/alpha_local_model_selection_20261020'
sys.path.insert(0,str(P))
from pileup.analysis.root_primary_time_recovery import scoring_time_recovery
adapter=P/'pileup/analysis/root_primary_time_recovery.py'
spec=importlib.util.spec_from_file_location('cached_native',R/'resume_with_read_cache.py')
io=importlib.util.module_from_spec(spec);spec.loader.exec_module(io)
m=json.loads((E/'manifest.json').read_text())
validation={}
for label in ('alpha_100000000_00','alpha_100000000_01'):
    row=next(r for r in m['runs'] if r['label']==label)
    stem=Path(row['root']).stem;local=Path('/tmp/gagg-alpha-local-read-cache')/(stem+'.root')
    sidecar=json.loads(local.with_name(stem+'-read-cache.json').read_text())
    assert io.digest(local)==sidecar['sha256']
    folder=O/('redundant-time-reference' if label.endswith('_00') else 'redundant-time-recovery')
    assert not folder.exists();checks={}
    with scoring_time_recovery(io.export_window,checks):
        meta=io.export_window(local,folder,expected={**row,'root_sha256':sidecar['sha256']})
    assert checks['one_to_one_complete_coverage']
    if label.endswith('_00'):
        keys=[]
        for kind in ('waveforms','truth'):
            with np.load(O/'inputs'/kind/(stem+'.npz')) as a,np.load(folder/kind/(stem+'.npz')) as b:
                assert set(a.files)==set(b.files)
                for key in a.files:
                    assert a[key].dtype==b[key].dtype
                    assert np.array_equal(a[key],b[key],equal_nan=a[key].dtype.kind in 'fc'),(kind,key)
                    keys.append(kind+'/'+key)
        validation.update(reference_root=row['root'],reference_checks=checks,
            all_arrays_exactly_equal=keys,adapter_source_sha256=io.digest(adapter),reference_passed=True)
        io.save(R/'redundant-time-validation.json',validation)
        print('REDUNDANT_TIME_REFERENCE_PASSED',len(keys),checks,flush=True)
    else:
        assert validation['reference_passed']
        for kind in ('waveforms','truth'):
            dest=O/'inputs'/kind/(stem+'.npz')
            if dest.exists():
                old=R/'before-redundant-time-recovery'/kind/dest.name
                old.parent.mkdir(parents=True,exist_ok=True);assert not old.exists();shutil.copyfile(dest,old)
            shutil.copyfile(folder/kind/(stem+'.npz'),dest)
        meta.update(root=str(P/row['root']),read_cache=str(local),scoring_field_recovery=checks,
            recovery_method='Exact unique stored x/y/energy float-bit join to duplicated WaveformEvent time')
        io.save(O/'inputs/metadata'/(stem+'.json'),meta)
        row.update(status='completed',root_sha256=meta['root_sha256'],root_bytes=(P/row['root']).stat().st_size,
            export_metadata=str((O/'inputs/metadata'/(stem+'.json')).relative_to(P)),scoring_field_recovery=checks)
        m['scoring_export_recovery']=dict(root=row['root'],method=meta['recovery_method'],
            adapter=str(adapter.relative_to(P)),adapter_sha256=io.digest(adapter),
            validation_record=str((R/'redundant-time-validation.json').relative_to(P)),fit_inputs_unchanged=True)
        io.save(E/'manifest.json',m)
        path=P/m['capacity_extension']['experiment']/'manifest.json';cm=json.loads(path.read_text())
        cm['runs']=[r for r in m['runs'] if r['rate_hz']==100_000_000];io.save(path,cm)
        validation.update(passed=True,recovered_root=row['root'],recovered_checks=checks,
            scoring_only=True,blind_voltage_unchanged=True)
        io.save(R/'redundant-time-validation.json',validation)
        print('REDUNDANT_TIME_EXPORT_RECOVERED',checks,flush=True)
