"""I/O-only adapter for frozen joint digitization, with per-case checkpoints."""
from pathlib import Path
import faulthandler
import importlib.util
import json
import os
import shutil
import sys
import tempfile
import time
import zipfile

PROJECT=next(p for p in Path(__file__).resolve().parents if (p/'g4/data').is_dir()
             and (p/'pileup/experiments/alpha_pair_resolution_20261022').is_dir())
EXPERIMENT=PROJECT/'pileup/experiments/alpha_pair_resolution_20261022'
for key in ('GAGG_RESOLUTION_PROJECT','GAGG_PAIR_PROJECT','GAGG_COMPONENT_PROJECT'):
    os.environ[key]=str(PROJECT)
sys.path.insert(0,str(EXPERIMENT/'replay-source'))
import numpy as np
from pileup.analysis import alpha_pair_resolution_data as data
ORIGINAL_CASE=data.joint_case


def initialize():
    faulthandler.enable();faulthandler.dump_traceback_later(120,repeat=True)
    cache=Path(os.environ['GAGG_RESOLUTION_LOCAL_CACHE'])
    print('WORKER_INIT',os.getpid(),flush=True)
    donors=[]
    for i in range(240):
        with np.load(cache/'photons'/f'{i:04}.npz') as f:
            donors.append({k:f[k].copy() for k in f.files})
    data._DONORS=donors;data._CONFIG=json.loads((cache/'config.json').read_text())
    validation=data.load(data.RECORD/'digitizer-optimization-validation.json')
    source=EXPERIMENT/'optimization-source/pileup/analysis/alpha_pair_digitizer_optimized.py'
    assert validation['passed'] and data.digest(source)==validation['optimized_source_sha256']
    spec=importlib.util.spec_from_file_location('pileup.analysis._optimized_digitizer',source)
    optimized=importlib.util.module_from_spec(spec);spec.loader.exec_module(optimized)
    from pileup.analysis import sipm_response
    sipm_response.simulate_avalanches=optimized.simulate_avalanches
    print('WORKER_READY',os.getpid(),flush=True)


def case_name(case):
    gap='none' if case['gap_ns'] is None else f"{case['gap_ns']:g}"
    return f"{case['kind']}_{case['condition']}_{gap}_{case['replicate']:03}"


def execute(case):
    cache=Path(os.environ['GAGG_RESOLUTION_LOCAL_CACHE'])
    name=case_name(case);path=cache/'cases'/(name+'.npz')
    if path.exists():
        try:
            with np.load(path) as f:
                return f['voltage'],f['charge'],f['time'],f['components'],json.loads(str(f['metadata']))
        except (OSError,ValueError,EOFError,zipfile.BadZipFile):
            path.rename(path.with_suffix('.interrupted'))
    if case['replicate']==0:print('CASE_START',name,flush=True)
    started=time.monotonic();row=ORIGINAL_CASE(case)
    with tempfile.NamedTemporaryFile(dir=cache/'cases',prefix=name+'.',suffix='.partial',delete=False) as temp:
        np.savez_compressed(temp,voltage=row[0],charge=row[1],time=row[2],components=row[3],
                            metadata=json.dumps(row[4],allow_nan=False))
        temporary=Path(temp.name)
    temporary.replace(path)
    if case['replicate']==0:print('CASE_DONE',name,round(time.monotonic()-started,3),flush=True)
    return row


def cached_records():
    return data.load(data.OUTPUT/'native-records.json')['test_single']


def main():
    faulthandler.enable();faulthandler.dump_traceback_later(600,repeat=True)
    record=data.RECORD/'physical-local-cache.json'
    if record.exists():
        cache=Path(data.load(record)['directory'])
    else:
        cache=Path(tempfile.mkdtemp(prefix='gagg-resolution-physical-'))
        (cache/'photons').mkdir();(cache/'cases').mkdir()
        data.save(record,dict(directory=str(cache),purpose='Regular local copies and case checkpoints; physics unchanged.'))
    records=cached_records()
    expected=data.load(data.OUTPUT/'photons/manifest.json')['files']
    print('CACHE_COPY_BEGIN',cache,flush=True)
    for i in range(240):
        source=data.OUTPUT/'photons'/f'{i:04}.npz';target=cache/'photons'/source.name
        if not target.exists():shutil.copyfile(source,target)
        assert data.digest(target)==expected[source.name]
        if (i+1)%40==0:print('CACHE_VERIFIED',i+1,flush=True)
    (cache/'config.json').write_text(json.dumps(data.load(PROJECT/records[0]['metadata'])['config']))
    os.environ['GAGG_RESOLUTION_LOCAL_CACHE']=str(cache)
    print('FROZEN_CHECK_BEGIN',flush=True);data.check_freeze();print('FROZEN_CHECK_DONE',flush=True)
    data.donor_cache=cached_records
    data._init_digitizer=initialize
    data.joint_case=execute
    data.physical_controls(workers=4)
    faulthandler.cancel_dump_traceback_later()
    data.save(data.RECORD/'physical-io-adapter.json',dict(passed=True,workers=4,
        adapter_sha256=data.digest(Path(__file__)),frozen_joint_case_unchanged=True,
        equivalent_digitizer_optimization=data.load(data.RECORD/'digitizer-optimization-validation.json'),
        local_cache=str(cache),checkpoints=len(list((cache/'cases').glob('*.npz'))),
        interrupted_original_stage='Exited130 via Ctrl-C after no progress beyond photon-cache preparation; no physical waveform dataset was saved.'))
    print('CACHED_PHYSICAL_COMPLETE',flush=True)


if __name__=='__main__':main()
