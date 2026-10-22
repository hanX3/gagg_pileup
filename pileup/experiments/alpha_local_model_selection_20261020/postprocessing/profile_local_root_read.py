"""Measure a regular local read cache; canonical production ROOT stays shared."""
from pathlib import Path
import json
import shutil
import sys
import time
import numpy as np

project=Path.cwd();exp=project/'pileup/experiments/alpha_local_model_selection_20261020'
out=project/'jupyter/artifacts/alpha_local_model_selection_20261020';rec=project/'jupyter/run_records/alpha_local_model_selection_20261020'
sys.path.insert(0,str(exp/'inference-source'))
from pileup.analysis.alpha_rate_data import export_window
from pileup.analysis.alpha_rate_campaign import digest,save
stem='20261020_17h41m54s';root=project/'g4/data'/f'{stem}.root'
row=next(r for r in json.loads((exp/'manifest.json').read_text())['runs'] if Path(r['root']).stem==stem)
assert row['status']=='completed'
cache=Path('/tmp/gagg-alpha-local-read-cache');cache.mkdir(parents=True,exist_ok=True)
local=cache/root.name;tick=time.monotonic()
if not local.exists():shutil.copyfile(root,local)
assert digest(local)==row['root_sha256']
copy_seconds=time.monotonic()-tick
print('LOCAL_CACHE_READY',round(copy_seconds,2),'seconds',flush=True)
tick=time.monotonic();meta=export_window(local,out/'read-cache-equivalence',expected=row)
seconds=time.monotonic()-tick
print('LOCAL_EXPORT_DONE',round(seconds,2),'seconds',flush=True)
checks=[]
for sub in ('waveforms','truth'):
    with np.load(out/'inputs'/sub/f'{stem}.npz') as old,np.load(out/'read-cache-equivalence'/sub/f'{stem}.npz') as new:
        assert set(old.files)==set(new.files)
        for key in old.files:
            a,b=old[key],new[key]
            assert a.dtype==b.dtype and np.array_equal(a,b,equal_nan=a.dtype.kind in 'fc'),key
            checks.append(sub+'/'+key)
save(rec/'read-cache-validation.json',dict(passed=True,canonical_root=str(root.relative_to(project)),
    local_cache=str(local),root_sha256=row['root_sha256'],copy_seconds=copy_seconds,
    local_export_seconds=seconds,all_arrays_exactly_equal=checks,
    note='Same frozen exporter and same ROOT bytes; no link or model change. Canonical ROOT retained under g4/data.'))
print('LOCAL_CACHE_EQUIVALENCE_PASSED',len(checks),'arrays',flush=True)
