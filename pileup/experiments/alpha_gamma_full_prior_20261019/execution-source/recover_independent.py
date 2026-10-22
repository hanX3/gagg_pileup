"""Resume only missing blind outputs after the outer process exited 143."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

project=Path(__file__).resolve().parents[3];record=Path(__file__).resolve().parent
folder=project/'pileup/experiments/alpha_gamma_full_prior_20261019'
work=Path('/tmp/gagg-full-prior-20261019')
methods=['baseline','unet','hybrid_none','hybrid_both','class_only']
runs=json.loads((folder/'manifest.json').read_text())['runs']
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,obj):path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
completed=[];pending=[];hashes={}
for row in runs:
    stem=Path(row['root']).stem
    files=[work/'test-outputs/fits'/m/f'{stem}.json' for m in methods]
    audit=work/'test-logs'/f'{stem}.file-access.json'
    if all(p.exists() for p in files) and audit.exists():
        assert json.loads(audit.read_text())['fit_completed']
        for path in files:
            json.loads(path.read_text());hashes[str(path.relative_to(work))]=sha(path)
        completed.append(dict(stem=stem,fit_completed=True,returncode=None,
            recovery='Existing five complete JSON results and successful file-access record; outer exit status not persisted'))
    else:
        assert not any(p.exists() for p in files) and not audit.exists()
        pending.append(stem)
assert len(completed)==14 and len(pending)==2
old=record/'interrupted-test-logs';old.mkdir(exist_ok=False)
for stem in pending:
    shutil.copy2(work/'test-logs'/f'{stem}.log',old/f'{stem}.log')
save(record/'independent-interruption.json',dict(observed_parent_exit_code=143,
    observed=datetime.now().astimezone().isoformat(),cause='Outer process termination; exact external cause not established',
    complete_windows=14,pending_stems=pending,existing_result_sha256=hashes,
    settings_sha256=sha(folder/'frozen-fit-settings.json'),settings_changed=False,roots_regenerated=False,
    recovery='Rerun only the two windows with no saved method results; keep same input, model, settings and blind-access guard'))
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
         PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(work/'mpl'),
         PYTHONPATH='/tmp/gagg-full-prior-ml-deps-20261019')
def run(stem):
    command=[sys.executable,str(work/'blind_test_runner.py'),'fit','--work',str(work),'--stem',stem]
    started=time.perf_counter();print('RESUME',stem,flush=True)
    with (work/'test-logs'/f'{stem}.log').open('w') as stream:
        result=subprocess.run(command,cwd=work,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=1800)
    row=dict(stem=stem,returncode=result.returncode,fit_completed=result.returncode==0,
             seconds=time.perf_counter()-started,command=command,finished=datetime.now().astimezone().isoformat())
    print(json.dumps(row),flush=True)
    if result.returncode:raise RuntimeError(stem)
    return row
with ThreadPoolExecutor(max_workers=2) as pool:
    resumed=list(pool.map(run,pending))
for relative,expected in hashes.items():assert sha(work/relative)==expected
save(record/'independent-fit-runs.json',completed+resumed)
command=[sys.executable,str(work/'run_test.py'),'report','--project',str(project),'--folder',str(folder),'--work',str(work)]
with (record/'independent-report.log').open('w') as stream:
    result=subprocess.run(command,cwd=project,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=600)
if result.returncode:raise RuntimeError(str(record/'independent-report.log'))
summary=json.loads((work/'test-outputs/summary.json').read_text())
print(json.dumps(dict(pooled=[{k:r[k] for k in ('method','particle','target_count','recovered_count','full_precision')} for r in summary['pooled']],
    diagnostics=summary['diagnostics'],frozen_screen=summary['frozen_screen'])),flush=True)
