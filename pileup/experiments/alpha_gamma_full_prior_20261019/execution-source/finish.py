"""Finalize only after all sixteen blind results and the separate report exist."""
from pathlib import Path
from datetime import datetime
import json
import os
import subprocess
import sys
import time

project=Path(__file__).resolve().parents[3];record=Path(__file__).resolve().parent
folder=project/'pileup/experiments/alpha_gamma_full_prior_20261019'
work=Path('/tmp/gagg-full-prior-20261019')
deadline=time.monotonic()+1200
print('WAITING_FOR_COMPLETE_INDEPENDENT_REPORT',flush=True)
while True:
    try:
        assert len(json.loads((record/'independent-fit-runs.json').read_text()))==16
        assert json.loads((work/'test-outputs/summary.json').read_text())['windows']==16
        events=json.loads((work/'test-outputs/events.json').read_text());assert events
        break
    except (OSError,ValueError,AssertionError):
        if time.monotonic()>deadline:raise TimeoutError('Independent run has not finished; existing outputs retained')
        time.sleep(2)
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
         PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(work/'mpl'))
def call(command,name):
    print('START',name,datetime.now().astimezone().isoformat(),flush=True)
    with (record/(name+'.log')).open('w') as stream:
        result=subprocess.run(command,cwd=project,env=env,stdout=stream,stderr=subprocess.STDOUT)
    print('DONE',name,result.returncode,flush=True)
    if result.returncode:raise RuntimeError(str(record/(name+'.log')))
call([sys.executable,str(folder/'postprocessing-source-v3/make_figures.py'),'--output',str(work/'test-outputs'),
      '--recovery-info',str(record/'independent-interruption.json')],'independent-figures')
call([sys.executable,str(record/'archive_test.py')],'archive-test')
call([sys.executable,str(record/'write_docs.py')],'write-docs')
pending=record/'validation.json'
assert not pending.exists()
pending.write_text(json.dumps(dict(status='pending',note='Full validation runs after Notebook outputs and links are generated'),indent=2)+'\n')
call([sys.executable,str(record/'write_notebook.py')],'write-notebook')
try:
    call([sys.executable,str(record/'validate.py')],'validation')
except BaseException:
    pending.write_text(json.dumps(dict(status='failed',note='See validation.log; results retained for repair'),indent=2)+'\n')
    raise
print('COMPLETE',flush=True)
