from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import time

project=Path(__file__).resolve().parents[3];record=Path(__file__).resolve().parent
folder=project/'pileup/experiments/alpha_gamma_full_prior_20261019'
work=Path('/tmp/gagg-full-prior-20261019')
assert json.loads((record/'inference-preflight.json').read_text())['all_five_methods_reproduced']
settings=json.loads((folder/'frozen-fit-settings.json').read_text())
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
         PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(work/'mpl'),
         PYTHONPATH='/tmp/gagg-full-prior-ml-deps-20261019')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,obj):path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
def call(command,log):
    print('STAGE',log.name,flush=True);started=time.perf_counter()
    with log.open('w') as stream:
        result=subprocess.run(command,cwd=project,env=env,stdout=stream,stderr=subprocess.STDOUT)
    print('FINISHED',log.name,'code',result.returncode,'seconds',round(time.perf_counter()-started,1),flush=True)
    if result.returncode:raise RuntimeError(str(log))

for stage in ('preflight','test_mixed'):
    call([sys.executable,str(folder/'run_native.py'),'--stage',stage,'--inputs',str(work/'test-inputs'),
          '--workers','4'],record/f'native-{stage}.log')
    manifest=json.loads((folder/'manifest.json').read_text())
    if stage=='preflight':
        assert len(manifest['runs'])==4 and manifest['preflight_validated']
        assert all(r['status']=='completed' and r['validation']['macro_verified'] and r['validation']['seeds_verified'] for r in manifest['runs'])
        save(record/'native-preflight.json',dict(roots=[r['root'] for r in manifest['runs']],
            seeds=[r['seed'] for r in manifest['runs']],all_valid=True,rows=manifest['runs']))
assert len(manifest['runs'])==16 and all(r['status']=='completed' for r in manifest['runs'])
assert all(datetime.fromisoformat(r['started'])>datetime.fromisoformat(settings['created']) for r in manifest['runs'])
(work/'test-cases').mkdir(exist_ok=False);logs=work/'test-logs';logs.mkdir(exist_ok=False)
for row in manifest['runs']:
    stem=Path(row['root']).stem
    wave=work/'test-inputs/waveforms'/f'{stem}.npz'
    save(work/'test-cases'/f'{stem}.json',dict(stem=stem,waveform=str(wave.relative_to(work)),waveform_sha256=sha(wave)))

guard='''from pathlib import Path
import json
import runpy
import sys
access=[]
def audit(event,args):
    if event=='open' and isinstance(args[0],(str,bytes)):
        path=str(args[0])
        if path.endswith('.root') or '/truth/' in path:
            raise RuntimeError('Blind fit attempted ROOT/truth access: '+path)
        access.append(path)
sys.addaudithook(audit)
stem=sys.argv[sys.argv.index('--stem')+1]
runpy.run_path(str(Path(__file__).with_name('run_test.py')),run_name='__main__')
path=Path(__file__).parent/'test-logs'/(stem+'.file-access.json')
path.write_text(json.dumps(dict(stem=stem,root_truth_guard_enabled=True,fit_completed=True,
    opened_files=sorted(set(access)),number_of_opens=len(access)),indent=2)+'\\n')
'''
(work/'blind_test_runner.py').write_text(guard)

def run(row):
    stem=Path(row['root']).stem
    command=[sys.executable,str(work/'blind_test_runner.py'),'fit','--work',str(work),'--stem',stem]
    started=time.perf_counter()
    with (logs/f'{stem}.log').open('w') as stream:
        result=subprocess.run(command,cwd=work,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=1800)
    value=dict(stem=stem,returncode=result.returncode,seconds=time.perf_counter()-started,
               command=command,finished=datetime.now().astimezone().isoformat())
    print(json.dumps(value),flush=True)
    if result.returncode:raise RuntimeError(str(logs/f'{stem}.log'))
    return value

with ThreadPoolExecutor(max_workers=4) as pool:
    rows=list(pool.map(run,manifest['runs']))
save(record/'independent-fit-runs.json',rows)
call([sys.executable,str(work/'run_test.py'),'report','--project',str(project),
      '--folder',str(folder),'--work',str(work)],record/'independent-report.log')
summary=json.loads((work/'test-outputs/summary.json').read_text())
print(json.dumps(dict(pooled=[{k:r[k] for k in ('method','particle','target_count','recovered_count','full_precision')} for r in summary['pooled']],
    diagnostics=summary['diagnostics'],frozen_screen=summary['frozen_screen'])),flush=True)
