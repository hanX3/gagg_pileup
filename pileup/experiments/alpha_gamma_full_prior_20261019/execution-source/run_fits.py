from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from datetime import datetime
import json
import os
import subprocess
import sys
import time
import numpy as np

project=Path(__file__).resolve().parents[3]
record=Path(__file__).resolve().parent
folder=project/'pileup/experiments/alpha_gamma_full_prior_20261019'
work=Path('/tmp/gagg-full-prior-20261019')
protocol=json.loads((folder/'protocol.json').read_text())
stems=sorted(Path(r['root']).stem for r in protocol['scoring_rows'])
logs=work/'logs';logs.mkdir(exist_ok=True)
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
         PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(work/'mpl'))

def save(path,obj):path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
def run(stem,controls=False):
    command=[sys.executable,str(work/'run_analysis.py'),'fit','--work',str(work),'--stem',stem]
    if controls:command+=['--controls']
    start=time.perf_counter()
    with (logs/f'{stem}.log').open('w') as log:
        result=subprocess.run(command,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1200)
    row=dict(stem=stem,returncode=result.returncode,seconds=time.perf_counter()-start,command=command,
             controls=controls,finished=datetime.now().astimezone().isoformat())
    print(json.dumps(row),flush=True)
    if result.returncode:raise RuntimeError(str(logs/f'{stem}.log'))
    return row

pilot=run(stems[0],True)
differences={}
for name,old in {'replay_none':'hybrid_unet_0','replay_both':'hybrid_unet_10'}.items():
    a=json.loads((work/'outputs/fits'/name/f'{stems[0]}.json').read_text())['fit']
    b=json.loads((project/protocol['parent_artifacts']/'results/validation'/old/f'{stems[0]}.json').read_text())['fit']
    maximum=0.
    for key in ('times_ns','charge_pe','weighted_cost_history'):
        np.testing.assert_allclose(a[key],b[key],rtol=1e-10,atol=1e-7)
        maximum=max(maximum,float(np.max(abs(np.asarray(a[key])-np.asarray(b[key])))))
    assert a['particle_name']==b['particle_name']
    assert a['all_components']['particle_name']==b['all_components']['particle_name']
    for key in ('times_ns','charge_pe','shape_charge_pe'):
        np.testing.assert_allclose(a['all_components'][key],b['all_components'][key],rtol=1e-10,atol=1e-7)
    assert a['pair_trials']==b['pair_trials']
    differences[name]=maximum
save(record/'pilot.json',dict(run=pilot,original_endpoints_max_abs_difference=differences))
with ThreadPoolExecutor(max_workers=4) as pool:
    rows=[pilot]+list(pool.map(run,stems[1:]))
save(record/'fit-runs.json',rows)
command=[sys.executable,str(work/'run_analysis.py'),'report','--project',str(project),
         '--folder',str(folder),'--output',str(work/'outputs')]
with (logs/'report.log').open('w') as log:
    result=subprocess.run(command,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=600)
if result.returncode:raise RuntimeError(str(logs/'report.log'))
print((logs/'report.log').read_text(),flush=True)
