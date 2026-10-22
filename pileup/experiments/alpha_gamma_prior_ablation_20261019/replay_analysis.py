"""Replay the frozen ablation from archived blind input bundles."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--workers',type=int,default=4)
p.add_argument('--stems',nargs='+',help='Optional diagnostic subset; aggregate scoring requires all 20 windows')
args=p.parse_args()
if not 1 <= args.workers <= 4:p.error('workers must be between 1 and 4')
folder=Path(__file__).resolve().parent
project=folder.parents[2]
output=args.output.resolve()
if output.exists():raise FileExistsError(output)
for filename in ('source-manifest.json','dependency-source-manifest.json'):
    manifest=json.loads((folder/filename).read_text())
    for name,expected in manifest['files'].items():
        assert hashlib.sha256((folder/name).read_bytes()).hexdigest()==expected,name
available=sorted(path.stem for path in (folder/'cases').glob('*.json'))
stems=sorted(args.stems or available)
assert stems and len(stems)==len(set(stems)) and set(stems)<=set(available)
work=Path(tempfile.mkdtemp(prefix='gagg-prior-replay-'))
print('LOCAL_WORK',work,flush=True)
for name in ('source','dependency-source'):
    shutil.copytree(folder/name,work/name)
shutil.copy2(folder/'calibration.json',work/'calibration.json')
(work/'cases').mkdir();(work/'logs').mkdir()
for stem in stems:
    shutil.copy2(folder/'cases'/f'{stem}.json',work/'cases'/f'{stem}.json')
    source=project/'jupyter/artifacts/alpha_gamma_prior_ablation_20261019/blind-inputs'/stem
    shutil.copytree(source,work/'inputs'/stem)
env=os.environ.copy()
env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
           PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(work/'mpl'))
runner=str(folder/'run_analysis.py')


def run(stem):
    cmd=[sys.executable,runner,'fit','--work',str(work),'--stem',stem]
    with (work/'logs'/f'{stem}.log').open('w') as log:
        result=subprocess.run(cmd,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180)
    if result.returncode:raise RuntimeError(f'{stem}: see {work}/logs/{stem}.log')
    return stem


with ThreadPoolExecutor(max_workers=args.workers) as pool:
    completed=list(pool.map(run,stems))
if stems==available:
    subprocess.run([sys.executable,runner,'report','--project',str(project),'--folder',str(folder),
                    '--output',str(work/'outputs')],cwd=work,env=env,check=True,timeout=300)
shutil.copytree(work/'outputs',output)
record=dict(windows=len(completed),full_twenty_window_replay=stems==available,stems=stems,
            work_directory=str(work),source_manifests={name:hashlib.sha256((folder/name).read_bytes()).hexdigest()
                for name in ('source-manifest.json','dependency-source-manifest.json')})
(output/'replay-manifest.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record),flush=True)
