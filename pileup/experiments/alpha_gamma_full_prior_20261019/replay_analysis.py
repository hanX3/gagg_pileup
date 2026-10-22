"""Replay frozen full-pipeline fits from the previous ablation's blind packets."""
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


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--stems',nargs='+',help='Optional diagnostic subset; only all 20 windows generate a report')
    args=parser.parse_args()
    if not 1 <= args.workers <= 4:parser.error('workers must be between 1 and 4')
    folder=Path(__file__).resolve().parent;project=folder.parents[2]
    output=args.output.resolve()
    if output.exists():raise FileExistsError(output)
    for filename in ('source-manifest.json','dependency-source-manifest.json'):
        for name,expected in json.loads((folder/filename).read_text())['files'].items():
            assert hashlib.sha256((folder/name).read_bytes()).hexdigest()==expected,name
    available=sorted(path.stem for path in (folder/'cases').glob('*.json'))
    stems=sorted(args.stems or available)
    assert stems and len(stems)==len(set(stems)) and set(stems)<=set(available)
    work=Path(tempfile.mkdtemp(prefix='gagg-full-prior-replay-'))
    print('LOCAL_WORK',work,flush=True)
    for name in ('source','dependency-source'):
        shutil.copytree(folder/name,work/name)
    for name in ('calibration.json','run_analysis.py'):
        shutil.copy2(folder/name,work/name)
    (work/'cases').mkdir();(work/'logs').mkdir()
    protocol=json.loads((folder/'protocol.json').read_text())
    packet_root=project/'jupyter/artifacts'/Path(protocol['ablation_parent']).name/'blind-inputs'
    for stem in stems:
        shutil.copy2(folder/'cases'/f'{stem}.json',work/'cases'/f'{stem}.json')
        shutil.copytree(packet_root/stem,work/'inputs'/stem)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
             PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(work/'mpl'))

    def run(stem):
        command=[sys.executable,str(work/'run_analysis.py'),'fit','--work',str(work),'--stem',stem]
        with (work/'logs'/f'{stem}.log').open('w') as log:
            result=subprocess.run(command,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1200)
        if result.returncode:raise RuntimeError(f'{stem}: see {work}/logs/{stem}.log')
        return stem

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        completed=list(pool.map(run,stems))
    if stems==available:
        subprocess.run([sys.executable,str(work/'run_analysis.py'),'report','--project',str(project),
                        '--folder',str(folder),'--output',str(work/'outputs')],cwd=work,env=env,check=True,timeout=600)
    shutil.copytree(work/'outputs',output)
    record=dict(windows=len(completed),full_twenty_window_replay=stems==available,stems=stems,
                work_directory=str(work),protocol_sha256=hashlib.sha256((folder/'protocol.json').read_bytes()).hexdigest())
    (output/'replay-manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)


if __name__=='__main__':main()
