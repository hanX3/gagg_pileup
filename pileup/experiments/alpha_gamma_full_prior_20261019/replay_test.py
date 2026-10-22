"""Replay five frozen methods on the archived independent native test inputs.

Run from the project root with analysis requirements and CPU PyTorch installed.
The original U-Net weights and local input archive are required. No generation,
training, tuning or overwrite is performed.
"""
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
    parser.add_argument('--models',type=Path,default=Path('jupyter/artifacts/alpha_gamma_full_prior_20261019/models/unet'))
    parser.add_argument('--inputs',type=Path,default=Path('jupyter/artifacts/alpha_gamma_full_prior_20261019/test-inputs'))
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--stems',nargs='+',help='Optional subset; report requires the complete 16-window cohort')
    args=parser.parse_args()
    if not 1<=args.workers<=4:parser.error('workers must be between 1 and 4')
    folder=Path(__file__).resolve().parent;project=folder.parents[2]
    if Path.cwd().resolve()!=project:parser.error('Run from the project root')
    output=args.output.resolve()
    if output.exists():raise FileExistsError(output)
    def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    for name in ('source-manifest.json','dependency-source-manifest.json','test-source-manifest.json'):
        for path,expected in json.loads((folder/name).read_text())['files'].items():
            assert sha(folder/path)==expected,path
    frozen=json.loads((folder/'frozen-fit-settings.json').read_text())
    assert sha(args.models/'model.pt')==frozen['neural_model_sha256']
    runs=json.loads((folder/'manifest.json').read_text())['runs']
    assert len(runs)==16 and all(r['status']=='completed' for r in runs)
    available=sorted(Path(r['root']).stem for r in runs)
    stems=sorted(args.stems or available)
    assert stems and len(set(stems))==len(stems) and set(stems)<=set(available)
    work=Path(tempfile.mkdtemp(prefix='gagg-full-prior-test-replay-'))
    print('LOCAL_WORK',work,flush=True)
    for name in ('source','dependency-source','test-source'):
        shutil.copytree(folder/name,work/name)
    for name in ('calibration.json','run_test.py'):
        shutil.copy2(folder/name,work/name)
    shutil.copy2(folder/'frozen-fit-settings.json',work/'test-settings.json')
    destination=work/'models/unet';destination.mkdir(parents=True)
    for name in ('model.pt','frozen-model.json','training-config.json'):
        shutil.copy2(args.models/name,destination/name)
    (work/'test-cases').mkdir();(work/'logs').mkdir()
    for stem in stems:
        for sub,extension in (('waveforms','npz'),('truth','npz'),('metadata','json')):
            destination=work/'test-inputs'/sub/f'{stem}.{extension}'
            destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(args.inputs/sub/f'{stem}.{extension}',destination)
        wave=work/'test-inputs/waveforms'/f'{stem}.npz'
        case=dict(stem=stem,waveform=str(wave.relative_to(work)),waveform_sha256=sha(wave))
        (work/'test-cases'/f'{stem}.json').write_text(json.dumps(case,indent=2)+'\n')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
             PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(work/'mpl'))
    def run(stem):
        command=[sys.executable,str(work/'run_test.py'),'fit','--work',str(work),'--stem',stem]
        with (work/'logs'/f'{stem}.log').open('w') as log:
            result=subprocess.run(command,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
        if result.returncode:raise RuntimeError(str(work/'logs'/f'{stem}.log'))
        return stem
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        completed=list(pool.map(run,stems))
    if stems==available:
        subprocess.run([sys.executable,str(work/'run_test.py'),'report','--project',str(project),
                        '--folder',str(folder),'--work',str(work)],cwd=work,env=env,check=True,timeout=600)
    shutil.copytree(work/'test-outputs',output)
    record=dict(windows=len(completed),stems=stems,full_sixteen_window_replay=stems==available,
                work_directory=str(work),settings_sha256=sha(folder/'frozen-fit-settings.json'))
    (output/'replay-manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)


if __name__=='__main__':main()
