"""Replay frozen methods into a new directory, retaining ROOT and model identity.

Run from project root with the analysis requirements and CPU PyTorch installed.
No training or tuning is performed; no existing result is overwritten.
"""
from pathlib import Path
import argparse,json,shutil,subprocess,sys,tempfile,os

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--models',type=Path,default=Path('jupyter/artifacts/alpha_gamma_ml_20261016'))
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    folder=Path(__file__).resolve().parent;project=folder.parents[2]
    assert Path.cwd().resolve()==project,'Run from project root'
    assert 1<=args.workers<=4
    output=args.output.resolve();relative_output=output.relative_to(project)
    output.mkdir(parents=True,exist_ok=False)
    cache=Path(tempfile.mkdtemp(prefix='gagg-mixed-replay-'))
    local_folder=cache/folder.relative_to(project);local_folder.mkdir(parents=True)
    local_output=cache/relative_output;local_output.mkdir(parents=True)
    (cache/'g4/data').mkdir(parents=True)
    shutil.copy2(project/'README.md',cache/'README.md')
    for part in ('ml-source-v2','postprocessing-source-v2'):
        shutil.copytree(folder/part,local_folder/part)
    for name in ('calibration.json','manifest.json','frozen-fit-settings.json','run_ml.py','run_postprocessing.py'):
        shutil.copy2(folder/name,local_folder/name)
    print('Machine-local replay cache:',cache,flush=True)
    sys.path.insert(0,str(folder/'ml-source-v2'))
    from pileup.analysis.alpha_rate_data import export_window
    from pileup.analysis.alpha_rate_campaign import digest
    frozen=json.loads((folder/'frozen-fit-settings.json').read_text())
    for kind in ('cnn','unet'):
        source=args.models/kind
        assert digest(source/'model.pt')==frozen['neural_models'][kind]['model_sha256']
        shutil.copytree(source,local_output/kind)
    manifest=json.loads((folder/'manifest.json').read_text())
    for row in manifest['runs']:
        if row['role']!='test_mixed':continue
        assert digest(project/row['root'])==row['root_sha256']
        export_window(project/row['root'],local_output/'inputs',allowed_particles=('alpha','gamma'))
    processes=[]
    for i in range(args.workers):
        command=[sys.executable,str(local_folder/'run_ml.py'),'pileup.analysis.alpha_gamma_ml_evaluate','evaluate',
                 '--folder',str(local_folder),'--output',str(local_output),'--role','test','--methods',*frozen['common_methods'],
                 '--shard',str(i),str(args.workers)]
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
        processes.append(subprocess.Popen(command,cwd=cache,env=env))
    codes=[p.wait() for p in processes]
    assert codes==[0]*args.workers,codes
    subprocess.run([sys.executable,str(local_folder/'run_postprocessing.py'),'pileup.analysis.alpha_gamma_ml_report',
        '--folder',str(local_folder),'--output',str(local_output),'--role','test'],cwd=cache,check=True,env=env)
    shutil.copytree(local_output,output,dirs_exist_ok=True)
    (output/'replay-cache.json').write_text(json.dumps(dict(cache=str(cache),
        policy='Machine-local compute, frozen models, native ROOTs unchanged; archived only after success'),indent=2))

if __name__=='__main__':main()
