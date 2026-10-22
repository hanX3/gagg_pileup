"""Freeze the new models, then make fresh isolated transport parents in G4."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
import threading
import time

import awkward as ak

from .alpha_pair_data import PROJECT,EXPERIMENT,OUTPUT,RECORD,load
from .alpha_rate_campaign import digest,save
from .alpha_pair_study import SEEDS,SINGLE_VARIANTS
from .alpha_pair_registry import ALL_PAIRS,load_model


def freeze():
    if (OUTPUT/'freeze.json').exists():
        raise FileExistsError('Already frozen')
    assert (OUTPUT/'fitting-selection.json').exists(), 'Select fitting penalty on validation first'
    assert (OUTPUT/'development-review.json').exists(), 'Complete development comparison before generating fresh tests'
    files={}
    for variant in SINGLE_VARIANTS+ALL_PAIRS:
        for seed in SEEDS:
            folder=OUTPUT/'models'/f'{variant}_{seed}'
            load_model(f'{variant}_{seed}')
            for name in ('model.pt','config.json','frozen.json','selection.json'):
                p=folder/name;files[str(p.relative_to(PROJECT))]=digest(p)
            run=load(folder/'config.json')
            for name in (('alpha_pair_ml.py','alpha_pair_data.py','alpha_pair_ordered.py') if variant.startswith('pair_ordered') else ('alpha_pair_ml.py','alpha_pair_data.py','alpha_pair_study.py')):
                key='pileup/analysis/'+name
                assert digest(PROJECT/key)==run['source_sha256'][key],key
    parent=PROJECT/'pileup/experiments/alpha_gamma_advanced_fit_20261020'
    shutil.copytree(parent/'source/brems/g4',EXPERIMENT/'source/brems/g4')
    shutil.copyfile(parent/'bin/gagg-brems',EXPERIMENT/'bin/gagg-brems')
    (EXPERIMENT/'bin/gagg-brems').chmod(0o755)
    for package in ('analysis','pileup/analysis'):
        for p in (PROJECT/package).glob('*.py'):
            q=EXPERIMENT/'inference-source'/p.relative_to(PROJECT)
            q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
    for package in ('analysis','pileup','pileup/analysis'):
        p=EXPERIMENT/'inference-source'/package/'__init__.py';p.parent.mkdir(parents=True,exist_ok=True);p.touch(exist_ok=True)
    (EXPERIMENT/'replay.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_PAIR_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
from pileup.analysis.alpha_pair_study import main
main()
''')
    plans=[]
    for energy in (1.5,2.5,5.486,8.):
        for x,y in ((0.,0.),(2.,0.),(4.,4.)):
            for rep in range(2):
                plans.append(dict(label=f'test_alpha_{energy:g}_{x:g}_{y:g}_{rep:02}',
                    particle='alpha',energy_MeV=energy,x_mm=x,y_mm=y,role='new_test_single',
                    gamma_model='brems',seed=262051000+len(plans),preflight=len(plans)==0))
    seeds={r['seed'] for r in plans}
    for p in (PROJECT/'pileup/experiments').glob('*/manifest.json'):
        if p==EXPERIMENT/'manifest.json':continue
        for r in load(p).get('planned_runs',[]):
            assert r.get('seed') not in seeds,p
    save(EXPERIMENT/'assets/test-plan.json',dict(planned_runs=plans))
    ancestor=load(parent/'manifest.json')
    save(EXPERIMENT/'manifest.json',dict(created=datetime.now().astimezone().isoformat(),
        purpose='Post-freeze independent single-alpha G4 transports, then controlled shared-cell photon pairs',
        planned_runs=plans,runs=[],preflight_validated=False,max_workers=4,window_us=20.,gate_us=[5.,15.],
        binary_sha256=digest(EXPERIMENT/'bin/gagg-brems'),
        environment=ancestor['parents']['brems']['environment'],adequate_dynamic_range=True,
        alpha_source='1.5/2.5/5.486/8 MeV primaries; realized Q target, not energy reconstruction',
        source_sha256={str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'source').rglob('*') if p.is_file()}))
    for p in (EXPERIMENT/'inference-source').rglob('*.py'):
        files[str(p.relative_to(PROJECT))]=digest(p)
    for p in (EXPERIMENT/'source').rglob('*'):
        if p.is_file():files[str(p.relative_to(PROJECT))]=digest(p)
    for role in ('train','validation'):
        for p in (OUTPUT/'photons'/role).rglob('*'):
            if p.is_file():files[str(p.relative_to(PROJECT))]=digest(p)
    for p in (EXPERIMENT/'settings.json',EXPERIMENT/'data-split.json',EXPERIMENT/'assets/test-plan.json',
              EXPERIMENT/'bin/gagg-brems',EXPERIMENT/'replay.py',OUTPUT/'protocol.json',OUTPUT/'single-selection.json',OUTPUT/'fitting-selection.json',
              OUTPUT/'train-pool.npz',OUTPUT/'validation-pool.npz',OUTPUT/'train-shared-pair.npz',OUTPUT/'validation-shared-pair.npz'):
        files[str(p.relative_to(PROJECT))]=digest(p)
    save(OUTPUT/'freeze.json',dict(frozen_at=datetime.now().astimezone().isoformat(),files=files,
        new_test_generation_not_started=True,test_tuning=False,models=14))
    print('FROZEN',len(files),'files;24 test ROOTs planned',flush=True)


def check_freeze():
    freeze=load(OUTPUT/'freeze.json')
    for path,sha in freeze['files'].items():
        assert digest(PROJECT/path)==sha,path
    return freeze


def native(workers=4):
    from .alpha_gamma_regularized_campaign import macro
    from .alpha_rate_data import export_window,open_root
    check_freeze()
    manifest=load(EXPERIMENT/'manifest.json')
    assert 1<=workers<=4
    cache=Path(tempfile.mkdtemp(prefix='gagg-pair-native-'))
    binary=cache/'gagg-brems';shutil.copyfile(EXPERIMENT/'bin/gagg-brems',binary);binary.chmod(0o755)
    env=os.environ.copy();env.update(manifest['environment']['g4_dataset_environment'])
    completed={r['label'] for r in manifest['runs']}
    plans=[r for r in manifest['planned_runs'] if r['label'] not in completed]
    launch_lock=threading.Lock();last_launch=[0.]

    def one(row):
        row=dict(row);root=PROJECT/row['root']
        path=EXPERIMENT/'macros'/(row['label']+'.mac');log=EXPERIMENT/'logs'/(row['label']+'.log')
        body=macro(EXPERIMENT,row,root).replace('/gun/position 0 0 0 mm',f"/gun/position {row['x_mm']} {row['y_mm']} 0 mm")
        path.write_text(body)
        row.update(macro=str(path.relative_to(PROJECT)),macro_sha256=digest(path),
            log=str(log.relative_to(PROJECT)),started=datetime.now().astimezone().isoformat(),
            sipm_seed=row['seed']+100000,readout_seed=row['seed']+200000)
        started=time.monotonic()
        with log.open('w') as stream:
            with launch_lock:
                delay=1.05-(time.monotonic()-last_launch[0])
                if delay>0.:time.sleep(delay)
                process=subprocess.Popen([str(binary),str(path)],cwd=PROJECT/'g4/build',env=env,stdout=stream,stderr=subprocess.STDOUT)
                last_launch[0]=time.monotonic()
            returncode=process.wait()
        if returncode:
            raise RuntimeError(f'G4 failed {row["label"]}: '+log.read_text()[-2000:])
        local=cache/root.name;shutil.copyfile(root,local)
        sha=digest(root);assert digest(local)==sha
        meta=export_window(local,OUTPUT/'native-inputs',allowed_particles=('alpha',))
        with open_root(local) as handle:
            info=ak.to_list(handle['RunInfo'].arrays(library='ak'))[0]
            assert handle['WaveformEvent'].num_entries==handle['SiPMEvent'].num_entries==1
            assert info['random_seed']==row['seed'] and info['macro_contents']==body
            assert info['sipm_avalanche_base_seed']==row['sipm_seed']
            assert info['sipm_readout_base_seed']==row['readout_seed']
        assert meta['root_sha256']==sha and meta['primaries']==1
        assert meta['macro_contents']==body
        meta['root']=str(root.relative_to(PROJECT));meta['io_note']='Verified regular local byte copy, no links'
        meta['actual_run_seeds']={k:int(info[k]) for k in ('random_seed','sipm_avalanche_base_seed','sipm_readout_base_seed')}
        save(OUTPUT/'native-inputs/metadata'/f'{root.stem}.json',meta)
        row.update(status='completed',returncode=0,elapsed_s=time.monotonic()-started,
            root_sha256=sha,root_bytes=root.stat().st_size,finished=datetime.now().astimezone().isoformat(),
            validation=dict(runinfo_macro_verified=True,seed_verified=True,charge_balance_error_C=meta['charge_balance_error_C']))
        # Keep the small original/copy identities in the record, not the temporary copy.
        local.unlink()
        return row

    def allocate(row):
        while True:
            name=datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss')+'.root'
            root=PROJECT/'g4/data'/name
            if not root.exists():break
            time.sleep(1.)
        row=dict(row,root=str(root.relative_to(PROJECT)))
        print('NATIVE_START',row['label'],row['root'],flush=True)
        time.sleep(1.05)
        return row

    if not manifest['preflight_validated']:
        row=one(allocate(plans.pop(0)))
        manifest['runs'].append(row);manifest['preflight_validated']=True
        save(EXPERIMENT/'manifest.json',manifest)
        print('NATIVE_PREFLIGHT_PASSED',row['elapsed_s'],flush=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(one,allocate(row)) for row in plans]
        for future in as_completed(futures):
            row=future.result();manifest['runs'].append(row);save(EXPERIMENT/'manifest.json',manifest)
            print('NATIVE_DONE',len(manifest['runs']),row['label'],round(row['elapsed_s'],1),flush=True)
    records=[]
    for r in manifest['runs']:
        stem=Path(r['root']).stem
        records.append(dict(root=r['root'],seed=r['seed'],particle='alpha',root_original_sha256=r['root_sha256'],
            source_manifest=str((EXPERIMENT/'manifest.json').relative_to(PROJECT)),
            waveform=str((OUTPUT/'native-inputs/waveforms'/f'{stem}.npz').relative_to(PROJECT)),
            truth=str((OUTPUT/'native-inputs/truth'/f'{stem}.npz').relative_to(PROJECT)),
            metadata=str((OUTPUT/'native-inputs/metadata'/f'{stem}.json').relative_to(PROJECT))))
    save(OUTPUT/'native-records.json',dict(test_single=records))
    data_module=__import__('pileup.analysis.alpha_pair_data',fromlist=['native_pool'])
    data_module.native_pool(records,OUTPUT/'test-pool.npz')
    print('FRESH_NATIVE_COMPLETE',len(records),flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('freeze','native'));p.add_argument('--workers',type=int,default=4)
    a=p.parse_args()
    if a.stage=='freeze':freeze()
    else:native(a.workers)


if __name__=='__main__':main()
