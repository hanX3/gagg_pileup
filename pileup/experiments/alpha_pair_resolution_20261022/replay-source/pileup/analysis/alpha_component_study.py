"""Orchestrate the frozen component-supervision ablation and fresh native test."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import json
from pathlib import Path
import platform
import shutil
import subprocess

import numpy as np
import torch

from .alpha_component_ml import (PROJECT,EXPERIMENT,OUTPUT,VARIANTS,SEEDS,
                                 native_component_cache, train, load_model)
from .alpha_rate_campaign import digest, save

RECORD = PROJECT/'jupyter/run_records/alpha_component_ml_20261021'


def load(path):
    return json.loads(Path(path).read_text())


def prepare():
    EXPERIMENT.mkdir(parents=True,exist_ok=False)
    OUTPUT.mkdir(parents=True,exist_ok=False)
    RECORD.mkdir(parents=True,exist_ok=False)
    for name in ('assets','bin','macros','logs'):
        (EXPERIMENT/name).mkdir()
    shutil.copyfile(PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb',RECORD/'notebook.before.ipynb')
    shutil.copyfile('/tmp/gagg-ml-component-initial-status.txt',RECORD/'git-status.before.txt')
    (RECORD/'working-changes.patch').write_bytes(subprocess.check_output(
        ['git','-c','core.filemode=false','diff','--binary'],cwd=PROJECT))
    parent=PROJECT/'pileup/experiments/alpha_ml_comparison_20261014'
    shutil.copyfile(parent/'calibration-settings.json',EXPERIMENT/'calibration-settings.json')
    refresh=PROJECT/'jupyter/artifacts/alpha_ml_refresh_20261021'
    records=load(parent/'data-split.json')['windows']
    for phase,role in [('common','train'),('confirmation','validation')]:
        for row in load(refresh/phase/'summary.json')['windows']:
            metadata=Path(row['waveform']).parent.parent/'metadata'/(Path(row['root']).stem+'.json')
            records.append(dict(label=f'{phase}_{row["seed"]}',root=row['root'],root_sha256=row['root_sha256'],
                seed=row['seed'],rate_hz=row['rate_hz'],split=role,waveforms=row['waveform'],
                truth=row['truth'],metadata=str(metadata),previous_role='Published test; now explicit development data'))
    identities=[r['root_sha256'] for r in records]
    assert len(identities)==len(set(identities))==75
    assert sum(r['split']=='train' for r in records)==49
    assert sum(r['split']=='validation' for r in records)==26
    for r in records:
        r['component_cache']=str((OUTPUT/'training-cache'/(Path(r['root']).stem+'.npz')).relative_to(PROJECT))
        r['input_hashes']={k:digest(PROJECT/r[k]) for k in ('waveforms','truth','metadata')}
    save(EXPERIMENT/'data-split.json',dict(windows=records,split_unit='Whole ROOT and all seeds before crops',
        train_count=49,validation_count=26,previous_tests='Previously inspected; no longer called independent tests'))
    plans=[]
    for rate,count in [(100000,8),(1000000,3),(10000000,2),(20000000,2),(40000000,2),
                       (60000000,2),(80000000,2),(100000000,2)]:
        for index in range(count):
            plans.append(dict(label=f'component_{rate:09d}_{index:02}',rate_hz=rate,
                              seed=262042000+len(plans),role='new_test',
                              preflight=rate in (1000000,10000000) and index==0))
    seeds={p['seed'] for p in plans}
    save(EXPERIMENT/'assets/test-plan.json',dict(planned_runs=plans))
    for p in (PROJECT/'pileup/experiments').glob('*/manifest.json'):
        d=load(p)
        assert not seeds.intersection(r.get('seed') for r in d.get('runs',[])+d.get('planned_runs',[])),p
    native_parent=PROJECT/'pileup/experiments/alpha_ml_refresh_20261021'
    old=load(native_parent/'manifest.json')
    shutil.copytree(native_parent/'source/g4',EXPERIMENT/'source/g4')
    shutil.copyfile(native_parent/'bin/gagg',EXPERIMENT/'bin/gagg')
    assert digest(EXPERIMENT/'bin/gagg')==old['binary_sha256']
    for name in ('baseline-settings.json','covariance.json','physical-selection.json'):
        shutil.copyfile(native_parent/'assets'/name,EXPERIMENT/'assets'/name)
    shutil.copytree(native_parent/'models',EXPERIMENT/'assets/reference-models')
    save(EXPERIMENT/'manifest.json',dict(created=datetime.now().astimezone().isoformat(),
        planned_runs=plans,runs=[],preflight_validated=False,max_workers=4,
        binary_sha256=old['binary_sha256'],environment=old['environment'],
        source_sha256={str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'source').rglob('*') if p.is_file()},
        source_model='Unchanged native shared-cell p-11B pure alpha, no gamma',
        window_us=20.,gate_us=[5.,15.],adequate_dynamic_range=True))
    save(OUTPUT/'protocol.json',dict(created=datetime.now().astimezone().isoformat(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        variants=list(VARIANTS),training_seeds=list(SEEDS),steps=2400,batch_size=16,
        primary_selection='Per model: validation macro charge-qualified F1, checkpoint every400 steps, threshold .15..95',
        secondary_selection='Validation extra-output budget from paired GRU-base selected checkpoint; no test threshold tuning',
        loss='Old BCE+2 logQ Huber+0.5 time Huber, plus component=1 and conditional-component-sum=0.1',
        component_target='Native ancestry voltage, original shared-cell history, linear electronics only',
        component_loss_support_ns=[-16.,1200.],component_loss_sampling_ns=2.,
        supervised_anchor_policy='Training only, nearest 2 ns bin; same-bin collisions excluded only from added loss',
        sum_loss='Sum of supervised central components only, not full observed voltage; weak/dark are not forced to zero',
        representation='Unchanged dense 3 heads and 8 ns peak separation; unknown count at inference',
        train_root_count=49,validation_root_count=26,fresh_test_root_count=len(plans),
        old_data_reuse='Old 23 tests now train, old 16 confirmations now validation; old 26/10 split preserved',
        new_amplitude_augmentation=False,contrastive_pretraining=False,
        fixed_emission_tau_ns=[95.,351.],no_truth_or_rate_at_inference=True,
        scoring='Q>=500 pe, max-cardinality time-only ±10 ns pairing, Q error<=20%',
        environment=dict(python=platform.python_version(),torch=torch.__version__,numpy=np.__version__),
        new_test_generation='Only after all8 models, thresholds, source and protocol are frozen'))
    print('PREPARED',len(records),'development ROOTs;',len(plans),'new test windows',flush=True)


def cache_one(record):
    settings=load(EXPERIMENT/'calibration-settings.json')
    dest=PROJECT/record['component_cache']
    result=native_component_cache(record,settings,dest)
    return dict(root=record['root'],seed=record['seed'],split=record['split'],rate_hz=record['rate_hz'],
                cache=record['component_cache'],cache_sha256=digest(dest),**result)


def cache(workers=2):
    if not 1<=workers<=4:raise ValueError('Use at most four cache processes')
    records=load(EXPERIMENT/'data-split.json')['windows']
    results=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(cache_one,r):r for r in records}
        for future in as_completed(futures):
            row=future.result();results.append(row)
            save(OUTPUT/'cache-progress.json',dict(completed=len(results),planned=len(records),windows=results))
            print('CACHE',len(results),row['split'],row['rate_hz'],row['events'],round(row['local_full_max_error_mV'],6),flush=True)
    save(OUTPUT/'cache-validation.json',dict(status='passed',windows=results,truth_used_at_inference=False))


def freeze():
    files={}
    for variant in VARIANTS:
        for seed in SEEDS:
            name=f'{variant}_{seed}'
            model,frozen=load_model(name)
            folder=OUTPUT/'models'/name
            for filename in ('model.pt','config.json','selection.json','frozen.json'):
                p=folder/filename;files[str(p.relative_to(PROJECT))]=digest(p)
    # Thresholds at a common validation extra-output budget are selected only
    # among the grid of the already selected checkpoint, never by new test data.
    matched={}
    for seed in SEEDS:
        base=load(OUTPUT/'models'/f'gru_base_{seed}'/'selection.json')['best']
        for variant in VARIANTS:
            name=f'{variant}_{seed}';d=load(OUTPUT/'models'/name/'selection.json')
            eligible=[c for c in d['grid'] if c['extra_count']<=base['extra_count']]
            if not eligible:
                chosen=min(d['grid'],key=lambda c:c['extra_count']);met=False
            else:
                chosen=max(eligible,key=lambda c:(c['objective'],c['threshold']));met=True
            matched[name]=dict(threshold=chosen['threshold'],validation_extra=chosen['extra_count'],
                               paired_base_budget=base['extra_count'],budget_met=met)
    save(OUTPUT/'matched-extra-thresholds.json',matched)
    for package in ('analysis','pileup/analysis'):
        for p in (PROJECT/package).glob('*.py'):
            q=EXPERIMENT/'inference-source'/p.relative_to(PROJECT)
            q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
    for package in ('analysis','pileup','pileup/analysis'):
        p=EXPERIMENT/'inference-source'/package/'__init__.py';p.parent.mkdir(parents=True,exist_ok=True);p.touch(exist_ok=True)
    for p in (EXPERIMENT/'inference-source').rglob('*.py'):
        files[str(p.relative_to(PROJECT))]=digest(p)
    for row in load(OUTPUT/'cache-validation.json')['windows']:
        assert digest(PROJECT/row['cache'])==row['cache_sha256']
    for p in [OUTPUT/'protocol.json',OUTPUT/'matched-extra-thresholds.json',OUTPUT/'cache-validation.json',EXPERIMENT/'data-split.json',
              EXPERIMENT/'calibration-settings.json',EXPERIMENT/'bin/gagg',EXPERIMENT/'replay.py',EXPERIMENT/'run_frozen.py',
              *[p for folder in ('assets','source') for p in (EXPERIMENT/folder).rglob('*') if p.is_file()]]:
        files[str(p.relative_to(PROJECT))]=digest(p)
    save(OUTPUT/'freeze.json',dict(frozen_at=datetime.now().astimezone().isoformat(),frozen_files=files,
                                 new_test_tuning=False,models=8))
    print('ALL_MODELS_FROZEN',len(files),flush=True)


def native(workers=2):
    from . import alpha_correlated_native as runner
    from .root_redundant_photon_fields import exact_redundant_fields
    assert load(EXPERIMENT/'manifest.json')['planned_runs']==load(EXPERIMENT/'assets/test-plan.json')['planned_runs']
    runner.PROJECT,runner.EXPERIMENT,runner.OUTPUT,runner.RECORD=PROJECT,EXPERIMENT,OUTPUT,RECORD
    original=runner.export_window
    def export(root,output,*,expected):
        try:
            return original(root,output,expected=expected)
        except Exception as error:
            if 'decompress' not in str(error).lower():raise
            checks=dict(original_error=repr(error))
            with exact_redundant_fields(original,checks):
                result=original(root,output,expected=expected)
            save(RECORD/(Path(root).stem+'-exact-fields.json'),checks)
            return result
    runner.export_window=export
    runner.run(workers)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','cache','train','freeze','native'])
    parser.add_argument('--variant',choices=VARIANTS)
    parser.add_argument('--seed',type=int,choices=SEEDS)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--steps',type=int,default=2400)
    args=parser.parse_args()
    if args.stage=='prepare':prepare()
    elif args.stage=='cache':cache(args.workers)
    elif args.stage=='train':train(args.variant,args.seed,args.steps,args.workers)
    elif args.stage=='freeze':freeze()
    else:native(args.workers)


if __name__=='__main__':main()
