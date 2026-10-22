"""Independent tests for calibration, unknown0/1/2 counts and actual components."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import time

import numpy as np
from scipy.signal import lfilter
import torch

from . import alpha_pair_data as data
from .alpha_pair_data import PROJECT,EXPERIMENT,OUTPUT,load
from .alpha_pair_registry import ALL_PAIRS,load_model as registry_load,predict as registry_predict
from .alpha_pair_fitting import PairFitter
from .alpha_pair_study import (SEEDS,SINGLE_VARIANTS,PAIR_VARIANTS,load_model,
                              batch_scores,aggregate,predict_batches)
from .alpha_pair_native import check_freeze
from .alpha_rate_campaign import digest,save
from .alpha_pileup_fit import score

PENALTIES=(.0001,.0003,.001,.003,.01,.03,.1,.3,1.,3.)


def pad_slots(batch):
    result=dict(batch)
    for k in ('active','charge_pe','times_ns'):
        result[k]=np.pad(batch[k],((0,0),(0,1)),constant_values=1. if k=='charge_pe' else 0.)
    result['components_mV']=np.pad(batch['components_mV'],((0,0),(0,1),(0,0)))
    return result


def validation_batch():
    linear=data.read_batch(OUTPUT/'assets/validation-pair.npz')
    shared=data.read_batch(OUTPUT/'validation-shared-pair.npz')
    return {k:np.concatenate((linear[k],shared[k])) for k in linear}


def fit_validation(workers=2):
    fitter=PairFitter(load(EXPERIMENT/'settings.json'))
    batch=validation_batch();destination=OUTPUT/'fit-validation-candidates.json'
    if destination.exists():
        rows=load(destination)['rows']
    else:
        rows=[None]*len(batch['active'])
        def one(i):
            c=fitter.candidates(batch['voltage_mV'][i])
            return i,[{k:v for k,v in x.items() if k!='model_mV'} for x in c]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for f in as_completed([pool.submit(one,i) for i in range(len(rows))]):
                i,rows[i]=f.result()
                done=sum(r is not None for r in rows)
                if done%64==0:print('FIT_VALIDATION',done,len(rows),flush=True)
        save(destination,dict(rows=rows,settings_sha256=digest(EXPERIMENT/'settings.json')))
    choices=[]
    simple=batch['active'].sum(1)<=1
    for penalty in PENALTIES:
        fits=[fitter.choose(r,penalty) for r in rows]
        scores=batch_scores(fits,batch);summary=aggregate(scores)
        extras=float(np.mean([scores[i]['unmatched_predictions']>0 for i in np.flatnonzero(simple)]))
        choices.append(dict(penalty=penalty,objective=summary['qualified_f1'],simple_frame_extra_fraction=extras,
                            eligible=extras<=.05,summary=summary))
    eligible=[c for c in choices if c['eligible']]
    chosen=max(eligible or choices,key=lambda c:(c['objective'] if c['eligible'] else -c['simple_frame_extra_fraction'],-c['penalty']))
    save(OUTPUT/'fitting-selection.json',dict(**chosen,grid=choices,
        selection='Validation only; same full local frame, unknownK in0/1/2; fixed tau; no truth at fit input',
        candidates_sha256=digest(destination)))
    print('FIT_SELECTED',chosen,flush=True)


def make_mean_test():
    rng=np.random.default_rng(262054500)
    sampler=data.FrameSampler(data.read_pool('test'),load(EXPERIMENT/'settings.json'))
    cases=[]
    for noise in ('white','AR1'):
        for phase in (0,1):
            for gap in data.GAPS:
                for ratio in data.RATIOS:
                    cases.append(dict(count=2,gap_ns=gap,ratio=ratio,q_pe=4000.,noise=noise,phase=phase))
            for q in (500.,1000.,4000.,8000.,16000.,40000.):
                cases.append(dict(count=1,q_pe=q,noise=noise,phase=phase))
            for j in range(8):cases.append(dict(count=0,noise=noise,phase=phase))
    values={k:[] for k in ('voltage_mV','active','charge_pe','times_ns','components_mV')}
    t=np.arange(data.SAMPLES)*2.+.25
    for c in cases:
        curves=np.zeros((2,data.SAMPLES));q=np.ones(2);ts=np.full(2,128.+.37*c['phase']);active=np.zeros(2)
        for j in range(c['count']):
            ts[j]=128.+.37*c['phase']+j*c.get('gap_ns',0.)
            q[j]=c['q_pe']*(c.get('ratio',1.) if j else 1.)
            curves[j]=sampler.mean_curve(q[j],ts[j]);active[j]=1.
        eps=rng.normal(0.,.1,data.SAMPLES)
        if c['noise']=='AR1':eps=lfilter([np.sqrt(1-.98**2)],[1.,-.98],eps)
        for key,value in zip(values,(10.+curves.sum(0)+eps,active,q,ts,curves)):
            values[key].append(value)
    return {k:np.asarray(v,np.float32) for k,v in values.items()},cases


def test_inputs():
    check_freeze()
    pool=data.read_pool('test')
    native,parents=data.single_frames(pool,seed=262054001,include_weak=True)
    data.save_batch(OUTPUT/'assets/test-native.npz',pad_slots(native),[dict(parent=i,domain='native-original-isolated') for i in parents])
    native_single,parents=data.single_frames(pool,seed=262054001,include_weak=False)
    data.save_batch(OUTPUT/'assets/test-native-single.npz',native_single,[dict(parent=i,known_single=True) for i in parents])
    mean,info=make_mean_test();data.save_batch(OUTPUT/'assets/test-mean.npz',mean,info)
    linear,info=data.FrameSampler(pool,load(EXPERIMENT/'settings.json')).sample(np.random.default_rng(262054002),256,mean_probability=0.)
    data.save_batch(OUTPUT/'assets/test-linear.npz',linear,info)
    from .alpha_pair_photons import create
    records=load(OUTPUT/'native-records.json')['test_single']
    create('test',count=256,seed=262054000,workers=2,records=records)
    shutil_module=__import__('shutil')
    shutil_module.copyfile(OUTPUT/'test-shared-pair.npz',OUTPUT/'assets/test-shared.npz')
    shutil_module.copyfile(OUTPUT/'test-shared-pair.json',OUTPUT/'assets/test-shared.json')
    save(OUTPUT/'test-inputs.json',dict(datasets={kind:dict(path=str((OUTPUT/f'assets/test-{kind}.npz').relative_to(PROJECT)),
         sha256=digest(OUTPUT/f'assets/test-{kind}.npz'),metadata_sha256=digest(OUTPUT/f'assets/test-{kind}.json'))
         for kind in ('native','native-single','mean','linear','shared')},
         root_parents=24,new_independent_transport_parents=True,after_freeze=True,
         independence='Repeated pair cases share24 new optical transports; report separately from24 native ROOTs'))
    print('TEST_INPUTS_READY',flush=True)


def evaluate_dataset(kind,*,workers=2):
    check_freeze();torch.set_num_threads(1)
    batch=data.read_batch(OUTPUT/f'assets/test-{kind}.npz')
    folder=OUTPUT/'test'/kind;folder.mkdir(parents=True,exist_ok=False)
    results=[dict(case=i,truth_slot_ids=np.flatnonzero(batch['active'][i]).tolist(),methods={}) for i in range(len(batch['active']))]
    models={}
    for variant in ALL_PAIRS:
        for seed in SEEDS:
            name=f'{variant}_{seed}';model,cfg,selected=registry_load(name)
            models[name]=(model,selected['threshold'])
    fitter=PairFitter(load(EXPERIMENT/'settings.json'))
    penalty=load(OUTPUT/'fitting-selection.json')['penalty']
    for name,(model,threshold) in models.items():
        outputs=[];fits=[];started=time.monotonic()
        for i in range(0,len(batch['active']),64):
            fit,arrays=registry_predict(model,batch['voltage_mV'][i:i+64],threshold)
            fits.extend(fit);outputs.append(arrays)
        arrays={k:np.concatenate([o[k] for o in outputs]) for k in outputs[0]}
        np.savez_compressed(folder/f'{name}-raw.npz',**arrays)
        seconds=(time.monotonic()-started)/len(fits)
        raw_scores=batch_scores(fits,batch)
        refined_curves=[]
        for i,(fit,s) in enumerate(zip(fits,raw_scores)):
            results[i]['methods'][name]=dict(fit=fit,score=s,seconds=seconds)
            tick=time.monotonic()
            refined=fitter.refine(batch['voltage_mV'][i],fit['times_ns'],radius=12.)
            curve=fitter.columns(np.arange(data.SAMPLES)*2.+.25,refined['times_ns']).T*np.asarray(refined['charge_pe'])[:,None]
            refined_curves.append(np.pad(curve,((0,2-len(curve)),(0,0))))
            rf={k:refined[k] for k in ('times_ns','charge_pe')}
            rs=batch_scores([rf],{k:v[i:i+1] for k,v in batch.items()})[0]
            results[i]['methods'][name+'_refit']=dict(fit=rf,score=rs,seconds=seconds+time.monotonic()-tick)
        np.savez_compressed(folder/f'{name}-refit.npz',components_mV=np.asarray(refined_curves,np.float32))
        print('TEST_MODEL',kind,name,aggregate(raw_scores),flush=True)
    def physical(i):
        tick=time.monotonic();candidates=fitter.candidates(batch['voltage_mV'][i]);elapsed=time.monotonic()-tick
        fit=fitter.choose(candidates,penalty)
        count=int(batch['active'][i].sum())
        # This separate diagnostic receives the true count, and is never mixed
        # into the blind main-method table or used to tune any test threshold.
        oracle=candidates[min(count,len(candidates)-1)]
        oracle={k:oracle[k] for k in ('times_ns','charge_pe')}
        return i,fit,oracle,elapsed,[{k:v for k,v in c.items() if k!='model_mV'} for c in candidates]
    candidates=[None]*len(results)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for f in as_completed([pool.submit(physical,i) for i in range(len(results))]):
            i,fit,oracle,elapsed,candidates[i]=f.result()
            for name,pred in (('enumerated_fit',fit),('known_K_fit_diagnostic',oracle)):
                s=batch_scores([pred],{k:v[i:i+1] for k,v in batch.items()})[0]
                results[i]['methods'][name]=dict(fit=pred,score=s,seconds=elapsed,
                    truth_inputs=['K'] if name.startswith('known_K') else [])
    # Strict pair tolerance is evaluation only; it prevents two times near the
    # midpoint from being called resolved merely because the default10ns gate
    # is wider than the true separation.
    for i,row in enumerate(results):
        active=batch['active'][i].astype(bool);times=batch['times_ns'][i,active]
        truth=dict(first_active_ns=times,emission_ns=times,charge_pe=batch['charge_pe'][i,active])
        tolerance=min(10.,float(np.diff(np.sort(times))[0])/4.) if len(times)==2 else 10.
        for name,result in row['methods'].items():
            result['strict_score']=score(result['fit'],truth,[0.,data.SAMPLES*2.],tolerance_ns=tolerance)
    save(folder/'results.json',dict(kind=kind,rows=results,candidates=candidates,
        input_sha256=digest(OUTPUT/f'assets/test-{kind}.npz'),after_freeze=True,
        inference_truth_inputs=[],oracle_separate=True,strict_pair_time_tolerance='min(10ns,true gap/4); evaluation only'))
    print('TEST_DATASET_DONE',kind,len(results),flush=True)


def single_test():
    check_freeze();torch.set_num_threads(1)
    batch=data.read_batch(OUTPUT/'assets/test-native-single.npz')
    rows=[]
    for variant in SINGLE_VARIANTS:
        for seed in SEEDS:
            name=f'{variant}_{seed}';model,_,_=load_model(name)
            fits=predict_batches(model,batch);scores=batch_scores(fits,batch)
            rows.append(dict(name=name,summary=aggregate(scores),fits=fits,scores=scores,
                known_K=1,targets=len(batch['active']),independent_parent_roots=len(batch['active']),
                scope='Qualified original native isolated roots only; not unknown-K pileup'))
    save(OUTPUT/'single-test.json',dict(rows=rows,selected_variant=load(OUTPUT/'single-selection.json')['variant']))
    print('SINGLE_TEST_DONE',[(r['name'],r['summary']['recovered_count'],r['summary']['target_count']) for r in rows],flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=('fit-validation','test-inputs','evaluate','single-test'))
    p.add_argument('--kind',choices=('native','mean','linear','shared'))
    p.add_argument('--workers',type=int,default=2)
    a=p.parse_args()
    if a.stage=='fit-validation':fit_validation(a.workers)
    elif a.stage=='test-inputs':test_inputs()
    elif a.stage=='single-test':single_test()
    else:evaluate_dataset(a.kind,workers=a.workers)


if __name__=='__main__':main()
