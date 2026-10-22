"""Reproducible single-to-pair experiment with validation-only selection."""
from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
import platform
import shutil
import subprocess
import time

import numpy as np
import torch

from . import alpha_pair_data as data
from .alpha_pair_data import PROJECT,EXPERIMENT,OUTPUT,RECORD,load
from .alpha_pair_ml import ComponentNet,set_loss,blind_predict,tensors
from .alpha_rate_campaign import digest,save

SEEDS=(262050001,262050002)
SINGLE_VARIANTS=('single_logq','single_relative','single_fixed')
PAIR_VARIANTS=('pair_direct','pair_physical')
THRESHOLDS=(.2,.35,.5,.65,.8,.9)


def prepare():
    EXPERIMENT.mkdir(exist_ok=False)
    OUTPUT.mkdir(exist_ok=False)
    RECORD.mkdir(exist_ok=False)
    shutil.copyfile(PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb',RECORD/'notebook.before.ipynb')
    shutil.copyfile(PROJECT/'docs/literature-use.md',RECORD/'literature-use.before.txt')
    for folder in ('assets','models','figures'):
        (OUTPUT/folder).mkdir()
    for folder in ('assets','macros','logs','bin'):
        (EXPERIMENT/folder).mkdir()
    source=PROJECT/'pileup/experiments/alpha_component_ml_20261021/calibration-settings.json'
    config=load(source)
    # Fixed, train-only global reference; actual median pool Q is inserted after
    # pool construction. This provisional value is not used in any training.
    config['component_fixed_norm']=1.
    save(EXPERIMENT/'settings.json',config)
    records=load(PROJECT/'jupyter/artifacts/alpha_gamma_advanced_fit_20261020/input-records.json')
    grouped={name:[r for r in records[role] if r['particle']=='alpha']
             for name,role in (('train','calibration'),('validation','development'))}
    train_roots={r['root'] for r in grouped['train']}
    assert not train_roots.intersection(r['root'] for r in grouped['validation'])
    save(EXPERIMENT/'data-split.json',dict(**grouped,split_before_augmentation=True,
        roles='Existing isolated calibration/development only; previously viewed tests not independent confirmation'))
    save(OUTPUT/'protocol.json',dict(created=datetime.now().astimezone().isoformat(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        seeds=list(SEEDS),single_variants=list(SINGLE_VARIANTS),pair_variants=list(PAIR_VARIANTS),
        single_steps=1200,pair_steps=2400,batch_size=32,validation_every=200,
        fixed_scintillation_taus_ns=[95.,351.],samples=data.SAMPLES,sampling_ns=2.,
        task='Known-single calibration, then unknown0/1/2 event reconstruction in complete local 2048 ns frames',
        single_training='75% unscaled actual native isolated voltage; 25% matched-mean software controls',
        pair_training='After200-step empty/single warmup:50% joint shared-cell photon redigitizations,25% native-component linear controls,25% matched mean',
        pair_augmentation='Native donors independently scaled0.5..2 and linearly resampled; not a shared-cell digitization',
        native_train_roots=len(grouped['train']),native_validation_roots=len(grouped['validation']),
        component_loss='Compare no component, true-component relative norm and fixed train-derived global norm; one other setting held common',
        single_selection='Validation known-single joint Q20%/time10ns recovery, native75% + independent mean25%; tie lower native median absolute Q error',
        pair_selection='Validation charge-qualified F1 with no more than5% zero/single frames having an unmatched output; threshold grid fixed',
        pair_pretraining='Copy encoder/fine/time-head from same-seed validation-selected single model; global count/partition heads new',
        single_residual_q_bound='exp(0.25*tanh(residual)); observed positive charge integral plus learned-time finite-tail correction',
        pair_head_charge='Direct: component-mask integrals. Physical: global charge fractions. Both use learned-time mean tail correction.',
        assignment='Minimum of exactly2 permutations during training only; no true K/time/Q or rate at inference',
        final_target='Q>=500pe; time-only one-to-one10ns; thenQerror<=20%',
        fresh_test='24 new native isolated parents after source/models/threshold freeze; independent software mean and linear-pair controls generated afterward',
        further_gate='Extend unknownK>=3/high-rate training only if controlled pair separation improves without sacrificing native isolated calibration; otherwise document bottleneck first',
        max_cpu_threads=4,environment=dict(python=platform.python_version(),torch=torch.__version__,numpy=np.__version__),
        dynamic_range='Adequate analog voltage range; no ADC rejection/quantization',
        scope='Realized avalanche-equivalent pe; not deposited energy or material lifetime inference'))
    status=subprocess.check_output(['git','-c','core.fileMode=false','status','--short'],cwd=PROJECT,text=True)
    (RECORD/'git-status-at-prepare.txt').write_text(status)
    print('PREPARED',len(grouped['train']),len(grouped['validation']),flush=True)


def pools():
    split=load(EXPERIMENT/'data-split.json')
    for role in ('train','validation'):
        data.native_pool(split[role],OUTPUT/f'{role}-pool.npz')
    config=load(EXPERIMENT/'settings.json')
    pool=data.read_pool('train')
    eligible=pool['charge_pe']>=500.
    qref=float(np.median(pool['charge_pe'][eligible]))
    from .alpha_component_ml import PhysicalDecoder
    decoder=PhysicalDecoder(config)
    with torch.no_grad():
        reference=decoder.response(torch.arange(data.SAMPLES)*2.+.25-128.)*qref
    config['component_fixed_norm']=float(reference.square().sum())
    config['component_norm_reference_q_pe']=qref
    save(EXPERIMENT/'settings.json',config)
    vpool=data.read_pool('validation')
    single,parents=data.single_frames(vpool,seed=262050500)
    data.save_batch(OUTPUT/'assets/validation-native-single.npz',single,[dict(parent=i) for i in parents])
    sampler=data.FrameSampler(vpool,config)
    mean,info=sampler.sample(np.random.default_rng(262050501),128,single=True,mean_probability=1.)
    data.save_batch(OUTPUT/'assets/validation-mean-single.npz',mean,info)
    pair,info=sampler.sample(np.random.default_rng(262050502),512)
    data.save_batch(OUTPUT/'assets/validation-pair.npz',pair,info)
    save(OUTPUT/'pool-summary.json',dict(train=len(pool['charge_pe']),validation=len(vpool['charge_pe']),
        train_qualified=int(eligible.sum()),validation_qualified=int((vpool['charge_pe']>=500.).sum()),
        train_q_p0_p50_p100=np.quantile(pool['charge_pe'][eligible],[0,.5,1]).tolist(),
        fixed_reference_q_pe=qref,fixed_component_norm=config['component_fixed_norm']))
    print('POOLS_READY',qref,flush=True)


def batch_scores(fits,batch):
    from .alpha_pileup_fit import score
    rows=[]
    for i,fit in enumerate(fits):
        active=batch['active'][i].astype(bool)
        times=batch['times_ns'][i,active]
        truth=dict(first_active_ns=times,emission_ns=times,charge_pe=batch['charge_pe'][i,active])
        rows.append(score(fit,truth,[0.,data.SAMPLES*2.]))
    return rows


def aggregate(rows):
    names=('target_count','matched_count','recovered_count','missed_count','unmatched_predictions','merge_flags','split_flags','reconstructed_count')
    out={k:sum(r[k] for r in rows) for k in names}
    out['q_failed']=out['matched_count']-out['recovered_count']
    out['qualified_f1']=2*out['recovered_count']/max(1,out['target_count']+out['reconstructed_count'])
    for field,prefix in (('relative_charge_errors','q_error'),('time_errors_ns','t_error_ns')):
        a=np.array([v for r in rows for v in r[field]])
        out[prefix+'_p16_p50_p84']=np.quantile(a,[.16,.5,.84]).tolist() if len(a) else None
        out[prefix+'_median_abs']=float(np.median(abs(a))) if len(a) else None
    return out


def predict_batches(model,batch,threshold=.5):
    fits=[]
    for i in range(0,len(batch['voltage_mV']),64):
        part,_=blind_predict(model,batch['voltage_mV'][i:i+64],threshold=threshold)
        fits.extend(part)
    return fits


def single_validation(model):
    results={}
    for kind in ('native','mean'):
        batch=data.read_batch(OUTPUT/f'assets/validation-{kind}-single.npz')
        fits=predict_batches(model,batch)
        summary=aggregate(batch_scores(fits,batch))
        results[kind]=summary
    objective=.75*results['native']['recovered_count']/max(1,results['native']['target_count'])+.25*results['mean']['recovered_count']/max(1,results['mean']['target_count'])
    return dict(objective=objective,by_domain=results)


def pair_validation(model):
    batch=data.read_batch(OUTPUT/'assets/validation-pair.npz')
    shared=data.read_batch(OUTPUT/'validation-shared-pair.npz')
    batch={k:np.concatenate((batch[k],shared[k])) for k in batch}
    # Run each voltage once, and decode the fixed threshold grid from slots.
    raw=[]
    with torch.inference_mode():
        for i in range(0,len(batch['voltage_mV']),64):
            p=model(torch.from_numpy(batch['voltage_mV'][i:i+64]))
            raw.append({k:v.numpy() for k,v in p.items() if k in ('presence_logits','times_ns','charge_pe')})
    z={k:np.concatenate([r[k] for r in raw]) for k in raw[0]}
    probability=1/(1+np.exp(-np.clip(z['presence_logits'],-40,40)))
    choices=[]
    null=batch['active'].sum(1)==0
    simple=batch['active'].sum(1)<=1
    for threshold in THRESHOLDS:
        fits=[]
        for i in range(len(probability)):
            keep=(probability[i]>=threshold)&(z['charge_pe'][i]>=500.)
            fits.append(dict(times_ns=z['times_ns'][i,keep].astype(float).tolist(),charge_pe=z['charge_pe'][i,keep].astype(float).tolist()))
        scores=batch_scores(fits,batch)
        summary=aggregate(scores)
        null_fp=float(np.mean([bool(fits[i]['times_ns']) for i in np.flatnonzero(null)]))
        simple_extra=float(np.mean([scores[i]['unmatched_predictions']>0 for i in np.flatnonzero(simple)]))
        choices.append(dict(threshold=threshold,objective=summary['qualified_f1'],null_frame_fp_fraction=null_fp,
                            simple_frame_extra_fraction=simple_extra,summary=summary,eligible=simple_extra<=.05))
    eligible=[c for c in choices if c['eligible']]
    best=max(eligible or choices,key=lambda c:(c['objective'] if c['eligible'] else -c['null_frame_fp_fraction'],c['threshold']))
    return dict(**best,grid=choices)


def train(variant,seed,*,steps=None,threads=2,destination=None):
    if variant not in SINGLE_VARIANTS+PAIR_VARIANTS or seed not in SEEDS or not 1<=threads<=4:
        raise ValueError('Undeclared model/seed/CPU budget')
    torch.set_num_threads(threads);torch.manual_seed(seed)
    rng=np.random.default_rng(seed)
    single=variant.startswith('single_')
    steps=steps or (1200 if single else 2400)
    folder=Path(destination) if destination else OUTPUT/'models'/f'{variant}_{seed}'
    folder.mkdir(parents=True,exist_ok=False)
    config=load(EXPERIMENT/'settings.json')
    pool=data.read_pool('train');sampler=data.FrameSampler(pool,config)
    shared=data.read_batch(OUTPUT/'train-shared-pair.npz') if not single else None
    model=ComponentNet(config,slots=1 if single else 2,head='direct' if variant.endswith('direct') else 'physical')
    pretrained=None
    if not single:
        chosen=load(OUTPUT/'single-selection.json')['variant']
        pretrained=OUTPUT/'models'/f'{chosen}_{seed}'/'model.pt'
        state=torch.load(pretrained,map_location='cpu',weights_only=True)
        for name in ('encoder','fine'):
            getattr(model,name).load_state_dict({k[len(name)+1:]:v for k,v in state.items() if k.startswith(name+'.')})
        with torch.no_grad():
            model.time_head.weight.copy_(state['time_head.weight'].repeat(2,1,1))
            model.time_head.bias.copy_(state['time_head.bias'].repeat(2))
        normalization=load(OUTPUT/'single-selection.json')['pair_component_normalization']
    else:
        normalization='relative' if variant.endswith('relative') else 'fixed'
    weight=0. if variant=='single_logq' else 1.
    run_config=dict(variant=variant,seed=seed,steps=steps,threads=threads,batch_size=32,
        parameters=sum(p.numel() for p in model.parameters()),normalization=normalization,component_weight=weight,
        single_known_K=single,pretrained=str(pretrained.relative_to(PROJECT)) if pretrained else None,
        pretrained_sha256=digest(pretrained) if pretrained else None,
        settings_sha256=digest(EXPERIMENT/'settings.json'),train_pool_sha256=digest(OUTPUT/'train-pool.npz'),
        source_sha256={str(p.relative_to(PROJECT)):digest(p) for p in (PROJECT/'pileup/analysis').glob('alpha_pair_*.py')})
    save(folder/'config.json',run_config)
    optimizer=torch.optim.AdamW(model.parameters(),lr=.001,weight_decay=.0001)
    history=[];best=None;started=time.monotonic()
    print('TRAIN_BEGIN',variant,seed,run_config['parameters'],flush=True)
    for step in range(1,steps+1):
        batch,_=sampler.sample(rng,32,single=single,mean_probability=.25 if single else .5)
        if not single and step<=200:
            # Start the count head on empty/single examples; labels remain
            # unchanged. This is a declared project curriculum, not a paper claim.
            batch,_=sampler.sample(rng,32,single=False,forced={'count':int(rng.integers(0,2)),'mean':False})
        if not single and step>200:
            ids=rng.integers(0,len(shared['active']),16)
            batch={k:np.concatenate((v[:16],shared[k][ids])) for k,v in batch.items()}
        target=tensors({k:v for k,v in batch.items() if k!='voltage_mV'})
        model.train();optimizer.zero_grad(set_to_none=True)
        pred=model(torch.from_numpy(batch['voltage_mV']))
        warmup=min(1.,max(0.,(step-100)/200.))
        loss,parts=set_loss(pred,target,config,normalization=normalization,
                            component_weight=weight*warmup,single=single)
        if not torch.isfinite(loss):
            raise FloatingPointError('Nonfinite loss')
        loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step()
        if step==1 or step%200==0 or step==steps:
            model.eval();selection=single_validation(model) if single else pair_validation(model)
            row=dict(step=step,loss=float(loss.detach()),parts={k:float(v.detach()) for k,v in parts.items()},
                     objective=selection['objective'],elapsed_s=time.monotonic()-started)
            history.append(row);save(folder/'history.json',history)
            native_error=selection['by_domain']['native']['q_error_median_abs'] if single else None
            key=(selection['objective'],-(native_error if native_error is not None else 1e6)) if single else (selection['objective'] if selection['eligible'] else -1.,selection['threshold'])
            if best is None or key>best:
                best=key;torch.save(model.state_dict(),folder/'model.pt')
                save(folder/'selection.json',dict(step=step,**selection))
            print('TRAIN_STEP',variant,seed,step,round(row['loss'],5),round(selection['objective'],4),round(row['elapsed_s'],1),flush=True)
    save(folder/'frozen.json',dict(frozen_at=datetime.now().astimezone().isoformat(),
        weights_sha256=digest(folder/'model.pt'),config_sha256=digest(folder/'config.json'),training_seconds=time.monotonic()-started))


def select_single():
    rows=[]
    for variant in SINGLE_VARIANTS:
        selections=[load(OUTPUT/'models'/f'{variant}_{s}'/'selection.json') for s in SEEDS]
        objective=float(np.mean([s['objective'] for s in selections]))
        rows.append(dict(variant=variant,objective=objective,selections=selections))
    choice=max(rows,key=lambda r:r['objective'])
    norm='relative' if choice['variant']=='single_relative' else 'fixed'
    save(OUTPUT/'single-selection.json',dict(variant=choice['variant'],pair_component_normalization=norm,choices=rows,
        selection='Validation only; mean objective over two training seeds'))
    print('SINGLE_SELECTED',choice['variant'],choice['objective'],flush=True)


def load_model(name):
    folder=OUTPUT/'models'/name
    cfg=load(folder/'config.json');frozen=load(folder/'frozen.json')
    assert digest(folder/'model.pt')==frozen['weights_sha256']
    assert digest(folder/'config.json')==frozen['config_sha256']
    single=cfg['single_known_K']
    model=ComponentNet(load(EXPERIMENT/'settings.json'),slots=1 if single else 2,
                       head='direct' if cfg['variant'].endswith('direct') else 'physical')
    model.load_state_dict(torch.load(folder/'model.pt',map_location='cpu',weights_only=True))
    return model.eval(),cfg,load(folder/'selection.json')


def sanity():
    torch.set_num_threads(2);torch.manual_seed(262050700)
    cfg=load(EXPERIMENT/'settings.json')
    sampler=data.FrameSampler(data.read_pool('train'),cfg)
    rng=np.random.default_rng(262050700)
    pieces=[]
    for count,gap in ((0,200.),(1,200.),(2,200.),(2,40.)):
        b,_=sampler.sample(rng,4,forced=dict(count=count,mean=True,gap_ns=gap,q_pe=4000.,ratio=1.))
        pieces.append(b)
    batch={k:np.concatenate([p[k] for p in pieces]) for k in pieces[0]}
    data.save_batch(OUTPUT/'assets/sanity-train.npz',batch,[dict(training_only=True)]*len(batch['active']))
    model=ComponentNet(cfg,slots=2,head='direct')
    optimizer=torch.optim.AdamW(model.parameters(),lr=.002)
    target=tensors({k:v for k,v in batch.items() if k!='voltage_mV'})
    for step in range(800):
        optimizer.zero_grad(set_to_none=True)
        pred=model(torch.from_numpy(batch['voltage_mV']))
        loss,_=set_loss(pred,target,cfg,component_weight=1.)
        loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step()
    fits=predict_batches(model.eval(),batch)
    rows=batch_scores(fits,batch)
    summary=aggregate(rows)
    save(OUTPUT/'sanity.json',dict(summary=summary,rows=rows,scope='Overfit fixed16 training cases; not independent accuracy',
        loss=float(loss.detach()),steps=800,passed=summary['recovered_count']==summary['target_count'] and summary['unmatched_predictions']==0))
    torch.save(model.state_dict(),OUTPUT/'assets/sanity-model.pt')
    print('SANITY',summary,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=('prepare','pools','sanity','train','select-single'))
    p.add_argument('--variant',choices=SINGLE_VARIANTS+PAIR_VARIANTS)
    p.add_argument('--seed',type=int,choices=SEEDS)
    p.add_argument('--steps',type=int)
    p.add_argument('--threads',type=int,default=2)
    p.add_argument('--destination',type=Path)
    a=p.parse_args()
    if a.stage=='prepare':prepare()
    elif a.stage=='pools':pools()
    elif a.stage=='sanity':sanity()
    elif a.stage=='select-single':select_single()
    else:train(a.variant,a.seed,steps=a.steps,threads=a.threads,destination=a.destination)


if __name__=='__main__':main()
