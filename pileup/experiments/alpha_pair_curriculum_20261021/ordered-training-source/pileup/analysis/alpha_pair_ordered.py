"""Ordered components, an explicit count classifier and global conditioning.

An additional development ablation after the unordered-slot validation failure;
no held-out test has been generated. Its changes are evaluated as a combined
method, not attributed solely to any one of count/order/global conditioning.
"""
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import time

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

from . import alpha_pair_data as data
from .alpha_pair_data import PROJECT,EXPERIMENT,OUTPUT,load
from .alpha_pair_ml import ComponentNet,tensors
from .alpha_pair_study import SEEDS,batch_scores,aggregate
from .alpha_rate_campaign import digest,save

VARIANTS=('pair_ordered_direct','pair_ordered_physical')
THRESHOLDS=(.2,.35,.5,.65,.8,.9,.95,.99)


class OrderedNet(ComponentNet):
    def __init__(self,settings,*,head):
        super().__init__(settings,slots=2,head=head)
        self.condition=nn.Linear(32*16+1,48)
        nn.init.zeros_(self.condition.weight);nn.init.zeros_(self.condition.bias)
        self.count_head=nn.Sequential(nn.Linear(32*16+1,64),nn.SiLU(),nn.Linear(64,3))

    def forward(self,voltage):
        y=voltage-self.baseline;positive=y.clamp_min(0.)
        derivative=F.pad((y[:,2:]-y[:,:-2])/4.,(1,1))
        smooth=F.avg_pool1d(y[:,None],33,stride=1,padding=16)[:,0]
        coord=(self.grid/self.grid[-1]*2-1).expand_as(y)
        x=torch.stack((torch.asinh(y/500.),torch.asinh(derivative/10.),torch.asinh(smooth/500.),coord),1)
        encoded=self.encoder(x)
        area=positive.sum(-1)*2./self.area_per_pe
        pooled=torch.cat((F.adaptive_avg_pool1d(encoded,16).flatten(1),torch.log1p(area/3000.)[:,None]),1)
        fine=self.fine(torch.cat((F.interpolate(encoded,size=self.samples,mode='linear',align_corners=False),x),1))
        scale,shift=self.condition(pooled).chunk(2,1)
        fine=fine*(1.+.1*torch.tanh(scale[:,:,None]))+.1*shift[:,:,None]
        time_logits=self.time_head(fine)
        times=(F.softmax(time_logits,-1)*self.grid).sum(-1)
        unit=self.decoder.response(self.grid[None,None,:]-times[:,:,None])
        fraction=(unit.sum(-1)*2./self.area_per_pe).clamp_min(.02)
        if self.head=='direct':
            curves=F.softmax(self.mask_head(fine),1)[:,:2]*positive[:,None]
            q=curves.sum(-1)*2./self.area_per_pe/fraction
        else:
            z=self.global_head(pooled)
            weights=F.softmax(z[:,2:5],1)
            q=area[:,None]*weights[:,:2]/fraction
            curves=q[:,:,None]*unit
        count_logits=self.count_head(pooled)
        cp=F.softmax(count_logits,-1)
        existence=torch.stack((cp[:,1]+cp[:,2],cp[:,2]),1).clamp(1e-6,1-1e-6)
        return dict(count_logits=count_logits,presence_logits=torch.logit(existence),time_logits=time_logits,
            times_ns=times,charge_pe=q.clamp_min(.01),components_mV=curves,area_fraction=fraction)


def ordered_loss(pred,target,settings,weight=1.):
    # Qualified targets first, then in time order. Unqualified physical signals
    # stay in the observed waveform/nuisance channel, never become zero input.
    order=torch.argsort(torch.where(target['active'].bool(),target['times_ns'],torch.full_like(target['times_ns'],1e6)),1)
    a=torch.gather(target['active'],1,order)
    q=torch.gather(target['charge_pe'],1,order)
    t=torch.gather(target['times_ns'],1,order)
    curves=torch.gather(target['components_mV'],1,order[:,:,None].expand_as(target['components_mV']))
    qloss=(F.smooth_l1_loss(torch.log(pred['charge_pe']/q.clamp_min(1.)),torch.zeros_like(q),reduction='none')*a).sum(1)
    grid=torch.arange(pred['time_logits'].shape[-1])*2.+.25
    density=torch.exp(-.5*((grid[None,None]-t[:,:,None])/2.)**2)
    density=density/density.sum(-1,keepdim=True).clamp_min(1e-20)
    tloss=(-(density*F.log_softmax(pred['time_logits'],-1)).sum(-1)*a).sum(1)
    closs=(pred['components_mV']-curves).square().sum((1,2))/settings['component_fixed_norm']
    count=F.cross_entropy(pred['count_logits'],a.sum(1).long(),reduction='none')
    loss=(count+6.*qloss+.5*tloss+weight*closs).mean()
    return loss,dict(count=count.mean(),logq=qloss.mean(),time_ce=tloss.mean(),component=closs.mean())


def decode(arrays,threshold):
    z={k:np.array(v,copy=True) for k,v in arrays.items()}
    logits=z['count_logits'];prob=np.exp(logits-logits.max(1,keepdims=True));prob/=prob.sum(1,keepdims=True)
    rawq=z.get('raw_charge_pe',z['charge_pe']).copy()
    rawt=z.get('raw_times_ns',z['times_ns']).copy()
    rawc=z.get('raw_components_mV',z['components_mV']).copy()
    z.update(raw_charge_pe=rawq,raw_times_ns=rawt,raw_components_mV=rawc)
    z['charge_pe']=rawq.copy();z['times_ns']=rawt.copy();z['components_mV']=rawc.copy()
    fits=[]
    for i,p in enumerate(prob):
        count=2 if p[2]>=threshold else 1 if p[1]+p[2]>=.5 else 0
        if count==1:
            # Conditional single hypothesis uses the whole learned signal,
            # so uncertainty between slots does not halve a real single Q.
            z['charge_pe'][i,0]=float(np.sum(rawq[i]*z['area_fraction'][i])/z['area_fraction'][i,0])
            z['components_mV'][i,0]=rawc[i].sum(0);z['components_mV'][i,1]=0.
        ids=np.arange(count)
        ids=ids[z['charge_pe'][i,ids]>=500.]
        ids=ids[np.argsort(z['times_ns'][i,ids])]
        fits.append(dict(status='predicted',times_ns=z['times_ns'][i,ids].astype(float).tolist(),
            charge_pe=z['charge_pe'][i,ids].astype(float).tolist(),slot_ids=ids.tolist(),predicted_K=count,
            count_probabilities=p.astype(float).tolist(),known_single=False,truth_inputs=[]))
    return fits,z


def predict(model,voltage,threshold):
    y=np.asarray(voltage,np.float32);one=y.ndim==1
    if one:y=y[None]
    with torch.inference_mode():z={k:v.numpy() for k,v in model(torch.from_numpy(y)).items()}
    fits,arrays=decode(z,threshold)
    return (fits[0] if one else fits),arrays


def validation(model):
    from .alpha_pair_evaluate import validation_batch
    batch=validation_batch();chunks=[]
    with torch.inference_mode():
        for i in range(0,len(batch['active']),64):chunks.append({k:v.numpy() for k,v in model(torch.from_numpy(batch['voltage_mV'][i:i+64])).items()})
    z={k:np.concatenate([x[k] for x in chunks]) for k in chunks[0]}
    simple=batch['active'].sum(1)<=1;choices=[]
    for threshold in THRESHOLDS:
        fits,_=decode(z,threshold);scores=batch_scores(fits,batch);summary=aggregate(scores)
        extra=float(np.mean([scores[i]['unmatched_predictions']>0 for i in np.flatnonzero(simple)]))
        choices.append(dict(threshold=threshold,objective=summary['qualified_f1'],simple_frame_extra_fraction=extra,eligible=extra<=.05,summary=summary))
    eligible=[x for x in choices if x['eligible']]
    chosen=max(eligible or choices,key=lambda x:(x['objective'] if x['eligible'] else -x['simple_frame_extra_fraction'],x['threshold']))
    return dict(**chosen,grid=choices)


def train(variant,seed,steps=2400):
    torch.set_num_threads(1);torch.manual_seed(seed);rng=np.random.default_rng(seed)
    settings=load(EXPERIMENT/'settings.json');folder=OUTPUT/'models'/f'{variant}_{seed}'
    folder.mkdir(exist_ok=False)
    model=OrderedNet(settings,head='direct' if variant.endswith('direct') else 'physical')
    chosen=load(OUTPUT/'single-selection.json')['variant'];pre=OUTPUT/'models'/f'{chosen}_{seed}'/'model.pt'
    state=torch.load(pre,map_location='cpu',weights_only=True)
    for name in ('encoder','fine'):
        getattr(model,name).load_state_dict({k[len(name)+1:]:v for k,v in state.items() if k.startswith(name+'.')})
    with torch.no_grad():
        model.time_head.weight.copy_(state['time_head.weight'].repeat(2,1,1));model.time_head.bias.copy_(state['time_head.bias'].repeat(2))
    config=dict(variant=variant,seed=seed,steps=steps,threads=1,batch_size=32,single_known_K=False,
        parameters=sum(p.numel() for p in model.parameters()),pretrained=str(pre.relative_to(PROJECT)),pretrained_sha256=digest(pre),
        source_sha256={str(p.relative_to(PROJECT)):digest(p) for p in (Path(__file__),PROJECT/'pileup/analysis/alpha_pair_ml.py',PROJECT/'pileup/analysis/alpha_pair_data.py')},
        model='Ordered targets + explicit0/1/2 classifier + global decoder conditioning + conditional single merge; combined ablation',
        reason='Development unordered-slot failure before any new test; not test-tuned')
    save(folder/'config.json',config)
    sampler=data.FrameSampler(data.read_pool('train'),settings);shared=data.read_batch(OUTPUT/'train-shared-pair.npz')
    optimizer=torch.optim.AdamW(model.parameters(),lr=.001,weight_decay=.0001)
    history=[];best=-1.;started=time.monotonic()
    for step in range(1,steps+1):
        batch,_=sampler.sample(rng,32,mean_probability=.5)
        if step<=200:batch,_=sampler.sample(rng,32,forced={'count':int(rng.integers(0,2)),'mean':False})
        else:
            ids=rng.integers(0,len(shared['active']),16);batch={k:np.concatenate((v[:16],shared[k][ids])) for k,v in batch.items()}
        target=tensors({k:v for k,v in batch.items() if k!='voltage_mV'})
        optimizer.zero_grad(set_to_none=True);model.train()
        loss,parts=ordered_loss(model(torch.from_numpy(batch['voltage_mV'])),target,settings,min(1.,max(0.,(step-100)/200.)))
        assert torch.isfinite(loss)
        loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step()
        if step==1 or step%200==0 or step==steps:
            model.eval();choice=validation(model)
            row=dict(step=step,loss=float(loss.detach()),parts={k:float(v.detach()) for k,v in parts.items()},objective=choice['objective'],elapsed_s=time.monotonic()-started)
            history.append(row);save(folder/'history.json',history)
            if choice['eligible'] and choice['objective']>best:
                best=choice['objective'];torch.save(model.state_dict(),folder/'model.pt');save(folder/'selection.json',dict(step=step,**choice))
            print('ORDERED_STEP',variant,seed,step,round(choice['objective'],4),round(row['elapsed_s'],1),flush=True)
    save(folder/'frozen.json',dict(frozen_at=datetime.now().astimezone().isoformat(),weights_sha256=digest(folder/'model.pt'),config_sha256=digest(folder/'config.json'),training_seconds=time.monotonic()-started))


def load_model(name):
    folder=OUTPUT/'models'/name;config=load(folder/'config.json');frozen=load(folder/'frozen.json')
    assert digest(folder/'model.pt')==frozen['weights_sha256'] and digest(folder/'config.json')==frozen['config_sha256']
    model=OrderedNet(load(EXPERIMENT/'settings.json'),head='direct' if config['variant'].endswith('direct') else 'physical')
    model.load_state_dict(torch.load(folder/'model.pt',map_location='cpu',weights_only=True))
    return model.eval(),config,load(folder/'selection.json')


def main():
    p=argparse.ArgumentParser();p.add_argument('--variant',choices=VARIANTS,required=True);p.add_argument('--seed',type=int,choices=SEEDS,required=True)
    a=p.parse_args();train(a.variant,a.seed)


if __name__=='__main__':main()
