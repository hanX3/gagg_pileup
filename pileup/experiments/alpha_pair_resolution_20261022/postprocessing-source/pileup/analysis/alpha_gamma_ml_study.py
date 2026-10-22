"""Grouped supervised training and blind mixed-event neural comparisons."""
from __future__ import annotations
import argparse
from collections import defaultdict
from dataclasses import asdict
from datetime import datetime
import importlib.metadata
import json
from pathlib import Path
import platform
import time
import numpy as np
import torch
from .alpha_rate_campaign import save,digest
from .alpha_gamma_ml import (MixedEventNet,NetworkConfig,voltage_features,
    supervised_targets,loss_function,predict_features,decode)
from .alpha_gamma_pileup_fit import score_mixed


def project_root():
    return next(p for p in [Path.cwd(),*Path.cwd().parents]
                if (p/'g4/data').is_dir() and (p/'README.md').is_file())


def summarize_scores(records):
    """Use separate truth and prediction gate denominators, including guards."""
    groups=defaultdict(lambda:dict(target=0,recovered=0,predicted=0,precise=0))
    for r in records:
        score,fit,truth,gate=r['score'],r['fit'],r['truth'],r['gate_ns']
        for particle in ('alpha','gamma'):
            key=(r['row']['gamma_model'],r['row']['rate_hz'],particle)
            g=groups[key];s=score['by_particle'][particle]
            g['target']+=s['target_count'];g['recovered']+=s['recovered_count']
            mask=np.array([name==particle and gate[0]<=t<gate[1]
                for name,t in zip(fit['particle_name'],fit['times_ns'])],bool)
            g['predicted']+=int(mask.sum())
            g['precise']+=int(sum(bool(mask[j]) and truth['particle_name'][i]==particle and
                abs(fit['charge_pe'][j]/truth['charge_pe'][i]-1)<=.2 for i,j in score['all_matched_pairs']))
    rows=[]
    for (model,rate,particle),g in sorted(groups.items()):
        recall=g['recovered']/g['target'] if g['target'] else 0.
        precision=g['precise']/g['predicted'] if g['predicted'] else 0.
        f1=2*recall*precision/(recall+precision) if recall+precision else 0.
        rows.append(dict(gamma_model=model,rate_hz=rate,particle=particle,**g,
                         recall=recall,precision=precision,full_f1=f1))
    populated=[g['full_f1'] for g in rows if g['target']]
    return dict(objective=float(np.mean(populated)) if populated else 0.,groups=rows)


def prepare_data(folder,output):
    root=project_root();folder,output=Path(folder).resolve(),Path(output).resolve()
    assert not (output/'data-split.json').exists()
    rows=json.loads((folder/'development-split.json').read_text())['windows']
    manifest=json.loads((folder/'manifest.json').read_text())
    new=[r for r in manifest['runs'] if r['role']=='train_rich']
    assert len(new)==12 and all(r['status']=='completed' for r in new)
    for r in new:
        stem=Path(r['root']).stem;base=output/'inputs'
        rows.append(dict(r,split='train',provenance_experiment=folder.name,
            **{k:str((base/sub/f'{stem}.{ext}').relative_to(root)) for k,sub,ext in
               [('waveform','waveforms','npz'),('truth','truth','npz'),('metadata','metadata','json')]}))
    assert len({r['root'] for r in rows})==len(rows)
    assert len({r['seed'] for r in rows})==len(rows)
    calibration=json.loads((folder/'calibration.json').read_text());cfg=NetworkConfig()
    counts=defaultdict(lambda:dict(windows=0,alpha250=0,gamma250=0,alpha500=0,gamma500=0,collisions=0))
    for r in rows:
        r['input_sha256']={k:digest(root/r[k]) for k in ('waveform','truth','metadata')}
        metadata=json.loads((root/r['metadata']).read_text())
        assert metadata['root_sha256']==r['root_sha256'] and metadata['seed']==r['seed']
        if r['split']=='regression':continue  # Open no regression labels during preparation/training.
        waveform=dict(np.load(root/r['waveform']));truth=dict(np.load(root/r['truth']))
        features=voltage_features(waveform['time_ns'],waveform['analog_mV'],calibration)
        labels,collision=supervised_targets(waveform['time_ns'],truth,cfg)
        destination=output/'features'/f"{Path(r['root']).stem}.npz";destination.parent.mkdir(exist_ok=True)
        np.savez_compressed(destination,features=features,labels=labels)
        r['features']=str(destination.relative_to(root));r['feature_sha256']=digest(destination)
        count=counts[r['split']];count['windows']+=1;count['collisions']+=sum(collision.values())
        for p in ('alpha','gamma'):
            for q in (250,500):
                keep=(truth['particle_name']==p)&np.isfinite(truth['first_active_ns'])&(truth['charge_pe']>=q)
                count[f'{p}{q}']+=int(keep.sum())
    save(output/'data-split.json',dict(created=datetime.now().astimezone().isoformat(),windows=rows,
        counts=dict(counts),calibration_sha256=digest(folder/'calibration.json'),config=asdict(cfg),
        split_policy='Whole native ROOT/seed before crops; regression held outside training/selection; new tests not generated',
        training_physics='Native shared SiPM cells; no linear combination of independently digitized inputs'))
    print('DATA READY',dict(counts),flush=True)


def load_items(rows,root):
    items=[]
    for r in rows:
        assert digest(root/r['features'])==r['feature_sha256']
        data=dict(np.load(root/r['features']));wave=dict(np.load(root/r['waveform']))
        truth=dict(np.load(root/r['truth']));meta=json.loads((root/r['metadata']).read_text())
        events=[]
        for p in ('alpha','gamma'):
            keep=(truth['particle_name']==p)&np.isfinite(truth['first_active_ns'])&(truth['charge_pe']>=250.)
            keep&=(truth['first_active_ns']>=2000)&(truth['first_active_ns']<18000)
            events.append(truth['first_active_ns'][keep])
        items.append(dict(row=r,t=wave['time_ns'],voltage=wave['analog_mV'],truth=truth,
                          gate_ns=meta['gate_ns'],events=events,**data))
    return items


def validation_choice(model,items,cfg):
    dense=[predict_features(model,r['features'],cfg) for r in items]
    choices=[]
    for a in (.15,.25,.4,.55,.7):
        for g in (.15,.25,.4,.55,.7):
            records=[]
            for item,prediction in zip(items,dense):
                fit=decode(item['t'],prediction,(a,g),cfg)
                records.append(dict(row=item['row'],fit=fit,truth=item['truth'],gate_ns=item['gate_ns'],
                    score=score_mixed(fit,item['truth'],item['gate_ns'])))
            choices.append(dict(thresholds=[a,g],**summarize_scores(records)))
    return max(choices,key=lambda r:(r['objective'],sum(r['thresholds']))),choices


def train(folder,output,kind,steps=2400,threads=2,batch_size=8):
    folder,output=Path(folder),Path(output);root=project_root();destination=output/kind
    destination.mkdir(exist_ok=False)
    torch.set_num_threads(threads);torch.manual_seed(261082001 if kind=='cnn' else 261082002)
    rng=np.random.default_rng(261082003);cfg=NetworkConfig()
    data=json.loads((output/'data-split.json').read_text())
    training=load_items([r for r in data['windows'] if r['split'] in ('train','train_single')],root)
    validation=load_items([r for r in data['windows'] if r['split']=='validation'],root)
    assert len(validation)==20
    mixed=[i for i,r in enumerate(training) if r['row']['split']=='train']
    singles=[i for i,r in enumerate(training) if r['row']['split']=='train_single']
    focused={(which,k):[i for i in pool if len(training[i]['events'][k])]
             for which,pool in [('mixed',mixed),('single',singles)] for k in (0,1)}
    model=MixedEventNet(kind,cfg);optimizer=torch.optim.AdamW(model.parameters(),lr=.001,weight_decay=.0001)
    specification=dict(kind=kind,config=asdict(cfg),steps=steps,batch_size=batch_size,threads=threads,
        optimizer='AdamW lr=0.001, weight_decay=0.0001; clipped gradient norm 5',
        training_seed=261082001 if kind=='cnn' else 261082002,crop_seed=261082003,
        sampling='85% native mixed / 15% isolated-history; 60% class-balanced event-centred crops, 40% uniform crops',
        loss='Per-class focal peak loss normalized by class positives + 2 masked logQ Huber + 0.5 masked timing Huber',
        selection='Class and gamma-spectrum/alpha-rate macro full recovery F1 on 20 entire disclosed development ROOTs',
        thresholds=[.15,.25,.4,.55,.7],inference_truth_inputs=[],source_rate_input=False,
        parameters=sum(p.numel() for p in model.parameters()),
        source_sha256={n:digest(Path(__file__).with_name(n)) for n in ('alpha_gamma_ml.py','alpha_gamma_ml_study.py')},
        environment=dict(python=platform.python_version(),torch=torch.__version__,numpy=np.__version__,
                         scipy=importlib.metadata.version('scipy'),cuda_available=torch.cuda.is_available()))
    save(destination/'training-config.json',specification)
    history=[];best=None;start=time.monotonic()
    print('TRAIN',kind,'windows',len(training),'validation',len(validation),'parameters',specification['parameters'],flush=True)
    for step in range(1,steps+1):
        xx=[];yy=[]
        for _ in range(batch_size):
            group='mixed' if rng.random()<.85 else 'single';pool=mixed if group=='mixed' else singles
            focus=rng.random()<.6;k=int(rng.integers(2))
            if focus and focused[(group,k)]:pool=focused[(group,k)]
            item=training[int(rng.choice(pool))];n=item['features'].shape[-1]
            if focus and len(item['events'][k]):
                event=float(rng.choice(item['events'][k]));position=int(rng.integers(50,cfg.core_samples-50))
                origin=int(round(event/2.))-cfg.guard_samples-position
            else:
                # The supervised core lies inside [5,15) us for uniform samples.
                origin=int(rng.integers(2500-cfg.guard_samples,7500-cfg.guard_samples-cfg.core_samples+1))
            origin=int(np.clip(origin,0,n-cfg.crop_samples))
            xx.append(item['features'][:,origin:origin+cfg.crop_samples])
            yy.append(item['labels'][:,:,origin:origin+cfg.crop_samples])
        model.train();optimizer.zero_grad(set_to_none=True)
        loss,terms=loss_function(model(torch.from_numpy(np.stack(xx))),torch.from_numpy(np.stack(yy)),cfg)
        if not torch.isfinite(loss):raise FloatingPointError('Nonfinite training objective')
        loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step()
        if step==1 or step%400==0 or step==steps:
            model.eval();choice,grid=validation_choice(model,validation,cfg)
            row=dict(step=step,loss=float(loss.detach()),loss_terms=terms.detach().tolist(),
                objective=choice['objective'],thresholds=choice['thresholds'],elapsed_s=time.monotonic()-start)
            history.append(row);save(destination/'training-history.json',history)
            save(destination/f'validation-{step}.json',grid)
            torch.save(model.state_dict(),destination/f'checkpoint-{step}.pt')
            if best is None or choice['objective']>best['objective']:
                best=dict(choice,step=step);torch.save(model.state_dict(),destination/'model.pt')
                save(destination/'validation-selection.json',best)
            print('STEP',kind,step,'loss',round(row['loss'],4),'val',round(row['objective'],4),
                  'thresholds',row['thresholds'],'seconds',round(row['elapsed_s'],1),flush=True)
    frozen=dict(created=datetime.now().astimezone().isoformat(),kind=kind,state='frozen_before_new_test_generation',
        model_sha256=digest(destination/'model.pt'),config_sha256=digest(destination/'training-config.json'),
        data_split_sha256=digest(output/'data-split.json'),calibration_sha256=digest(folder/'calibration.json'),
        selected=best,training_seconds=time.monotonic()-start,source_sha256=specification['source_sha256'])
    save(destination/'frozen-model.json',frozen);print('FROZEN',kind,best['step'],best['objective'],flush=True)


def load_model(path):
    path=Path(path);frozen=json.loads((path/'frozen-model.json').read_text())
    assert digest(path/'model.pt')==frozen['model_sha256']
    specification=json.loads((path/'training-config.json').read_text());cfg=NetworkConfig(**specification['config'])
    model=MixedEventNet(frozen['kind'],cfg)
    model.load_state_dict(torch.load(path/'model.pt',map_location='cpu',weights_only=True));model.eval()
    return model,cfg,frozen


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=('prepare','train'))
    p.add_argument('--folder',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--kind',choices=('cnn','unet'));p.add_argument('--steps',type=int,default=2400)
    p.add_argument('--threads',type=int,default=2);p.add_argument('--batch-size',type=int,default=8)
    a=p.parse_args()
    if a.stage=='prepare':prepare_data(a.folder,a.output)
    else:train(a.folder,a.output,a.kind,a.steps,a.threads,a.batch_size)


if __name__=='__main__':main()
