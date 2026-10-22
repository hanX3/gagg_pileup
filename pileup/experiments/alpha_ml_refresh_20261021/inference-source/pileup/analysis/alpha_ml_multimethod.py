"""Train and freeze U-Net/Transformer/hybrid, then compare on common ROOTs.

Run from project root. Only validation ROOTs select model/threshold/prior;
previous 17 tests are a public regression benchmark. New tests are not opened
for model selection. Truth is opened after blind inference when evaluating.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from datetime import datetime
from dataclasses import replace,asdict
import json
from pathlib import Path
import platform
import time
import numpy as np
import torch

from .alpha_rate_campaign import digest,save
from .alpha_full_signal_study import make_template,summarize
from .alpha_pileup_fit import FitConfig,fit_waveform,score
from .alpha_pileup_ml import EventNet,waveform_features,dense_prediction,decode_prediction,load_pilot,loss_function
from .alpha_ml_models import make_model,windowed_prediction
from .alpha_pileup_hybrid import HybridConfig,hybrid_waveform

METHODS=('fit','cnn','unet','transformer','cnn_physical')
OLD=Path('jupyter/artifacts/alpha_ml_comparison_20261014_162000')


def source_hashes():
    folder=Path(__file__).parent
    return {name:digest(folder/name) for name in ('alpha_ml_multimethod.py','alpha_ml_models.py',
            'alpha_pileup_hybrid.py','alpha_pileup_ml.py','alpha_pileup_fit.py','alpha_full_signal_study.py')}


def choice(model,validation):
    dense=[windowed_prediction(model,item['features']) for item in validation]
    options=[]
    for threshold in (.15,.25,.35,.45,.55,.65,.75,.85):
        groups=defaultdict(lambda:np.zeros(3));rows=[]
        for item,prediction in zip(validation,dense):
            result=decode_prediction(item['t'],prediction,threshold)
            s=score(result,item['truth'],item['gate_ns'])
            groups[item['record']['rate_hz']]+=[s['target_count'],s['reconstructed_count'],s['recovered_count']]
            rows.append(dict(label=item['record']['label'],rate_hz=item['record']['rate_hz'],fit=result,score=s))
        objective=float(np.mean([2*v[2]/max(v[0]+v[1],1.) for v in groups.values()]))
        options.append(dict(threshold=threshold,objective=objective,windows=rows))
    return max(options,key=lambda r:(r['objective'],r['threshold']))


def train(kind,experiment,output,steps=2400,threads=2):
    if not 1<=threads<=4 or steps<1:raise ValueError('Invalid training budget')
    destination=output/'models'/kind;destination.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(threads);torch.manual_seed(2610141620)
    rng=np.random.default_rng(2610141620)
    settings=json.loads((experiment/'calibration-settings.json').read_text())
    split=json.loads((experiment/'data-split.json').read_text())
    items=[load_pilot(r,settings) for r in split['windows']]
    training=[r for r in items if r['record']['split']=='train']
    validation=[r for r in items if r['record']['split']=='validation']
    assert len(training)==26 and len(validation)==10
    assert not {r['record']['root_sha256'] for r in training}&{r['record']['root_sha256'] for r in validation}
    by_rate=defaultdict(list)
    for item in training:by_rate[item['record']['rate_hz']].append(item)
    rates=sorted(by_rate);model=make_model(kind)
    config=dict(architecture=kind,parameters=sum(p.numel() for p in model.parameters()),
        steps=steps,batch_size=16,seed=2610141620,optimizer='AdamW',learning_rate=.001,weight_decay=.0001,
        gradient_clip_norm=5.,cpu_threads=threads,features='Same fixed five voltage-derived features as original CNN',
        crop_samples=1024,core_samples=500,left_guard_samples=262,
        inference='Sliding 1024-point windows, emit each central 500-point core once; no absolute time/rate/truth inputs',
        loss='Identical weighted soft BCE + 2x masked log-charge Huber + 0.5x time offset Huber',
        split_sha256=digest(experiment/'data-split.json'),calibration_sha256=digest(experiment/'calibration-settings.json'),
        model_selection='Validation-only macro charge-qualified F1; global threshold grid .15..85',
        environment=dict(python=platform.python_version(),torch=torch.__version__,numpy=np.__version__),
        sources=source_hashes(),paper_reproduction=False)
    save(destination/'config.json',config)
    optimizer=torch.optim.AdamW(model.parameters(),lr=.001,weight_decay=.0001)
    history=[];best=None;started=time.monotonic()
    print('TRAIN',kind,'parameters',config['parameters'],flush=True)
    for step in range(1,steps+1):
        x,labels=[],[]
        for _ in range(16):
            group=by_rate[rates[int(rng.integers(len(rates)))]]
            item=group[int(rng.integers(len(group)))];start=int(rng.integers(2250,6739))
            x.append(item['features'][:,start:start+1024]);labels.append(item['targets'][:,start:start+1024])
        model.train();optimizer.zero_grad(set_to_none=True)
        loss,parts=loss_function(model(torch.from_numpy(np.stack(x))),torch.from_numpy(np.stack(labels)))
        if not torch.isfinite(loss):raise FloatingPointError('Nonfinite training loss')
        loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step()
        if step==1 or step%200==0 or step==steps:
            model.eval();selected=choice(model,validation)
            row=dict(step=step,loss=float(loss.detach()),loss_terms=[float(v.detach()) for v in parts],
                     validation_objective=selected['objective'],threshold=selected['threshold'],elapsed_s=time.monotonic()-started)
            history.append(row);save(destination/'history.json',history)
            if best is None or selected['objective']>best['objective']:
                best=dict(selected,step=step);torch.save(model.state_dict(),destination/'model.pt')
                save(destination/'validation-selection.json',best)
            print('STEP',kind,step,'loss',round(row['loss'],4),'val',round(selected['objective'],4),
                  'threshold',selected['threshold'],'seconds',round(row['elapsed_s'],1),flush=True)
    frozen=dict(state='frozen_before_test',created=datetime.now().astimezone().isoformat(),architecture=kind,
        weights_sha256=digest(destination/'model.pt'),config_sha256=digest(destination/'config.json'),
        threshold=best['threshold'],step=best['step'],validation_objective=best['objective'],
        training_seconds=time.monotonic()-started,parameters=config['parameters'],sources=source_hashes())
    save(destination/'frozen.json',frozen);print('MODEL FROZEN',kind,best['step'],best['threshold'],flush=True)


def load_network(kind,output):
    if kind=='cnn':
        frozen=json.loads((output/'models/cnn/origin-frozen-model.json').read_text())
        model=EventNet()
    else:
        frozen=json.loads((output/f'models/{kind}/frozen.json').read_text());model=make_model(kind)
    assert digest(output/f'models/{kind}/model.pt')==frozen['weights_sha256']
    model.load_state_dict(torch.load(output/f'models/{kind}/model.pt',map_location='cpu',weights_only=True))
    model.eval();return model,frozen


def select_hybrid(experiment,output):
    torch.set_num_threads(1)
    settings=json.loads((experiment/'calibration-settings.json').read_text())
    split=json.loads((experiment/'data-split.json').read_text())
    validation=[r for r in split['windows'] if r['split']=='validation']
    model,frozen=load_network('cnn',output)
    selected=None;trials=[]
    for strength in (0.,.1,1.):
        config=HybridConfig(threshold=frozen['threshold'],prior_strength=strength)
        rows=[];groups=defaultdict(lambda:np.zeros(3))
        for record in validation:
            with np.load(record['waveforms']) as f:t,y=f['time_ns'],f['analog_mV']
            dense=dense_prediction(model,waveform_features(t,y,settings))
            started=time.monotonic();fit=hybrid_waveform(t,y,dense,settings,config);fit.pop('model_mV')
            with np.load(record['truth']) as f:truth={k:f[k] for k in ('first_active_ns','emission_ns','charge_pe')}
            metadata=json.loads(Path(record['metadata']).read_text())
            s=score(fit,truth,metadata['gate_ns'])
            groups[record['rate_hz']]+=[s['target_count'],s['reconstructed_count'],s['recovered_count']]
            rows.append(dict(label=record['label'],rate_hz=record['rate_hz'],fit=fit,score=s,seconds=time.monotonic()-started))
            print('HYBRID VALIDATION',strength,record['label'],s['recovered_count'],'/',s['target_count'],flush=True)
        objective=float(np.mean([2*v[2]/max(v[0]+v[1],1.) for v in groups.values()]))
        trial=dict(config=asdict(config),validation_objective=objective,windows=rows);trials.append(trial)
        save(output/'hybrid-validation-trials.json',trials)
        if selected is None or objective>selected['validation_objective']:selected=trial
    save(output/'hybrid-config.json',dict(selected['config'],selection='Validation-only macro charge-qualified F1',
         validation_objective=selected['validation_objective'],cnn_weights_sha256=frozen['weights_sha256']))
    methods={k:json.loads((output/f'models/{k}/frozen.json').read_text()) for k in ('unet','transformer')}
    methods['cnn']=frozen
    save(output/'frozen-methods.json',dict(state='all_methods_frozen_before_test',created=datetime.now().astimezone().isoformat(),
        methods=methods,hybrid_config_sha256=digest(output/'hybrid-config.json'),sources=source_hashes(),
        calibration_sha256=digest(experiment/'calibration-settings.json'),split_sha256=digest(experiment/'data-split.json'),
        test_tuning=False,previous_17_tests='Regression only'))
    print('ALL METHODS FROZEN','hybrid strength',selected['config']['prior_strength'],flush=True)


def evaluate(experiment,output,data,benchmark):
    torch.set_num_threads(1)
    frozen=json.loads((output/'frozen-methods.json').read_text())
    assert source_hashes()==frozen['sources']
    assert digest(output/'hybrid-config.json')==frozen['hybrid_config_sha256']
    settings=json.loads((experiment/'calibration-settings.json').read_text())
    assert digest(experiment/'calibration-settings.json')==frozen['calibration_sha256']
    networks={k:load_network(k,output) for k in ('cnn','unet','transformer')}
    raw_config=json.loads((output/'hybrid-config.json').read_text())
    hybrid_config=HybridConfig(**{k:raw_config[k] for k in asdict(HybridConfig())})
    if benchmark=='regression':
        previous=json.loads((OLD/'summary.json').read_text());runs=previous['windows']
        cache={r['root_sha256']:r for r in runs}
    else:
        runs=json.loads((experiment/'manifest.json').read_text())['runs'];cache={}
        assert len(runs)==22 and all(r['status']=='completed' for r in runs)
    destination=output/benchmark;destination.mkdir(exist_ok=True)
    (destination/'windows').mkdir(exist_ok=True)
    roots=[];training={r['root_sha256'] for r in json.loads((experiment/'data-split.json').read_text())['windows']}
    for run in runs:
        stem=Path(run['root']).stem;result_path=destination/'windows'/f'{stem}.json'
        assert run['root_sha256'] not in training
        if result_path.exists():
            row=json.loads(result_path.read_text());assert row['root_sha256']==run['root_sha256'];roots.append(row);continue
        with np.load(data/'waveforms'/f'{stem}.npz') as f:t,y=f['time_ns'],f['analog_mV']
        started=time.monotonic();features=waveform_features(t,y,settings);feature_seconds=time.monotonic()-started
        entries={};cnn_dense=None
        for kind in ('cnn','unet','transformer'):
            model,checkpoint=networks[kind];started=time.monotonic()
            dense=dense_prediction(model,features) if kind=='cnn' else windowed_prediction(model,features)
            fit=decode_prediction(t,dense,checkpoint['threshold'])
            entries[kind]=dict(fit=fit,seconds=time.monotonic()-started+feature_seconds)
            np.savez_compressed(destination/'windows'/f'{stem}_{kind}.npz',**dense)
            if kind=='cnn':cnn_dense=dense
        started=time.monotonic();fit=hybrid_waveform(t,y,cnn_dense,settings,hybrid_config)
        np.savez_compressed(destination/'windows'/f'{stem}_cnn_physical.npz',model_mV=fit.pop('model_mV'))
        entries['cnn_physical']=dict(fit=fit,seconds=time.monotonic()-started+entries['cnn']['seconds'])
        if run['root_sha256'] in cache:
            saved=cache[run['root_sha256']]
            fit=saved['fit'];seconds=saved['fit_seconds']
            assert entries['cnn']['fit']==saved['neural_fit']
            entries['fit']=dict(fit=fit,seconds=seconds,cached_from_first_round=True)
        else:
            started=time.monotonic();fit=fit_waveform(t,y,make_template(settings),FitConfig(**settings['fit_config']))
            np.savez_compressed(destination/'windows'/f'{stem}_fit.npz',model_mV=fit.pop('model_mV'))
            entries['fit']=dict(fit=fit,seconds=time.monotonic()-started,cached_from_first_round=False)
        # Test ancestry and source rate/seed metadata are not inputs to any method.
        with np.load(data/'truth'/f'{stem}.npz') as f:truth={k:f[k] for k in ('first_active_ns','emission_ns','charge_pe')}
        meta=json.loads((data/'metadata'/f'{stem}.json').read_text())
        for entry in entries.values():entry['score']=score(entry['fit'],truth,meta['gate_ns'])
        row=dict(label=run['label'],root=run['root'],root_sha256=run['root_sha256'],seed=run['seed'],
            rate_hz=run['rate_hz'],gate_ns=meta['gate_ns'],waveform_path=str(data/'waveforms'/f'{stem}.npz'),
            truth_path=str(data/'truth'/f'{stem}.npz'),waveform_sha256=digest(data/'waveforms'/f'{stem}.npz'),
            truth_sha256=digest(data/'truth'/f'{stem}.npz'),methods=entries)
        save(result_path,row);roots.append(row)
        save(destination/'progress.json',dict(completed=len(roots),planned=len(runs)))
        print('COMPARE',benchmark,run['label'],{k:v['score']['recovered_count'] for k,v in entries.items()},
              '/',entries['fit']['score']['target_count'],flush=True)
    rates={k:summarize([dict(rate_hz=r['rate_hz'],**r['methods'][k]) for r in roots]) for k in METHODS}
    save(destination/'summary.json',dict(state='completed',benchmark=benchmark,windows=roots,rates=rates,
        method_order=list(METHODS),frozen_methods_sha256=digest(output/'frozen-methods.json'),
        adequate_dynamic_range=True,inference_truth_fields=[],energy_calibration=False))
    print('COMPARISON COMPLETE',benchmark,len(roots),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=('train','select-hybrid','evaluate'))
    p.add_argument('--experiment',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--model',choices=('unet','transformer'));p.add_argument('--steps',type=int,default=2400)
    p.add_argument('--threads',type=int,default=2);p.add_argument('--data',type=Path)
    p.add_argument('--benchmark',choices=('regression','new_test'),default='new_test')
    args=p.parse_args()
    if args.stage=='train':
        if args.model is None:p.error('--model required')
        train(args.model,args.experiment,args.output,args.steps,args.threads)
    elif args.stage=='select-hybrid':select_hybrid(args.experiment,args.output)
    else:
        if args.data is None:p.error('--data required')
        evaluate(args.experiment,args.output,args.data,args.benchmark)


if __name__=='__main__':main()
