"""Blind evaluation, software controls and score verification for component ML."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
import json
from pathlib import Path
import time

import numpy as np
import torch
from scipy.signal import lfilter

from .alpha_component_ml import (PROJECT,EXPERIMENT,OUTPUT,VARIANTS,SEEDS,
                                 load_model,blind_predict)
from .alpha_component_study import load
from .alpha_correlated_fit import CorrelatedConfig,ResidualCalibration,fit_correlated
from .alpha_correlated_gru import GRUEvent
from .alpha_full_signal_study import make_template,summarize
from .alpha_ml_models import make_model,windowed_prediction
from .alpha_pileup_ml import EventNet,waveform_features,dense_prediction,decode_prediction
from .alpha_pileup_fit import FitConfig,fit_waveform,score
from .alpha_rate_campaign import save,digest


def configuration():
    settings=load(EXPERIMENT/'assets/baseline-settings.json')
    selected=load(EXPERIMENT/'assets/physical-selection.json')
    covariance=ResidualCalibration(**load(EXPERIMENT/'assets/covariance.json')['parameters'])
    return settings,CorrelatedConfig(**selected['config']),covariance


def networks():
    torch.set_num_threads(1)
    result={}
    matched=load(OUTPUT/'matched-extra-thresholds.json')
    for variant in VARIANTS:
        for seed in SEEDS:
            name=f'{variant}_{seed}';model,frozen=load_model(name)
            result[name]=(model,frozen['threshold'],matched[name]['threshold'])
    source=EXPERIMENT/'assets/reference-models'
    for kind in ('cnn','unet','transformer','gru'):
        folder=source/kind;frozen=load(folder/'frozen-model.json')
        assert digest(folder/'model.pt')==frozen['weights_sha256']
        model=EventNet() if kind=='cnn' else GRUEvent() if kind=='gru' else make_model(kind)
        model.load_state_dict(torch.load(folder/'model.pt',map_location='cpu',weights_only=True))
        result['old_'+kind]=(model.eval(),frozen['threshold'],None)
    return result


def infer(t,y,nets,*,physical=True):
    """Complete blind boundary; all references to truth/rate are outside it."""
    settings,correlated,covariance=configuration()
    template=make_template(settings);base=FitConfig(**settings['fit_config'])
    tick=time.monotonic();features=waveform_features(t,y,settings);feature_seconds=time.monotonic()-tick
    fits,elapsed,denses={}, {}, {}
    for name,(model,threshold,matched_threshold) in nets.items():
        tick=time.monotonic()
        dense=dense_prediction(model,features) if name=='old_cnn' else windowed_prediction(model,features)
        fits[name]=decode_prediction(t,dense,threshold)
        elapsed[name]=time.monotonic()-tick+feature_seconds
        denses[name]=dense
        if matched_threshold is not None:
            fits[name+'_matched_extra']=decode_prediction(t,dense,matched_threshold)
            elapsed[name+'_matched_extra']=elapsed[name]
    if physical:
        tick=time.monotonic();initial=fit_waveform(t,y,template,base);pilot_seconds=time.monotonic()-tick
        fits['sse']=initial;elapsed['sse']=pilot_seconds
        tick=time.monotonic()
        fits['wls']=fit_correlated(t,y,template,base,covariance,correlated,initial=initial)
        elapsed['wls']=pilot_seconds+time.monotonic()-tick
    return fits,elapsed,denses


def check_freeze():
    frozen=load(OUTPUT/'freeze.json')
    for name,sha in frozen['frozen_files'].items():
        assert digest(PROJECT/name)==sha,name
    return frozen


def evaluate(workers=2):
    frozen=check_freeze();nets=networks()
    manifest=load(EXPERIMENT/'manifest.json')
    assert len(manifest['runs'])==len(manifest['planned_runs'])==load(OUTPUT/'protocol.json')['fresh_test_root_count']
    assert all(r['status']=='completed' for r in manifest['runs'])
    for name in ('windows','models','dense'):(OUTPUT/'test'/name).mkdir(parents=True,exist_ok=True)
    known={r['root_sha256'] for r in load(EXPERIMENT/'data-split.json')['windows']}
    def one(row):
        stem=Path(row['root']).stem;destination=OUTPUT/'test/windows'/(stem+'.json')
        if destination.exists():return load(destination)
        assert row['root_sha256'] not in known
        wp=OUTPUT/'inputs/waveforms'/(stem+'.npz');tp=OUTPUT/'inputs/truth'/(stem+'.npz')
        with np.load(wp) as f:t,y=f['time_ns'].copy(),f['analog_mV'].copy()
        fits,seconds,denses=infer(t,y,nets)
        # Open scoring truth only after all blind methods have completed.
        with np.load(tp) as f:truth={k:f[k] for k in ('first_active_ns','emission_ns','charge_pe')}
        methods={}
        for name,fit in fits.items():
            if 'model_mV' in fit:
                np.savez_compressed(OUTPUT/'test/models'/f'{stem}_{name}.npz',model_mV=fit.pop('model_mV'))
            methods[name]=dict(fit=fit,score=score(fit,truth,[5000.,15000.]),seconds=seconds[name])
        for name,dense in denses.items():np.savez_compressed(OUTPUT/'test/dense'/f'{stem}_{name}.npz',**dense)
        result=dict(root=row['root'],root_sha256=row['root_sha256'],seed=row['seed'],rate_hz=row['rate_hz'],
            waveform=str(wp.relative_to(PROJECT)),truth=str(tp.relative_to(PROJECT)),gate_ns=[5000.,15000.],
            waveform_sha256=digest(wp),truth_sha256=digest(tp),methods=methods,
            inference_truth_inputs=[],rate_input_at_inference=False)
        save(destination,result)
        print('TEST_DONE',row['rate_hz'],row['seed'],[(n,m['score']['recovered_count']) for n,m in methods.items() if not n.endswith('matched_extra')],flush=True)
        return result
    with ThreadPoolExecutor(max_workers=workers) as pool:
        rows=[f.result() for f in as_completed([pool.submit(one,r) for r in manifest['runs']])]
    rows.sort(key=lambda r:(r['rate_hz'],r['seed']))
    names=list(rows[0]['methods'])
    rates={name:summarize([dict(rate_hz=r['rate_hz'],fit=r['methods'][name]['fit'],score=r['methods'][name]['score']) for r in rows]) for name in names}
    result=dict(state='completed',windows=rows,rates=rates,frozen_at=frozen['frozen_at'],method_order=names)
    save(OUTPUT/'test-summary.json',result)
    return result


def controls():
    check_freeze();nets=networks();settings,_,_=configuration()
    template=make_template(settings);base=FitConfig(**settings['fit_config'])
    t=np.arange(10000)*2.+.25;rng=np.random.default_rng(262043000)
    cases=[('empty',[],[])]
    cases += [(f'single_{q:g}',[8000.3],[q]) for q in (500.,1000.,4000.,8000.,16000.,40000.)]
    for gap in (4.,8.,10.,12.,20.,40.,80.,200.):
        for ratio in (.125,1.,4.):
            cases.append((f'double_{gap:g}ns_{ratio:g}',[8000.3,8000.3+gap],[4000.,4000.*ratio]))
    rows=[]
    for noise in ('white','correlated'):
        for label,times,charges in cases:
            times,charges=np.asarray(times),np.asarray(charges)
            signal=template(t[:,None]-times[None,:])@(charges*base.area_per_pe_mV_ns*template.unit_area_peak)
            eps=rng.normal(0,.1,len(t))
            if noise=='correlated':eps=lfilter([np.sqrt(1-.98**2)],[1.,-.98],eps)
            y=base.baseline_mV+signal+eps
            fits,seconds,_=infer(t,y,nets,physical=True)
            truth=dict(first_active_ns=times,emission_ns=times,charge_pe=charges)
            rows.append(dict(label=label,noise=noise,truth={k:v.tolist() for k,v in truth.items()},
                methods={name:dict(fit=fit,score=score(fit,truth,[5000.,15000.]),seconds=seconds[name]) for name,fit in fits.items()}))
        print('CONTROLS',noise,len(cases),flush=True)
    model_folder=OUTPUT/'control-models';model_folder.mkdir(exist_ok=True)
    for row in rows:
        for name,result in row['methods'].items():
            if 'model_mV' in result['fit']:
                path=model_folder/f'{row["noise"]}_{row["label"]}_{name}.npz'
                np.savez_compressed(path,model_mV=result['fit'].pop('model_mV'))
                result['saved_model']=str(path.relative_to(PROJECT))
    save(OUTPUT/'controls.json',dict(cases=rows,seed=262043000,
        scope='Matched mean template with independent synthetic white/AR1 noise; not native shared-cell separation limit'))


def diagnostics():
    summary=load(OUTPUT/'test-summary.json')
    names=summary['method_order'];table=[];null=[];neighbors=[]
    for name in names:
        scores=[r['methods'][name]['score'] for r in summary['windows']]
        keys=('target_count','recovered_count','missed_count','matched_count','unmatched_predictions','merge_flags','split_flags','reconstructed_count')
        row={key:sum(s[key] for s in scores) for key in keys}
        row.update(method=name,matched_charge_failed=row['matched_count']-row['recovered_count'],
                   exposure_us=sum(s['duration_ns'] for s in scores)/1000.)
        dq=np.array([v for s in scores for v in s['relative_charge_errors']])
        row['charge_p16_p50_p84']=np.quantile(dq,[.16,.5,.84]).tolist() if len(dq) else None
        table.append(row)
    for row in summary['windows']:
        with np.load(PROJECT/row['truth']) as f:truth={k:f[k] for k in ('first_active_ns','emission_ns','charge_pe')}
        valid=np.flatnonzero(np.isfinite(truth['first_active_ns'])&(truth['charge_pe']>=500.))
        valid=valid[np.argsort(truth['first_active_ns'][valid])]
        gap=np.full(len(valid),np.inf)
        if len(valid)>1:
            dd=np.diff(truth['first_active_ns'][valid]);gap[1:]=np.minimum(gap[1:],dd);gap[:-1]=np.minimum(gap[:-1],dd)
        by_id={int(i):float(g) for i,g in zip(valid,gap)}
        for name in names:
            fit=row['methods'][name]['fit'];times=np.asarray(fit['times_ns']);q=np.asarray(fit['charge_pe'])
            s=row['methods'][name]['score'];matched=dict(s['matched_pairs']);groups=np.zeros((8,3),dtype=int)
            for i in s['target_indices']:
                b=int(np.searchsorted([4.,8.,12.,20.,40.,80.,200.],by_id[i],side='right'))
                outcome=0 if i not in matched else (2 if abs(q[matched[i]]/truth['charge_pe'][i]-1)<=.2 else 1)
                groups[b,outcome]+=1
            neighbors.append(dict(rate_hz=row['rate_hz'],seed=row['seed'],method=name,groups=groups.tolist(),
                definition='Nearest other Q>=500 finite-time alpha, including guard events; post-inference grouping only'))
            if row['rate_hz']<40000000:continue
            gate=(times>=5000.)&(times<15000.)
            # Identical to prior study: permute reported central output times,
            # preserving output counts and Q distribution; diagnostic only.
            original=dict(times_ns=times[gate].tolist(),charge_pe=q[gate].tolist())
            actual=score(original,truth,[5000.,15000.])['recovered_count']
            shifted=[]
            for shift in (137.,353.,911.,1703.,2999.):
                pred=dict(times_ns=(5000.+(times[gate]-5000.+shift)%10000.).tolist(),charge_pe=q[gate].tolist())
                shifted.append(score(pred,truth,[5000.,15000.])['recovered_count'])
            null.append(dict(rate_hz=row['rate_hz'],seed=row['seed'],method=name,actual=actual,shifted=shifted,
                             target_count=row['methods'][name]['score']['target_count']))
    save(OUTPUT/'report-tables.json',dict(totals=table,rates=summary['rates'],timing_null=null,neighbor_groups=neighbors))
    print('DIAGNOSTICS_DONE',len(table),len(null),flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('evaluate','controls','diagnostics'))
    p.add_argument('--workers',type=int,default=2);a=p.parse_args()
    if a.stage=='evaluate':evaluate(a.workers)
    elif a.stage=='controls':controls()
    else:diagnostics()


if __name__=='__main__':main()
