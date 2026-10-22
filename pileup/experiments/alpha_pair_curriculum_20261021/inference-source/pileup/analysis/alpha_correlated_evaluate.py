"""Frozen common-window evaluation and independent software/null controls."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
from pathlib import Path
import time

import numpy as np

from .alpha_correlated_fit import CorrelatedConfig,ResidualCalibration,InnovationsMetric,optimize,fit_correlated
from .alpha_correlated_study import PROJECT,EXPERIMENT,OUTPUT,OLD,settings,pilot,fit_row
from .alpha_full_signal_study import make_template,summarize
from .alpha_local_fit import LocalConfig,fit_local
from .alpha_pileup_fit import FitConfig,score
from .alpha_rate_campaign import digest,save


def frozen():
    result=json.loads((OUTPUT/'freeze.json').read_text())
    for path,expected in result['frozen_files'].items():assert digest(PROJECT/path)==expected,path
    return result


def physical_window(row,selection):
    folder=OUTPUT/'test/windows';folder.mkdir(parents=True,exist_ok=True);stem=Path(row['root']).stem;path=folder/(stem+'.json')
    if path.exists():return json.loads(path.read_text())
    t,y,truth_path,gate,shift,waveform,initial,baseline_seconds=pilot(row)
    cfg=LocalConfig(**json.loads((OLD/'freeze.json').read_text())['selected']['config'])
    d=settings();template=make_template(d);base=FitConfig(**d['fit_config'])
    tick=time.monotonic();local=fit_local(t,y,template,base,cfg,initial=initial);local_seconds=time.monotonic()-tick
    with np.load(truth_path) as z:truth={k:z[k].copy() for k in ('first_active_ns','emission_ns','charge_pe')}
    methods={}
    for name,fit,seconds in [('baseline',initial,baseline_seconds),('local',local,baseline_seconds+local_seconds)]:
        np.savez_compressed(folder/(stem+'_'+name+'.npz'),model_mV=fit['model_mV'])
        methods[name]=dict(fit={k:v for k,v in fit.items() if k!='model_mV'},score=score(fit,truth,gate),seconds=seconds)
    for name,selected in selection.items():
        record=fit_row(row,CorrelatedConfig(**selected['config']),phase='test',tag=name)
        methods[name]={k:record[k] for k in ('fit','score','seconds')}
    result=dict(root=row['root'],root_sha256=row['root_sha256'],seed=row['seed'],rate_hz=row['rate_hz'],
        waveform=str(waveform.relative_to(PROJECT)),truth=str(truth_path.relative_to(PROJECT)),gate_ns=gate,
        waveform_sha256=digest(waveform),truth_sha256=digest(truth_path),methods=methods,inference_truth_inputs=[])
    save(path,result);print('TEST_PHYSICAL',stem,[(k,v['score']['recovered_count']) for k,v in methods.items()],flush=True)
    return result


def physical(workers=2):
    f=frozen();manifest=json.loads((EXPERIMENT/'manifest.json').read_text())
    assert len(manifest['runs'])==23 and all(r['status']=='completed' and 'export_metadata' in r for r in manifest['runs'])
    with ThreadPoolExecutor(max_workers=workers) as pool:
        tasks=[pool.submit(physical_window,r,f['selected']) for r in manifest['runs']]
        rows=[p.result() for p in as_completed(tasks)]
    summarize_test(rows)


def neural():
    f=frozen()
    import torch
    from .alpha_correlated_gru import load
    from .alpha_pileup_ml import EventNet,waveform_features,dense_prediction,decode_prediction
    from .alpha_ml_models import windowed_prediction
    torch.set_num_threads(1);model,gru,gru_settings=load()
    old=EXPERIMENT/'reference-cnn'
    cnn_frozen=json.loads((old/'frozen-model.json').read_text());cnn_settings=json.loads((old/'analysis_config.json').read_text())
    assert digest(old/'model.pt')==cnn_frozen['weights_sha256']
    cnn=EventNet();cnn.load_state_dict(torch.load(old/'model.pt',map_location='cpu',weights_only=True));cnn.eval()
    save(OUTPUT/'cnn-reference.json',dict(directory=str(old.relative_to(PROJECT)),frozen=cnn_frozen,
        config_sha256=digest(old/'analysis_config.json'),retrained=False))
    rows=[]
    for path in sorted((OUTPUT/'test/windows').glob('*.json')):
        row=json.loads(path.read_text())
        if all(k in row['methods'] for k in ('cnn','gru')):rows.append(row);continue
        with np.load(PROJECT/row['waveform']) as z:t,y=z['time_ns'].copy(),z['analog_mV'].copy()
        predictions={}
        for name,network,configuration,threshold in [('cnn',cnn,cnn_settings,cnn_frozen['threshold']),('gru',model,gru_settings,gru['threshold'])]:
            tick=time.monotonic();features=waveform_features(t,y,configuration)
            dense=dense_prediction(network,features) if name=='cnn' else windowed_prediction(network,features)
            predictions[name]=(decode_prediction(t,dense,threshold),time.monotonic()-tick)
        # Neural predictions are complete before opening scoring truth.
        with np.load(PROJECT/row['truth']) as z:truth={k:z[k].copy() for k in ('first_active_ns','emission_ns','charge_pe')}
        for name,(fit,elapsed) in predictions.items():row['methods'][name]=dict(fit=fit,score=score(fit,truth,row['gate_ns']),seconds=elapsed)
        save(path,row);rows.append(row)
        print('TEST_NEURAL',Path(row['root']).stem,[(k,row['methods'][k]['score']['recovered_count']) for k in ('cnn','gru')],flush=True)
    assert len(rows)==23
    summarize_test(rows)


def summarize_test(rows):
    methods=sorted(set.intersection(*(set(r['methods']) for r in rows)))
    rates={name:summarize([dict(rate_hz=r['rate_hz'],fit=r['methods'][name]['fit'],score=r['methods'][name]['score']) for r in rows]) for name in methods}
    save(OUTPUT/'test-summary.json',dict(windows=sorted(rows,key=lambda r:(r['rate_hz'],r['seed'])),rates=rates,
        test_windows=len(rows),frozen_at=json.loads((OUTPUT/'freeze.json').read_text())['frozen_at'],
        adequate_dynamic_range=True,inference_truth_inputs=[],energy_calibration=False,
        cohort='New native seeds generated after common freeze; high rates have only two windows'))


def controls():
    f=frozen();cal=ResidualCalibration(**json.loads((OUTPUT/'calibration/covariance.json').read_text())['parameters'])
    d=settings();template=make_template(d);base=FitConfig(**d['fit_config'])
    old_config=LocalConfig(**json.loads((OLD/'freeze.json').read_text())['selected']['config'])
    rng=np.random.default_rng(262026700);t=np.arange(0.,6000.,2.)+.137
    plans=[('empty',[],[],0.,0.)]*12
    plans += [('single',[1000.3],[q],0.,0.) for q in (500.,1000.,4000.,16000.) for _ in range(3)]
    plans += [('double',[1000.3,1000.3+gap],[4000.,4000.*ratio],gap,ratio)
              for gap in (4.,8.,12.,20.,40.,80.,200.) for ratio in (.125,.25,.5,1.,2.,4.) for _ in range(2)]
    from .alpha_pileup_fit import fit_waveform
    rows=[]
    for noise in ('white','correlated'):
        for i,(kind,onsets,charges,gap,ratio) in enumerate(plans):
            identifier=f'{noise}_{i:03}';path=OUTPUT/'controls'/(identifier+'.json')
            a=np.asarray(charges)*base.area_per_pe_mV_ns*template.unit_area_peak
            mu=template(t[:,None]-np.asarray(onsets)[None,:])@a
            white=rng.normal(0,.1,len(t))
            if noise=='correlated':
                latent=np.empty(len(t));latent[0]=rng.normal()
                draws=rng.normal(size=len(t)-1)
                for j in range(1,len(t)):latent[j]=cal.rho*latent[j-1]+np.sqrt(1-cal.rho**2)*draws[j-1]
                white+=np.sqrt(cal.dark_variance_mV2+cal.shot_coefficient_mV*np.maximum(mu,0))*latent
            y=10+mu+white
            # Generate the entire deterministic stream even when resuming.
            if path.exists():rows.append(json.loads(path.read_text()));continue
            initial=fit_waveform(t,y,template,base)
            fits={'baseline':initial,'local':fit_local(t,y,template,base,old_config,initial=initial)}
            for method,chosen in f['selected'].items():fits[method]=fit_correlated(t,y,template,base,cal,CorrelatedConfig(**chosen['config']),initial=initial)
            truth=dict(first_active_ns=np.asarray(onsets),emission_ns=np.asarray(onsets),charge_pe=np.asarray(charges))
            row=dict(id=identifier,noise=noise,kind=kind,gap_ns=gap,ratio=ratio,
                scores={k:score(v,truth,[500.,4000.]) for k,v in fits.items()},
                global_evaluations=sum(r.get('evaluations',0) for r in fits['global'].get('optimizer_trials',[]) if r['kind']=='differential_evolution'),
                known_count_input=False,known_times_input=False)
            np.savez_compressed(OUTPUT/'controls'/(identifier+'.npz'),time_ns=t,analog_mV=y,onsets_ns=onsets,charge_pe=charges)
            save(path,row);rows.append(row)
        print('CONTROLS_DONE',noise,flush=True)
    save(OUTPUT/'controls-summary.json',dict(seed=262026700,rows=rows,
        limitation='Matched mean-template and synthetic noise controls; no native photon statistics/shared-cell response, not detector resolution'))


def shifts():
    result=json.loads((OUTPUT/'test-summary.json').read_text());rows=[]
    shifts_ns=[1379.,2933.,4871.,7039.,8911.]
    for row in result['windows']:
        with np.load(PROJECT/row['truth']) as z:truth={k:z[k].copy() for k in ('first_active_ns','emission_ns','charge_pe')}
        for method,data in row['methods'].items():
            for delta in shifts_ns:
                fit=dict(data['fit']);times=np.asarray(fit['times_ns'],float);central=(times>=5000)&(times<15000)
                shifted=times.copy();shifted[central]=5000+(shifted[central]-5000+delta)%10000
                fit['times_ns']=shifted.tolist()
                rows.append(dict(root=row['root'],rate_hz=row['rate_hz'],method=method,shift_ns=delta,score=score(fit,truth,[5000.,15000.])))
    save(OUTPUT/'time-shift-controls.json',dict(shifts_ns=shifts_ns,rows=rows,
        definition='Circularly shift central-window predictions and preserve their charges/counts; guard predictions stay unchanged. Null control, not exact false-discovery probability.'))


def optimizer_controls():
    """Explicit known-K diagnostics, including three-pulse nonlinear searches.

    These are conditional tests, excluded from blind event-efficiency tables.
    Both searches retain identical fixed trigger-region priors and starting
    grids; DE additionally spends its declared population-search budget.
    """
    f=frozen();d=settings();template=make_template(d);base=FitConfig(**d['fit_config'])
    cal=ResidualCalibration(**json.loads((OUTPUT/'calibration/covariance.json').read_text())['parameters'])
    cfg=CorrelatedConfig(**f['selected']['global']['config'])
    t=np.arange(0.,3000.,2.)+.137;rng=np.random.default_rng(262026701);rows=[]
    from .alpha_pileup_fit import fit_waveform
    for k in (2,3):
        for gap in (4.,8.,12.,20.,40.,80.):
            for ratio in (.125,1.,4.):
                true_t=1000.3+np.arange(k)*gap;true_q=np.array([4000.]+[4000.*ratio]*(k-1))
                a=true_q*base.area_per_pe_mV_ns*template.unit_area_peak
                mu=template(t[:,None]-true_t)@a
                latent=np.empty(len(t));latent[0]=rng.normal()
                for j in range(1,len(t)):latent[j]=cal.rho*latent[j-1]+np.sqrt(1-cal.rho**2)*rng.normal()
                y=10+mu+np.sqrt(cal.dark_variance_mV2+cal.shot_coefficient_mV*np.maximum(mu,0))*latent+rng.normal(0,.1,len(t))
                initial=fit_waveform(t,y,template,base)
                metric=InnovationsMetric(initial['model_mV']-10,cal)
                seeds=[np.linspace(970.+shift,1150.+shift,k) for shift in (-10.,0.,10.)]
                results={}
                for name,global_rng in [('slsqp',None),('population_plus_slsqp',np.random.default_rng(262026702+len(rows)))]:
                    tick=time.monotonic();fit,trials=optimize(t,y-10,seeds,template,940.,1200.,metric,cfg,global_rng=global_rng)
                    results[name]=dict(times_ns=fit['times'].tolist(),charge_pe=(fit['amplitude']/(base.area_per_pe_mV_ns*template.unit_area_peak)).tolist(),
                        cost=fit['cost'],seconds=time.monotonic()-tick,evaluations=sum(r['evaluations'] for r in trials),
                        time_errors_ns=(fit['times']-true_t).tolist(),relative_charge_errors=(fit['amplitude']/a-1).tolist())
                rows.append(dict(known_count=k,gap_ns=gap,ratio=ratio,known_times_input=False,known_charge_input=False,results=results))
    save(OUTPUT/'optimizer-controls.json',dict(seed=262026701,rows=rows,
        known_count_input=True,known_times_input=False,counts_are_not_blind=True,
        limitation='Known two/three components, matched template and synthetic correlated noise; fixed trigger-region prior; DE uses additional evaluations. Not native unknown-K efficiency.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['physical','neural','controls','shifts','optimizer_controls']);parser.add_argument('--workers',type=int,default=2);a=parser.parse_args()
    if a.stage=='physical':physical(a.workers)
    else:globals()[a.stage]()
