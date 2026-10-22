"""Frozen independent-window inference and separate scoring for the full hybrid."""
from collections import defaultdict
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
import argparse
import json
import time
import numpy as np
from pileup.analysis.alpha_gamma_full_prior_fit import FullPriorConfig, fit_full_prior
from pileup.analysis.alpha_gamma_model_selection_fit import ModelSelectionConfig, fit_model_selection
from pileup.analysis.alpha_gamma_prior_study import digest, save, measure
from pileup.analysis.alpha_gamma_pileup_fit import score_mixed

METHODS = ['baseline','unet','hybrid_none','hybrid_both','class_only']
CONFIGS = dict(hybrid_none=FullPriorConfig(charge_strength=0.,class_strength=0.),
               hybrid_both=FullPriorConfig(charge_strength=10.,class_strength=10.),
               class_only=FullPriorConfig())


def infer(time_ns, voltage_mV, calibration, model, network_config, thresholds):
    """No truth, source conditions, known count or metadata enter this interface."""
    from pileup.analysis.alpha_gamma_ml import voltage_features, predict_features, decode
    started=time.perf_counter()
    dense=predict_features(model,voltage_features(time_ns,voltage_mV,calibration),network_config)
    neural=decode(time_ns,dense,thresholds,network_config)
    fits={'unet':neural};seconds={'unet':time.perf_counter()-started}
    started=time.perf_counter()
    fits['baseline']=fit_model_selection(time_ns,voltage_mV,calibration,ModelSelectionConfig(addition_cost=1.))
    seconds['baseline']=time.perf_counter()-started
    for name,config in CONFIGS.items():
        started=time.perf_counter()
        fits[name]=fit_full_prior(time_ns,voltage_mV,calibration,neural,dense,fits['baseline'],config)
        seconds[name]=time.perf_counter()-started
    return fits,dense,seconds


def fit_case(work,stem):
    import torch
    from pileup.analysis.alpha_gamma_ml_study import load_model
    work=Path(work);torch.set_num_threads(1)
    case=json.loads((work/'test-cases'/f'{stem}.json').read_text())
    settings=json.loads((work/'test-settings.json').read_text())
    path=work/case['waveform']
    assert digest(path)==case['waveform_sha256']
    assert digest(work/'calibration.json')==settings['calibration_sha256']
    model,cfg,frozen=load_model(work/'models/unet')
    assert frozen['model_sha256']==settings['neural_model_sha256']
    assert frozen['selected']['thresholds']==settings['neural_thresholds']
    assert {m:asdict(c) for m,c in CONFIGS.items()}==settings['hybrid_configs']
    wave=dict(np.load(path));t,y=wave['time_ns'],wave['analog_mV']
    fits,dense,seconds=infer(t,y,json.loads((work/'calibration.json').read_text()),model,cfg,settings['neural_thresholds'])
    densepath=work/'test-outputs/dense'/f'{stem}.npz';densepath.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(densepath,**dense)
    for method in METHODS:
        destination=work/'test-outputs/fits'/method/f'{stem}.json'
        if destination.exists():raise FileExistsError(destination)
        fitted=dict(fits[method]);model_voltage=fitted.pop('model_mV',None)
        destination.parent.mkdir(parents=True,exist_ok=True)
        if model_voltage is not None:
            np.savez_compressed(destination.with_suffix('.npz'),time_ns=t,model_mV=model_voltage)
        save(destination,dict(stem=stem,method=method,fit=fitted,fit_seconds=seconds[method],
            pipeline_seconds=seconds[method]+(seconds['baseline']+seconds['unet'] if method in CONFIGS else 0.),
            waveform_sha256=case['waveform_sha256'],calibration_sha256=settings['calibration_sha256'],
            neural_model_sha256=frozen['model_sha256'],settings_sha256=digest(work/'test-settings.json'),
            fit_truth_inputs=[],source_rate_input=False,truth_file_access=False))
    print(json.dumps(dict(stem=stem,seconds=seconds,reported={m:len(fits[m]['times_ns']) for m in METHODS})),flush=True)


def report(project,folder,work):
    project,folder,work=map(Path,(project,folder,work))
    manifest=json.loads((folder/'manifest.json').read_text())
    settings=json.loads((folder/'frozen-fit-settings.json').read_text())
    assert digest(work/'test-settings.json')==digest(folder/'frozen-fit-settings.json')
    rows=manifest['runs'];assert len(rows)==16 and all(r['status']=='completed' for r in rows)
    output=work/'test-outputs'
    assert all((output/'fits'/m/(Path(r['root']).stem+'.json')).exists() for r in rows for m in METHODS)
    by_method=defaultdict(list);strata=defaultdict(list);events=[];numerics=[]
    for row in rows:
        stem=Path(row['root']).stem;strata[(row['gamma_model'],row['rate_hz'])].append(stem)
        meta=json.loads((work/'test-inputs/metadata'/f'{stem}.json').read_text())
        assert meta['root_sha256']==row['root_sha256'] and meta['seed']==row['seed']
        gate=meta['gate_ns']
        with np.load(work/'test-inputs/truth'/f'{stem}.npz') as data:
            truth={k:data[k] for k in ('first_active_ns','emission_ns','particle_name','charge_pe')}
        wave=dict(np.load(work/'test-inputs/waveforms'/f'{stem}.npz'))
        for method in METHODS:
            record=json.loads((output/'fits'/method/f'{stem}.json').read_text())
            assert record['waveform_sha256']==digest(work/'test-inputs/waveforms'/f'{stem}.npz')
            assert record['settings_sha256']==digest(folder/'frozen-fit-settings.json')
            fitted=record['fit'];score=score_mixed(fitted,truth,gate)
            entry=dict(row=row,fit=fitted,truth=truth,score=score,gate_ns=gate)
            by_method[method].append(entry)
            save(output/'scores'/method/f'{stem}.json',dict(row=row,method=method,score=score,gate_ns=gate))
            if method in CONFIGS:
                history=fitted.get('weighted_cost_history',[])
                numerics.append(dict(stem=stem,method=method,status=fitted['status'],
                    cost_nonincreasing=bool(np.all(np.diff(history)<=1e-6)),
                    pair_trials=fitted.get('pair_trials',0),optimizer_status_counts=fitted.get('optimizer_status_counts',{}),
                    additions=len(fitted.get('residual_additions',[]))))
        physical=np.flatnonzero(np.isfinite(truth['first_active_ns'])&(truth['charge_pe']>=50.))
        for i in by_method['baseline'][-1]['score']['target_indices']:
            onset,charge=float(truth['first_active_ns'][i]),float(truth['charge_pe'][i])
            others=physical[physical!=i]
            gap=float(np.min(abs(truth['first_active_ns'][others]-onset))) if len(others) else None
            mask=(wave['time_ns']>=onset-24.)&(wave['time_ns']<onset-4.)
            history=max(0.,float(wave['analog_mV'][mask].mean()-settings['baseline_mV'])) if mask.any() else 0.
            event=dict(stem=stem,truth_index=i,particle=str(truth['particle_name'][i]),true_charge_pe=charge,
                preceding_voltage_mV=history,nearest_other_q50_ns=gap,strata=dict(
                    charge='500-1000' if charge<1000 else '1000-2000' if charge<2000 else '>=2000',
                    history='<50' if history<50 else '50-200' if history<200 else '>=200',
                    spacing='<50' if gap is not None and gap<50 else '50-150' if gap is not None and gap<150 else
                            '150-500' if gap is not None and gap<500 else '>=500'),outcomes={})
            for method in METHODS:
                result=by_method[method][-1];pair=dict(result['score']['matched_pairs']);prediction=None
                outcome='missed_at_report_threshold'
                if i in pair:
                    j=pair[i];fit=result['fit'];kind,q=fit['particle_name'][j],fit['charge_pe'][j]
                    outcome='wrong_type' if kind!=event['particle'] else 'bad_charge' if abs(q/charge-1)>.2 else 'recovered'
                    prediction=dict(index=j,time_ns=fit['times_ns'][j],particle=kind,charge_pe=q)
                event['outcomes'][method]=dict(outcome=outcome,prediction=prediction)
            events.append(event)
    pooled=[dict(method=m,particle=p,**measure(by_method[m],p)) for m in METHODS for p in ('alpha','gamma')]
    groups=[]
    for method in METHODS:
        for model,rate in sorted(strata):
            subset=[r for r in by_method[method] if (r['row']['gamma_model'],r['row']['rate_hz'])==(model,rate)]
            groups.extend(dict(method=method,gamma_model=model,rate_hz=rate,particle=p,windows=len(subset),
                               **measure(subset,p)) for p in ('alpha','gamma'))
    indexed={m:{Path(r['row']['root']).stem:r for r in by_method[m]} for m in METHODS}
    rng=np.random.default_rng(settings['test_bootstrap_seed'])
    draws=[[str(s) for cell in strata.values() for s in rng.choice(cell,len(cell),replace=True)] for _ in range(2000)]
    paired=[]
    for reference in ('baseline','hybrid_none','hybrid_both'):
        for particle in ('alpha','gamma'):
            n={s:indexed[reference][s]['score']['by_particle'][particle]['target_count'] for s in indexed[reference]}
            delta={s:indexed['class_only'][s]['score']['by_particle'][particle]['recovered_count']-
                     indexed[reference][s]['score']['by_particle'][particle]['recovered_count'] for s in indexed[reference]}
            samples=[sum(delta[s] for s in draw)/sum(n[s] for s in draw) for draw in draws if sum(n[s] for s in draw)]
            paired.append(dict(a='class_only',b=reference,particle=particle,delta_count=sum(delta.values()),
                delta_recovery=sum(delta.values())/sum(n.values()),paired_window_p025_p50_p975=np.quantile(samples,[.025,.5,.975]).tolist()))
    diagnostics={m:{k:sum(r['score'][k] for r in by_method[m]) for k in
        ('target_count','matched_count','missed_count','unmatched_predictions','split_flags','merge_flags')} for m in METHODS}
    bins=defaultdict(lambda:defaultdict(int))
    for event in events:
        for method,item in event['outcomes'].items():
            for feature,interval in event['strata'].items():
                bins[(method,event['particle'],feature,interval)]['n']+=1
                bins[(method,event['particle'],feature,interval)][item['outcome']]+=1
    checks=[];gate_rule=settings['development_gate']
    for particle in ('alpha','gamma'):
        a=next(r for r in pooled if r['method']=='class_only' and r['particle']==particle)
        b=next(r for r in pooled if r['method']=='baseline' and r['particle']==particle)
        checks.extend([dict(check=particle+'_pooled_recovery',passed=a['recovered_count']>=b['recovered_count'],
                            candidate=a['recovered_count'],reference=b['recovered_count']),
                       dict(check=particle+'_precision_margin',passed=a['full_precision']>=b['full_precision']-gate_rule['precision_margin_absolute'],
                            candidate=a['full_precision'],reference=b['full_precision'])])
    extra_limit=diagnostics['baseline']['unmatched_predictions']*gate_rule['extra_multiplier']+gate_rule['extra_offset']
    checks.append(dict(check='extra_predictions_budget',passed=diagnostics['class_only']['unmatched_predictions']<=extra_limit,
                       candidate=diagnostics['class_only']['unmatched_predictions'],limit=extra_limit))
    for a in (r for r in groups if r['method']=='class_only' and r['target_count']>=gate_rule['cell_minimum_targets']):
        b=next(r for r in groups if r['method']=='baseline' and (r['particle'],r['gamma_model'],r['rate_hz'])==(a['particle'],a['gamma_model'],a['rate_hz']))
        checks.append(dict(check=f"cell_{a['gamma_model']}_{a['rate_hz']}_{a['particle']}",
            passed=a['full_recovery']>=b['full_recovery']-gate_rule['cell_recovery_margin_absolute'],
            candidate=a['full_recovery'],reference=b['full_recovery']))
    summary=dict(created=datetime.now().astimezone().isoformat(),windows=16,methods=METHODS,pooled=pooled,groups=groups,
        paired_contrasts=paired,diagnostics=diagnostics,numerical_diagnostics=numerics,
        event_strata=[dict(method=m,particle=p,feature=f,bin=b,**v) for (m,p,f,b),v in sorted(bins.items())],
        frozen_screen=dict(passed=all(c['passed'] for c in checks),checks=checks),
        settings_sha256=digest(folder/'frozen-fit-settings.json'),
        scope='16 native independent windows, settings and seeds fixed before generation, all fits before scoring',
        target='Q>=500 pe; match<=10 ns; correct class; |dQ/Q|<=20%; emission [5,15) us',
        limits=['Only four windows per source condition; limited gamma statistics',
                'Window bootstrap excludes detector/model systematic errors','Uncalibrated class scores and neural-dependent initial shape Q',
                'Charge is not deposited energy','No tuning or retraining after test generation'])
    save(output/'summary.json',summary);save(output/'events.json',events)
    print(json.dumps(dict(pooled=pooled,frozen_screen=summary['frozen_screen'])),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='action',required=True)
    fit=sub.add_parser('fit');fit.add_argument('--work',type=Path,required=True);fit.add_argument('--stem',required=True)
    scoring=sub.add_parser('report')
    for name in ('project','folder','work'):scoring.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    if args.action=='fit':fit_case(args.work,args.stem)
    else:report(args.project,args.folder,args.work)


if __name__=='__main__':main()
