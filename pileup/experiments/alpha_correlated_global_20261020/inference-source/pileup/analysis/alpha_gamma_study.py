"""Run reproducible alpha/gamma fitting checks and score frozen native campaigns."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime
import importlib.metadata
import json
from pathlib import Path
import platform
import time

import numpy as np
import awkward as ak

from .alpha_gamma_pileup_fit import TemplateBank, MixedFitConfig, fit_waveform, score_mixed
from .alpha_rate_campaign import digest, save
from .alpha_rate_data import export_window, open_root


def calibration(path):
    data=json.loads(Path(path).read_text());cfg=MixedFitConfig(**data['fit_config'])
    bank=TemplateBank(data['templates']['alpha']['readout'],
                      alpha_slow=data['templates']['alpha']['slow_fraction'],
                      gamma_slow=data['templates']['gamma']['slow_fraction'])
    return data,bank,cfg


def validate_root(row,inputs):
    root=Path(row['root']);macro=Path(row['macro'])
    assert digest(root)==row['root_sha256']
    assert digest(macro)==row['macro_sha256']
    with open_root(root) as f:
        assert all(f[k].num_entries==1 for k in ('RunInfo','WaveformEvent','SiPMEvent'))
        info=ak.to_list(f['RunInfo'].arrays())[0]
        assert info['macro_contents']==macro.read_text()
        assert info['random_seed']==row['seed']
        assert info['sipm_avalanche_base_seed']==row['sipm_seed']
        assert info['sipm_readout_base_seed']==row['readout_seed']
        assert info['alpha_p11b_rate_hz']==row['rate_hz']
        assert info['t_length_ps']==20000000
        p=f['PrimaryPhoton'].arrays(['primary_particle_name','primary_energy_MeV'],library='ak')
        names=np.asarray(ak.to_list(p['primary_particle_name']))
        energies=ak.to_numpy(p['primary_energy_MeV'])[names=='gamma']
        assert len(energies)>0 and np.all(energies>0)
        if row['gamma_model']=='mono662':assert np.allclose(energies,.661657,rtol=0,atol=1e-7)
        else:assert energies.min()>=.010-1e-7 and energies.max()<=5.+1e-7
    meta=export_window(root,inputs,expected=row,allowed_particles=('alpha','gamma'))
    truth=np.load(Path(inputs)/'truth'/f'{root.stem}.npz')
    counts={}
    for name in ('alpha','gamma'):
        selected=truth['particle_name']==name
        counts[name]=dict(injected=int(selected.sum()),nonzero=int(np.sum(selected&(truth['charge_pe']>0))),
                         above_500_pe=int(np.sum(selected&(truth['charge_pe']>=500))))
    return dict(root=str(root),actual_macro_verified=True,seed_verified=True,gamma_energies_verified=True,
                charge_conserved=True,counts=counts,metadata=meta)


def analyze_window(row,settings,inputs,output):
    output=Path(output);inputs=Path(inputs)
    for sub in ('windows','reconstructions'): (output/sub).mkdir(parents=True,exist_ok=True)
    stem=Path(row['root']).stem;destination=output/'windows'/f'{stem}.json'
    settings_hash=digest(settings);code_hash=digest(Path(__file__).with_name('alpha_gamma_pileup_fit.py'))
    if destination.exists():
        old=json.loads(destination.read_text())
        if old['settings_sha256']!=settings_hash or old['fitter_sha256']!=code_hash:
            raise ValueError('Existing result uses different settings/code; use a new output directory')
        return old
    checked=validate_root(row,inputs)
    _,bank,cfg=calibration(settings)
    waveform=np.load(inputs/'waveforms'/f'{stem}.npz')
    started=time.perf_counter()
    fit=fit_waveform(waveform['time_ns'],waveform['analog_mV'],bank,cfg)
    elapsed=time.perf_counter()-started
    # Truth is deliberately loaded only after reconstruction is finished.
    truth=dict(np.load(inputs/'truth'/f'{stem}.npz'))
    gate=checked['metadata']['gate_ns'];scored=score_mixed(fit,truth,gate)
    np.savez_compressed(output/'reconstructions'/f'{stem}.npz',time_ns=waveform['time_ns'],model_mV=fit.pop('model_mV'),
                        times_ns=fit['times_ns'],charge_pe=fit['charge_pe'],particle_name=fit['particle_name'])
    result=dict(row=row,settings_sha256=settings_hash,fitter_sha256=code_hash,
                fit_seconds=elapsed,fit=fit,score=scored,validation={k:v for k,v in checked.items() if k!='metadata'},
                source_rates_hz=checked['metadata']['source_rates_hz'],
                waveform_sha256=digest(inputs/'waveforms'/f'{stem}.npz'),
                truth_sha256=digest(inputs/'truth'/f'{stem}.npz'))
    save(destination,result)
    print('FIT',row['label'],'targets',scored['target_count'],'recovered',scored['charge_and_type_recovered_count'],
          'extra',scored['unmatched_predictions'],f'{elapsed:.2f}s',flush=True)
    return result


def controlled(settings,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    _,bank,cfg=calibration(settings);t=np.arange(.25,12000.,2.)
    rows=[];seed=261025000
    # These are ideal template + electronics-noise controls. Native shared-cell
    # simulations are assessed separately and never made by adding reset pulses.
    specifications=[('empty',[],[],[],None,None)]
    for particle in (0,1):
        for charge in (700.,2000.,6000.):
            specifications.append(('single',[5000.],[particle],[charge],None,None))
    for gap in (0.,10.,20.,40.,80.,200.,1000.):
        for ratio in (.25,1.,2.):
            for order in ((0,1),(1,0)):
                charges=[4000.,4000.*ratio] if order[0]==0 else [4000.*ratio,4000.]
                specifications.append(('double',[5000.,5000.+gap],list(order),charges,gap,ratio))
    for index,(kind,times,types,charges,gap,ratio) in enumerate(specifications):
        for repetition in range(3):
            rng=np.random.default_rng(seed+3*index+repetition)
            values=bank.matrix(t,times,types,cfg)@charges+cfg.baseline_mV+rng.normal(0.,.1,len(t))
            fit=fit_waveform(t,values,bank,cfg)
            truth=dict(first_active_ns=np.asarray(times),emission_ns=np.asarray(times),charge_pe=np.asarray(charges),
                       particle_name=bank.names[np.asarray(types,dtype=int)])
            result=score_mixed(fit,truth,[4000.,8000.])
            rows.append(dict(kind=kind,gap_ns=gap,gamma_to_alpha_charge_ratio=ratio,
                seed=seed+3*index+repetition,truth_times_ns=times,truth_particle_name=truth['particle_name'].tolist(),
                truth_charge_pe=charges,fit={k:v for k,v in fit.items() if k!='model_mV'},score=result))
        if index%10==0:print('CONTROL',index+1,'/',len(specifications),flush=True)
    summary=dict(input_model='Calibrated ideal template sum + independent 0.1 mV electronics noise; no photon/SiPM fluctuations',
        repeats_per_case=3,zero_gap_policy='Blind baseline merges onsets closer than 8 ns; coincident components need not be recoverable',
        calibration_sha256=digest(settings),fitter_sha256=digest(Path(__file__).with_name('alpha_gamma_pileup_fit.py')),rows=rows)
    save(output/'controlled.json',summary);return summary


def summarize(rows,output):
    output=Path(output);groups=[];rng=np.random.default_rng(261026000)
    count_keys=('target_count','central_injected','central_nonzero_charge','matched_count','correct_type_count',
                'recovered_count','missed_count','wrong_type_count')
    for model in ('brems','mono662'):
        for rate in (1000000,10000000):
            subset=[r for r in rows if r['row']['gamma_model']==model and r['row']['rate_hz']==rate]
            if not subset:continue
            for particle in ('alpha','gamma'):
                records=[r['score']['by_particle'][particle] for r in subset]
                counts={k:sum(r[k] for r in records) for k in count_keys}
                q=np.array([v for r in records for v in r['relative_charge_errors']])
                dt=np.array([v for r in records for v in r['time_errors_ns']])
                n=counts['target_count'];matched=counts['matched_count']
                arrays={k:np.array([r[k] for r in records]) for k in count_keys}
                draws=rng.integers(0,len(subset),size=(3000,len(subset)))
                intervals={}
                for name,num,den in [('detection','matched_count','target_count'),('classification','correct_type_count','matched_count'),
                                     ('recovery','recovered_count','target_count')]:
                    denominator=arrays[den][draws].sum(axis=1);numerator=arrays[num][draws].sum(axis=1)
                    valid=denominator>0
                    intervals[name]=np.quantile(numerator[valid]/denominator[valid],[.025,.975]).tolist() if valid.any() else None
                group=dict(gamma_model=model,alpha_injection_MHz=rate/1e6,gamma_injection_MHz=10.,particle=particle,
                    windows=len(subset),effective_time_us=10*len(subset),**counts,
                    signal_rate_MHz=n/(10*len(subset)),detection_efficiency=matched/n if n else None,
                    conditional_classification=counts['correct_type_count']/matched if matched else None,
                    full_recovery_efficiency=counts['recovered_count']/n if n else None,
                    mean_relative_charge_error=float(q.mean()) if len(q) else None,
                    std_relative_charge_error=float(q.std(ddof=1)) if len(q)>1 else None,
                    charge_error_percentiles=np.quantile(q,[.16,.5,.84]).tolist() if len(q) else None,
                    time_bias_ns=float(dt.mean()) if len(dt) else None,
                    time_rmse_ns=float(np.sqrt(np.mean(dt**2))) if len(dt) else None,
                    unmatched_predictions_class=sum(r['score']['unmatched_predictions_by_class'][particle] for r in subset),
                    cluster_bootstrap_95=intervals)
                groups.append(group)
    summary=dict(created=datetime.now().astimezone().isoformat(),groups=groups,windows=len(rows),
        input='Native full-window shared-cell SiPM analog_mV; fixed 10 mV baseline; no ADC range cut',
        success='Q_truth>=500 pe, one-to-one class-blind time match within 10 ns, correct class and charge error<=20%',
        charge_error_statistics='Conditional on time matching, including incorrectly classified events',
        uncertainty='3000 resamples of whole independent windows within each gamma/rate condition; descriptive pilot intervals',
        environment=dict(python=platform.python_version(),**{k:importlib.metadata.version(k) for k in ('numpy','scipy','uproot','matplotlib')}),
        totals={k:sum(r['score'][k] for r in rows) for k in ('target_count','matched_count','missed_count',
            'correct_type_count','charge_and_type_recovered_count','unmatched_predictions','split_flags','merge_flags')},
        fit_seconds=sum(r['fit_seconds'] for r in rows))
    save(output/'summary.json',summary)
    save(output/'window-results.json',rows)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings',type=Path,required=True)
    parser.add_argument('--inputs',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--stage',choices=('checks','production','controlled'),required=True)
    parser.add_argument('--experiments',nargs='+',type=Path,default=[])
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    if args.stage=='controlled':controlled(args.settings,args.output);return
    rows=[]
    for folder in args.experiments:
        manifest=json.loads((folder/'manifest.json').read_text())
        plans=manifest['checks' if args.stage=='checks' else 'runs']
        assert plans and all(r['status']=='completed' for r in plans)
        if args.stage=='production':assert len(plans)==len(manifest['planned_runs'])
        for row in plans:rows.append(analyze_window(row,args.settings,args.inputs,args.output))
    summary=summarize(rows,args.output)
    print('SUMMARY',args.output/'summary.json',summary['totals'],flush=True)


if __name__=='__main__':main()
