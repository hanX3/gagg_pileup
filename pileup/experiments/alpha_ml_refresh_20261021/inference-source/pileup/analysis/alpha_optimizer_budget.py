"""Supplemental known-K search diagnostic with equal variable-projection work.

This post-freeze diagnostic uses NEW synthetic noise, not new native selection.
Both strategies spend exactly the same number of actual whitened NNLS profile
evaluations, including their common coarse starts and polishing. All candidate
charges/times are inferred. K and a fixed trigger region are explicitly given.
"""
from dataclasses import replace
from datetime import datetime
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import LinearConstraint, differential_evolution, minimize

from .alpha_correlated_fit import InnovationsMetric, ResidualCalibration, projected
from .alpha_correlated_study import OUTPUT, EXPERIMENT
from .alpha_full_signal_study import make_template
from .alpha_pileup_fit import FitConfig, fit_waveform
from .alpha_rate_campaign import digest, save


class BudgetStop(Exception):
    pass


class Objective:
    def __init__(self, t, values, template, metric, budget, lower=940., upper=1200., separation=3.):
        self.t, self.values, self.template, self.metric = t, values, template, metric
        self.budget = self.stop_at = budget
        self.lower, self.upper, self.separation = lower, upper, separation
        self.calls = self.invalid_calls = 0; self.best = None

    def evaluate(self, times):
        times = np.asarray(times, float)
        if np.any(times < self.lower-1e-6) or np.any(times > self.upper+1e-6) or np.any(np.diff(times) < self.separation-1e-5):
            self.invalid_calls += 1
            return 1e9, np.zeros(len(times))
        if self.calls >= self.stop_at:
            raise BudgetStop()
        value, gradient, amplitude = projected(self.t, self.values, times, self.template, self.metric)
        self.calls += 1
        if self.best is None or value < self.best['cost']:
            self.best = dict(cost=value, times_ns=times.copy(), amplitude=amplitude.copy())
        return value, gradient


def local(objective, seed, scale):
    seed=np.asarray(seed,float);k=len(seed);width=objective.upper-objective.lower
    difference=np.eye(k,k,1)[:-1]-np.eye(k)[:-1]
    constraints=(LinearConstraint(difference,objective.separation/width,np.inf),) if k>1 else ()
    def function(z):
        value,gradient=objective.evaluate(objective.lower+width*z)
        return value/scale,gradient*width/scale
    result=minimize(function,(seed-objective.lower)/width,jac=True,method='SLSQP',
        bounds=[(0.,1.)]*k,constraints=constraints,options=dict(maxiter=35,ftol=1e-9))
    return bool(result.success)


def search(t,values,template,metric,k,*,seed,budget=500,population=False):
    objective=Objective(t,values,template,metric,budget)
    rng=np.random.default_rng(seed);tick=time.monotonic();trials=0;de_calls=0
    try:
        # Identical three fixed coarse starts, without true onsets/charges.
        starts=[np.linspace(970.+shift,1150.+shift,k) for shift in (-10.,0.,10.)]
        scale=max(objective.evaluate(starts[0])[0],1.)
        for start in starts:
            local(objective,start,scale);trials+=1
        common_calls=objective.calls
        if population and objective.calls<budget:
            # Spend 70% of remaining profile work on DE, then polish and fill
            # with local restarts. Invalid time vectors do not spend NNLS work.
            objective.stop_at=objective.calls+max(1,int(.7*(budget-objective.calls)))
            before=objective.calls;width=objective.upper-objective.lower
            def function(z):
                return objective.evaluate(objective.lower+width*np.sort(z))[0]/scale
            initial=np.clip((objective.best['times_ns']-objective.lower)/width,1e-12,1-1e-12)
            try:
                differential_evolution(function,[(0.,1.)]*k,popsize=6,maxiter=10000,
                    rng=rng,polish=False,tol=0.,atol=0.,workers=1,updating='immediate',x0=initial)
            except BudgetStop:
                pass
            de_calls=objective.calls-before
            objective.stop_at=budget
            local(objective,objective.best['times_ns'],scale);trials+=1
        while objective.calls<budget:
            start=np.sort(rng.uniform(objective.lower,objective.upper,k))
            if np.min(np.diff(start))<objective.separation:continue
            local(objective,start,scale);trials+=1
    except BudgetStop:
        pass
    assert objective.calls==budget and objective.best is not None
    result={k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in objective.best.items()}
    result.update(profile_evaluations=objective.calls,invalid_time_evaluations=objective.invalid_calls,
        de_profile_evaluations=de_calls,completed_local_starts=trials,seconds=time.monotonic()-tick)
    return result


def main():
    folder=OUTPUT/'optimizer-budget';folder.mkdir(exist_ok=True)
    protocol_path=folder/'protocol.json'
    protocol=dict(created=datetime.now().astimezone().isoformat(),
        scope='Supplement added after seeing the earlier coarse-start known-K diagnostic; new noise, not native test tuning',
        noise_seed=262027700,optimizer_seeds=[262027711,262027712,262027713],
        known_counts=[2,3],gaps_ns=[4,8,12,20,40,80],ratios=[.125,1.,4.],
        profile_evaluations_per_strategy=500,common_starts='Three fixed 970..1150 ns grids, shifts -10/0/+10 ns',
        search_bounds_ns=[940.,1200.],minimum_separation_ns=3.,
        strategies=['random_multistart_slsqp','de_then_slsqp_and_restart_fill'],
        de_fraction_of_remaining_profile_budget=.7,known_count_input=True,known_times_or_charges_input=False,
        native_inference_or_selection_changed=False,source_sha256=digest(Path(__file__)),
        covariance_sha256=digest(OUTPUT/'calibration/covariance.json'),
        baseline_settings_sha256=digest(OUTPUT/'baseline-settings.json'),
        equality_scope='Actual whitened NNLS profile evaluations; not equal wall time, invalid-vector overhead differs')
    if protocol_path.exists():
        previous=json.loads(protocol_path.read_text());assert previous['source_sha256']==protocol['source_sha256'];protocol=previous
    else:save(protocol_path,protocol)
    settings=json.loads((OUTPUT/'baseline-settings.json').read_text());template=make_template(settings);base=FitConfig(**settings['fit_config'])
    cal=ResidualCalibration(**json.loads((OUTPUT/'calibration/covariance.json').read_text())['parameters'])
    t=np.arange(0.,3000.,2.)+.137;rng=np.random.default_rng(protocol['noise_seed']);rows=[];case=0
    for k in protocol['known_counts']:
        for gap in protocol['gaps_ns']:
            for ratio in protocol['ratios']:
                true_t=1000.3+np.arange(k)*gap;true_q=np.array([4000.]+[4000.*ratio]*(k-1))
                amplitudes=true_q*base.area_per_pe_mV_ns*template.unit_area_peak
                mu=template(t[:,None]-true_t)@amplitudes
                latent=np.empty(len(t));latent[0]=rng.normal()
                for j in range(1,len(t)):latent[j]=cal.rho*latent[j-1]+np.sqrt(1-cal.rho**2)*rng.normal()
                y=10+mu+np.sqrt(cal.dark_variance_mV2+cal.shot_coefficient_mV*np.maximum(mu,0))*latent+rng.normal(0,.1,len(t))
                pilot=fit_waveform(t,y,template,base);metric=InnovationsMetric(pilot['model_mV']-10,cal)
                np.savez_compressed(folder/f'case-{case:02d}.npz',time_ns=t,analog_mV=y,true_onsets_ns=true_t,true_charge_pe=true_q)
                for optimizer_seed in protocol['optimizer_seeds']:
                    path=folder/f'case-{case:02d}-seed-{optimizer_seed}.json'
                    if path.exists():rows.append(json.loads(path.read_text()));continue
                    predictions={}
                    for name,population in zip(protocol['strategies'],(False,True)):
                        predictions[name]=search(t,y-10,template,metric,k,seed=optimizer_seed+case*97,population=population,
                            budget=protocol['profile_evaluations_per_strategy'])
                    # Labels enter evaluation only, after both searches above.
                    results={}
                    for name,prediction in predictions.items():
                        charge=np.asarray(prediction.pop('amplitude'))/(base.area_per_pe_mV_ns*template.unit_area_peak)
                        dt=np.asarray(prediction['times_ns'])-true_t;dq=charge/true_q-1
                        results[name]=dict(**prediction,charge_pe=charge.tolist(),time_errors_ns=dt.tolist(),
                            relative_charge_errors=dq.tolist(),all_components_complete=bool(np.all(abs(dt)<=10)&np.all(abs(dq)<=.2)))
                    row=dict(case=case,known_count=k,gap_ns=gap,ratio=ratio,optimizer_seed=optimizer_seed,results=results)
                    save(path,row);rows.append(row)
                case+=1
                print('BUDGET_DIAGNOSTIC_CASE',case,'/36',flush=True)
    summary=[]
    for k in protocol['known_counts']:
        selected=[r for r in rows if r['known_count']==k]
        for name in protocol['strategies']:
            summary.append(dict(known_count=k,strategy=name,independent_noise_cases=18,optimizer_seeds_per_case=3,
                runs=len(selected),all_components_complete=sum(r['results'][name]['all_components_complete'] for r in selected),
                total_profile_evaluations=sum(r['results'][name]['profile_evaluations'] for r in selected),
                summed_seconds=sum(r['results'][name]['seconds'] for r in selected)))
    assert len(rows)==108
    save(OUTPUT/'optimizer-budget-controls.json',dict(protocol=protocol,rows=rows,summary=summary,
        limitation='Known K and trigger region, matched mean and synthetic noise. Three optimizer seeds share each noise record; 108 runs are not 108 independent noise samples. No native method/threshold changes.'))
    snapshot=EXPERIMENT/'supplemental-search-source';snapshot.mkdir(exist_ok=True)
    import shutil
    destination=snapshot/Path(__file__).name
    if not destination.exists():shutil.copyfile(Path(__file__),destination)
    assert digest(destination)==protocol['source_sha256']
    print('EQUAL_PROFILE_BUDGET_COMPLETE',summary,flush=True)


if __name__=='__main__':main()
