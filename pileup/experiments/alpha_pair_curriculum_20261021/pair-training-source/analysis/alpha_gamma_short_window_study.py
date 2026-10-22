"""Calibrate, develop and evaluate a frozen isolated-prefix reconstruction."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from .alpha_gamma_short_window import (
    METHODS, PrefixCalibration, PrefixConfig, build_calibration, find_trigger, fit_prefix,
)

PROJECT = Path(__file__).resolve().parents[1]
POSTS = (50.,75.,100.,150.,200.,250.,300.,400.,500.,750.,1000.)
PRES = (20.,50.,100.,200.)
SOURCES = (
    ('pileup/experiments/alpha_gamma_regularized_fit_20261015',
     'jupyter/artifacts/alpha_gamma_regularized_20261015/inputs','calibration','validation_single'),
    ('pileup/experiments/alpha_gamma_context_calibration_20261016',
     'jupyter/artifacts/alpha_gamma_context_calibration_20261016/inputs','calibration_single','validation_single'),
)


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def save(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def records_from_sources():
    groups = {'calibration':[],'development':[]}
    for folder, inputs, training_role, development_role in SOURCES:
        folder,inputs=PROJECT/folder,PROJECT/inputs
        manifest=json.loads((folder/'manifest.json').read_text())
        for row in manifest['runs']:
            role = 'calibration' if row['role']==training_role else 'development' if row['role']==development_role else None
            if role is None:
                continue
            stem=Path(row['root']).stem
            meta=inputs/'metadata'/f'{stem}.json'
            assert row['status']=='completed' and json.loads(meta.read_text())['root_sha256']==row['root_sha256']
            groups[role].append({'root':row['root'],'seed':row['seed'],'particle':row['particle'],
                'waveform':str((inputs/'waveforms'/f'{stem}.npz').relative_to(PROJECT)),
                'truth':str((inputs/'truth'/f'{stem}.npz').relative_to(PROJECT)),
                'metadata':str(meta.relative_to(PROJECT)), 'source_manifest':str((folder/'manifest.json').relative_to(PROJECT)),
                'root_original_sha256':row['root_sha256']})
    roots=[r['root'] for group in groups.values() for r in group]
    assert len(roots)==len(set(roots))
    assert not set(r['seed'] for r in groups['calibration']) & set(r['seed'] for r in groups['development'])
    return groups


def prepare(output):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    groups=records_from_sources()
    for group in groups.values():
        for row in group:
            row['cached_sha256']={key:digest(PROJECT/row[key]) for key in ('waveform','truth','metadata')}
    meta=json.loads((PROJECT/groups['calibration'][0]['metadata']).read_text())
    readout=meta['config']['readout']
    for row in [r for group in groups.values() for r in group]:
        actual=json.loads((PROJECT/row['metadata']).read_text())['config']['readout']
        assert actual==readout, row['root']
    cfg=PrefixConfig()
    save(output/'input-records.json',groups)
    save(output/'protocol.json',{
        'created':datetime.now().astimezone().isoformat(),'post_ns':POSTS,'pre_ns':PRES,'methods':METHODS,
        'primary_method':'template_gls','primary_adaptive_pre_ns':50.,'charge_success_tolerance':.05,
        'requested_point_success_fraction':.95,'matching_trigger_tolerance_ns':60.,
        'calibration_rows':len(groups['calibration']),'development_rows':len(groups['development']),
        'input':'Full native analog voltage, adequate range; common new input for ALL compared methods',
        'scope':'At most one isolated signal, unknown type/time/charge; no pulse-count estimation or pileup performance claim',
        'target':'Realized photon-descendant avalanche charge excluding dark chains, in full-cell pe; not source/deposited energy',
        'no_future':'Trigger uses only a pre_ns lookback and independent fixed threshold; fit only [trigger-pre,trigger+post]',
        'calibration_onsets':'First active photon and class/Q labels allowed ONLY in independent calibration',
        'shape':'Charge-binned native mean waveform; interpolation at a voltage-only pilot Q; not a truth type input',
        'covariance':'Same pilot-conditioned C for both class hypotheses; covariance restricted BEFORE Cholesky; no Gaussian probability claim',
        'prefix_lookback_note':'Onset search covers min(Tpre,60 ns); trigger confirmation delay is included in Tpre; readout/filter is native and adds no new software delay',
        'adaptive':'Development per-event worst-window standardized error and residual guard; empirical rule, no sequential coverage guarantee',
        'precision_note':'5%/95% is an engineering target evaluated per charge band; low-charge conditions may not reach it even at long windows',
        'known_limits':['Calibration includes several source positions, no dense-history calibration',
                        'Fixed95 physics model keeps original tau inputs; native templates independently constrain full effective shape',
                        'Population covariance and shape interpolation approximate nonstationary SiPM statistics',
                        'New test will be created only after configuration/code freeze; small sample cannot certify a universal 95% bound']})
    start=time.perf_counter()
    cal=build_calibration(groups['calibration'],readout,13.618501389,output/'calibration',cfg)
    print('CALIBRATED',cal.meta['used_training_rows'],cal.meta['classes'],'seconds',round(time.perf_counter()-start,2),flush=True)
    return groups


def _worker(task):
    row,calibration_folder,output,pre_values,post_values=task
    cal=PrefixCalibration.load(calibration_folder)
    with np.load(PROJECT/row['waveform']) as waveform:
        t,y=waveform['time_ns'],waveform['analog_mV']
        results=[];started=time.perf_counter()
        for pre in pre_values:
            trigger=find_trigger(t,y,pre,cal.config.threshold_mV)
            if trigger is None:
                results.extend({'detected':False,'pre_ns':pre,'post_ns':post,'methods':{}} for post in post_values)
            else:
                for post in post_values:
                    results.append(fit_prefix(t,y,cal,pre_ns=pre,post_ns=post,trigger=trigger))
    # Truth is opened after all reconstructions of this event have completed.
    with np.load(PROJECT/row['truth']) as truth:
        q=float(truth['charge_pe'][0]);name=str(truth['particle_name'][0]);onset=float(truth['first_active_ns'][0])
        score_truth={'charge_pe':q,'particle':name,'onset_first_active_ns':onset if math.isfinite(onset) else None,
                     'edep_MeV':float(truth['edep_MeV'][0])}
    for fit in results:
        matched=bool(fit['detected'] and math.isfinite(onset) and abs(fit['trigger']['time_ns']-onset)<=60.)
        fit['evaluation']={'matched_signal':matched,'target_above_500pe':q>=500.,
                           'no_deposited_signal':q==0.}
        for method,r in fit['methods'].items():
            r['evaluation']={'relative_charge_error':r['charge_pe']/q-1 if q>0 and matched else None,
                'correct_type':r['particle']==name if r['particle'] is not None and matched else False,
                'onset_error_ns':fit['trigger']['time_ns']+r['onset_relative_ns']-onset if matched else None}
    result={'row':row,'truth_for_evaluation_only':score_truth,'prefixes':results,
            'seconds':time.perf_counter()-started,'inference_truth_inputs':[],
            'fitter_sha256':digest(PROJECT/'analysis/alpha_gamma_short_window.py'),
            'study_sha256':digest(Path(__file__)),
            'calibration_sha256':digest(Path(calibration_folder)/'calibration-arrays.npz')}
    destination=Path(output)/(Path(row['root']).stem+'.json')
    save(destination,result)
    return row['root'],result['seconds']


def scan(rows,calibration,output,workers=4,pre_values=PRES,post_values=POSTS):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    for row in rows:
        previous=output/(Path(row['root']).stem+'.json')
        if previous.exists():
            stored=json.loads(previous.read_text())
            assert stored['fitter_sha256']==digest(PROJECT/'analysis/alpha_gamma_short_window.py')
            assert stored['study_sha256']==digest(Path(__file__))
            assert stored['calibration_sha256']==digest(Path(calibration)/'calibration-arrays.npz')
    missing=[row for row in rows if not (output/(Path(row['root']).stem+'.json')).exists()]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        tasks=[pool.submit(_worker,(row,str(calibration),str(output),pre_values,post_values)) for row in missing]
        for task in as_completed(tasks):
            root,seconds=task.result()
            print('SCANNED',root,round(seconds,2),'s',flush=True)


def load_results(folder):
    return [json.loads(p.read_text()) for p in sorted(Path(folder).glob('*.json'))]


def wilson(k,n):
    if not n:return None
    z=1.959963984540054;p=k/n;denom=1+z*z/n
    center=(p+z*z/(2*n))/denom
    radius=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/denom
    return [center-radius,center+radius]


def summarize(results):
    groups=[]
    for pre in PRES:
        for post in POSTS:
            for method in METHODS:
                for particle in ('alpha','gamma'):
                    events=[r for r in results if r['truth_for_evaluation_only']['particle']==particle]
                    for band,low,high in [('all_above_500',500.,np.inf),('500_1000',500.,1000.),('1000_3000',1000.,3000.),('above_3000',3000.,np.inf)]:
                        selected=[r for r in events if low<=r['truth_for_evaluation_only']['charge_pe']<high]
                        errors=[];types=0;success=0;joint=0;detected=0;bounds=0;dt=[]
                        for event in selected:
                            fit=next(f for f in event['prefixes'] if f['pre_ns']==pre and f['post_ns']==post)
                            if not fit['evaluation']['matched_signal']:continue
                            detected+=1;r=fit['methods'][method];ev=r['evaluation'];error=ev['relative_charge_error']
                            errors.append(error);dt.append(ev['onset_error_ns']);types+=ev['correct_type']
                            success+=abs(error)<=.05;joint+=abs(error)<=.05 and ev['correct_type']
                            bounds+=r['onset_at_search_bound']
                        groups.append({'pre_ns':pre,'post_ns':post,'method':method,'particle':particle,'band':band,
                            'targets':len(selected),'detected':detected,'charge_5pct_success':success,
                            'charge_5pct_fraction':success/len(selected) if selected else None,
                            'charge_5pct_wilson95':wilson(success,len(selected)),
                            'correct_type':types if method.startswith('template') else None,
                            'joint_charge_type_success':joint if method.startswith('template') else None,
                            'onset_bound_count':bounds,
                            'charge_bias':float(np.mean(errors)) if errors else None,
                            'charge_error_sd':float(np.std(errors,ddof=1)) if len(errors)>1 else None,
                            'charge_rmse':float(np.sqrt(np.mean(np.square(errors)))) if errors else None,
                            'abs_charge_error_p90':float(np.quantile(abs(np.array(errors)),.90)) if errors else None,
                            'time_rmse_ns':float(np.sqrt(np.mean(np.square(dt)))) if dt else None})
    return {'groups':groups,'independent_events':len(results),
            'truth_counts':{name:{'injected':sum(r['truth_for_evaluation_only']['particle']==name for r in results),
                'above_500':sum(r['truth_for_evaluation_only']['particle']==name and r['truth_for_evaluation_only']['charge_pe']>=500 for r in results),
                'zero_charge':sum(r['truth_for_evaluation_only']['particle']==name and r['truth_for_evaluation_only']['charge_pe']==0 for r in results)} for name in ('alpha','gamma')}}


def freeze(output):
    output=Path(output);results=load_results(output/'development')
    summary=summarize(results);save(output/'development-summary.json',summary)
    multipliers={};residual_limits={}
    for particle in ('alpha','gamma'):
        ratios=[];scores=[]
        for event in results:
            truth=event['truth_for_evaluation_only']
            if truth['particle']!=particle or truth['charge_pe']<500:continue
            event_ratios=[];event_scores=[]
            for fit in event['prefixes']:
                if fit['pre_ns']!=50. or not fit['evaluation']['matched_signal']:continue
                r=fit['methods']['template_gls'];sigma=r['local_charge_sigma_pe']
                if sigma is not None and sigma>0:
                    event_ratios.append(abs(r['charge_pe']-truth['charge_pe'])/sigma)
                    event_scores.append(r['objective']/max(1,fit['sample_count']-3))
            if event_ratios:ratios.append(max(event_ratios));scores.append(max(event_scores))
        multipliers[particle]=max(1.96,float(np.quantile(ratios,.95))) if ratios else None
        residual_limits[particle]=float(np.quantile(scores,.95)) if scores else None
    assert all(v is not None for v in multipliers.values())
    code_files=['analysis/alpha_gamma_short_window.py','analysis/alpha_gamma_short_window_study.py',
                'analysis/alpha_gamma_waveform_fit.py','analysis/test_alpha_gamma_short_window.py']
    settings={'created':datetime.now().astimezone().isoformat(),'pre_ns':50.,'post_ns':POSTS,
        'primary_method':'template_gls','adaptive_relative_error_tolerance':.05,
        'adaptive_multiplier':max(multipliers.values()),'per_class_development_multiplier':multipliers,
        'adaptive_residual_limit':max(residual_limits.values()),
        'guard_status':'Empirical development guard; no calibrated confidence/coverage guarantee',
        'development_summary_sha256':digest(output/'development-summary.json'),
        'calibration_sha256':digest(output/'calibration/calibration.json'),
        'calibration_arrays_sha256':digest(output/'calibration/calibration-arrays.npz'),
        'source_sha256':{p:digest(PROJECT/p) for p in code_files}}
    save(output/'frozen-settings.json',settings)
    print('FROZEN',settings['adaptive_multiplier'],settings['adaptive_residual_limit'],flush=True)
    return settings


def adaptive(results,settings):
    rows=[]
    for event in results:
        truth=event['truth_for_evaluation_only'];chosen=None
        for post in settings['post_ns']:
            fit=next(f for f in event['prefixes'] if f['pre_ns']==settings['pre_ns'] and f['post_ns']==post)
            if not fit['detected']:continue
            r=fit['methods'][settings['primary_method']];sigma=r['local_charge_sigma_pe']
            if sigma is None or r['charge_pe']<=0 or r['onset_at_search_bound']:continue
            bound=settings['adaptive_multiplier']*sigma/r['charge_pe']
            if bound<=.05 and r['objective']/max(1,fit['sample_count']-3)<=settings['adaptive_residual_limit']:
                chosen={'post_ns':post,'charge_pe':r['charge_pe'],'particle':r['particle'],
                        'predicted_relative_guard':bound,'matched':fit['evaluation']['matched_signal']}
                break
        rows.append({'root':event['row']['root'],'truth_for_evaluation_only':truth,'chosen':chosen})
    counts={}
    for particle in ('alpha','gamma'):
        events=[r for r in rows if r['truth_for_evaluation_only']['particle']==particle and r['truth_for_evaluation_only']['charge_pe']>=500]
        chosen=[r for r in events if r['chosen'] is not None]
        ok=[r for r in chosen if r['chosen']['matched'] and abs(r['chosen']['charge_pe']/r['truth_for_evaluation_only']['charge_pe']-1)<=.05]
        counts[particle]={'targets':len(events),'accepted':len(chosen),'charge_success':len(ok),
            'charge_success_fraction':len(ok)/len(events) if events else None,
            'conditional_charge_success':len(ok)/len(chosen) if chosen else None,
            'joint_charge_type_success':sum(r['chosen']['particle']==particle for r in ok),
            'median_post_ns':float(np.median([r['chosen']['post_ns'] for r in chosen])) if chosen else None,
            'success_wilson95':wilson(len(ok),len(events))}
    return {'rows':rows,'counts':counts,'scope':'Empirical development-calibrated stopping; abstentions remain failures in unconditional success'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=('prepare','development','freeze','summarize','test'))
    p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=4)
    p.add_argument('--test-records',type=Path)
    a=p.parse_args();output=a.output
    if a.stage=='prepare':prepare(output)
    elif a.stage=='development':
        rows=json.loads((output/'input-records.json').read_text())['development']
        scan(rows,output/'calibration',output/'development',a.workers)
    elif a.stage=='freeze':freeze(output)
    elif a.stage=='test':
        settings=json.loads((output/'frozen-settings.json').read_text())
        for path,expected in settings['source_sha256'].items():assert digest(PROJECT/path)==expected
        assert digest(output/'calibration/calibration.json')==settings['calibration_sha256']
        assert digest(output/'calibration/calibration-arrays.npz')==settings['calibration_arrays_sha256']
        rows=json.loads(a.test_records.read_text());scan(rows,output/'calibration',output/'test',a.workers)
    else:
        result=load_results(output/'test');settings=json.loads((output/'frozen-settings.json').read_text())
        save(output/'test-summary.json',summarize(result));save(output/'adaptive-test.json',adaptive(result,settings))


if __name__=='__main__':main()
