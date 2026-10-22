"""Development/freeze/test protocol for correlated fitting, run at project root."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from dataclasses import asdict,replace
from datetime import datetime
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess
import time

import numpy as np
from scipy.optimize import nnls

from .alpha_correlated_fit import CorrelatedConfig,ResidualCalibration,InnovationsMetric,optimize,fit_correlated
from .alpha_full_signal_study import make_template,summarize
from .alpha_local_fit import LocalConfig,fit_local,optimize_local
from .alpha_pileup_fit import FitConfig,fit_waveform,score
from .alpha_rate_campaign import digest,save

PROJECT=Path(__file__).resolve().parents[2]
OLD=PROJECT/'jupyter/artifacts/alpha_local_model_selection_20261020'
EXPERIMENT=PROJECT/'pileup/experiments/alpha_correlated_global_20261020'
OUTPUT=PROJECT/'jupyter/artifacts/alpha_correlated_global_20261020'
RECORD=PROJECT/'jupyter/run_records/alpha_correlated_global_20261020'
METHODS=('diagonal','correlated','global')


def settings():return json.loads((OUTPUT/'baseline-settings.json').read_text())


def prepare():
    assert not (OUTPUT/'protocol.json').exists(), 'A prepared experiment must be preserved'
    for p in (OUTPUT,EXPERIMENT,RECORD):p.mkdir(parents=True,exist_ok=True)
    for sub in ('pilots','calibration','development','test','inputs','controls'): (OUTPUT/sub).mkdir(exist_ok=True)
    for name in ('jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb','jupyter/Run_Catalog.ipynb','docs/literature-use.md','README.md','pileup/README.md'):
        backup=RECORD/('before__'+name.replace('/','__')+'.before')
        if backup.exists():assert digest(backup)==digest(PROJECT/name)
        else:shutil.copyfile(PROJECT/name,backup)
    shutil.copyfile(OLD/'baseline-settings.json',OUTPUT/'baseline-settings.json')
    old=json.loads((OLD/'input-records.json').read_text())
    # ROOT-level calibration and development identities fixed before fitting.
    save(OUTPUT/'input-records.json',dict(calibration_single=old['calibration_single'],
        development_single=old['development_single'],development_pileup=old['development_pileup'],
        published_previous_test='All 23 previous local-fit test windows are exposed; not used as new test.'))
    save(OUTPUT/'protocol.json',dict(created=datetime.now().astimezone().isoformat(),
        hypotheses=['Full correlated covariance vs its diagonal, identical observed-voltage pilot',
            'Budgeted local differential evolution plus SLSQP vs SLSQP; K/time/Q remain unknown',
            'Light GRU vs previous frozen CNN and physical fits, independent ROOT splits'],
        covariance='sigma_e^2 I + diag(s) AR1(rho) diag(s); s^2=v_dark+eta max(pilot,0); frozen per record',
        covariance_estimation='Empirical residual moments in known-single calibration; includes template mismatch, not microscopic truth covariance',
        penalty_multipliers=[.5,1.,2.],
        selection='Minimize extra predictions in native development singles; then maximize macro complete-recovery F1 at 1/10/40 MHz; higher multiplier wins ties',
        optimizer_budget=asdict(CorrelatedConfig(global_calls=12)),
        test_plan={'100000':8,'1000000':3,'10000000':2,'20000000':2,'40000000':2,'60000000':2,'80000000':2,'100000000':2},
        software_controls='Empty/single/double matched-template white noise and independently sampled correlated noise; not native double-pulse resolution',
        adequate_dynamic_range=True,double_scintillation_decay_ns=[95.,351.],
        scoring=dict(charge_threshold_pe=500,time_tolerance_ns=10,relative_charge_tolerance=.2),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        environment=dict(python=platform.python_version(),platform=platform.platform(),
            packages={p:importlib.metadata.version(p) for p in ('numpy','scipy','uproot','awkward','matplotlib')})))
    (RECORD/'git-status-before.txt').write_text(subprocess.check_output(['git','-c','core.fileMode=false','status','--short','--untracked-files=no'],cwd=PROJECT,text=True))
    (RECORD/'working-changes.patch').write_bytes(subprocess.check_output(['git','-c','core.fileMode=false','diff','--binary'],cwd=PROJECT))
    print('PREPARED',flush=True)


def load_row(row,single=False):
    stem=Path(row['root']).stem
    if single:
        waveform=PROJECT/row['waveform'];truth=PROJECT/row['truth'];gate=[0.,6000.];shift=3000.
    else:
        source=PROJECT/row.get('source',str((OUTPUT/'inputs').relative_to(PROJECT)))
        waveform=source/'waveforms'/f'{stem}.npz';truth=source/'truth'/f'{stem}.npz';gate=[5000.,15000.];shift=0.
    with np.load(waveform) as z:t,y=z['time_ns'].copy(),z['analog_mV'].copy()
    if single:
        select=(t>=3000)&(t<9000);t,y=t[select]-3000,y[select]
    return t,y,truth,gate,shift,waveform


def pilot(row,single=False):
    t,y,truth,gate,shift,waveform=load_row(row,single)
    stem=Path(row['root']).stem+('_single' if single else '')
    path=OUTPUT/'pilots'/(stem+'.json');array_path=path.with_suffix('.npz')
    if path.exists():
        cached=json.loads(path.read_text());assert cached['voltage_sha256']==digest(waveform)
        initial=cached['fit'];initial['model_mV']=np.load(array_path)['model_mV'].copy()
        return t,y,truth,gate,shift,waveform,initial,cached['seconds']
    old=OLD/'baseline-cache'/(stem+'.json')
    if old.exists():
        cached=json.loads(old.read_text());initial=cached['fit'];elapsed=cached['seconds']
        initial['model_mV']=np.load(old.with_suffix('.npz'))['model_mV'].copy()
    else:
        d=settings();tick=time.monotonic();initial=fit_waveform(t,y,make_template(d),FitConfig(**d['fit_config']));elapsed=time.monotonic()-tick
    assert initial['model_mV'].shape==y.shape
    template=make_template(settings())
    expected=10+template(t[:,None]-np.asarray(initial['times_ns'])[None,:])@np.asarray(initial.get('peak_mV',[]))
    np.testing.assert_allclose(initial['model_mV'],expected,atol=1e-7)
    np.savez_compressed(array_path,model_mV=initial['model_mV'])
    save(path,dict(fit={k:v for k,v in initial.items() if k!='model_mV'},seconds=elapsed,
        voltage_sha256=digest(waveform),baseline_settings_sha256=digest(OUTPUT/'baseline-settings.json'),
        cache_origin=str(old.relative_to(PROJECT)) if old.exists() else 'fresh observed-voltage pilot'))
    return t,y,truth,gate,shift,waveform,initial,elapsed


def calibrate():
    data=json.loads((OUTPUT/'input-records.json').read_text());d=settings();template=make_template(d);base=FitConfig(**d['fit_config'])
    samples=[];diagnostics=[]
    for row in data['calibration_single']:
        t,y,_,_,_,_,initial,_=pilot(row,True)
        if not initial['times_ns']:continue
        center=initial['times_ns'][int(np.argmax(initial['charge_pe']))];mask=(t>=center-24-4)&(t<=center+24+1400)
        one,_=optimize_local(t[mask],y[mask]-10,[[center]],template,max(0,center-24),min(t[-1],center+24),LocalConfig())
        mu=template(t[:,None]-one['times'])@one['amplitude'];residual=y-10-mu
        samples.append((t,mu,residual,one))
        diagnostics.append(dict(root=row['root'],fitted_charge_pe=float(one['amplitude'].sum()/(base.area_per_pe_mV_ns*template.unit_area_peak))))
    sigma=.1
    # Equal waveform lengths and one row per single event; no test labels.
    mu=np.concatenate([s[1] for s in samples]);residual=np.concatenate([s[2] for s in samples])
    design=np.c_[np.ones(len(mu)),np.maximum(mu,0.)];weight=1/np.sqrt(np.maximum(mu,1.))
    variance=nnls(design*weight[:,None],(residual**2-sigma**2)*weight)[0]
    v0,eta=map(float,variance);eta=max(eta,1e-8)
    numerator=denominator=0.;empirical=[]
    for t,m,r,one in samples:
        s=np.sqrt(v0+eta*np.maximum(m,0.))
        numerator+=float(r[1:]@r[:-1]);denominator+=float(s[1:]@s[:-1])
        empirical.append(dict(mean_mV=m.tolist(),residual_mV=r.tolist()))
    rho=float(np.clip(numerator/max(denominator,1e-20),0.,.995))
    cal=ResidualCalibration(sigma,v0,eta,rho)
    save(OUTPUT/'calibration/covariance.json',dict(parameters=asdict(cal),single_windows=len(samples),
        rows=diagnostics,fit_inputs='Known-single identity, measured voltage; no charge/time truth used',
        covariance_scope='Empirical calibrated residual model; AR1 structure assumed, shared-cell covariance not proven'))
    np.savez_compressed(OUTPUT/'calibration/residuals.npz',mean=np.stack([s[1] for s in samples]),residual=np.stack([s[2] for s in samples]))
    for method in METHODS:
        cfg=CorrelatedConfig(diagonal=method=='diagonal',global_calls=12 if method=='global' else 0)
        gains=[]
        for i,(t,mu,r,one) in enumerate(samples):
            center=float(one['times'][0]);lo,hi=max(0,center-24),min(t[-1],center+24)
            mask=(t>=lo-4)&(t<=hi+cfg.tail_ns);values=mu+r
            metric=InnovationsMetric(mu,cal,diagonal=cfg.diagonal);local=metric.subset(int(np.flatnonzero(mask)[0]),int(np.flatnonzero(mask)[-1])+1)
            first,_=optimize(t[mask],values[mask],[[center]],template,lo,hi,local,cfg)
            seeds=[]
            for gap in cfg.split_gaps_ns:seeds.extend(([center-gap/2,center+gap/2],[center,center+gap],[center-gap,center]))
            second,_=optimize(t[mask],values[mask],seeds,template,lo,hi,local,cfg,
                global_rng=np.random.default_rng(cfg.global_seed+i) if cfg.global_calls else None)
            gain=0.
            if second is not None and np.min(second['amplitude']/(base.area_per_pe_mV_ns*template.unit_area_peak))>=cfg.component_charge_pe:
                r1=metric.apply(values-template(t[:,None]-first['times'])@first['amplitude'])
                r2=metric.apply(values-template(t[:,None]-second['times'])@second['amplitude'])
                gain=max(0.,float(r1@r1-r2@r2))
            q=max(500.,float(first['amplitude'].sum()/(base.area_per_pe_mV_ns*template.unit_area_peak)))
            gains.append(dict(index=i,gain=gain,fitted_charge_pe=q,gain_per_pe=gain/q))
        save(OUTPUT/f'calibration/penalty-{method}.json',dict(method=method,rows=gains,
            coefficient=max([g['gain_per_pe'] for g in gains]+[1e-8]),
            rule='Maximum calibration-single full-record split gain / fitted charge; empirical, not a p-value'))
        print('CALIBRATION',method,'cov',asdict(cal),'penalty',max(g['gain_per_pe'] for g in gains),flush=True)


def fit_row(row,config,*,single=False,phase='development',tag='correlated'):
    folder=OUTPUT/phase/tag;folder.mkdir(exist_ok=True,parents=True);stem=Path(row['root']).stem
    path=folder/(stem+'.json')
    if path.exists():return json.loads(path.read_text())
    t,y,truth_path,gate,shift,waveform,initial,baseline_seconds=pilot(row,single)
    d=settings();template=make_template(d);base=FitConfig(**d['fit_config'])
    cal=ResidualCalibration(**json.loads((OUTPUT/'calibration/covariance.json').read_text())['parameters'])
    tick=time.monotonic();fitted=fit_correlated(t,y,template,base,cal,config,initial=initial);refine_seconds=time.monotonic()-tick
    # Truth is first opened only after inference above.
    with np.load(truth_path) as z:truth={k:z[k].copy() for k in ('first_active_ns','emission_ns','charge_pe')}
    truth['first_active_ns']-=shift;truth['emission_ns']-=shift
    np.savez_compressed(folder/(stem+'.npz'),model_mV=fitted.pop('model_mV'))
    result=dict(root=row['root'],seed=row.get('seed'),rate_hz=row.get('rate_hz',0),single=single,
        fit=fitted,score=score(fitted,truth,gate),seconds=baseline_seconds+refine_seconds,
        refinement_seconds=refine_seconds,waveform=str(waveform.relative_to(PROJECT)),truth=str(truth_path.relative_to(PROJECT)),
        input_sha256=digest(waveform),config_sha256=digest(OUTPUT/'baseline-settings.json'),
        inference_truth_inputs=[],gate_ns=gate,truth_shift_ns=shift)
    save(path,result)
    z=result['score'];print('FIT',phase,tag,stem,z['recovered_count'],z['target_count'],z['unmatched_predictions'],round(refine_seconds,2),flush=True)
    return result


def develop(workers=2):
    data=json.loads((OUTPUT/'input-records.json').read_text());selected={};candidates=[]
    for method in METHODS:
        coefficient=json.loads((OUTPUT/f'calibration/penalty-{method}.json').read_text())['coefficient']
        for factor in (.5,1.,2.):
            cfg=CorrelatedConfig(penalty_per_pe=coefficient,penalty_multiplier=factor,
                diagonal=method=='diagonal',global_calls=12 if method=='global' else 0)
            tag=f'{method}_{factor:g}'
            tasks=[(r,True) for r in data['development_single']]+[(r,False) for r in data['development_pileup']]
            with ProcessPoolExecutor(max_workers=workers) as pool:
                fs=[pool.submit(fit_row,row,cfg,single=single,tag=tag) for row,single in tasks]
                rows=[f.result() for f in as_completed(fs)]
            extra=sum(r['score']['unmatched_predictions'] for r in rows if r['single'])
            f1=[2*r['score']['recovered_count']/max(1,r['score']['target_count']+r['score']['reconstructed_count']) for r in rows if not r['single']]
            candidates.append(dict(method=method,factor=factor,tag=tag,config=asdict(cfg),single_extras=extra,macro_f1=float(np.mean(f1))))
        selected[method]=max([c for c in candidates if c['method']==method],key=lambda c:(-c['single_extras'],c['macro_f1'],c['factor']))
        save(OUTPUT/'development-summary.json',dict(candidates=candidates,selected=selected))
        print('SELECTED',method,selected[method],flush=True)


def freeze():
    chosen=json.loads((OUTPUT/'development-summary.json').read_text())['selected'];assert set(chosen)==set(METHODS)
    assert (OUTPUT/'gru/frozen-model.json').exists()
    inference=EXPERIMENT/'inference-source';assert not inference.exists()
    cnn_origin=PROJECT/'jupyter/artifacts/alpha_ml_comparison_20261014_162000'
    cnn_reference=EXPERIMENT/'reference-cnn';cnn_reference.mkdir()
    for name in ('model.pt','frozen-model.json','analysis_config.json'):
        shutil.copyfile(cnn_origin/name,cnn_reference/name)
    assets=EXPERIMENT/'frozen-analysis';assets.mkdir()
    for name in ('baseline-settings.json','protocol.json','input-records.json','development-summary.json'):
        shutil.copyfile(OUTPUT/name,assets/name)
    shutil.copytree(OUTPUT/'calibration',assets/'calibration')
    shutil.copytree(OUTPUT/'gru',assets/'gru')
    (assets/'local-reference').mkdir()
    shutil.copyfile(OLD/'freeze.json',assets/'local-reference/freeze.json')
    originals=list((PROJECT/'pileup/analysis').glob('*.py'))+list((PROJECT/'analysis').glob('*.py'))
    for path in originals:
        destination=inference/path.relative_to(PROJECT);destination.parent.mkdir(exist_ok=True,parents=True);shutil.copyfile(path,destination)
    for package in ('analysis','pileup','pileup/analysis'):(inference/package/'__init__.py').touch(exist_ok=True)
    files=originals+list(inference.rglob('*.py'))+list((EXPERIMENT/'source').rglob('*'))+[EXPERIMENT/'bin/gagg']
    files+=[OUTPUT/'baseline-settings.json',OUTPUT/'protocol.json',OUTPUT/'input-records.json',OUTPUT/'development-summary.json']
    files+=list((OUTPUT/'calibration').glob('*.json'))+list((OUTPUT/'gru').glob('*.json'))+[OUTPUT/'gru/model.pt']
    files+=list(cnn_reference.glob('*'))
    files+=list(assets.rglob('*'))+[EXPERIMENT/'replay_fit.py']
    files=[p for p in files if p.is_file()]
    frozen=dict(frozen_at=datetime.now().astimezone().isoformat(),selected=chosen,
        frozen_files={str(p.relative_to(PROJECT)):digest(p) for p in files})
    save(OUTPUT/'freeze.json',frozen);save(EXPERIMENT/'frozen-fit-settings.json',frozen)
    print('FROZEN',frozen['frozen_at'],flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','calibrate','develop','freeze']);parser.add_argument('--workers',type=int,default=2)
    a=parser.parse_args();assert 1<=a.workers<=4
    if a.stage=='develop':develop(a.workers)
    else:globals()[a.stage]()


if __name__=='__main__':main()
