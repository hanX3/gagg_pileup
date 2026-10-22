"""Calibrate isolated native responses and compare blind fits on separate data."""
from __future__ import annotations
import argparse
from dataclasses import asdict
from datetime import datetime
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import least_squares

from .alpha_pileup_fit import AlphaTemplate
from .alpha_gamma_pileup_fit import fit_waveform, score_mixed
from .alpha_gamma_study import calibration as baseline_calibration
from .alpha_gamma_regularized_fit import CalibratedBank, FrozenWhitening, RegularizedConfig, fit_regularized
from . import alpha_gamma_pileup_fit as baseline_engine
from .alpha_gamma_nnls import reduced_nnls
from .alpha_rate_campaign import digest, save


def calibrate(folder,inputs,output):
    folder,inputs,output=map(Path,(folder,inputs,output))
    if output.exists():raise FileExistsError(output)
    manifest=json.loads((folder/'manifest.json').read_text())
    rows=[r for r in manifest['runs'] if r['role']=='calibration']
    assert len(rows)==len([r for r in manifest['planned_runs'] if r['role']=='calibration'])
    assert all(r['status']=='completed' for r in rows)
    old_path=folder/'baseline-settings.json'
    if not old_path.exists():old_path=Path('jupyter/run_records/alpha_gamma_fit_20261015_093500/calibration-settings.json')
    old=json.loads(old_path.read_text())
    readout=old['templates']['alpha']['readout'];area=old['fit_config']['area_per_pe_mV_ns']
    basis=[AlphaTemplate(readout,f) for f in (0.,1.)]
    samples=[];pre=[];individual=[]
    for row in rows:
        stem=Path(row['root']).stem
        wave=dict(np.load(inputs/'waveforms'/f'{stem}.npz'))
        truth=dict(np.load(inputs/'truth'/f'{stem}.npz'))
        assert len(truth['charge_pe'])==1
        t,y=wave['time_ns'],wave['analog_mV'];pre.append(y[(t>500)&(t<4000)])
        q=float(truth['charge_pe'][0]);onset=float(truth['first_active_ns'][0])
        item=dict(root=row['root'],seed=row['seed'],particle=row['particle'],source_energy_MeV=row['energy_MeV'],
            charge_pe=q,first_active_ns=onset if np.isfinite(onset) else None,
            role='independent isolated calibration; truth used only to calibrate',used=bool(q>=50 and np.isfinite(onset)))
        individual.append(item)
        if not item['used']:continue
        keep=(t>=onset-20)&(t<=onset+3000)
        tx=t[keep]-onset;values=y[keep]-10.
        def model(z):
            f,delay=z
            return q*area*((1-f)*basis[0](tx-delay)*basis[0].unit_area_peak+
                           f*basis[1](tx-delay)*basis[1].unit_area_peak)
        opt=least_squares(lambda z:model(z)-values,[.7 if row['particle']=='alpha' else .21,0.],
            bounds=([.02,-10.],[.98,10.]),x_scale=[.2,2.],max_nfev=100)
        item.update(slow_fraction=float(opt.x[0]),delay_ns=float(opt.x[1]),solver_success=bool(opt.success))
        samples.append((item,t,y,q,onset))
    classes={}
    for particle in ('alpha','gamma'):
        group=[s[0] for s in samples if s[0]['particle']==particle]
        assert len(group)>=6,(particle,len(group))
        q=np.array([s['charge_pe'] for s in group]);fractions=np.array([s['slow_fraction'] for s in group])
        X=np.c_[np.ones(len(q)),np.log(q/1000.)]
        # Shot fluctuations shrink with charge. Cap weights so a single bright
        # pulse cannot determine the entire energy dependence.
        w=np.sqrt(np.minimum(q,4000.))
        coeff=np.linalg.lstsq(X*w[:,None],fractions*w,rcond=None)[0]
        classes[particle]=dict(slow_intercept=float(coeff[0]),slow_log_slope=float(np.clip(coeff[1],-.08,.08)),
            delay_ns=float(np.average([s['delay_ns'] for s in group],weights=w*w)),
            charge_range_pe=[float(q.min()),float(q.max())],used_samples=len(group),
            rms_fraction_residual=float(np.sqrt(np.mean((fractions-X@coeff)**2))))
    # Electronics/dark baseline from pre-source data; robust differencing avoids
    # infrequent dark bursts determining the white electronics noise scale.
    pre=np.concatenate(pre);sigma=np.median(abs(np.diff(pre)-np.median(np.diff(pre))))/.6744897501960817/np.sqrt(2)
    cal=dict(created=datetime.now().astimezone().isoformat(),readout=readout,area_per_pe_mV_ns=area,
        baseline_mV=10.,classes=classes,source_roles=['calibration'],isolated_rows=individual,
        model='Calibrated log-charge-dependent slow fraction, fixed 95/351 ns basis and mean onset delay',
        physics_changed=False,extrapolation='Clamp shape charge to calibration range; charge estimate itself not clamped',
        baseline_policy='Known 10 mV electronics setting; variance estimated from independent pre-source samples',
        limitations=['Central single-particle calibration only; no position/occupancy calibration',
                    'AR(1) and variance proportional to observed signal are approximations',
                    'Charge means realized full-cell equivalent pe, not deposited energy'])
    bank=CalibratedBank(cal);rs=[];mus=[]
    for item,t,y,q,onset in samples:
        k=0 if item['particle']=='alpha' else 1
        mu=q*bank.column(t,onset,k,q);mask=(t>=onset)&(t<onset+2000)
        rs.append((y-10.-mu)[mask]);mus.append(mu[mask])
    r=np.concatenate(rs);mu=np.concatenate(mus)
    beta=max(1e-6,float(np.sum(r*r-sigma*sigma)/np.sum(mu)))
    z=[a/np.sqrt(sigma*sigma+beta*b) for a,b in zip(rs,mus)]
    numerator=sum(float(a[1:]@a[:-1]) for a in z)
    denominator=sum(float(a[:-1]@a[:-1]) for a in z)
    rho=float(np.clip(numerator/denominator,0.,.995))
    cal['noise']=dict(baseline_variance_mV2=float(sigma*sigma),variance_per_mV=beta,
        shot_rho_at_2ns=rho,estimated_from='Population-template residuals of calibration events, 0–2000 ns after first active photon',
        covariance='Heteroscedastic AR(1) factor; fixed weights from smoothed observed voltage per fit')
    cal['manifest_sha256']=digest(folder/'manifest.json')
    cal['waveform_sha256']={Path(r['root']).stem:digest(inputs/'waveforms'/f"{Path(r['root']).stem}.npz") for r in rows}
    cal['truth_sha256']={Path(r['root']).stem:digest(inputs/'truth'/f"{Path(r['root']).stem}.npz") for r in rows}
    save(output,cal)
    print('CALIBRATION',classes,cal['noise'],flush=True)
    return cal


def evaluate(folder,inputs,calibration_file,output,role,penalties=(3.,8.),baseline=True):
    folder,inputs,calibration_file,output=map(Path,(folder,inputs,calibration_file,output))
    output.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((folder/'manifest.json').read_text())
    rows=[r for r in manifest['runs'] if r['role']==role]
    assert rows and len(rows)==len([r for r in manifest['planned_runs'] if r['role']==role])
    cal=json.loads(calibration_file.read_text())
    # Same algebraically equivalent NNLS reduction used by the previous formal
    # baseline campaign; no regularization or changed baseline objective.
    baseline_engine.nnls=reduced_nnls
    old_path=folder/'baseline-settings.json'
    if not old_path.exists():old_path=Path('jupyter/run_records/alpha_gamma_fit_20261015_093500/calibration-settings.json')
    _,oldbank,oldcfg=baseline_calibration(old_path)
    source=Path(__file__).with_name('alpha_gamma_regularized_fit.py');code_hash=digest(source)
    if role=='test_mixed':
        frozen=json.loads((folder/'frozen-fit-settings.json').read_text())
        assert code_hash==frozen['fitter_sha256'] and digest(calibration_file)==frozen['calibration_sha256']
        penalties=(frozen['config']['selection_penalty'],)
    for row in rows:
        stem=Path(row['root']).stem;wave=dict(np.load(inputs/'waveforms'/f'{stem}.npz'))
        for name,penalty in ([('baseline',None)] if baseline else [])+[(f'regularized_{p:g}',p) for p in penalties]:
            destination=output/name/f'{stem}.json'
            if destination.exists():
                stored=json.loads(destination.read_text())
                assert stored['fitter_sha256']==code_hash and stored['calibration_sha256']==digest(calibration_file)
                continue
            start=time.perf_counter()
            if penalty is None:fit=fit_waveform(wave['time_ns'],wave['analog_mV'],oldbank,oldcfg)
            else:fit=fit_regularized(wave['time_ns'],wave['analog_mV'],cal,RegularizedConfig(selection_penalty=penalty))
            seconds=time.perf_counter()-start
            # Labels are opened after reconstruction; they do not enter inference.
            truth=dict(np.load(inputs/'truth'/f'{stem}.npz'))
            meta=json.loads((inputs/'metadata'/f'{stem}.json').read_text())
            scored=score_mixed(fit,truth,meta['gate_ns'])
            model=fit.pop('model_mV');destination.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(destination.with_suffix('.npz'),time_ns=wave['time_ns'],model_mV=model)
            result=dict(row=row,method=name,fit=fit,score=scored,fit_seconds=seconds,
                fitter_sha256=code_hash,baseline_fitter_sha256=digest(Path(__file__).with_name('alpha_gamma_pileup_fit.py')),
                baseline_linear_solver='Checked equivalent Gram/Cholesky reduction with direct NNLS fallback',
                nnls_helper_sha256=digest(Path(__file__).with_name('alpha_gamma_nnls.py')),
                calibration_sha256=digest(calibration_file),waveform_sha256=digest(inputs/'waveforms'/f'{stem}.npz'),
                truth_sha256=digest(inputs/'truth'/f'{stem}.npz'),fit_truth_inputs=[],gate_ns=meta['gate_ns'])
            save(destination,result)
            print('FIT',role,row['label'],name,{k:(v['recovered_count'],v['target_count']) for k,v in scored['by_particle'].items()},
                'extra',scored['unmatched_predictions'],'time',round(seconds,2),'components',fit.get('selected_components'),flush=True)


def controlled(calibration_file,output,penalty):
    output=Path(output)
    if output.exists():raise FileExistsError(output)
    cal=json.loads(Path(calibration_file).read_text());bank=CalibratedBank(cal)
    t=np.arange(.25,12000.,2.);rows=[]
    for gap in (10.,20.,40.,80.,200.):
        for pair in ((0,0),(0,1),(1,0),(1,1)):
            times=np.array([5000.,5000.+gap]);q=np.array([3000.,1000.])
            seed=261039000+len(rows)
            y=10.+bank.matrix(t,times,pair,q)@q+np.random.default_rng(seed).normal(0.,.1,len(t))
            truth=dict(first_active_ns=times,emission_ns=times,charge_pe=q,particle_name=bank.names[list(pair)])
            for stride in (1,2):
                fit=fit_regularized(t,y,cal,RegularizedConfig(selection_penalty=penalty,fit_stride=stride))
                fit.pop('model_mV')
                rows.append(dict(gap_ns=gap,pair=bank.names[list(pair)].tolist(),seed=seed,fit_stride=stride,
                    truth_charge_pe=q.tolist(),fit=fit,score=score_mixed(fit,truth,[4000.,8000.])))
    record=dict(rows=rows,calibration_sha256=digest(calibration_file),
        fitter_sha256=digest(Path(__file__).with_name('alpha_gamma_regularized_fit.py')),
        purpose='Numerical controls and 2/4 ns fit-grid comparison; not independent physical performance',
        input='Same calibrated template as inversion + 0.1 mV Gaussian electronics noise; no photon/shared-cell fluctuations')
    save(output,record);print('CONTROLS',len(rows),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--folder',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True)
    p.add_argument('--calibration',type=Path,required=True);p.add_argument('--output',type=Path)
    p.add_argument('--stage',choices=('calibrate','controlled','validation_single','validation_mixed','test_mixed'),required=True)
    p.add_argument('--penalties',nargs='+',type=float,default=[3.,8.]);p.add_argument('--no-baseline',action='store_true')
    a=p.parse_args()
    if a.stage=='calibrate':calibrate(a.folder,a.inputs,a.calibration)
    elif a.stage=='controlled':controlled(a.calibration,a.output,a.penalties[0])
    else:evaluate(a.folder,a.inputs,a.calibration,a.output,a.stage,a.penalties,not a.no_baseline)


if __name__=='__main__':main()
