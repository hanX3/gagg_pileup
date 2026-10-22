"""Blind fixed-double-decay fits with frozen correlated residual covariance.

C[i,j] = s[i] s[j] rho**abs(i-j) + sigma_e**2 delta[i,j], where
s**2 = dark_variance + shot_coefficient * max(observed-pilot mean, 0).
The scalar latent AR(1) innovations factor applies the full inverse covariance
without forming an N x N matrix. It is an empirical residual model, not a
microscopic SiPM covariance. C is fixed before any K/time/charge comparisons.
All discovery uses measured voltage; truth and source rate are not arguments.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import time

import numpy as np
from scipy.linalg.lapack import dgtsv
from scipy.optimize import LinearConstraint, differential_evolution, minimize, nnls

from .alpha_local_fit import matrix_and_slope
from .alpha_pileup_fit import FitConfig, candidate_peaks, fit_waveform


@dataclass(frozen=True)
class ResidualCalibration:
    electronic_sigma_mV: float = .1
    dark_variance_mV2: float = 0.
    shot_coefficient_mV: float = 1.
    rho: float = .9


class InnovationsMetric:
    """W with W C W.T = I for the stated finite-record covariance.

    The stationary unit-variance latent state is used at the first sample.
    Bidiagonal solves implement the nonstationary Kalman mean recurrence in
    compiled LAPACK. Identical C/W is used for every proposal in one record.
    """
    def __init__(self, pilot_mean_mV, calibration: ResidualCalibration, *, diagonal=False):
        mu=np.asarray(pilot_mean_mV,float)
        if mu.ndim!=1 or not len(mu) or not np.isfinite(mu).all():
            raise ValueError('Finite one-dimensional pilot required')
        if calibration.electronic_sigma_mV<=0 or calibration.dark_variance_mV2<0 or calibration.shot_coefficient_mV<0 or not 0<=calibration.rho<1:
            raise ValueError('Invalid covariance calibration')
        self.calibration=calibration
        self.pilot=mu.copy()
        self.diagonal=diagonal
        self.s=np.sqrt(calibration.dark_variance_mV2+calibration.shot_coefficient_mV*np.maximum(mu,0))
        self.rho=0. if diagonal else calibration.rho
        n=len(mu);variance=1.;f=np.empty(n);gain=np.empty(n)
        white=calibration.electronic_sigma_mV**2
        for i in range(n):
            f[i]=self.s[i]**2*variance+white
            gain[i]=variance*self.s[i]/f[i]
            posterior=variance*white/f[i]
            variance=self.rho**2*posterior+1-self.rho**2
        self.sqrt_f=np.sqrt(f)
        self.band=np.zeros((2,n));self.band[0]=1.
        self.band[1,:-1]=-self.rho*(1-gain[:-1]*self.s[:-1])
        self.diagonal_entries=np.ones(n);self.upper_entries=np.zeros(n-1)
        self.drive=self.rho*gain[:-1]
        self.logdet=float(np.log(f).sum())

    def apply(self, values):
        x=np.asarray(values,float)
        if x.ndim not in (1,2) or x.shape[0]!=len(self.s):
            raise ValueError('Whitening input has wrong dimensions')
        matrix=x[:,None] if x.ndim==1 else x
        if self.rho==0 or len(self.s)==1:
            answer=matrix/self.sqrt_f[:,None]
            return answer[:,0] if x.ndim==1 else answer
        if matrix.shape[1]==0:return matrix.copy()
        rhs=np.zeros_like(matrix)
        rhs[1:]=self.drive[:,None]*matrix[:-1]
        # This bidiagonal system is a tridiagonal one with a zero upper band.
        # Direct LAPACK avoids the generic band solver and array-API overhead.
        predicted,info=dgtsv(self.band[1,:-1],self.diagonal_entries,self.upper_entries,rhs,overwrite_b=1)[3:]
        if info:raise np.linalg.LinAlgError(f'Innovation recurrence failed: {info}')
        answer=(matrix-self.s[:,None]*predicted)/self.sqrt_f[:,None]
        return answer[:,0] if x.ndim==1 else answer

    def subset(self, first, last):
        # Marginal covariance of the local proposal window. Final acceptance
        # always uses the original full-record metric, including exterior data.
        return InnovationsMetric(self.pilot[first:last],self.calibration,diagonal=self.diagonal)


@dataclass(frozen=True)
class CorrelatedConfig:
    report_charge_pe: float = 500.
    component_charge_pe: float = 100.
    min_separation_ns: float = 3.
    pair_distance_ns: float = 200.
    move_ns: float = 24.
    tail_ns: float = 1400.
    passes: int = 2
    max_iterations: int = 35
    split_gaps_ns: tuple = (6.,16.,40.)
    max_additions: int = 160
    residual_budget: int = 32
    penalty_per_pe: float = .01
    penalty_multiplier: float = 1.
    diagonal: bool = False
    global_calls: int = 0
    global_popsize: int = 6
    global_iterations: int = 12
    global_seed: int = 262025007


def projected(t, target, times, template, metric):
    h,slope=matrix_and_slope(template,t[:,None]-np.asarray(times)[None,:])
    # One compiled whitening solve for the design, its derivatives and data.
    k=len(times);z=metric.apply(np.column_stack((h,slope,target)))
    design,derivative,data=z[:,:k],z[:,k:2*k],z[:,-1]
    amplitude=nnls(design,data,maxiter=max(100,20*k))[0]
    residual=design@amplitude-data
    return float(residual@residual),-2*amplitude*(derivative.T@residual),amplitude


def optimize(t,target,seeds,template,lower,upper,metric,config,*,global_rng=None):
    best=None;trials=[]
    for seed in seeds:
        seed=np.sort(np.asarray(seed,float));k=len(seed)
        if not k or upper-lower<(k-1)*config.min_separation_ns:continue
        seed=np.clip(seed,lower,upper)
        if k>1 and np.min(np.diff(seed))<config.min_separation_ns-1e-7:continue
        origin=float(np.mean(seed));scale=max(projected(t,target,seed,template,metric)[0],1.)
        def objective(z):
            value,gradient,_=projected(t,target,origin+10*z,template,metric)
            return value/scale,gradient*10/scale
        constraints=()
        if k>1:
            difference=np.eye(k,k,1)[:-1]-np.eye(k)[:-1]
            constraints=(LinearConstraint(difference,config.min_separation_ns/10,np.inf),)
        result=minimize(objective,(seed-origin)/10,jac=True,method='SLSQP',
            bounds=[((lower-origin)/10,(upper-origin)/10)]*k,constraints=constraints,
            options=dict(maxiter=config.max_iterations,ftol=1e-9))
        for times in (seed,origin+10*result.x):
            if np.any(times<lower-1e-6) or np.any(times>upper+1e-6):continue
            if k>1 and np.min(np.diff(times))<config.min_separation_ns-1e-5:continue
            cost,_,amplitude=projected(t,target,times,template,metric)
            if best is None or cost<best['cost']:
                best=dict(cost=cost,times=times.copy(),amplitude=amplitude,
                          optimizer_success=bool(result.success))
        trials.append(dict(kind='SLSQP',success=bool(result.success),iterations=int(result.nit),evaluations=int(result.nfev)))
    if global_rng is not None and best is not None and len(best['times'])<=3:
        k=len(best['times']);scale=max(best['cost'],1.)
        def objective_global(z):
            ordered=lower+(upper-lower)*np.sort(z)
            if k>1 and np.min(np.diff(ordered))<config.min_separation_ns:
                return 1e3+float(np.maximum(0,config.min_separation_ns-np.diff(ordered)).sum())
            return projected(t,target,ordered,template,metric)[0]/scale
        # Search on [0,1]: SciPy's internal conversion of large absolute times
        # can otherwise round an exactly-on-bound x0 outside its unit box.
        normalized=np.clip((best['times']-lower)/(upper-lower),1e-12,1-1e-12)
        result=differential_evolution(objective_global,[(0.,1.)]*k,
            maxiter=config.global_iterations,popsize=config.global_popsize,
            rng=global_rng,polish=False,tol=0.,atol=0.,workers=1,updating='immediate',x0=normalized)
        proposal=lower+(upper-lower)*np.sort(result.x)
        if k==1 or np.min(np.diff(proposal))>=config.min_separation_ns:
            refined,extra=optimize(t,target,[proposal],template,lower,upper,metric,config)
            trials+=extra
            if refined is not None and refined['cost']<best['cost']:best=refined
        trials.append(dict(kind='differential_evolution',success=bool(result.success),iterations=int(result.nit),evaluations=int(result.nfev)))
    return best,trials


def fit_correlated(time_ns,voltage_mV,template,baseline: FitConfig,calibration: ResidualCalibration,
                   config=CorrelatedConfig(),*,initial=None):
    t,y=np.asarray(time_ns,float),np.asarray(voltage_mV,float)
    if t.ndim!=1 or y.shape!=t.shape or len(t)<10 or not np.isfinite(y).all() or not np.isfinite(t).all():
        raise ValueError('Invalid measured waveform')
    if not np.allclose(np.diff(t),t[1]-t[0],atol=1e-8,rtol=0) or abs(t[0])>t[1]-t[0]:
        raise ValueError('Full uniformly sampled record required')
    if config.penalty_multiplier<=0 or config.min_separation_ns<=0:raise ValueError('Invalid fit configuration')
    tick=time.monotonic()
    if initial is None:initial=fit_waveform(t,y,template,baseline)
    if initial['status']!='fitted':return {**initial,'correlated_status':'initial_fit_failed'}
    target=y-baseline.baseline_mV
    # Only the unchanged observed-voltage pilot sets C. Do not update it with
    # a candidate's residual/charge: K comparisons must use identical weights.
    metric=InnovationsMetric(np.asarray(initial['model_mV'])-baseline.baseline_mV,calibration,diagonal=config.diagonal)
    times=np.asarray(initial['times_ns'],float).copy();amplitudes=np.asarray(initial['peak_mV'],float).copy()
    h=template(t[:,None]-times);wh=metric.apply(h);wy=metric.apply(target)
    if len(times):amplitudes=nnls(wh,wy,maxiter=max(100,20*len(times)))[0]
    residual=target-h@amplitudes;wr=wy-wh@amplitudes
    peak_per_pe=baseline.area_per_pe_mV_ns*template.unit_area_peak
    changes=[];diagnostics=[];trials=[];additions=0;global_remaining=config.global_calls
    rng=np.random.default_rng(config.global_seed);global_bins=set()

    def propose(ix,extra=(),seeds=None):
        nonlocal global_remaining
        anchors=np.r_[times[ix],extra]
        lower=max(t[0],float(anchors.min())-config.move_ns)
        upper=min(t[-1],float(anchors.max())+config.move_ns)
        outside=np.delete(times,ix);earlier=outside[outside<anchors.min()];later=outside[outside>anchors.max()]
        if len(earlier):lower=max(lower,earlier.max()+config.min_separation_ns)
        if len(later):upper=min(upper,later.min()-config.min_separation_ns)
        first=int(np.searchsorted(t,lower-4));last=int(np.searchsorted(t,upper+config.tail_ns,side='right'))
        values=(residual+h[:,ix]@amplitudes[ix])[first:last]
        # Distribute the finite budget across the acquired record; spending
        # every call on the first few guard-region pairs would miss its centre.
        slot=min(config.global_calls-1,int((float(np.mean(anchors))-t[0])/max(t[-1]-t[0],1)*config.global_calls)) if config.global_calls else -1
        use_global=global_remaining>0 and slot not in global_bins and seeds is not None and 1<len(seeds[0])<=3
        if use_global:
            global_remaining-=1;global_bins.add(slot)
        answer,attempts=optimize(t[first:last],values,seeds or [times[ix]],template,lower,upper,
            metric.subset(first,last),config,global_rng=rng if use_global else None)
        trials.extend(attempts)
        return answer

    def update(ix,proposal,kind,reference_charge=0.):
        nonlocal times,amplitudes,h,wh,residual,wr,additions
        if proposal is None:return False
        tt,aa=proposal['times'],proposal['amplitude']
        if len(aa) and np.min(aa/peak_per_pe)<config.component_charge_pe:return False
        new_h=template(t[:,None]-tt);new_wh=metric.apply(new_h)
        new_wr=wr+wh[:,ix]@amplitudes[ix]-new_wh@aa
        gain=float(wr@wr-new_wr@new_wr);delta_k=len(tt)-len(ix)
        qref=max(config.report_charge_pe,reference_charge,float(amplitudes[ix].sum()/peak_per_pe))
        penalty=config.penalty_multiplier*(1.+config.penalty_per_pe*qref)
        entry=dict(kind=kind,gain=gain,delta_k=delta_k,penalty=penalty,accepted=False)
        diagnostics.append(entry)
        if gain<=penalty*delta_k+1e-7:return False
        keep=np.ones(len(times),bool);keep[ix]=False
        if len(tt) and len(times[keep]) and np.min(abs(times[keep,None]-tt))<config.min_separation_ns:return False
        entry['accepted']=True
        old_times=times[ix].copy()
        if delta_k==0:
            times[ix]=tt;amplitudes[ix]=aa;h[:,ix]=new_h;wh[:,ix]=new_wh
        else:
            times=np.r_[times[keep],tt];amplitudes=np.r_[amplitudes[keep],aa]
            h=np.column_stack((h[:,keep],new_h));wh=np.column_stack((wh[:,keep],new_wh))
            order=np.argsort(times);times,amplitudes,h,wh=times[order],amplitudes[order],h[:,order],wh[:,order]
        wr=new_wr;residual=target-h@amplitudes;additions+=max(delta_k,0)
        changes.append(dict(kind=kind,old_times_ns=old_times.tolist(),new_times_ns=tt.tolist(),gain=gain,penalty=penalty,delta_k=delta_k))
        return True

    for iteration in range(config.passes):
        for j in range(max(0,len(times)-1)):
            if times[j+1]-times[j]>config.pair_distance_ns:continue
            ix=np.array([j,j+1]);seeds=[times[ix]]+[times[ix]+np.array([s,-s]) for s in (-4.,4.)]
            update(ix,propose(ix,seeds=seeds),'same_k')
        if len(times):
            amplitudes=nnls(wh,wy,maxiter=max(100,20*len(times)))[0]
            wr=wy-wh@amplitudes;residual=target-h@amplitudes
        from dataclasses import replace
        discovery=replace(baseline,min_charge_pe=config.report_charge_pe/2,min_separation_ns=config.min_separation_ns)
        extra,_=candidate_peaks(t,residual,template,discovery)
        ranked=[]
        for onset in extra:
            if len(times) and np.min(abs(times-onset))<config.min_separation_ns:continue
            col=metric.apply(template(t-onset));ranked.append((max(float(col@wr),0)**2/max(float(col@col),1e-30),onset))
        for _,onset in sorted(ranked,reverse=True)[:config.residual_budget]:
            if additions>=config.max_additions:break
            if len(times) and np.min(abs(times-onset))<config.min_separation_ns:continue
            ix=np.array([int(np.argmin(abs(times-onset)))],int) if len(times) else np.array([],int)
            update(ix,propose(ix,[onset],[np.sort(np.r_[times[ix],onset])]),'residual_add')
        for anchor in times.copy():
            if additions>=config.max_additions or not len(times):break
            j=int(np.argmin(abs(times-anchor)))
            if abs(times[j]-anchor)>config.move_ns:continue
            ix=np.array([j]);seeds=[]
            for gap in config.split_gaps_ns:
                seeds.extend(([times[j]-gap/2,times[j]+gap/2],[times[j],times[j]+gap],[times[j]-gap,times[j]]))
            update(ix,propose(ix,seeds=seeds),'split')
        j=0
        while j<len(times):
            ix=np.array([j]);empty=dict(times=np.array([]),amplitude=np.array([]))
            if update(ix,empty,'remove'):continue
            if j+1<len(times) and times[j+1]-times[j]<=config.pair_distance_ns:
                ix=np.array([j,j+1]);center=float(np.average(times[ix],weights=np.maximum(amplitudes[ix],1e-20)))
                if update(ix,propose(ix,seeds=[[center]]),'merge'):continue
            j+=1
        if len(times):
            amplitudes=nnls(wh,wy,maxiter=max(100,20*len(times)))[0]
            wr=wy-wh@amplitudes;residual=target-h@amplitudes
    charge=amplitudes/peak_per_pe;keep=charge>=config.report_charge_pe
    return dict(status='fitted',times_ns=times[keep].tolist(),peak_mV=amplitudes[keep].tolist(),
        charge_pe=charge[keep].tolist(),area_mV_ns=(amplitudes[keep]/template.unit_area_peak).tolist(),
        all_components=dict(times_ns=times.tolist(),peak_mV=amplitudes.tolist(),charge_pe=charge.tolist()),
        model_mV=target-residual+baseline.baseline_mV,rms_mV=float(np.sqrt(np.mean(residual**2))),
        whitened_cost=float(wr@wr),logdet_covariance=metric.logdet,local_changes=changes,
        selection_diagnostics=diagnostics,optimizer_trials=trials,
        global_calls_used=config.global_calls-global_remaining,refinement_seconds=time.monotonic()-tick,
        global_time_bins_used=sorted(global_bins),
        config=asdict(config),covariance=asdict(calibration),covariance_pilot='Fixed unchanged observed-voltage baseline fit',
        penalty_note='Empirical single-pulse calibrated full-record whitened gain, not a chi-square test',
        optimization_note='Local marginal-covariance proposals; accepted on full-record fixed-covariance objective',
        fit_truth_inputs=[])
