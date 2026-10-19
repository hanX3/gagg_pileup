"""Neural event candidates followed by genuine physical-template fitting.

No truth, injected rate or known pulse multiplicity is accepted. A configurable
amplitude prior comes only from the waveform-derived neural prediction. Zero
prior is ordinary nonnegative physical fitting; positive prior regularizes
ill-conditioned overlapping components. Times float around neural candidates.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
import numpy as np
from scipy.optimize import minimize,minimize_scalar,nnls
from scipy.signal import find_peaks
from .alpha_pileup_fit import FitConfig
from .alpha_full_signal_study import make_template


@dataclass(frozen=True)
class HybridConfig:
    threshold: float=.65
    min_neural_charge_pe: float=250.
    prior_strength: float=0.
    time_refine_ns: float=6.
    time_passes: int=1
    min_separation_ns: float=8.
    max_optimizer_iterations: int=400


def merge_candidates(times,amplitudes,separation):
    order=np.argsort(times);times=np.asarray(times)[order];amps=np.asarray(amplitudes)[order]
    groups=np.split(np.arange(len(times)),np.flatnonzero(np.diff(times)>=separation)+1)
    if not len(times):return times,amps
    return np.array([np.average(times[g],weights=np.maximum(amps[g],1e-12)) for g in groups]),np.array([amps[g].sum() for g in groups])


def amplitude_fit(H,values,prior,initial,strength,max_iterations):
    """Solve a convex, nonnegative physical data + amplitude-prior objective.

    Prior penalty is lambda*sum_j ||H_j||^2*(a_j-a_neural_j)^2. Small
    matrices use exact active-set NNLS; larger ones use analytic Gram gradients
    and bounded L-BFGS-B. Solver termination/KKT quality is always recorded.
    """
    n=H.shape[1]
    if not n:return np.array([]),dict(converged=True,iterations=0,kkt_relative=0.)
    norms=np.sqrt(np.sum(H*H,axis=0));safe=np.maximum(norms,1e-12)
    if n<=64:
        penalty=np.sqrt(strength)*safe
        matrix=np.vstack((H,np.diag(penalty))) if strength else H
        target=np.r_[values,penalty*prior] if strength else values
        amps=nnls(matrix,target,maxiter=max(1000,20*n))[0]
        gradient=H.T@(H@amps-values)+strength*norms**2*(amps-prior)
        projected=np.where(amps>1e-9,gradient,np.minimum(gradient,0.))
        quality=float(np.max(abs(projected))/max(1.,np.max(abs(H.T@values))))
        return amps,dict(converged=True,iterations=None,kkt_relative=quality,solver='active_set_nnls')
    B=H/safe
    G=B.T@B;b=B.T@values;up=prior*safe;constant=.5*float(values@values)
    def objective(u):
        d=u-up
        value=.5*float(u@(G@u))-float(b@u)+constant+.5*strength*float(d@d)
        gradient=G@u-b+strength*d
        return value,gradient
    result=minimize(objective,np.maximum(initial,0.)*safe,jac=True,method='L-BFGS-B',
        bounds=[(0.,None)]*n,options={'maxiter':max_iterations,'ftol':1e-11,'gtol':1e-5,'maxcor':20,'maxls':40})
    _,gradient=objective(result.x)
    projected=np.where(result.x>1e-9,gradient,np.minimum(gradient,0.))
    quality=float(np.max(abs(projected))/max(1.,np.max(abs(b))))
    return result.x/safe,dict(converged=bool(result.success),iterations=int(result.nit),
        kkt_relative=quality,solver='bounded_gram_lbfgs',termination=str(result.message))


def fit_candidates(t,voltage,times,neural_charges,settings,config):
    """Fit measured voltage using unknown-count neural candidates, not truth."""
    cfg=FitConfig(**settings['fit_config']);template=make_template(settings)
    t=np.asarray(t,float);values=np.asarray(voltage,float)-cfg.baseline_mV
    times=np.asarray(times,float);prior=np.asarray(neural_charges,float)*cfg.area_per_pe_mV_ns*template.unit_area_peak
    if values.shape!=t.shape or times.shape!=prior.shape or not np.isfinite(values).all():
        raise ValueError('Invalid measured input or neural candidates')
    if config.prior_strength<0 or config.time_refine_ns<0:
        raise ValueError('Invalid physical-fit configuration')
    valid=np.isfinite(times)&np.isfinite(prior)&(times>=0)&(times<=t[-1])&(prior>0)
    times,prior=merge_candidates(times[valid],prior[valid],config.min_separation_ns)
    initial_count=len(times);start=times.copy()
    H=template(t[:,None]-times[None,:])
    amps,solver=amplitude_fit(H,values,prior,prior,config.prior_strength,config.max_optimizer_iterations)
    solvers=[solver];residual=values-H@amps
    history=[dict(data_sse_mV2=float(residual@residual),
                  prior_penalty=float(config.prior_strength*np.sum(np.sum(H*H,axis=0)*(amps-prior)**2)))]
    for _ in range(config.time_passes):
        for j in range(len(times)):
            if amps[j]<.05*cfg.min_charge_pe*cfg.area_per_pe_mV_ns*template.unit_area_peak:
                continue
            mask=(t>=start[j]-config.time_refine_ns-4)&(t<=start[j]+6000+config.time_refine_ns)
            tx=t[mask];target=residual[mask]+H[mask,j]*amps[j]
            gamma=config.prior_strength*float(H[:,j]@H[:,j])
            def evaluate(time,return_amplitude=False):
                h=template(tx-time);denominator=float(h@h)+gamma
                amplitude=max(0.,(float(h@target)+gamma*prior[j])/denominator) if denominator else 0.
                error=target-amplitude*h
                return amplitude if return_amplitude else float(error@error)+gamma*(amplitude-prior[j])**2
            lo,hi=max(0.,start[j]-config.time_refine_ns),min(t[-1],start[j]+config.time_refine_ns)
            if hi<=lo:continue
            trial=minimize_scalar(evaluate,bounds=(lo,hi),method='bounded',options={'xatol':.05,'maxiter':25})
            old_cost=evaluate(times[j])
            new_time=trial.x if trial.fun<old_cost else times[j]
            new_amplitude=evaluate(new_time,True)
            residual+=H[:,j]*amps[j]
            times[j],amps[j]=new_time,new_amplitude
            H[:,j]=template(t-times[j]);residual-=H[:,j]*amps[j]
        amps,solver=amplitude_fit(H,values,prior,amps,config.prior_strength,config.max_optimizer_iterations)
        solvers.append(solver);residual=values-H@amps
        history.append(dict(data_sse_mV2=float(residual@residual),
                      prior_penalty=float(config.prior_strength*np.sum(np.sum(H*H,axis=0)*(amps-prior)**2))))
    order=np.argsort(times);times,amps,prior=times[order],amps[order],prior[order]
    if len(times)>1 and np.any(np.diff(times)<config.min_separation_ns):
        groups=np.split(np.arange(len(times)),np.flatnonzero(np.diff(times)>=config.min_separation_ns)+1)
        merged=np.array([np.average(times[g],weights=np.maximum(amps[g],1e-12)) for g in groups])
        prior=np.array([prior[g].sum() for g in groups]);initial=np.array([amps[g].sum() for g in groups])
        times=merged;H=template(t[:,None]-times)
        amps,solver=amplitude_fit(H,values,prior,initial,config.prior_strength,config.max_optimizer_iterations);solvers.append(solver)
    minimum=cfg.min_charge_pe*cfg.area_per_pe_mV_ns*template.unit_area_peak
    keep=amps>=minimum
    if not np.all(keep):
        times,prior,amps=times[keep],prior[keep],amps[keep]
        H=template(t[:,None]-times)
        amps,solver=amplitude_fit(H,values,prior,amps,config.prior_strength,config.max_optimizer_iterations);solvers.append(solver)
    keep=amps>=minimum;times,amps=times[keep],amps[keep]
    model=cfg.baseline_mV+template(t[:,None]-times)@amps
    return dict(status='hybrid_fitted',times_ns=times.tolist(),peak_mV=amps.tolist(),
        charge_pe=(amps/template.unit_area_peak/cfg.area_per_pe_mV_ns).tolist(),
        area_mV_ns=(amps/template.unit_area_peak).tolist(),model_mV=model,
        rms_mV=float(np.sqrt(np.mean((voltage-model)**2))),candidates=initial_count,
        fit_config=asdict(config),solver_runs=solvers,cost_history=history,
        limitations='Cannot discover an event absent from every neural candidate; fixed mean template and finite optimizer budget')


def hybrid_waveform(t,voltage,dense,settings,config=HybridConfig()):
    peaks,_=find_peaks(dense['probability'],height=config.threshold,distance=4)
    peaks=peaks[dense['charge_pe'][peaks]>=config.min_neural_charge_pe]
    times=t[peaks]+dense['time_offset_ns'][peaks]
    return fit_candidates(t,voltage,times,dense['charge_pe'][peaks],settings,config)
