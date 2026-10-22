"""Calibrated, whitened sparse selection followed by discrete pulse refitting.

This is a project implementation inspired by the recorded literature, not a
paper replica. Reconstruction accepts voltage and frozen calibration only.
Selection weights use smoothed voltage; final weights use a voltage-only pilot
model and remain fixed during refitting. This is feasible weighted least
squares, not a parameter-dependent likelihood.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict, replace
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize, minimize_scalar, nnls

from .alpha_pileup_fit import AlphaTemplate, FitConfig
from .alpha_gamma_pileup_fit import TemplateBank, discover
from .alpha_gamma_nnls import reduced_nnls


@dataclass(frozen=True)
class RegularizedConfig:
    selection_penalty: float = 5.
    discovery_charge_pe: float = 80.
    report_charge_pe: float = 500.
    min_separation_ns: float = 8.
    time_refine_ns: float = 12.
    max_candidates: int = 700
    max_components: int = 400
    refit_passes: int = 2
    pair_distance_ns: float = 100.
    fit_stride: int = 2


class CalibratedBank:
    names=np.array(['alpha','gamma'])

    def __init__(self, calibration):
        self.calibration=calibration
        self.area=float(calibration['area_per_pe_mV_ns'])
        self.basis=[AlphaTemplate(calibration['readout'],slow_fraction=f) for f in (0.,1.)]

    def column(self, t, onset, kind, shape_charge):
        p=self.calibration['classes'][self.names[kind]]
        charge=np.clip(max(shape_charge,1.),*p['charge_range_pe'])
        slow=np.clip(p['slow_intercept']+p['slow_log_slope']*np.log(charge/1000.),.02,.98)
        x=t-onset-p['delay_ns']
        return self.area*((1-slow)*self.basis[0](x)*self.basis[0].unit_area_peak+
                          slow*self.basis[1](x)*self.basis[1].unit_area_peak)

    def matrix(self,t,times,types,shape_charges):
        return np.column_stack([self.column(t,x,int(k),q) for x,k,q in zip(times,types,shape_charges)]) if len(times) else np.empty((len(t),0))

    def discovery_bank(self):
        p=self.calibration['classes']
        return TemplateBank(self.calibration['readout'],alpha_slow=p['alpha']['slow_intercept'],
                            gamma_slow=p['gamma']['slow_intercept'])


class FrozenWhitening:
    """Invert a fixed, calibrated heteroscedastic AR(1) covariance factor."""
    def __init__(self, values, dt, noise, enabled=True, smooth_ns=80.):
        if not enabled:
            self.sigma=np.ones(len(values));self.rho=np.zeros(len(values)-1)
            return
        smooth=np.maximum(gaussian_filter1d(values,smooth_ns/dt) if smooth_ns else values,0.)
        shot=noise['variance_per_mV']*smooth
        variance=noise['baseline_variance_mV2']+shot
        self.sigma=np.sqrt(np.maximum(variance,1e-8))
        # Readout noise is nearly white; correlated shot noise dominates in pulses.
        rho=noise['shot_rho_at_2ns']**(dt/2.)
        self.rho=np.clip(rho*np.sqrt(shot[:-1]*shot[1:]/(variance[:-1]*variance[1:])),0.,.995)

    def apply(self, a):
        scale=self.sigma if a.ndim==1 else self.sigma[:,None]
        z=a/scale
        result=z.copy()
        rho=self.rho if a.ndim==1 else self.rho[:,None]
        result[1:]=(z[1:]-rho*z[:-1])/np.sqrt(1-rho*rho)
        return result


def positive_lasso(matrix, target, penalty, maxiter=1500):
    """Unit-norm columns, nonnegative L1 regression, explicit KKT diagnostic."""
    norms=np.linalg.norm(matrix,axis=0)
    valid=norms>1e-12
    if not np.all(valid):raise ValueError('Zero discovery column')
    design=matrix/norms
    gram=design.T@design;rhs=design.T@target
    def objective(x):
        grad=gram@x-rhs+penalty
        return .5*x@(gram@x)-rhs@x+penalty*x.sum(),grad
    solution=minimize(objective,np.zeros(len(rhs)),jac=True,method='L-BFGS-B',
        bounds=[(0.,None)]*len(rhs),options={'maxiter':maxiter,'ftol':1e-12,'gtol':1e-5,'maxls':40,'maxcor':20})
    x=solution.x;gradient=gram@x-rhs+penalty
    projected=np.where(x>1e-7,gradient,np.minimum(gradient,0.))
    return x/norms,dict(success=bool(solution.success),iterations=int(solution.nit),
        message=str(solution.message),projected_gradient_max=float(np.max(abs(projected))),
        objective=float(solution.fun),penalty=float(penalty))


def fit_regularized(time_ns, voltage_mV, calibration, config=RegularizedConfig(), *, weighted=True):
    full_t=np.asarray(time_ns,float);full_y=np.asarray(voltage_mV,float)
    if full_t.ndim!=1 or full_y.shape!=full_t.shape or len(full_t)<10 or not np.isfinite(full_y).all():
        raise ValueError('Invalid voltage record')
    if not np.allclose(np.diff(full_t),full_t[1]-full_t[0],atol=1e-8,rtol=0):raise ValueError('Uniform samples required')
    if not isinstance(config.fit_stride,int) or config.fit_stride<1:raise ValueError('Positive integer stride required')
    bank=CalibratedBank(calibration);baseline=float(calibration['baseline_mV'])
    t=full_t[::config.fit_stride];values=full_y[::config.fit_stride]-baseline
    W=FrozenWhitening(values,t[1]-t[0],calibration['noise'],weighted)
    target=W.apply(values)
    cfg=FitConfig(baseline_mV=baseline,min_charge_pe=config.discovery_charge_pe,
        area_per_pe_mV_ns=bank.area,min_separation_ns=config.min_separation_ns,
        max_candidates=config.max_candidates)
    # Discovery retains the original 2 ns samples; sparse/refit uses a declared stride.
    initial=discover(full_t,full_y-baseline,bank.discovery_bank(),cfg)
    if len(initial)>config.max_candidates:
        return _empty(full_t,full_y,baseline,'candidate_limit',len(initial))
    if not len(initial):return _empty(full_t,full_y,baseline,'fitted',0)
    dual_times=np.repeat(initial,2);dual_types=np.tile([0,1],len(initial))
    dual_shapes=np.full(len(dual_times),1000.)
    design=W.apply(bank.matrix(t,dual_times,dual_types,dual_shapes))
    sparse,selection=positive_lasso(design,target,config.selection_penalty)
    sparse=sparse.reshape(-1,2);q=sparse.sum(axis=1)
    keep=q>max(1.,config.discovery_charge_pe*.08)
    times=initial[keep];types=np.argmax(sparse[keep],axis=1);q=q[keep]
    if not len(times):return _empty(full_t,full_y,baseline,'fitted',len(initial),selection)
    if len(times)>config.max_components:return _empty(full_t,full_y,baseline,'component_limit',len(initial),selection)
    starts=times.copy();shapes=np.maximum(q,1.)
    H=W.apply(bank.matrix(t,times,types,shapes))
    q=reduced_nnls(H,target,maxiter=max(100,20*len(times)))[0]
    # Using noisy troughs themselves as inverse-variance weights can pull
    # charges down. Freeze the refit covariance from a smoother pilot mean,
    # built solely from the voltage and the selected support.
    pilot=bank.matrix(t,times,types,shapes)@q
    W=FrozenWhitening(pilot,t[1]-t[0],calibration['noise'],weighted,smooth_ns=0.)
    target=W.apply(values)
    H=W.apply(bank.matrix(t,times,types,shapes))
    q=reduced_nnls(H,target,maxiter=max(100,20*len(times)))[0]
    residual=target-H@q;history=[float(residual@residual)]
    # Conditional amplitude refits remove the selection L1 shrinkage. Weak
    # components remain in the model; the reporting cut is applied only at end.
    for iteration in range(config.refit_passes):
        # Charge-dependent shape is frozen within this pass; update only if
        # the jointly refitted weighted residual improves.
        trial_shapes=np.maximum(q,1.)
        trial_H=W.apply(bank.matrix(t,times,types,trial_shapes))
        trial_q=reduced_nnls(trial_H,target,maxiter=max(100,20*len(times)))[0]
        trial_r=target-trial_H@trial_q
        if trial_r@trial_r < residual@residual:
            H,q,shapes,residual=trial_H,trial_q,trial_shapes,trial_r
        for j in range(len(times)):
            local_target=residual+H[:,j]*q[j]
            def evaluate(x,k):
                h=W.apply(bank.column(t,x,k,shapes[j]));norm=float(h@h)
                amp=max(0.,float(h@local_target)/norm) if norm else 0.
                r=local_target-amp*h
                return float(r@r),amp,h
            best=(float(residual@residual),times[j],types[j],q[j],H[:,j])
            for k in (0,1):
                lo=max(0.,starts[j]-config.time_refine_ns)
                hi=min(full_t[-1],starts[j]+config.time_refine_ns)
                # Avoid coalescing adjacent events during continuous refinement.
                if j>0:lo=max(lo,times[j-1]+config.min_separation_ns)
                if j+1<len(times):hi=min(hi,times[j+1]-config.min_separation_ns)
                if hi<=lo:continue
                opt=minimize_scalar(lambda x:evaluate(x,k)[0],bounds=(lo,hi),method='bounded',
                    options={'xatol':.1,'maxiter':20})
                cost,amp,h=evaluate(opt.x,k)
                if cost<best[0]:best=(cost,float(opt.x),k,amp,h)
            _,times[j],types[j],q[j],H[:,j]=best
            residual=local_target-H[:,j]*q[j]
        # Adjacent labels can exchange charge; compare four assignments with
        # both amplitudes refitted, never use known truth classes or pairings.
        for a in range(len(times)-1):
            b=a+1
            if times[b]-times[a]>config.pair_distance_ns:continue
            ix=[a,b];local_target=residual+H[:,ix]@q[ix]
            best=float(residual@residual)
            for pair in ((0,0),(0,1),(1,0),(1,1)):
                if np.array_equal(types[ix],pair):continue
                h=W.apply(bank.matrix(t,times[ix],pair,shapes[ix]));amp=nnls(h,local_target,maxiter=100)[0]
                r=local_target-h@amp;cost=float(r@r)
                if cost<best:best=cost;types[ix]=pair;H[:,ix]=h;q[ix]=amp;residual=r
        q=reduced_nnls(H,target,maxiter=max(100,20*len(times)))[0]
        residual=target-H@q;history.append(float(residual@residual))
    active=q>1e-6;times,types,q,shapes=times[active],types[active],q[active],shapes[active]
    model=bank.matrix(full_t,times,types,shapes)@q+baseline
    report=q>=config.report_charge_pe
    all_components=dict(times_ns=times.tolist(),particle_name=bank.names[types].tolist(),
                        charge_pe=q.tolist(),shape_charge_pe=shapes.tolist())
    return dict(status='fitted',times_ns=times[report].tolist(),particle_name=bank.names[types[report]].tolist(),
        charge_pe=q[report].tolist(),area_mV_ns=(q[report]*bank.area).tolist(),
        all_components=all_components,model_mV=model,rms_mV=float(np.sqrt(np.mean((full_y-model)**2))),
        initial_candidates=len(initial),selected_components=len(times),
        nuisance_component_count=int(np.sum(~report)),model_includes_nuisance_components=True,
        selection=selection,weighted_cost_history=history,config=asdict(config),weighted=weighted,
        objective='Selection uses 80 ns smoothed voltage weights; refit uses fixed voltage-only pilot-mean AR(1) weights; L1 then conditional unpenalized refit',
        fit_truth_inputs=[])


def _empty(t,y,baseline,status,candidates,selection=None):
    return dict(status=status,times_ns=[],particle_name=[],charge_pe=[],area_mV_ns=[],
        all_components=dict(times_ns=[],particle_name=[],charge_pe=[],shape_charge_pe=[]),
        model_mV=np.full_like(t,baseline),rms_mV=float(np.sqrt(np.mean((y-baseline)**2))),
        initial_candidates=candidates,selected_components=0,nuisance_component_count=0,
        model_includes_nuisance_components=True,selection=selection,weighted_cost_history=[],fit_truth_inputs=[])
