"""Blind pure-alpha template fits and separate truth scoring.

The fitter accepts only sampled voltages, their time grid and a frozen template.
It never accepts source rate, pulse count, event times, energies or truth labels.
This is a project baseline, not a reproduction of a paper's complete algorithm.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import numpy as np
from scipy.fft import rfft, irfft, rfftfreq, next_fast_len
from scipy.optimize import nnls, minimize_scalar, linear_sum_assignment, least_squares
from scipy.signal import find_peaks

from analysis.alpha_gamma_waveform_fit import sampled_component


@dataclass(frozen=True)
class FitConfig:
    baseline_mV: float = 10.
    min_charge_pe: float = 500.
    area_per_pe_mV_ns: float = 13.618501389
    deconvolution_sigma_ns: float = 8.
    deconvolution_ridge_fraction: float = .0001
    candidate_fraction: float = .65
    min_separation_ns: float = 8.
    time_refine_ns: float = 12.
    max_passes: int = 5
    residual_rounds: int = 1
    max_candidates: int = 500


class AlphaTemplate:
    def __init__(self, readout, slow_fraction=.740166, tau_fast_ns=95., tau_slow_ns=351.):
        self.readout = dict(readout)
        self.slow_fraction = float(slow_fraction)
        self.tau_fast_ns, self.tau_slow_ns = tau_fast_ns, tau_slow_ns
        self.grid = np.arange(-4., 6000.001, .25)
        h = ((1-slow_fraction)*sampled_component(self.grid,0.,tau_fast_ns,readout)
             + slow_fraction*sampled_component(self.grid,0.,tau_slow_ns,readout))
        self.unit_area_peak = float(h.max())
        self.shape = h/self.unit_area_peak
        self.slope = np.gradient(self.shape,.25)

    def __call__(self, relative_ns):
        return np.interp(relative_ns,self.grid,self.shape,left=0.,right=0.)

    def metadata(self):
        return dict(slow_fraction=self.slow_fraction, tau_fast_ns=self.tau_fast_ns,
                    tau_slow_ns=self.tau_slow_ns, area_per_peak_mV_ns=1/self.unit_area_peak,
                    readout=self.readout, interpolation_ns=.25, support_ns=6000.)


def candidate_peaks(t, residual, template, config):
    """Regularized inverse filter used for discovery only; amplitude is refit."""
    dt = float(t[1]-t[0]); n = len(t)
    nfft = next_fast_len(2*n)
    h = template(np.arange(n)*dt+t[0])
    H = rfft(h,nfft)
    gaussian = np.exp(-.5*(2*np.pi*rfftfreq(nfft,dt)*config.deconvolution_sigma_ns)**2)
    inverse = gaussian*H.conj()/(abs(H)**2+(config.deconvolution_ridge_fraction*abs(H).max())**2)
    resolution = irfft(H*inverse,nfft)
    shifted = irfft(rfft(residual,nfft)*inverse,nfft)[:n]/resolution.max()
    min_peak = config.min_charge_pe*config.area_per_pe_mV_ns*template.unit_area_peak
    peaks,_ = find_peaks(shifted, height=min_peak*config.candidate_fraction,
                        prominence=min_peak*config.candidate_fraction*.65,
                        distance=max(1,int(np.ceil(config.min_separation_ns/dt))))
    return peaks*dt, shifted


def refit(t, values, initial_times, template, config):
    """Nonnegative amplitudes and floating times by cyclic variable projection.

    Each time step minimizes the residual after adding back that component;
    an unpenalized joint NNLS solve removes discovery-filter amplitude bias.
    The calibrated baseline is fixed; drift sensitivity needs a later study.
    """
    times=np.sort(np.asarray(initial_times,float))
    if len(times)==0: return times,np.array([]),np.zeros_like(values),[]
    start=times.copy()
    H=template(t[:,None]-times[None,:])
    amps=nnls(H,values,maxiter=max(100,10*len(times)))[0]
    residual=values-H@amps
    costs=[float(residual@residual)]
    for _ in range(config.max_passes):
        for j in range(len(times)):
            mask=(t>=start[j]-config.time_refine_ns-4.)&(t<=start[j]+6000.+config.time_refine_ns)
            tx=t[mask]; target=residual[mask]+H[mask,j]*amps[j]
            def evaluate(time, return_amplitude=False):
                v=template(tx-time)
                denominator=float(v@v)
                amplitude=max(0.,float(v@target)/denominator) if denominator else 0.
                error=target-amplitude*v
                return amplitude if return_amplitude else float(error@error)
            trial=minimize_scalar(evaluate,bounds=(max(0.,start[j]-config.time_refine_ns),
                                                  min(t[-1],start[j]+config.time_refine_ns)),
                                  method='bounded', options={'xatol':.02,'maxiter':35})
            # Bounded scalar minimization is not guaranteed to improve an old solution.
            if trial.fun < evaluate(times[j]):
                residual += H[:,j]*amps[j]
                times[j]=trial.x; amps[j]=evaluate(trial.x,True)
                H[:,j]=template(t-times[j]); residual -= H[:,j]*amps[j]
        amps=nnls(H,values,maxiter=max(100,10*len(times)))[0]
        residual=values-H@amps
        costs.append(float(residual@residual))
        if costs[-2]-costs[-1] <= 1e-5*max(1.,costs[-2]): break
    order=np.argsort(times)
    if len(times)<=64:
        n=len(times)
        def fun(z):
            return template(t[:,None]-z[n:][None,:])@z[:n]-values
        def jac(z):
            u=t[:,None]-z[n:][None,:]
            H=template(u)
            dh=np.interp(u,template.grid,template.slope,left=0.,right=0.)
            return np.c_[H,-dh*z[:n]]
        solution=least_squares(fun,np.r_[amps,times],jac=jac,
            bounds=(np.r_[np.zeros(n),np.maximum(0.,start-config.time_refine_ns)],
                    np.r_[np.full(n,np.inf),np.minimum(t[-1],start+config.time_refine_ns)]),
            x_scale=np.r_[np.full(n,30.),np.full(n,5.)],max_nfev=60,ftol=1e-8,xtol=1e-7,gtol=1e-7)
        if 2*solution.cost < costs[-1]:
            amps,times=solution.x[:n],solution.x[n:]
            H=template(t[:,None]-times[None,:])
            costs.append(float(2*solution.cost))
            order=np.argsort(times)
    return times[order],amps[order],H@amps,costs


def fit_waveform(t_ns, voltage_mV, template, config=FitConfig(), *, rail_mask=None):
    """Fit the supplied total voltage without an intrinsic ADC range limit.

    Passing rail_mask explicitly opts into the historical whole-record rejection
    policy. Studies assuming adequate dynamic range supply complete sampled
    voltage without that mask; clipped ADC values cannot recreate this input.
    """
    t,y=np.asarray(t_ns,float),np.asarray(voltage_mV,float)
    if t.ndim!=1 or y.shape!=t.shape or len(t)<10 or not np.isfinite(y).all():
        raise ValueError('Invalid input waveform')
    if not np.allclose(np.diff(t),t[1]-t[0],atol=1e-8,rtol=0):
        raise ValueError('Uniform samples required')
    if abs(t[0])>t[1]-t[0]: raise ValueError('Expected a full record beginning near zero')
    if rail_mask is not None and np.asarray(rail_mask,bool).any():
        return dict(status='adc_rail_rejected',times_ns=[],peak_mV=[],area_mV_ns=[],
                    charge_pe=[],model_mV=np.full_like(y,config.baseline_mV),
                    candidates=0,cost_history=[],rail_fraction=float(np.mean(rail_mask)))
    target=y-config.baseline_mV
    times,_=candidate_peaks(t,target,template,config)
    count_initial=len(times)
    if len(times)>config.max_candidates:
        return dict(status='candidate_limit',times_ns=[],peak_mV=[],area_mV_ns=[],
                    charge_pe=[],model_mV=np.full_like(y,config.baseline_mV),
                    candidates=len(times),cost_history=[],rail_fraction=0.)
    def clean_refit(initial):
        fitted,amplitudes,prediction,history=refit(t,target,initial,template,config)
        # Collinear columns can otherwise represent a single pulse twice.
        while len(fitted)>1 and np.any(np.diff(fitted)<config.min_separation_ns):
            groups=np.split(np.arange(len(fitted)),np.flatnonzero(np.diff(fitted)>=config.min_separation_ns)+1)
            merged=[float(np.average(fitted[g],weights=np.maximum(amplitudes[g],1e-12))) for g in groups]
            fitted,amplitudes,prediction,last=refit(t,target,merged,template,config)
            history+=last
        return fitted,amplitudes,prediction,history
    times,amps,model,costs=clean_refit(times)
    for _ in range(config.residual_rounds):
        more,_=candidate_peaks(t,target-model,template,config)
        more=np.array([x for x in more if len(times)==0 or np.min(abs(times-x))>=config.min_separation_ns])
        if not len(more) or len(times)+len(more)>config.max_candidates: break
        times,amps,model,history=clean_refit(np.r_[times,more])
        costs+=history
    minimum=config.min_charge_pe*config.area_per_pe_mV_ns*template.unit_area_peak
    keep=amps>=minimum
    times,amps,model,last=clean_refit(times[keep])
    # The final selection is evaluated again after debiasing; do not return subthreshold components.
    keep=amps>=minimum
    times,amps=times[keep],amps[keep]
    model=template(t[:,None]-times[None,:])@amps if len(times) else np.zeros_like(t)
    area=amps/template.unit_area_peak
    return dict(status='fitted',times_ns=times.tolist(),peak_mV=amps.tolist(),
                area_mV_ns=area.tolist(),charge_pe=(area/config.area_per_pe_mV_ns).tolist(),
                model_mV=model+config.baseline_mV,candidates=count_initial,
                cost_history=costs+last,rail_fraction=0.,
                rms_mV=float(np.sqrt(np.mean((target-model)**2))))


def score(fit, truth, gate_ns, *, tolerance_ns=10., charge_threshold_pe=500.):
    """One-to-one time matching; charges are not used to choose assignments.

    Time/charge error statistics are conditional on matching. Recovery also
    requires charge error <=20%. A proximity graph separately counts unresolved
    groups (merge) and multiple candidates near one truth (split); these are
    operational flags, not a complete causal attribution at very high rate.
    """
    tt=np.asarray(truth['first_active_ns'],float)
    qt=np.asarray(truth['charge_pe'],float)
    central=(np.asarray(truth['emission_ns'])>=gate_ns[0])&(np.asarray(truth['emission_ns'])<gate_ns[1])
    detectable=(qt>=charge_threshold_pe)&np.isfinite(tt)
    target=central&detectable
    # Match against every signal including guards; filter denominators only afterwards.
    ix=np.flatnonzero(detectable); tx=tt[ix]
    tp=np.asarray(fit['times_ns']); qp=np.asarray(fit['charge_pe'])
    close=abs(tx[:,None]-tp[None,:])<=tolerance_ns
    n,m=len(tx),len(tp)
    # Dummies cost one; all valid time costs together are less than one,
    # so maximum matching cardinality dominates the tie-breaking time cost.
    cost=np.ones((n+m,n+m))*3.
    if n and m: cost[:n,:m]=np.where(close,abs(tx[:,None]-tp[None,:])/tolerance_ns/(n+m+1),3.)
    cost[:n,m:]=1.; cost[n:,:m]=1.; cost[n:,m:]=0.
    a,b=linear_sum_assignment(cost) if n+m else (np.array([],int),np.array([],int))
    pairs=[(int(ix[i]),int(j)) for i,j in zip(a,b) if i<n and j<m and close[i,j]]
    gate_pred=(tp>=gate_ns[0])&(tp<gate_ns[1])
    # Proximity gate uses first-arrival time, so report the small flight-time boundary convention.
    used={j for i,j in pairs}
    matched=[(i,j) for i,j in pairs if target[i]]
    dt=np.array([tp[j]-tt[i] for i,j in matched])
    dq=np.array([qp[j]/qt[i]-1 for i,j in matched])
    used_truth={i for i,j in pairs}
    # Flag a merge only when an unmatched truth lies near an already assigned
    # candidate, and a split only when an extra candidate lies near an assigned
    # truth. Perfectly recovered close pairs must not be labelled failures.
    merge_candidates={j for row,i in enumerate(ix) if target[i] and i not in used_truth
                      for j in used if close[row,j]}
    split_truth={i for row,i in enumerate(ix) if target[i] and i in used_truth
                 if any(close[row,j] for j in range(m) if j not in used)}
    merge=len(merge_candidates);split=len(split_truth)
    return dict(target_count=int(target.sum()),central_injected=int(central.sum()),
                central_nonzero_charge=int(np.sum(central&(qt>0))),
                reconstructed_count=int(gate_pred.sum()),matched_count=len(matched),
                recovered_count=int(np.sum(abs(dq)<=.20)),missed_count=int(target.sum())-len(matched),
                unmatched_predictions=int(sum(bool(gate_pred[j]) and j not in used for j in range(m))),
                merge_flags=merge,split_flags=split,time_errors_ns=dt.tolist(),
                relative_charge_errors=dq.tolist(),matched_pairs=matched,all_matched_pairs=pairs,
                target_indices=np.flatnonzero(target).tolist(),tolerance_ns=tolerance_ns,
                duration_ns=float(gate_ns[1]-gate_ns[0]))
