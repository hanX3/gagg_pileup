"""Regularized within-class pulse shapes, initialized from voltage alone.

Each event has a nonnegative charge and a slow fraction within three calibrated
RMS deviations of its class mean. The empirical quadratic penalty uses a fixed
voltage-pilot charge scale. It is not a calibrated Bayesian class probability.
Times, noise weights and charge-dependent shape centers remain fixed.
"""
from dataclasses import asdict, dataclass
import numpy as np
from scipy.optimize import nnls
from .alpha_gamma_regularized_fit import CalibratedBank, FrozenWhitening
from .alpha_gamma_model_selection_fit import fit_model_selection, ModelSelectionConfig
from .alpha_gamma_nnls import reduced_nnls


@dataclass(frozen=True)
class ShapeConfig:
    sigma_scale: float = 1.
    fraction_bound_sigma: float = 3.
    pilot_charge_floor_pe: float = 250.
    passes: int = 3
    pair_distance_ns: float = 250.
    fit_stride: int = 2


def fraction_center(calibration, kind, charge):
    p=calibration['classes'][('alpha','gamma')[int(kind)]]
    q=np.clip(max(charge,1.),*p['charge_range_pe'])
    return float(np.clip(p['slow_intercept']+p['slow_log_slope']*np.log(q/1000.),.02,.98))


def shape_matrix(bank,t,times,kinds,fractions):
    columns=[]
    for onset,kind,slow in zip(times,kinds,fractions):
        delay=bank.calibration['classes'][bank.names[int(kind)]]['delay_ns']
        x=t-onset-delay
        columns.append(bank.area*((1-slow)*bank.basis[0](x)*bank.basis[0].unit_area_peak+
                                  slow*bank.basis[1](x)*bank.basis[1].unit_area_peak))
    return np.column_stack(columns) if columns else np.empty((len(t),0))


def fit_shape_profile(time_ns,voltage_mV,calibration,config=ShapeConfig(),*,initial=None,
                      center_grid=None,sigma_grid=None):
    if config.sigma_scale<=0 or config.fraction_bound_sigma<=0 or config.pilot_charge_floor_pe<=0:
        raise ValueError('Positive shape regularization scales required')
    full_t=np.asarray(time_ns,float);full_y=np.asarray(voltage_mV,float)
    if initial is None:
        initial=fit_model_selection(full_t,full_y,calibration,ModelSelectionConfig(addition_cost=1.,fit_stride=config.fit_stride))
    if initial['status']!='fitted' or not initial['all_components']['times_ns']:return initial
    bank=CalibratedBank(calibration);c=initial['all_components'];times=np.array(c['times_ns'])
    kinds=np.array([0 if p=='alpha' else 1 for p in c['particle_name']]);q0=np.array(c['charge_pe']);shapes=np.array(c['shape_charge_pe']);n=len(times)
    t=full_t[::config.fit_stride];y=full_y[::config.fit_stride]-calibration['baseline_mV']
    W=FrozenWhitening(bank.matrix(t,times,kinds,shapes)@q0,t[1]-t[0],calibration['noise'],smooth_ns=0.)
    centers=(np.array([[fraction_center(calibration,k,q) for k in (0,1)] for q in shapes])
             if center_grid is None else np.asarray(center_grid,float))
    sigmas=config.sigma_scale*(np.broadcast_to(
        np.array([calibration['classes'][p]['rms_fraction_residual'] for p in bank.names]),(n,2))
        if sigma_grid is None else np.asarray(sigma_grid,float))
    if centers.shape!=(n,2) or sigmas.shape!=(n,2) or not np.isfinite(centers).all() or not np.isfinite(sigmas).all():
        raise ValueError('Shape parameter grids must be finite (components, 2) arrays')
    if np.any(sigmas<=0) or np.any(centers<.02) or np.any(centers>.98):
        raise ValueError('Positive shape RMS and centers in [.02,.98] required')
    bounds=np.clip(centers[:,:,None]+config.fraction_bound_sigma*sigmas[:,:,None]*np.array([-1.,1.]),.02,.98)
    design=W.apply(shape_matrix(bank,t,np.repeat(times,4),np.tile([0,0,1,1],n),bounds.ravel()))
    penalty=np.zeros((n,4*n))
    for j in range(n):
        penalty[j,4*j:4*j+4]=((bounds[j]-centers[j,:,None])/(sigmas[j,:,None]*max(q0[j],config.pilot_charge_floor_pe))).ravel()
    dictionary=np.vstack([design,penalty]);target=np.r_[W.apply(y),np.zeros(n)]
    columns=(4*np.arange(n)[:,None]+2*kinds[:,None]+np.arange(2)).ravel()
    # Represent the pilot center exactly, even when fraction limits are clipped.
    selected_bounds=bounds[np.arange(n),kinds];selected_centers=centers[np.arange(n),kinds]
    upper_weight=(selected_centers-selected_bounds[:,0])/np.diff(selected_bounds,axis=1)[:,0]
    amps=np.column_stack([q0*(1-upper_weight),q0*upper_weight]).ravel()
    H=dictionary[:,columns].copy();residual=target-H@amps;history=[float(residual@residual)];changes=[]
    for iteration in range(config.passes):
        before=history[-1]
        for j in np.argsort(amps.reshape(-1,2).sum(axis=1))[::-1]:
            ix=np.arange(2*j,2*j+2);local=residual+H[:,ix]@amps[ix];best=float(residual@residual)
            for kind in (0,1):
                cc=np.arange(4*j+2*kind,4*j+2*kind+2);h=dictionary[:,cc]
                q=nnls(h,local,maxiter=100)[0];r=local-h@q;cost=float(r@r)
                if cost<best-1e-9*max(1.,best):
                    best=cost;columns[ix]=cc;H[:,ix]=h;amps[ix]=q;residual=r
        for j in range(n-1):
            if times[j+1]-times[j]>config.pair_distance_ns or max(amps.reshape(-1,2).sum(axis=1)[j:j+2])<250.:continue
            ix=np.arange(2*j,2*j+4);local=residual+H[:,ix]@amps[ix];best=float(residual@residual)
            for a,b in ((0,0),(0,1),(1,0),(1,1)):
                cc=np.r_[np.arange(4*j+2*a,4*j+2*a+2),np.arange(4*(j+1)+2*b,4*(j+1)+2*b+2)]
                h=dictionary[:,cc];q=nnls(h,local,maxiter=100)[0];r=local-h@q;cost=float(r@r)
                if cost<best-1e-9*max(1.,best):
                    best=cost;columns[ix]=cc;H[:,ix]=h;amps[ix]=q;residual=r
        amps=reduced_nnls(H,target,maxiter=max(200,40*n))[0];residual=target-H@amps
        history.append(float(residual@residual));changes.append(int(np.sum((columns.reshape(-1,2)[:,0]%4)//2!=kinds)))
        if before-history[-1]<=1e-8*max(1.,before):break
    kinds=(columns.reshape(-1,2)[:,0]%4)//2;q=amps.reshape(-1,2).sum(axis=1)
    fractions=np.divide(np.sum(amps.reshape(-1,2)*bounds[np.arange(n),kinds],axis=1),q,out=centers[np.arange(n),kinds].copy(),where=q>0)
    data_cost=float(residual[:len(t)]@residual[:len(t)]);shape_cost=float(residual[len(t):]@residual[len(t):])
    keep=q>1e-6;times,kinds,q,shapes,fractions=times[keep],kinds[keep],q[keep],shapes[keep],fractions[keep]
    model=calibration['baseline_mV']+shape_matrix(bank,full_t,times,kinds,fractions)@q;report=q>=500.
    return dict(status='fitted',times_ns=times[report].tolist(),particle_name=bank.names[kinds[report]].tolist(),charge_pe=q[report].tolist(),
        area_mV_ns=(q[report]*bank.area).tolist(),model_mV=model,
        all_components=dict(times_ns=times.tolist(),particle_name=bank.names[kinds].tolist(),charge_pe=q.tolist(),shape_charge_pe=shapes.tolist(),slow_fraction=fractions.tolist()),
        rms_mV=float(np.sqrt(np.mean((full_y-model)**2))),config=asdict(config),fit_truth_inputs=[],selection=initial['selection'],
        objective_history=history,data_cost=data_cost,shape_penalty_cost=shape_cost,type_changes_per_pass=changes,
        model_includes_nuisance_components=True,nuisance_component_count=int(np.sum(~report)),selected_components=len(times),
        objective='Fixed voltage-pilot times/covariance; class-specific bounded shape variation with isolated-calibration RMS regularization',
        shape_uncertainty='Empirical isolated calibration RMS includes finite-photon and fit variation; no calibrated class posterior')
