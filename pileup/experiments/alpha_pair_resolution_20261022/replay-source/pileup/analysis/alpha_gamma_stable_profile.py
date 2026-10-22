"""Profile all charges while stabilizing each voltage-derived leading edge.

The regularizer constrains a pilot rate jump, not the particle class or charge.
Its strength is chosen on development data; it is not an independent Bayesian
prior, calibrated posterior, or new physical detector parameter.
"""
from dataclasses import dataclass, asdict
import numpy as np
from .alpha_gamma_regularized_fit import CalibratedBank, FrozenWhitening
from .alpha_gamma_model_selection_fit import fit_model_selection, ModelSelectionConfig
from .alpha_gamma_profile_fit import ProfileSystem


@dataclass(frozen=True)
class StableConfig:
    relative_jump_scale: float = .5
    absolute_jump_scale_pe_per_ns: float = .5
    passes: int = 3
    pair_distance_ns: float = 250.
    fit_stride: int = 2


def jump_coefficient(calibration,kind,shape):
    p=calibration['classes'][('alpha','gamma')[int(kind)]]
    q=np.clip(max(shape,1.),*p['charge_range_pe'])
    slow=np.clip(p['slow_intercept']+p['slow_log_slope']*np.log(q/1000.),.02,.98)
    return (1-slow)/95.+slow/351.


def fit_stable(time_ns,voltage_mV,calibration,config=StableConfig(),*,initial=None):
    if config.relative_jump_scale<=0 or config.absolute_jump_scale_pe_per_ns<=0:
        raise ValueError('Positive stabilizer scales required')
    full_t=np.asarray(time_ns,float);full_y=np.asarray(voltage_mV,float)
    if initial is None:
        initial=fit_model_selection(full_t,full_y,calibration,ModelSelectionConfig(addition_cost=1.,fit_stride=config.fit_stride))
    if initial['status']!='fitted' or not initial['all_components']['times_ns']:return initial
    bank=CalibratedBank(calibration);c=initial['all_components'];times=np.array(c['times_ns'])
    kinds=np.array([0 if k=='alpha' else 1 for k in c['particle_name']]);q0=np.array(c['charge_pe']);shapes=np.array(c['shape_charge_pe'])
    t=full_t[::config.fit_stride];y=full_y[::config.fit_stride]-calibration['baseline_mV']
    pilot=bank.matrix(t,times,kinds,shapes)@q0
    W=FrozenWhitening(pilot,t[1]-t[0],calibration['noise'],smooth_ns=0.)
    target=W.apply(y);n=len(times)
    design=W.apply(bank.matrix(t,np.repeat(times,2),np.tile([0,1],n),np.repeat(shapes,2)))
    coefficients=np.array([[jump_coefficient(calibration,k,s) for k in (0,1)] for s in shapes])
    jumps=q0*coefficients[np.arange(n),kinds]
    scales=np.hypot(config.relative_jump_scale*jumps,config.absolute_jump_scale_pe_per_ns)
    regularizer=np.zeros((n,2*n))
    regularizer[np.arange(n),2*np.arange(n)]=coefficients[:,0]/scales
    regularizer[np.arange(n),2*np.arange(n)+1]=coefficients[:,1]/scales
    system=ProfileSystem(np.vstack([design,regularizer]),np.r_[target,jumps/scales])
    columns=2*np.arange(n)+kinds;charges,cost=system.solve(columns);history=[cost];changes=[]
    for iteration in range(config.passes):
        old_cost=cost
        for j in np.argsort(charges)[::-1]:
            trial=columns.copy();trial[j]^=1
            q,value=system.solve(trial)
            if value<cost-1e-8*max(1.,cost):
                changes.append(dict(pass_index=iteration,indices=[int(j)],decrease=cost-value))
                columns,charges,cost=trial,q,value
        for j in range(n-1):
            if times[j+1]-times[j]>config.pair_distance_ns or max(charges[j:j+2])<250.:continue
            best=None
            for a,b in ((0,0),(0,1),(1,0),(1,1)):
                trial=columns.copy();trial[j:j+2]=[2*j+a,2*(j+1)+b]
                if np.array_equal(trial,columns):continue
                q,value=system.solve(trial)
                if value<cost-1e-8*max(1.,cost):best=(trial,q,value);cost=value
            if best is not None:
                columns,charges,cost=best;changes.append(dict(pass_index=iteration,indices=[j,j+1]))
        history.append(cost)
        if old_cost-cost<=1e-8*max(1.,old_cost):break
    kinds=columns%2;data_cost=float(np.sum((design[:,columns]@charges-target)**2))
    prior_cost=float(np.sum((regularizer[:,columns]@charges-jumps/scales)**2))
    keep=charges>1e-6;times,kinds,charges,shapes=times[keep],kinds[keep],charges[keep],shapes[keep]
    model=calibration['baseline_mV']+bank.matrix(full_t,times,kinds,shapes)@charges;report=charges>=500.
    return dict(status='fitted',times_ns=times[report].tolist(),particle_name=bank.names[kinds[report]].tolist(),charge_pe=charges[report].tolist(),
        area_mV_ns=(charges[report]*bank.area).tolist(),model_mV=model,
        all_components=dict(times_ns=times.tolist(),particle_name=bank.names[kinds].tolist(),charge_pe=charges.tolist(),shape_charge_pe=shapes.tolist()),
        rms_mV=float(np.sqrt(np.mean((full_y-model)**2))),config=asdict(config),fit_truth_inputs=[],selection=initial['selection'],
        objective_history=history,data_cost=data_cost,stabilizer_cost=prior_cost,type_changes=changes,
        profile_nnls_calls=system.calls,model_includes_nuisance_components=True,nuisance_component_count=int(np.sum(~report)),
        selected_components=len(times),objective='Fixed voltage-derived times and covariance; full-charge profile with pilot leading-rate regularizer',
        stabilizer='Same-data empirical regularization, not an independently calibrated probabilistic prior')
