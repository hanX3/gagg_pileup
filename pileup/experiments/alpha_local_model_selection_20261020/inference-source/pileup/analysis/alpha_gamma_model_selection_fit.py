"""Refine sparse candidates by residual proposals and joint time/type searches.

Separate experimental implementation; old frozen fits are unchanged. Receives
only voltage and calibration. Selection and refinement diagnostics are saved.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.optimize import nnls, minimize_scalar, least_squares
from .alpha_gamma_regularized_fit import (CalibratedBank, FrozenWhitening,
    RegularizedConfig, positive_lasso, _empty)
from .alpha_gamma_pileup_fit import discover
from .alpha_pileup_fit import FitConfig
from .alpha_gamma_nnls import reduced_nnls


def strict_positive_lasso(matrix, target, penalty, tolerance=1e-5):
    """Polish the active quadratic subproblem until its projected KKT check.

    The L1 objective is unchanged. No truth or regularizing ridge is added.
    A least-squares solve handles nearly dependent active dictionary columns.
    """
    charge, initial = positive_lasso(matrix, target, penalty)
    norms=np.linalg.norm(matrix,axis=0);a=matrix/norms
    gram=a.T@a;rhs=a.T@target-penalty
    x=charge*norms;active=x>1e-7
    iterations=0
    for iterations in range(300):
        gradient=gram@x-rhs
        projected=np.where(x>1e-7,gradient,np.minimum(gradient,0.))
        violation=float(np.max(abs(projected)))
        if violation<=tolerance:break
        inactive=np.flatnonzero(~active)
        if len(inactive) and np.min(gradient[inactive]) < -tolerance:
            active[inactive[np.argmin(gradient[inactive])]]=True
        for _ in range(len(x)+1):
            ix=np.flatnonzero(active)
            if not len(ix):break
            proposed=np.zeros_like(x)
            proposed[ix]=np.linalg.lstsq(gram[np.ix_(ix,ix)],rhs[ix],rcond=1e-12)[0]
            negative=ix[proposed[ix]<=0.]
            if not len(negative):x=proposed;break
            ratios=x[negative]/np.maximum(x[negative]-proposed[negative],1e-30)
            step=float(np.clip(np.min(ratios),0.,1.))
            x+=step*(proposed-x)
            drop=active&(x<=1e-8);active[drop]=False;x[drop]=0.
    gradient=gram@x-rhs
    projected=np.where(x>1e-7,gradient,np.minimum(gradient,0.))
    violation=float(np.max(abs(projected)))
    details=dict(initial,initial_projected_gradient_max=initial['projected_gradient_max'],
        projected_gradient_max=violation,success=violation<=tolerance,
        active_set_iterations=iterations+1,kkt_tolerance=tolerance,
        message='L-BFGS-B followed by active quadratic subproblem polishing',
        objective=float(.5*x@(gram@x)-rhs@x))
    return x/norms,details


@dataclass(frozen=True)
class ModelSelectionConfig:
    selection_penalty: float = 3.
    addition_cost: float = 4.
    pair_distance_ns: float = 500.
    joint_passes: int = 2
    max_additions: int = 12
    proposal_budget: int = 16
    pair_min_charge_pe: float = 250.
    fit_stride: int = 2


def _joint_pair(t, target, W, bank, times, types, charges, shapes, lower, upper):
    """Compare four classes, jointly optimizing two times and two charges.

    Shape-charge remains fixed within a trial, as in the preceding fitter.
    The second time is represented by first time plus a positive >=8 ns gap.
    """
    times=np.array(times,float);charges=np.array(charges,float)
    gap=times[1]-times[0];best=None;trials=0
    scale=bank.area
    for pair in ((0,0),(0,1),(1,0),(1,1)):
        H=W.apply(bank.matrix(t,times,pair,shapes))
        q=nnls(H,target,maxiter=100)[0]
        def calc(z, derivatives=False):
            q1,q2,start,delta=z;onsets=[start,start+delta]
            cols=[];slopes=[]
            for onset,kind,shape in zip(onsets,pair,shapes):
                p=bank.calibration['classes'][bank.names[kind]]
                sq=np.clip(max(shape,1.),*p['charge_range_pe'])
                slow=np.clip(p['slow_intercept']+p['slow_log_slope']*np.log(sq/1000.),.02,.98)
                x=t-onset-p['delay_ns']
                cols.append(bank.column(t,onset,kind,shape))
                if derivatives:
                    fast,tail=bank.basis
                    slopes.append(-scale*((1-slow)*np.interp(x,fast.grid,fast.slope,left=0.,right=0.)*fast.unit_area_peak+
                        slow*np.interp(x,tail.grid,tail.slope,left=0.,right=0.)*tail.unit_area_peak))
            if derivatives:
                return W.apply(np.column_stack([cols[0],cols[1],q1*slopes[0]+q2*slopes[1],q2*slopes[1]]))
            return W.apply(cols[0]*q1+cols[1]*q2)-target
        bounds=([0.,0.,lower,8.],[np.inf,np.inf,min(upper-8.,times[0]+16.),max(8.01,min(upper-lower,gap+40.))])
        if bounds[1][2]<=bounds[0][2]:continue
        z=np.r_[q,times[0],max(8.01,gap)]
        z=np.maximum(z,np.array(bounds[0])+1e-8);z=np.minimum(z,np.array(bounds[1])-1e-8)
        opt=least_squares(calc,z,jac=lambda z:calc(z,True),bounds=bounds,
            x_scale=[1000.,1000.,5.,5.],max_nfev=35,ftol=1e-8,xtol=1e-7,gtol=1e-7)
        tt=np.array([opt.x[2],opt.x[2]+opt.x[3]])
        trials+=1
        if tt[1]>upper:continue
        cost=float(2*opt.cost)
        if best is None or cost<best['cost']:
            best=dict(cost=cost,times=tt,types=np.array(pair),charges=opt.x[:2],
                      matrix=W.apply(bank.matrix(t,tt,pair,shapes)),success=bool(opt.success))
    return best,trials


def fit_model_selection(time_ns,voltage_mV,calibration,config=ModelSelectionConfig()):
    """Blind reconstruction plus cost-penalized residual proposals."""
    initial=_initial_fit(time_ns,voltage_mV,calibration,RegularizedConfig(
        selection_penalty=config.selection_penalty,fit_stride=config.fit_stride))
    if initial['status']!='fitted' or not initial['all_components']['times_ns']:return initial
    bank=CalibratedBank(calibration);full_t=np.asarray(time_ns,float);full_y=np.asarray(voltage_mV,float)
    t=full_t[::config.fit_stride];values=full_y[::config.fit_stride]-calibration['baseline_mV']
    components=initial['all_components'];times=np.array(components['times_ns']);q=np.array(components['charge_pe'])
    types=np.array([0 if s=='alpha' else 1 for s in components['particle_name']]);shapes=np.array(components['shape_charge_pe'])
    W=FrozenWhitening(initial['model_mV'][::config.fit_stride]-calibration['baseline_mV'],
                      t[1]-t[0],calibration['noise'],smooth_ns=0.)
    target=W.apply(values);H=W.apply(bank.matrix(t,times,types,shapes));residual=target-H@q
    costs=[float(residual@residual)];pair_trials=0;additions=[]
    cfg=FitConfig(min_charge_pe=80.,area_per_pe_mV_ns=bank.area)
    for iteration in range(config.joint_passes):
        # Near-neighbour choice is based only on reconstructed events.
        pairs=[(a,a+1) for a in range(len(times)-1) if times[a+1]-times[a]<=config.pair_distance_ns and
               max(q[a:a+2])>=config.pair_min_charge_pe]
        for a,b in pairs:
            ix=[a,b];local= residual+H[:,ix]@q[ix]
            lower=max(0.,times[a]-16.,times[a-1]+8. if a else 0.)
            upper=min(t[-1],times[b]+40.,times[b+1]-8. if b+1<len(times) else t[-1])
            best,trials=_joint_pair(t,local,W,bank,times[ix],types[ix],q[ix],shapes[ix],lower,upper)
            pair_trials+=trials
            if best is not None and best['cost']<residual@residual:
                times[ix]=best['times'];types[ix]=best['types'];q[ix]=best['charges'];H[:,ix]=best['matrix']
                residual=local-H[:,ix]@q[ix]
        q=reduced_nnls(H,target,maxiter=max(100,20*len(times)))[0];residual=target-H@q
        # Model shape may update with charge only when the same fixed weighted
        # objective decreases. Failed updates leave the previous model intact.
        trial_shapes=np.maximum(q,1.);trial_H=W.apply(bank.matrix(t,times,types,trial_shapes))
        trial_q=reduced_nnls(trial_H,target,maxiter=max(100,20*len(times)))[0];trial_r=target-trial_H@trial_q
        if trial_r@trial_r < residual@residual:H,q,shapes,residual=trial_H,trial_q,trial_shapes,trial_r
        costs.append(float(residual@residual))
        if len(additions)>=config.max_additions:continue
        full_model=bank.matrix(full_t,times,types,shapes)@q
        proposals=discover(full_t,full_y-calibration['baseline_mV']-full_model,bank.discovery_bank(),cfg)
        proposals=[x for x in proposals if np.min(abs(times-x))>=8.]
        ranked=[]
        for x in proposals:
            score=0.
            for k in (0,1):
                h=W.apply(bank.column(t,x,k,500.));score=max(score,max(float(h@residual),0.)**2/max(float(h@h),1e-20))
            ranked.append((score,x))
        for _,x in sorted(ranked,reverse=True)[:config.proposal_budget]:
            if len(additions)>=config.max_additions:break
            if np.min(abs(times-x))<8.:continue
            j=int(np.argmin(abs(times-x)))
            if abs(times[j]-x)>config.pair_distance_ns:continue
            if np.any((times>min(x,times[j]))&(times<max(x,times[j]))):continue
            tt=np.sort([x,times[j]]);existing=int(times[j]>x)
            qs=np.zeros(2);qs[existing]=q[j];ss=np.full(2,500.);ss[existing]=shapes[j]
            local=residual+H[:,j]*q[j]
            left=times[times<tt[0]];right=times[times>tt[1]]
            lo=max(0.,tt[0]-16.,left[-1]+8. if len(left) else 0.)
            hi=min(t[-1],tt[1]+40.,right[0]-8. if len(right) else t[-1])
            best,trials=_joint_pair(t,local,W,bank,tt,[0,0],qs,ss,lo,hi);pair_trials+=trials
            if best is None:continue
            gain=float(residual@residual-best['cost'])
            if gain<=config.addition_cost or min(best['charges'])<50.:continue
            times=np.r_[np.delete(times,j),best['times']];types=np.r_[np.delete(types,j),best['types']]
            q=np.r_[np.delete(q,j),best['charges']];shapes=np.r_[np.delete(shapes,j),ss]
            order=np.argsort(times);times,types,q,shapes=times[order],types[order],q[order],shapes[order]
            H=W.apply(bank.matrix(t,times,types,shapes));residual=target-H@q
            additions.append(dict(proposal_time_ns=float(x),weighted_improvement=gain,pass_index=iteration))
        costs.append(float(residual@residual))
    active=q>1e-6;times,types,q,shapes=times[active],types[active],q[active],shapes[active]
    model=calibration['baseline_mV']+bank.matrix(full_t,times,types,shapes)@q
    report=q>=500.
    return dict(status='fitted',times_ns=times[report].tolist(),particle_name=bank.names[types[report]].tolist(),
        charge_pe=q[report].tolist(),area_mV_ns=(q[report]*bank.area).tolist(),model_mV=model,
        all_components=dict(times_ns=times.tolist(),particle_name=bank.names[types].tolist(),charge_pe=q.tolist(),shape_charge_pe=shapes.tolist()),
        rms_mV=float(np.sqrt(np.mean((full_y-model)**2))),selection=initial['selection'],
        initial_reported_count=len(initial['times_ns']),initial_components=initial['selected_components'],
        selected_components=len(times),nuisance_component_count=int(np.sum(~report)),
        model_includes_nuisance_components=True,weighted_cost_history=costs,pair_trials=pair_trials,
        residual_additions=additions,config=asdict(config),fit_truth_inputs=[],
        objective='Fixed pilot-mean AR(1), strict L1 initialization, joint pairs and residual proposals with fixed cost penalty')

def _initial_fit(time_ns, voltage_mV, calibration, config=RegularizedConfig(), *, weighted=True):
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
    sparse,selection=strict_positive_lasso(design,target,config.selection_penalty)
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

