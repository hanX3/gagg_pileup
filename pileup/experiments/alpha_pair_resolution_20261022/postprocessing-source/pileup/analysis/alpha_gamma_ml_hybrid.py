"""Neural candidates and soft priors with a voltage-only joint physical refit.

Unknown classes remain discrete hypotheses; NN labels are never fixed. Candidate
support includes the existing blind fit and new neural events. A residual pass
can add events absent from both. The soft-prior strength is validation-selected.
"""
from dataclasses import dataclass,asdict
import numpy as np
from scipy.optimize import least_squares,minimize_scalar
from .alpha_gamma_regularized_fit import CalibratedBank,FrozenWhitening
from .alpha_gamma_nnls import reduced_nnls
from .alpha_gamma_pileup_fit import discover
from .alpha_pileup_fit import FitConfig


@dataclass(frozen=True)
class MixedHybridConfig:
    strength: float=10.
    passes: int=2
    time_range_ns: float=8.
    min_separation_ns: float=4.
    pair_distance_ns: float=100.
    max_additions: int=8
    addition_cost: float=9.
    max_components: int=400
    fit_stride: int=2


def _prior_qfit(H,target,weight,prior):
    if np.any(weight):
        return reduced_nnls(np.vstack([H,np.diag(np.sqrt(weight))]),
                            np.r_[target,np.sqrt(weight)*prior],maxiter=max(100,20*H.shape[1]))[0]
    return reduced_nnls(H,target,maxiter=max(100,20*H.shape[1]))[0]


def fit_hybrid(time_ns,voltage_mV,calibration,neural,dense,baseline,config=MixedHybridConfig()):
    if config.strength<0:raise ValueError('Nonnegative prior strength required')
    bank=CalibratedBank(calibration);full_t=np.asarray(time_ns,float);full_y=np.asarray(voltage_mV,float)
    initial=baseline.get('all_components',baseline);nn=neural['all_components']
    entries=[]
    # Neural coordinates are used where the two candidates refer to one onset;
    # nearby distinct candidates retain a declared >=4 ns output capacity.
    for j,(t,p,q) in enumerate(zip(nn['times_ns'],nn['particle_name'],nn['charge_pe'])):
        entries.append((2.,t,0 if p=='alpha' else 1,q,j))
    for t,p,q in zip(initial['times_ns'],initial['particle_name'],initial['charge_pe']):
        entries.append((1.,t,0 if p=='alpha' else 1,q,-1))
    kept=[]
    for row in sorted(entries,reverse=True):
        if row[3]<40:continue
        if all(abs(row[1]-x[1])>=config.min_separation_ns for x in kept):kept.append(row)
    kept.sort(key=lambda r:r[1])
    if len(kept)>config.max_components:
        return dict(status='component_limit',times_ns=[],charge_pe=[],particle_name=[],
                    model_mV=np.full_like(full_y,calibration['baseline_mV']),fit_truth_inputs=[])
    if not kept:
        return dict(status='fitted',times_ns=[],charge_pe=[],particle_name=[],all_components=dict(times_ns=[],charge_pe=[],particle_name=[]),
                    model_mV=np.full_like(full_y,calibration['baseline_mV']),fit_truth_inputs=[],config=asdict(config))
    times=np.array([r[1] for r in kept]);types=np.array([r[2] for r in kept]);q=np.array([r[3] for r in kept]);starts=times.copy()
    shapes=np.maximum(q,1.);prior=np.zeros(len(q));weight=np.zeros(len(q));class_cost=np.zeros((len(q),2))
    nt=np.array(nn['times_ns'])
    for j,t in enumerate(times):
        if not len(nt):continue
        m=int(np.argmin(abs(nt-t)))
        if abs(nt[m]-t)>8.:continue
        prior[j]=nn['charge_pe'][m]
        weight[j]=config.strength/(.35*max(500.,prior[j]))**2
        index=int(np.clip(np.rint(nt[m]/2.),0,len(full_t)-1))
        scores=np.maximum(dense['event_score'][:,index],1e-4);prob=scores/scores.sum()
        class_cost[j]=-2*config.strength*np.log(prob)
    t=full_t[::config.fit_stride];values=full_y[::config.fit_stride]-calibration['baseline_mV']
    pilot=np.asarray(baseline.get('model_mV',full_y))[::config.fit_stride]-calibration['baseline_mV']
    W=FrozenWhitening(pilot,t[1]-t[0],calibration['noise'],smooth_ns=0.)
    target=W.apply(values);H=W.apply(bank.matrix(t,times,types,shapes))
    q=_prior_qfit(H,target,weight,prior);residual=target-H@q
    def total_cost(r,qq,kk):return float(r@r+np.sum(weight*(qq-prior)**2)+class_cost[np.arange(len(kk)),kk].sum())
    costs=[total_cost(residual,q,types)];pair_trials=0;stopping=[]
    for iteration in range(config.passes):
        for j in range(len(times)):
            local=residual+H[:,j]*q[j]
            lo=max(t[0],starts[j]-config.time_range_ns,times[j-1]+config.min_separation_ns if j else t[0])
            hi=min(t[-1],starts[j]+config.time_range_ns,times[j+1]-config.min_separation_ns if j+1<len(times) else t[-1])
            if hi<=lo:continue
            def evaluate(x,k):
                h=W.apply(bank.column(t,x,k,shapes[j]));den=float(h@h)+weight[j]
                amp=max(0.,(float(h@local)+weight[j]*prior[j])/max(den,1e-20));r=local-h*amp
                return float(r@r+weight[j]*(amp-prior[j])**2+class_cost[j,k]),amp,h
            best=(*evaluate(times[j],types[j])[:2],times[j],types[j],H[:,j])
            for k in (0,1):
                opt=minimize_scalar(lambda x:evaluate(x,k)[0],bounds=(lo,hi),method='bounded',options={'xatol':.15,'maxiter':18})
                cost,amp,h=evaluate(opt.x,k)
                if cost<best[0]:best=(cost,amp,float(opt.x),k,h)
            _,q[j],times[j],types[j],H[:,j]=best;residual=local-H[:,j]*q[j]
        # Jointly float two amplitudes/times for every close pair and all classes.
        for a in range(len(times)-1):
            b=a+1;gap=times[b]-times[a]
            if gap>config.pair_distance_ns:continue
            ix=np.array([a,b]);local=residual+H[:,ix]@q[ix]
            lo=max(t[0],starts[a]-config.time_range_ns,times[a-1]+4. if a else t[0])
            hi=min(t[-1],starts[b]+config.time_range_ns,times[b+1]-4. if b+1<len(times) else t[-1])
            start_upper=min(starts[a]+config.time_range_ns,hi-4.01)
            if start_upper<=lo or hi-lo<=4.01:continue
            current=float(residual@residual+np.sum(weight[ix]*(q[ix]-prior[ix])**2)+class_cost[ix,types[ix]].sum())
            best=None
            for pair in ((0,0),(0,1),(1,0),(1,1)):
                def fun(z):
                    onsets=[z[2],z[2]+z[3]];h=W.apply(bank.matrix(t,onsets,pair,shapes[ix]))
                    return np.r_[h@z[:2]-local,np.sqrt(weight[ix])*(z[:2]-prior[ix])]
                def jac(z):
                    onsets=[z[2],z[2]+z[3]];h=W.apply(bank.matrix(t,onsets,pair,shapes[ix]))
                    derivatives=[]
                    for x,k,s in zip(onsets,pair,shapes[ix]):
                        derivatives.append(W.apply((bank.column(t,x+.1,k,s)-bank.column(t,x-.1,k,s))/.2))
                    matrix=np.column_stack([h,z[0]*derivatives[0]+z[1]*derivatives[1],z[1]*derivatives[1]])
                    return np.vstack([matrix,np.c_[np.diag(np.sqrt(weight[ix])),np.zeros((2,2))]])
                lower=np.array([0.,0.,lo,4.]);upper=np.array([np.inf,np.inf,start_upper,min(gap+16.,hi-lo)])
                z=np.minimum(np.maximum(np.r_[q[ix],times[a],gap],lower+1e-7),upper-1e-7)
                opt=least_squares(fun,z,jac=jac,bounds=(lower,upper),x_scale=[1000.,1000.,5.,5.],
                    max_nfev=20,ftol=1e-7,xtol=1e-6,gtol=1e-6)
                pair_trials+=1;stopping.append(int(opt.status))
                if opt.x[2]+opt.x[3]>hi:continue
                cost=float(2*opt.cost+class_cost[ix,list(pair)].sum())
                if cost<current:
                    current=cost;best=(pair,opt.x)
            if best is not None:
                pair,z=best;times[ix]=[z[2],z[2]+z[3]];types[ix]=pair;q[ix]=z[:2]
                H[:,ix]=W.apply(bank.matrix(t,times[ix],types[ix],shapes[ix]));residual=local-H[:,ix]@q[ix]
        q=_prior_qfit(H,target,weight,prior);residual=target-H@q;costs.append(total_cost(residual,q,types))
    additions=[]
    model=bank.matrix(full_t,times,types,shapes)@q
    cfg=FitConfig(area_per_pe_mV_ns=bank.area,min_charge_pe=100.,min_separation_ns=4.)
    proposals=discover(full_t,full_y-calibration['baseline_mV']-model,bank.discovery_bank(),cfg)
    ranked=[]
    for x in proposals:
        if min(abs(times-x))<8.:continue
        for k in (0,1):
            h=W.apply(bank.column(t,x,k,500.));amp=max(0.,float(h@residual)/max(float(h@h),1e-20))
            gain=float(residual@residual-(residual-h*amp)@(residual-h*amp))
            ranked.append((gain,float(x),k,amp,h))
    for _,x,k,_,h in sorted(ranked,key=lambda r:-r[0]):
        if len(additions)>=config.max_additions or min(abs(times-x))<8.:continue
        amp=max(0.,float(h@residual)/max(float(h@h),1e-20));rr=residual-h*amp
        gain=float(residual@residual-rr@rr)
        if gain<=config.addition_cost or amp<100.:continue
        times=np.r_[times,x];types=np.r_[types,k];q=np.r_[q,amp];shapes=np.r_[shapes,500.]
        H=np.column_stack([H,h]);prior=np.r_[prior,0.];weight=np.r_[weight,0.];class_cost=np.vstack([class_cost,[0.,0.]])
        residual=rr;additions.append(dict(time_ns=x,weighted_improvement=gain))
    q=_prior_qfit(H,target,weight,prior);residual=target-H@q;costs.append(total_cost(residual,q,types))
    order=np.argsort(times);times,types,q,shapes=times[order],types[order],q[order],shapes[order]
    model=bank.matrix(full_t,times,types,shapes)@q+calibration['baseline_mV'];active=q>1e-6;report=q>=500.
    return dict(status='fitted',times_ns=times[report].tolist(),particle_name=bank.names[types[report]].tolist(),charge_pe=q[report].tolist(),
        model_mV=model,all_components=dict(times_ns=times[active].tolist(),particle_name=bank.names[types[active]].tolist(),
        charge_pe=q[active].tolist(),shape_charge_pe=shapes[active].tolist()),
        config=asdict(config),weighted_cost_history=costs,pair_trials=pair_trials,optimizer_status_counts={str(s):stopping.count(s) for s in sorted(set(stopping))},
        residual_additions=additions,initial_neural_count=len(nn['times_ns']),initial_fit_count=len(initial['times_ns']),
        rms_mV=float(np.sqrt(np.mean((full_y-model)**2))),fit_truth_inputs=[],
        objective='Fixed whitened physical residual + validation-selected soft neural charge/class priors; labels remain variable')
