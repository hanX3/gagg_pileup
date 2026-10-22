"""Blind unknown-count, unknown-type alpha/gamma voltage template fitting.

Only the sample grid, total voltage and frozen calibration enter reconstruction.
The separate scoring function reads labels afterwards. This is a project
baseline with discrete template selection, not a literature algorithm replica.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.optimize import least_squares, minimize_scalar, nnls

from .alpha_pileup_fit import AlphaTemplate, FitConfig, candidate_peaks, score


@dataclass(frozen=True)
class MixedFitConfig(FitConfig):
    time_refine_ns: float = 20.
    max_passes: int = 6
    residual_rounds: int = 1


class TemplateBank:
    def __init__(self, readout, *, alpha_slow=.7414997984399583, gamma_slow=.20473191192272483):
        self.templates=[AlphaTemplate(readout,alpha_slow),AlphaTemplate(readout,gamma_slow)]
        self.names=np.array(['alpha','gamma'])

    def matrix(self,t,times,types,config):
        return np.column_stack([self.templates[k](t-x)*self.templates[k].unit_area_peak*
                                config.area_per_pe_mV_ns for x,k in zip(times,types)]) if len(times) else np.empty((len(t),0))

    def metadata(self):
        return {name:h.metadata() for name,h in zip(self.names,self.templates)}


def discover(t,values,bank,config):
    candidates=[]
    for template in bank.templates:
        times,inverse=candidate_peaks(t,values,template,config)
        for time in times:
            index=int(round(time/(t[1]-t[0])))
            q=inverse[index]/(template.unit_area_peak*config.area_per_pe_mV_ns)
            candidates.append((time,q))
    # Inverting two templates can find the same onset twice. Keep the strongest
    # charge-equivalent discovery, then let floating-time fits select the type.
    kept=[]
    for time,q in sorted(candidates,key=lambda v:-v[1]):
        if all(abs(time-other)>=config.min_separation_ns for other in kept):kept.append(time)
    return np.sort(kept)


def refit_mixed(t,values,initial_times,bank,config):
    times=np.sort(np.asarray(initial_times,float)); starts=times.copy();n=len(times)
    if not n:return times,np.array([],int),np.array([]),np.zeros_like(values),[]
    # A relaxed two-template NNLS initialization lets overlapping neighbours
    # exchange charge before the discrete class decision. A local template
    # swap at fixed neighbour amplitudes can otherwise trap even exact doubles.
    relaxed=bank.matrix(t,np.repeat(times,2),np.tile([0,1],n),config)
    dual_q=nnls(relaxed,values,maxiter=max(200,20*n))[0].reshape(n,2)
    types=np.argmax(dual_q,axis=1)
    H=bank.matrix(t,times,types,config);q=nnls(H,values,maxiter=max(100,10*n))[0]
    residual=values-H@q;costs=[float(residual@residual)]
    for _ in range(config.max_passes):
        relaxed=bank.matrix(t,np.repeat(times,2),np.tile([0,1],n),config)
        proposed=np.argmax(nnls(relaxed,values,maxiter=max(200,20*n))[0].reshape(n,2),axis=1)
        if np.any(proposed!=types):
            trial_H=bank.matrix(t,times,proposed,config)
            trial_q=nnls(trial_H,values,maxiter=max(100,10*n))[0]
            trial_residual=values-trial_H@trial_q
            if trial_residual@trial_residual < residual@residual:
                types,H,q,residual=proposed,trial_H,trial_q,trial_residual
        for j in range(n):
            mask=(t>=starts[j]-config.time_refine_ns-4)&(t<=starts[j]+6000+config.time_refine_ns)
            tx=t[mask];target=residual[mask]+H[mask,j]*q[j]
            def evaluate(x,k):
                template=bank.templates[k]
                v=template(tx-x)*template.unit_area_peak*config.area_per_pe_mV_ns
                denom=float(v@v);amp=max(0.,float(v@target)/denom) if denom else 0.
                error=target-amp*v
                return float(error@error),amp
            old=evaluate(times[j],types[j]);best=(old[0],times[j],types[j],old[1])
            for k in range(2):
                result=minimize_scalar(lambda x:evaluate(x,k)[0],
                    bounds=(max(0.,starts[j]-config.time_refine_ns),min(t[-1],starts[j]+config.time_refine_ns)),
                    method='bounded',options={'xatol':.02,'maxiter':35})
                if result.fun<best[0]:best=(result.fun,result.x,k,evaluate(result.x,k)[1])
            if best[0]<old[0]:
                residual+=H[:,j]*q[j];times[j],types[j],q[j]=best[1:]
                H[:,j]=bank.matrix(t,[times[j]],[types[j]],config)[:,0];residual-=H[:,j]*q[j]
        q=nnls(H,values,maxiter=max(100,10*n))[0];residual=values-H@q
        # A close alpha/gamma pair can be stuck with exchanged classes. Test
        # the four pair assignments while refitting both amplitudes together.
        # Neighbour selection uses fitted times, never true pairing or types.
        order=np.argsort(times)
        for a,b in zip(order[:-1],order[1:]):
            if times[b]-times[a]>500.:continue
            index=[a,b];target=residual+H[:,index]@q[index]
            best_cost=float(residual@residual);best=None
            for pair in ((0,0),(0,1),(1,0),(1,1)):
                if np.array_equal(types[index],pair):continue
                trial_H=bank.matrix(t,times[index],pair,config)
                trial_q=nnls(trial_H,target,maxiter=100)[0];error=target-trial_H@trial_q
                cost=float(error@error)
                if cost<best_cost:
                    best_cost=cost;best=(pair,trial_H,trial_q,error)
            if best is not None:
                pair,trial_H,trial_q,residual=best
                types[index]=pair;H[:,index]=trial_H;q[index]=trial_q
        q=nnls(H,values,maxiter=max(100,10*n))[0];residual=values-H@q;costs.append(float(residual@residual))
        if costs[-2]-costs[-1]<=1e-5*max(1.,costs[-2]):break
    if n<=64:
        def fun(z):return bank.matrix(t,z[n:],types,config)@z[:n]-values
        def jac(z):
            H=bank.matrix(t,z[n:],types,config);derivatives=[]
            for x,k,amp in zip(z[n:],types,z[:n]):
                h=bank.templates[k]
                dh=np.interp(t-x,h.grid,h.slope,left=0.,right=0.)*h.unit_area_peak*config.area_per_pe_mV_ns
                derivatives.append(-dh*amp)
            return np.c_[H,np.column_stack(derivatives)]
        solution=least_squares(fun,np.r_[q,times],jac=jac,
            bounds=(np.r_[np.zeros(n),np.maximum(0.,starts-config.time_refine_ns)],
                    np.r_[np.full(n,np.inf),np.minimum(t[-1],starts+config.time_refine_ns)]),
            x_scale=np.r_[np.full(n,1000.),np.full(n,5.)],max_nfev=60,ftol=1e-8,xtol=1e-7,gtol=1e-7)
        if 2*solution.cost<costs[-1]:q,times=solution.x[:n],solution.x[n:];costs.append(float(2*solution.cost))
    # Below 100 ns, changing both classes and both times together is needed:
    # fixed-time label swaps may still prefer the wrong local minimum.
    H=bank.matrix(t,times,types,config);residual=values-H@q
    order=np.argsort(times)
    for a,b in zip(order[:-1],order[1:]):
        if times[b]-times[a]>100.:continue
        index=[a,b]
        mask=(t>=min(starts[index])-config.time_refine_ns-4)&(t<=max(starts[index])+6000+config.time_refine_ns)
        tx=t[mask];target=residual[mask]+H[mask][:,index]@q[index]
        old=float(residual[mask]@residual[mask]);best=None
        for pair in ((0,0),(0,1),(1,0),(1,1)):
            if np.array_equal(types[index],pair):continue
            pair=np.asarray(pair)
            initial_H=bank.matrix(tx,times[index],pair,config)
            initial_q=nnls(initial_H,target,maxiter=100)[0]
            def fun_pair(z):return bank.matrix(tx,z[2:],pair,config)@z[:2]-target
            def jac_pair(z):
                matrix=bank.matrix(tx,z[2:],pair,config);derivatives=[]
                for x,k,amp in zip(z[2:],pair,z[:2]):
                    h=bank.templates[k]
                    dh=np.interp(tx-x,h.grid,h.slope,left=0.,right=0.)*h.unit_area_peak*config.area_per_pe_mV_ns
                    derivatives.append(-dh*amp)
                return np.c_[matrix,np.column_stack(derivatives)]
            solution=least_squares(fun_pair,np.r_[initial_q,times[index]],jac=jac_pair,
                bounds=(np.r_[0.,0.,np.maximum(0.,starts[index]-config.time_refine_ns)],
                        np.r_[np.inf,np.inf,np.minimum(t[-1],starts[index]+config.time_refine_ns)]),
                x_scale=[1000.,1000.,5.,5.],max_nfev=50,ftol=1e-8,xtol=1e-7,gtol=1e-7)
            if 2*solution.cost<old:
                old=2*solution.cost;best=(pair.copy(),solution.x.copy())
        if best is not None:
            pair,z=best;residual+=H[:,index]@q[index]
            types[index]=pair;q[index]=z[:2];times[index]=z[2:]
            H[:,index]=bank.matrix(t,times[index],pair,config);residual-=H[:,index]@q[index]
    q=nnls(H,values,maxiter=max(100,10*n))[0];costs.append(float(np.sum((values-H@q)**2)))
    order=np.argsort(times)
    return times[order],types[order],q[order],bank.matrix(t,times,types,config)@q,costs


def fit_waveform(time_ns,voltage_mV,bank,config=MixedFitConfig()):
    t,y=np.asarray(time_ns,float),np.asarray(voltage_mV,float)
    if t.ndim!=1 or y.shape!=t.shape or len(t)<10 or not np.isfinite(y).all():raise ValueError('Invalid waveform')
    if not np.allclose(np.diff(t),t[1]-t[0],atol=1e-8,rtol=0):raise ValueError('Uniform samples required')
    if abs(t[0])>t[1]-t[0]:raise ValueError('Full record starting near zero required')
    values=y-config.baseline_mV;initial=discover(t,values,bank,config)
    def clean(initial):
        times,k,q,model,history=refit_mixed(t,values,initial,bank,config)
        while len(times)>1 and np.any(np.diff(times)<config.min_separation_ns):
            groups=np.split(np.arange(len(times)),np.flatnonzero(np.diff(times)>=config.min_separation_ns)+1)
            merged=[float(np.average(times[g],weights=np.maximum(q[g],1e-12))) for g in groups]
            times,k,q,model,last=refit_mixed(t,values,merged,bank,config);history+=last
        return times,k,q,model,history
    if len(initial)>config.max_candidates:
        return dict(status='candidate_limit',times_ns=[],particle_name=[],charge_pe=[],model_mV=np.full_like(y,config.baseline_mV))
    times,k,q,model,history=clean(initial)
    for _ in range(config.residual_rounds):
        more=discover(t,values-model,bank,config)
        more=[x for x in more if len(times)==0 or np.min(abs(times-x))>=config.min_separation_ns]
        if not more or len(more)+len(times)>config.max_candidates:break
        times,k,q,model,last=clean(np.r_[times,more]);history+=last
    times,k,q,model,last=clean(times[q>=config.min_charge_pe]);history+=last
    keep=q>=config.min_charge_pe;times,k,q=times[keep],k[keep],q[keep]
    model=bank.matrix(t,times,k,config)@q
    return dict(status='fitted',times_ns=times.tolist(),particle_name=bank.names[k].tolist(),charge_pe=q.tolist(),
        area_mV_ns=(q*config.area_per_pe_mV_ns).tolist(),model_mV=model+config.baseline_mV,
        candidates=len(initial),cost_history=history,rms_mV=float(np.sqrt(np.mean((values-model)**2))))


def score_mixed(fit,truth,gate_ns, *, tolerance_ns=10.,charge_threshold_pe=500.):
    result=score(fit,truth,gate_ns,tolerance_ns=tolerance_ns,charge_threshold_pe=charge_threshold_pe)
    names=np.asarray(truth['particle_name']);pred=np.asarray(fit['particle_name']);q=np.asarray(truth['charge_pe'])
    pairs=result['matched_pairs'];targets=result['target_indices'];by_truth={i:j for i,j in pairs}
    counts={};confusion={}
    for name in ('alpha','gamma'):
        indices=[i for i in targets if names[i]==name];mp=[(i,by_truth[i]) for i in indices if i in by_truth]
        correct=[(i,j) for i,j in mp if pred[j]==name]
        dq=[fit['charge_pe'][j]/q[i]-1 for i,j in mp]
        dt=[fit['times_ns'][j]-truth['first_active_ns'][i] for i,j in mp]
        recovered=sum(abs(fit['charge_pe'][j]/q[i]-1)<=.2 for i,j in correct)
        central=(truth['emission_ns']>=gate_ns[0])&(truth['emission_ns']<gate_ns[1])&(names==name)
        counts[name]=dict(target_count=len(indices),central_injected=int(central.sum()),
            central_nonzero_charge=int(np.sum(central&(q>0))),matched_count=len(mp),correct_type_count=len(correct),
            recovered_count=int(recovered),missed_count=len(indices)-len(mp),wrong_type_count=len(mp)-len(correct),
            relative_charge_errors=[float(v) for v in dq],time_errors_ns=[float(v) for v in dt])
        confusion[name]={p:int(sum(pred[j]==p for i,j in mp)) for p in ('alpha','gamma')}
        confusion[name]['missed']=len(indices)-len(mp)
    result['by_particle']=counts;result['confusion']=confusion
    result['correct_type_count']=sum(v['correct_type_count'] for v in counts.values())
    result['charge_and_type_recovered_count']=sum(v['recovered_count'] for v in counts.values())
    used={j for i,j in result['all_matched_pairs']}
    extras=[j for j,x in enumerate(fit['times_ns']) if gate_ns[0]<=x<gate_ns[1] and j not in used]
    result['unmatched_prediction_indices']=extras
    result['unmatched_predictions_by_class']={p:int(sum(pred[j]==p for j in extras)) for p in ('alpha','gamma')}
    result['success_definition']=f'Time within {tolerance_ns:g} ns, correct particle class, realized charge within 20%; Q>={charge_threshold_pe:g} pe truths'
    return result
