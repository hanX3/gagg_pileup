"""Frozen observation-only inference and separately invoked oracle diagnostics."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import multiprocessing as mp
import time

import numpy as np
from scipy.optimize import nnls

from .alpha_pair_fitting import PairFitter
from .alpha_pair_resolution_data import (
    EXPERIMENT, OUTPUT, MODEL_NAMES, GRID, settings, load, save, digest, check_freeze)

_FITTER = None


def _init_fitter():
    global _FITTER
    _FITTER = PairFitter(settings())


def strip_fit(fit):
    return {k:v for k,v in fit.items() if k!='model_mV'}


def _candidates(task):
    index, voltage = task
    start = time.monotonic()
    candidates = _FITTER.candidates(voltage)
    return dict(index=index, candidates=[strip_fit(x) for x in candidates],
                seconds=time.monotonic()-start)


def _refine(task):
    index, voltage, initial = task
    start = time.monotonic()
    result = _FITTER.refine(voltage, initial, radius=12.)
    return dict(index=index, fit=strip_fit(result), seconds=time.monotonic()-start)


def _load_own_model(name):
    import torch
    from .alpha_pair_ordered import OrderedNet
    folder=EXPERIMENT/'models'/name
    frozen=load(folder/'frozen.json')
    assert digest(folder/'model.pt')==frozen['weights_sha256']
    assert digest(folder/'config.json')==frozen['config_sha256']
    model=OrderedNet(settings(), head='physical')
    model.load_state_dict(torch.load(folder/'model.pt',map_location='cpu',weights_only=True))
    return model.eval(), load(folder/'selection.json')['threshold']


def predict(name, workers=4):
    """This entry opens observations, calibration and frozen weights only."""
    check_freeze()
    assert 1<=workers<=4
    observation=OUTPUT/'assets'/(name+'.waveforms.npz')
    with np.load(observation, allow_pickle=False) as f:
        assert f.files==['voltage_mV'], 'Observation file must contain no labels.'
        voltage=f['voltage_mV'].astype(np.float64)
    assert voltage.ndim==2 and voltage.shape[1]==1024 and np.isfinite(voltage).all()
    destination=OUTPUT/'predictions'/name
    destination.mkdir(parents=True,exist_ok=True)
    header=dict(observation=str(observation), waveform_sha256=digest(observation),
                freeze_sha256=digest(EXPERIMENT/'freeze.json'), truth_inputs=[],
                events=len(voltage), models=list(MODEL_NAMES), workers=workers,
                all_K_hypotheses='0/1/2 generated from voltage before any truth is opened')
    if (destination/'inputs.json').exists():
        assert load(destination/'inputs.json')==header
    else: save(destination/'inputs.json',header)
    path=destination/'candidates.jsonl'
    done={}
    if path.exists():
        for line in path.read_text().splitlines():
            row=json.loads(line);assert row['index'] not in done
            done[row['index']]=row
    with ProcessPoolExecutor(max_workers=workers,initializer=_init_fitter,
                             mp_context=mp.get_context('spawn')) as pool:
        tasks={pool.submit(_candidates,(i,y)):i for i,y in enumerate(voltage) if i not in done}
        with path.open('a') as stream:
            for future in as_completed(tasks):
                row=future.result();done[row['index']]=row
                stream.write(json.dumps(row,allow_nan=False)+'\n')
                if len(done)%100==0:
                    stream.flush();print('FIT',name,len(done),'/',len(voltage),flush=True)
        assert len(done)==len(voltage)
        # Import neural code only after the independent fitting workers start.
        import torch
        from .alpha_pair_ordered import predict as neural_predict
        torch.set_num_threads(1)
        for model_name in MODEL_NAMES:
            path=destination/(model_name+'.json')
            if path.exists():
                old=load(path)
                assert old['waveform_sha256']==header['waveform_sha256']
                continue
            model,threshold=_load_own_model(model_name)
            fits=[]; arrays=[]; start=time.monotonic()
            for first in range(0,len(voltage),64):
                part,raw=neural_predict(model,voltage[first:first+64].astype(np.float32),threshold)
                fits.extend(part)
                arrays.append({k:raw[k] for k in ('count_logits','raw_charge_pe','raw_times_ns','area_fraction')})
            neural_seconds=time.monotonic()-start
            np.savez_compressed(destination/(model_name+'.npz'),
                                **{k:np.concatenate([a[k] for a in arrays]) for k in arrays[0]})
            refined=[None]*len(fits)
            futures={pool.submit(_refine,(i,voltage[i],f['times_ns'])):i for i,f in enumerate(fits)}
            for future in as_completed(futures):
                row=future.result();refined[row['index']]=row
            save(path,dict(fits=fits,refined=refined,threshold=threshold,
                           neural_seconds=neural_seconds,truth_inputs=[],
                           waveform_sha256=header['waveform_sha256'],
                           weights_sha256=digest(EXPERIMENT/'models'/model_name/'model.pt')))
            print('NEURAL',name,model_name,'done',flush=True)
    save(destination/'complete.json',dict(events=len(voltage),truth_inputs=[],
         candidate_file_sha256=digest(destination/'candidates.jsonl'),
         waveform_sha256=header['waveform_sha256']))


def whiten(x, sigma, rho):
    x=np.asarray(x,float)
    result=x.copy()
    if rho:
        result[1:]=(x[1:]-rho*x[:-1])/np.sqrt(1.-rho*rho)
    return result/sigma


def diagnostics(name):
    """Known true K/t and electronic-noise covariance: never a blind result."""
    check_freeze()
    assert (OUTPUT/'predictions'/name/'complete.json').exists(), 'Blind inference must finish first.'
    path=OUTPUT/'diagnostics'/(name+'.json')
    assert not path.exists()
    with np.load(OUTPUT/'assets'/(name+'.waveforms.npz')) as f:
        voltage=f['voltage_mV'].astype(float)
    with np.load(OUTPUT/'assets'/(name+'.truth.npz')) as f:
        q=f['charge_pe']; times=f['times_ns']; active=f['active']
    cases=load(OUTPUT/'assets'/(name+'.json'))['cases']
    fitter=PairFitter(settings()); rows=[]
    for i,(y,c) in enumerate(zip(voltage,cases)):
        ids=np.flatnonzero(active[i]); tt=times[i,ids]; qq=q[i,ids]
        row=dict(index=i,true_K=len(ids),truth_inputs=['K','times','Q_for_Fisher_only','noise_covariance'],
                 validity='Exact fixed-mean/Gaussian-noise model only. Shared/linear physics has model fluctuations not represented by C.')
        if not len(ids):
            row.update(known_times_fit=dict(times_ns=[],charge_pe=[]),rank=0)
            rows.append(row);continue
        H=fitter.columns(GRID,tt)
        sigma=float(c.get('noise_rms',.1));rho_noise=float(c.get('ar1',0.))
        A=whiten(H,sigma,rho_noise);b=whiten(y-fitter.baseline,sigma,rho_noise)
        rank=int(np.linalg.matrix_rank(A))
        row.update(rank=rank, condition=float(np.linalg.cond(A)) if rank==len(ids) else None)
        if rank==len(ids):
            fitted=nnls(A,b,maxiter=100)[0]
            covariance=np.linalg.inv(A.T@A)
            row['known_times_fit']=dict(times_ns=tt.tolist(),charge_pe=fitted.tolist())
            row['sigma_q_known_time_pe']=np.sqrt(np.maximum(0.,np.diag(covariance))).tolist()
        else:
            row['known_times_fit']=None
        if len(ids)==2:
            correlation=float(A[:,0]@A[:,1]/np.sqrt((A[:,0]@A[:,0])*(A[:,1]@A[:,1])))
            row['template_correlation']=correlation
            row['variance_inflation']=float(1./(1.-correlation**2)) if 1.-correlation**2>1e-14 else None
            derivative=np.interp(GRID[:,None]-tt[None],fitter.grid,fitter.slope,left=0.,right=0.)
            # Relative-Q coordinates and ns time coordinates avoid mixing pe
            # and ns scales when testing the Fisher matrix rank.
            J=np.column_stack((H*qq[None],-derivative*qq[None]))
            J=whiten(J,sigma,rho_noise)
            scales=np.linalg.norm(J,axis=0)
            normalized=J/np.maximum(scales,1e-30)[None]
            gram=normalized.T@normalized
            if np.linalg.matrix_rank(gram)==4:
                covariance=np.linalg.inv(gram)/scales[:,None]/scales[None]
                diagonal=np.diag(covariance)
                if np.all(diagonal>=0):
                    row['sigma_relative_q_unknown_times']=np.sqrt(diagonal[:2]).tolist()
                    row['sigma_time_unknown_times_ns']=np.sqrt(diagonal[2:]).tolist()
            row['joint_information_rank']=int(np.linalg.matrix_rank(gram))
        rows.append(row)
    save(path,dict(rows=rows,oracle=True,truth_inputs='Explicitly listed per row',
                   waveform_sha256=digest(OUTPUT/'assets'/(name+'.waveforms.npz')),
                   definition='Known true times: nonnegative GLS. Fisher diagnostics assume exact mean response and electronic noise only.'))
    print('DIAGNOSTICS',name,len(rows),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=('predict','diagnostics'))
    parser.add_argument('--dataset',required=True,choices=('mean','shared','linear','shared-noise'))
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    if args.stage=='predict':predict(args.dataset,args.workers)
    else:diagnostics(args.dataset)


if __name__=='__main__':
    main()
