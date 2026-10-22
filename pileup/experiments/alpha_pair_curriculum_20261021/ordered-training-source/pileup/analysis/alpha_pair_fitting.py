"""Blind zero/one/two fixed-response fits for the local component experiment."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize,nnls
from scipy.signal import find_peaks

from .alpha_full_signal_study import make_template


@dataclass
class PairFitter:
    settings: dict

    def __post_init__(self):
        template=make_template(self.settings)
        self.grid=np.arange(-24.,6000.01,.25)
        self.kernel=template(self.grid)*self.settings['fit_config']['area_per_pe_mV_ns']*template.unit_area_peak
        self.slope=np.gradient(self.kernel,.25)
        self.area_per_pe=self.settings['fit_config']['area_per_pe_mV_ns']
        self.baseline=self.settings['fit_config']['baseline_mV']

    def columns(self,t,times):
        u=t[:,None]-np.asarray(times)[None]
        return np.interp(u,self.grid,self.kernel,left=0.,right=0.)

    def refine(self,voltage,initial,*,radius=None):
        """Variable projection: NNLS Q for every trial t, fixed double tau.

        Initial times must come from observed voltage or a blind neural output.
        The known-count variant is explicitly marked by its caller.
        """
        y=np.asarray(voltage,float)-self.baseline;t=np.arange(len(y))*2.+.25
        initial=np.asarray(initial,float)
        if not len(initial):
            return dict(times_ns=[],charge_pe=[],cost=float(y@y),model_mV=np.zeros_like(y))
        def objective(times):
            H=self.columns(t,times)
            q=nnls(H,y,maxiter=100)[0]
            r=H@q-y
            slope=np.interp(t[:,None]-times[None],self.grid,self.slope,left=0.,right=0.)
            gradient=-2.*np.sum(slope*q[None]*r[:,None],axis=0)
            return float(r@r),gradient
        bounds=[(0.,float(t[-1]))]*len(initial) if radius is None else [(max(0.,v-radius),min(float(t[-1]),v+radius)) for v in initial]
        fit=minimize(objective,initial,method='L-BFGS-B',jac=True,bounds=bounds,
                     options=dict(maxiter=100,ftol=1e-12,gtol=1e-6,maxls=30))
        times=fit.x if fit.fun<=objective(initial)[0] else initial
        H=self.columns(t,times);q=nnls(H,y,maxiter=100)[0]
        order=np.argsort(times)
        return dict(times_ns=times[order].astype(float).tolist(),charge_pe=q[order].astype(float).tolist(),
                    cost=float(np.sum((H@q-y)**2)),model_mV=H@q,converged=bool(fit.success))

    def candidates(self,voltage):
        """Estimate model0/1/2 from the observation; no K/t/Q truth accepted."""
        y=np.asarray(voltage,float)-self.baseline;t=np.arange(len(y))*2.+.25
        zero=dict(times_ns=[],charge_pe=[],cost=float(y@y),model_mV=np.zeros_like(y))
        peak=float(y.max())
        if peak<=1.:
            return [zero]
        crossed=np.flatnonzero(y>max(.5,.03*peak))
        onset=max(0.,float(t[crossed[0]])-3.) if len(crossed) else float(t[y.argmax()])
        one=self.refine(voltage,[onset])
        derivative=gaussian_filter1d(np.gradient(y,2.),1.)
        peaks,_=find_peaks(derivative,distance=2,prominence=max(.05,.1*derivative.max()))
        peaks=peaks[np.argsort(derivative[peaks])[-4:]]
        starts=[[one['times_ns'][0],one['times_ns'][0]+gap] for gap in (4.,12.,40.,100.,200.)]
        starts += [[max(0.,float(t[i])-4.),max(0.,float(t[j])-4.)] for k,i in enumerate(peaks) for j in peaks[k+1:] if i!=j]
        fits=[self.refine(voltage,np.clip(s,0.,t[-1])) for s in starts]
        two=min(fits,key=lambda f:f['cost'])
        return [zero,one,two]

    @staticmethod
    def choose(candidates,penalty):
        norm=max(candidates[0]['cost'],1.)
        eligible=[c for c in candidates if all(q>=500. for q in c['charge_pe'])]
        chosen=min(eligible,key=lambda c:c['cost']/norm+penalty*len(c['times_ns']))
        return {k:chosen[k] for k in ('times_ns','charge_pe')}
