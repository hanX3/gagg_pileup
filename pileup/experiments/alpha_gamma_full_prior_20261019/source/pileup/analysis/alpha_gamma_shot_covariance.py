"""Fixed compound-shot covariance through the documented linear readout.

This is a project approximation: independent arrivals and a scalar excess
variance calibrated from isolated signals. Shared-cell correlations and
parameter-dependent likelihood terms are not claimed to be modeled exactly.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.linalg import cholesky_banded, solve_banded
from scipy.signal import lfilter
from .alpha_gamma_regularized_fit import CalibratedBank
from .alpha_rate_campaign import digest, save


def arrival_intensity(t, times, kinds, charges, shapes, calibration):
    """Expected charge-equivalent arrivals/ns from a voltage-derived pilot."""
    result=np.full(len(t),.002)  # Known configured 2 MHz dark primary rate.
    for onset,kind,charge,shape in zip(times,kinds,charges,shapes):
        p=calibration['classes'][('alpha','gamma')[int(kind)]]
        shape=np.clip(max(shape,1.),*p['charge_range_pe'])
        f=np.clip(p['slow_intercept']+p['slow_log_slope']*np.log(shape/1000.),.02,.98)
        u=t-onset-p['delay_ns'];positive=u>=0
        result[positive]+=charge*((1-f)*np.exp(-u[positive]/95.)/95.+f*np.exp(-u[positive]/351.)/351.)
    return result


def spe(t, readout, area):
    """Unit-charge voltage for Exp(decay)*Exp(rise)*Exp(lowpass), rise=lowpass."""
    r,d,low=(float(readout[k]) for k in ('rise_ns','decay_ns','lowpass_ns'))
    if not np.isclose(r,low) or d<=r:
        raise ValueError('Current shot model requires the frozen repeated-rise readout')
    t=np.asarray(t,float);u=np.maximum(t,0.)
    value=d/(d-r)**2*(np.exp(-u/d)-np.exp(-u/r))-u/(r*(d-r))*np.exp(-u/r)
    return np.where(t>=0.,np.maximum(value,0.),0.)*area


def covariance_bands(t, intensity, calibration, *, scale=1., support_ns=240., quadrature=12):
    """Finite SPE support yields a positive semidefinite shot covariance."""
    dt=float(t[1]-t[0]);n=len(t);length=min(n,int(np.ceil(support_ns/dt))+1)
    # Each sub-bin arrival is propagated through the same finite ADC aperture.
    offsets=(np.arange(quadrature)+.5)*dt/quadrature
    aperture=calibration['readout']['aperture_ns']
    grid=np.arange(length)[:,None]*dt+offsets[None,:]
    h=.5*(spe(grid-aperture/4,calibration['readout'],calibration['area_per_pe_mV_ns'])+
          spe(grid+aperture/4,calibration['readout'],calibration['area_per_pe_mV_ns']))
    covariance=np.zeros((length,n))
    for lag in range(length):
        kernel=np.mean(h[lag:]*h[:length-lag],axis=1)*dt*scale
        covariance[lag,:n-lag]=lfilter(kernel,[1.],intensity)[:n-lag]
    covariance[0]+=calibration['noise']['baseline_variance_mV2']
    return covariance


class ShotWhitening:
    def __init__(self,t,times,kinds,charges,shapes,calibration):
        intensity=arrival_intensity(np.asarray(t)-(t[1]-t[0])/2,times,kinds,charges,shapes,calibration)
        settings=calibration['physical_noise']
        covariance=covariance_bands(t,intensity,calibration,scale=settings['excess_variance_scale'],
                                   support_ns=settings['spe_support_ns'],quadrature=settings['quadrature'])
        self.lower=cholesky_banded(covariance,lower=True,check_finite=False)
        self.bandwidth=self.lower.shape[0]-1
        self.diagonal=covariance[0].copy()

    def apply(self,values):
        return solve_banded((self.bandwidth,0),self.lower,values,check_finite=False)


def calibrate(folder, inputs, parent_calibration, output):
    folder,inputs,parent_calibration,output=map(Path,(folder,inputs,parent_calibration,output))
    if output.exists():raise FileExistsError(output)
    calibration=json.loads(parent_calibration.read_text());bank=CalibratedBank(calibration)
    numerator=denominator=0.;records=[]
    for row in calibration['isolated_rows']:
        if not row['used']:continue
        stem=Path(row['root']).stem;wp=inputs/'waveforms'/f'{stem}.npz';tp=inputs/'truth'/f'{stem}.npz'
        wave=dict(np.load(wp));truth=dict(np.load(tp));q=float(truth['charge_pe'][0]);onset=float(truth['first_active_ns'][0])
        kind=int(row['particle']=='gamma');t=wave['time_ns'][::2];y=wave['analog_mV'][::2]-calibration['baseline_mV']
        model=bank.matrix(t,[onset],[kind],[q])[:,0]*q
        intensity=arrival_intensity(t-2.,[onset],[kind],[q],[q],calibration)
        diagonal=covariance_bands(t,intensity,calibration)[0]
        mask=(t>=onset)&(t<onset+2000.)
        numerator+=float(np.sum((y[mask]-model[mask])**2-calibration['noise']['baseline_variance_mV2']))
        denominator+=float(np.sum(diagonal[mask]-calibration['noise']['baseline_variance_mV2']))
        records.append(dict(root=row['root'],waveform_sha256=digest(wp),truth_sha256=digest(tp)))
    calibration['parent_calibration_sha256']=digest(parent_calibration)
    calibration['physical_noise']=dict(excess_variance_scale=max(.01,numerator/denominator),spe_support_ns=240.,quadrature=12,
        calibrated_samples=len(records),data=records,source_sha256=digest(Path(__file__)),
        approximation='Fixed independent charge-arrival covariance through the linear readout; one isolated-calibrated variance multiplier',
        truth_use='Independent isolated calibration only; blind covariance uses voltage-derived pilot parameters')
    save(output,calibration)
    print('SHOT CALIBRATION',calibration['physical_noise']['excess_variance_scale'],len(records),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('folder','inputs','parent_calibration','output'):
        parser.add_argument('--'+name.replace('_','-'),type=Path,required=True)
    args=parser.parse_args();calibrate(args.folder,args.inputs,args.parent_calibration,args.output)
