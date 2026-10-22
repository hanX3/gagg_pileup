"""Grouped native single-pulse pools and explicitly linear pair controls."""
from __future__ import annotations

from dataclasses import replace
import json
import os
from pathlib import Path

import numpy as np
from scipy.signal import lfilter

from .alpha_full_signal_study import make_template
from .alpha_rate_campaign import digest, save
from .alpha_pair_ml import SAMPLES
from .sipm_readout import ReadoutConfig, synthesize_readout
from .sipm_response import AvalancheResult, ELECTRON_CHARGE_C

PROJECT = Path(os.environ.get('GAGG_PAIR_PROJECT',Path(__file__).resolve().parents[2])).resolve()
EXPERIMENT = PROJECT/'pileup/experiments/alpha_pair_curriculum_20261021'
OUTPUT = PROJECT/'jupyter/artifacts/alpha_pair_curriculum_20261021'
RECORD = PROJECT/'jupyter/run_records/alpha_pair_curriculum_20261021'
POOL_SAMPLES = SAMPLES+256
GAPS = (4., 8., 10., 12., 20., 40., 80., 200.)
RATIOS = (.125, .25, .5, 1., 2., 4., 8.)


def load(path):
    return json.loads(Path(path).read_text())


def native_pool(records, destination):
    """Retain original native observed voltage; replay only linear electronics.

    A pool row is one original ROOT, not an independently simulated waveform
    every time it is cropped. Dark chains and electronics noise remain in the
    observed input. The component target excludes dark ancestry.
    """
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    values = {k:[] for k in ('observed', 'component', 'background', 'charge_pe', 'times_ns')}
    provenance = []
    for j,r in enumerate(records):
        meta = load(PROJECT/r['metadata'])
        if meta['primaries'] != 1:
            raise ValueError('The isolated pool requires exactly one primary')
        with np.load(PROJECT/r['waveform']) as f:
            t,y = f['time_ns'].copy(),f['analog_mV'].copy()
        with np.load(PROJECT/r['truth']) as f:
            truth = {key:f[key].copy() for key in f.files}
        if list(truth['particle_name']) != ['alpha']:
            raise ValueError('Pure isolated alpha pool only')
        q = float(truth['charge_pe'][0]); first = float(truth['first_active_ns'][0])
        center = int(np.rint((first-t[0])/2.)) if np.isfinite(first) else 2512
        start = center-256
        if start < 0 or start+POOL_SAMPLES > len(t):
            raise ValueError('Missing isolated pre/post samples')
        q0 = meta['config']['sensor']['gain_electrons']*ELECTRON_CHARGE_C
        owners = truth['avalanche_owner']
        block = np.empty(int((owners==0).sum()), dtype=[('time_ns','f8'),('charge_C','f8')])
        block['time_ns'] = truth['avalanche_time_ns'][owners==0]
        block['charge_C'] = truth['avalanche_charge_C'][owners==0]
        charge_error = abs(float(block['charge_C'].sum()/q0-q))
        assert charge_error < 1e-7
        readout = replace(ReadoutConfig(**meta['config']['readout']),baseline_mV=0.,noise_rms_mV=0.)
        origin = t[start]-readout.aperture_ns/2
        assert not len(block) or block['time_ns'].min() >= origin
        block['time_ns'] -= origin
        # Full-Q above includes every recorded avalanche. A local voltage
        # target only replays avalanches inside its finite readout interval.
        block=block[block['time_ns']<POOL_SAMPLES*2.]
        signal = synthesize_readout(AvalancheResult(block,POOL_SAMPLES*2.,q0,np.array([]),{}),readout,0)
        component = signal.arrays['signal_mV'][:POOL_SAMPLES]
        assert len(component)==POOL_SAMPLES
        assert not np.isfinite(first) or np.all(t[:POOL_SAMPLES] < first)
        values['observed'].append(y[start:start+POOL_SAMPLES])
        values['component'].append(component)
        values['background'].append(y[:POOL_SAMPLES]-meta['config']['readout']['baseline_mV'])
        values['charge_pe'].append(q)
        values['times_ns'].append(first-origin if np.isfinite(first) else 512.)
        provenance.append(dict(**r,input_hashes={k:digest(PROJECT/r[k]) for k in ('waveform','truth','metadata')},
                               first_active_ns=first if np.isfinite(first) else None,charge_pe=q,
                               native_crop_first_index=start,charge_error_pe=charge_error))
        if (j+1)%20==0:
            print('POOL',destination.stem,j+1,len(records),flush=True)
    arrays={k:np.asarray(v,np.float32) for k,v in values.items()}
    np.savez_compressed(destination,**arrays)
    save(destination.with_suffix('.json'),dict(records=provenance,pool_sha256=digest(destination),
        truth_use='Training targets / separate scoring only',
        scope='Native original isolated voltage and actual light-ancestry component; shared-cell state unaltered',
        crop_grid_ns=[.25,2.],parent_split_before_augmentation=True))
    return arrays


def read_pool(split):
    with np.load(OUTPUT/f'{split}-pool.npz') as f:
        return {k:f[k].copy() for k in f.files}


def settings():
    return load(EXPERIMENT/'settings.json')


def single_frames(pool, *, seed, repeats=1, include_weak=False):
    """Original native local voltage, with integer crop jitter, not summed donors."""
    rng=np.random.default_rng(seed)
    eligible=np.arange(len(pool['charge_pe'])) if include_weak else np.flatnonzero(pool['charge_pe']>=500.)
    batch={k:[] for k in ('voltage_mV','active','charge_pe','times_ns','components_mV')}
    parents=[]
    for i in eligible:
        for repetition in range(repeats):
            start=int(rng.integers(96,224)) if repeats>1 else 192
            active=pool['charge_pe'][i]>=500.
            batch['voltage_mV'].append(pool['observed'][i,start:start+SAMPLES])
            batch['active'].append([float(active)])
            batch['charge_pe'].append([pool['charge_pe'][i]])
            batch['times_ns'].append([pool['times_ns'][i]-start*2.])
            batch['components_mV'].append([pool['component'][i,start:start+SAMPLES] if active else np.zeros(SAMPLES)])
            parents.append(int(i))
    return {k:np.asarray(v,np.float32) for k,v in batch.items()},parents


class FrameSampler:
    """Generate single training frames and linear 0/1/2 component controls.

    Pair controls are sums after linear electronics, with one nuisance/noise
    record. They are not native shared-microcell two-event digitizations.
    Any scaling/interpolation is explicitly a software augmentation.
    """
    def __init__(self,pool,config):
        self.pool=pool
        self.qualified=np.flatnonzero(pool['charge_pe']>=500.)
        if not len(self.qualified):
            raise ValueError('No qualified alpha components')
        self.t=np.arange(SAMPLES)*2.+.25
        self.pool_t=np.arange(POOL_SAMPLES)*2.+.25
        template=make_template(config)
        self.kernel_t=np.arange(-24.,6000.01,.25)
        self.kernel=template(self.kernel_t)*config['fit_config']['area_per_pe_mV_ns']*template.unit_area_peak
        self.baseline=config['fit_config']['baseline_mV']

    def mean_curve(self,q,t):
        return q*np.interp(self.t-t,self.kernel_t,self.kernel,left=0.,right=0.)

    def sample(self,rng,size,*,single=False,mean_probability=.25,forced=None):
        out={k:[] for k in ('voltage_mV','active','charge_pe','times_ns','components_mV')}
        info=[]
        for _ in range(size):
            slots=1 if single else 2
            count=1 if single else int(rng.choice((0,1,2),p=(.2,.3,.5)))
            mean=rng.random()<mean_probability
            if forced is not None:
                count=forced['count'];mean=forced.get('mean',True)
            first=float(rng.uniform(64.,320.))
            gap=float(rng.choice(GAPS)) if rng.random()<.6 else float(np.exp(rng.uniform(np.log(4),np.log(200))))
            if forced is not None:
                gap=float(forced.get('gap_ns',gap))
            curves=np.zeros((slots,SAMPLES),np.float32)
            q=np.ones(slots,np.float32)
            times=np.full(slots,first,np.float32)
            active=np.zeros(slots,np.float32)
            parents=[]
            if single and not mean:
                i=int(rng.choice(self.qualified));start=int(rng.integers(96,224))
                y=self.pool['observed'][i,start:start+SAMPLES].copy()
                q[0]=self.pool['charge_pe'][i]
                times[0]=self.pool['times_ns'][i]-2.*start
                curves[0]=self.pool['component'][i,start:start+SAMPLES]
                active[0]=1.;parents=[i]
            else:
                ni=int(rng.integers(len(self.pool['charge_pe'])))
                ns=int(rng.integers(0,257))
                noise=self.pool['background'][ni,ns:ns+SAMPLES].copy()
                if mean:
                    noise=rng.normal(0.,.1,SAMPLES)
                    if rng.random()<.5:
                        noise=lfilter([np.sqrt(1-.98**2)],[1.,-.98],noise)
                desired0=float(np.exp(rng.uniform(np.log(500.),np.log(32000.))))
                if forced is not None:
                    desired0=float(forced.get('q_pe',desired0))
                for j in range(count):
                    times[j]=first+j*gap
                    if mean:
                        q[j]=max(500.,desired0 if j==0 else desired0*float(forced.get('ratio',1.) if forced else np.exp(rng.uniform(np.log(.125),np.log(8.)))))
                        curves[j]=self.mean_curve(q[j],times[j])
                    else:
                        i=int(rng.choice(self.qualified));parents.append(i)
                        scale=float(np.exp(rng.uniform(np.log(.5),np.log(2.))))
                        scale=max(scale,500./self.pool['charge_pe'][i])
                        q[j]=self.pool['charge_pe'][i]*scale
                        curves[j]=scale*np.interp(self.t-times[j]+self.pool['times_ns'][i],
                                                 self.pool_t,self.pool['component'][i],left=0.,right=0.)
                    active[j]=1.
                y=self.baseline+curves.sum(0)+noise
            for key,val in zip(out,(y,active,q,times,curves)):
                out[key].append(val)
            info.append(dict(count=count,mean=bool(mean),parents=parents,gap_ns=gap if count==2 else None))
        return {k:np.asarray(v,np.float32) for k,v in out.items()},info


def save_batch(path,batch,info):
    np.savez_compressed(path,**batch)
    save(Path(path).with_suffix('.json'),dict(cases=info,arrays_sha256=digest(path)))


def read_batch(path):
    with np.load(path) as f:
        return {k:f[k].copy() for k in f.files}
