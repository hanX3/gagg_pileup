"""Class-balanced dense alpha/gamma event reconstruction from measured voltage.

This is a project-specific multitask CNN/U-Net, not a paper reproduction.
Training targets and blind inference are separate interfaces. No event count,
source rate, particle label, charge or true onset is accepted by inference.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import torch
from torch import nn
from torch.nn import functional as F
from .alpha_pileup_fit import FitConfig, candidate_peaks
from .alpha_gamma_regularized_fit import CalibratedBank

NAMES=np.array(['alpha','gamma'])


@dataclass(frozen=True)
class NetworkConfig:
    crop_samples: int = 2048
    guard_samples: int = 600
    channels: int = 24
    label_charge_pe: float = 250.
    report_charge_pe: float = 500.
    charge_scale_pe: float = 1000.
    min_separation_ns: float = 4.

    @property
    def core_samples(self):
        return self.crop_samples-2*self.guard_samples


class DoubleConv(nn.Module):
    def __init__(self,a,b):
        super().__init__()
        self.layers=nn.Sequential(nn.Conv1d(a,b,5,padding=2),nn.SiLU(),
                                  nn.Conv1d(b,b,5,padding=2),nn.SiLU())
    def forward(self,x):return self.layers(x)


class MixedEventNet(nn.Module):
    """Two independent class event maps, log-charge maps and time offsets.

    CNN: 1029-sample receptive field; 600-sample guards exceed its radius.
    U-Net: five downsampling stages and skip connections. Both use eight
    observable channels and emit (batch, class, [event, logQ, offset], sample).
    """
    def __init__(self,kind='cnn',config=NetworkConfig()):
        super().__init__();self.kind=kind;self.config=config
        if kind=='cnn':
            c=config.channels;self.stem=nn.Conv1d(8,c,9,padding=4)
            self.blocks=nn.ModuleList([nn.Sequential(nn.Conv1d(c,c,5,padding=2*d,dilation=d),
                nn.SiLU(),nn.Conv1d(c,c,1)) for d in (1,2,4,8,16,32,64,128)])
            self.head=nn.Conv1d(c,6,1)
        elif kind=='unet':
            sizes=(12,16,24,32,48,64)
            self.encoders=nn.ModuleList([DoubleConv(8 if i==0 else sizes[i-1],s) for i,s in enumerate(sizes)])
            self.decoders=nn.ModuleList([DoubleConv(sizes[i]+sizes[i-1],sizes[i-1]) for i in range(5,0,-1)])
            self.head=nn.Conv1d(sizes[0],6,1)
        else:raise ValueError('Expected cnn or unet')
        with torch.no_grad():
            self.head.bias.zero_();self.head.bias[[0,3]]=-2.2

    def forward(self,x):
        if x.ndim!=3 or x.shape[1]!=8:raise ValueError('Expected eight voltage features')
        if self.kind=='cnn':
            x=F.silu(self.stem(x))
            for block in self.blocks:x=F.silu(x+.3*block(x))
        else:
            skips=[]
            for i,encoder in enumerate(self.encoders):
                if i:x=F.avg_pool1d(x,2)
                x=encoder(x);skips.append(x)
            for decoder,skip in zip(self.decoders,reversed(skips[:-1])):
                x=F.interpolate(x,size=skip.shape[-1],mode='linear',align_corners=False)
                x=decoder(torch.cat([x,skip],dim=1))
        return self.head(x).reshape(x.shape[0],2,3,x.shape[-1])


def voltage_features(time_ns,voltage_mV,calibration):
    """Fixed transforms of the full observed record; no truth normalization."""
    t,y=np.asarray(time_ns,float),np.asarray(voltage_mV,float)
    if t.ndim!=1 or y.shape!=t.shape or len(t)<2 or not np.isfinite(y).all():raise ValueError('Invalid voltage')
    if not np.allclose(np.diff(t),2.,atol=1e-8,rtol=0):raise ValueError('Model requires 2 ns sampling')
    bank=CalibratedBank(calibration);v=y-calibration['baseline_mV']
    inv=[]
    for sigma in (8.,2.):
        cfg=FitConfig(area_per_pe_mV_ns=bank.area,deconvolution_sigma_ns=sigma)
        for h in bank.discovery_bank().templates:inv.append(candidate_peaks(t,v,h,cfg)[1]/200.)
    smooth=gaussian_filter1d(v,50.)
    features=np.array([v/500.,np.gradient(v,2.)/10.,*inv,smooth/500.,
                       v/(50.+np.maximum(smooth,0.))],dtype=np.float32)
    return np.arcsinh(features)


def supervised_targets(t,truth,config=NetworkConfig()):
    """Per-class peaks with nearest-event regression labels, q >= 250 pe.

    Same-class onsets in one sample bin cannot be separately represented;
    their count is retained as a capacity diagnostic, never dropped at scoring.
    """
    t=np.asarray(t,float);n=len(t);labels=np.zeros((2,4,n),np.float32);collisions={}
    names=np.asarray(truth['particle_name']);tt=np.asarray(truth['first_active_ns']);q=np.asarray(truth['charge_pe'])
    for k,name in enumerate(NAMES):
        keep=(names==name)&np.isfinite(tt)&(q>=config.label_charge_pe)&(tt>=t[0])&(tt<=t[-1])
        times,charge=tt[keep],q[keep];order=np.argsort(times);times,charge=times[order],charge[order]
        bins=np.clip(np.rint((times-t[0])/2.).astype(int),0,n-1)
        collisions[name]=len(bins)-len(np.unique(bins))
        if not len(times):continue
        j=np.searchsorted(times,t);a=np.maximum(j-1,0);b=np.minimum(j,len(times)-1)
        ix=np.where(abs(t-times[a])<=abs(t-times[b]),a,b)
        labels[k,0]=np.exp(-.5*(np.arange(n)-bins[ix])**2)
        labels[k,1]=np.log(charge[ix]/config.charge_scale_pe)
        labels[k,2]=(times[ix]-t)/4.
        labels[k,3]=(abs(t-times[ix])<=4.).astype(float)
        # Nearest-event ties must not remove an occupied grid center.
        labels[k,0,bins]=1.
    return labels,collisions


def loss_function(prediction,labels,config=NetworkConfig()):
    a,b=config.guard_samples,config.crop_samples-config.guard_samples
    z,y=prediction[:,:,:,a:b],labels[:,:,:,a:b]
    heat,logq,offset,mask=y.unbind(2)
    p=torch.sigmoid(z[:,:,0]).clamp(1e-6,1-1e-6)
    positive=(heat>=1.-1e-6).to(p.dtype)
    negative=(1.-positive)*(1.-heat)**4
    norm=positive.sum((0,2)).clamp(min=1.)
    event=(-torch.log(p)*(1-p)**2*positive-torch.log1p(-p)*p**2*negative).sum((0,2))/norm
    den=mask.sum((0,2)).clamp(min=1.)
    charge=(F.smooth_l1_loss(z[:,:,1],logq,reduction='none')*mask).sum((0,2))/den
    timing=(F.smooth_l1_loss(torch.tanh(z[:,:,2]),offset,reduction='none')*mask).sum((0,2))/den
    terms=torch.stack([event.mean(),charge.mean(),timing.mean()])
    return terms[0]+2.*terms[1]+.5*terms[2],terms


def predict_features(model,features,config=NetworkConfig(),batch_size=4):
    if features.ndim!=2 or features.shape[0]!=8 or not np.isfinite(features).all():raise ValueError('Invalid features')
    n=features.shape[-1];core=config.core_samples;g=config.guard_samples
    starts=list(range(0,n,core));padded=np.pad(features,((0,0),(g,config.crop_samples)),mode='constant')
    output=np.empty((2,3,n),np.float32)
    model.eval()
    with torch.inference_mode():
        for i in range(0,len(starts),batch_size):
            ss=starts[i:i+batch_size]
            x=np.stack([padded[:,s:s+config.crop_samples] for s in ss])
            z=model(torch.from_numpy(x)).cpu().numpy()
            for j,s in enumerate(ss):
                length=min(core,n-s);output[:,:,s:s+length]=z[j,:,:,g:g+length]
    return dict(event_score=1/(1+np.exp(-np.clip(output[:,0],-30,30))),
                charge_pe=config.charge_scale_pe*np.exp(np.clip(output[:,1],-6,6)),
                time_offset_ns=4*np.tanh(output[:,2]))


def decode(t,dense,thresholds=(.5,.5),config=NetworkConfig()):
    t=np.asarray(t);candidates=[]
    for k in range(2):
        peaks=find_peaks(dense['event_score'][k],height=thresholds[k],distance=2)[0]
        for i in peaks:
            q=float(dense['charge_pe'][k,i])
            if q<config.label_charge_pe:continue
            candidates.append((float(dense['event_score'][k,i]),float(t[i]+dense['time_offset_ns'][k,i]),k,q))
    kept=[]
    for entry in sorted(candidates,reverse=True):
        if all(abs(entry[1]-old[1])>=config.min_separation_ns for old in kept):kept.append(entry)
    kept.sort(key=lambda r:r[1])
    def packed(rows):
        return dict(times_ns=[r[1] for r in rows],particle_name=[str(NAMES[r[2]]) for r in rows],
                    charge_pe=[r[3] for r in rows],event_score=[r[0] for r in rows])
    return dict(status='predicted',**packed([r for r in kept if r[3]>=config.report_charge_pe]),
                all_components=packed(kept),thresholds=list(map(float,thresholds)),
                input_policy='complete_observed_voltage',fit_truth_inputs=[],
                minimum_output_spacing_ns=config.min_separation_ns,
                scores_are_calibrated_probabilities=False,config=asdict(config))


def predict_waveform(model,t,voltage,calibration,thresholds,config=NetworkConfig()):
    return decode(t,predict_features(model,voltage_features(t,voltage,calibration),config),thresholds,config)
