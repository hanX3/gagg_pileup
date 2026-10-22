"""Historical fixed-range ADC analysis of the pure-alpha count-rate pilot.

Run from project root with BLAS threads limited to one. ROOT acquisition is a
separate explicit command. This analysis never launches G4 or modifies a ROOT.
For the current adequate-dynamic-range comparison, use alpha_full_signal_study;
this historical entry retains the original whole-record ADC rejection policy.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path
import json
import time
import numpy as np

from .alpha_rate_campaign import digest,save,PROJECT
from .alpha_rate_data import export_window
from .alpha_pileup_fit import AlphaTemplate,FitConfig,fit_waveform,score

CALIBRATION=Path('jupyter/artifacts/alpha_gamma_fit_20261013_165844/summary.json')


def setup(output):
    s=json.loads(CALIBRATION.read_text())['results']['alpha']
    p=s['voltage']['adc']['two_fixed']['params']
    h=AlphaTemplate(s['recorded_config']['readout'],p['slow_fraction'],p['tau_fast_ns'],p['tau_slow_ns'])
    cfg=FitConfig()
    metadata=dict(template=h.metadata(),fit_config=asdict(cfg),calibration=str(CALIBRATION),
                  calibration_sha256=digest(CALIBRATION),calibration_root=s['metadata']['root'],
                  threshold_note='500 full-avalanche-equivalent charges; amplitude is realized charge, not deposited energy',
                  true_time='First photon inside the active 6x6 mm plane; source-emission/flight delay is not reconstructed',
                  scoring=dict(time_tolerance_ns=10.,relative_charge_tolerance=.20),
                  fit_truth_inputs=[],rate_input_to_fit=False)
    Path(output).mkdir(parents=True,exist_ok=True)
    save(Path(output)/'analysis_config.json',metadata)
    return h,cfg


def controls(output,h,cfg):
    """Template-model controls; Gaussian electronic noise and actual ADC bins.

    They contain neither new G4 photons nor shared-cell saturation and are not
    estimates of the physical detector's double-pulse resolving time.
    """
    rng=np.random.default_rng(261013900)
    t=np.arange(3000)*2+.25; lsb=1000/4096
    rows=[]
    plans=[('empty',None,None,k) for k in range(24)]
    plans += [('single',None,None,k) for k in range(24)]
    plans += [('double',gap,ratio,k) for gap in (4.,8.,12.,16.,20.,30.,50.,80.,160.)
              for ratio in (.25,.5,1.) for k in range(8)]
    for index,(kind,gap,ratio,rep) in enumerate(plans):
        t0=1000.+rng.uniform(0,2)
        times=np.array([] if kind=='empty' else [t0] if kind=='single' else [t0,t0+gap])
        charges=np.array([] if kind=='empty' else [4000.] if kind=='single' else [4000.,4000.*ratio])
        exact=10.+h(t[:,None]-times)@(charges*cfg.area_per_pe_mV_ns*h.unit_area_peak)
        analog=exact+rng.normal(0,.1,len(t))
        adc=(np.clip(np.floor(analog/lsb),0,4095)+.5)*lsb
        f=fit_waveform(t,adc,h,cfg)
        truth=dict(first_active_ns=times,emission_ns=times,charge_pe=charges)
        sc=score(f,truth,[500.,4000.])
        rows.append(dict(kind=kind,gap_ns=gap,ratio=ratio,repeat=rep,
                         exact_recovery=sc['recovered_count']==len(times) and sc['reconstructed_count']==len(times),
                         **{k:v for k,v in sc.items() if k not in ('target_indices','matched_pairs')}))
        if index%36==0:print('controls',index,'/',len(plans),flush=True)
    save(Path(output)/'controls.json',dict(seed=261013900,model='Frozen template + 0.1 mV Gaussian noise + 12-bit ADC',
        physics_fluctuations=False,rows=rows))
    return rows


def summarize(rows):
    result=[]
    for rate in sorted({r['rate_hz'] for r in rows}):
        subset=[r for r in rows if r['rate_hz']==rate]
        scores=[r['adc_score'] for r in subset]
        totals={k:int(sum(s[k] for s in scores)) for k in (
            'target_count','central_injected','central_nonzero_charge','reconstructed_count','matched_count',
            'recovered_count','missed_count','unmatched_predictions','merge_flags','split_flags')}
        dt=np.array([x for s in scores for x in s['time_errors_ns']])
        dq=np.array([x for s in scores for x in s['relative_charge_errors']])
        def percentiles(a):return [float(x) for x in np.percentile(a,[16,50,84])] if len(a) else None
        result.append(dict(rate_hz=rate,windows=len(subset),gate_time_ns=sum(s['duration_ns'] for s in scores),
            adc_rejected_windows=sum(r['adc_fit']['status']!='fitted' for r in subset),
            clipped_fraction=sum(r['central_clipped_samples'] for r in subset)/sum(r['central_samples'] for r in subset),
            recovery_fraction=totals['recovered_count']/totals['target_count'] if totals['target_count'] else None,
            match_fraction=totals['matched_count']/totals['target_count'] if totals['target_count'] else None,
            charge_relative_p16_p50_p84=percentiles(dq),time_error_ns_p16_p50_p84=percentiles(dt),**totals))
    return result


def analyze(campaign,output,h,cfg):
    output=Path(output);(output/'fits').mkdir(exist_ok=True)
    manifest=json.loads((Path(campaign)/'manifest.json').read_text())
    if any(r['status']!='completed' for r in manifest['runs']):raise ValueError('Acquisition still running or failed')
    rows=[]
    # The pilot was inspected during development and is excluded from the rate table.
    runs=[r for r in manifest['runs'] if r['label']!='pilot_1e6']
    for run in runs:
        stem=Path(run['root']).stem
        print('export/fit',run['label'],flush=True)
        meta=export_window(run['root'],output,expected=run)
        assert meta['root_sha256']==run['root_sha256']
        wave=np.load(output/'waveforms'/f'{stem}.npz')
        begin=time.monotonic()
        fit=fit_waveform(wave['time_ns'],wave['adc_mV'],h,cfg,rail_mask=wave['adc_rail_mask'])
        elapsed=time.monotonic()-begin
        model=fit.pop('model_mV')
        np.savez_compressed(output/'fits'/f'{stem}.npz',model_mV=model)
        # Open labels only after the blind fit has completed.
        truth=np.load(output/'truth'/f'{stem}.npz')
        sc=score(fit,truth,meta['gate_ns'])
        analog=None;asc=None
        if fit['status']!='fitted':
            analog=fit_waveform(wave['time_ns'],wave['analog_mV'],h,cfg)
            amodel=analog.pop('model_mV')
            np.savez_compressed(output/'fits'/f'{stem}_analog_diagnostic.npz',model_mV=amodel)
            asc=score(analog,truth,meta['gate_ns'])
        gate=(wave['time_ns']>=meta['gate_ns'][0])&(wave['time_ns']<meta['gate_ns'][1])
        result=dict(label=run['label'],root=run['root'],root_sha256=meta['root_sha256'],rate_hz=run['rate_hz'],
                    adc_fit=fit,adc_score=sc,analog_fit=analog,analog_score=asc,fit_seconds=elapsed,
                    gate_ns=meta['gate_ns'],central_samples=int(gate.sum()),
                    central_clipped_samples=int(wave['adc_rail_mask'][gate].sum()),
                    analog_peak_mV=float(wave['analog_mV'].max()),
                    tolerance_scan={str(tol):{k:v for k,v in score(fit,truth,meta['gate_ns'],tolerance_ns=tol).items()
                                             if k in ('matched_count','recovered_count','missed_count','merge_flags','split_flags')}
                                    for tol in (5.,10.,20.)})
        save(output/'fits'/f'{stem}.json',result);rows.append(result)
        save(output/'results.partial.json',rows)
        print('result',run['label'],fit['status'],sc['recovered_count'],'/',sc['target_count'],flush=True)
    summary=dict(campaign=str(campaign),campaign_manifest_sha256=digest(Path(campaign)/'manifest.json'),
        pilot_excluded='pilot_1e6',rates=summarize(rows),windows=rows,
        caveats=['Exploratory independent windows, not a final held-out efficiency measurement',
                 'ADC policy rejects the full record on any rail contact; local recovery from clipped windows is not attempted',
                 'Analog diagnostic is pre-ADC simulation information, not available from a clipped ADC record',
                 'Charge is realized shared-cell output including photon-induced crosstalk; no deposited-energy calibration',
                 'Known pure-alpha family and calibration shape; unknown count/time/amplitudes'],
        source_sha256={str(p):digest(p) for p in map(Path,[
            'pileup/analysis/alpha_pileup_fit.py','pileup/analysis/alpha_rate_data.py',
            'pileup/analysis/alpha_rate_study.py','analysis/alpha_gamma_waveform_fit.py'])})
    save(output/'summary.json',summary)
    return summary


def figures(output,h,cfg):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output=Path(output);s=json.loads((output/'summary.json').read_text())
    rows=s['windows']; colors=['#2576b7','#d47b10','#31845c','#b33e52']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':130})
    fig,axes=plt.subplots(4,1,figsize=(12,10),sharex=True,layout='constrained')
    for ax,rate,color in zip(axes,s['rates'],colors):
        row=next(x for x in rows if x['rate_hz']==rate['rate_hz']);stem=Path(row['root']).stem
        a=np.load(output/'waveforms'/f'{stem}.npz');tr=np.load(output/'truth'/f'{stem}.npz');f=np.load(output/'fits'/f'{stem}.npz')
        ax.plot(a['time_ns']/1000,tr['signal_mV']+10,color='#8bc39e',lw=1,label='True pre-noise voltage')
        ax.plot(a['time_ns']/1000,a['adc_mV'],color='#333333',lw=.8,label='ADC bin centers')
        if row['adc_fit']['status']=='fitted':ax.plot(a['time_ns']/1000,f['model_mV'],color=color,lw=1.3,label='Blind fit')
        ax.axhline(1000,color='#b33e52',ls=':',lw=.7)
        ax.text(.015,.95,f"Injection {rate['rate_hz']:.0e} Hz | {row['adc_fit']['status']}\nRecovered {row['adc_score']['recovered_count']}/{row['adc_score']['target_count']} above 500 p.e.",transform=ax.transAxes,va='top',fontsize=9,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
        ax.set_ylabel('Voltage (mV)');ax.set_xlim(5,15);ax.grid(alpha=.15)
    axes[0].legend(loc='upper right',fontsize=8,ncol=3);axes[-1].set_xlabel('Record time (us)')
    fig.savefig(output/'rate_waveforms.png');plt.close(fig)

    control=json.loads((output/'controls.json').read_text())['rows']
    gaps=sorted({x['gap_ns'] for x in control if x['kind']=='double'})
    fig,ax=plt.subplots(figsize=(9,4.5),layout='constrained')
    for ratio,color in zip((.25,.5,1.),colors):
        values=[np.mean([x['exact_recovery'] for x in control if x['gap_ns']==gap and x['ratio']==ratio]) for gap in gaps]
        ax.plot(gaps,np.asarray(values)*100,'o-',color=color,label=f'Q2/Q1={ratio:g}')
    ax.set_xscale('log');ax.set_xticks(gaps,labels=[f'{g:g}' for g in gaps]);ax.set_ylim(-5,105)
    ax.set_xlabel('True pulse separation (ns)');ax.set_ylabel('Both pulses recovered, no extras (%)')
    ax.set_title('Template controls: 8 noise/phase trials per point; no photon/SiPM fluctuations')
    ax.legend();ax.grid(alpha=.2);fig.savefig(output/'double_pulse_controls.png');plt.close(fig)

    fig,ax=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    x=np.arange(4)
    rates=s['rates'];labels=['1e5','1e6','1e7','1e8']
    ax[0,0].bar(x,[100*r['recovery_fraction'] if r['recovery_fraction'] is not None else 0 for r in rates],color=colors)
    for i,r in enumerate(rates):ax[0,0].text(i,5,f"{r['recovered_count']}/{r['target_count']}\n{r['windows']} windows",ha='center',va='bottom',fontsize=8)
    ax[0,0].set_ylabel('ADC recovery including rejections (%)');ax[0,0].set_ylim(0,110)
    ax[0,1].bar(x,[100*r['clipped_fraction'] for r in rates],color=colors);ax[0,1].set_ylabel('Central ADC rail samples (%)')
    for i,r in enumerate(rates):
        for a,key,scale in [(ax[1,0],'charge_relative_p16_p50_p84',100.),(ax[1,1],'time_error_ns_p16_p50_p84',1.)]:
            if r[key] is not None:
                lo,mid,hi=np.array(r[key])*scale;a.errorbar(i,mid,[[mid-lo],[hi-mid]],fmt='o',color=colors[i],capsize=4)
    ax[1,0].set_ylabel('Matched charge error: median, P16-P84 (%)');ax[1,1].set_ylabel('Matched timing error: median, P16-P84 (ns)')
    for a in ax.ravel():a.set_xticks(x,labels);a.set_xlabel('Source alpha injection rate (Hz)');a.grid(axis='y',alpha=.15)
    for a in ax[1]:a.axhline(0,color='gray',ls=':',lw=.8)
    fig.savefig(output/'rate_metrics.png');plt.close(fig)

    # Select the first predeclared 1e6-Hz window; truth only selects a display crop.
    row=next(x for x in rows if x['rate_hz']==1e6)
    stem=Path(row['root']).stem;tr=np.load(output/'truth'/f'{stem}.npz');wave=np.load(output/'waveforms'/f'{stem}.npz')
    fit=row['adc_fit'];tim=tr['first_active_ns'];idx=np.array(row['adc_score']['target_indices'],int)
    if len(idx):
        chosen=idx[np.argmin(np.diff(tim[idx]))] if len(idx)>1 else idx[0]
        lo=max(5000.,tim[chosen]-120.);hi=min(15000.,lo+1200.)
        fig,axes=plt.subplots(3,1,figsize=(11,8),sharex=True,layout='constrained')
        tt=wave['time_ns'];model=np.load(output/'fits'/f'{stem}.npz')['model_mV']
        axes[0].plot(tt,wave['adc_mV'],color='.3',lw=1,label='Observed ADC');axes[0].plot(tt,model,color='#b33e52',lw=1.3,label='Blind total fit')
        axes[0].legend();axes[0].set_ylabel('Total (mV)')
        from .sipm_response import AvalancheResult,AVALANCHE_DTYPE
        from .sipm_readout import ReadoutConfig,synthesize_readout
        rc=ReadoutConfig(**h.readout)
        pairs={i:j for i,j in row['adc_score']['matched_pairs']}
        for k,i in enumerate(np.flatnonzero((tim>lo-1000)&(tim<hi)&(tr['charge_pe']>=500))):
            color=plt.get_cmap('tab10')(k%10);selected=tr['avalanche_owner']==i
            av=np.zeros(int(selected.sum()),dtype=AVALANCHE_DTYPE)
            av['time_ns']=tr['avalanche_time_ns'][selected];av['charge_C']=tr['avalanche_charge_C'][selected]
            component=synthesize_readout(AvalancheResult(av,20000.,1.602176634e-19*1.7e6,np.array([]),{}),rc,261013990).arrays['signal_mV']
            axes[1].plot(tt,component,'--',color=color,lw=1.2,label=f'True pulse {i}')
            if i in pairs:
                j=pairs[i];axes[1].plot(tt,h(tt-fit['times_ns'][j])*fit['peak_mV'][j],color=color,lw=1.)
        axes[1].set_ylabel('Components (mV)');axes[1].set_title('Dashed: actual attributed avalanches; solid: fitted components',fontsize=10)
        axes[1].legend(ncol=4,fontsize=7)
        axes[2].plot(tt,wave['adc_mV']-model,color='.3',lw=.8);axes[2].axhline(0,color='gray',ls=':');axes[2].set_ylabel('Residual (mV)');axes[2].set_xlabel('Record time (ns)')
        for a in axes:a.set_xlim(lo,hi);a.grid(alpha=.15)
        fig.savefig(output/'individual_pulse_recovery.png');plt.close(fig)

    fig,ax=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
    for r,color in zip(rates,colors):
        true=[];reco=[];photon_q=[]
        for row in rows:
            if row['rate_hz']!=r['rate_hz']:continue
            tr=np.load(output/'truth'/f"{Path(row['root']).stem}.npz")
            for i,j in row['adc_score']['matched_pairs']:
                true.append(tr['charge_pe'][i]);reco.append(row['adc_fit']['charge_pe'][j])
            gate=(tr['emission_ns']>=5000)&(tr['emission_ns']<15000)&(tr['active_photons']>=100)
            photon_q.extend((tr['charge_pe'][gate]/tr['active_photons'][gate]).tolist())
        if true:ax[0].scatter(true,reco,s=14,alpha=.65,color=color,label=f"{r['rate_hz']:.0e} Hz")
        if photon_q:ax[1].scatter(np.full(len(photon_q),np.log10(r['rate_hz'])),photon_q,s=9,alpha=.25,color=color)
    ax[0].plot([0,15000],[0,15000],'k:',lw=.8);ax[0].set_xlabel('True realized charge (equivalent p.e.)');ax[0].set_ylabel('Reconstructed charge (equivalent p.e.)');ax[0].legend(fontsize=8)
    ax[1].set_xticks([5,6,7,8],labels);ax[1].set_xlabel('Source alpha injection rate (Hz)');ax[1].set_ylabel('True charge / active-plane photons');ax[1].set_title('Truth diagnostic: includes recovery and crosstalk',fontsize=10)
    fig.savefig(output/'charge_comparison.png');plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--campaign',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--controls-only',action='store_true')
    parser.add_argument('--figures-only',action='store_true');args=parser.parse_args()
    h,cfg=setup(args.output)
    if args.figures_only:figures(args.output,h,cfg)
    elif args.controls_only:controls(args.output,h,cfg)
    else:
        if not (args.output/'controls.json').exists():controls(args.output,h,cfg)
        analyze(args.campaign,args.output,h,cfg);figures(args.output,h,cfg)
