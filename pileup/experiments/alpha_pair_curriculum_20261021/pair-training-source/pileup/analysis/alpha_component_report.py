"""Saved-data figures for the component-supervision ablation; no model tuning."""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

import numpy as np

from .alpha_component_ml import PROJECT,EXPERIMENT,OUTPUT,VARIANTS,SEEDS
from .alpha_component_study import load
from .alpha_component_evaluate import configuration
from .alpha_full_signal_study import make_template
from .alpha_pileup_fit import score
from .alpha_rate_campaign import save,digest

LABELS=dict(gru_base='GRU / original loss',gru_component='GRU / component loss',
            unetgru_base='U-Net-GRU / original loss',unetgru_component='U-Net-GRU / component loss',
            old_gru='Previous frozen GRU',wls='Current WLS')
COLORS=dict(gru_base='#4477aa',gru_component='#ee7733',unetgru_base='#228833',
            unetgru_component='#cc3377',old_gru='#888888',wls='#111111')
RATES=[100000,1000000,10000000,20000000,40000000,60000000,80000000,100000000]


def plotting():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.18,'savefig.dpi':160,
                         'legend.loc':'upper right'})
    return plt


def finish(fig,name,number,title):
    fig.suptitle(f'Figure 0004-{number}: {title}',fontsize=13)
    folder=OUTPUT/'figures';folder.mkdir(exist_ok=True)
    fig.savefig(folder/(name+'.png'));fig.savefig(folder/(name+'.svg'))
    import matplotlib.pyplot as plt
    plt.close(fig)


def select_gallery_model():
    choices=[]
    for variant in ('gru_component','unetgru_component'):
        for seed in SEEDS:
            name=f'{variant}_{seed}';frozen=load(OUTPUT/'models'/name/'frozen.json')
            choices.append((frozen['validation_objective'],name))
    # The new test scores are deliberately not opened for this choice.
    selected=max(choices)[1]
    save(OUTPUT/'gallery-selection.json',dict(model=selected,criterion='Highest validation macro F1 among four component-loss runs; no test-based selection',
        candidates=[dict(validation_objective=s,name=n) for s,n in choices]))
    return selected


def figures():
    plt=plotting();summary=load(OUTPUT/'test-summary.json');report=load(OUTPUT/'report-tables.json')
    selected=select_gallery_model()
    # Cache mismatch at true parameters is a reference error, not a lower bound.
    split=load(EXPERIMENT/'data-split.json')['windows'];floor=defaultdict(list)
    for record in split:
        with np.load(PROJECT/record['component_cache']) as data:
            floor[record['rate_hz']].extend(data['mean_decoder_floor'].tolist())
    checked=load(OUTPUT/'cache-validation.json')['windows']
    fig,ax=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    ax[0].plot(np.array(RATES)/1e6,[100*np.median(floor[r]) for r in RATES],'-o',label='Median')
    ax[0].plot(np.array(RATES)/1e6,[100*np.quantile(floor[r],.9) for r in RATES],'--s',label='90th percentile')
    ax[0].set(xscale='log',xlabel='Source injection (MHz)',ylabel='Normalized squared waveform mismatch (%)',title='(a) Mean decoder at true Q/time; not an optimized bound')
    ax[0].legend()
    fractions=[]
    for rate in RATES:
        rows=[r for r in checked if r['rate_hz']==rate]
        excluded=sum(r['same_bin_excluded_from_component_loss'] for r in rows);included=sum(r['events'] for r in rows)
        fractions.append(100*excluded/max(excluded+included,1))
    ax[1].plot(np.array(RATES)/1e6,fractions,'-o',color='#aa3377')
    ax[1].set(xscale='log',xlabel='Source injection (MHz)',ylabel='Added-loss exclusions (%)',title='(b) Same-bin labels; no test events excluded')
    finish(fig,'reference_mismatch','23-1','Training/development target checks')
    fig,axes=plt.subplots(1,2,figsize=(13,4.5),layout='constrained')
    for variant in VARIANTS:
        for j,seed in enumerate(SEEDS):
            history=load(OUTPUT/'models'/f'{variant}_{seed}'/'history.json')
            x=[r['step'] for r in history]
            axes[0].plot(x,[r['validation_objective'] for r in history],ls='-' if j==0 else '--',color=COLORS[variant],label=LABELS[variant] if j==0 else None)
            axes[1].plot(x,[r['base_loss'] for r in history],ls='-' if j==0 else '--',color=COLORS[variant])
    axes[0].set(xlabel='Training steps',ylabel='Validation macro charge-qualified F1',title='(a) Same ROOT split; solid/dashed = two training seeds')
    axes[1].set(xlabel='Training steps',ylabel='Original supervised objective',title='(b) Shared loss scale; added terms not mixed into this axis')
    axes[0].legend(fontsize=8)
    finish(fig,'training','24-1','Paired training runs and checkpoint selection')

    def rate_arrays(method,key):
        names=[f'{method}_{seed}' for seed in SEEDS] if method in VARIANTS else [method]
        arrays=[]
        for name in names:
            rows=summary['rates'][name];values=[]
            for r in rows:
                if key=='recovery':v=r['recovered_count']/r['target_count'] if r['target_count'] else np.nan
                elif key=='extras':v=r['unmatched_predictions']/(r['gate_time_ns']/10000.)
                elif key=='charge_failed':v=r['matched_charge_failed']/r['target_count'] if r['target_count'] else np.nan
                else:v=r['missed_count']/r['target_count'] if r['target_count'] else np.nan
                values.append(v)
            arrays.append(values)
        return np.array(arrays)
    methods=[*VARIANTS,'old_gru','wls'];x=np.array(RATES)/1e6
    fig,ax=plt.subplots(figsize=(10,5),layout='constrained')
    for method in methods:
        values=rate_arrays(method,'recovery')*100;mean=np.nanmean(values,axis=0)
        ax.plot(x,mean,'-o',color=COLORS[method],label=LABELS[method])
        if len(values)==2:ax.fill_between(x,np.nanmin(values,axis=0),np.nanmax(values,axis=0),color=COLORS[method],alpha=.12)
    ax.set(xscale='log',xlabel='Source injection (MHz)',ylabel='Complete recovery (%)',ylim=(-2,103))
    ax.legend(fontsize=8,ncol=2)
    finish(fig,'rate_recovery','25-1',f'Fresh {len(summary["windows"])} windows; shaded range is two training seeds, not a CI')
    fig,axes=plt.subplots(1,3,figsize=(17,4.5),layout='constrained')
    for method in methods:
        for ax,key in zip(axes,('extras','charge_failed','missed')):
            values=rate_arrays(method,key)*(1 if key=='extras' else 100)
            ax.plot(x,np.nanmean(values,axis=0),'-o',color=COLORS[method],label=LABELS[method])
    for ax in axes:ax.set_xscale('log');ax.set_xlabel('Source injection (MHz)')
    axes[0].set(ylabel='Extra reports per 10 us',title='(a) Unmatched output')
    axes[1].set(ylabel='Fraction of targets (%)',title='(b) Time matched, Q failed')
    axes[2].set(ylabel='Fraction of targets (%)',title='(c) No time match')
    axes[0].legend(fontsize=7)
    finish(fig,'error_budget','25-2','Recovery must be read together with errors and extra outputs')
    pair_rows=[]
    for seed in SEEDS:
        for architecture in ('gru','unetgru'):
            base=f'{architecture}_base_{seed}';new=f'{architecture}_component_{seed}'
            gain=loss=0
            for row in summary['windows']:
                with np.load(PROJECT/row['truth']) as f:q=f['charge_pe']
                good=[]
                for name in (base,new):
                    m=row['methods'][name];fit=m['fit'];s=m['score']
                    good.append({i for i,j in s['matched_pairs'] if abs(fit['charge_pe'][j]/q[i]-1)<=.2})
                gain+=len(good[1]-good[0]);loss+=len(good[0]-good[1])
            pair_rows.append(dict(seed=seed,architecture=architecture,gained=gain,lost=loss,net=gain-loss))
    fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained');p=np.arange(len(pair_rows))
    ax.bar(p-.18,[r['gained'] for r in pair_rows],.36,label='New complete events',color='#228833')
    ax.bar(p+.18,[-r['lost'] for r in pair_rows],.36,label='Lost complete events',color='#cc6677')
    ax.set_xticks(p,[f'{r["architecture"]}, seed {j//2+1}' for j,r in enumerate(pair_rows)])
    ax.set(ylabel='Paired true-event count',title='Added component loss vs same architecture/seed base loss')
    ax.axhline(0,color='.3',lw=.7);ax.legend()
    finish(fig,'paired_events','25-3','Same-event gains and losses on identical fresh ROOTs')
    save(OUTPUT/'paired-events.json',pair_rows)

    controls=load(OUTPUT/'controls.json')['cases'];gaps=[4,8,10,12,20,40,80,200];ratios=[.125,1.,4.]
    control_methods=[selected,'old_gru','wls']
    fig,axes=plt.subplots(2,3,figsize=(15,7),layout='constrained')
    for i,noise in enumerate(('white','correlated')):
        for j,method in enumerate(control_methods):
            values=np.empty((3,8))
            for a,ratio in enumerate(ratios):
                for b,gap in enumerate(gaps):
                    row=next(r for r in controls if r['noise']==noise and r['label']==f'double_{gap:g}ns_{ratio:g}')
                    values[a,b]=row['methods'][method]['score']['recovered_count']
            im=axes[i,j].imshow(values,vmin=0,vmax=2,aspect='auto',cmap='viridis')
            axes[i,j].set_xticks(range(8),gaps);axes[i,j].set_yticks(range(3),ratios)
            axes[i,j].set(xlabel='True gap (ns)',ylabel='Q2/Q1',title=f'{noise}: {method}')
            for a in range(3):
                for b in range(8):axes[i,j].text(b,a,str(int(values[a,b])),ha='center',va='center',color='white' if values[a,b]<1.5 else 'black')
    fig.colorbar(im,ax=axes,label='Complete components / 2',shrink=.7)
    finish(fig,'software_pairs','26-1','Matched-template software pairs; not a native detector resolution limit')
    fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for ax,method in zip(axes,control_methods):
        rows=[r for r in report['timing_null'] if r['method']==method]
        rates=sorted({r['rate_hz'] for r in rows});actual=[];shifted=[]
        for rate in rates:
            rs=[r for r in rows if r['rate_hz']==rate];den=sum(r['target_count'] for r in rs)
            actual.append(100*sum(r['actual'] for r in rs)/den)
            shifted.append(100*sum(np.mean(r['shifted']) for r in rs)/den)
        ax.plot(np.array(rates)/1e6,actual,'-o',label='Actual times')
        ax.plot(np.array(rates)/1e6,shifted,'--s',label='Five shifted controls')
        ax.set(xlabel='Source injection (MHz)',ylabel='Complete matching (%)',title=method);ax.legend(fontsize=8)
    finish(fig,'timing_null','26-2','Time-shift diagnostic; not an exact FDR or efficiency correction')
    fig,ax=plt.subplots(figsize=(11,4.5),layout='constrained')
    for method in methods:
        names=[f'{method}_{seed}' for seed in SEEDS] if method in VARIANTS else [method]
        values=[]
        for name in names:
            groups=sum((np.asarray(r['groups']) for r in report['neighbor_groups'] if r['method']==name),start=np.zeros((8,3),dtype=int))
            count=groups.sum(1);values.append(np.divide(groups[:,2],count,out=np.full(8,np.nan),where=count>0))
        ax.plot(range(8),100*np.nanmean(values,axis=0),'-o',color=COLORS[method],label=LABELS[method])
    ax.set_xticks(range(8),['<4','4–8','8–12','12–20','20–40','40–80','80–200','>=200'])
    ax.set(xlabel='Nearest other evaluable alpha gap (ns)',ylabel='Complete recovery (%)',ylim=(-2,103),title='Post-inference truth grouping; repeated training seeds are averaged')
    ax.legend(fontsize=8,ncol=2)
    finish(fig,'native_neighbors','26-3','Native close-pulse diagnostic on all fresh windows')
    print('FIGURES_DONE',selected,flush=True)


def gallery():
    from .sipm_readout import ReadoutConfig,synthesize_readout
    from .sipm_response import AvalancheResult,ELECTRON_CHARGE_C
    from matplotlib.lines import Line2D
    plt=plotting();selected=load(OUTPUT/'gallery-selection.json')['model']
    summary=load(OUTPUT/'test-summary.json');settings,_,_=configuration();template=make_template(settings)
    baseline=settings['fit_config']['baseline_mV'];scale=settings['fit_config']['area_per_pe_mV_ns']*template.unit_area_peak
    manifests=[]
    for number,rate in enumerate(RATES,1):
        row=min((r for r in summary['windows'] if r['rate_hz']==rate),key=lambda r:r['seed'])
        stem=Path(row['root']).stem;folder=OUTPUT/'gallery'/f'rate_{rate/1e6:g}MHz'
        folder.mkdir(parents=True,exist_ok=False)
        with np.load(PROJECT/row['waveform']) as f:t,y=f['time_ns'].copy(),f['analog_mV'].copy()
        with np.load(PROJECT/row['truth']) as f:
            truth={k:f[k] for k in ('first_active_ns','emission_ns','charge_pe','avalanche_owner','avalanche_time_ns','avalanche_charge_C','signal_mV')}
        metadata=load(OUTPUT/'inputs/metadata'/(stem+'.json'));duration=metadata['duration_ns']
        q0=metadata['config']['sensor']['gain_electrons']*ELECTRON_CHARGE_C
        readout=replace(ReadoutConfig(**metadata['config']['readout']),baseline_mV=0.,noise_rms_mV=0.)
        owners=truth['avalanche_owner'];order=np.argsort(owners,kind='stable')
        records=np.empty(len(order),dtype=[('time_ns','f8'),('charge_C','f8')])
        records['time_ns']=truth['avalanche_time_ns'][order];records['charge_C']=truth['avalanche_charge_C'][order]
        bounds=np.r_[0,np.cumsum(np.bincount(owners+1,minlength=len(truth['charge_pe'])+1))]
        def voltage(block):return synthesize_readout(AvalancheResult(block,duration,q0,np.array([]),{}),readout,0).arrays['signal_mV']
        ids=np.flatnonzero(truth['charge_pe']>0);actual=np.array([voltage(records[bounds[i+1]:bounds[i+2]]) for i in ids])
        if not len(ids):actual=np.empty((0,len(t)))
        dark=voltage(records[bounds[0]:bounds[1]])
        add_error=float(np.max(abs(actual.sum(0)+dark-truth['signal_mV'])))
        q_error=max([abs(records[bounds[i+1]:bounds[i+2]]['charge_C'].sum()/q0-truth['charge_pe'][i]) for i in ids],default=0.)
        if add_error>1e-8 or q_error>1e-7:raise ValueError('Native ancestry decomposition failed')
        m=row['methods'][selected];fit=m['fit'];rescored=score(fit,truth,row['gate_ns'])
        assert json_compatible(rescored)==m['score']
        times=np.array(fit['times_ns']);charges=np.array(fit['charge_pe'])
        predicted=template(t[None,:]-times[:,None])*(charges*scale)[:,None]
        total=baseline+predicted.sum(0)
        lookup={j:i for i,j in rescored['all_matched_pairs']}
        def color(i):return plt.get_cmap('hsv')((i*.61803398875)%1)
        fig,axes=plt.subplots(3,1,figsize=(17,10),sharex=True,layout='constrained',gridspec_kw={'height_ratios':[1.,1.3,.65]})
        axes[0].plot(t/1000,y,color='.3',lw=.6,label='Observed native voltage')
        axes[0].plot(t/1000,total,color='#cc3377',lw=.8,label='Sum of reported neural components + baseline')
        with np.load(OUTPUT/'test/models'/f'{stem}_wls.npz') as f:wls=f['model_mV'].copy()
        axes[0].plot(t/1000,wls,color='#228833',lw=.7,alpha=.65,label='Current WLS total, including internal terms')
        axes[0].set(ylabel='Total (mV)',title='(a) Complete 0–20 us record; neural total omits unreported weak/dark terms');axes[0].legend(fontsize=8)
        for k,i in enumerate(ids):axes[1].plot(t/1000,actual[k],ls='--',lw=.55,color=color(int(i)) if truth['charge_pe'][i]>=500 else '.7',alpha=.7)
        axes[1].plot(t/1000,dark,color='.55',ls='--',lw=.4)
        for j in range(len(times)):axes[1].plot(t/1000,predicted[j],color=color(lookup[j]) if j in lookup else '#d000d0',lw=.6,ls='-' if j in lookup else ':')
        axes[1].set(ylabel='Component (mV)',title='(b) Every Q>0 native alpha and dark chain; every reported neural component')
        axes[1].legend(handles=[Line2D([],[],color='.3',ls='--',label='Actual alpha'),Line2D([],[],color='.3',label='Time-matched prediction: same color'),Line2D([],[],color='#d000d0',ls=':',label='Extra report'),Line2D([],[],color='.7',ls='--',label='Weak alpha / dark')],fontsize=8,ncol=4)
        axes[2].plot(t/1000,y-total,color='.35',lw=.6)
        axes[2].set(xlabel='Recorded time (us)',ylabel='Residual (mV)',title='(c) Observed minus reported-neural sum; contains missing/nuisance signals')
        for ax in axes:ax.set_xlim(0,20);ax.axvspan(0,5,color='.9',alpha=.45);ax.axvspan(15,20,color='.9',alpha=.45)
        figure=f'27-{number}';fig.suptitle(f'Figure 0004-{figure}: {rate/1e6:g} MHz, first predeclared new ROOT, seed {row["seed"]}\n{selected}; central score {rescored["recovered_count"]}/{rescored["target_count"]}, extras {rescored["unmatched_predictions"]}',fontsize=13)
        fig.savefig(folder/'full_components.png');fig.savefig(folder/'full_components.svg');plt.close(fig)
        np.savez_compressed(folder/'components.npz',time_ns=t,observed_mV=y,actual_truth_ids=ids,actual_components_mV=actual,dark_mV=dark,predicted_times_ns=times,predicted_charge_pe=charges,predicted_components_mV=predicted,neural_total_mV=total,wls_model_mV=wls)
        manifest=dict(figure_id=figure,rate_hz=rate,seed=row['seed'],root=row['root'],root_sha256=row['root_sha256'],model=selected,
            samples=len(t),interval_ns=[0.,20000.],gate_ns=row['gate_ns'],native_components=len(ids),reported_components=len(times),
            native_additivity_max_error_mV=add_error,ancestry_charge_max_error_pe=float(q_error),score=rescored,
            truth_use='Post-inference decomposition, matching colors and unchanged scoring only',
            selection='First predeclared seed at each rate; model selected on validation only',
            neural_total_scope='Reported Q>=500 components only, no hidden dark or weak-truth correction')
        save(folder/'manifest.json',manifest);manifests.append(manifest)
        print('GALLERY',figure,rate,flush=True)
    save(OUTPUT/'gallery/index.json',dict(figures=manifests))


def json_compatible(value):
    import json
    return json.loads(json.dumps(value))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['figures','gallery']);a=p.parse_args()
    figures() if a.stage=='figures' else gallery()
