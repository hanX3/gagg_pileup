"""Figures and tables for the alpha/gamma blind-fitting pilot (saved results)."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .alpha_rate_campaign import save

COLORS={'alpha':'#cf5c36','gamma':'#267ba6'}
MODELS={'brems':'Brems gamma','mono662':'661.657 keV gamma'}


def figure_save(fig,output,name):
    fig.savefig(output/(name+'.png'),dpi=160,bbox_inches='tight')
    fig.savefig(output/(name+'.pdf'),bbox_inches='tight');plt.close(fig)


def make_report(output,inputs):
    output,inputs=Path(output),Path(inputs);summary=json.loads((output/'summary.json').read_text())
    rows=json.loads((output/'window-results.json').read_text());groups=summary['groups']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    columns=['gamma_model','alpha_injection_MHz','gamma_injection_MHz','particle','windows','central_injected',
             'central_nonzero_charge','target_count','matched_count','correct_type_count','recovered_count','missed_count',
             'wrong_type_count','signal_rate_MHz','detection_efficiency','conditional_classification','full_recovery_efficiency',
             'mean_relative_charge_error','std_relative_charge_error','time_rmse_ns','unmatched_predictions_class']
    with (output/'metrics.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader()
        writer.writerows({k:g[k] for k in columns} for g in groups)
    fig,axes=plt.subplots(2,2,figsize=(11,8),sharey=True)
    for ax,(model,rate) in zip(axes.flat,[(m,r) for m in MODELS for r in (1.,10.)]):
        subset=[g for g in groups if g['gamma_model']==model and g['alpha_injection_MHz']==rate]
        for j,(metric,ci,label) in enumerate([('detection_efficiency','detection','Detected within 10 ns'),
                ('conditional_classification','classification','Correct type | detected'),
                ('full_recovery_efficiency','recovery','Time + type + charge')]):
            positions=np.arange(2)+(j-1)*.24;values=[100*g[metric] if g[metric] is not None else 0 for g in subset]
            ax.bar(positions,values,width=.23,label=label,alpha=.85)
            for x,y,g in zip(positions,values,subset):
                interval=g['cluster_bootstrap_95'][ci]
                if interval:ax.plot([x,x],100*np.asarray(interval),color='black',lw=1)
                ax.text(x,y+2,f'{y:.0f}',ha='center',fontsize=8)
        ax.set_xticks([0,1],[f"{g['particle']} (N={g['target_count']})" for g in subset])
        ax.set_title(f'{MODELS[model]}, alpha injection {rate:g} MHz');ax.set_ylim(0,113);ax.set_ylabel('Percent')
        ax.grid(axis='y',alpha=.2)
    axes[0,0].legend(fontsize=8,loc='lower left')
    fig.suptitle('Native mixed waveforms: gamma injection 10 MHz; Q truth >= 500 pe\nBars use different stated denominators; lines: whole-window bootstrap 95% intervals')
    fig.tight_layout();figure_save(fig,output,'mixed_efficiencies')

    fig,axes=plt.subplots(2,2,figsize=(10,7))
    for ax,(model,rate) in zip(axes.flat,[(m,r) for m in MODELS for r in (1.,10.)]):
        subset=[r for r in rows if r['row']['gamma_model']==model and r['row']['rate_hz']==rate*1e6]
        matrix=np.array([[sum(r['score']['confusion'][truth][pred] for r in subset) for pred in ('alpha','gamma','missed')]
                          for truth in ('alpha','gamma')])
        denom=matrix.sum(axis=1,keepdims=True);normalized=np.divide(matrix,denom,out=np.zeros_like(matrix,dtype=float),where=denom>0)
        ax.imshow(normalized,vmin=0,vmax=1,cmap='Blues',aspect='auto')
        for i in range(2):
            for j in range(3):ax.text(j,i,f'{matrix[i,j]}\n{100*normalized[i,j]:.1f}%',ha='center',va='center',
                                     color='white' if normalized[i,j]>.6 else 'black')
        ax.set_xticks(range(3),['fit alpha','fit gamma','missed']);ax.set_yticks(range(2),['true alpha','true gamma'])
        ax.set_title(f'{MODELS[model]}, alpha {rate:g} MHz')
    fig.suptitle('Class-blind time matching, then particle classification\nDenominator: all central truths with Q >= 500 pe; misses retained')
    fig.tight_layout();figure_save(fig,output,'mixed_confusion')

    details=[]
    for row in rows:
        stem=Path(row['row']['root']).stem;truth=dict(np.load(inputs/'truth'/f'{stem}.npz'));fit=row['fit']
        match={i:j for i,j in row['score']['matched_pairs']}
        signal=np.flatnonzero((truth['charge_pe']>=500)&np.isfinite(truth['first_active_ns']))
        for i in row['score']['target_indices']:
            others=signal[signal!=i];gap=float(np.min(abs(truth['first_active_ns'][others]-truth['first_active_ns'][i]))) if len(others) else None
            other_type=others[truth['particle_name'][others]!=truth['particle_name'][i]]
            cross_gap=float(np.min(abs(truth['first_active_ns'][other_type]-truth['first_active_ns'][i]))) if len(other_type) else None
            j=match.get(i);error=fit['charge_pe'][j]/truth['charge_pe'][i]-1 if j is not None else None
            correct=j is not None and fit['particle_name'][j]==truth['particle_name'][i]
            details.append(dict(root=row['row']['root'],gamma_model=row['row']['gamma_model'],alpha_injection_MHz=row['row']['rate_hz']/1e6,
                particle=str(truth['particle_name'][i]),truth_index=int(i),charge_pe=float(truth['charge_pe'][i]),
                nearest_signal_ns=gap,nearest_other_type_ns=cross_gap,matched=j is not None,correct_type=bool(correct),
                relative_charge_error=float(error) if error is not None else None,recovered=bool(correct and abs(error)<=.2)))
    save(output/'event-diagnostics.json',details)
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    bins=np.array([0.,20.,50.,100.,300.,1000.,np.inf]);labels=['<20','20–50','50–100','100–300','300–1000','>=1000']
    for ax,(model,rate) in zip(axes.flat,[(m,r) for m in MODELS for r in (1.,10.)]):
        for particle in ('alpha','gamma'):
            subset=[d for d in details if d['gamma_model']==model and d['alpha_injection_MHz']==rate and d['particle']==particle]
            counts=[];good=[]
            for lower,upper in zip(bins[:-1],bins[1:]):
                here=[d for d in subset if lower<=(d['nearest_signal_ns'] if d['nearest_signal_ns'] is not None else 1e100)<upper]
                counts.append(len(here));good.append(sum(d['recovered'] for d in here))
            y=[100*g/n if n else np.nan for g,n in zip(good,counts)]
            ax.plot(range(len(labels)),y,'o-',color=COLORS[particle],label=particle)
            for x,value,n in zip(range(len(labels)),y,counts):
                if n:ax.annotate(f'n={n}',(x,value),xytext=(0,5 if particle=='alpha' else -12),textcoords='offset points',ha='center',fontsize=7)
        ax.set_title(f'{MODELS[model]}, alpha {rate:g} MHz');ax.set_xticks(range(len(labels)),labels,rotation=25)
        ax.set_ylim(-5,110);ax.set_ylabel('Time + type + charge recovery (%)');ax.set_xlabel('Nearest Q>=500 pe signal (ns)');ax.grid(alpha=.2)
    axes[0,0].legend();fig.tight_layout();figure_save(fig,output,'mixed_gap_dependence')

    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for ax,(model,rate) in zip(axes.flat,[(m,r) for m in MODELS for r in (1.,10.)]):
        row=next(r for r in rows if r['row']['gamma_model']==model and r['row']['rate_hz']==rate*1e6)
        stem=Path(row['row']['root']).stem;wave=np.load(inputs/'waveforms'/f'{stem}.npz')
        truth=np.load(inputs/'truth'/f'{stem}.npz');reco=np.load(output/'reconstructions'/f'{stem}.npz')
        mask=(wave['time_ns']>=5500)&(wave['time_ns']<=7500)
        ax.plot(wave['time_ns'][mask]/1000,wave['analog_mV'][mask],color='.6',lw=1,label='observed')
        ax.plot(reco['time_ns'][mask]/1000,reco['model_mV'][mask],color='black',lw=1,label='fit')
        ymin,ymax=ax.get_ylim();height=ymax-ymin
        for particle in ('alpha','gamma'):
            use=(truth['particle_name']==particle)&(truth['charge_pe']>=500)&(truth['first_active_ns']>=5500)&(truth['first_active_ns']<=7500)
            ax.scatter(truth['first_active_ns'][use]/1000,np.full(use.sum(),ymin-.08*height),marker='|',s=70,color=COLORS[particle],label=f'true {particle}')
            use=(reco['particle_name']==particle)&(reco['times_ns']>=5500)&(reco['times_ns']<=7500)
            ax.scatter(reco['times_ns'][use]/1000,np.full(use.sum(),ymin-.16*height),marker='^',s=16,color=COLORS[particle])
        ax.set_title(f'{MODELS[model]}, alpha {rate:g} MHz: first planned window');ax.set_xlabel('Time (us)');ax.set_ylabel('Voltage (mV)')
    axes[0,0].legend(fontsize=7,ncol=2);fig.suptitle('Fixed 5.5–7.5 us excerpt; truth ticks above reconstruction triangles')
    fig.tight_layout();figure_save(fig,output,'mixed_waveform_examples')

    fig,axes=plt.subplots(1,2,figsize=(12,4))
    conditions=[(m,r) for m in MODELS for r in (1.,10.)]
    positions=[];charge_data=[];time_data=[];labels=[];colors=[]
    for i,(model,rate) in enumerate(conditions):
        for j,particle in enumerate(('alpha','gamma')):
            selection=[r['score']['by_particle'][particle] for r in rows if r['row']['gamma_model']==model and r['row']['rate_hz']==rate*1e6]
            charge_data.append([100*v for r in selection for v in r['relative_charge_errors']]);time_data.append([v for r in selection for v in r['time_errors_ns']])
            positions.append(i*3+j);labels.append(f'{model}\n{rate:g}M {particle}');colors.append(COLORS[particle])
    for ax,data,title in zip(axes,(charge_data,time_data),('Charge error (%) | time matched','Time error (ns) | matched')):
        boxes=ax.boxplot(data,positions=positions,tick_labels=labels,patch_artist=True,showfliers=False)
        for patch,color in zip(boxes['boxes'],colors):patch.set_facecolor(color);patch.set_alpha(.65)
        ax.axhline(0,color='.4',lw=.7);ax.set_title(title);ax.tick_params(axis='x',labelsize=7);ax.grid(axis='y',alpha=.2)
    fig.suptitle('Conditional error distributions; boxes: 25–75%, lines: median, whiskers: 1.5 IQR\nIncludes misclassified matches; omitted outliers remain in CSV/JSON statistics')
    fig.tight_layout();figure_save(fig,output,'mixed_conditional_errors')

    control=json.loads((output/'controlled/controlled.json').read_text());cr=control['rows']
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    gaps=[0.,10.,20.,40.,80.,200.,1000.];ratios=[.25,1.,2.]
    for ax,order in zip(axes,[['alpha','gamma'],['gamma','alpha']]):
        counts=np.zeros((3,7),int);totals=np.zeros_like(counts)
        for i,ratio in enumerate(ratios):
            for j,gap in enumerate(gaps):
                subset=[r for r in cr if r['kind']=='double' and r['truth_particle_name']==order and r['gap_ns']==gap and r['gamma_to_alpha_charge_ratio']==ratio]
                totals[i,j]=len(subset);counts[i,j]=sum(r['score']['charge_and_type_recovered_count']==2 and r['score']['unmatched_predictions']==0 for r in subset)
        ax.imshow(counts/np.maximum(totals,1),vmin=0,vmax=1,cmap='Blues',aspect='auto')
        for i in range(3):
            for j in range(7):ax.text(j,i,f'{counts[i,j]}/{totals[i,j]}',ha='center',va='center',color='white' if counts[i,j]/max(totals[i,j],1)>.6 else 'black')
        ax.set_xticks(range(7),[f'{g:g}' for g in gaps]);ax.set_yticks(range(3),[str(r) for r in ratios]);ax.set_xlabel('Onset separation (ns)')
        ax.set_ylabel('Gamma / alpha charge');ax.set_title(' then '.join(order))
    fig.suptitle('Ideal template + 0.1 mV noise controls: both pulses recovered, no extra component\n3 noise realizations per cell; not native photon/SiPM performance')
    fig.tight_layout();figure_save(fig,output,'mixed_controlled_doubles')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True)
    a=p.parse_args();make_report(a.output,a.inputs)
