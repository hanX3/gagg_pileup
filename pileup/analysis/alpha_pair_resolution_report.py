"""Evaluation of saved blind predictions, with distinct oracle diagnostics."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
import json
from pathlib import Path

import numpy as np

from .alpha_pair_fitting import PairFitter
from .alpha_pileup_fit import score
from .alpha_pair_resolution_data import (
    PROJECT, EXPERIMENT, OUTPUT, MODEL_NAMES, GAPS, GRID,
    protocol, settings, load, save, digest)

DOMAINS=('mean','shared','linear','shared-noise')
METHODS=('blind_fit','ml_1','hybrid_1','ml_2','hybrid_2','known_K_fit','known_K_t_GLS')
BLIND=METHODS[:5]
COLORS=dict(blind_fit='#1f77b4',ml_1='#168c82',hybrid_1='#d62728',
            ml_2='#9757bb',hybrid_2='#e58b25',known_K_fit='#8855aa',known_K_t_GLS='#222222')
LABELS=dict(blind_fit='Blind fit',ml_1='ML seed 1',hybrid_1='ML + fit seed 1',
            ml_2='ML seed 2',hybrid_2='ML + fit seed 2',
            known_K_fit='Known K (diagnostic)',known_K_t_GLS='Known K,t GLS (diagnostic)')


def wilson(success,n,z=1.959963984540054):
    if not n:return [None,None]
    p=success/n;d=1.+z*z/n
    center=(p+z*z/(2*n))/d
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return [float(max(0.,center-half)),float(min(1.,center+half))]


def quantiles(values):
    return np.percentile(values,[16,50,84]).astype(float).tolist() if len(values) else None


def predictions(domain):
    folder=OUTPUT/'predictions'/domain
    complete=load(folder/'complete.json')
    assert complete['waveform_sha256']==digest(OUTPUT/'assets'/(domain+'.waveforms.npz'))
    candidates=[None]*complete['events']
    for line in (folder/'candidates.jsonl').read_text().splitlines():
        row=json.loads(line);assert candidates[row['index']] is None
        candidates[row['index']]=row['candidates']
    assert all(c is not None for c in candidates)
    penalty=protocol()['fitting_penalty']
    fits={'blind_fit':[PairFitter.choose(c,penalty) for c in candidates]}
    for j,name in enumerate(MODEL_NAMES,1):
        neural=load(folder/(name+'.json'))
        fits[f'ml_{j}']=neural['fits']
        fits[f'hybrid_{j}']=[r['fit'] for r in neural['refined']]
    return fits,candidates


def native_false_conditions(domain,condition):
    if domain=='mean':return [condition],[condition]
    if domain in ('linear','shared-noise'):return ['all'],['mid']
    return ['all'],list(dict.fromkeys(condition.split('_')))


def make_report():
    report=dict(curves=[],false_controls=[],delta90=[],coincident=[],inputs={})
    event_rows=[]
    for domain in DOMAINS:
        info=load(OUTPUT/'assets'/(domain+'.json'))['cases']
        with np.load(OUTPUT/'assets'/(domain+'.truth.npz')) as f:
            q=f['charge_pe'].copy();t=f['times_ns'].copy();active=f['active'].copy()
        fits,candidates=predictions(domain)
        diagnostic=load(OUTPUT/'diagnostics'/(domain+'.json'))['rows']
        fits['known_K_fit']=[c[min(int(a.sum()),len(c)-1)] for c,a in zip(candidates,active)]
        fits['known_K_t_GLS']=[r.get('known_times_fit') for r in diagnostic]
        buckets=defaultdict(list);false_buckets=defaultdict(list);zero_buckets=defaultdict(list)
        for i,c in enumerate(info):
            ids=np.flatnonzero(active[i]);true_times=t[i,ids];true_q=q[i,ids]
            truth=dict(first_active_ns=true_times,emission_ns=true_times,charge_pe=true_q)
            true_gap=float(np.diff(np.sort(true_times))[0]) if len(ids)==2 else None
            for method in METHODS:
                fit=fits[method][i]
                if fit is None:
                    if c['kind']=='pair' and c['gap_ns']==0.:continue
                    raise ValueError((domain,i,method,'unexpected unavailable diagnostic'))
                normal=score(fit,truth,[0.,2048.])
                strict=score(fit,truth,[0.,2048.],tolerance_ns=min(10.,true_gap/4.)) if true_gap and true_gap>0 else None
                row=dict(domain=domain,condition=c['condition'],kind=c['kind'],index=i,
                         gap_ns=c['gap_ns'],method=method,oracle=method not in BLIND,
                         target_count=len(ids),predicted_count=len(fit['times_ns']),
                         returned_slots_below_500_pe=sum(float(v)<500. for v in fit['charge_pe']),
                         ordinary_recovered=normal['recovered_count'],
                         ordinary_matched=normal['matched_count'],
                         extra=normal['unmatched_predictions'],missed=normal['missed_count'],
                         split_flags=normal['split_flags'],merge_flags=normal['merge_flags'],
                         q_errors=normal['relative_charge_errors'],time_errors=normal['time_errors_ns'],
                         strict_recovered=strict['recovered_count'] if strict else None,
                         strict_matched=strict['matched_count'] if strict else None,
                         cluster=c['cluster'],true_q=true_q.tolist(),
                         sum_q_relative_error=float((sum(fit['charge_pe'])/true_q.sum())-1.) if true_q.sum()>0 else None)
                if len(ids)==2 and strict:
                    row['pair_success']=strict['recovered_count']==2 and normal['unmatched_predictions']==0
                    row['ordinary_pair_success']=normal['recovered_count']==2 and normal['unmatched_predictions']==0
                    if row['predicted_count']<2:row['failure_category']='count_low'
                    elif row['predicted_count']>2:row['failure_category']='count_high'
                    elif strict['matched_count']<2:row['failure_category']='time'
                    elif strict['recovered_count']<2:row['failure_category']='charge'
                    else:row['failure_category']='success'
                else:row['pair_success']=None
                event_rows.append(row)
                if c['kind']=='pair' and c['gap_ns']>0:
                    buckets[(c['condition'],c['gap_ns'],method)].append(row)
                elif c['kind']=='pair':
                    zero_buckets[(c['condition'],method)].append(row)
                else:
                    false_buckets[(c['condition'],c['kind'],method)].append(row)
        for (condition,gap,method),rows in buckets.items():
            eligible=[r for r in rows if r['target_count']==2]
            assert len({r['cluster'] for r in rows})==len(rows), 'Repeated parents within one condition/gap.'
            n=len(eligible);s=sum(r['pair_success'] for r in eligible)
            qerrors=[x for r in rows for x in r['q_errors']]
            terrors=[x for r in rows for x in r['time_errors']]
            fq=defaultdict(int)
            for r in eligible:fq[r['failure_category']]+=1
            qs=np.array([r['true_q'] for r in eligible])
            report['curves'].append(dict(domain=domain,condition=condition,gap_ns=gap,method=method,
                generated_cases=len(rows),eligible_pairs=n,success=s,efficiency=s/n if n else None,
                wilson95=wilson(s,n),ordinary_success=sum(r['ordinary_pair_success'] for r in eligible),
                target_events=sum(r['target_count'] for r in rows),
                recovered_events=sum(r['ordinary_recovered'] for r in rows),
                strict_recovered_events=sum(r['strict_recovered'] or 0 for r in rows),
                extra_outputs=sum(r['extra'] for r in rows),missed_events=sum(r['missed'] for r in rows),
                returned_slots_below_500_pe=sum(r['returned_slots_below_500_pe'] for r in rows),
                split_flags=sum(r['split_flags'] for r in rows),merge_flags=sum(r['merge_flags'] for r in rows),
                q_error_p16_p50_p84=quantiles(qerrors),time_error_p16_p50_p84_ns=quantiles(terrors),
                q_error_std=float(np.std(qerrors,ddof=1)) if len(qerrors)>1 else None,
                q1_p16_p50_p84_pe=quantiles(qs[:,0]) if n else None,
                q2_p16_p50_p84_pe=quantiles(qs[:,1]) if n else None,
                ratio_q2_q1_p16_p50_p84=quantiles(qs[:,1]/qs[:,0]) if n else None,
                failure_counts=dict(fq),oracle=method not in BLIND))
        for (condition,kind,method),rows in false_buckets.items():
            n=len(rows);f=sum(r['extra']>0 for r in rows)
            report['false_controls'].append(dict(domain=domain,condition=condition,kind=kind,method=method,
                frames=n,false_output_frames=f,false_output_fraction=f/n,wilson95=wilson(f,n),
                target_events=sum(r['target_count'] for r in rows),
                recovered_events=sum(r['ordinary_recovered'] for r in rows),
                missed_events=sum(r['missed'] for r in rows),extra_outputs=sum(r['extra'] for r in rows),
                returned_slots_below_500_pe=sum(r['returned_slots_below_500_pe'] for r in rows),
                split_flags=sum(r['split_flags'] for r in rows),merge_flags=sum(r['merge_flags'] for r in rows),
                q_error_p16_p50_p84=quantiles([x for r in rows for x in r['q_errors']]),
                time_error_p16_p50_p84_ns=quantiles([x for r in rows for x in r['time_errors']]),
                oracle=method not in BLIND))
        for (condition,method),rows in zero_buckets.items():
            report['coincident'].append(dict(domain=domain,condition=condition,method=method,
                cases=len(rows),pair_success=None,
                predicted_two=sum(r['predicted_count']==2 for r in rows),
                sum_q_error_p16_p50_p84=quantiles([r['sum_q_relative_error'] for r in rows]),
                explanation='Same-class exact coincidence: individual correctness not scored.'))
        report['inputs'][domain]=dict(waveform_sha256=digest(OUTPUT/'assets'/(domain+'.waveforms.npz')),
                                     truth_sha256=digest(OUTPUT/'assets'/(domain+'.truth.npz')),
                                     predicted_before_oracle=True,events=len(info))
    for domain,condition,method in sorted({(r['domain'],r['condition'],r['method']) for r in report['curves']}):
        curve=sorted([r for r in report['curves'] if (r['domain'],r['condition'],r['method'])==(domain,condition,method)],
                     key=lambda r:r['gap_ns'])
        possible=[r['gap_ns'] for i,r in enumerate(curve)
                  if all(x['efficiency'] is not None and x['efficiency']>=.9 for x in curve[i:])]
        supported=[r['gap_ns'] for i,r in enumerate(curve)
                   if all(x['wilson95'][0] is not None and x['wilson95'][0]>=.9 for x in curve[i:])]
        half=[r['gap_ns'] for i,r in enumerate(curve)
              if all(x['efficiency'] is not None and x['efficiency']>=.5 for x in curve[i:])]
        empty,single=native_false_conditions(domain,condition)
        controls=[r for r in report['false_controls'] if r['domain']==domain and r['method']==method
                  and ((r['kind']=='empty' and r['condition'] in empty) or
                       (r['kind']=='single' and r['condition'] in single))]
        assert len(controls)==len(empty)+len(single)
        budget=all(r['false_output_fraction']<=.05 for r in controls)
        confident=all(r['wilson95'][1]<=.05 for r in controls)
        report['delta90'].append(dict(domain=domain,condition=condition,method=method,
            delta90_point_ns=min(possible) if possible else None,
            delta50_point_ns=min(half) if half else None,
            delta90_pointwise95_ns=min(supported) if supported else None,
            false_budget_observed_pass=budget,false_budget_pointwise95_pass=confident,
            qualified_delta90_point_ns=min(possible) if possible and budget and method in BLIND else None,
            oracle=method not in BLIND,
            caution='Pointwise intervals, not simultaneous confidence; inherited networks trained on4–200ns,400/800ns extrapolate.'))
    report.update(protocol_sha256=digest(EXPERIMENT/'protocol.json'),
        freeze_sha256=digest(EXPERIMENT/'freeze.json'),
        confidence='Wilson95 per condition/gap with independent pairs; the same parents recur across gaps/conditions. No pooled independent-exposure claim.',
        strict_pair_definition='Both Q errors<=20%; one-to-one time errors<=min(10ns,true gap/4); no unmatched output. Delta0 excluded.',
        blind_methods=list(BLIND),oracle_methods=list(METHODS[5:]),
        interpretation='KnownK is one specific solver, not an upper bound; electronic-noise Fisher estimates omit native shape/shot fluctuations.',
        no_retraining=True,no_test_threshold_tuning=True)
    with gzip.open(OUTPUT/'score-rows.jsonl.gz','wt',encoding='utf-8') as stream:
        for row in event_rows:stream.write(json.dumps(row,ensure_ascii=False,allow_nan=False)+'\n')
    save(OUTPUT/'resolution-report.json',report)
    for key,filename in (('curves','curve-summary.csv'),('false_controls','false-controls.csv'),
                         ('delta90','resolution-thresholds.csv')):
        with (OUTPUT/filename).open('w',newline='',encoding='utf-8-sig') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(report[key][0]))
            writer.writeheader()
            writer.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v
                              for k,v in row.items()} for row in report[key])
    print('REPORT',len(event_rows),'scores;',len(report['curves']),'curve points',flush=True)
    return report


def _curve(ax,report,domain,condition,methods,uncertainty=False):
    for method in methods:
        rows=sorted([r for r in report['curves'] if r['domain']==domain and r['condition']==condition and r['method']==method],
                    key=lambda r:r['gap_ns'])
        x=np.array([r['gap_ns'] for r in rows]);y=np.array([r['efficiency'] for r in rows],float)
        style='--' if method not in BLIND else (':' if method.startswith('ml_') else '-')
        ax.plot(x,100*y,style,marker='.',ms=5,lw=1.4,color=COLORS[method],label=LABELS[method])
        if uncertainty:
            lo=np.array([r['wilson95'][0] for r in rows]);hi=np.array([r['wilson95'][1] for r in rows])
            ax.fill_between(x,100*lo,100*hi,color=COLORS[method],alpha=.09)
    ax.axhline(90,color='black',ls=':',lw=.8)
    ax.axvspan(200,800,color='gray',alpha=.08)
    ax.set_xscale('symlog',linthresh=20)
    ax.set_xticks([2,8,16,40,80,200,800],labels=['2','8','16','40','80','200','800'])
    ax.set(xlim=(0,830),ylim=(0,103),xlabel='Separation (ns)',ylabel='Strict pair recovery (%)')


def figures(report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9,'axes.grid':True,'grid.alpha':.2,'savefig.dpi':160})
    out=OUTPUT/'figures';out.mkdir(exist_ok=True)
    main=('blind_fit','hybrid_1','hybrid_2','known_K_fit','known_K_t_GLS')
    fig,axes=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
    _curve(axes[0],report,'mean','reference',main,uncertainty=True)
    axes[0].set_title('(a) Fixed mean: Q1=Q2=4000 pe, noise0.1 mV')
    info=load(OUTPUT/'assets/mean.json')['cases'];diag=load(OUTPUT/'diagnostics/mean.json')['rows']
    selected=[i for i,c in enumerate(info) if c['kind']=='pair' and c['condition']=='reference' and c['replicate']==0 and c['gap_ns']>0]
    rows=[(info[i],diag[i]) for i in selected]
    # The original frozen diagnostic uses a smoothed grid gradient. At short
    # gaps its small derivative error can be amplified by the inverse Fisher
    # matrix. Evaluate the actual interpolated mean model by finite differences
    # here; keep and plot the original approximation for transparency.
    with np.load(OUTPUT/'assets/mean.truth.npz') as f:
        times=f['times_ns'];charges=f['charge_pe']
    fitter=PairFitter(settings());finite_sd=[];sensitivity=[]
    for i in selected:
        tt=times[i];qq=charges[i];eps=1e-5
        H=fitter.columns(GRID,tt)
        slope=(fitter.columns(GRID,tt-eps)-fitter.columns(GRID,tt+eps))/(2*eps)
        J=np.column_stack((H*qq[None],-slope*qq[None]))/.1
        scales=np.linalg.norm(J,axis=0);normalized=J/scales[None]
        gram=normalized.T@normalized
        covariance=np.linalg.inv(gram)/scales[:,None]/scales[None]
        sd=np.sqrt(np.diag(covariance));finite_sd.append(float(max(sd[:2])))
        original=np.array(diag[i]['sigma_relative_q_unknown_times']+diag[i]['sigma_time_unknown_times_ns'])
        sensitivity.append(dict(index=i,gap_ns=info[i]['gap_ns'],
            finite_difference_step_ns=eps,finite_difference_standard_deviations=sd.tolist(),
            stored_grid_gradient_standard_deviations=original.tolist(),
            maximum_relative_difference=float(np.max(np.abs(original/sd-1.)))))
    save(OUTPUT/'information-sensitivity.json',dict(rows=sensitivity,
        affects='Conditional mean-model information plot only; no change to frozen blind fits or scores.',
        scope='Reference Q4000/4000, 0.1mV white noise, replicate0; local known-K information, not a global limit.'))
    x=[c['gap_ns'] for c,d in rows]
    axes[1].plot(x,[d['template_correlation'] for c,d in rows],'.-',color='#228866')
    axes[1].set(xscale='log',xlabel='Separation (ns)',ylabel='Whitened template correlation',
                title='(b) Known times: template similarity')
    a=[100*max(d.get('sigma_q_known_time_pe',[np.nan]))/4000 for c,d in rows]
    b=[100*max(d.get('sigma_relative_q_unknown_times',[np.nan])) for c,d in rows]
    axes[2].loglog(x,a,'.-',label='Known times GLS')
    axes[2].loglog(x,100*np.asarray(finite_sd),'.-',label='Free times: numerical-derivative Fisher')
    axes[2].loglog(x,b,':',color='gray',label='Stored grid-gradient approximation')
    axes[2].axhline(20,color='gray',ls=':',label='20% charge tolerance')
    axes[2].set(xlabel='Separation (ns)',ylabel='Relative Q standard deviation (%)',
                title='(c) Exact mean / electronic noise only')
    axes[2].legend(fontsize=8);axes[0].legend(fontsize=7,loc='lower right')
    fig.suptitle('Figure 0004-37-1: Algorithm recovery versus conditional information diagnostics')
    fig.savefig(out/'mean_information.png');plt.close(fig)

    fig,axes=plt.subplots(2,4,figsize=(16,8),layout='constrained')
    for i,c in enumerate(protocol()['mean_conditions']):
        ax=axes.flat[i];_curve(ax,report,'mean',c['name'],main)
        ax.set_title(f"({chr(97+i)}) {c['name']}\nQ={c['q1']:g}/{c['q2']:g} pe, noise={c['noise_rms']:g} mV")
    axes.flat[7].axis('off')
    handles,labels=axes.flat[0].get_legend_handles_labels()
    axes.flat[7].legend(handles,labels,loc='center',fontsize=9)
    axes.flat[7].text(.02,.15,'40 paired noise realizations / point\nGray region: ML gap extrapolation\nDashed curves use true K or true K,t',transform=axes.flat[7].transAxes)
    fig.suptitle('Figure 0004-37-2: Signal strength, pulse order and noise in fixed-mean controls')
    fig.savefig(out/'mean_conditions.png');plt.close(fig)

    fig,axes=plt.subplots(2,3,figsize=(15,8),layout='constrained')
    for i,c in enumerate(protocol()['shared_conditions']):
        ax=axes.flat[i];_curve(ax,report,'shared',c,main,uncertainty=True)
        ax.set_title(f"({chr(97+i)}) Shared cells: {c}")
    ax=axes.flat[5]
    for domain,color in (('linear','#9467bd'),('shared','#1f77b4'),('shared-noise','#d62728')):
        condition='mid_mid'
        rows=sorted([r for r in report['curves'] if r['domain']==domain and r['condition']==condition and r['method']=='blind_fit'],
                    key=lambda r:r['gap_ns'])
        ax.plot([r['gap_ns'] for r in rows],[100*r['efficiency'] for r in rows],'.-',label=domain,color=color)
    ax.set(xscale='symlog',xlim=(0,830),ylim=(0,103),xlabel='Separation (ns)',ylabel='Strict pair recovery (%)',
           title='(f) Blind fit: same optical parents\nDifferent digitization / noise')
    ax.set_xscale('symlog',linthresh=20);ax.axhline(90,color='black',ls=':',lw=.8);ax.legend(fontsize=8)
    handles,labels=axes.flat[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=3,fontsize=8)
    fig.suptitle('Figure 0004-38-1: Fresh optical parents and shared microcells; 40 independent pairs / condition')
    fig.savefig(out/'shared_conditions.png');plt.close(fig)

    fig,axes=plt.subplots(2,2,figsize=(14,8),layout='constrained')
    categories=('success','count_low','time','charge','count_high')
    colors=('#29975d','#d9534f','#eead42','#8b68b6','#555555')
    for i,(condition,method) in enumerate((('mid_mid','blind_fit'),('mid_mid','hybrid_1'),
                                         ('low_high','blind_fit'),('low_high','hybrid_1'))):
        ax=axes.flat[i]
        rows=sorted([r for r in report['curves'] if r['domain']=='shared' and r['condition']==condition and
                     r['method']==method and r['gap_ns']<=200],key=lambda r:r['gap_ns'])
        x=np.arange(len(rows));bottom=np.zeros(len(rows))
        for cat,color in zip(categories,colors):
            height=np.array([100*r['failure_counts'].get(cat,0)/r['eligible_pairs'] for r in rows])
            ax.bar(x,height,bottom=bottom,label=cat,color=color,width=.85);bottom+=height
        ax.set_xticks(x,[f"{r['gap_ns']:g}" for r in rows],rotation=45)
        ax.set(ylim=(0,100),xlabel='Separation (ns)',ylabel='Pair fraction (%)',
               title=f"({chr(97+i)}) {condition} / {LABELS[method]}")
    handles,labels=axes.flat[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=5)
    fig.suptitle('Figure 0004-38-2: Exclusive failure categories; count, then timing, then charge')
    fig.savefig(out/'failure_categories.png');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(14,5),layout='constrained')
    controls=report['false_controls'];kinds=[('all','empty'),('low','single'),('mid','single'),('high','single')]
    for j,m in enumerate(BLIND):
        rows=[next(r for r in controls if r['domain']=='shared' and r['method']==m and r['kind']==kind and r['condition']==cond)
              for cond,kind in kinds]
        y=np.array([r['false_output_fraction'] for r in rows])*100
        lo=np.array([r['wilson95'][0] for r in rows])*100;hi=np.array([r['wilson95'][1] for r in rows])*100
        x=np.arange(len(rows))+(j-2)*.13
        axes[0].errorbar(x,y,yerr=np.maximum(0.,np.vstack((y-lo,hi-y))),fmt='.',capsize=2,color=COLORS[m],label=LABELS[m])
        singles=rows[1:];eff=np.array([r['recovered_events']/r['target_events'] for r in singles])*100
        axes[1].plot(np.arange(3)+(j-2)*.06,eff,'.-',color=COLORS[m],label=LABELS[m])
    axes[0].axhline(5,color='black',ls=':',label='5% budget')
    axes[0].set(xticks=range(4),xticklabels=['Empty','Low single','Mid single','High single'],
                ylabel='Frames with extra output (%)',title='(a) Point estimates and pointwise Wilson95')
    axes[1].set(xticks=range(3),xticklabels=['Low single','Mid single','High single'],
                ylabel='Q20% / time10ns recovery (%)',ylim=(0,103),title='(b) Preserve isolated-event recovery')
    axes[1].axhline(90,color='black',ls=':',lw=.8)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=3,fontsize=8)
    fig.suptitle('Figure 0004-38-3: Fresh empty/single controls; frozen thresholds')
    fig.savefig(out/'false_outputs.png');plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(15,5),layout='constrained')
    for i,condition in enumerate(('mid_mid','high_low','low_low')):
        _curve(axes[i],report,'shared',condition,BLIND)
        axes[i].set_title(f'({chr(97+i)}) {condition}: raw ML versus refit')
    fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=3,fontsize=8)
    fig.suptitle('Figure 0004-38-4: A physical-template refit can lose correct neural reconstructions')
    fig.savefig(out/'ml_refit_comparison.png');plt.close(fig)
    examples(report)
    print('FIGURES',7,flush=True)


def examples(report):
    import matplotlib.pyplot as plt
    info=load(OUTPUT/'assets/shared.json')['cases'];fits,_=predictions('shared')
    with np.load(OUTPUT/'assets/shared.waveforms.npz') as f:y=f['voltage_mV']
    with np.load(OUTPUT/'assets/shared.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    with np.load(OUTPUT/'assets/shared.pair-components.npz') as f:components=f['components_mV']
    fitter=PairFitter(settings())
    fig,axes=plt.subplots(3,2,figsize=(15,10),layout='constrained')
    selections=[]
    for col,gap in enumerate((4.,80.)):
        index=next(i for i,c in enumerate(info) if c['kind']=='pair' and c['condition']=='mid_mid' and c['gap_ns']==gap and c['replicate']==0)
        fit=fits['hybrid_1'][index]
        H=fitter.columns(GRID,fit['times_ns'])
        reconstructed=(H*np.asarray(fit['charge_pe'])[None]).T
        total=reconstructed.sum(0) if len(reconstructed) else np.zeros(1024)
        axes[0,col].plot(GRID,y[index],color='black',lw=.8,label='Observed full local frame')
        axes[0,col].plot(GRID,10.+total,color='#d62728',lw=1,label='Hybrid fit sum')
        axes[0,col].set_title(f"({chr(97+col)}) Same parent pair, separation {gap:g} ns")
        for j in range(2):
            axes[1,col].plot(GRID,components[index,j],color=('#1f77b4','#ff7f0e')[j],lw=1,
                             label=f'True {j+1}: Q={q[index,j]:.0f} pe')
        for j,curve in enumerate(reconstructed):
            axes[1,col].plot(GRID,curve,'--',lw=1,label=f'Reconstructed {j+1}: Q={fit["charge_pe"][j]:.0f} pe')
        axes[1,col].set_title(f"({chr(99+col)}) True and reconstructed individual components")
        axes[2,col].plot(GRID,y[index]-10.-total,lw=.8,color='#995533')
        axes[2,col].set_title(f"({chr(101+col)}) Residual; true times {t[index,0]:.2f}/{t[index,1]:.2f} ns")
        for row in range(3):
            axes[row,col].set(xlim=(0,2048),xlabel='Time in local frame (ns)',ylabel='Voltage (mV)')
            if row<2:axes[row,col].legend(fontsize=7)
        selections.append(dict(index=index,gap_ns=gap,parents=info[index]['parents'],true_times_ns=t[index].tolist(),
                               true_Q_pe=q[index].tolist(),fit=fit,selection='Predefined mid_mid replicate0, gaps4/80; not selected by success.'))
    fig.suptitle('Figure 0004-39-1: Complete 2.048 us pile-up, individual components and residuals')
    fig.savefig(OUTPUT/'figures/full_frame_examples.png');plt.close(fig)
    save(OUTPUT/'examples.json',selections)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=('tables','figures','all'))
    args=parser.parse_args()
    report=make_report() if args.stage in ('tables','all') else load(OUTPUT/'resolution-report.json')
    if args.stage in ('figures','all'):figures(report)


if __name__=='__main__':
    main()
