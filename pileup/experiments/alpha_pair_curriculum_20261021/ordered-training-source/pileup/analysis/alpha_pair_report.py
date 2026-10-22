"""Tables/figures from saved single/pair results, without training or inference."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import alpha_pair_data as data
from .alpha_pair_data import PROJECT,EXPERIMENT,OUTPUT,load,read_batch,GAPS,RATIOS
from .alpha_pair_study import SEEDS,SINGLE_VARIANTS,PAIR_VARIANTS,aggregate
from .alpha_pair_fitting import PairFitter
from .alpha_rate_campaign import save


def tables():
    rows=[];by_count=[];pair_gaps=[]
    for kind in ('native','mean','linear','shared'):
        result=load(OUTPUT/'test'/kind/'results.json')
        batch=read_batch(OUTPUT/f'assets/test-{kind}.npz')
        info=load(OUTPUT/f'assets/test-{kind}.json')['cases']
        methods=list(result['rows'][0]['methods'])
        for method in methods:
            scores=[r['methods'][method]['score'] for r in result['rows']]
            strict=[r['methods'][method]['strict_score'] for r in result['rows']]
            row=dict(domain=kind,method=method,oracle=method.startswith('known_K'),**aggregate(scores))
            row['strict_recovered_count']=sum(r['recovered_count'] for r in strict)
            row['cases']=len(scores)
            row['mean_seconds']=float(np.mean([r['methods'][method]['seconds'] for r in result['rows']]))
            rows.append(row)
            for count in (0,1,2):
                ids=np.flatnonzero(batch['active'].sum(1)==count)
                selected=[scores[i] for i in ids]
                by_count.append(dict(domain=kind,method=method,count=count,cases=len(ids),
                    both_complete=sum(scores[i]['recovered_count']==2 for i in ids) if count==2 else None,
                    strict_both_complete=sum(strict[i]['recovered_count']==2 for i in ids) if count==2 else None,
                    any_extra=sum(scores[i]['unmatched_predictions']>0 for i in ids),**aggregate(selected)))
            for gap in GAPS:
                ids=[i for i,x in enumerate(info) if batch['active'][i].sum()==2 and x.get('gap_ns')==gap]
                if ids:
                    pair_gaps.append(dict(domain=kind,method=method,gap_ns=gap,cases=len(ids),
                        both_complete=sum(scores[i]['recovered_count']==2 for i in ids),
                        strict_both_complete=sum(strict[i]['recovered_count']==2 for i in ids),
                        extras=sum(scores[i]['unmatched_predictions'] for i in ids)))
    singles=load(OUTPUT/'single-test.json')['rows']
    single_rows=[]
    for r in singles:
        selected=load(OUTPUT/'models'/r['name']/'selection.json')
        single_rows.append(dict(method=r['name'],best_step=selected['step'],validation_objective=selected['objective'],
            validation_native_median_abs_q=selected['by_domain']['native']['q_error_median_abs'],**r['summary']))
    report=dict(summary=rows,by_count=by_count,pair_gaps=pair_gaps,singles=single_rows,
        scope='Native original isolated, fixed mean, linear donors and joint offline shared-cell controls are distinct',
        oracle='known_K_fit_diagnostic is separated from blind comparisons',
        strict='True-gap/4 capped at10ns, evaluation only; not applied to inference',
        source_records='24 fresh native isolated transport parents; repeated controls are not independent ROOTs')
    save(OUTPUT/'report-tables.json',report)
    return report


def figures(report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9,'axes.grid':True,'grid.alpha':.18,'savefig.dpi':150})
    out=OUTPUT/'figures';out.mkdir(exist_ok=True)
    colors={'single_logq':'#3366aa','single_relative':'#bb5555','single_fixed':'#228866',
            'pair_direct':'#3366aa','pair_physical':'#bb5555'}
    fig,axes=plt.subplots(1,3,figsize=(14,3.7),layout='constrained')
    for v in SINGLE_VARIANTS:
        for j,s in enumerate(SEEDS):
            history=load(OUTPUT/'models'/f'{v}_{s}'/'history.json')
            step=[h['step'] for h in history]
            style='-' if j==0 else '--'
            label=v.replace('single_','')+f' / seed{j+1}'
            axes[0].plot(step,[h['objective'] for h in history],style,color=colors[v],label=label)
            axes[1].plot(step,[h['parts']['logq'] for h in history],style,color=colors[v])
            axes[2].plot(step,[h['parts']['component'] for h in history],style,color=colors[v])
    for ax,title in zip(axes,('(a) Validation joint recovery objective','(b) Sampled training log-Q loss','(c) Sampled component loss')):
        ax.set(title=title,xlabel='Training step')
    axes[1].set_yscale('log');axes[2].set_yscale('log');axes[0].legend(fontsize=7,loc='lower right')
    fig.suptitle('Figure 0004-30-1: Single-pulse loss ablation, two seeds')
    fig.savefig(out/'single_training.png');plt.close(fig)

    singles=load(OUTPUT/'single-test.json')['rows'];batch=read_batch(OUTPUT/'assets/test-native-single.npz')
    trueq=batch['charge_pe'][:,0];truet=batch['times_ns'][:,0]
    fig,axes=plt.subplots(1,2,figsize=(12,3.7),layout='constrained')
    for r in singles:
        v=r['name'].rsplit('_',1)[0];j=SEEDS.index(int(r['name'].rsplit('_',1)[1]))
        q=np.array([f['charge_pe'][0] for f in r['fits']]);t=np.array([f['times_ns'][0] for f in r['fits']])
        for ax,value in zip(axes,(100*(q/trueq-1),t-truet)):
            ax.scatter(trueq,value,s=24,marker='o' if j==0 else 'x',color=colors[v],alpha=.75,label=v.replace('single_','')+f' / seed{j+1}')
    axes[0].axhline(20,color='gray',ls=':');axes[0].axhline(-20,color='gray',ls=':')
    axes[1].axhline(10,color='gray',ls=':');axes[1].axhline(-10,color='gray',ls=':')
    for ax in axes:ax.set(xlabel='True realized avalanche charge (pe)',xscale='log')
    axes[0].set(title='(a) Known-single Q calibration',ylabel='Relative Q error (%)')
    axes[1].set(title='(b) First-active-photon timing',ylabel='Time error (ns)')
    axes[0].legend(fontsize=7,ncol=2,loc='upper right')
    fig.suptitle('Figure 0004-30-2: Fresh native isolated controls; count supplied as one')
    fig.savefig(out/'single_calibration.png');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(12,3.7),layout='constrained')
    for v in PAIR_VARIANTS:
        for j,s in enumerate(SEEDS):
            h=load(OUTPUT/'models'/f'{v}_{s}'/'history.json');x=[a['step'] for a in h]
            label=v.replace('pair_','')+f' / seed{j+1}';style='-' if j==0 else '--'
            axes[0].plot(x,[a['objective'] for a in h],style,color=colors[v],label=label)
            axes[1].plot(x,[a['loss'] for a in h],style,color=colors[v])
    axes[0].set(title='(a) Validation charge-qualified F1',xlabel='Training step')
    axes[1].set(title='(b) Sampled training loss',xlabel='Training step',yscale='log')
    axes[0].legend(fontsize=8,loc='lower right');fig.suptitle('Figure 0004-31-1: Unknown 0/1/2 component training')
    fig.savefig(out/'pair_training.png');plt.close(fig)

    result=load(OUTPUT/'test/mean/results.json');info=load(OUTPUT/'assets/test-mean.json')['cases']
    methods=('pair_direct_262050001','pair_physical_262050001','enumerated_fit')
    fig,axes=plt.subplots(2,3,figsize=(14,7),layout='constrained',sharex=True,sharey=True)
    for c,m in enumerate(methods):
        for r,metric in enumerate(('score','strict_score')):
            grid=np.zeros((len(RATIOS),len(GAPS)))
            for iy,ratio in enumerate(RATIOS):
                for ix,gap in enumerate(GAPS):
                    ids=[i for i,x in enumerate(info) if x['count']==2 and x['gap_ns']==gap and x['ratio']==ratio]
                    grid[iy,ix]=np.mean([result['rows'][i]['methods'][m][metric]['recovered_count']==2 for i in ids])
            im=axes[r,c].imshow(grid,vmin=0,vmax=1,cmap='viridis',aspect='auto',origin='lower')
            axes[r,c].set(title=f'({chr(97+r*3+c)}) {m.replace("_262050001", " / seed1")} / {"10ns" if r==0 else "strict"}',
                          xticks=np.arange(len(GAPS)),xticklabels=[f'{g:g}' for g in GAPS],
                          yticks=np.arange(len(RATIOS)),yticklabels=[f'{q:g}' for q in RATIOS],xlabel='True gap (ns)',ylabel='Q2/Q1')
    fig.colorbar(im,ax=axes.ravel().tolist(),label='Both events completely recovered / four controls per cell',shrink=.8)
    fig.suptitle('Figure 0004-32-1: Independent fixed-mean controls; ordinary vs strict timing')
    fig.savefig(out/'mean_gap_ratio.png');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(12,3.8),layout='constrained')
    lines=('pair_direct_262050001','pair_direct_262050002','pair_physical_262050001','pair_physical_262050002','enumerated_fit')
    for method in lines:
        rows=[r for r in report['pair_gaps'] if r['domain']=='shared' and r['method']==method]
        rows.sort(key=lambda r:r['gap_ns'])
        x=[r['gap_ns'] for r in rows]
        for ax,key in zip(axes,('both_complete','strict_both_complete')):
            ax.plot(x,[r[key]/r['cases'] for r in rows],marker='o',label=method.replace('_26205000',' / seed'))
    for ax in axes:ax.set(xlabel='Injected first-active gap (ns)',xscale='log',ylim=(-.02,1.02),ylabel='Both complete / qualified pairs')
    axes[0].set_title('(a) Ordinary time tolerance 10ns');axes[1].set_title('(b) Strict min(10ns, gap/4)')
    axes[0].legend(fontsize=7,loc='upper left');fig.suptitle('Figure 0004-32-2: Shared-cell redigitizations of fresh optical parents')
    fig.savefig(out/'shared_gap_success.png');plt.close(fig)

    # Fixed data-order examples, never selected by reconstruction success.
    shared=read_batch(OUTPUT/'assets/test-shared.npz');sinfo=load(OUTPUT/'assets/test-shared.json')['cases']
    sr=load(OUTPUT/'test/shared/results.json')['rows']
    native=read_batch(OUTPUT/'assets/test-native.npz');nr=load(OUTPUT/'test/native/results.json')['rows']
    near=next(i for i,x in enumerate(sinfo) if shared['active'][i].sum()==2 and x.get('gap_ns',999)<=12.)
    far=next(i for i,x in enumerate(sinfo) if shared['active'][i].sum()==2 and x.get('gap_ns',0)>=80.)
    nid=int(np.flatnonzero(native['active'].sum(1)==1)[0])
    mid=next(i for i,x in enumerate(info) if x.get('gap_ns')==4. and x.get('ratio')==1.)
    examples=[('native',nid,'Native isolated'),('shared',far,'Shared cells, separated pair'),('shared',near,'Shared cells, near pair'),('mean',mid,'Mean, 4ns equal pair')]
    fig,axes=plt.subplots(4,1,figsize=(15,11),layout='constrained')
    fitter=PairFitter(load(EXPERIMENT/'settings.json'));t=np.arange(data.SAMPLES)*2.+.25
    preferred='pair_direct_262050001'
    for ax,(kind,i,label) in zip(axes,examples):
        b=read_batch(OUTPUT/f'assets/test-{kind}.npz');r=load(OUTPUT/'test'/kind/'results.json')['rows'][i]
        with np.load(OUTPUT/'test'/kind/f'{preferred}-raw.npz') as f:c=f['components_mV'][i]
        fit=r['methods'][preferred]['fit'];baseline=r['methods']['enumerated_fit']['fit']
        ax.plot(t/1000,b['voltage_mV'][i]-10.,color='.6',lw=.8,label='Observed - baseline')
        for j in np.flatnonzero(b['active'][i]):ax.plot(t/1000,b['components_mV'][i,j],ls='--',lw=1.1,label=f'True component {j+1}')
        for j in fit['slot_ids']:ax.plot(t/1000,c[j],lw=.9,label=f'Learned component {j+1}')
        hc=fitter.columns(t,baseline['times_ns']).T*np.asarray(baseline['charge_pe'])[:,None]
        ax.plot(t/1000,hc.sum(0) if len(hc) else np.zeros_like(t),color='k',ls=':',lw=.9,label='Blind fit total')
        s=r['methods'][preferred]['score'];strict=r['methods'][preferred]['strict_score']
        ax.set(title=f'({chr(97+list(axes).index(ax))}) {label}, case{i}; learned complete={s["recovered_count"]}/{s["target_count"]}, strict={strict["recovered_count"]}',ylabel='Voltage (mV)',xlabel='Local time (us)',xlim=(0,data.SAMPLES*.002))
        ax.legend(fontsize=7,ncol=5,loc='upper right')
    fig.suptitle('Figure 0004-33-1: Complete local frames and all reported learned components')
    fig.savefig(out/'full_examples.png');plt.close(fig)
    save(OUTPUT/'example-selection.json',dict(method=preferred,examples=[dict(domain=k,case=i,label=l) for k,i,l in examples],selection='First in data order within predefined physical categories; no success filter'))


def main():
    r=tables();figures(r)
    print('REPORT_READY',len(r['summary']),'rows;6 figures',flush=True)


if __name__=='__main__':main()
