"""Read-only summaries and figures for the frozen native prefix experiment."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from .alpha_gamma_short_window_study import METHODS, POSTS, PRES, adaptive, load_results, save, summarize


def selected_group(summary, method, particle, pre, post, band='all_above_500'):
    return next(r for r in summary['groups'] if r['method']==method and r['particle']==particle
                and r['pre_ns']==pre and r['post_ns']==post and r['band']==band)


def make_report(output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    output=Path(output)
    results=load_results(output/'test')
    summary=summarize(results);save(output/'test-summary.json',summary)
    settings=json.loads((output/'frozen-settings.json').read_text())
    adaptive_result=adaptive(results,settings);save(output/'adaptive-test.json',adaptive_result)
    colors={'physics_ols':'#737373','physics_gls':'#ce7b22','template_ols':'#2875aa','template_gls':'#b64046'}
    labels={'physics_ols':'Fixed 95/351, OLS','physics_gls':'Fixed 95/351, GLS',
            'template_ols':'Native class templates, OLS','template_gls':'Native class templates, GLS'}
    figures=[]
    with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
        fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True,sharey=True,layout='constrained')
        for row,pre in enumerate((20.,50.)):
            for col,particle in enumerate(('alpha','gamma')):
                ax=axes[row,col]
                for method in METHODS:
                    group=[selected_group(summary,method,particle,pre,t) for t in POSTS]
                    fraction=[100*r['charge_5pct_fraction'] if r['charge_5pct_fraction'] is not None else np.nan for r in group]
                    ax.plot(POSTS,fraction,'o-',ms=4,color=colors[method],label=labels[method])
                    if method=='template_gls' and all(r['charge_5pct_wilson95'] for r in group):
                        band=np.array([r['charge_5pct_wilson95'] for r in group])*100
                        ax.fill_between(POSTS,band[:,0],band[:,1],color=colors[method],alpha=.1,label='GLS Wilson 95%')
                n=group[0]['targets']
                ax.axhline(95,color='0.4',ls='--',lw=1)
                ax.set(title=f'{particle}: {n} targets; pre-trigger {pre:g} ns',xscale='log',ylim=(-3,104))
                ax.set_xticks([50,100,200,300,500,1000],labels=['50','100','200','300','500','1000'])
                ax.grid(alpha=.2)
                if col==0:ax.set_ylabel('Unconditional charge-within-5% success (%)')
                if row==1:ax.set_xlabel('T after observed trigger (ns)')
        axes[0,0].legend(fontsize=7,loc='lower right')
        fig.suptitle('New native isolated test events; Q >= 500 pe; same full analog input for every method')
        for ext in ('png','pdf'):fig.savefig(output/f'prefix_success.{ext}',dpi=170)
        plt.close(fig);figures.append('prefix_success.png')

        fig,axes=plt.subplots(2,2,figsize=(12,7),sharex=True,layout='constrained')
        for col,particle in enumerate(('alpha','gamma')):
            for method in METHODS:
                rows=[selected_group(summary,method,particle,50.,t) for t in POSTS]
                axes[0,col].plot(POSTS,[100*r['charge_bias'] if r['charge_bias'] is not None else np.nan for r in rows],
                                 'o-',ms=3,label=labels[method],color=colors[method])
                axes[1,col].plot(POSTS,[100*r['abs_charge_error_p90'] if r['abs_charge_error_p90'] is not None else np.nan for r in rows],
                                 'o-',ms=3,label=labels[method],color=colors[method])
            axes[0,col].axhspan(-5,5,color='0.85',alpha=.5)
            axes[0,col].axhline(0,color='0.5',lw=.8)
            axes[1,col].axhline(5,color='0.5',lw=.8,ls='--')
            axes[0,col].set_title(particle)
            for row in range(2):
                axes[row,col].set(xscale='log');axes[row,col].grid(alpha=.2)
            axes[1,col].set_xlabel('T after observed trigger (ns)')
        axes[0,0].set_ylabel('Mean relative charge error (%)')
        axes[1,0].set_ylabel('90th percentile absolute error (%)')
        axes[0,0].legend(fontsize=7)
        fig.suptitle('Matched-signal error distribution, pre-trigger 50 ns; misses counted separately in success plots')
        for ext in ('png','pdf'):fig.savefig(output/f'prefix_errors.{ext}',dpi=170)
        plt.close(fig);figures.append('prefix_errors.png')

        fig,axes=plt.subplots(1,2,figsize=(13,4),layout='constrained')
        for ax,particle in zip(axes,('alpha','gamma')):
            grid=np.array([[selected_group(summary,'template_gls',particle,p,t)['charge_5pct_fraction'] for t in POSTS] for p in PRES],float)
            im=ax.imshow(grid*100,origin='lower',aspect='auto',vmin=0,vmax=100,cmap='viridis')
            for i,p in enumerate(PRES):
                for j,t in enumerate(POSTS):
                    r=selected_group(summary,'template_gls',particle,p,t)
                    ax.text(j,i,f"{r['charge_5pct_success']}/{r['targets']}",ha='center',va='center',fontsize=7,color='white' if grid[i,j]<.65 else 'black')
            ax.set(xticks=range(len(POSTS)),xticklabels=[str(int(t)) for t in POSTS],
                   yticks=range(len(PRES)),yticklabels=[str(int(t)) for t in PRES],
                   xlabel='T after observed trigger (ns)',ylabel='Pre-trigger span (ns)',title=particle)
        fig.colorbar(im,ax=axes,label='Charge-within-5% success (%)',shrink=.9)
        fig.suptitle('Native class-template GLS: every tested window; cells are counts, not independent samples')
        for ext in ('png','pdf'):fig.savefig(output/f'prefix_window_map.{ext}',dpi=170)
        plt.close(fig);figures.append('prefix_window_map.png')

    # Describe all grid points; do not retune the frozen method on these curves.
    earliest=[]
    for method in METHODS:
        for pre in PRES:
            passing=[]
            for post in POSTS:
                rows=[selected_group(summary,method,particle,pre,t) for particle in ('alpha','gamma') for t in POSTS if t>=post]
                if all(r['targets']>0 and r['charge_5pct_fraction']>=.95 for r in rows):passing.append(post)
            earliest.append({'method':method,'pre_ns':pre,'earliest_post_ns':min(passing) if passing else None,
                             'criterion':'Exploratory point fraction >=95% per class at this and all longer tested windows, not certified population threshold'})
    false_triggers={}
    for pre in PRES:
        zeros=[r for r in results if r['truth_for_evaluation_only']['charge_pe']==0]
        false_triggers[str(int(pre))]={'zero_charge_windows':len(zeros),
            'triggered_windows':sum(next(f for f in r['prefixes'] if f['pre_ns']==pre)['detected'] for r in zeros),
            'predicted_above_500_at_300ns':sum(
                bool(f['detected'] and f['methods']['template_gls']['charge_pe']>=500)
                for r in zeros for f in r['prefixes'] if f['pre_ns']==pre and f['post_ns']==300.)}
    secondary_tolerance=[]
    for pre in (20.,50.):
        for post in POSTS:
            for particle in ('alpha','gamma'):
                events=[r for r in results if r['truth_for_evaluation_only']['particle']==particle
                        and r['truth_for_evaluation_only']['charge_pe']>=500]
                successes=0;correct=0
                for event in events:
                    fit=next(f for f in event['prefixes'] if f['pre_ns']==pre and f['post_ns']==post)
                    if not fit['evaluation']['matched_signal']:continue
                    r=fit['methods']['template_gls'];ok=abs(r['evaluation']['relative_charge_error'])<=.10
                    successes+=ok;correct+=ok and r['evaluation']['correct_type']
                secondary_tolerance.append({'pre_ns':pre,'post_ns':post,'particle':particle,'targets':len(events),
                    'charge_10pct_success':successes,'joint_charge10_type_success':correct,
                    'scope':'Secondary descriptive tolerance; not the primary predeclared5% target'})
    details=[]
    for event in results:
        truth=event['truth_for_evaluation_only']
        if truth['charge_pe']<500:continue
        for post in (100.,300.):
            fit=next(f for f in event['prefixes'] if f['pre_ns']==50. and f['post_ns']==post)
            if fit['detected']:
                model=fit['methods']['template_gls']
                details.append({'root':event['row']['root'],'particle':truth['particle'],'true_charge_pe':truth['charge_pe'],
                    'post_ns':post,'predicted_charge_pe':model['charge_pe'],'predicted_type':model['particle'],
                    'relative_charge_error':model['evaluation']['relative_charge_error'],
                    'correct_type':model['evaluation']['correct_type'],'matched':fit['evaluation']['matched_signal'],
                    'onset_bound':model['onset_at_search_bound']})
    report={'events':summary['independent_events'],'truth_counts':summary['truth_counts'],
        'earliest_grid_points':earliest,'false_triggers_on_no_signal':false_triggers,
        'secondary_charge_10pct':secondary_tolerance,
        'adaptive':adaptive_result['counts'],'figures':figures,'diagnostic_rows':details,
        'development_choice':json.loads((output/'development-window-choice.json').read_text()),
        'notes':['All errors refer to realized light-descendant avalanche Q, not source/deposited energy',
                 'Four methods share causal analog input; historical two-event ADC scan is not a matched performance comparison',
                 'Wilson bands are per-condition finite-event intervals; no multiple-window or sequential coverage guarantee',
                 'At most one pulse is assumed; no claim of mixed/high-rate or free-lifetime recovery',
                 'T=75 ns uses samples through74 ns on a2 ns grid; span is data span, not CPU runtime or source-onset delay']}
    save(output/'report.json',report)
    print(json.dumps({'events':report['events'],'counts':report['truth_counts'],'adaptive':report['adaptive']},ensure_ascii=False))
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();make_report(a.output)


if __name__=='__main__':main()
