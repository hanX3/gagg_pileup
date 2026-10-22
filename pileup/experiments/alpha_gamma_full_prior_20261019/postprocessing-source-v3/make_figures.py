"""Render saved scores and outcome transitions; no reconstruction or tuning."""
from collections import Counter
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--recovery-info',type=Path)
args=p.parse_args();output=args.output
summary=json.loads((output/'summary.json').read_text())
events=json.loads((output/'events.json').read_text())
names=dict(baseline='Physical\nfit',unet='U-Net',hybrid_none='Hybrid\nno penalty',
           hybrid_both='Hybrid\nQ + class',class_only='Hybrid\nclass only')
methods=summary['methods']
fig,axes=plt.subplots(1,2,figsize=(11.5,4.2))
for ax,particle in zip(axes,('alpha','gamma')):
    rows=[next(r for r in summary['pooled'] if r['method']==m and r['particle']==particle) for m in methods]
    x=np.arange(len(rows))
    ax.bar(x-.18,[100*r['full_recovery'] for r in rows],.36,label='Full recovery')
    ax.bar(x+.18,[100*r['full_precision'] for r in rows],.36,label='Full precision')
    ax.set(xticks=x,xticklabels=[names[m] for m in methods],title=particle,ylim=(0,100),ylabel='Percent')
    ax.grid(axis='y',alpha=.2)
axes[0].legend();fig.suptitle(f"{summary['windows']} new independent native windows; all methods frozen")
fig.tight_layout()
for ext in ('png','pdf'):fig.savefig(output/f'independent-test-comparison.{ext}',dpi=175,bbox_inches='tight')
plt.close(fig)

fig,axes=plt.subplots(2,2,figsize=(10,7),sharey=True)
display_methods=['baseline','hybrid_both','class_only']
colors=['#71839d','#d59439','#3b9c75']
for ax,(model,rate) in zip(axes.flat,[('brems',1000000),('brems',10000000),('mono662',1000000),('mono662',10000000)]):
    for j,(method,color) in enumerate(zip(display_methods,colors)):
        rows=[next(r for r in summary['groups'] if (r['method'],r['gamma_model'],r['rate_hz'],r['particle'])==(method,model,rate,particle)) for particle in ('alpha','gamma')]
        x=np.arange(2)+(j-1)*.25
        bars=ax.bar(x,[100*r['full_recovery'] if r['full_recovery'] is not None else 0 for r in rows],.24,color=color,label=names[method].replace('\n',' '))
        for bar,row in zip(bars,rows):
            ax.text(bar.get_x()+bar.get_width()/2,bar.get_height()+1.2,f"{row['recovered_count']}/{row['target_count']}",
                    ha='center',va='bottom',fontsize=8)
    ax.set(title=f'{model}; alpha source {rate/1e6:g} MHz',xticks=[0,1],xticklabels=['alpha','gamma'],ylim=(0,105),ylabel='Full recovery / %')
    ax.grid(axis='y',alpha=.2)
handles,labels=axes[0,0].get_legend_handles_labels()
fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=9)
fig.suptitle('Independent test: four windows per source condition; gamma source 10 MHz')
fig.tight_layout(rect=(0,.05,1,.97))
for ext in ('png','pdf'):fig.savefig(output/f'independent-test-conditions.{ext}',dpi=175,bbox_inches='tight')
plt.close(fig)

transitions=[]
for reference in ('baseline','hybrid_none','hybrid_both'):
    for particle in ('alpha','gamma'):
        subset=[e for e in events if e['particle']==particle]
        lost=[e for e in subset if e['outcomes'][reference]['outcome']=='recovered' and e['outcomes']['class_only']['outcome']!='recovered']
        gained=[e for e in subset if e['outcomes'][reference]['outcome']!='recovered' and e['outcomes']['class_only']['outcome']=='recovered']
        transitions.append(dict(reference=reference,particle=particle,lost=len(lost),gained=len(gained),
            lost_to=dict(Counter(e['outcomes']['class_only']['outcome'] for e in lost)),
            gained_from=dict(Counter(e['outcomes'][reference]['outcome'] for e in gained)),
            lost_events=[dict(stem=e['stem'],truth_index=e['truth_index']) for e in lost],
            gained_events=[dict(stem=e['stem'],truth_index=e['truth_index']) for e in gained]))
groups=[('weak_gamma_500_1000',[e for e in events if e['particle']=='gamma' and e['true_charge_pe']<1000.]),
        ('gamma_history_ge200',[e for e in events if e['particle']=='gamma' and e['preceding_voltage_mV']>=200.]),
        ('gamma_nearest_q50_lt50',[e for e in events if e['particle']=='gamma' and e['nearest_other_q50_ns'] is not None and e['nearest_other_q50_ns']<50.])]
strata=[dict(group=name,n=len(subset),outcomes={m:dict(Counter(e['outcomes'][m]['outcome'] for e in subset)) for m in methods}) for name,subset in groups]
result=dict(scope='Post-fit diagnostic on the new independent test; not used for retuning',transitions=transitions,strata=strata)
(output/'transition-diagnostics.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
timing=[];grouped_timing=[]
resumed=set(json.loads(args.recovery_info.read_text())['pending_stems']) if args.recovery_info else set()
def timing_row(method,rows):
    fits=[row['fit_seconds'] for row in rows];pipelines=[row['pipeline_seconds'] for row in rows]
    return dict(method=method,windows=len(rows),fit_median_s=float(np.median(fits)),
        pipeline_median_s=float(np.median(pipelines)),pipeline_p95_s=float(np.quantile(pipelines,.95)),
        pipeline_min_s=min(pipelines),pipeline_max_s=max(pipelines))
for method in methods:
    rows=[json.loads(path.read_text()) for path in sorted((output/'fits'/method).glob('*.json'))]
    timing.append(timing_row(method,rows))
    for label,subset in [('initial_4_workers',[r for r in rows if r['stem'] not in resumed]),
                         ('resumed_2_workers',[r for r in rows if r['stem'] in resumed])]:
        if subset:grouped_timing.append(dict(execution_group=label,**timing_row(method,subset)))
(output/'timing-summary.json').write_text(json.dumps(dict(rows=timing,groups=grouped_timing,
    context=('14 complete windows at four workers; two interrupted windows resumed at two workers; one BLAS/Torch thread each; 20 us input' if resumed else 'Observed wall time; concurrency not supplied; check execution logs'),
    limits='Includes host contention; excludes process/model startup and file export; not a hardware throughput benchmark'),indent=2)+'\n')
print(json.dumps(dict(figures=['independent-test-comparison','independent-test-conditions'],
    transitions=[{k:r[k] for k in ('reference','particle','lost','gained')} for r in transitions])),flush=True)
