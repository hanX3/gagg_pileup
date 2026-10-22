"""Common performance plots and time-shift diagnosis for five methods."""
from pathlib import Path
import json
import numpy as np
from .alpha_pileup_fit import score
from .alpha_rate_campaign import save

LABELS={'fit':'Blind template fit','cnn':'Frozen CNN','unet':'Small U-Net',
        'transformer':'Local Transformer','cnn_physical':'CNN + physical refit'}


def timing_null(folder,repetitions=20):
    folder=Path(folder);s=json.loads((folder/'summary.json').read_text())
    rows=[]
    for method in s['method_order']:
        # Identical circular shifts for every method on each ROOT.
        rng=np.random.default_rng(261016900);groups={r['rate_hz']:[] for r in s['windows']}
        for row in s['windows']:
            with np.load(row['truth_path']) as f:truth={k:f[k] for k in ('first_active_ns','emission_ns','charge_pe')}
            fitted=row['methods'][method]['fit'];times=np.asarray(fitted['times_ns']);q=np.asarray(fitted['charge_pe'])
            lo,hi=row['gate_ns'];central=(times>=lo)&(times<hi);trials=[]
            for _ in range(repetitions):
                shifted=times.copy();shift=float(rng.uniform(1000.,hi-lo-1000.))
                shifted[central]=lo+(times[central]-lo+shift)%(hi-lo);order=np.argsort(shifted)
                result=score(dict(times_ns=shifted[order],charge_pe=q[order]),truth,row['gate_ns'])
                trials.append(dict(target_count=result['target_count'],recovered_count=result['recovered_count']))
            groups[row['rate_hz']].append(trials)
        for rate,windows in sorted(groups.items()):
            targets=sum(w[0]['target_count'] for w in windows)
            values=np.array([sum(w[j]['recovered_count'] for w in windows) for j in range(repetitions)])
            rows.append(dict(method=method,rate_hz=rate,target_count=targets,
                mean_recovery_fraction=float(values.mean()/targets) if targets else None,
                trial_p16_p50_p84=np.percentile(values/targets,[16,50,84]).tolist() if targets else None))
        print('TIMING NULL',s['benchmark'],method,flush=True)
    save(folder/'timing-null.json',dict(seed=261016900,repetitions=repetitions,rates=rows,
        used_for_model_selection=False,interpretation='20 circular shifts of the same ROOT are a chance-coincidence diagnostic, not independent-root confidence intervals'))


def figures(output,benchmark='new_test'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output=Path(output);folder=output/benchmark
    s=json.loads((folder/'summary.json').read_text());methods=s['method_order']
    rates=s['rates']['fit'];x=np.arange(len(rates));labels=[f"{r['rate_hz']/1e6:g}" for r in rates]
    plt.rcParams.update({'font.size':10,'savefig.dpi':170})
    fig,axes=plt.subplots(2,2,figsize=(13,8.5),constrained_layout=True)
    for method in methods:
        data=s['rates'][method]
        axes[0,0].plot(x,[100*r['recovery_fraction'] if r['recovery_fraction'] is not None else np.nan for r in data],'o-',label=LABELS[method])
        axes[0,1].plot(x,[100*r['matched_count']/r['target_count'] if r['target_count'] else np.nan for r in data],'o-',label=LABELS[method])
        axes[1,0].plot(x,[100*r['unmatched_predictions']/r['reconstructed_count'] if r['reconstructed_count'] else 0 for r in data],'o-',label=LABELS[method])
        axes[1,1].plot(x,[100*r['charge_relative_p16_p50_p84'][1] if r['charge_relative_p16_p50_p84'] else np.nan for r in data],'o-',label=LABELS[method])
    for ax,title in zip(axes.flat,['Time + charge success (%)','Time-matched targets (%)',
                                  'Unmatched / predictions (%)','Conditional charge error: median (%)']):
        ax.set_title(title);ax.set_xticks(x,labels);ax.set_xlabel('Source injection rate (MHz)');ax.grid(alpha=.25);ax.legend(fontsize=8)
    for ax in (axes[0,0],axes[0,1],axes[1,0]):ax.set_ylim(-2,105)
    axes[1,1].axhline(0,color='0.4',lw=.8)
    fig.suptitle(f"{benchmark}: common independent ROOTs, q >= 500 pe, |dt| <= 10 ns, |dq/q| <= 20%")
    for suffix in ('png','pdf'):fig.savefig(folder/f'multimethod_rates.{suffix}')
    plt.close(fig)
    if (folder/'timing-null.json').exists():
        null=json.loads((folder/'timing-null.json').read_text())
        lookup={(r['method'],r['rate_hz']):r for r in null['rates']}
        fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
        for method in methods:
            actual=np.array([r['recovery_fraction'] if r['recovery_fraction'] is not None else np.nan for r in s['rates'][method]])
            chance=np.array([lookup[method,r['rate_hz']]['mean_recovery_fraction']
                             if lookup[method,r['rate_hz']]['mean_recovery_fraction'] is not None else np.nan for r in rates],dtype=float)
            axes[0].plot(x,100*chance,'o-',label=LABELS[method])
            axes[1].plot(x,100*(actual-chance),'o-',label=LABELS[method])
        axes[0].set_title('Time-shifted success (%)');axes[1].set_title('Actual minus shifted success (percentage points)')
        for ax in axes:ax.set_xticks(x,labels);ax.set_xlabel('Source injection rate (MHz)');ax.grid(alpha=.25);ax.legend(fontsize=8)
        for suffix in ('png','pdf'):fig.savefig(folder/f'multimethod_timing_null.{suffix}')
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,3.8),constrained_layout=True)
    for kind in ('unet','transformer'):
        history=json.loads((output/f'models/{kind}/history.json').read_text())
        axes[0].plot([r['step'] for r in history],[r['loss'] for r in history],label=LABELS[kind])
        axes[1].plot([r['step'] for r in history],[r['validation_objective'] for r in history],label=LABELS[kind])
    axes[0].set_ylabel('Training batch loss');axes[1].set_ylabel('Validation macro charge-qualified F1')
    for ax in axes:ax.set_xlabel('Training step');ax.grid(alpha=.25);ax.legend()
    for suffix in ('png','pdf'):fig.savefig(output/f'multimethod_training.{suffix}')
    plt.close(fig)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path)
    p.add_argument('--benchmark',choices=('new_test','regression'),default='new_test')
    a=p.parse_args();timing_null(a.output/a.benchmark);figures(a.output,a.benchmark)
