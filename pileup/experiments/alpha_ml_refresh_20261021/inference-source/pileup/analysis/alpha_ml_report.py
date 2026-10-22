"""Plots and transparent failure metrics for the frozen ML-vs-fit experiment."""
from pathlib import Path
import json
import numpy as np

from .alpha_pileup_fit import score
from .alpha_rate_campaign import save


def timing_null(output, repetitions=20):
    """Destroy central event timing while retaining predicted charges/counts.

    Post-test diagnostic only, never used to tune the model. High-rate time
    coincidences and charge priors can otherwise inflate the recovery metric.
    The random shifts are not independent ROOTs and do not give statistical CI.
    """
    output=Path(output)
    summary=json.loads((output/'summary.json').read_text())
    rng=np.random.default_rng(261015822)
    results=[]
    for method in ('fit','neural_fit'):
        by_rate={r['rate_hz']:[] for r in summary['windows']}
        for row in summary['windows']:
            with np.load(row['truth_path']) as f:
                truth={k:f[k] for k in ('first_active_ns','emission_ns','charge_pe')}
            fitted=row[method]
            times=np.asarray(fitted['times_ns'],float)
            charges=np.asarray(fitted['charge_pe'],float)
            lo,hi=row['gate_ns'];central=(times>=lo)&(times<hi)
            trials=[]
            for _ in range(repetitions):
                shifted=times.copy()
                shift=float(rng.uniform(1000.,hi-lo-1000.))
                shifted[central]=lo+(times[central]-lo+shift)%(hi-lo)
                order=np.argsort(shifted)
                s=score(dict(times_ns=shifted[order],charge_pe=charges[order]),truth,row['gate_ns'])
                trials.append({k:s[k] for k in ('target_count','matched_count','recovered_count','reconstructed_count')})
            by_rate[row['rate_hz']].append(trials)
        for rate,windows in sorted(by_rate.items()):
            targets=sum(w[0]['target_count'] for w in windows)
            totals=np.array([sum(w[k]['recovered_count'] for w in windows) for k in range(repetitions)])
            results.append(dict(method=method,rate_hz=rate,windows=len(windows),target_count=targets,
                mean_recovered_count=float(totals.mean()),
                mean_recovery_fraction=float(totals.mean()/targets) if targets else None,
                trial_p16_p50_p84=([float(v) for v in np.percentile(totals/targets,[16,50,84])] if targets else None)))
    save(output/'timing-null.json',dict(seed=261015822,repetitions=repetitions,rates=results,
        interpretation='Circular shift central predictions by 1–9 us, retaining count/charges/guards; diagnostic chance-coincidence/charge-prior baseline, not independent-root confidence intervals',
        used_for_model_selection=False))
    return results


def figures(output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output=Path(output)
    summary=json.loads((output/'summary.json').read_text())
    plt.rcParams.update({'font.size':10,'savefig.dpi':170})
    x=np.arange(len(summary['fit_rates']))
    labels=[f"{r['rate_hz']/1e6:g}" for r in summary['fit_rates']]
    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for key,label,color in [('fit_rates','Blind template fit','C0'),('neural_rates','Dense CNN','C1')]:
        rates=summary[key]
        def fraction(num,den):
            return [r[num]/r[den] if r[den] else np.nan for r in rates]
        axes[0,0].plot(x,100*np.array([r['recovery_fraction'] if r['recovery_fraction'] is not None else np.nan for r in rates]),'o-',label=label,color=color)
        axes[0,1].plot(x,100*np.array(fraction('matched_count','target_count')),'o-',label=label,color=color)
        axes[1,0].plot(x,100*np.array(fraction('unmatched_predictions','reconstructed_count')),'o-',label=label,color=color)
        med=np.array([r['charge_relative_p16_p50_p84'] if r['charge_relative_p16_p50_p84'] else [np.nan]*3 for r in rates])
        axes[1,1].plot(x,100*med[:,1],'o-',label=label,color=color)
        axes[1,1].fill_between(x,100*med[:,0],100*med[:,2],alpha=.14,color=color)
    for ax,title in zip(axes.flat,['Time + charge recovery (%)','Time-matched target fraction (%)',
                                  'Unmatched / predicted events (%)','Matched charge error: median, P16–P84 (%)']):
        ax.set_title(title);ax.set_xticks(x,labels);ax.set_xlabel('Source injection rate (MHz)')
        ax.grid(alpha=.25);ax.legend()
    axes[0,0].set_ylim(-2,105);axes[0,1].set_ylim(-2,105);axes[1,0].set_ylim(-2,105)
    axes[1,1].axhline(0,color='0.4',lw=.8)
    fig.suptitle('New independent ROOTs; adequate dynamic range; q >= 500 pe, |dt| <= 10 ns, |dq/q| <= 20%')
    for suffix in ('png','pdf'):fig.savefig(output/f'ml_vs_fit_rates.{suffix}')
    plt.close(fig)
    history=json.loads((output/'training-history.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(10,3.5),constrained_layout=True)
    axes[0].plot([r['step'] for r in history],[r['loss'] for r in history],'o-')
    axes[0].set_ylabel('Training batch loss');axes[0].set_xlabel('Training step')
    axes[1].plot([r['step'] for r in history],[r['validation_objective'] for r in history],'o-')
    axes[1].axvline(summary['frozen']['step'],color='C1',ls='--',label='Frozen checkpoint')
    axes[1].set_ylabel('Validation macro charge-qualified F1');axes[1].set_xlabel('Training step');axes[1].legend()
    for ax in axes:ax.grid(alpha=.25)
    for suffix in ('png','pdf'):fig.savefig(output/f'ml_training.{suffix}')
    plt.close(fig)
    # Select examples by declared rate and label, never by performance.
    examples=[min([r for r in summary['windows'] if r['rate_hz']==rate],key=lambda r:r['label'])
              for rate in (1e6,1e7,2e7,8e7)]
    fig,axes=plt.subplots(4,1,figsize=(12,11),constrained_layout=True)
    for ax,row in zip(axes,examples):
        with np.load(row['waveform_path']) as f:t,y=f['time_ns'],f['analog_mV']
        with np.load(row['truth_path']) as f:
            qt,tt=f['charge_pe'],f['first_active_ns']
        stem=Path(row['root']).stem
        with np.load(output/'fit'/f'{stem}.npz') as f: fitted=f['model_mV']
        view=(t>=7000)&(t<8500)
        ax.plot(t[view],y[view],lw=1,color='0.4',label='Noisy complete voltage')
        ax.plot(t[view],fitted[view],lw=1,color='C0',label='Template fit voltage')
        ax.set_xlim(7000,8500);ax.set_ylabel('Voltage (mV)')
        ax.set_title(f"{row['rate_hz']/1e6:g} MHz, {row['label']}")
        second=ax.twinx()
        visible_charge=[]
        for times,charges,color,marker,label in [
            (tt[qt>=500],qt[qt>=500],'k','|','Truth realized charge'),
            (row['fit']['times_ns'],row['fit']['charge_pe'],'C0','o','Fit events'),
            (row['neural_fit']['times_ns'],row['neural_fit']['charge_pe'],'C1','x','CNN events')]:
            times,charges=np.asarray(times),np.asarray(charges)
            mask=(times>=7000)&(times<8500)
            visible_charge.extend(charges[mask].tolist())
            second.scatter(times[mask],charges[mask],s=28,color=color,marker=marker,label=label)
        second.set_ylabel('Event charge (pe)')
        second.set_ylim(0.,max(1000.,max(visible_charge,default=float(qt.max(initial=0)))*1.12))
        if not visible_charge:
            ax.text(.98,.06,'No threshold events in this fixed time view',transform=ax.transAxes,
                    ha='right',fontsize=9,color='0.3')
        ax.grid(alpha=.2)
        if row is examples[0]:
            ax.legend(loc='upper left',fontsize=8);second.legend(loc='upper right',fontsize=8)
    axes[-1].set_xlabel('Sample / arrival time (ns)')
    for suffix in ('png','pdf'):fig.savefig(output/f'ml_vs_fit_examples.{suffix}')
    plt.close(fig)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    figures(p.parse_args().output)
