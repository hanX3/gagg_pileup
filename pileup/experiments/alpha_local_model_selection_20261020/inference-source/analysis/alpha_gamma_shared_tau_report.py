"""Describe all frozen-test outcomes without choosing parameters on the test."""
import json
from pathlib import Path

import numpy as np

from .alpha_gamma_shared_tau import ProfileProblem
from .alpha_gamma_shared_tau_study import METHODS, POSTS, load_results
from .alpha_gamma_short_window_study import save

LABELS = dict(fixed_tau='Fixed 95/351 ns', shared_tau='Shared calibration taus',
              constrained_tau='Calibration-constrained taus', free_tau='Free event taus',
              template_gls='Calibrated full-template GLS')
COLORS = dict(zip(METHODS, ('#777777', '#2466a8', '#ba7b18', '#b34567', '#25855f')))


def group(summary, method, particle, post):
    return next(g for g in summary['groups'] if g['method'] == method and g['particle'] == particle
                and g['post_ns'] == post and g.get('band', 'all_above_500') == 'all_above_500')


def earliest(summary, method, metric):
    for i, post in enumerate(POSTS):
        groups = [group(summary, method, c, t) for c in ('alpha', 'gamma') for t in POSTS[i:]]
        if all(g['targets'] and g[metric]/g['targets'] >= .95 for g in groups):
            return post
    return None


def run(output):
    output = Path(output).resolve()
    summary = json.loads((output/'test-summary.json').read_text())
    shared = json.loads((output/'shared-fit.json').read_text())
    frozen = json.loads((output/'freeze.json').read_text())
    matched = json.loads((output/'matched-model-summary.json').read_text())
    results = load_results(output/'test')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.titlesize': 11})

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    for j, c in enumerate(('alpha', 'gamma')):
        for k, metric in enumerate(('charge_5pct_success', 'joint_charge_type_success')):
            ax = axes[k, j]
            for method in METHODS:
                groups = [group(summary, method, c, t) for t in POSTS]
                ax.plot(POSTS, [100*g[metric]/g['targets'] if g['targets'] else np.nan for g in groups],
                        'o-', ms=4, color=COLORS[method], label=LABELS[method])
            ax.axhline(95, color='black', ls=':', lw=1)
            ax.set(title=f"({chr(97+2*k+j)}) {c}: "+('charge within5%' if k == 0 else 'charge within5% and correct type'),
                   ylabel='All target events / %', xlabel='Time after observed trigger / ns', ylim=(-3, 103))
            ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=8, loc='lower right')
    fig.suptitle('Figure12-1. Independent native test; pre-trigger20 ns, Q >=500 pe')
    fig.savefig(output/'shared_tau_success.png', dpi=170); plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    for j, c in enumerate(('alpha', 'gamma')):
        selected = [r for r in results if r['truth_for_evaluation_only']['particle'] == c
                    and r['truth_for_evaluation_only']['charge_pe'] >= 500.]
        for k, (field, physical, learned) in enumerate((('tau_fast_ns', 95., shared['taus_ns'][0]),
                                                      ('tau_slow_ns', 351., shared['taus_ns'][1]))):
            ax = axes[k, j]
            ax.axhline(physical, color='gray', ls='--', label='G4 emission constant')
            ax.axhline(learned, color=COLORS['shared_tau'], ls=':', label='Independent shared fit')
            for method in ('constrained_tau', 'free_tau'):
                values = []
                for post in POSTS:
                    fits = [next(p for p in r['prefixes'] if p['post_ns'] == post) for r in selected]
                    sample = [p['methods'][method][field] for p in fits if p['evaluation']['matched_signal']]
                    values.append(np.percentile(sample, [25, 50, 75]) if sample else [np.nan]*3)
                values = np.array(values)
                ax.plot(POSTS, values[:, 1], 'o-', color=COLORS[method], label=LABELS[method])
                ax.fill_between(POSTS, values[:, 0], values[:, 2], color=COLORS[method], alpha=.15)
            ax.set(title=f'({chr(97+2*k+j)}) {c}: '+('fast' if k == 0 else 'slow'),
                   ylabel='Effective response lifetime / ns', xlabel='Time after observed trigger / ns')
            ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle('Figure12-2. Native lifetime median and interquartile range; model mismatch retained')
    fig.savefig(output/'shared_tau_lifetimes.png', dpi=170); plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    for j, c in enumerate(('alpha', 'gamma')):
        for k, metric in enumerate(('charge_5pct_success', 'both_tau_10pct_success')):
            ax = axes[k, j]
            methods = METHODS[:-1] if k == 0 else ('constrained_tau', 'free_tau')
            for method in methods:
                groups = [group(matched, method, c, t) for t in POSTS]
                ax.plot(POSTS, [100*g[metric]/g['targets'] for g in groups], 'o-',
                        color=COLORS[method], label=LABELS[method])
            ax.axhline(95, color='black', ls=':', lw=1)
            ax.set(title=f'({chr(97+2*k+j)}) {c}: '+('charge within5%' if k == 0 else 'both event taus within10%'),
                   ylabel='Independent matched events / %', xlabel='Time after preset observation origin / ns', ylim=(-3, 103))
            ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=8, loc='lower right')
    axes[1, 0].legend(fontsize=8, loc='lower right')
    fig.suptitle('Figure12-3. Exactly matched mean + white electronics only; excludes native shot fluctuations')
    fig.savefig(output/'shared_tau_matched.png', dpi=170); plt.close(fig)

    project = output.parents[2]
    calibration_meta = json.loads((project/'jupyter/artifacts/alpha_gamma_advanced_fit_20261020/calibration/calibration.json').read_text())
    chosen = {}
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    for j, c in enumerate(('alpha', 'gamma')):
        r = next(r for r in results if r['truth_for_evaluation_only']['particle'] == c
                 and r['truth_for_evaluation_only']['charge_pe'] >= 500.
                 and all(p['evaluation']['matched_signal'] for p in r['prefixes']))
        chosen[c] = r['row']['root']
        with np.load(project/r['row']['waveform']) as wave:
            t, y = wave['time_ns'], wave['analog_mV']
        for k, post in enumerate((100., 300.)):
            prefix = next(p for p in r['prefixes'] if p['post_ns'] == post)
            x = t-prefix['trigger']['time_ns']; use = (x >= -20.) & (x <= 1000.)
            x, values = x[use], y[use]
            ax = axes[k, j]
            ax.plot(x, values, color='gray', lw=.7, alpha=.4, label='Observed (tail for display only)')
            ax.plot(x[x <= post], values[x <= post], color='black', lw=1., label='Fitted prefix')
            for method in METHODS[:-1]:
                fit = prefix['methods'][method]
                p = ProfileProblem(x, values, calibration_meta['readout'], calibration_meta['area_per_pe_mV_ns'], 10., .2, .1)
                prediction = fit['baseline_mV']+p.columns((fit['tau_fast_ns'], fit['tau_slow_ns']),
                                                         fit['onset_relative_ns'])@fit['amplitudes_pe']
                ax.plot(x[x <= post], prediction[x <= post], color=COLORS[method], lw=1.1, label=LABELS[method])
                ax.plot(x[x > post], prediction[x > post], color=COLORS[method], ls='--', lw=1.)
            ax.axvline(post, color='black', ls=':', lw=.8)
            ax.set(title=f'({chr(97+2*k+j)}) {c}, post{post:g} ns', xlabel='Time relative to observed trigger / ns',
                   ylabel='Analog voltage / mV', xlim=(-20, 1000))
            ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle('Figure12-4. Prefix fits and extrapolation; dashed tails never used in inference')
    fig.savefig(output/'shared_tau_examples.png', dpi=170); plt.close(fig)

    rng = np.random.default_rng(262022001)
    paired = []
    for c in ('alpha', 'gamma'):
        selected = [r for r in results if r['truth_for_evaluation_only']['particle'] == c
                    and r['truth_for_evaluation_only']['charge_pe'] >= 500.]
        for post in (150., 300.):
            for method in ('shared_tau', 'constrained_tau', 'free_tau'):
                for control in ('fixed_tau', 'template_gls'):
                    values = []
                    for r in selected:
                        p = next(p for p in r['prefixes'] if p['post_ns'] == post)
                        def success(m):
                            return int(p['evaluation']['matched_signal'] and abs(p['methods'][m]['evaluation']['relative_charge_error']) <= .05
                                       and p['methods'][m]['evaluation']['correct_type'])
                        values.append(success(method)-success(control))
                    draws = np.mean(rng.choice(values, size=(2000, len(values)), replace=True), axis=1)
                    paired.append(dict(particle=c, post_ns=post, method=method, control=control,
                        difference=float(np.mean(values)), paired_bootstrap95=np.percentile(draws, [2.5, 97.5]).tolist()))
    report = dict(taus_ns=shared['taus_ns'], log_tau_covariance=shared['log_tau_covariance'],
        selected_prior_inflation=frozen['prior_inflation'], native_test=summary['truth_counts'],
        earliest_native=[dict(method=m, charge5_post_ns=earliest(summary, m, 'charge_5pct_success'),
                             joint_post_ns=earliest(summary, m, 'joint_charge_type_success')) for m in METHODS],
        earliest_matched=[dict(method=m, charge5_post_ns=earliest(matched, m, 'charge_5pct_success'),
                              both_tau10_post_ns=earliest(matched, m, 'both_tau_10pct_success')) for m in METHODS[:-1]],
        paired_comparisons=paired, example_roots=chosen,
        figures=['shared_tau_success.png', 'shared_tau_lifetimes.png', 'shared_tau_matched.png', 'shared_tau_examples.png'],
        interpretation='Shared/constrained taus use independent calibration. Agreement with95/351 is microscopic only in matched model; native fitted taus are effective response parameters. No universal minimum or blind-pileup result.')
    save(output/'report.json', report)


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args(); run(a.output)
