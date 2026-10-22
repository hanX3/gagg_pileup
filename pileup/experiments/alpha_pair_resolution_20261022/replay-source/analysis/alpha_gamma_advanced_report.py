"""Read-only reporting and explicit single-model failure controls.

No method/window is selected here. The choice is read from the development
freeze, and every native test method/window remains in the saved tables.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .alpha_gamma_advanced_fit import AdvancedCalibration, METHODS, fit_advanced_prefix
from .alpha_gamma_advanced_study import OLD, POSTS, PRE, group, load_results, summary
from .alpha_gamma_short_window import PrefixCalibration
from .alpha_gamma_short_window_study import save

PROJECT = Path(__file__).resolve().parents[1]
LABELS = {
    'legacy_gls': 'Original calibration + GLS',
    'expanded_gls': 'Expanded calibration + GLS',
    'conditional_mle': 'Conditional Gaussian',
    'pca8_mle': 'PCA rank 8',
    'pca16_mle': 'PCA rank 16',
    'bootstrap_mle': 'Conditional + calibration error',
    'pca8_bootstrap_mle': 'PCA 8 + calibration error',
}
COLORS = dict(zip(METHODS, ('#737373', '#d17c23', '#007a9d', '#8b4ca1', '#598d35', '#c13e52', '#92734d')))


def failure_controls(output):
    """Mean-additive doubles test assumption failure, not native pileup recovery."""
    output = Path(output).resolve()
    frozen = json.loads((output/'freeze.json').read_text())
    methods = tuple(dict.fromkeys(('legacy_gls', 'expanded_gls', frozen['primary_method'])))
    calibration = AdvancedCalibration.load(output/'calibration')
    legacy = PrefixCalibration.load(PROJECT/OLD/'calibration')
    rng = np.random.default_rng(262020919)
    t = np.arange(0., 20000., 2.)
    rows = []
    for first, second in (('alpha', 'gamma'), ('gamma', 'alpha')):
        for gap in (50., 100., 200., 400.):
            for ratio in (.25, 1., 4.):
                q1, q2, start = 2000., 2000.*ratio, 5000.
                voltage = (calibration.meta['baseline_mV']
                           + q1*calibration.template(t-start, first, q1)
                           + q2*calibration.template(t-start-gap, second, q2)
                           + rng.normal(0., calibration.meta['electronics_sigma_mV'], len(t)))
                for post in (150., 300.):
                    fitted = fit_advanced_prefix(t, voltage, calibration, legacy,
                                                  pre_ns=PRE, post_ns=post, methods=methods)
                    assert fitted['detected']
                    for method, result in fitted['methods'].items():
                        result['relative_first_charge_error'] = result['charge_pe']/q1-1
                        result['first_type_correct'] = result['particle'] == first
                    rows.append({'first': first, 'second': second, 'gap_ns': gap,
                                 'first_charge_pe': q1, 'second_charge_pe': q2, 'ratio': ratio,
                                 'post_ns': post,
                                 'second_signal_can_enter_prefix': fitted['last_observed_ns'] >= start+gap,
                                 'fit': fitted})
    counts = []
    for method in methods:
        for post in (150., 300.):
            for gap in (50., 100., 200., 400.):
                selected = [r for r in rows if r['post_ns'] == post and r['gap_ns'] == gap]
                fits = [r['fit']['methods'][method] for r in selected]
                counts.append({'method': method, 'post_ns': post, 'gap_ns': gap, 'prefixes': len(fits),
                               'first_charge_5pct_success': sum(abs(f['relative_first_charge_error']) <= .05 for f in fits),
                               'first_joint_success': sum(abs(f['relative_first_charge_error']) <= .05
                                                          and f['first_type_correct'] for f in fits),
                               'abs_first_error_p90': float(np.quantile([abs(f['relative_first_charge_error']) for f in fits], .9)),
                               'second_enters': sum(r['second_signal_can_enter_prefix'] for r in selected)})
    result = {'generator': 'Add two independent calibrated mean voltages and white electronics noise',
              'native_shared_cell_double': False, 'pulse_count_fitted': False,
              'scored_quantity': 'First pulse Q/type only; second pulse is not reconstructed',
              'limits': 'No shot-noise realization, native nonlinear cell sharing, rate scan or blind double recovery',
              'seed': 262020919, 'independent_synthetic_waveforms': 24,
              'rows': rows, 'counts': counts}
    save(output/'controlled-doubles.json', result)
    return result


def make_report(output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    output = Path(output).resolve()
    settings = json.loads((output/'freeze.json').read_text())
    results = load_results(output/'test')
    s = summary(results)
    assert len(results) == 168 and s == json.loads((output/'test-summary.json').read_text())
    primary = settings['primary_method']
    controls = json.loads((output/'controlled-doubles.json').read_text())
    cal = json.loads((output/'calibration/calibration.json').read_text())
    old = json.loads((PROJECT/OLD/'calibration/calibration.json').read_text())
    figures = []
    highlights = tuple(dict.fromkeys(('legacy_gls', 'expanded_gls', primary)))

    def style(method):
        return dict(color=COLORS[method], label=LABELS[method]+(' [selected]' if method == primary else ''),
                    lw=2.4 if method == primary else 1.6 if method in highlights else 1.,
                    alpha=1. if method in highlights else .6,
                    ls='-' if method in highlights else '--', marker='o', ms=3)

    with plt.rc_context({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False}):
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True, layout='constrained')
        for col, particle in enumerate(('alpha', 'gamma')):
            for method in METHODS:
                groups = [group(s, method, particle, post) for post in POSTS]
                for row, key in enumerate(('charge_5pct_fraction', 'joint_fraction')):
                    axes[row, col].plot(POSTS, [100*g[key] if g[key] is not None else np.nan for g in groups], **style(method))
                if method == primary and all(g['joint_wilson95'] is not None for g in groups):
                    band = np.array([g['joint_wilson95'] for g in groups], float)*100
                    axes[1, col].fill_between(POSTS, band[:, 0], band[:, 1], color=COLORS[method], alpha=.1)
            axes[0, col].set_title(f"{particle}: {groups[0]['targets']} targets, Q >= 500 pe")
            for row in (0, 1):
                axes[row, col].axhline(95., color='.4', lw=.8, ls=':')
                axes[row, col].set(ylim=(-3, 103), xticks=POSTS)
                axes[row, col].grid(alpha=.2)
            axes[1, col].set_xlabel('T after observed trigger (ns); pre-trigger 20 ns')
        axes[0, 0].set_ylabel('Charge within 5% / all targets (%)')
        axes[1, 0].set_ylabel('Charge within 5% AND correct type (%)')
        axes[0, 0].legend(fontsize=7, loc='lower right')
        fig.suptitle('Independent native isolated test; method selected on development only')
        fig.savefig(output/'advanced_success.png', dpi=170)
        plt.close(fig); figures.append('advanced_success.png')

        fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex='col', layout='constrained')
        for col, particle in enumerate(('alpha', 'gamma')):
            for method in highlights:
                groups = [group(s, method, particle, post) for post in POSTS]
                center = np.array([g['charge_bias'] for g in groups], float)*100
                spread = np.array([g['charge_error_sd'] for g in groups], float)*100
                axes[0, col].errorbar(POSTS, center, yerr=spread, capsize=3, **style(method))
                axes[1, col].plot(POSTS, [100*g['interval_coverage_all_targets'] if g['targets'] else np.nan for g in groups], **style(method))
            axes[0, col].axhspan(-5., 5., color='.9', alpha=.5)
            axes[0, col].axhline(0., color='.4', lw=.8)
            axes[0, col].set_title(particle)
            axes[1, col].axhline(95., color='.4', lw=.8, ls=':')
            axes[1, col].set(ylim=(-3, 103))
            axes[1, col].set_xlabel('T after observed trigger (ns)')
            for ax in axes[:, col]:
                ax.set_xticks(POSTS); ax.grid(alpha=.2)
        axes[0, 0].set_ylabel('Matched charge error: mean +/- SD (%)')
        axes[1, 0].set_ylabel('Local nominal 95% interval covers Q / all targets (%)')
        axes[0, 0].legend(fontsize=7)
        fig.suptitle('Spread is NOT an error bar on the mean; local intervals have no coverage guarantee')
        fig.savefig(output/'advanced_error_coverage.png', dpi=170)
        plt.close(fig); figures.append('advanced_error_coverage.png')

        fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True, layout='constrained')
        for ax, post in zip(axes, (150., 300.)):
            for method in highlights:
                groups = [r for r in controls['counts'] if r['method'] == method and r['post_ns'] == post]
                ax.plot([g['gap_ns'] for g in groups], [100*g['first_joint_success']/g['prefixes'] for g in groups], **style(method))
            ax.set(title=f'T after trigger = {post:g} ns', xlabel='True second-minus-first onset (ns)',
                   xticks=(50, 100, 200, 400), ylim=(-3, 103))
            ax.grid(alpha=.2)
        axes[0].set_ylabel('First pulse Q within 5% AND correct type (%)')
        axes[0].legend(fontsize=7)
        fig.suptitle('Synthetic mean-additive doubles: assumption-failure control, not native pileup recovery')
        fig.savefig(output/'advanced_double_control.png', dpi=170)
        plt.close(fig); figures.append('advanced_double_control.png')

    earliest = []
    for method in METHODS:
        passing = [post for i, post in enumerate(POSTS) if all(
            group(s, method, particle, longer)['targets'] > 0
            and group(s, method, particle, longer)['joint_fraction'] >= .95
            for particle in ('alpha', 'gamma') for longer in POSTS[i:])]
        earliest.append({'method': method, 'post_ns': passing[0] if passing else None,
                         'scope': 'Exploratory test grid description only; no re-selection or population certificate'})
    diagnostic = []
    for post in POSTS:
        for method in METHODS:
            for particle in ('alpha', 'gamma'):
                g = group(s, method, particle, post)
                diagnostic.append({key: g[key] for key in (
                    'post_ns', 'method', 'particle', 'targets', 'detected', 'correct_type',
                    'charge_5pct_success', 'joint_charge_type_success', 'charge_bias', 'charge_error_sd',
                    'charge_rmse', 'abs_charge_error_p90', 'interval_valid_count', 'interval_covered_count',
                    'interval_coverage_all_targets', 'optimizer_failure_count', 'onset_bound_count', 'charge_bound_count')}
                    | {'missed': g['targets']-g['detected'], 'wrong_type': g['detected']-g['correct_type']})
    report = {
        'primary_method': primary, 'primary_short_goal_ns': settings['primary_short_goal_ns'],
        'development_selected_95pct_point_window_ns': settings['selected_95pct_point_window_ns'],
        'native_test_events': len(results), 'truth_counts': s['truth_counts'],
        'calibration': {'old_total': old['training_rows'], 'old_used': old['used_training_rows'],
                        'old_classes': old['classes'], 'expanded_total': cal['training_rows'],
                        'expanded_used': cal['used_training_rows'], 'expanded_classes': cal['classes']},
        'exploratory_earliest_test_grid': earliest, 'diagnostic_rows': diagnostic,
        'paired_comparisons': json.loads((output/'paired-comparisons.json').read_text()),
        'false_triggers': s['false_triggers'], 'controlled_doubles': controls['counts'],
        'figures': figures,
        'runtime': {'total_worker_seconds_all_methods_all_windows': sum(r['seconds'] for r in results),
                    'median_worker_seconds_per_event': float(np.median([r['seconds'] for r in results])),
                    'scope': 'Offline Python benchmark including all seven methods/five windows; not an online throughput limit'},
        'limits': s['limits']+[
            'Same G4/SiPM model for calibration and test; sparse energies/positions, not hardware validation',
            'PCA modes replace correlated covariance, not additional unrestricted pulse parameters',
            'Bootstrap estimates finite bin-mean error only; bins/covariance/systematics are not fully resampled',
            'Mean-template uncertainty is shared across events; eventwise intervals are not independent',
            'At most one pulse, fixed double-emission physics; no free-tau minimum window or source-rate capability'],
    }
    save(output/'report.json', report)
    print('REPORTED_ADVANCED', primary, s['truth_counts'], flush=True)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=('controls', 'report'))
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    (failure_controls if args.stage == 'controls' else make_report)(args.output)


if __name__ == '__main__':
    main()
