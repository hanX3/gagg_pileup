"""Figures for Notebook 0004, from saved refresh scores and fitted models."""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from . import alpha_ml_refresh as study
from .alpha_full_signal_study import make_template
from .alpha_pileup_fit import score
from .alpha_rate_campaign import digest, save


LABELS = dict(baseline='Original SSE', diagonal='Current WLS', cnn='CNN', unet='U-Net',
              transformer='Transformer', gru='GRU', cnn_physical='Old CNN + physics',
              cnn_wls='CNN + residual WLS', local='Local-K SSE', correlated='Correlated GLS', **{'global':'GLS + DE'})
COLORS = dict(zip(study.METHODS, ('#8c564b', '#2864b4', '#ff7f0e', '#2ca02c', '#9467bd', '#17a0a0', '#888888', '#cf243e')))


def plotting():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': .15, 'savefig.dpi': 170})
    return plt


def finish(fig, name):
    folder = study.OUTPUT / 'figures'
    folder.mkdir(exist_ok=True)
    fig.savefig(folder / (name + '.png'), bbox_inches='tight')
    fig.savefig(folder / (name + '.svg'), bbox_inches='tight')
    plotting().close(fig)


def rates_figure(report):
    plt = plotting()
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True, layout='constrained')
    for column, phase in enumerate(('common', 'confirmation')):
        block = report['benchmarks'][phase]
        for name in study.METHODS:
            rows = sorted((r for r in block['rates'] if r['method'] == name), key=lambda r: r['rate_hz'])
            x = np.array([r['rate_hz'] / 1.e6 for r in rows])
            width = 2. if name == 'cnn_wls' else 1.
            axes[0, column].plot(x, [100 * r['recovery_fraction'] if r['recovery_fraction'] is not None else np.nan for r in rows],
                                 '-o', color=COLORS[name], lw=width, ms=3, label=LABELS[name])
            axes[1, column].plot(x, [100 * r['unmatched_predictions'] / r['target_count'] if r['target_count'] else np.nan for r in rows],
                                 '-o', color=COLORS[name], lw=width, ms=3)
        axes[0, column].set_title(f'({"a" if column == 0 else "b"}) {phase}: {block["windows"]} native windows')
        axes[1, column].set_title(f'({"c" if column == 0 else "d"}) Extra reported outputs / true targets')
        axes[1, column].set_xlabel('Source-injected alpha rate (MHz)')
        for ax in axes[:, column]:
            ax.set_xscale('log')
            ax.set_xticks([.1, 1, 10, 20, 40, 60, 80, 100], ['0.1', '1', '10', '20', '40', '60', '80', '100'],
                          rotation=55, ha='right')
    axes[0, 0].set_ylabel('Time AND charge qualified recovery (%)')
    axes[1, 0].set_ylabel('Extra outputs per 100 targets')
    axes[0, 0].legend(fontsize=8, ncol=2)
    fig.suptitle('Figure 0004-17-1: identical-window comparison; public re-evaluation and new confirmation kept separate')
    finish(fig, 'comparison')


def errors_figure(report):
    plt = plotting()
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), layout='constrained')
    block = report['benchmarks']['confirmation']
    for name in study.METHODS:
        rows = sorted((r for r in block['rates'] if r['method'] == name), key=lambda r: r['rate_hz'])
        x = [r['rate_hz'] / 1.e6 for r in rows]
        p = [r['charge_relative_p16_p50_p84'] for r in rows]
        axes[0, 0].plot(x, [100 * q[1] if q else np.nan for q in p], '-o', ms=3, color=COLORS[name], label=LABELS[name])
        axes[0, 1].plot(x, [50 * (q[2] - q[0]) if q else np.nan for q in p], '-o', ms=3, color=COLORS[name])
    for ax in axes[0]:
        ax.set_xscale('log')
        ax.set_xlabel('Source-injected rate (MHz)')
    axes[0, 0].set(title='(a) Conditional charge median bias', ylabel='Median (Qfit/Qtrue - 1), %')
    axes[0, 0].legend(fontsize=7, ncol=2)
    axes[0, 1].set(title='(b) Conditional central 68% charge half-width', ylabel='(P84 - P16) / 2, %')
    axes[0, 0].axhline(0., color='.4', ls=':')
    gains = []
    for rate in sorted({r['rate_hz'] for r in block['paired_vs_wls']}):
        rows = [r for r in block['paired_vs_wls'] if r['method'] == 'cnn_wls' and r['rate_hz'] == rate]
        gains.append((rate, sum(r['method_only'] for r in rows), sum(r['wls_only'] for r in rows)))
    x = np.arange(len(gains))
    axes[1, 0].bar(x - .18, [r[1] for r in gains], .36, label='Newly recovered vs WLS', color='#2a9d65')
    axes[1, 0].bar(x + .18, [r[2] for r in gains], .36, label='Lost vs WLS', color='#d04a4a')
    axes[1, 0].set_xticks(x, [f'{r[0] / 1.e6:g}' for r in gains])
    axes[1, 0].set(title='(c) Paired event gains and losses: CNN + residual WLS', xlabel='Source-injected rate (MHz)', ylabel='Events')
    axes[1, 0].legend(fontsize=8)
    timing = study.load(study.OUTPUT / 'confirmation/summary.json')['windows']
    timing_rates = sorted({r['rate_hz'] for r in timing})
    for name in study.METHODS:
        values = [float(np.median([r['methods'][name]['seconds'] for r in timing if r['rate_hz'] == rate])) for rate in timing_rates]
        axes[1, 1].plot(np.array(timing_rates) / 1.e6, values, '-o', ms=3, color=COLORS[name])
    axes[1, 1].set(xscale='log', yscale='log', xlabel='Source-injected rate (MHz)',
                   ylabel='Seconds / 20 us record (median at rate)', title='(d) Offline CPU time by rate, concurrent load')
    fig.suptitle('Figure 0004-17-2: matched-only charge errors, paired event recovery, and complete-record runtime')
    finish(fig, 'errors_and_pairs')


def controls_figure():
    plt = plotting()
    s = study.load(study.OUTPUT / 'controls/summary.json')
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), layout='constrained')
    for column, noise in enumerate(('white', 'correlated')):
        for name in study.METHODS:
            xs, ys = [], []
            for gap in (4., 8., 10., 12., 20., 40., 80., 200.):
                rows = [r['methods'][name]['score'] for r in s['rows'] if r['noise'] == noise and r['kind'] == 'double' and r['gap_ns'] == gap]
                xs.append(gap)
                ys.append(100 * sum(r['recovered_count'] for r in rows) / sum(r['target_count'] for r in rows))
            axes[0, column].plot(xs, ys, '-o', color=COLORS[name], ms=3, label=LABELS[name])
        axes[0, column].set(title=f'({"a" if column == 0 else "b"}) Matched mean template: {noise} noise', xscale='log',
                            xlabel='True separation (ns)', ylabel='Time AND charge qualified recovery (%)')
    axes[0, 0].legend(fontsize=7, ncol=2)
    values = []
    for name in study.METHODS:
        rows = [r['methods'][name]['score'] for r in s['rows'] if r['kind'] == 'single' and r['charges_pe'] == [40000.]]
        values.append(sum(r['recovered_count'] for r in rows) / len(rows))
    axes[1, 0].barh(range(len(values)), values, color=[COLORS[m] for m in study.METHODS])
    axes[1, 0].set_yticks(range(len(values)), [LABELS[m] for m in study.METHODS])
    axes[1, 0].set(title='(c) 40000 pe single pulse: above training range', xlabel='Fraction of two noise cases recovered', xlim=(0, 1.1))
    confirmation = study.load(study.OUTPUT / 'confirmation/summary.json')
    rows = sorted(confirmation['windows'], key=lambda r: (r['rate_hz'], r['seed']))
    x = np.arange(len(rows))
    axes[1, 1].bar(x-.18, [r['methods']['cnn_wls']['fit'].get('neural_proposals', {}).get('added_neural_components', 0) for r in rows], .36, label='Added CNN seeds')
    axes[1, 1].bar(x+.18, [r['methods']['cnn_wls']['fit'].get('neural_proposals', {}).get('excluded_edge_peaks', 0) for r in rows], .36, label='Excluded unsupported edge peaks', alpha=.8)
    axes[1, 1].set_xticks(x, [f'{r["rate_hz"] / 1.e6:g}' for r in rows], rotation=60)
    axes[1, 1].set(title='(d) Confirmation: candidate construction', xlabel='Predeclared windows, source rate (MHz)', ylabel='Candidates')
    axes[1, 1].legend(fontsize=7)
    fig.suptitle(f'Figure 0004-18-1: {len(s["rows"])} software controls; no native photon/shared-cell statistics')
    finish(fig, 'controls')


def null_figure(report):
    plt = plotting()
    block = report['benchmarks']['confirmation']
    fig, axes = plt.subplots(2, 2, figsize=(14, 8), layout='constrained')
    x = np.arange(len(study.METHODS))
    for ax, rate, panel in zip(axes.flat, (40000000, 60000000, 80000000, 100000000), 'abcd'):
        actual = [next(r for r in block['rates'] if r['method'] == m and r['rate_hz'] == rate)['recovery_fraction'] for m in study.METHODS]
        null = [next(r for r in block['timing_null'] if r['method'] == m and r['rate_hz'] == rate)['mean_recovery_fraction'] for m in study.METHODS]
        ax.bar(x - .18, np.array(actual) * 100, .36, label='Actual waveform match', color='#3876af')
        ax.bar(x + .18, np.array(null) * 100, .36, label='Five circular shifts (same window)', color='#ed9d36')
        ax.set_xticks(x, [LABELS[m] for m in study.METHODS], rotation=65, ha='right', fontsize=7)
        ax.set(title=f'({panel}) {rate / 1.e6:g} MHz, one new window', ylabel='Time AND charge qualified matching (%)')
    axes[0, 0].legend(fontsize=8)
    fig.suptitle('Figure 0004-18-2: accidental-match reference; no efficiency subtraction or exact false-discovery interpretation')
    finish(fig, 'timing_null')


def low_rate_figure(report):
    plt = plotting()
    block = report['benchmarks']['confirmation']
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), layout='constrained')
    for ax, rate, panel in zip(axes, (100000, 1000000, 10000000), 'abc'):
        rows = [next(r for r in block['rates'] if r['method'] == m and r['rate_hz'] == rate) for m in study.METHODS]
        y = np.arange(len(rows))
        recovered = np.array([r['recovered_count'] for r in rows])
        failed = np.array([r['matched_charge_failed'] for r in rows])
        missed = np.array([r['missed_count'] for r in rows])
        ax.barh(y, recovered, color='#4a9f65', label='Recovered')
        ax.barh(y, failed, left=recovered, color='#e3a43b', label='Matched, charge failed')
        ax.barh(y, missed, left=recovered + failed, color='#ba5260', label='No time match')
        ax.set_yticks(y, [LABELS[m] for m in study.METHODS])
        ax.set(title=f'({panel}) {rate / 1.e6:g} MHz, {rows[0]["target_count"]} targets', xlabel='Target events (new confirmation)')
    axes[0].legend(fontsize=8)
    fig.suptitle('Figure 0004-19-1: low-rate failures separated into detection and charge regression; small denominators retained')
    finish(fig, 'low_rates')


def figures():
    r = study.load(study.OUTPUT / 'report-tables.json')
    rates_figure(r)
    errors_figure(r)
    controls_figure()
    null_figure(r)
    low_rate_figure(r)
    print('FIGURES_READY', flush=True)


def gallery():
    from .sipm_readout import ReadoutConfig, synthesize_readout
    from .sipm_response import AvalancheResult, ELECTRON_CHARGE_C
    plt = plotting()
    settings, _, _ = study.configuration()
    template, baseline = make_template(settings), settings['fit_config']['baseline_mV']
    summary = study.load(study.OUTPUT / 'confirmation/summary.json')
    manifests = []
    for number, rate in enumerate(sorted({r['rate_hz'] for r in summary['windows']}), 1):
        row = min((r for r in summary['windows'] if r['rate_hz'] == rate), key=lambda r: r['seed'])
        stem, figure = Path(row['root']).stem, f'20-{number}'
        folder = study.OUTPUT / 'gallery' / f'rate_{rate / 1.e6:g}MHz'
        if (folder / 'manifest.json').exists():
            manifests.append(study.load(folder / 'manifest.json'))
            continue
        folder.mkdir(parents=True, exist_ok=False)
        with np.load(study.PROJECT / row['waveform']) as z:
            t, observed = z['time_ns'].copy(), z['analog_mV'].copy()
        with np.load(study.PROJECT / row['truth']) as z:
            truth = {k: z[k].copy() for k in ('first_active_ns', 'emission_ns', 'charge_pe', 'avalanche_owner',
                     'avalanche_time_ns', 'avalanche_charge_C', 'signal_mV')}
        metadata = study.load(study.OUTPUT / 'inputs/metadata' / (stem + '.json'))
        q0 = metadata['config']['sensor']['gain_electrons'] * ELECTRON_CHARGE_C
        readout = replace(ReadoutConfig(**metadata['config']['readout']), baseline_mV=0., noise_rms_mV=0.)
        duration = metadata['duration_ns']
        owner = truth['avalanche_owner']
        order = np.argsort(owner, kind='stable')
        records = np.empty(len(order), dtype=[('time_ns', 'f8'), ('charge_C', 'f8')])
        records['time_ns'] = truth['avalanche_time_ns'][order]
        records['charge_C'] = truth['avalanche_charge_C'][order]
        bounds = np.r_[0, np.cumsum(np.bincount(owner + 1, minlength=len(truth['charge_pe']) + 1))]
        del order

        def voltage(block):
            out = synthesize_readout(AvalancheResult(block, duration, q0, np.array([]), {}), readout, 0)
            np.testing.assert_array_equal(out.arrays['time_ns'], t)
            return out.arrays['signal_mV']

        actual_indices = np.flatnonzero(truth['charge_pe'] > 0.)
        actual = np.empty((len(actual_indices), len(t)))
        ancestry_error = 0.
        for k, i in enumerate(actual_indices):
            block = records[bounds[i + 1]:bounds[i + 2]]
            actual[k] = voltage(block)
            ancestry_error = max(ancestry_error, abs(float(block['charge_C'].sum() / q0 - truth['charge_pe'][i])))
        dark = voltage(records[bounds[0]:bounds[1]])
        native_error = float(np.max(abs(actual.sum(axis=0) + dark - truth['signal_mV'])))
        assert native_error < 1.e-8 and ancestry_error < 1.e-7
        models = {}
        for name in ('cnn_wls', 'diagonal'):
            with np.load(study.OUTPUT / 'confirmation/models' / (stem + '_' + name + '.npz')) as z:
                models[name] = z['model_mV'].copy()
        data = row['methods']['cnn_wls']
        fit = data['fit']
        full = fit.get('all_components', fit)
        times, peaks = np.asarray(full['times_ns']), np.asarray(full['peak_mV'])
        components = template(t[None, :] - times[:, None]) * peaks[:, None]
        model_error = float(np.max(abs(baseline + components.sum(axis=0) - models['cnn_wls'])))
        assert model_error < 1.e-8
        assert json.loads(json.dumps(score(fit, truth, row['gate_ns']))) == data['score']
        report_index = {(tt, aa): j for j, (tt, aa) in enumerate(zip(fit['times_ns'], fit['peak_mV']))}
        fitted_reports = np.array([report_index.get((tt, aa), -1) for tt, aa in zip(times, peaks)])
        matched = {j: i for i, j in data['score']['all_matched_pairs']}
        fitted_truth = np.array([matched.get(int(j), -1) for j in fitted_reports])
        arrays = dict(time_ns=t, observed_mV=observed, actual_components_mV=actual,
                      actual_truth_indices=actual_indices, dark_mV=dark, fitted_components_mV=components,
                      fitted_reported_indices=fitted_reports, fitted_truth_indices=fitted_truth,
                      hybrid_model_mV=models['cnn_wls'], wls_model_mV=models['diagonal'])
        np.savez_compressed(folder / 'components.npz', **arrays)
        fig, axes = plt.subplots(3, 1, figsize=(17, 10), sharex=True, layout='constrained',
                                 gridspec_kw={'height_ratios': [1., 1.3, .65]})
        x = t / 1000.
        axes[0].plot(x, observed, color='.25', lw=.7, label='Measured voltage')
        axes[0].plot(x, models['diagonal'], color=COLORS['diagonal'], lw=.9, label='Current WLS')
        axes[0].plot(x, models['cnn_wls'], color=COLORS['cnn_wls'], lw=.9, label='CNN + residual WLS')
        axes[0].set(title='(a) Complete 0-20 us observed record and blind total fits', ylabel='Total (mV)')
        axes[0].legend(ncol=3, fontsize=8)

        def color(i):
            return plt.get_cmap('hsv')((int(i) * .618033988749895) % 1.)

        for k, i in enumerate(actual_indices):
            axes[1].plot(x, actual[k], color=color(i) if truth['charge_pe'][i] >= 500. else '.7', ls='--', lw=.55, alpha=.8)
        axes[1].plot(x, dark, color='.45', ls='--', lw=.45)
        for k, (j, i) in enumerate(zip(fitted_reports, fitted_truth)):
            c, style = ('.65', ':') if j < 0 else (color(i), '-') if i >= 0 else ('#cf0084', ':')
            axes[1].plot(x, components[k], color=c, ls=style, lw=.6, alpha=.8)
        from matplotlib.lines import Line2D
        axes[1].legend(handles=[Line2D([], [], color='.3', ls='--', label='Actual alpha (same-color match)'),
                               Line2D([], [], color='.3', label='CNN + residual WLS component'),
                               Line2D([], [], color='#cf0084', ls=':', label='Unmatched reported fit'),
                               Line2D([], [], color='.65', ls=':', label='Internal subthreshold term')], fontsize=8, ncol=2)
        axes[1].set(title=f'(b) All {len(actual_indices)} nonzero alpha contributions and all {len(times)} hybrid fit terms', ylabel='Components (mV)')
        low, high = axes[1].get_ylim()
        axes[1].set_ylim(low, high * 1.18 if high > 0. else 1.)
        for name in ('diagonal', 'cnn_wls'):
            axes[2].plot(x, observed - models[name], color=COLORS[name], lw=.65, label=LABELS[name])
        axes[2].axhline(0., color='.4', ls=':')
        axes[2].set(title='(c) Measured minus total fit', xlabel='Record time (us)', ylabel='Residual (mV)')
        axes[2].legend(fontsize=8, ncol=2)
        for ax in axes:
            ax.axvspan(0., 5., color='.9', alpha=.5, zorder=-10)
            ax.axvspan(15., 20., color='.9', alpha=.5, zorder=-10)
            ax.set_xlim(0., 20.)
        fig.suptitle(f'Figure 0004-{figure}: {rate / 1.e6:g} MHz, first predeclared new window, seed {row["seed"]}\n'
                     'Every recorded sample and all alpha/fit components; performance scored only in central 5-15 us', fontsize=12)
        fig.savefig(folder / 'full_components.png')
        fig.savefig(folder / 'full_components.svg')
        plt.close(fig)
        manifest = dict(figure_id=figure, rate_hz=rate, root=row['root'], seed=row['seed'],
                        sample_count=len(t), actual_nonzero_alpha_count=len(actual_indices), internal_fit_count=len(times),
                        gate_ns=row['gate_ns'], display_ns=[0., 20000.], score=data['score'],
                        validation=dict(actual_sum_max_error_mV=native_error, ancestry_charge_max_error_pe=ancestry_error,
                                        model_sum_max_error_mV=model_error, score_reproduced=True),
                        waveform_sha256=digest(study.PROJECT / row['waveform']), truth_sha256=digest(study.PROJECT / row['truth']),
                        selection='First predeclared seed; retain empty record, no outcome-based selection',
                        truth_use='Post-inference display and scoring only')
        save(folder / 'manifest.json', manifest)
        manifests.append(manifest)
        print('GALLERY', figure, rate / 1.e6, len(actual_indices), len(times), flush=True)
    save(study.OUTPUT / 'gallery/index.json', dict(phase='confirmation', figures=manifests))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['figures', 'gallery'])
    args = parser.parse_args()
    globals()[args.stage]()
