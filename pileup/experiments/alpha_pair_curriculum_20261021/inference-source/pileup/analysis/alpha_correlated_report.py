"""Postprocessing only: frozen correlated-fit study (no inference or selection).

Run from the project root after all predeclared evaluations have completed.
Saved figures carry explicit Notebook 0003 section-19 identifiers. Reading this
module/its figures does not require PyTorch or access to the original ROOTs.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np

from .alpha_correlated_fit import InnovationsMetric, ResidualCalibration
from .alpha_correlated_study import PROJECT, OUTPUT
from .alpha_rate_campaign import save


METHODS = ('baseline', 'local', 'diagonal', 'correlated', 'global', 'cnn', 'gru')
LABELS = {'baseline': 'Original SSE', 'local': 'Local-K SSE', 'diagonal': 'Diagonal WLS',
          'correlated': 'Correlated GLS', 'global': 'GLS + local DE', 'cnn': 'Frozen CNN', 'gru': 'Light GRU'}
COLORS = dict(zip(METHODS, ('#777777', '#1f77b4', '#17becf', '#d62728', '#9467bd', '#ff7f0e', '#2ca02c')))
COUNTS = ('target_count', 'central_injected', 'central_nonzero_charge', 'reconstructed_count',
          'matched_count', 'recovered_count', 'missed_count', 'unmatched_predictions',
          'merge_flags', 'split_flags')


def good_indices(scored):
    return {pair[0] for pair, error in zip(scored['matched_pairs'], scored['relative_charge_errors'])
            if abs(error) <= .2}


def statistics(values):
    a = np.asarray(values, float)
    return dict(n=len(a), mean=float(a.mean()) if len(a) else None,
                sd=float(a.std(ddof=1)) if len(a) > 1 else None,
                p16_p50_p84=np.percentile(a, [16, 50, 84]).tolist() if len(a) else None)


def tables():
    test = json.loads((OUTPUT/'test-summary.json').read_text())
    assert test['test_windows'] == 23 and set(test['rates']) == set(METHODS)
    controls = json.loads((OUTPUT/'controls-summary.json').read_text())
    optimizer = json.loads((OUTPUT/'optimizer-controls.json').read_text())
    shifted = json.loads((OUTPUT/'time-shift-controls.json').read_text())
    rates = sorted({r['rate_hz'] for r in test['windows']})
    native = []
    for rate in rates:
        windows = [r for r in test['windows'] if r['rate_hz'] == rate]
        for name in METHODS:
            scores = [w['methods'][name]['score'] for w in windows]
            totals = {k: sum(s[k] for s in scores) for k in COUNTS}
            good, lost, common = 0, 0, 0
            for w in windows:
                reference = good_indices(w['methods']['local']['score'])
                candidate = good_indices(w['methods'][name]['score'])
                good += len(candidate-reference); lost += len(reference-candidate)
                common += len(reference & candidate)
            null = []
            for delta in shifted['shifts_ns']:
                sc = [r['score'] for r in shifted['rows'] if r['rate_hz'] == rate
                      and r['method'] == name and r['shift_ns'] == delta]
                assert len(sc) == len(windows)
                n = sum(s['target_count'] for s in sc)
                null.append(sum(s['recovered_count'] for s in sc)/n if n else 0.)
            native.append(dict(rate_hz=rate, method=name, windows=len(windows), **totals,
                fit_status_counts=dict(Counter(w['methods'][name]['fit']['status'] for w in windows)),
                recovery_fraction=totals['recovered_count']/totals['target_count'] if totals['target_count'] else None,
                charge=statistics([e for s in scores for e in s['relative_charge_errors']]),
                time_ns=statistics([e for s in scores for e in s['time_errors_ns']]),
                seconds=statistics([w['methods'][name]['seconds'] for w in windows]),
                gained_vs_local=good, lost_vs_local=lost, common_vs_local=common,
                shifted_fractions=null, shifted_fraction_mean=float(np.mean(null))))
    aggregate = []
    for name in METHODS:
        rows = [r for r in native if r['method'] == name]
        counts = {k: sum(r[k] for r in rows) for k in COUNTS}
        aggregate.append(dict(method=name, **counts,
            recovery_fraction=counts['recovered_count']/max(1, counts['target_count']),
            summed_inference_seconds=sum(w['methods'][name]['seconds'] for w in test['windows']),
            gained_vs_local=sum(r['gained_vs_local'] for r in rows),
            lost_vs_local=sum(r['lost_vs_local'] for r in rows)))
    software = []
    for noise in ('white', 'correlated'):
        for kind in ('empty', 'single', 'double'):
            for name in METHODS[:5]:
                rows = [r for r in controls['rows'] if r['noise'] == noise and r['kind'] == kind]
                totals = {k: sum(r['scores'][name][k] for r in rows) for k in COUNTS}
                software.append(dict(noise=noise, kind=kind, method=name, windows=len(rows), **totals))
    conditional = []
    for k in (2, 3):
        rows = [r for r in optimizer['rows'] if r['known_count'] == k]
        for name in ('slsqp', 'population_plus_slsqp'):
            complete = [np.all(abs(np.asarray(r['results'][name]['time_errors_ns'])) <= 10)
                        and np.all(abs(np.asarray(r['results'][name]['relative_charge_errors'])) <= .2) for r in rows]
            conditional.append(dict(known_count=k, method=name, cases=len(rows), all_components_complete=int(sum(complete)),
                total_evaluations=sum(r['results'][name]['evaluations'] for r in rows),
                total_seconds=sum(r['results'][name]['seconds'] for r in rows)))
    diagnostics = dict(native=native, aggregate=aggregate, software=software, known_count=conditional,
        candidate_limits=dict(initial_candidate_guard=500, maximum_accepted_additions=160,
            initial_candidate_limit_windows=sum(w['methods']['baseline']['fit']['status']=='candidate_limit' for w in test['windows']),
            initial_candidate_counts=[dict(root=w['root'], rate_hz=w['rate_hz'], candidates=w['methods']['baseline']['fit']['candidates']) for w in test['windows']]),
        known_count_objective_improved=sum(r['results']['population_plus_slsqp']['cost'] < r['results']['slsqp']['cost']-1e-6
                                            for r in optimizer['rows']),
        interpretation='Descriptive fixed-cohort comparison; no new-test model selection, confidence interval, energy calibration or native resolution claim')
    save(OUTPUT/'report-tables.json', diagnostics)
    return test, controls, optimizer, diagnostics


def plot_style():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': .22,
                         'figure.dpi': 110, 'savefig.dpi': 150})
    return plt


def finish(fig, number, title, filename):
    fig.suptitle(f'Figure 0003-19-{number}: {title}', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, .96))
    fig.savefig(OUTPUT/filename)
    import matplotlib.pyplot as plt
    plt.close(fig)


def calibration_figure():
    plt = plot_style()
    cal = ResidualCalibration(**json.loads((OUTPUT/'calibration/covariance.json').read_text())['parameters'])
    with np.load(OUTPUT/'calibration/residuals.npz') as z:
        mu, residual = z['mean'].copy(), z['residual'].copy()
    variance = cal.electronic_sigma_mV**2 + cal.dark_variance_mV2 + cal.shot_coefficient_mV*np.maximum(mu, 0)
    standardized = residual/np.sqrt(variance)
    white = np.stack([InnovationsMetric(m, cal).apply(r) for m, r in zip(mu, residual)])
    # Global residual moments, matching calibration's zero-mean model; no
    # per-record mean subtraction that could hide baseline mismatch.
    def acf(x):
        return np.array([np.sum(x*x) if lag == 0 else np.sum(x[:, lag:]*x[:, :-lag])
                         for lag in range(61)])/np.sum(x*x)
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    flat_mu, r2 = mu.ravel(), residual.ravel()**2
    edges = np.unique(np.quantile(flat_mu, np.linspace(0, 1, 45)))
    centers, moments = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        chosen = (flat_mu >= lo) & (flat_mu <= hi)
        if chosen.sum() >= 20:
            centers.append(float(flat_mu[chosen].mean())); moments.append(float(r2[chosen].mean()))
    ax[0].plot(centers, moments, 'o', ms=3, label='Single-window residual second moment')
    xx = np.linspace(0, flat_mu.max(), 200)
    ax[0].plot(xx, cal.electronic_sigma_mV**2+cal.dark_variance_mV2+cal.shot_coefficient_mV*xx,
               label='Frozen linear variance model')
    ax[0].set(xlabel='Refitted single-pulse mean (mV)', ylabel='Residual second moment (mV$^2$)', title='(a) Calibration only')
    ax[0].legend(fontsize=7)
    lags = np.arange(61)
    scaled_latent = np.sqrt(np.maximum(variance-cal.electronic_sigma_mV**2, 0)/variance)
    predicted_acf = np.array([1. if lag == 0 else cal.rho**lag*np.sum(
        scaled_latent[:, lag:]*scaled_latent[:, :-lag])/scaled_latent.size for lag in lags])
    ax[1].plot(2*lags, acf(standardized), label='Variance-standardized residual')
    ax[1].plot(2*lags, predicted_acf, '--', label='Frozen total covariance prediction')
    ax[1].set(xlabel='Lag (ns)', ylabel='Normalized uncentered correlation', title='(b) Correlation is an approximation')
    ax[1].legend(fontsize=7)
    ax[2].plot(2*lags[1:], acf(white)[1:], label='After full innovation whitening')
    ax[2].axhline(0, color='k', lw=.7)
    ax[2].set(xlabel='Lag (ns)', ylabel='Normalized uncentered correlation', title='(c) Remaining model mismatch')
    ax[2].legend(fontsize=7)
    save(OUTPUT/'calibration/residual-diagnostics.json', dict(
        raw_standardized_acf_lag1=float(acf(standardized)[1]),
        predicted_standardized_acf_lag1=float(predicted_acf[1]),
        whitened_acf_lag1=float(acf(white)[1]),
        whitened_mean=float(white.mean()), whitened_second_moment=float(np.mean(white**2)),
        note='Calibration diagnostic only. AR1 and linear variance are approximations; no test-dependent covariance refit.'))
    finish(fig, 1, 'Empirical single-pulse residual covariance', 'covariance_calibration.png')


def native_figures(report):
    plt = plot_style(); rows = report['native']
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    for name in METHODS:
        selected = [r for r in rows if r['method'] == name]
        x = np.array([r['rate_hz']/1e6 for r in selected])
        series = [100*np.array([r['recovery_fraction'] or 0 for r in selected]),
                  100*np.array([r['unmatched_predictions']/max(1, r['reconstructed_count']) for r in selected]),
                  100*np.array([r['missed_count']/max(1, r['target_count']) for r in selected]),
                  100*np.array([(r['matched_count']-r['recovered_count'])/max(1, r['matched_count']) for r in selected])]
        for a, y in zip(ax.ravel(), series):
            a.plot(x, y, '-o', color=COLORS[name], label=LABELS[name], ms=4, lw=1.3)
    titles = ['(a) Complete recovery: time and charge', '(b) Unmatched outputs / central outputs',
              '(c) No time match / target events', '(d) Charge failure / time-matched targets']
    for a, title in zip(ax.ravel(), titles):
        a.set(xscale='log', xlabel='Injected alpha source rate (MHz)', ylabel='Percent', title=title)
        a.set_ylim(-2, 102)
    ax[0, 0].legend(ncol=2, fontsize=8)
    finish(fig, 2, 'New native test windows, seven frozen methods', 'native_recovery.png')

    fig, ax = plt.subplots(3, 2, figsize=(14, 13))
    for name in METHODS:
        selected = [r for r in rows if r['method'] == name]
        x = np.array([r['rate_hz']/1e6 for r in selected])
        med = [100*r['charge']['p16_p50_p84'][1] if r['charge']['n'] else np.nan for r in selected]
        sd = [100*r['charge']['sd'] if r['charge']['sd'] is not None else np.nan for r in selected]
        time_sd = [r['time_ns']['sd'] if r['time_ns']['sd'] is not None else np.nan for r in selected]
        seconds = [r['seconds']['p16_p50_p84'][1] for r in selected]
        for a, y in zip(ax.ravel()[:4], (med, sd, time_sd, seconds)):
            a.plot(x, y, '-o', color=COLORS[name], label=LABELS[name], ms=4)
    for a, title, ylabel in zip(ax.ravel()[:4],
        ('(a) Conditional charge bias', '(b) Conditional charge spread', '(c) Conditional time spread', '(d) Full 20 us inference'),
        ('Median charge error (%)', 'Sample SD charge error (%)', 'Sample SD time error (ns)', 'Median seconds / window')):
        a.set(xscale='log', xlabel='Injected alpha source rate (MHz)', ylabel=ylabel, title=title)
    ax[1, 1].set_yscale('log'); ax[0, 0].legend(ncol=2, fontsize=7)
    totals = report['aggregate']; xx = np.arange(len(METHODS)); width = .38
    ax[2, 0].bar(xx-width/2, [r['gained_vs_local'] for r in totals], width, label='Newly complete vs Local-K SSE', color='#2ca02c')
    ax[2, 0].bar(xx+width/2, [-r['lost_vs_local'] for r in totals], width, label='Lost vs Local-K SSE', color='#d62728')
    ax[2, 0].set(xticks=xx, xticklabels=[LABELS[m] for m in METHODS], ylabel='Target count', title='(e) Paired event identity, all 23 windows')
    ax[2, 0].tick_params(axis='x', rotation=30); ax[2, 0].legend(fontsize=7)
    selected = [next(r for r in rows if r['rate_hz'] == 1e8 and r['method'] == m) for m in METHODS]
    ax[2, 1].bar(xx-width/2, [100*r['recovery_fraction'] for r in selected], width, label='Actual', color='#1f77b4')
    ax[2, 1].bar(xx+width/2, [100*r['shifted_fraction_mean'] for r in selected], width, label='Mean of 5 circular shifts', color='#aaaaaa')
    ax[2, 1].set(xticks=xx, xticklabels=[LABELS[m] for m in METHODS], ylabel='Complete matches (%)', title='(f) 100 MHz chance-matching reference')
    ax[2, 1].tick_params(axis='x', rotation=30); ax[2, 1].legend(fontsize=7)
    finish(fig, 3, 'Conditional errors, computation and matched-event controls', 'native_diagnostics.png')


def controls_figure(controls, optimizer):
    plt = plot_style(); fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    for a, noise in zip(ax[0], ('white', 'correlated')):
        for name in METHODS[:5]:
            gaps, efficiency = [], []
            for gap in (4., 8., 12., 20., 40., 80., 200.):
                scores = [r['scores'][name] for r in controls['rows'] if r['noise'] == noise and r['kind'] == 'double' and r['gap_ns'] == gap]
                gaps.append(gap); efficiency.append(100*sum(s['recovered_count'] for s in scores)/sum(s['target_count'] for s in scores))
            a.plot(gaps, efficiency, '-o', ms=4, color=COLORS[name], label=LABELS[name])
        a.set(xscale='log', xlabel='Injected double separation (ns)', ylabel='Complete recovery (%)', ylim=(-2, 102),
              title=f'({"a" if noise == "white" else "b"}) Blind, matched template, {noise} noise')
        a.legend(fontsize=7)
    markers = {2: 'o', 3: '^'}
    for k in (2, 3):
        rows = [r for r in optimizer['rows'] if r['known_count'] == k]
        old = np.array([r['results']['slsqp']['cost'] for r in rows])
        new = np.array([r['results']['population_plus_slsqp']['cost'] for r in rows])
        ax[1, 0].scatter(old, new, marker=markers[k], s=28, label=f'Known K={k}')
    limits = ax[1, 0].get_xlim(); ax[1, 0].plot(limits, limits, 'k--', lw=.8)
    ax[1, 0].set(xlabel='Local-search weighted residual cost', ylabel='DE + local-search weighted residual cost',
                 title='(c) Known-count diagnostic: lower cost is not enough')
    ax[1, 0].legend()
    for offset, name in enumerate(('slsqp', 'population_plus_slsqp')):
        values = []
        for k in (2, 3):
            rows = [r for r in optimizer['rows'] if r['known_count'] == k]
            values.append(sum(np.all(abs(np.asarray(r['results'][name]['time_errors_ns'])) <= 10)
                              and np.all(abs(np.asarray(r['results'][name]['relative_charge_errors'])) <= .2) for r in rows))
        ax[1, 1].bar(np.arange(2)+(offset-.5)*.35, values, .35, label=name.replace('_', ' '))
    ax[1, 1].set(xticks=[0, 1], xticklabels=['Known K=2', 'Known K=3'], ylim=(0, 19),
                 ylabel='All components complete / 18 cases', title='(d) Conditional software diagnostic; extra DE compute')
    ax[1, 1].legend(fontsize=7)
    finish(fig, 4, 'Software controls: blind double pulses and known-count searches', 'software_controls.png')


def examples_figure(test):
    plt = plot_style(); fig, ax = plt.subplots(2, 2, figsize=(14, 9))
    for i, rate in enumerate((10000000, 40000000)):
        row = min((r for r in test['windows'] if r['rate_hz'] == rate), key=lambda r: r['seed'])
        stem = Path(row['root']).stem
        with np.load(PROJECT/row['waveform']) as z: t, y = z['time_ns'].copy(), z['analog_mV'].copy()
        with np.load(PROJECT/row['truth']) as z:
            tt, qt = z['first_active_ns'].copy(), z['charge_pe'].copy()
        select = (t >= 5000) & (t <= 5800)
        ax[i, 0].plot(t[select]/1000, y[select], color='k', lw=1, alpha=.6, label='Observed voltage')
        for name in ('baseline', 'local', 'correlated', 'global'):
            model_path = OUTPUT/'test/windows'/(stem+'_'+name+'.npz') if name in ('baseline', 'local') else OUTPUT/'test'/name/(stem+'.npz')
            with np.load(model_path) as z: model = z['model_mV'].copy()
            ax[i, 0].plot(t[select]/1000, model[select], color=COLORS[name], lw=1, label=LABELS[name])
        ax[i, 0].set(xlabel='Record time (us)', ylabel='Analog voltage (mV)', title=f'({chr(97+2*i)}) {rate/1e6:g} MHz, first predeclared test')
        ax[i, 0].legend(fontsize=7, ncol=2)
        select = (tt >= 5000) & (tt <= 5800) & (qt >= 500)
        ax[i, 1].scatter(tt[select]/1000, qt[select], marker='x', color='k', s=35, label='Scoring truth, Q >= 500 pe')
        for j, name in enumerate(METHODS):
            fit = row['methods'][name]['fit']; times = np.asarray(fit['times_ns']); charge = np.asarray(fit['charge_pe'])
            select = (times >= 5000) & (times <= 5800)
            style = dict(color=COLORS[name], s=26, label=LABELS[name])
            if j in (0, 1, 2, 4, 5):
                style.update(facecolors='none', edgecolors=COLORS[name])
            ax[i, 1].scatter(times[select]/1000, charge[select],
                             marker=('o', '^', 's', '+', 'v', 'd', '*')[j], **style)
        ax[i, 1].set(xlabel='First active photon / fitted onset (us)', ylabel='Realized / fitted charge (pe)',
                    title=f'({chr(98+2*i)}) {rate/1e6:g} MHz; truth only after inference')
        ax[i, 0].set_xlim(5., 5.8)
        ax[i, 1].set_xlim(5., 5.8)
        ax[i, 1].legend(fontsize=6, ncol=2)
    finish(fig, 5, 'Fixed 5.0-5.8 us examples, no success-based selection', 'native_examples.png')


def main():
    test, controls, optimizer, report = tables()
    calibration_figure(); native_figures(report); controls_figure(controls, optimizer); examples_figure(test)
    print('REPORT', {r['method']: (r['recovered_count'], r['target_count'], r['unmatched_predictions']) for r in report['aggregate']})


if __name__ == '__main__':
    main()
