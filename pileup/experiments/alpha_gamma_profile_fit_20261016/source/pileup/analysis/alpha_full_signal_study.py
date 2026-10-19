"""Blind rate comparison assuming adequate acquisition dynamic range.

Run from the project root. Every declared window uses its saved, noisy,
pre-quantization total voltage. ADC codes and clipping flags are not inputs.
Existing ROOT files, acquisition configurations and historical results are read-only.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from datetime import datetime
import importlib.metadata
import json
from pathlib import Path
import platform
import time

import numpy as np

from .alpha_pileup_fit import AlphaTemplate, FitConfig, fit_waveform, score
from .alpha_rate_campaign import digest, save


def make_template(settings):
    t = settings['template']
    return AlphaTemplate(t['readout'], t['slow_fraction'], t['tau_fast_ns'], t['tau_slow_ns'])


def fit_window(task):
    source, original, settings, destination = task
    source, destination = Path(source), Path(destination)
    stem = Path(original['root']).stem
    paths = {name: source/name/f'{stem}.{suffix}' for name, suffix in
             [('waveforms', 'npz'), ('truth', 'npz'), ('metadata', 'json')]}
    metadata = json.loads(paths['metadata'].read_text())
    if metadata['root_sha256'] != original['root_sha256']:
        raise ValueError(f'ROOT identity mismatch: {stem}')
    with np.load(paths['waveforms']) as waveform:
        t = waveform['time_ns']
        voltage = waveform['analog_mV']
    template, config = make_template(settings), FitConfig(**settings['fit_config'])
    started = time.monotonic()
    fitted = fit_waveform(t, voltage, template, config)
    seconds = time.monotonic() - started
    model = fitted.pop('model_mV')
    np.savez_compressed(destination/'fits'/f'{stem}.npz', model_mV=model)

    # Pulse ancestry, count, times and charges are opened only after blind fitting.
    with np.load(paths['truth']) as labels:
        truth = {k: labels[k] for k in ('first_active_ns', 'emission_ns', 'charge_pe')}
    scored = score(fitted, truth, metadata['gate_ns'])
    result = dict(label=original['label'], root=original['root'],
                  root_sha256=original['root_sha256'], rate_hz=original['rate_hz'],
                  source=str(source), input_field='analog_mV', gate_ns=metadata['gate_ns'],
                  seed=metadata['seed'], fit=fitted, score=scored, fit_seconds=seconds,
                  input_samples=len(t), input_min_mV=float(voltage.min()),
                  input_max_mV=float(voltage.max()),
                  input_sha256={str(p): digest(p) for p in paths.values()})
    save(destination/'fits'/f'{stem}.json', result)
    print(f"FULL SIGNAL {original['label']}: {fitted['status']}, "
          f"{scored['recovered_count']}/{scored['target_count']}, {seconds:.1f} s", flush=True)
    return result


def summarize(rows):
    rates = []
    for rate in sorted({r['rate_hz'] for r in rows}):
        selected = [r for r in rows if r['rate_hz'] == rate]
        scores = [r['score'] for r in selected]
        totals = {key: sum(s[key] for s in scores) for key in (
            'target_count', 'central_injected', 'central_nonzero_charge', 'reconstructed_count',
            'matched_count', 'recovered_count', 'missed_count', 'unmatched_predictions',
            'merge_flags', 'split_flags')}
        def percentiles(key):
            values = [v for s in scores for v in s[key]]
            return np.percentile(values, [16, 50, 84]).tolist() if values else None
        exposure = sum(s['duration_ns'] for s in scores)
        rates.append(dict(rate_hz=rate, windows=len(selected), gate_time_ns=exposure,
            fit_status_counts=dict(Counter(r['fit']['status'] for r in selected)),
            recovery_fraction=totals['recovered_count']/totals['target_count'] if totals['target_count'] else None,
            matched_charge_failed=totals['matched_count']-totals['recovered_count'],
            charge_relative_p16_p50_p84=percentiles('relative_charge_errors'),
            time_error_ns_p16_p50_p84=percentiles('time_errors_ns'),
            charge_producing_rate_hz=totals['central_nonzero_charge']/(exposure*1e-9),
            target_rate_hz=totals['target_count']/(exposure*1e-9), **totals))
    return rates


def controls(settings, output):
    """Check full-voltage empty/single/double controls, including peaks above 1 V.

    These use the known template with electronic noise, not new detector events.
    Truth is used to construct and score controls, never supplied to the fitter.
    """
    template, config = make_template(settings), FitConfig(**settings['fit_config'])
    rng = np.random.default_rng(2610141430)
    t = np.arange(3000)*2. + .25
    rows = []
    for name, times, charges in [
        ('empty', [], []), ('isolated', [1000.3], [4000.]),
        ('isolated_above_old_range', [1000.3], [40000.]),
        ('double_above_old_range', [1000.3, 1080.7], [40000., 20000.])]:
        times, charges = np.asarray(times), np.asarray(charges)
        amplitude = charges*config.area_per_pe_mV_ns*template.unit_area_peak
        voltage = config.baseline_mV + template(t[:, None]-times)@amplitude
        voltage += rng.normal(0., .1, len(t))
        fitted = fit_waveform(t, voltage, template, config)
        fitted.pop('model_mV')
        scored = score(fitted, dict(first_active_ns=times, emission_ns=times, charge_pe=charges), [500., 4000.])
        passed = (fitted['status'] == 'fitted' and scored['recovered_count'] == len(times)
                  and scored['reconstructed_count'] == len(times))
        if not passed:
            raise AssertionError(f'Full-voltage control failed: {name}')
        rows.append(dict(name=name, peak_mV=float(voltage.max()), passed=passed, fit=fitted, score=scored))
    save(Path(output)/'controls.json', dict(seed=2610141430, rows=rows,
        model='Frozen alpha template plus 0.1 mV Gaussian noise; no clipping or quantization',
        limitation='Software controls, not physical double-pulse resolution or detector calibration'))


def figures(output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output = Path(output)
    summary = json.loads((output/'summary.json').read_text())
    rates, rows = summary['rates'], summary['windows']
    labels = [f"{r['rate_hz']/1e6:g}" for r in rates]
    x = np.arange(len(rates))
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    def savefig(fig, name):
        for ext in ('png', 'pdf'):
            fig.savefig(output/f'{name}.{ext}', dpi=170)
        plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), layout='constrained')
    axes[0].plot(x, [100*r['recovery_fraction'] for r in rates], 'o-', color='#176994')
    for i, r in enumerate(rates):
        axes[0].annotate(f"{r['recovered_count']}/{r['target_count']}",
            (i, 100*r['recovery_fraction']), xytext=(0, 10), textcoords='offset points', ha='center', fontsize=8)
        if r['charge_relative_p16_p50_p84']:
            lo, mid, hi = np.array(r['charge_relative_p16_p50_p84'])*100
            axes[1].errorbar(i, mid, [[mid-lo], [hi-mid]], fmt='o', capsize=4, color='#a14c32')
    axes[0].set_ylim(-5, 116)
    axes[0].set_ylabel('Event recovery (%)')
    axes[0].set_title('Time error <= 10 ns and charge error <= 20%', fontsize=10)
    axes[1].axhline(0, ls=':', color='.5', lw=.8)
    axes[1].set_ylabel('Matched charge error: median, P16-P84 (%)')
    axes[1].set_title('Conditional on a time match; not energy resolution', fontsize=10)
    for ax in axes:
        ax.set_xticks(x, labels)
        ax.set_xlabel('Source alpha injection rate (MHz)')
        ax.grid(alpha=.18)
    fig.suptitle('Adequate dynamic range: all 36 windows use complete noisy voltage', fontsize=12)
    savefig(fig, 'full_signal_rates')

    selected = [r for r in rates if 2e7 <= r['rate_hz'] <= 8e7]
    fig, axes = plt.subplots(len(selected), 1, figsize=(11, 10), sharex=True, layout='constrained', squeeze=False)
    for ax, rate in zip(axes[:, 0], selected):
        row = next(r for r in rows if r['rate_hz'] == rate['rate_hz'])
        stem = Path(row['root']).stem
        with np.load(Path(row['source'])/'waveforms'/f'{stem}.npz') as w:
            t, voltage = w['time_ns']/1000, w['analog_mV']
        with np.load(output/'fits'/f'{stem}.npz') as f:
            model = f['model_mV']
        ax.plot(t, voltage, color='.3', lw=.9, label='Complete sampled voltage')
        ax.plot(t, model, color='#b05239', lw=.9, alpha=.85, label='Blind total fit')
        ax.text(.01, .97, f"{rate['rate_hz']/1e6:g} MHz | first predeclared window | "
            f"recovered {row['score']['recovered_count']}/{row['score']['target_count']}",
            transform=ax.transAxes, va='top', fontsize=9, bbox=dict(facecolor='white', alpha=.9, edgecolor='none'))
        ax.set_ylabel('Voltage (mV)'); ax.set_xlim(5, 15); ax.grid(alpha=.15)
        ax.margins(y=.2)
    axes[0, 0].legend(loc='lower right', fontsize=8)
    axes[-1, 0].set_xlabel('Record time (us)')
    savefig(fig, 'full_signal_waveforms')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, nargs='+', required=True, help='Existing analysis result folders')
    parser.add_argument('--output', type=Path, required=True, help='New output folder')
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        parser.error('workers must be between 1 and 8')
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/'fits').mkdir()
    tasks, input_records, roots = [], [], set()
    settings = json.loads((args.inputs[0]/'analysis_config.json').read_text())
    for source in args.inputs:
        other = json.loads((source/'analysis_config.json').read_text())
        for key in ('template', 'fit_config', 'calibration_sha256', 'scoring'):
            if settings[key] != other[key]:
                raise ValueError(f'Incompatible {key}: {source}')
        previous = json.loads((source/'summary.json').read_text())
        for row in previous['windows']:
            if row['root_sha256'] in roots:
                raise ValueError('Duplicated physical window')
            roots.add(row['root_sha256'])
            tasks.append((str(source), {k: row[k] for k in ('label', 'root', 'root_sha256', 'rate_hz')},
                          settings, str(args.output)))
        input_records.append(dict(folder=str(source), windows=len(previous['windows']),
            summary_sha256=digest(source/'summary.json'), config_sha256=digest(source/'analysis_config.json')))
    config = dict(settings, input_policy='complete_sampled_voltage', input_field='analog_mV',
        dynamic_range_assumption='Adequate range; no ADC clipping, rejection or quantization in this comparison',
        retained_response='Actual shared-cell SiPM history, linear electronics, sample aperture and electronic noise',
        historical_template_note='Calibration metadata retains original ADC settings; they do not limit input voltage',
        fit_truth_inputs=[], rate_input_to_fit=False, created=datetime.now().astimezone().isoformat(),
        environment=dict(python=platform.python_version(), platform=platform.platform(),
            **{k: importlib.metadata.version(k) for k in ('numpy', 'scipy', 'matplotlib', 'uproot')}))
    save(args.output/'analysis_config.json', config)
    controls(settings, args.output)
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(fit_window, tasks):
            rows.append(result)
            save(args.output/'progress.json', dict(completed=len(rows), planned=len(tasks)))
    summary = dict(state='completed', input_policy=config['input_policy'],
        windows=rows, rates=summarize(rows), inputs=input_records,
        caveats=['All declared independent windows included; no ADC-conditioned selection',
                 'Adequate-range simulation benchmark without quantization, not a particular 16-bit ADC validation',
                 'Known pure-alpha family; pulse count, times and charge unknown to the fitter',
                 'Realized output charge, not deposited-energy reconstruction',
                 'Exploratory small samples and a fixed pulse width; not a universal rate limit'])
    save(args.output/'summary.json', summary)
    figures(args.output)
    print('COMPLETE', args.output/'summary.json', flush=True)


if __name__ == '__main__':
    main()
