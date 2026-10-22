"""Post-test voltage-residual and fixed-region diagnostics; no model selection."""
from pathlib import Path
import json
import numpy as np
from .alpha_full_signal_study import make_template
from .alpha_rate_campaign import save
from .alpha_ml_multimethod_report import LABELS


def diagnostics(output, experiment, benchmark='new_test'):
    output, experiment = Path(output), Path(experiment)
    folder = output / benchmark
    summary = json.loads((folder / 'summary.json').read_text())
    settings = json.loads((experiment / 'calibration-settings.json').read_text())
    template = make_template(settings)
    baseline = settings['fit_config']['baseline_mV']
    charge_scale = settings['fit_config']['area_per_pe_mV_ns'] * template.unit_area_peak
    rows = []
    for window in summary['windows']:
        with np.load(window['waveform_path']) as wave:
            t, y = wave['time_ns'], wave['analog_mV']
        gate = (t >= window['gate_ns'][0]) & (t < window['gate_ns'][1])
        for method in summary['method_order']:
            fitted = window['methods'][method]['fit']
            times, charges = np.asarray(fitted['times_ns']), np.asarray(fitted['charge_pe'])
            predicted = baseline + template(t[:, None] - times) @ (charges * charge_scale)
            central_predictions = (times >= window['gate_ns'][0]) & (times < window['gate_ns'][1])
            rows.append(dict(label=window['label'], rate_hz=window['rate_hz'], method=method,
                template_voltage_rms_mV=float(np.sqrt(np.mean((y - predicted) ** 2))),
                central_template_voltage_rms_mV=float(np.sqrt(np.mean((y[gate] - predicted[gate]) ** 2))),
                maximum_guard_prediction_charge_pe=float(charges[~central_predictions].max(initial=0.)),
                maximum_central_prediction_charge_pe=float(charges[central_predictions].max(initial=0.)),
                recovered_count=window['methods'][method]['score']['recovered_count'],
                is_neural_template_projection=method in ('cnn', 'unet', 'transformer')))
    rates = []
    for rate in sorted({r['rate_hz'] for r in rows}):
        for method in summary['method_order']:
            group = [r for r in rows if r['rate_hz'] == rate and r['method'] == method]
            rates.append(dict(rate_hz=rate, method=method,
                median_template_voltage_rms_mV=float(np.median([r['template_voltage_rms_mV'] for r in group])),
                median_central_template_voltage_rms_mV=float(np.median([r['central_template_voltage_rms_mV'] for r in group])),
                maximum_guard_prediction_charge_pe=max(r['maximum_guard_prediction_charge_pe'] for r in group),
                maximum_central_prediction_charge_pe=max(r['maximum_central_prediction_charge_pe'] for r in group)))
    save(folder / 'voltage-diagnostics.json', dict(windows=rows, rates=rates,
        used_for_model_selection=False,
        interpretation='Neural outputs projected through the same fixed mean template for diagnosis only; '
                       'a small voltage residual is not an event-identity or charge-recovery guarantee. '
                       'Report both full-record and central RMS because predictions outside the trained '
                       'central region can dominate full-record errors.'))
    _examples(folder, summary, settings, template)
    print('POST-TEST VOLTAGE DIAGNOSTICS', benchmark, len(rows), flush=True)


def _examples(folder, summary, settings, template):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    # Rates, first window, and display region fixed before looking for successes.
    selected = [next(w for w in summary['windows'] if w['rate_hz'] == rate)
                for rate in (20000000, 80000000)]
    lo, hi = 6000., 6600.
    baseline = settings['fit_config']['baseline_mV']
    scale = settings['fit_config']['area_per_pe_mV_ns'] * template.unit_area_peak
    fig, axes = plt.subplots(2, 2, figsize=(13, 7.5), constrained_layout=True)
    colors = dict(fit='tab:red', cnn='tab:blue', unet='tab:green',
                  transformer='tab:orange', cnn_physical='tab:purple')
    for column, window in enumerate(selected):
        with np.load(window['waveform_path']) as wave:
            t, y = wave['time_ns'], wave['analog_mV']
        region = (t >= lo) & (t <= hi)
        axes[0, column].plot(t[region] / 1000., y[region], color='black', lw=1., label='Measured voltage')
        for method in ('fit', 'cnn_physical'):
            fitted = window['methods'][method]['fit']
            model = baseline + template(t[region, None] - np.asarray(fitted['times_ns'])) @ (
                np.asarray(fitted['charge_pe']) * scale)
            axes[0, column].plot(t[region] / 1000., model, color=colors[method], label=LABELS[method])
        # Truth is used only in this post-test comparison figure.
        with np.load(window['truth_path']) as truth:
            tt, qt = truth['first_active_ns'], truth['charge_pe']
        target = np.isfinite(tt) & (qt >= 500.) & (tt >= lo) & (tt <= hi)
        axes[1, column].scatter(tt[target] / 1000., np.zeros(target.sum()), marker='|', s=110,
                                color='0.35', label='Truth q >= 500 pe')
        methods = summary['method_order']
        for level, method in enumerate(methods, start=1):
            tp = np.asarray(window['methods'][method]['fit']['times_ns'])
            tp = tp[(tp >= lo) & (tp <= hi)]
            axes[1, column].scatter(tp / 1000., np.full(len(tp), level), marker='|', s=110,
                                    color=colors[method])
        axes[1, column].set_yticks(range(6), ['Truth'] + [LABELS[k] for k in methods])
        axes[1, column].invert_yaxis()
        axes[1, column].set_title('Event times (charge scored separately)')
        axes[0, column].set_title(f"{window['rate_hz']/1e6:g} MHz source rate, {window['label']}")
        axes[0, column].set_ylabel('Voltage (mV)')
        axes[0, column].legend(fontsize=8)
        for row in range(2):
            axes[row, column].set_xlim(lo / 1000., hi / 1000.)
            axes[row, column].set_xlabel('Time (us)')
            axes[row, column].grid(alpha=.2)
    fig.suptitle(f"{summary['benchmark']}: fixed 6.0-6.6 us regions; sufficient acquisition range")
    for suffix in ('png', 'pdf'):
        fig.savefig(folder / f'multimethod_examples.{suffix}', dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--experiment', type=Path, required=True)
    parser.add_argument('--benchmark', choices=('regression', 'new_test'), default='new_test')
    arguments = parser.parse_args()
    diagnostics(arguments.output, arguments.experiment, arguments.benchmark)
