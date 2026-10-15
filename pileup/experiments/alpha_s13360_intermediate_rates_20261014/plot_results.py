"""Static figures comparing the new intermediate rates with the frozen pilot."""
from pathlib import Path
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT))
from pileup.analysis.alpha_rate_campaign import digest, save

out = PROJECT / 'jupyter/artifacts/alpha_intermediate_rates_20261014_114000'
old = PROJECT / 'jupyter/artifacts/alpha_rate_fit_20261013_211541/summary.json'
width_file = PROJECT / 'jupyter/artifacts/pulse_duration_rate_20261014_105829/summary.json'
new = json.loads((out / 'summary.json').read_text())
assert new['state'] == 'completed'
prior = json.loads(old.read_text())
widths = json.loads(width_file.read_text())['widths']
rates = sorted(prior['rates'] + new['rates'], key=lambda x:x['rate_hz'])
rows = prior['windows'] + new['windows']
plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False,
                     'figure.dpi':140, 'savefig.dpi':180})

def write_figure(fig, name):
    for suffix in ('png', 'pdf'):
        fig.savefig(out / f'{name}.{suffix}')
    plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout='constrained')
x = np.array([r['rate_hz'] / 1e6 for r in rates])
axes[0].plot(x, [100*r['recovery_fraction'] for r in rates], 'o-', color='#176994', label='Pooled event recovery')
for r in rows:
    sc = r['adc_score']
    if sc['target_count']:
        axes[0].scatter(r['rate_hz']/1e6, 100*sc['recovered_count']/sc['target_count'],
                        s=23, facecolors='none', edgecolors='#176994', alpha=.45)
    axes[1].scatter(r['rate_hz']/1e6, 100*r['central_clipped_samples']/r['central_samples'],
                    s=23, facecolors='none', edgecolors='#ad4c30', alpha=.45)
axes[1].plot(x, [100*r['clipped_fraction'] for r in rates], 's-', color='#ad4c30', label='Pooled central rail samples')
for ax in axes:
    ax.set_xscale('log'); ax.set_xticks(x, labels=[f'{v:g}' for v in x])
    plt.setp(ax.get_xticklabels(), rotation=55, ha='right', fontsize=9)
    ax.set_ylim(-4, 110); ax.set_xlabel('Source alpha injection rate (MHz)'); ax.grid(alpha=.18)
    ax.legend(fontsize=8)
axes[0].set_ylabel('ADC recovery, including rejected windows (%)')
axes[1].set_ylabel('ADC rail fraction in [5, 15) us (%)')
axes[0].set_title('Timing within 10 ns and charge error within 20%', fontsize=10)
axes[1].set_title('1000 mV ADC ceiling; any rail rejects whole record', fontsize=10)
fig.suptitle('Fixed GAGG alpha template and readout; open markers = individual windows', fontsize=11)
write_figure(fig, 'rate_comparison')

fig, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True, layout='constrained')
for ax, rate in zip(axes, new['rates']):
    row = next(r for r in new['windows'] if r['rate_hz'] == rate['rate_hz'])
    stem = Path(row['root']).stem
    with np.load(out/'waveforms'/f'{stem}.npz') as wave:
        ax.plot(wave['time_ns']/1000, wave['analog_mV'], color='#4d9679', lw=1, label='Pre-ADC simulation voltage')
        ax.plot(wave['time_ns']/1000, wave['adc_mV'], color='#333333', lw=.7, label='ADC')
        if row['adc_fit']['status'] == 'fitted':
            with np.load(out/'fits'/f'{stem}.npz') as fit:
                ax.plot(wave['time_ns']/1000, fit['model_mV'], color='#176994', lw=1, label='Blind ADC fit')
    ax.axhline(1000, color='#ad4c30', ls=':', lw=.9)
    ax.text(.01, .97, f"{rate['rate_hz']/1e6:g} MHz | first predeclared window | {row['adc_fit']['status']}\n"
            f"ADC recovered {row['adc_score']['recovered_count']}/{row['adc_score']['target_count']}",
            transform=ax.transAxes, va='top', fontsize=9,
            bbox=dict(facecolor='white', edgecolor='none', alpha=.85))
    ax.set_ylabel('Voltage (mV)'); ax.grid(alpha=.15); ax.set_xlim(5,15)
axes[0].legend(loc='upper right', fontsize=8, ncol=3)
axes[-1].set_xlabel('Record time (us)')
write_figure(fig, 'intermediate_waveforms')

derived = []
for r in new['rates']:
    duration_s = r['gate_time_ns'] * 1e-9
    signal_rate = r['central_nonzero_charge'] / duration_s
    derived.append(dict(source_rate_hz=r['rate_hz'], window_count=r['windows'], exposure_s=duration_s,
        realized_injection_rate_hz=r['central_injected']/duration_s,
        charge_producing_rate_hz=signal_rate, target_rate_hz=r['target_count']/duration_s,
        nominal_injection_mean_interval_ns=1e9/r['rate_hz'],
        charge_rate_reciprocal_ns=1e9/signal_rate,
        occupancy_fwhm=signal_rate*widths['fwhm_ns']*1e-9,
        occupancy_w99=signal_rate*widths['area_99_duration_ns']*1e-9))
save(out/'derived_rate_metrics.json', dict(rows=derived, template_widths=widths,
    meaning='Observed charge-producing rate times fixed isolated-pulse width; descriptive occupancy, not measured resolving time or failure probability',
    inputs_sha256={str(old.relative_to(PROJECT)):digest(old),str(width_file.relative_to(PROJECT)):digest(width_file),
                  str((out/'summary.json').relative_to(PROJECT)):digest(out/'summary.json')}))
print('Saved two figures (PNG/PDF) and duration/rate metrics', flush=True)
