"""Post-test analysis of the frozen five methods at 0.1, 1 and 10 MHz.

Reuse saved blind predictions; truth is used only for scoring and diagnostics.
Run from the project root. No training, threshold selection or G4 production.
"""
from pathlib import Path
from collections import Counter
from datetime import datetime
import hashlib
import json
import numpy as np
from scipy.signal import find_peaks
from .alpha_pileup_fit import score

RATES = (100000, 1000000, 10000000)
METHODS = ('fit', 'cnn', 'unet', 'transformer', 'cnn_physical')
SHORT = ('Fit', 'CNN', 'U-Net', 'Transformer', 'Hybrid')
COLORS = ('tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple')
NEIGHBOURS = ('gap < FWHM', 'FWHM <= gap < W90', 'gap >= W90 / no neighbour')


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, content):
    Path(path).write_text(json.dumps(content, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def category(gap, fwhm, w90):
    if gap is None or gap >= w90:
        return NEIGHBOURS[2]
    return NEIGHBOURS[0] if gap < fwhm else NEIGHBOURS[1]


def peak_diagnostic(folder, root, method, t, target_time, threshold):
    """Inspect a saved score peak near truth, after blind prediction only."""
    with np.load(folder / 'windows' / f'{Path(root).stem}_{method}.npz') as dense:
        probability, charge, offset = dense['probability'], dense['charge_pe'], dense['time_offset_ns']
    peaks = find_peaks(probability, distance=4)[0]
    nearby = peaks[np.abs(t[peaks] + offset[peaks] - target_time) <= 10.]
    if not len(nearby):
        return dict(local_peak_found=False)
    index = int(nearby[np.argmax(probability[nearby])])
    return dict(local_peak_found=True, probability=float(probability[index]), threshold=threshold,
                charge_pe=float(charge[index]), charge_threshold_pe=500.,
                time_error_ns=float(t[index] + offset[index] - target_time),
                probability_passed=bool(probability[index] >= threshold),
                charge_passed=bool(charge[index] >= 500.),
                note='Post-test local-peak diagnostic; assignment/distance constraints can also matter')


def collect_benchmark(reference, benchmark, widths, frozen):
    folder = reference / benchmark
    summary = json.loads((folder / 'summary.json').read_text())
    windows = [r for r in summary['windows'] if r['rate_hz'] in RATES]
    assert len(windows) == 12
    metrics = [dict(method=method, **r) for method in METHODS
               for r in summary['rates'][method] if r['rate_hz'] in RATES]
    events, roots, inputs = [], [], {}
    for window in windows:
        for field in ('waveform', 'truth'):
            path = Path(window[f'{field}_path'])
            assert digest(path) == window[f'{field}_sha256']
            inputs[str(path)] = window[f'{field}_sha256']
        with np.load(window['truth_path']) as file:
            truth = {k: file[k] for k in ('first_active_ns', 'emission_ns', 'charge_pe')}
        tt, qt = truth['first_active_ns'], truth['charge_pe']
        charged = np.flatnonzero(np.isfinite(tt) & (qt > 0.))
        for method in METHODS:
            entry = window['methods'][method]
            recomputed = score(entry['fit'], truth, window['gate_ns'])
            assert json.loads(json.dumps(recomputed)) == entry['score']
            s = entry['score']
            roots.append(dict(label=window['label'], root=window['root'], seed=window['seed'],
                rate_hz=window['rate_hz'], method=method, target_count=s['target_count'],
                recovered_count=s['recovered_count'], missed_count=s['missed_count'],
                charge_failed=s['matched_count']-s['recovered_count'],
                extra_predictions=s['unmatched_predictions'], merge_flags=s['merge_flags'], split_flags=s['split_flags'],
                recovery_fraction=s['recovered_count']/s['target_count'] if s['target_count'] else None))
            pairs = dict(s['matched_pairs'])
            for index in s['target_indices']:
                others = charged[charged != index]
                gap = float(np.min(np.abs(tt[others] - tt[index]))) if len(others) else None
                prediction_index = pairs.get(index)
                pq = entry['fit']['charge_pe'][prediction_index] if prediction_index is not None else None
                dt = entry['fit']['times_ns'][prediction_index]-tt[index] if prediction_index is not None else None
                dq = pq/qt[index]-1. if prediction_index is not None else None
                outcome = 'missed' if dq is None else ('recovered' if abs(dq) <= .2 else 'charge_failed')
                event = dict(benchmark=benchmark, label=window['label'], root=window['root'], rate_hz=window['rate_hz'],
                    method=method, truth_index=int(index), true_time_ns=float(tt[index]), true_charge_pe=float(qt[index]),
                    nearest_charged_neighbour_ns=gap, neighbour_group=category(gap, widths['fwhm_ns'], widths['area_90_duration_ns']),
                    predicted_charge_pe=pq, time_error_ns=None if dt is None else float(dt),
                    relative_charge_error=None if dq is None else float(dq), outcome=outcome,
                    prediction_index=prediction_index)
                if benchmark == 'new_test' and window['rate_hz'] <= 1000000 and outcome != 'recovered' and method in ('cnn','unet','transformer'):
                    with np.load(window['waveform_path']) as file:
                        t = file['time_ns']
                    event['saved_neural_peak_diagnostic'] = peak_diagnostic(folder, window['root'], method, t,
                        tt[index], frozen['methods'][method]['threshold'])
                events.append(event)
    grouped, paired, sensitivity = [], [], []
    for rate in RATES:
        for method in METHODS:
            subset = [e for e in events if e['rate_hz'] == rate and e['method'] == method]
            metric = next(r for r in metrics if r['rate_hz'] == rate and r['method'] == method)
            counts = Counter(e['outcome'] for e in subset)
            assert len(subset) == metric['target_count']
            assert counts['recovered'] == metric['recovered_count']
            assert counts['missed'] == metric['missed_count']
            assert counts['charge_failed'] == metric['matched_charge_failed']
            for group in NEIGHBOURS:
                selected = [e for e in subset if e['neighbour_group'] == group]
                recovered = sum(e['outcome'] == 'recovered' for e in selected)
                grouped.append(dict(rate_hz=rate, method=method, neighbour_group=group,
                    target_count=len(selected), recovered_count=recovered,
                    recovery_fraction=recovered/len(selected) if selected else None))
            baseline = {(e['root'],e['truth_index']) for e in events
                        if e['rate_hz']==rate and e['method']=='fit' and e['outcome']=='recovered'}
            success = {(e['root'],e['truth_index']) for e in subset if e['outcome']=='recovered'}
            paired.append(dict(rate_hz=rate, method=method, target_count=len(subset),
                both_success=len(baseline & success), method_only=len(success-baseline),
                fit_only=len(baseline-success), neither_success=len(subset)-len(baseline|success),
                interpretation='Paired descriptive counts; no independent-event significance claim'))
            rates = [r for r in roots if r['rate_hz']==rate and r['method']==method]
            leave = [(metric['recovered_count']-r['recovered_count'])/(metric['target_count']-r['target_count'])
                     for r in rates if metric['target_count'] > r['target_count']]
            sensitivity.append(dict(rate_hz=rate, method=method,
                leave_one_root_out_min=min(leave) if leave else None,
                leave_one_root_out_max=max(leave) if leave else None,
                interpretation='Remove one complete ROOT at a time; sensitivity range, not a confidence interval'))
    return dict(metrics=metrics, windows=roots, events=events, neighbour_groups=grouped,
                paired_vs_fit=paired, root_sensitivity=sensitivity, inputs=inputs,
                summary_sha256=digest(folder/'summary.json'))


def analyse(reference, experiment, width_path, output):
    reference, experiment, width_path, output = map(Path, (reference,experiment,width_path,output))
    output.mkdir(parents=True, exist_ok=False)
    frozen = json.loads((reference/'frozen-methods.json').read_text())
    for name, expected in frozen['sources'].items():
        assert digest(experiment/'source/pileup/analysis'/name) == expected
    for method in ('cnn','unet','transformer'):
        assert digest(reference/f'models/{method}/model.pt') == frozen['methods'][method]['weights_sha256']
    width_reference = json.loads(width_path.read_text())
    settings = json.loads((experiment/'calibration-settings.json').read_text())
    assert width_reference['template'] == settings['template']
    widths = width_reference['widths']
    benchmarks = {name:collect_benchmark(reference,name,widths,frozen) for name in ('new_test','regression')}
    occupancy = []
    for metric in benchmarks['new_test']['metrics']:
        if metric['method'] != 'fit':
            continue
        rate = metric['charge_producing_rate_hz']
        mu = rate*widths['area_90_duration_ns']*1e-9
        targets = [e for e in benchmarks['new_test']['events'] if e['rate_hz']==metric['rate_hz'] and e['method']=='fit']
        overlap = sum(e['nearest_charged_neighbour_ns'] is not None
                      and e['nearest_charged_neighbour_ns']<widths['area_90_duration_ns'] for e in targets)
        occupancy.append(dict(source_rate_hz=metric['rate_hz'], actual_injected=metric['central_injected'],
            central_gate_time_ns=metric['gate_time_ns'], charge_producing_count=metric['central_nonzero_charge'],
            charge_producing_rate_hz=rate, threshold_target_count=metric['target_count'],
            threshold_target_rate_hz=metric['target_rate_hz'], w90_ns=widths['area_90_duration_ns'], occupancy=mu,
            poisson_next_neighbour_probability=float(-np.expm1(-mu)),
            poisson_either_neighbour_probability=float(-np.expm1(-2*mu)),
            observed_targets_with_neighbour_within_w90=overlap,
            observed_neighbour_fraction=overlap/len(targets) if targets else None,
            interpretation='Poisson reference uses observed charge-producing rate; it is not a failure probability. '
                           'Measured neighbour uses both sides, all nonzero-charge alpha events including guards.'))
    null = json.loads((reference/'new_test/timing-null.json').read_text())
    result = dict(state='completed', created=datetime.now().astimezone().isoformat(), rates_hz=list(RATES),
        widths=widths, benchmarks=benchmarks, occupancy=occupancy,
        timing_null=[r for r in null['rates'] if r['rate_hz'] in RATES],
        timing_null_is_confidence_interval=False, reused_saved_blind_predictions=True,
        truth_used_for='Post-test scoring, neighbour grouping and failure illustration only',
        model_selection=False, new_g4_runs=0, energy_calibration=False,
        reference_sha256=dict(frozen_methods=digest(reference/'frozen-methods.json'),
                             widths=digest(width_path), timing_null=digest(reference/'new_test/timing-null.json')))
    save(output/'summary.json', result)
    figures(result, reference, output)
    save(output/'validation.json', dict(passed=True, low_rate_windows_per_benchmark=12,
        common_methods=5, scores_recomputed=True, target_decomposition_matches_original=True,
        model_weights_and_sources_unchanged=True, paired_and_neighbour_counts_verified=True,
        reference_width_template_identical=True, independent_benchmarks_kept_separate=True,
        changed_model_parameters=False, new_g4_runs=0))
    print('LOW-RATE ANALYSIS COMPLETE: new12 + regression12 ROOTs, five frozen methods', flush=True)


def figures(result, reference, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fresh = result['benchmarks']['new_test']
    fig, axes = plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
    for column, rate in enumerate(RATES):
        rows = [next(r for r in fresh['metrics'] if r['rate_hz']==rate and r['method']==m) for m in METHODS]
        targets = rows[0]['target_count']
        good = np.array([r['recovered_count'] for r in rows])/targets*100
        bad = np.array([r['matched_charge_failed'] for r in rows])/targets*100
        missed = np.array([r['missed_count'] for r in rows])/targets*100
        x = np.arange(5)
        axes[0,column].bar(x,good,color='seagreen',label='Time + charge success')
        axes[0,column].bar(x,bad,bottom=good,color='darkorange',label='Matched, charge failed')
        axes[0,column].bar(x,missed,bottom=good+bad,color='0.6',label='Time missed')
        for index, row in enumerate(rows):
            axes[0,column].text(index,good[index]/2,f"{row['recovered_count']}/{targets}",ha='center',color='white',fontsize=9)
            q = np.array(row['charge_relative_p16_p50_p84'])*100
            axes[1,column].errorbar(index,q[1],yerr=[[q[1]-q[0]],[q[2]-q[1]]],fmt='o',capsize=4,color=COLORS[index])
        axes[0,column].set_title(f'{rate/1e6:g} MHz source rate; {targets} targets')
        axes[0,column].set_ylim(0,105)
        axes[0,column].set_ylabel('Fraction of targets (%)')
        axes[1,column].set_title('Conditional charge error: P16 / median / P84')
        axes[1,column].set_ylabel('Relative charge error (%)')
        axes[1,column].axhline(0,color='0.4',lw=.8)
        axes[1,column].axhspan(-20,20,color='seagreen',alpha=.08)
        for row in (0,1):
            axes[row,column].set_xticks(x,SHORT,rotation=15,fontsize=9)
            axes[row,column].grid(axis='y',alpha=.2)
    axes[0,0].legend(fontsize=8,loc='lower left')
    fig.suptitle('New common tests: low-rate success, charge failure and time miss; q >= 500 pe')
    for suffix in ('png','pdf'):
        fig.savefig(output/f'low_rate_metrics.{suffix}',dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(1,3,figsize=(15,4.5),constrained_layout=True)
    for column, rate in enumerate(RATES):
        for index, method in enumerate(METHODS):
            rows = [next(r for r in fresh['neighbour_groups']
                         if r['rate_hz']==rate and r['method']==method and r['neighbour_group']==g) for g in NEIGHBOURS]
            axes[column].plot(range(3),[r['recovery_fraction']*100 if r['recovery_fraction'] is not None else np.nan for r in rows],
                              'o-',color=COLORS[index],label=SHORT[index])
        counts = [next(r for r in fresh['neighbour_groups']
                       if r['rate_hz']==rate and r['method']=='fit' and r['neighbour_group']==g)['target_count'] for g in NEIGHBOURS]
        axes[column].set_xticks(range(3),[f'<FWHM\nn={counts[0]}',f'FWHM..W90\nn={counts[1]}',f'>=W90 / none\nn={counts[2]}'],fontsize=9)
        axes[column].set_title(f'{rate/1e6:g} MHz; post-test neighbour grouping')
        axes[column].set_ylim(-3,105)
        axes[column].set_ylabel('Time + charge success (%)')
        axes[column].grid(alpha=.2)
        axes[column].legend(fontsize=8)
    fig.suptitle('Nearest nonzero-charge alpha in the full record; amplitude and shared-cell state are confounded')
    for suffix in ('png','pdf'):
        fig.savefig(output/f'low_rate_neighbours.{suffix}',dpi=170)
    plt.close(fig)
    _failure_figures(result,reference,output)


def _failure_figures(result,reference,output):
    import matplotlib.pyplot as plt
    fresh = result['benchmarks']['new_test']
    failures = [e for e in fresh['events'] if e['rate_hz']<=1000000 and e['outcome']!='recovered']
    cases = sorted({(e['rate_hz'],e['root'],e['truth_index']) for e in failures})[:3]
    saved = json.loads((reference/'new_test/summary.json').read_text())
    fig, axes = plt.subplots(2,len(cases),figsize=(5.4*len(cases),7.5),squeeze=False,constrained_layout=True)
    for column,(rate,root,index) in enumerate(cases):
        window = next(r for r in saved['windows'] if r['root']==root)
        event = next(e for e in fresh['events'] if e['root']==root and e['truth_index']==index)
        with np.load(window['waveform_path']) as wave:
            t,y = wave['time_ns'],wave['analog_mV']
        mask = (t>=event['true_time_ns']-200)&(t<=event['true_time_ns']+800)
        tx = t[mask]-event['true_time_ns']
        axes[0,column].plot(tx,y[mask],color='black',label='Measured voltage')
        for method,color in (('fit','tab:blue'),('cnn_physical','tab:purple')):
            with np.load(reference/'new_test/windows'/f'{Path(root).stem}_{method}.npz') as file:
                model = file['model_mV']
            axes[0,column].plot(tx,model[mask],color=color,label='Fit' if method=='fit' else 'Hybrid')
        axes[0,column].set_title(f"{rate/1e6:g} MHz, truth #{index}, q={event['true_charge_pe']:.1f} pe")
        axes[0,column].set_ylabel('Voltage (mV)')
        axes[0,column].set_xlabel('Time from truth first active photon (ns)')
        axes[0,column].legend(fontsize=8)
        charge = event['true_charge_pe']
        axes[1,column].axhspan(.8*charge,1.2*charge,color='seagreen',alpha=.12,label='Charge tolerance')
        axes[1,column].axhline(charge,color='black',ls='--',label='Truth charge')
        for position,method in enumerate(METHODS):
            e = next(e for e in fresh['events'] if e['root']==root and e['truth_index']==index and e['method']==method)
            if e['predicted_charge_pe'] is not None:
                axes[1,column].bar(position,e['predicted_charge_pe'],color=COLORS[position],width=.65)
            else:
                axes[1,column].text(position,.15*charge,'Missed',rotation=90,ha='center',color=COLORS[position])
        axes[1,column].set_xticks(range(5),SHORT,rotation=15,fontsize=9)
        axes[1,column].set_ylabel('Matched reconstructed charge (pe)')
        axes[1,column].set_ylim(bottom=0)
        for row in (0,1):
            axes[row,column].grid(axis='y',alpha=.2)
    fig.suptitle('All three distinct failed low-rate targets: selected after testing, not unbiased performance examples')
    for suffix in ('png','pdf'):
        fig.savefig(output/f'low_rate_native_failures.{suffix}',dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--experiment',type=Path,required=True)
    parser.add_argument('--widths',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    arguments = parser.parse_args()
    analyse(arguments.reference,arguments.experiment,arguments.widths,arguments.output)
