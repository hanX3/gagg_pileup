"""Frozen-method controls and post-inference diagnostics for the ML refresh."""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime
import json
from pathlib import Path

import numpy as np

from . import alpha_ml_refresh as study
from .alpha_full_signal_study import make_template
from .alpha_pileup_fit import FitConfig, score
from .alpha_rate_campaign import digest, save


def controls():
    study.check_frozen()
    settings, _, covariance = study.configuration()
    template, base = make_template(settings), FitConfig(**settings['fit_config'])
    nets = study.networks()
    folder = study.OUTPUT / 'controls'
    for name in ('waveforms', 'truth', 'results'):
        (folder / name).mkdir(exist_ok=True)
    t = np.arange(0., 6000., 2.)
    plans = [('empty', [], [], 0., 0.)] * 2
    plans += [('single', [1000.3], [q], 0., 0.) for q in (500., 1000., 4000., 16000., 40000.)]
    plans += [('double', [1000.3, 1000.3 + gap], [4000., 4000. * ratio], gap, ratio)
              for gap in (4., 8., 10., 12., 20., 40., 80., 200.) for ratio in (.125, 1., 4.)]
    plans += [('double', [1000.3, 1010.3], [4000., 2000.], 10., .5),
              ('edge_context', [1000.3, 5790.3], [1000., 4000.], 0., 0.)]
    rng = np.random.default_rng(262030800)
    rows = []
    for noise in ('white', 'correlated'):
        for i, (kind, onsets, charges, gap, ratio) in enumerate(plans):
            identifier = f'software_{noise}_{i:03}'
            output = folder / 'results' / (identifier + '.json')
            signal = template(t[:, None] - np.asarray(onsets)[None, :]) @ (
                np.asarray(charges) * base.area_per_pe_mV_ns * template.unit_area_peak)
            residual = rng.normal(0., covariance.electronic_sigma_mV, len(t))
            if noise == 'correlated':
                latent = np.empty(len(t))
                latent[0] = rng.normal()
                draws = rng.normal(size=len(t) - 1)
                for j in range(1, len(t)):
                    latent[j] = covariance.rho * latent[j - 1] + np.sqrt(1. - covariance.rho ** 2) * draws[j - 1]
                residual += np.sqrt(covariance.dark_variance_mV2 + covariance.shot_coefficient_mV * np.maximum(signal, 0.)) * latent
            if output.exists():
                rows.append(study.load(output))
                continue
            waveform, truth_path = folder / 'waveforms' / (identifier + '.npz'), folder / 'truth' / (identifier + '.npz')
            np.savez_compressed(waveform, time_ns=t, analog_mV=base.baseline_mV + signal + residual)
            np.savez_compressed(truth_path, first_active_ns=onsets, emission_ns=onsets, charge_pe=charges)
            record = dict(root='software/' + identifier, waveform=str(waveform.relative_to(study.PROJECT)),
                          truth=str(truth_path.relative_to(study.PROJECT)), gate_ns=[500., 4000.])
            fits, timing, _, _ = study.inference(record, nets, 'controls')
            truth = dict(first_active_ns=np.asarray(onsets), emission_ns=np.asarray(onsets), charge_pe=np.asarray(charges))
            methods = {name: dict(fit={k: v for k, v in fit.items() if k != 'model_mV'},
                                  score=score(fit, truth, record['gate_ns']), seconds=timing[name])
                       for name, fit in fits.items()}
            result = dict(identifier=identifier, noise=noise, kind=kind, gap_ns=gap, ratio=ratio,
                          onsets_ns=onsets, charges_pe=charges, waveform=record['waveform'],
                          methods=methods, known_count_input=False, known_time_input=False)
            save(output, result)
            rows.append(result)
            print('CONTROL', identifier, kind, flush=True)
    save(folder / 'summary.json', dict(seed=262030800, windows=len(rows), rows=rows,
         limitation='Matched mean template plus independently sampled synthetic noise; no native photon/shared-cell statistics, not detector resolution'))


def recovered_indices(data):
    return {int(i) for (i, _), dq in zip(data['score']['matched_pairs'], data['score']['relative_charge_errors'])
            if abs(dq) <= .2}


def report():
    study.check_frozen()
    result = dict(created=datetime.now().astimezone().isoformat(), benchmarks={})
    for phase in ('validation', 'common', 'confirmation'):
        s = study.load(study.OUTPUT / phase / 'summary.json')
        totals, rates, paired, windows, events = [], [], [], [], []
        for name in s['method_order']:
            rows = s['rates'][name]
            values = {k: sum(r[k] for r in rows) for k in
                      ('target_count', 'recovered_count', 'matched_count', 'missed_count', 'matched_charge_failed',
                       'reconstructed_count', 'unmatched_predictions', 'merge_flags', 'split_flags')}
            dq = [v for w in s['windows'] for v in w['methods'][name]['score']['relative_charge_errors']]
            dt = [v for w in s['windows'] for v in w['methods'][name]['score']['time_errors_ns']]
            timing = [w['methods'][name]['seconds'] for w in s['windows']]
            totals.append(dict(method=name, **values,
                recovery_fraction=values['recovered_count'] / max(values['target_count'], 1),
                charge_relative_p16_p50_p84=np.percentile(dq, [16, 50, 84]).tolist() if dq else None,
                time_error_ns_p16_p50_p84=np.percentile(dt, [16, 50, 84]).tolist() if dt else None,
                median_seconds=float(np.median(timing)), total_seconds=float(sum(timing))))
            rates.extend(dict(method=name, **r) for r in rows)
        for w in s['windows']:
            saved = study.load(study.OUTPUT / phase / 'windows' / (Path(w['root']).stem + '.json'))
            assert saved['waveform_sha256'] == digest(study.PROJECT / w['waveform'])
            with np.load(study.PROJECT / w['truth']) as z:
                truth = {k: z[k].copy() for k in ('first_active_ns', 'emission_ns', 'charge_pe')}
            baseline = recovered_indices(w['methods']['diagonal'])
            tt, qt = truth['first_active_ns'], truth['charge_pe']
            all_charged = np.flatnonzero(np.isfinite(tt) & (qt > 0.))
            for name, data in w['methods'].items():
                rebuilt = saved['methods'][name]['score']
                assert data['score'] == rebuilt
                s0 = data['score']
                assert s0['target_count'] == s0['recovered_count'] + s0['missed_count'] + s0['matched_count'] - s0['recovered_count']
                successful = recovered_indices(data)
                paired.append(dict(root=w['root'], rate_hz=w['rate_hz'], method=name,
                                   both=len(successful & baseline), method_only=len(successful - baseline),
                                   wls_only=len(baseline - successful)))
                pairs = {i: j for i, j in s0['matched_pairs']}
                for i in s0['target_indices']:
                    neighbours = all_charged[all_charged != i]
                    gap = float(np.min(abs(tt[neighbours] - tt[i]))) if len(neighbours) else None
                    j = pairs.get(i)
                    dq = None if j is None else float(data['fit']['charge_pe'][j] / qt[i] - 1.)
                    events.append(dict(root=w['root'], rate_hz=w['rate_hz'], method=name,
                        truth_index=i, true_charge_pe=float(qt[i]), first_active_ns=float(tt[i]),
                        nearest_charged_neighbour_ns=gap, relative_charge_error=dq,
                        outcome='missed' if j is None else 'recovered' if abs(dq) <= .2 else 'charge_failed'))
                windows.append(dict(root=w['root'], seed=w['seed'], rate_hz=w['rate_hz'], method=name,
                    target_count=s0['target_count'], recovered_count=s0['recovered_count'],
                    missed_count=s0['missed_count'], unmatched_predictions=s0['unmatched_predictions'],
                    matched_charge_failed=s0['matched_count'] - s0['recovered_count'],
                    recovery_fraction=s0['recovered_count'] / s0['target_count'] if s0['target_count'] else None,
                    neural_proposals=data['fit'].get('neural_proposals')))
        null = []
        if (study.OUTPUT / phase / 'timing-null.json').exists():
            groups = defaultdict(lambda: [0, 0])
            for r in study.load(study.OUTPUT / phase / 'timing-null.json')['rows']:
                a = groups[r['method'], r['rate_hz']]
                a[0] += r['score']['recovered_count']
                a[1] += r['score']['target_count']
            null = [dict(method=m, rate_hz=r, mean_recovery_fraction=v[0] / v[1] if v[1] else None)
                    for (m, r), v in groups.items()]
        result['benchmarks'][phase] = dict(windows=len(s['windows']), method_order=s['method_order'],
             totals=totals, rates=rates, window_metrics=windows, paired_vs_wls=paired, events=events, timing_null=null)
    result.update(energy_calibration=False, inference_truth_inputs=[], test_tuning=False,
                  cohort_interpretation='Public common and new confirmation reported separately; validation is not test performance')
    save(study.OUTPUT / 'report-tables.json', result)
    print('REPORT_READY', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['controls', 'report'])
    args = parser.parse_args()
    globals()[args.stage]()
