"""Blind mixed-template fit with separate modelling and reporting thresholds.

This pilot retains fitted weak pulses as nuisance contributions to the waveform.
It does not supply their truth labels, times, charges or counts to reconstruction.
The original baseline and its immutable campaign snapshots remain unchanged.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime
import importlib.metadata
import json
from pathlib import Path
import platform
import time

import numpy as np

from .alpha_gamma_pileup_fit import fit_waveform, score_mixed
from .alpha_gamma_study import calibration
from .alpha_rate_campaign import digest, save


def fit_with_weak_components(time_ns, voltage_mV, bank, config, *, internal_threshold_pe=200.,
                             report_threshold_pe=500.):
    if not 0 < internal_threshold_pe <= report_threshold_pe:
        raise ValueError('Require 0 < internal threshold <= reporting threshold')
    internal_config = replace(config, min_charge_pe=float(internal_threshold_pe))
    full = fit_waveform(time_ns, voltage_mV, bank, internal_config)
    # Keep the full model for the residual: the reporting cut never deletes its
    # weak contributions or re-fits the strong amplitudes without those terms.
    full_model = full.pop('model_mV')
    keep = np.asarray(full['charge_pe']) >= report_threshold_pe
    reported = dict(full)
    for key in ('times_ns', 'particle_name', 'charge_pe', 'area_mV_ns'):
        if key in full:
            reported[key] = np.asarray(full[key])[keep].tolist()
    reported['nuisance_component_count'] = int(np.sum(~keep))
    reported['model_includes_nuisance_components'] = True
    return reported, full, full_model


def run(settings, inputs, original_result, output, *, internal_threshold_pe=200.):
    output, inputs = Path(output), Path(inputs)
    destination = output / 'result.json'
    if destination.exists():
        raise FileExistsError('Existing pilot result; choose a new output path')
    original = json.loads(Path(original_result).read_text())
    row = original['row']
    stem = Path(row['root']).stem
    _, bank, config = calibration(settings)
    with np.load(inputs/'waveforms'/f'{stem}.npz') as a:
        waveform = dict(a)
    start = time.perf_counter()
    reported, full, model = fit_with_weak_components(
        waveform['time_ns'], waveform['analog_mV'], bank, config,
        internal_threshold_pe=internal_threshold_pe, report_threshold_pe=500.)
    elapsed = time.perf_counter()-start
    # No truth file is opened until the reconstruction above has finished.
    with np.load(inputs/'truth'/f'{stem}.npz') as a:
        truth = dict(a)
    meta = json.loads((inputs/'metadata'/f'{stem}.json').read_text())
    scored = score_mixed(reported, truth, meta['gate_ns'])
    source_files = [Path(__file__), Path(__file__).with_name('alpha_gamma_pileup_fit.py'),
                    Path(__file__).with_name('alpha_pileup_fit.py')]
    result = dict(created=datetime.now().astimezone().isoformat(), row=row,
                  role='Blind weak-component pilot; preflight tuning if input is a preflight window',
                  fit_truth_inputs=[], internal_threshold_pe=internal_threshold_pe,
                  reporting_threshold_pe=500., evaluation_threshold_pe=500.,
                  original_result=str(original_result), original_result_sha256=digest(original_result),
                  settings_sha256=digest(settings), sources={str(p): digest(p) for p in source_files},
                  waveform_sha256=digest(inputs/'waveforms'/f'{stem}.npz'),
                  truth_sha256=digest(inputs/'truth'/f'{stem}.npz'), fit_seconds=elapsed,
                  environment=dict(python=platform.python_version(),
                                   numpy=importlib.metadata.version('numpy'),
                                   scipy=importlib.metadata.version('scipy')),
                  fit=reported, all_components=full, score=scored,
                  original_score=original['score'], original_rms_mV=original['fit']['rms_mV'])
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'reconstruction.npz', time_ns=waveform['time_ns'],
                        model_including_weak_mV=model,
                        times_ns=full['times_ns'], charge_pe=full['charge_pe'],
                        particle_name=full['particle_name'])
    save(destination, result)
    print('WEAKFIT', row['label'], internal_threshold_pe,
          {p: (v['recovered_count'], v['target_count']) for p, v in scored['by_particle'].items()},
          'extras', scored['unmatched_predictions'], 'seconds', round(elapsed, 2), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('settings', 'inputs', 'original-result', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--internal-threshold-pe', type=float, default=200.)
    args = parser.parse_args()
    run(args.settings, args.inputs, args.original_result, args.output,
        internal_threshold_pe=args.internal_threshold_pe)
