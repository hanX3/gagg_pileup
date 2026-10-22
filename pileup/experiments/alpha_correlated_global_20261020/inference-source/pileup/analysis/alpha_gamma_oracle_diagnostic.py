"""Truth-assisted onset/type diagnostics; never a blind reconstruction result.

Run on the two preflight windows, outside the 32-window formal benchmark.
The observed analog waveform and original templates are unchanged. Knowing
all nonzero-charge arrival times includes the weak population as nuisance
components. Compare with knowing only the above-threshold arrival times.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import importlib.metadata
import json
from pathlib import Path
import platform
import time

import numpy as np
from scipy.optimize import nnls

from .alpha_gamma_pileup_fit import score_mixed
from .alpha_gamma_study import calibration
from .alpha_rate_campaign import digest, save


def oracle_fit(t, voltage, truth, bank, config, *, include_weak, known_types):
    """Deliberately consumes truth: this function must not enter a blind fit."""
    selected = np.flatnonzero(
        (truth['charge_pe'] > 0 if include_weak else truth['charge_pe'] >= 500.)
        & np.isfinite(truth['first_active_ns'])
    )
    times = truth['first_active_ns'][selected]
    order = np.argsort(times)
    selected, times = selected[order], times[order]
    values = voltage - config.baseline_mV
    if known_types:
        types = (truth['particle_name'][selected] == 'gamma').astype(int)
        matrix = bank.matrix(t, times, types, config)
        charge = nnls(matrix, values, maxiter=max(1000, 30 * len(times)))[0]
        model = matrix @ charge
        assignment = 'True particle types; only amplitudes estimated'
    else:
        matrix = bank.matrix(t, np.repeat(times, 2), np.tile([0, 1], len(times)), config)
        dual = nnls(matrix, values, maxiter=max(1000, 60 * len(times)))[0].reshape(-1, 2)
        model = matrix @ dual.ravel()
        charge = dual.sum(axis=1)
        types = np.argmax(dual, axis=1)
        assignment = ('Relaxed alpha+gamma coefficients at each true onset; sum gives charge, '
                      'larger coefficient gives label. Not an exact discrete type optimum.')
    keep = charge >= 500.
    fit = dict(times_ns=times[keep].tolist(), particle_name=bank.names[types[keep]].tolist(),
               charge_pe=charge[keep].tolist(), status='truth_assisted')
    # The model includes every nuisance component, even those below the output
    # threshold. Thus its residual is not a reconstruction from reported peaks.
    result = dict(fit=fit, candidate_count=len(times), include_weak=include_weak,
                  known_types=known_types, assignment=assignment,
                  all_component_rms_mV=float(np.sqrt(np.mean((values-model)**2))),
                  fitted_truth_indices=selected.tolist(), all_charge_pe=charge.tolist(),
                  all_fitted_types=bank.names[types].tolist(),
                  all_component_model_includes_below_output_threshold=True)
    assert np.isfinite(model).all() and np.all(charge >= 0)
    # Numerical optimum check for the continuous NNLS problem (not the labels).
    coefficients = charge if known_types else dual.ravel()
    gradient = matrix.T @ (model-values)
    scale = max(1., float(np.max(np.abs(matrix.T @ values))))
    violation = max(float(np.max(np.maximum(-gradient, 0))),
                    float(np.max(np.abs(gradient[coefficients > 1e-6])))
                    if np.any(coefficients > 1e-6) else 0.) / scale
    assert violation < 1e-7, violation
    result['nnls_relative_kkt_violation'] = violation
    return result


def run(settings, inputs, originals, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError('Choose a new diagnostic output; preserve old runs')
    inputs = Path(inputs)
    _, bank, config = calibration(settings)
    control_t = np.arange(.25, 6000., 2.)
    control_truth = dict(charge_pe=np.array([2000., 800., 120.]),
                         first_active_ns=np.array([1000., 1030., 1700.]),
                         particle_name=np.array(['alpha', 'gamma', 'gamma']))
    control_y = (bank.matrix(control_t, control_truth['first_active_ns'], [0, 1, 1], config)
                 @ control_truth['charge_pe'] + config.baseline_mV)
    control_checks = []
    for known in (True, False):
        control = oracle_fit(control_t, control_y, control_truth, bank, config,
                             include_weak=True, known_types=known)
        error = float(np.max(np.abs(np.array(control['all_charge_pe'])-control_truth['charge_pe'])))
        assert error < 1e-5 and control['all_component_rms_mV'] < 1e-7
        assert control['all_fitted_types'] == control_truth['particle_name'].tolist()
        control_checks.append(dict(known_types=known, max_charge_error_pe=error,
                                   rms_mV=control['all_component_rms_mV']))
    rows = []
    for original_path in originals:
        original = json.loads(Path(original_path).read_text())
        row = original['row']
        # Preflight role must be explicit in the campaign row or its label.
        assert 'preflight' in row['label'] or 'check' in row['label'], row['label']
        stem = Path(row['root']).stem
        with np.load(inputs/'waveforms'/f'{stem}.npz') as a:
            waveform = dict(a)
        with np.load(inputs/'truth'/f'{stem}.npz') as a:
            truth = dict(a)
        meta = json.loads((inputs/'metadata'/f'{stem}.json').read_text())
        cases = []
        for include_weak, known_types in ((False, True), (True, True), (True, False)):
            start = time.perf_counter()
            case = oracle_fit(waveform['time_ns'], waveform['analog_mV'], truth,
                              bank, config, include_weak=include_weak, known_types=known_types)
            case['score'] = score_mixed(case['fit'], truth, meta['gate_ns'])
            case['elapsed_s'] = time.perf_counter()-start
            cases.append(case)
            print(row['label'], 'weak', include_weak, 'types', known_types,
                  {p: (v['recovered_count'], v['target_count'])
                   for p, v in case['score']['by_particle'].items()},
                  'rms', round(case['all_component_rms_mV'], 3), flush=True)
        assert cases[1]['all_component_rms_mV'] <= cases[0]['all_component_rms_mV'] + 1e-8
        assert cases[2]['all_component_rms_mV'] <= cases[1]['all_component_rms_mV'] + 1e-8
        rows.append(dict(row=row, original_result=str(original_path),
                         original_sha256=digest(original_path),
                         waveform_sha256=digest(inputs/'waveforms'/f'{stem}.npz'),
                         truth_sha256=digest(inputs/'truth'/f'{stem}.npz'),
                         original_score=original['score'], original_fit_rms_mV=original['fit']['rms_mV'],
                         cases=cases))
    source_files = [Path(__file__), Path(__file__).with_name('alpha_gamma_pileup_fit.py'),
                    Path(__file__).with_name('alpha_pileup_fit.py')]
    result = dict(created=datetime.now().astimezone().isoformat(),
                  purpose='Truth-assisted failure diagnosis on preflight windows; not blind performance',
                  settings=str(settings), settings_sha256=digest(settings),
                  sources={str(p): digest(p) for p in source_files},
                  environment=dict(python=platform.python_version(),
                                   numpy=importlib.metadata.version('numpy'),
                                   scipy=importlib.metadata.version('scipy')),
                  truth_inputs=['nonzero/above-threshold event count', 'first_active_ns',
                                'particle_name when known_types=True'],
                  unchanged='Observed analog waveform, electronics noise, microcell history and original templates',
                  output_threshold_pe=500., evaluation_threshold_pe=500.,
                  exact_template_double_plus_weak_controls=control_checks, rows=rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    save(output, result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--originals', type=Path, nargs='+', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.settings, args.inputs, args.originals, args.output)
