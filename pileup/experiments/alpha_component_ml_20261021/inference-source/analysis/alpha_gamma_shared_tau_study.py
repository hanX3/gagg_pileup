"""Calibrate shared lifetimes, choose constraints on development, freeze/test.

Run from the project root; all saved paths are project relative. Prior public
test waveforms are deliberately absent from this experiment's input records.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from functools import lru_cache
import json
from pathlib import Path
import shutil
import time

import numpy as np
from scipy.linalg import cholesky

from .alpha_gamma_shared_tau import (
    ProfileProblem, build_fraction_classifier, fit_joint_shared, fit_problem,
    fit_shared_prefix,
)
from .alpha_gamma_short_window import PrefixCalibration, find_trigger, fit_prefix
from .alpha_gamma_short_window_study import digest, save, wilson

PROJECT = Path(__file__).resolve().parents[1]
OLD = 'jupyter/artifacts/alpha_gamma_advanced_fit_20261020'
METHODS = ('fixed_tau', 'shared_tau', 'constrained_tau', 'free_tau', 'template_gls')
POSTS = (50., 75., 100., 150., 200., 300., 500., 750., 1000.)
INFLATIONS = (1., 3., 10.)


def calibration_problem(row, calibration):
    with np.load(PROJECT/row['waveform']) as wave:
        t, y = wave['time_ns'], wave['analog_mV']
        trigger = find_trigger(t, y, 20., calibration.config.threshold_mV)
        if trigger is None:
            return None, None
        base = fit_prefix(t, y, calibration, pre_ns=20., post_ns=1000., trigger=trigger,
                          methods=('physics_ols',))['methods']['physics_ols']
        mask = (t >= trigger['time_ns']-20.-1e-7) & (t <= trigger['time_ns']+1000.+1e-7)
        x, values = t[mask]-trigger['time_ns'], y[mask]
        C = calibration.covariance(x, base['onset_relative_ns'], max(1., base['charge_pe']))
        return ProfileProblem(x, values, calibration.readout, calibration.meta['area_per_pe_mV_ns'],
                              calibration.meta['baseline_mV'], .2,
                              cholesky(C, lower=True, check_finite=False)), trigger


@lru_cache(maxsize=2)
def models(project, output):
    project, output = Path(project), Path(output)
    calibration = PrefixCalibration.load(project/OLD/'calibration')
    shared = json.loads((output/'shared-fit.json').read_text())
    classifier = (json.loads((output/'classifier.json').read_text())
                  if (output/'classifier.json').exists() else None)
    return calibration, shared, classifier


def calibration_worker(task):
    project, output, row = task
    global PROJECT
    PROJECT = Path(project)
    cal, shared, _ = models(project, output)
    problem, _ = calibration_problem(row, cal)
    fits = {}
    for method in ('fixed_tau', 'shared_tau', 'free_tau'):
        fits[method] = fit_problem(problem, method, shared, initial=fits.get('shared_tau'))
    for value in INFLATIONS:
        fits[f'constrained_tau_{value:g}'] = fit_problem(
            problem, 'constrained_tau', shared, value, fits['free_tau'])
    result = {'root': row['root'], 'q': row['calibration_charge_pe'],
              'particle': row['particle'], 'fits': fits,
              'fitter_sha256': digest(Path(__file__).with_name('alpha_gamma_shared_tau.py')),
              'shared_sha256': digest(Path(output)/'shared-fit.json')}
    save(Path(output)/'calibration-fits'/(Path(row['root']).stem+'.json'), result)
    return row['root']


def calibrate(output, workers=4):
    output = Path(output).resolve()
    rows = json.loads((output/'input-records.json').read_text())['calibration']
    cal = PrefixCalibration.load(PROJECT/OLD/'calibration')
    retained = []
    for row in rows:
        for key, expected in row['cached_sha256'].items():
            assert digest(PROJECT/row[key]) == expected, row[key]
        with np.load(PROJECT/row['truth']) as truth:
            q = float(truth['charge_pe'][0])
            onset = float(truth['first_active_ns'][0])
        if q >= 250.:
            retained.append({**row, 'calibration_charge_pe': q, 'calibration_first_active_ns': onset})
    if not (output/'shared-fit.json').exists():
        problems, selected = [], []
        for row in retained:
            problem, trigger = calibration_problem(row, cal)
            if problem is None or abs(trigger['time_ns']-row['calibration_first_active_ns']) > 60.:
                continue
            problems.append(problem); selected.append(row)
        save(output/'shared-calibration-inputs.json', selected)
        shared = fit_joint_shared(problems)
        shared.update(fitter_sha256=digest(Path(__file__).with_name('alpha_gamma_shared_tau.py')),
                      covariance_source=OLD+'/calibration', scope='Effective double-decay/readout model, not a microscopic lifetime measurement')
        save(output/'shared-fit.json', shared)
        del problems
    selected = json.loads((output/'shared-calibration-inputs.json').read_text())
    current_hash = digest(Path(__file__).with_name('alpha_gamma_shared_tau.py'))
    assert json.loads((output/'shared-fit.json').read_text())['fitter_sha256'] == current_hash
    tasks = []
    for row in selected:
        p = output/'calibration-fits'/(Path(row['root']).stem+'.json')
        if p.exists():
            old = json.loads(p.read_text())
            assert old['fitter_sha256'] == current_hash and old['shared_sha256'] == digest(output/'shared-fit.json')
        else:
            tasks.append((str(PROJECT), str(output), row))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, result in enumerate(as_completed([pool.submit(calibration_worker, t) for t in tasks]), 1):
            print('CALIBRATION_FIT', i, '/', len(tasks), result.result(), flush=True)
    fits = [json.loads((output/'calibration-fits'/(Path(r['root']).stem+'.json')).read_text()) for r in selected]
    save(output/'classifier.json', build_fraction_classifier(fits))
    print('CALIBRATION_COMPLETE', len(fits), flush=True)


def event_worker(task):
    project, output, row, destination, inflation, candidates = task
    global PROJECT
    PROJECT = Path(project)
    output = Path(output)
    cal, shared, classifier = models(project, str(output))
    started = time.perf_counter()
    with np.load(PROJECT/row['waveform']) as wave:
        t, y = wave['time_ns'], wave['analog_mV']
        trigger = find_trigger(t, y, 20., cal.config.threshold_mV)
        prefixes = [fit_shared_prefix(t, y, cal, shared, post_ns=post, inflation=inflation,
                     classifier=classifier, trigger=trigger, candidate_inflations=candidates)
                    for post in POSTS] if trigger else [dict(detected=False, pre_ns=20., post_ns=post, methods={}) for post in POSTS]
    # Independent truth is opened only after every method/prefix is complete.
    with np.load(PROJECT/row['truth']) as truth:
        q = float(truth['charge_pe'][0]); onset = float(truth['first_active_ns'][0])
        particle = str(truth['particle_name'][0]); energy = float(truth['edep_MeV'][0])
    for prefix in prefixes:
        matched = bool(prefix['detected'] and np.isfinite(onset) and abs(prefix['trigger']['time_ns']-onset) <= 60.)
        prefix['evaluation'] = dict(matched_signal=matched, target_above_500pe=q >= 500., no_signal=q == 0.)
        for name, fit in prefix['methods'].items():
            error = fit['charge_pe']/q-1 if matched and q > 0 else None
            interval = fit.get('charge_interval95_pe')
            fit['evaluation'] = dict(relative_charge_error=error, correct_type=bool(matched and fit['particle'] == particle),
                onset_error_ns=prefix['trigger']['time_ns']+fit['onset_relative_ns']-onset if matched else None,
                interval_covers_true_Q=bool(matched and interval and interval[0] <= q <= interval[1]),
                tau_fast_relative_to_G4=fit['tau_fast_ns']/95.-1 if 'tau_fast_ns' in fit else None,
                tau_slow_relative_to_G4=fit['tau_slow_ns']/351.-1 if 'tau_slow_ns' in fit else None)
    result = dict(row=row, prefixes=prefixes, truth_for_evaluation_only=dict(charge_pe=q, particle=particle,
        first_active_ns=onset if np.isfinite(onset) else None, edep_MeV=energy),
        seconds=time.perf_counter()-started, inference_truth_inputs=[],
        fitter_sha256=digest(Path(__file__).with_name('alpha_gamma_shared_tau.py')),
        study_sha256=digest(__file__), shared_sha256=digest(output/'shared-fit.json'),
        classifier_sha256=digest(output/'classifier.json'))
    save(Path(destination)/(Path(row['root']).stem+'.json'), result)
    return row['root'], result['seconds']


def scan(rows, output, stage, workers=4):
    output = Path(output).resolve()
    destination = output/stage
    destination.mkdir(exist_ok=True)
    inflation = 1. if stage == 'development' else json.loads((output/'freeze.json').read_text())['prior_inflation']
    candidates = INFLATIONS if stage == 'development' else None
    tasks = []
    for row in rows:
        p = destination/(Path(row['root']).stem+'.json')
        if p.exists():
            previous = json.loads(p.read_text())
            for key, path in (('fitter_sha256', Path(__file__).with_name('alpha_gamma_shared_tau.py')),
                              ('study_sha256', Path(__file__)), ('shared_sha256', output/'shared-fit.json'),
                              ('classifier_sha256', output/'classifier.json')):
                assert previous[key] == digest(path), path
            assert previous['row'] == row
            continue
        tasks.append((str(PROJECT), str(output), row, str(destination), inflation, candidates))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, future in enumerate(as_completed([pool.submit(event_worker, t) for t in tasks]), 1):
            root, seconds = future.result()
            print(stage.upper(), i, '/', len(tasks), root, round(seconds, 2), 's', flush=True)


def load_results(folder):
    return [json.loads(p.read_text()) for p in sorted(Path(folder).glob('*.json'))]


def summary(results, methods):
    groups = []
    for post in POSTS:
        for method in methods:
            for particle in ('alpha', 'gamma'):
                for band, lo, hi in (('all_above_500', 500., np.inf), ('500_1000', 500., 1000.),
                                     ('1000_3000', 1000., 3000.), ('above_3000', 3000., np.inf)):
                    events = [r for r in results if r['truth_for_evaluation_only']['particle'] == particle
                              and lo <= r['truth_for_evaluation_only']['charge_pe'] < hi]
                    counts = dict(detected=0, correct_type=0, charge_5pct_success=0, joint_charge_type_success=0,
                                  tau_both_within_10pct_G4=0, tau_bound_count=0, onset_bound_count=0,
                                  optimizer_failure_count=0, interval_valid_count=0, interval_covered_count=0)
                    errors, fasts, slows, tau_errors, timing = [], [], [], [], []
                    for event in events:
                        prefix = next(p for p in event['prefixes'] if p['post_ns'] == post)
                        if not prefix['evaluation']['matched_signal']:
                            continue
                        fit = prefix['methods'][method]; ev = fit['evaluation']
                        error = ev['relative_charge_error']
                        errors.append(error); timing.append(ev['onset_error_ns'])
                        counts['detected'] += 1; counts['correct_type'] += ev['correct_type']
                        counts['charge_5pct_success'] += abs(error) <= .05
                        counts['joint_charge_type_success'] += abs(error) <= .05 and ev['correct_type']
                        counts['tau_bound_count'] += fit.get('tau_at_search_bound', False)
                        counts['onset_bound_count'] += fit.get('onset_at_search_bound', False)
                        counts['optimizer_failure_count'] += not fit.get('optimizer_converged', True)
                        counts['interval_valid_count'] += fit.get('charge_interval95_pe') is not None
                        counts['interval_covered_count'] += ev['interval_covers_true_Q']
                        if 'tau_fast_ns' in fit:
                            fasts.append(fit['tau_fast_ns']); slows.append(fit['tau_slow_ns'])
                            err = [ev['tau_fast_relative_to_G4'], ev['tau_slow_relative_to_G4']]
                            tau_errors.append(err)
                            counts['tau_both_within_10pct_G4'] += bool(np.max(abs(np.array(err))) <= .1)
                    n = len(events)
                    groups.append(dict(pre_ns=20., post_ns=post, method=method, particle=particle, band=band,
                        targets=n, **counts, missed=n-counts['detected'], wrong_type=counts['detected']-counts['correct_type'],
                        joint_fraction=counts['joint_charge_type_success']/n if n else None,
                        charge_5pct_fraction=counts['charge_5pct_success']/n if n else None,
                        joint_wilson95=wilson(counts['joint_charge_type_success'], n),
                        charge_bias=float(np.mean(errors)) if errors else None,
                        charge_error_sd=float(np.std(errors, ddof=1)) if len(errors) > 1 else None,
                        charge_rmse=float(np.sqrt(np.mean(np.square(errors)))) if errors else None,
                        onset_rmse_ns=float(np.sqrt(np.mean(np.square(timing)))) if timing else None,
                        tau_fast_median_ns=float(np.median(fasts)) if fasts else None,
                        tau_slow_median_ns=float(np.median(slows)) if slows else None,
                        tau_fast_rmse_vs_G4=float(np.sqrt(np.mean(np.array(tau_errors)[:, 0]**2))) if tau_errors else None,
                        tau_slow_rmse_vs_G4=float(np.sqrt(np.mean(np.array(tau_errors)[:, 1]**2))) if tau_errors else None,
                        interval_coverage=counts['interval_covered_count']/n if n else None))
    zero = [r for r in results if r['truth_for_evaluation_only']['charge_pe'] == 0.]
    false = []
    for post in POSTS:
        for method in methods:
            prefixes = [next(p for p in r['prefixes'] if p['post_ns'] == post) for r in zero]
            false.append(dict(post_ns=post, method=method, zero_signal_windows=len(zero),
                false_triggers=sum(p['detected'] for p in prefixes),
                false_Q_above_500=sum(p['detected'] and p['methods'][method]['charge_pe'] >= 500. for p in prefixes)))
    return dict(groups=groups, false_triggers=false, independent_windows=len(results),
                truth_counts={c: dict(injected=sum(r['truth_for_evaluation_only']['particle'] == c for r in results),
                    targets=sum(r['truth_for_evaluation_only']['particle'] == c and r['truth_for_evaluation_only']['charge_pe'] >= 500. for r in results))
                              for c in ('alpha', 'gamma')},
                scope='Isolated pulse; Q is avalanche charge, not energy; microscopic tau agreement is a model-mismatch diagnostic in native data')


def freeze(output, experiment):
    output, experiment = Path(output).resolve(), Path(experiment).resolve()
    assert not (output/'freeze.json').exists()
    rows = json.loads((output/'input-records.json').read_text())
    results = load_results(output/'development')
    assert len(results) == len(rows['development']) == 129
    all_methods = METHODS+tuple(f'constrained_tau_{v:g}' for v in INFLATIONS)
    s = summary(results, all_methods)
    save(output/'development-summary.json', s)
    ranked = []
    for value in INFLATIONS:
        groups = [g for g in s['groups'] if g['method'] == f'constrained_tau_{value:g}'
                  and g['band'] == 'all_above_500' and g['post_ns'] in (100., 150., 200., 300.)]
        ranked.append(dict(inflation=value, mean_joint_success=float(np.mean([g['joint_fraction'] for g in groups])),
                           mean_charge_rmse=float(np.mean([g['charge_rmse'] for g in groups]))))
    ranked.sort(key=lambda r: (-r['mean_joint_success'], r['mean_charge_rmse'], r['inflation']))
    files = list((PROJECT/'analysis').glob('*.py'))+list((PROJECT/'pileup/analysis').glob('*.py'))
    files += [PROJECT/'requirements-analysis.lock.txt']
    for p in files:
        target = experiment/'inference-source'/p.relative_to(PROJECT)
        target.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(p, target)
    for package in ('analysis', 'pileup', 'pileup/analysis'):
        p = experiment/'inference-source'/package/'__init__.py'
        if not p.exists():
            p.write_text('')
    selected = [output/n for n in ('protocol.json', 'shared-fit.json', 'classifier.json', 'input-records.json')]
    selected += [PROJECT/OLD/'calibration'/n for n in ('calibration.json', 'calibration-arrays.npz')]
    selected += list((experiment/'inference-source').rglob('*.py'))
    settings = dict(created_utc=datetime.now(timezone.utc).isoformat(), prior_inflation=ranked[0]['inflation'],
                    development_ranking=ranked, methods=METHODS, post_ns=POSTS,
                    primary_goal='Each class >=95% joint Q5/type point fraction; descriptive finite sample only',
                    frozen_files={str(p.relative_to(PROJECT)): digest(p) for p in selected})
    save(output/'freeze.json', settings); save(experiment/'frozen-fit-settings.json', settings)
    print('FROZEN_PRIOR', ranked, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=('calibrate', 'development', 'freeze', 'test', 'summarize'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--experiment', type=Path)
    p.add_argument('--workers', type=int, default=4)
    args = p.parse_args(); output = args.output.resolve()
    if args.stage == 'calibrate':
        calibrate(output, args.workers)
    elif args.stage == 'freeze':
        freeze(output, args.experiment)
    elif args.stage == 'development':
        scan(json.loads((output/'input-records.json').read_text())['development'], output, 'development', args.workers)
    elif args.stage == 'test':
        settings = json.loads((output/'freeze.json').read_text())
        for relative, expected in settings['frozen_files'].items():
            assert digest(PROJECT/relative) == expected, relative
        scan(json.loads((output/'native-records.json').read_text())['test_single'], output, 'test', args.workers)
    else:
        results = load_results(output/'test')
        assert len(results) == 168
        save(output/'test-summary.json', summary(results, METHODS))


if __name__ == '__main__':
    main()
