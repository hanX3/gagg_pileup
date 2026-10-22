"""Independent calibration/development/frozen-test comparison of prefix fits."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from functools import lru_cache
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import shutil
import time

import numpy as np

from .alpha_gamma_advanced_fit import (
    METHODS, AdvancedCalibration, AdvancedConfig, build_advanced_calibration, fit_advanced_prefix,
)
from .alpha_gamma_short_window import PrefixCalibration, find_trigger
from .alpha_gamma_short_window_study import digest, save, wilson

PROJECT = Path(__file__).resolve().parents[1]
OLD = 'jupyter/artifacts/alpha_gamma_short_window_implementation_20261019'
POSTS = (100., 150., 200., 300., 500.)
PRE = 20.


def checked_records(records, readout=None):
    if readout is None:
        readout = json.loads((PROJECT/records[0]['metadata']).read_text())['config']['readout']
    for row in records:
        meta = json.loads((PROJECT/row['metadata']).read_text())
        assert meta['config']['readout'] == readout
        assert meta['root_sha256'] == row['root_original_sha256']
        row['cached_sha256'] = {key: digest(PROJECT/row[key]) for key in ('waveform', 'truth', 'metadata')}
    return records, readout


def calibrate(output):
    """Calibration can finish while independent development is still simulated."""
    output = Path(output).resolve()
    old = json.loads((PROJECT/OLD/'input-records.json').read_text())
    fresh = json.loads((output/'native-records.json').read_text())
    assert len(fresh['calibration_single']) == 168
    records, readout = checked_records(old['calibration']+fresh['calibration_single'])
    paths = ('analysis/alpha_gamma_advanced_fit.py', 'analysis/alpha_gamma_short_window.py')
    provenance = {'records': records, 'readout': readout,
                  'builder_sha256': {p: digest(PROJECT/p) for p in paths}}
    if (output/'calibration').exists():
        assert json.loads((output/'calibration-inputs.json').read_text()) == provenance
        return AdvancedCalibration.load(output/'calibration')
    save(output/'calibration-inputs.json', provenance)
    started = time.perf_counter()
    cal = build_advanced_calibration(records, readout, output/'calibration',
                                     config=AdvancedConfig(), project=PROJECT)
    save(output/'calibration-summary.json', {
        'training_rows': cal.meta['training_rows'], 'used_training_rows': cal.meta['used_training_rows'],
        'classes': cal.meta['classes'], 'bootstrap_repetitions': 32,
        'old_used_training_rows': 70, 'old_gamma_bin_counts': [6, 6, 6],
        'seconds': time.perf_counter()-started,
        'covariance_interpolation': 'Interpolate correlated PSD component; interpolate independent variance onto observed diagonal',
    })
    print('ADVANCED_CALIBRATION', cal.meta['used_training_rows'], cal.meta['classes'], flush=True)
    return cal


def prepare(output):
    output = Path(output).resolve()
    assert (output/'protocol.json').is_file()
    old = json.loads((PROJECT/OLD/'input-records.json').read_text())
    fresh = json.loads((output/'native-records.json').read_text())
    assert len(fresh['calibration_single']) == 168 and len(fresh['development_single']) == 84
    groups = {'calibration': old['calibration']+fresh['calibration_single'],
              'development': old['development']+fresh['development_single']}
    roots = [r['root'] for group in groups.values() for r in group]
    seeds = [r['seed'] for group in groups.values() for r in group]
    assert len(roots) == len(set(roots)) and len(seeds) == len(set(seeds))
    _, readout = checked_records(groups['calibration'])
    checked_records(groups['development'], readout)
    save(output/'input-records.json', groups)
    calibrate(output)


def ready_development(output):
    """Use completed/validated native rows while other simulations finish."""
    output = Path(output).resolve()
    plan = json.loads((output/'production-plan.json').read_text())
    folder = PROJECT/plan['native_experiment']
    manifest = json.loads((folder/'manifest.json').read_text())
    rows = json.loads((PROJECT/OLD/'input-records.json').read_text())['development']
    for row in manifest['runs']:
        if row['role'] != 'development_single' or row['status'] != 'completed' or not row.get('validation'):
            continue
        stem = Path(row['root']).stem
        rows.append({'root': row['root'], 'seed': row['seed'], 'particle': row['particle'],
            'source_manifest': str((folder/'manifest.json').relative_to(PROJECT)),
            'root_original_sha256': row['root_sha256'],
            **{name: str((output/'native-inputs'/sub/(stem+suffix)).relative_to(PROJECT))
               for name, sub, suffix in (('waveform','waveforms','.npz'),('truth','truth','.npz'),
                                        ('metadata','metadata','.json'))}})
    checked_records(rows)
    save(output/'development-progress-inputs.json', rows)
    return rows


@lru_cache(maxsize=2)
def _worker_calibrations(project, calibration_folder):
    """Load immutable calibration once per process, without altering inference."""
    folder = Path(calibration_folder)
    return (AdvancedCalibration.load(folder), PrefixCalibration.load(Path(project)/OLD/'calibration'),
            {name: digest(folder/name) for name in (
                'calibration.json', 'calibration-arrays.npz', 'advanced-arrays.npz', 'advanced-calibration.json')},
            digest(Path(__file__).with_name('alpha_gamma_advanced_fit.py')), digest(Path(__file__)))


def _worker(task):
    project, row, calibration_folder, destination, methods, posts = task
    project = Path(project)
    cal, legacy, calibration_hashes, fitter_hash, study_hash = _worker_calibrations(str(project), calibration_folder)
    started = time.perf_counter()
    with np.load(project/row['waveform']) as wave:
        t, y = wave['time_ns'], wave['analog_mV']
        trigger = find_trigger(t, y, PRE, cal.config.threshold_mV)
        prefixes = ([fit_advanced_prefix(t, y, cal, legacy, pre_ns=PRE, post_ns=post,
                                         methods=methods, trigger=trigger) for post in posts]
                    if trigger is not None else
                    [{'detected': False, 'pre_ns': PRE, 'post_ns': post, 'methods': {}} for post in posts])
    # Only after ALL requested inference windows/methods are complete:
    with np.load(project/row['truth']) as truth:
        q = float(truth['charge_pe'][0]); name = str(truth['particle_name'][0])
        onset = float(truth['first_active_ns'][0])
        value = {'charge_pe': q, 'particle': name,
                 'first_active_ns': onset if np.isfinite(onset) else None,
                 'edep_MeV': float(truth['edep_MeV'][0])}
    for prefix in prefixes:
        matched = bool(prefix['detected'] and np.isfinite(onset)
                       and abs(prefix['trigger']['time_ns']-onset) <= 60.)
        prefix['evaluation'] = {'matched_signal': matched, 'target_above_500pe': q >= 500.,
                                'no_signal': q == 0.}
        for method, fitted in prefix['methods'].items():
            interval = fitted['charge_interval95_pe']
            fitted['evaluation'] = {
                'relative_charge_error': fitted['charge_pe']/q-1 if matched and q > 0 else None,
                'correct_type': bool(fitted['particle'] == name and matched),
                'onset_error_ns': prefix['trigger']['time_ns']+fitted['onset_relative_ns']-onset if matched else None,
                'interval_covers_true_Q': bool(interval is not None and matched and interval[0] <= q <= interval[1]),
            }
    result = {'row': row, 'truth_for_evaluation_only': value, 'prefixes': prefixes,
              'seconds': time.perf_counter()-started, 'inference_truth_inputs': [],
              'fitter_sha256': fitter_hash, 'study_sha256': study_hash,
              'calibration_sha256': calibration_hashes}
    save(Path(destination)/(Path(row['root']).stem+'.json'), result)
    return row['root'], result['seconds']


def scan(rows, output, destination, workers=4, methods=METHODS, posts=POSTS):
    output, destination = Path(output).resolve(), Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    missing = []
    expected = digest(Path(__file__).with_name('alpha_gamma_advanced_fit.py'))
    calibration_hashes = {name: digest(output/'calibration'/name) for name in (
        'calibration.json', 'calibration-arrays.npz', 'advanced-arrays.npz', 'advanced-calibration.json')}
    for row in rows:
        for key, value in row.get('cached_sha256', {}).items():
            assert digest(PROJECT/row[key]) == value, row[key]
        path = destination/(Path(row['root']).stem+'.json')
        if path.exists():
            previous = json.loads(path.read_text())
            assert previous['fitter_sha256'] == expected
            assert previous['study_sha256'] == digest(Path(__file__))
            assert previous['calibration_sha256'] == calibration_hashes
            assert previous['row'] == row
            assert [p['post_ns'] for p in previous['prefixes']] == list(posts)
            assert all(not prefix['detected'] or set(prefix['methods']) == set(methods)
                       for prefix in previous['prefixes'])
        else:
            missing.append(row)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_worker, (str(PROJECT), row, str(output/'calibration'),
                                         str(destination), methods, posts)) for row in missing]
        for i, future in enumerate(as_completed(futures), 1):
            root, seconds = future.result()
            print('FITTED', i, '/', len(missing), root, round(seconds, 2), 's', flush=True)


def load_results(folder):
    return [json.loads(path.read_text()) for path in sorted(Path(folder).glob('*.json'))]


def summary(results, methods=METHODS, posts=POSTS):
    groups = []
    for post in posts:
        for method in methods:
            for particle in ('alpha', 'gamma'):
                events = [r for r in results if r['truth_for_evaluation_only']['particle'] == particle]
                for band, lower, upper in (('all_above_500', 500., np.inf), ('500_1000', 500., 1000.),
                                           ('1000_3000', 1000., 3000.), ('above_3000', 3000., np.inf)):
                    selected = [r for r in events if lower <= r['truth_for_evaluation_only']['charge_pe'] < upper]
                    errors, time_errors = [], []
                    counts = dict(detected=0, correct_type=0, charge_5pct_success=0, joint_charge_type_success=0,
                                  onset_bound_count=0, charge_bound_count=0, optimizer_failure_count=0,
                                  interval_valid_count=0, interval_covered_count=0)
                    for event in selected:
                        prefix = next(p for p in event['prefixes'] if p['post_ns'] == post)
                        if not prefix['evaluation']['matched_signal']:
                            continue
                        fit = prefix['methods'][method]; ev = fit['evaluation']
                        error = ev['relative_charge_error']
                        errors.append(error); time_errors.append(ev['onset_error_ns'])
                        counts['detected'] += 1
                        counts['correct_type'] += ev['correct_type']
                        counts['charge_5pct_success'] += abs(error) <= .05
                        counts['joint_charge_type_success'] += abs(error) <= .05 and ev['correct_type']
                        counts['onset_bound_count'] += fit['onset_at_search_bound']
                        counts['charge_bound_count'] += fit.get('charge_at_search_bound', False)
                        counts['optimizer_failure_count'] += fit['optimizer_converged'] is False
                        counts['interval_valid_count'] += fit['charge_interval95_pe'] is not None
                        counts['interval_covered_count'] += ev['interval_covers_true_Q']
                    n = len(selected)
                    groups.append({'pre_ns': PRE, 'post_ns': post, 'method': method,
                                   'particle': particle, 'band': band, 'targets': n, **counts,
                                   'charge_5pct_fraction': counts['charge_5pct_success']/n if n else None,
                                   'joint_fraction': counts['joint_charge_type_success']/n if n else None,
                                   'joint_wilson95': wilson(counts['joint_charge_type_success'], n),
                                   'charge_bias': float(np.mean(errors)) if errors else None,
                                   'charge_error_sd': float(np.std(errors, ddof=1)) if len(errors)>1 else None,
                                   'charge_rmse': float(np.sqrt(np.mean(np.square(errors)))) if errors else None,
                                   'abs_charge_error_p90': float(np.quantile(np.abs(errors), .9)) if errors else None,
                                   'time_rmse_ns': float(np.sqrt(np.mean(np.square(time_errors)))) if time_errors else None,
                                   'interval_coverage_all_targets': counts['interval_covered_count']/n if n else None,
                                   'interval_coverage_valid': counts['interval_covered_count']/counts['interval_valid_count']
                                    if counts['interval_valid_count'] else None})
    false = []
    for post in posts:
        for method in methods:
            zero = [r for r in results if r['truth_for_evaluation_only']['charge_pe'] == 0.]
            triggers = 0; above = 0
            for event in zero:
                prefix = next(p for p in event['prefixes'] if p['post_ns'] == post)
                triggers += prefix['detected']
                above += bool(prefix['detected'] and prefix['methods'][method]['charge_pe'] >= 500.)
            false.append({'method': method, 'post_ns': post, 'zero_signal_windows': len(zero),
                          'false_triggers': triggers, 'false_Q_above_500': above})
    return {'groups': groups, 'false_triggers': false, 'independent_events': len(results),
            'truth_counts': {particle: {
                'injected': sum(r['truth_for_evaluation_only']['particle'] == particle for r in results),
                'above_500': sum(r['truth_for_evaluation_only']['particle'] == particle
                                  and r['truth_for_evaluation_only']['charge_pe'] >= 500. for r in results),
                'zero_charge': sum(r['truth_for_evaluation_only']['particle'] == particle
                                    and r['truth_for_evaluation_only']['charge_pe'] == 0. for r in results),
            } for particle in ('alpha', 'gamma')},
            'limits': ['At most one pulse assumed; no pulse-count search or blind pileup result',
                       'Charge pe is not deposited/source energy or energy resolution',
                       'Bias/SD/quantiles are conditional on matching; successes/coverage retain all targets',
                       'Intervals are local/model-conditional; coverage is observed, not a calibration guarantee']}


def group(s, method, particle, post):
    return next(g for g in s['groups'] if g['method'] == method and g['particle'] == particle
                and g['post_ns'] == post and g['band'] == 'all_above_500')


def freeze(output, experiment):
    output, experiment = Path(output).resolve(), Path(experiment).resolve()
    assert not (output/'freeze.json').exists()
    dev = load_results(output/'development')
    groups = json.loads((output/'input-records.json').read_text())
    assert len(dev) == len(groups['development'])
    s = summary(dev); save(output/'development-summary.json', s)
    ranking = []
    for method in METHODS[1:]:
        selected = [group(s, method, particle, post) for particle in ('alpha', 'gamma')
                    for post in (100., 150., 200., 300.)]
        ranking.append({'method': method,
                        'class_and_window_equal_joint_score': float(np.mean([g['joint_fraction'] for g in selected])),
                        'mean_charge_rmse': float(np.mean([g['charge_rmse'] for g in selected]))})
    ranking.sort(key=lambda r: (-r['class_and_window_equal_joint_score'], r['mean_charge_rmse'], METHODS.index(r['method'])))
    primary = ranking[0]['method']
    passing = [post for i, post in enumerate(POSTS) if all(
        group(s, primary, particle, longer)['joint_fraction'] >= .95
        for particle in ('alpha', 'gamma') for longer in POSTS[i:])]
    selected_window = passing[0] if passing else None
    snapshot = experiment/'inference-source'
    snapshot.mkdir(exist_ok=False)
    for area in ('analysis', 'pileup/analysis'):
        for source in (PROJECT/area).glob('*.py'):
            destination = snapshot/source.relative_to(PROJECT)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
    (snapshot/'analysis/__init__.py').write_text('')
    (snapshot/'pileup/__init__.py').write_text('')
    shutil.copy2(PROJECT/'requirements-analysis.lock.txt', snapshot/'requirements-analysis.lock.txt')
    code = ['analysis/alpha_gamma_advanced_fit.py', 'analysis/alpha_gamma_advanced_study.py',
            'analysis/alpha_gamma_advanced_native.py', 'analysis/test_alpha_gamma_advanced_fit.py',
            'analysis/alpha_gamma_advanced_native_local.py',
            'analysis/alpha_gamma_short_window.py', 'analysis/alpha_gamma_short_window_study.py',
            'analysis/alpha_gamma_waveform_fit.py']
    paths = [PROJECT/p for p in code]
    paths += [output/name for name in ('protocol.json', 'input-records.json', 'calibration-inputs.json', 'development-summary.json')]
    paths += [p for p in (output/'calibration').iterdir() if p.is_file()]
    paths += [PROJECT/OLD/'calibration'/name for name in ('calibration.json', 'calibration-arrays.npz')]
    for p in code:
        assert digest(PROJECT/p) == digest(snapshot/p)
    frozen = {
        'created': datetime.now().astimezone().isoformat(), 'selected_on_development_only': True,
        'primary_method': primary, 'primary_short_goal_ns': 150.,
        'selected_95pct_point_window_ns': selected_window,
        'selection_ranking': ranking, 'pre_ns': PRE, 'posts_ns': POSTS,
        'calibration_records': len(groups['calibration']), 'development_records': len(dev),
        'calibration_sha256': {p.name: digest(p) for p in (output/'calibration').iterdir() if p.is_file()},
        'new_native_test_generated': 0,
        'development_root_files': [r['root'] for r in groups['development']],
        'environment': {'python': platform.python_version(), 'platform': platform.platform(),
                        **{n: importlib.metadata.version(n) for n in ('numpy', 'scipy', 'uproot', 'awkward')}},
        'frozen_files': {str(p.relative_to(PROJECT)): digest(p) for p in paths},
        'snapshot_files': {str(p.relative_to(experiment)): digest(p) for p in snapshot.rglob('*') if p.is_file()},
    }
    manifest = json.loads((experiment/'manifest.json').read_text())
    assert not any(r['role'] == 'test_single' for r in manifest['runs'])
    save(output/'freeze.json', frozen)
    (experiment/'analyze_frozen.py').write_text('''from pathlib import Path
import sys
folder = Path(__file__).resolve().parent
sys.path.insert(0, str(folder / 'inference-source'))
from analysis import alpha_gamma_advanced_study as module
module.PROJECT = folder.parents[2]
module.main()
''')
    (experiment/'run_test_local.py').write_text('''from pathlib import Path
import argparse
import sys
folder = Path(__file__).resolve().parent
sys.path.insert(0, str(folder / 'inference-source'))
from analysis import alpha_gamma_advanced_native_local as module
module.PROJECT = folder.parents[2]
p = argparse.ArgumentParser()
p.add_argument('--workers', type=int, default=4)
a = p.parse_args()
module.run(folder, workers=a.workers)
''')
    print('FROZEN_ADVANCED', primary, 'selected post', selected_window, flush=True)
    return frozen


def paired_comparison(results, primary, post=150.):
    rows = []; rng = np.random.default_rng(262020909)
    for particle in ('alpha', 'gamma'):
        events = [r for r in results if r['truth_for_evaluation_only']['particle'] == particle
                  and r['truth_for_evaluation_only']['charge_pe'] >= 500.]
        for control in ('legacy_gls', 'expanded_gls'):
            differences = []
            for event in events:
                p = next(p for p in event['prefixes'] if p['post_ns'] == post)
                def success(method):
                    if not p['evaluation']['matched_signal']:
                        return False
                    ev = p['methods'][method]['evaluation']
                    return abs(ev['relative_charge_error']) <= .05 and ev['correct_type']
                differences.append(int(success(primary))-int(success(control)))
            a = np.asarray(differences, float)
            distribution = a[rng.integers(0, len(a), size=(2000, len(a)))].mean(axis=1) if len(a) else []
            rows.append({'particle': particle, 'post_ns': post, 'primary': primary, 'control': control,
                         'independent_pairs': len(a), 'joint_success_difference': float(a.mean()) if len(a) else None,
                         'paired_bootstrap95': np.quantile(distribution, [.025, .975]).tolist() if len(a) else None})
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=('calibrate', 'prepare', 'development-old', 'development-ready', 'development', 'freeze', 'test', 'summarize'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--experiment', type=Path)
    p.add_argument('--workers', type=int, default=4)
    args = p.parse_args(); output = args.output.resolve()
    if args.stage == 'calibrate':
        calibrate(output)
    elif args.stage == 'prepare':
        prepare(output)
    elif args.stage == 'development-old':
        rows = json.loads((PROJECT/OLD/'input-records.json').read_text())['development']
        checked_records(rows)
        scan(rows, output, output/'development', workers=args.workers)
    elif args.stage == 'development-ready':
        scan(ready_development(output), output, output/'development', workers=args.workers)
    elif args.stage == 'development':
        rows = json.loads((output/'input-records.json').read_text())['development']
        scan(rows, output, output/'development', workers=args.workers)
    elif args.stage == 'freeze':
        freeze(output, args.experiment)
    elif args.stage == 'test':
        settings = json.loads((output/'freeze.json').read_text())
        for relative, expected in settings['frozen_files'].items():
            assert digest(PROJECT/relative) == expected, relative
        for name, expected in settings['calibration_sha256'].items():
            assert digest(output/'calibration'/name) == expected, name
        rows = json.loads((output/'native-records.json').read_text())['test_single']
        assert len(rows) == 168
        training = json.loads((output/'input-records.json').read_text())
        assert not {r['seed'] for r in rows}.intersection(r['seed'] for group in training.values() for r in group)
        assert not {r['root'] for r in rows}.intersection(r['root'] for group in training.values() for r in group)
        scan(rows, output, output/'test', workers=args.workers)
    else:
        results = load_results(output/'test'); settings = json.loads((output/'freeze.json').read_text())
        assert len(results) == 168
        save(output/'test-summary.json', summary(results))
        save(output/'paired-comparisons.json', paired_comparison(results, settings['primary_method']))
        print('SUMMARIZED', len(results), flush=True)


if __name__ == '__main__':
    main()
