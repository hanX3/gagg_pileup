"""Repeat the pure-alpha ML comparison with current fitting and fresh seeds.

Run from the project root. Existing weights remain frozen. The public 23-window
benchmark is reported separately from newly acquired confirmation windows.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from datetime import datetime
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import time

import numpy as np

from .alpha_correlated_fit import CorrelatedConfig, ResidualCalibration, fit_correlated
from .alpha_full_signal_study import make_template, summarize
from .alpha_neural_refit import NeuralProposalConfig, fit_neural_wls
from .alpha_pileup_fit import FitConfig, fit_waveform, score
from .alpha_rate_campaign import digest, save


PROJECT = Path(os.environ.get('GAGG_ML_REFRESH_PROJECT', Path(__file__).resolve().parents[2])).resolve()
EXPERIMENT = PROJECT / 'pileup/experiments/alpha_ml_refresh_20261021'
OUTPUT = PROJECT / 'jupyter/artifacts/alpha_ml_refresh_20261021'
RECORD = PROJECT / 'jupyter/run_records/alpha_ml_refresh_20261021'
PARENT = PROJECT / 'jupyter/artifacts/alpha_correlated_global_20261020'
LEGACY = PROJECT / 'jupyter/artifacts/alpha_ml_multimethod_20261014_180000'
NETWORKS = ('cnn', 'unet', 'transformer', 'gru')
METHODS = ('baseline', 'diagonal', *NETWORKS, 'cnn_physical', 'cnn_wls')
REFERENCE_METHODS = ('local', 'correlated', 'global')


def load(path):
    return json.loads(Path(path).read_text())


def configuration():
    settings = load(EXPERIMENT / 'assets/baseline-settings.json')
    physical = CorrelatedConfig(**load(EXPERIMENT / 'assets/physical-selection.json')['config'])
    covariance = ResidualCalibration(**load(EXPERIMENT / 'assets/covariance.json')['parameters'])
    return settings, physical, covariance


def prepare():
    import torch
    assert not EXPERIMENT.exists() or {p.name for p in EXPERIMENT.iterdir()} == {'replay.py'}
    EXPERIMENT.mkdir(parents=True, exist_ok=True)
    for folder in (OUTPUT, RECORD):
        folder.mkdir(parents=True, exist_ok=False)
    for name in ('assets', 'models', 'macros', 'logs', 'bin'):
        (EXPERIMENT / name).mkdir()
    for name in ('pilots', 'validation', 'common', 'confirmation', 'controls', 'inputs'):
        (OUTPUT / name).mkdir()
    for name in ('jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb',
                 'jupyter/Run_Catalog.ipynb', 'docs/literature-use.md', 'README.md', 'pileup/README.md'):
        shutil.copyfile(PROJECT / name, RECORD / ('before__' + name.replace('/', '__')))
    for name in ('cnn', 'unet', 'transformer'):
        shutil.copytree(LEGACY / 'models' / name, EXPERIMENT / 'models' / name)
        origin = 'origin-frozen-model.json' if name == 'cnn' else 'frozen.json'
        shutil.copyfile(EXPERIMENT / 'models' / name / origin,
                        EXPERIMENT / 'models' / name / 'frozen-model.json')
    shutil.copytree(PARENT / 'gru', EXPERIMENT / 'models/gru')
    for origin, name in [
        (PARENT / 'baseline-settings.json', 'baseline-settings.json'),
        (PARENT / 'calibration/covariance.json', 'covariance.json'),
        (LEGACY / 'hybrid-config.json', 'legacy-hybrid-config.json'),
        (PROJECT / 'pileup/experiments/alpha_ml_multimethod_20261014/data-split.json', 'data-split.json'),
        (PARENT / 'freeze.json', 'parent-freeze.json'),
        (LEGACY / 'frozen-methods.json', 'legacy-freeze.json'),
    ]:
        shutil.copyfile(origin, EXPERIMENT / 'assets' / name)
    selected = load(PARENT / 'freeze.json')['selected']['diagonal']
    save(EXPERIMENT / 'assets/physical-selection.json', selected)
    cnn = load(EXPERIMENT / 'models/cnn/frozen-model.json')
    save(EXPERIMENT / 'assets/neural-proposals.json', asdict(NeuralProposalConfig(threshold=cnn['threshold'])))
    old = load(PROJECT / 'pileup/experiments/alpha_correlated_global_20261020/manifest.json')
    capacity = PROJECT / 'pileup/experiments/alpha_local_capacity_20261020'
    cap = load(capacity / 'manifest.json')
    assert cap['equivalence_passed']
    shutil.copytree(capacity / 'source/g4', EXPERIMENT / 'source/g4')
    shutil.copyfile(capacity / 'bin/gagg', EXPERIMENT / 'bin/gagg')
    assert digest(EXPERIMENT / 'bin/gagg') == cap['binary_sha256']
    plans = []
    for rate, count in [(100000, 8), (1000000, 2), (10000000, 1), (20000000, 1),
                        (40000000, 1), (60000000, 1), (80000000, 1), (100000000, 1)]:
        for index in range(count):
            plans.append(dict(label=f'mlrefresh_{rate:09d}_{index:02}', rate_hz=rate,
                              seed=262030000 + len(plans), role='confirmation',
                              preflight=rate in (1000000, 10000000) and index == 0))
    seeds = {p['seed'] for p in plans}
    for path in (PROJECT / 'pileup/experiments').glob('*/manifest.json'):
        data = load(path)
        for row in data.get('runs', []) + data.get('planned_runs', []):
            assert row.get('seed') not in seeds, path
    save(EXPERIMENT / 'manifest.json', dict(
        created=datetime.now().astimezone().isoformat(), planned_runs=plans, runs=[],
        preflight_validated=False, max_workers=2, binary_sha256=cap['binary_sha256'],
        environment=old['environment'], parent=str(capacity.relative_to(PROJECT)),
        source_sha256={str(p.relative_to(EXPERIMENT)): digest(p)
                       for p in (EXPERIMENT / 'source').rglob('*') if p.is_file()},
        source_model='Unchanged p-11B effective single-alpha spectrum; no gamma',
        window_us=20., gate_us=[5., 15.], adequate_dynamic_range=True,
        root_rule='One serial G4 event per ROOT, seconds-only filename, no overwrite'))
    protocol = dict(
        created=datetime.now().astimezone().isoformat(),
        hypotheses=['Compare four frozen neural architectures and current physical fitting on identical windows',
                    'Guarded CNN + observed-pilot union, no neural-Q prior, residual add/split/merge WLS'],
        methods=list(METHODS), public_reference_methods=list(REFERENCE_METHODS),
        network_training='No retraining; preserve exact CNN/U-Net/Transformer/GRU weights and thresholds',
        physical_configuration='Inherit diagonal WLS and calibrated covariance/penalty from Notebook 0003',
        new_parameters='No tuning: CNN threshold inherited, 262-sample context from training, 250 pe proposal gate from old hybrid',
        validation='Original ten validation ROOTs; diagnostics only, no test-based hyperparameter changes',
        public_common='Existing 23 published windows; re-evaluation, not a second independent replicate',
        confirmation='16 new seeds after complete method/config/weight/source freeze',
        scoring=dict(charge_threshold_pe=500., time_tolerance_ns=10., relative_charge_tolerance=.2),
        gate_ns=[5000., 15000.], adequate_dynamic_range=True, fixed_scintillation_decay_ns=[95., 351.],
        truth_input_at_inference=[], rate_input_at_inference=False, energy_calibration=False,
        confirmation_high_rate_windows_per_rate=1, confirmation_is_small_sample=True,
        git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=PROJECT, text=True).strip(),
        environment=dict(python=platform.python_version(), platform=platform.platform(), torch=torch.__version__,
                         packages={p: importlib.metadata.version(p) for p in ('numpy', 'scipy', 'uproot', 'awkward', 'matplotlib')}))
    save(OUTPUT / 'protocol.json', protocol)
    (RECORD / 'working-changes.patch').write_bytes(subprocess.check_output(
        ['git', '-c', 'core.fileMode=false', 'diff', '--binary'], cwd=PROJECT))
    # Actual executable source is kept independently of future workspace edits.
    for package in ('analysis', 'pileup/analysis'):
        for path in (PROJECT / package).glob('*.py'):
            target = EXPERIMENT / 'inference-source' / path.relative_to(PROJECT)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)
    for package in ('analysis', 'pileup', 'pileup/analysis'):
        p = EXPERIMENT / 'inference-source' / package / '__init__.py'
        p.parent.mkdir(parents=True, exist_ok=True)
        p.touch(exist_ok=True)
    print('PREPARED', len(plans), flush=True)


def freeze():
    assert not (OUTPUT / 'freeze.json').exists()
    files = [p for sub in ('assets', 'models', 'inference-source', 'source')
             for p in (EXPERIMENT / sub).rglob('*') if p.is_file()]
    files += [EXPERIMENT / 'bin/gagg', EXPERIMENT / 'replay.py', OUTPUT / 'protocol.json']
    save(OUTPUT / 'freeze.json', dict(frozen_at=datetime.now().astimezone().isoformat(),
         frozen_files={str(p.relative_to(PROJECT)): digest(p) for p in files},
         new_test_tuning=False, network_retraining=False))
    print('FROZEN', flush=True)


def check_frozen():
    f = load(OUTPUT / 'freeze.json')
    for path, expected in f['frozen_files'].items():
        assert digest(PROJECT / path) == expected, path
    return f


def native(workers=2):
    check_frozen()
    from . import alpha_correlated_native as runner
    from .root_redundant_photon_fields import exact_redundant_fields
    runner.PROJECT, runner.EXPERIMENT, runner.OUTPUT, runner.RECORD = PROJECT, EXPERIMENT, OUTPUT, RECORD
    original = runner.export_window

    def export(root, output, *, expected):
        try:
            return original(root, output, expected=expected)
        except Exception as error:
            if 'decompress' not in str(error).lower():
                raise
            checks = dict(original_error=repr(error), exact_redundant_fields_only=True)
            with exact_redundant_fields(original, checks):
                result = original(root, output, expected=expected)
            path = RECORD / (Path(root).stem + '-exact-fields.json')
            save(path, checks)
            expected['redundant_fields_recovery'] = str(path.relative_to(PROJECT))
            return result

    runner.export_window = export
    runner.run(workers)


def networks():
    import torch
    from .alpha_pileup_ml import EventNet
    from .alpha_ml_models import make_model
    from .alpha_correlated_gru import GRUEvent
    torch.set_num_threads(1)
    result = {}
    for name in NETWORKS:
        folder = EXPERIMENT / 'models' / name
        frozen = load(folder / 'frozen-model.json')
        assert digest(folder / 'model.pt') == frozen['weights_sha256']
        model = EventNet() if name == 'cnn' else GRUEvent() if name == 'gru' else make_model(name)
        model.load_state_dict(torch.load(folder / 'model.pt', map_location='cpu', weights_only=True))
        model.eval()
        result[name] = (model, float(frozen['threshold']))
    return result


def records(phase):
    if phase == 'common':
        return load(PARENT / 'test-summary.json')['windows']
    if phase == 'validation':
        result = []
        for r in load(EXPERIMENT / 'assets/data-split.json')['windows']:
            if r['split'] != 'validation':
                continue
            result.append(dict(root=r['root'], root_sha256=r['root_sha256'], seed=r['seed'],
                rate_hz=r['rate_hz'], waveform=r['waveforms'], truth=r['truth'],
                gate_ns=load(PROJECT / r['metadata'])['gate_ns']))
        return result
    manifest = load(EXPERIMENT / 'manifest.json')
    assert len(manifest['runs']) == 16 and all(r['status'] == 'completed' for r in manifest['runs'])
    return [dict(root=r['root'], root_sha256=r['root_sha256'], seed=r['seed'], rate_hz=r['rate_hz'],
                 waveform=str((OUTPUT / 'inputs/waveforms' / (Path(r['root']).stem + '.npz')).relative_to(PROJECT)),
                 truth=str((OUTPUT / 'inputs/truth' / (Path(r['root']).stem + '.npz')).relative_to(PROJECT)),
                 gate_ns=[5000., 15000.]) for r in manifest['runs']]


def measured(row):
    path = PROJECT / row['waveform']
    with np.load(path) as data:
        return data['time_ns'].copy(), data['analog_mV'].copy()


def pilot(row, t, y):
    stem = Path(row['root']).stem
    cached = OUTPUT / 'pilots' / (stem + '.json')
    waveform = PROJECT / row['waveform']
    if cached.exists():
        result = load(cached)
        assert result['voltage_sha256'] == digest(waveform)
        result['fit']['model_mV'] = np.load(cached.with_suffix('.npz'))['model_mV'].copy()
        return result['fit'], result['seconds'], result['origin']
    old = PARENT / 'pilots' / (stem + '.json')
    if old.exists() and load(old)['voltage_sha256'] == digest(waveform):
        result = load(old)
        initial = result['fit']
        initial['model_mV'] = np.load(old.with_suffix('.npz'))['model_mV'].copy()
        elapsed, origin = result['seconds'], str(old.relative_to(PROJECT))
    else:
        settings, _, _ = configuration()
        tick = time.monotonic()
        initial = fit_waveform(t, y, make_template(settings), FitConfig(**settings['fit_config']))
        elapsed, origin = time.monotonic() - tick, 'fresh measured-voltage inference'
    np.savez_compressed(cached.with_suffix('.npz'), model_mV=initial['model_mV'])
    save(cached, dict(fit={k: v for k, v in initial.items() if k != 'model_mV'}, seconds=elapsed,
                     voltage_sha256=digest(waveform), origin=origin))
    return initial, elapsed, origin


def inference(row, nets, phase):
    from .alpha_pileup_ml import waveform_features, dense_prediction, decode_prediction
    from .alpha_ml_models import windowed_prediction
    from .alpha_pileup_hybrid import HybridConfig, hybrid_waveform
    t, y = measured(row)
    settings, physical, covariance = configuration()
    template, base = make_template(settings), FitConfig(**settings['fit_config'])
    tick = time.monotonic()
    features = waveform_features(t, y, settings)
    feature_seconds = time.monotonic() - tick
    fits, seconds, dense = {}, {}, {}
    for name, (model, threshold) in nets.items():
        tick = time.monotonic()
        prediction = dense_prediction(model, features) if name == 'cnn' else windowed_prediction(model, features)
        fits[name] = decode_prediction(t, prediction, threshold)
        seconds[name] = feature_seconds + time.monotonic() - tick
        dense[name] = prediction
    initial, elapsed, origin = pilot(row, t, y)
    fits['baseline'], seconds['baseline'] = initial, elapsed
    if phase == 'common':
        for name in ('diagonal', *REFERENCE_METHODS):
            old = row['methods'][name]
            fits[name], seconds[name] = dict(old['fit']), old['seconds']
            source = (PARENT / 'test/windows' / (Path(row['root']).stem + '_local.npz')
                      if name == 'local' else PARENT / 'test' / name / (Path(row['root']).stem + '.npz'))
            fits[name]['model_mV'] = np.load(source)['model_mV'].copy()
    else:
        tick = time.monotonic()
        fits['diagonal'] = fit_correlated(t, y, template, base, covariance, physical, initial=initial)
        seconds['diagonal'] = elapsed + time.monotonic() - tick
    legacy_config = load(EXPERIMENT / 'assets/legacy-hybrid-config.json')
    old_config = HybridConfig(**{k: v for k, v in legacy_config.items() if k in HybridConfig.__dataclass_fields__})
    tick = time.monotonic()
    fits['cnn_physical'] = hybrid_waveform(t, y, dense['cnn'], settings, old_config)
    seconds['cnn_physical'] = seconds['cnn'] + time.monotonic() - tick
    tick = time.monotonic()
    proposals = NeuralProposalConfig(**load(EXPERIMENT / 'assets/neural-proposals.json'))
    fits['cnn_wls'] = fit_neural_wls(t, y, dense['cnn'], template, base, covariance, physical, initial, proposals)
    seconds['cnn_wls'] = elapsed + seconds['cnn'] + time.monotonic() - tick
    return fits, seconds, dense, origin


def evaluate(phase, workers=1):
    check_frozen()
    nets = networks()
    folder = OUTPUT / phase
    for name in ('windows', 'models', 'dense'):
        (folder / name).mkdir(exist_ok=True)

    def one(row):
        stem = Path(row['root']).stem
        destination = folder / 'windows' / (stem + '.json')
        if destination.exists():
            result = load(destination)
            assert result['waveform_sha256'] == digest(PROJECT / row['waveform'])
            return result
        fits, elapsed, dense, origin = inference(row, nets, phase)
        # Independent scoring labels are not opened until every blind fit ends.
        with np.load(PROJECT / row['truth']) as data:
            truth = {k: data[k].copy() for k in ('first_active_ns', 'emission_ns', 'charge_pe')}
        methods = {}
        for name, fit in fits.items():
            s = json.loads(json.dumps(score(fit, truth, row['gate_ns'])))
            if phase == 'common' and name in row['methods']:
                assert s == row['methods'][name]['score'], (stem, name, 'original score changed')
            if 'model_mV' in fit:
                np.savez_compressed(folder / 'models' / (stem + '_' + name + '.npz'), model_mV=fit['model_mV'])
            methods[name] = dict(fit={k: v for k, v in fit.items() if k != 'model_mV'}, score=s, seconds=elapsed[name])
        for name, prediction in dense.items():
            np.savez_compressed(folder / 'dense' / (stem + '_' + name + '.npz'), **prediction)
        result = {k: v for k, v in row.items() if k != 'methods'}
        result.update(methods=methods, waveform_sha256=digest(PROJECT / row['waveform']),
                      truth_sha256=digest(PROJECT / row['truth']), pilot_origin=origin,
                      inference_truth_inputs=[], rate_input_at_inference=False, phase=phase)
        save(destination, result)
        print('EVALUATED', phase, stem, [(m, methods[m]['score']['recovered_count']) for m in METHODS], flush=True)
        return result

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(one, row) for row in records(phase)]
        rows = [future.result() for future in as_completed(futures)]
    methods = [m for m in (*METHODS, *REFERENCE_METHODS) if all(m in r['methods'] for r in rows)]
    rates = {m: summarize([dict(rate_hz=r['rate_hz'], fit=r['methods'][m]['fit'], score=r['methods'][m]['score'])
                            for r in rows]) for m in methods}
    result = dict(phase=phase, windows=sorted(rows, key=lambda r: (r['rate_hz'], r['seed'])),
                  rates=rates, method_order=methods, adequate_dynamic_range=True, energy_calibration=False,
                  inference_truth_inputs=[], frozen_at=load(OUTPUT / 'freeze.json')['frozen_at'])
    save(folder / 'summary.json', result)
    print('SUMMARY', phase, {m: [sum(r[k] for r in rates[m]) for k in ('recovered_count', 'target_count', 'unmatched_predictions')]
                             for m in methods}, flush=True)


def shifts(phase):
    study = load(OUTPUT / phase / 'summary.json')
    deltas = [1379., 2933., 4871., 7039., 8911.]
    rows = []
    for row in study['windows']:
        with np.load(PROJECT / row['truth']) as z:
            truth = {k: z[k].copy() for k in ('first_active_ns', 'emission_ns', 'charge_pe')}
        for name, data in row['methods'].items():
            for delta in deltas:
                fit = dict(data['fit'])
                times = np.asarray(fit['times_ns'], float)
                central = (times >= 5000.) & (times < 15000.)
                times[central] = 5000. + (times[central] - 5000. + delta) % 10000.
                fit['times_ns'] = times.tolist()
                rows.append(dict(root=row['root'], rate_hz=row['rate_hz'], method=name,
                                 shift_ns=delta, score=score(fit, truth, row['gate_ns'])))
    save(OUTPUT / phase / 'timing-null.json', dict(shifts_ns=deltas, rows=rows,
         interpretation='Same-window accidental-match reference; not exact FDR, not independent replicas, not subtracted efficiency'))
    print('SHIFTS', phase, len(rows), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['prepare', 'freeze', 'native', 'validation', 'common', 'confirmation', 'shifts'])
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--phase', choices=['common', 'confirmation'], default='confirmation')
    args = parser.parse_args()
    if args.stage in ('prepare', 'freeze'):
        globals()[args.stage]()
    elif args.stage == 'native':
        native(args.workers)
    elif args.stage == 'shifts':
        shifts(args.phase)
    else:
        evaluate(args.stage, args.workers)


if __name__ == '__main__':
    main()
