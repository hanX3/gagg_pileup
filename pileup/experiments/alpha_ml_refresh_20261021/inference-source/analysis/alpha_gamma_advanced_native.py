"""Native isolated-pulse calibration/development/test data for advanced fits.

Each serial Geant4 process writes one ROOT. Run starts are at least one second
apart, and test production requires an immutable inference/calibration freeze.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess

from pileup.analysis import alpha_gamma_regularized_campaign as native
from .alpha_gamma_short_window_study import digest, save

PROJECT = Path(__file__).resolve().parents[1]


def prepare(folder, output):
    folder, output = Path(folder).resolve(), Path(output).resolve()
    if folder.exists():
        assert not (folder / 'manifest.json').exists(), 'Preparation already completed'
        assert not list((folder / 'macros').glob('*.mac')), 'Production already started'
    folder.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    assert not (output / 'production-plan.json').exists()
    parent = PROJECT / 'pileup/experiments/alpha_gamma_regularized_fit_20261015'
    previous = json.loads((parent / 'manifest.json').read_text())
    assert digest(parent / 'bin/gagg-brems') == previous['parents']['brems']['binary_sha256']
    for sub in ('macros', 'logs', 'bin', 'source/analysis', 'source/pileup/analysis'):
        (folder / sub).mkdir(parents=True, exist_ok=True)
    shutil.copytree(parent / 'source/brems/g4', folder / 'source/brems/g4', dirs_exist_ok=True)
    shutil.copy2(parent / 'bin/gagg-brems', folder / 'bin/gagg-brems')
    files = list((PROJECT / 'analysis').glob('*.py'))
    files += [PROJECT / 'pileup/analysis' / n for n in (
        '__init__.py', 'alpha_gamma_regularized_campaign.py', 'alpha_rate_campaign.py',
        'alpha_rate_data.py', 'sipm_response.py', 'sipm_readout.py',
    )]
    files.append(PROJECT / 'requirements-analysis.lock.txt')
    for source in files:
        target = folder / 'source' / source.relative_to(PROJECT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    (folder / 'source/pileup/__init__.py').write_text('')
    (folder / 'source/analysis/__init__.py').write_text('')
    plans = []
    for role, alpha_repeats, gamma_repeats in (
        ('calibration_single', 2, 12), ('development_single', 1, 6), ('test_single', 2, 12),
    ):
        for particle, energies, repetitions in (
            ('alpha', (1.5, 2.5, 5.486, 8.), alpha_repeats),
            ('gamma', (.1, .2, .661657, 1.5), gamma_repeats),
        ):
            for energy in energies:
                for x, y in ((0., 0.), (4., 0.), (4., 4.)):
                    for repetition in range(repetitions):
                        plans.append({
                            'label': f'{role}_{particle}_{energy:g}_{x:g}_{y:g}_{repetition:02}',
                            'role': role, 'particle': particle, 'energy_MeV': energy,
                            'x_mm': x, 'y_mm': y, 'gamma_model': 'brems',
                            'rate_hz': 0, 'gamma_rate_hz': 0, 'seed': 262020000 + len(plans),
                        })
    for particle, energy in (('alpha', 1.5), ('alpha', 8.), ('gamma', .1), ('gamma', .661657)):
        next(p for p in plans if p['role'] == 'calibration_single' and p['particle'] == particle
             and p['energy_MeV'] == energy and p['x_mm'] == p['y_mm'] == 0.)['preflight'] = True
    seeds = {p['seed'] for p in plans}
    for path in (PROJECT / 'pileup/experiments').glob('*/manifest.json'):
        existing = json.loads(path.read_text())
        assert not seeds.intersection(r['seed'] for r in existing.get('planned_runs', []) if 'seed' in r), path
    save(folder / 'manifest.json', {
        'created': datetime.now().astimezone().isoformat(),
        'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=PROJECT, text=True).strip(),
        'purpose': 'New native calibration/development plus post-freeze isolated test for advanced conditional fits',
        'parents': {'brems': previous['parents']['brems']},
        'planned_runs': plans, 'runs': [], 'preflight_validated': False,
        'window_us': 20., 'gate_us': [5., 15.], 'max_workers': 8,
        'analysis_output': str(output.relative_to(PROJECT)),
        'environment': {'python': platform.python_version(), 'platform': platform.platform(),
                        **{n: importlib.metadata.version(n) for n in ('numpy', 'scipy', 'uproot', 'awkward')}},
        'input_policy': 'Complete analog voltage, adequate dynamic range; retain all gamma zero-signal windows',
        'rates': 'One primary per ROOT, not a Poisson source-rate experiment; rate fields unused in single mode',
        'source_sha256': {str(p.relative_to(folder)): digest(p) for p in (folder / 'source').rglob('*') if p.is_file()},
    })
    (folder / 'git-status-before.txt').write_text(subprocess.check_output(
        ['git', '-c', 'core.fileMode=false', 'status', '--short'], cwd=PROJECT, text=True))
    (folder / 'working-changes.patch').write_bytes(subprocess.check_output(
        ['git', '-c', 'core.fileMode=false', 'diff', '--binary'], cwd=PROJECT))
    (folder / 'run_frozen.py').write_text('''from pathlib import Path
import sys
folder = Path(__file__).resolve().parent
sys.path.insert(0, str(folder / 'source'))
from analysis import alpha_gamma_advanced_native as module
module.PROJECT = folder.parents[2]
module.main()
''')
    save(output / 'production-plan.json', {
        'native_experiment': str(folder.relative_to(PROJECT)),
        'planned_counts': {role: {particle: sum(p['role'] == role and p['particle'] == particle for p in plans)
                                  for particle in ('alpha', 'gamma')}
                           for role in ('calibration_single', 'development_single', 'test_single')},
        'test_policy': 'Do not generate test primaries until code, calibration, methods and selection are frozen',
    })
    print('PREPARED_NATIVE', len(plans), flush=True)


def run(folder, stage, workers):
    folder = Path(folder).resolve()
    manifest = json.loads((folder / 'manifest.json').read_text())
    output = PROJECT / manifest['analysis_output']
    if stage == 'test_single':
        freeze = json.loads((output / 'freeze.json').read_text())
        for relative, expected in freeze['frozen_files'].items():
            assert digest(PROJECT / relative) == expected, relative
        assert (folder / 'inference-source/analysis/alpha_gamma_advanced_fit.py').is_file()
        if not (folder / 'frozen-fit-settings.json').exists():
            save(folder / 'frozen-fit-settings.json', freeze)
    native.PROJECT = PROJECT
    original_macro = native.macro
    native.macro = lambda f, row, destination: original_macro(f, row, destination).replace(
        '/gun/position 0 0 0 mm', f"/gun/position {row['x_mm']} {row['y_mm']} 0 mm")
    native.run(folder, stage, output / 'native-inputs', workers=workers)
    manifest = json.loads((folder / 'manifest.json').read_text())
    groups = {role: [] for role in ('calibration_single', 'development_single', 'test_single')}
    for row in manifest['runs']:
        if row['status'] != 'completed':
            continue
        stem = Path(row['root']).stem
        groups[row['role']].append({
            'root': row['root'], 'seed': row['seed'], 'particle': row['particle'],
            'source_manifest': str((folder / 'manifest.json').relative_to(PROJECT)),
            'root_original_sha256': row['root_sha256'],
            **{name: str((output / 'native-inputs' / sub / (stem + suffix)).relative_to(PROJECT))
               for name, sub, suffix in (('waveform', 'waveforms', '.npz'), ('truth', 'truth', '.npz'),
                                        ('metadata', 'metadata', '.json'))},
        })
    save(output / 'native-records.json', groups)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=('prepare', 'preflight', 'train', 'test_single'))
    p.add_argument('--folder', type=Path, required=True)
    p.add_argument('--output', type=Path)
    p.add_argument('--workers', type=int, default=8)
    args = p.parse_args()
    if args.stage == 'prepare':
        prepare(args.folder, args.output)
    elif args.stage == 'train':
        run(args.folder, 'calibration_single', args.workers)
        run(args.folder, 'development_single', args.workers)
    else:
        run(args.folder, args.stage, args.workers)


if __name__ == '__main__':
    main()
