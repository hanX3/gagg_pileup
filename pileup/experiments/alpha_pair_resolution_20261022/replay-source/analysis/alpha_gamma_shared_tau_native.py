"""Independent post-freeze native test for lifetime-sharing fits.

The unchanged G4 binary/configuration are snapshotted, executed serially per
process, with regular local binary/log caches and ROOT files on the shared disk.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import time

from pileup.analysis import alpha_gamma_regularized_campaign as native
from .alpha_gamma_short_window_study import digest, save

PROJECT = Path(__file__).resolve().parents[1]


def prepare(folder, output):
    folder, output = Path(folder).resolve(), Path(output).resolve()
    assert not (folder/'manifest.json').exists()
    parent = PROJECT/'pileup/experiments/alpha_gamma_advanced_fit_20261020'
    previous = json.loads((parent/'manifest.json').read_text())
    for sub in ('bin', 'macros', 'logs'):
        (folder/sub).mkdir(parents=True, exist_ok=True)
    shutil.copytree(parent/'source/brems/g4', folder/'source/brems/g4')
    shutil.copy2(parent/'bin/gagg-brems', folder/'bin/gagg-brems')
    plans = []
    for particle, energies, repetitions in (
        ('alpha', (1.5, 2.5, 5.486, 8.), 2), ('gamma', (.1, .2, .661657, 1.5), 12)):
        for energy in energies:
            for x, y in ((0., 0.), (4., 0.), (4., 4.)):
                for i in range(repetitions):
                    plans.append(dict(label=f'test_{particle}_{energy:g}_{x:g}_{y:g}_{i:02}',
                        role='test_single', particle=particle, energy_MeV=energy,
                        x_mm=x, y_mm=y, gamma_model='brems', rate_hz=0, gamma_rate_hz=0,
                        seed=262021000+len(plans)))
    for particle, energy in (('alpha', 1.5), ('alpha', 8.), ('gamma', .1), ('gamma', .661657)):
        next(r for r in plans if r['particle'] == particle and r['energy_MeV'] == energy
             and r['x_mm'] == r['y_mm'] == 0.)['preflight'] = True
    seeds = {r['seed'] for r in plans}
    for path in (PROJECT/'pileup/experiments').glob('*/manifest.json'):
        previous_manifest = json.loads(path.read_text())
        assert not seeds.intersection(r['seed'] for r in previous_manifest.get('planned_runs', []) if 'seed' in r), path
    manifest = dict(created=datetime.now().astimezone().isoformat(),
        git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=PROJECT, text=True).strip(),
        purpose='Post-freeze independent native isolated test for fixed/shared/constrained/free double lifetimes',
        parents={'brems': previous['parents']['brems']}, planned_runs=plans, runs=[],
        window_us=20., gate_us=[5., 15.], max_workers=4, preflight_validated=False,
        analysis_output=str(output.relative_to(PROJECT)),
        environment=dict(python=platform.python_version(), platform=platform.platform(),
                         packages={n: importlib.metadata.version(n) for n in ('numpy', 'scipy', 'uproot', 'awkward')}),
        input_policy='Full analog voltage, adequate dynamic range; no clipping rejection; all zero-signal gamma retained',
        rate_policy='One primary per ROOT, unused rate fields are not a source-rate experiment',
        binary_sha256=digest(folder/'bin/gagg-brems'),
        source_sha256={str(p.relative_to(folder)): digest(p) for p in (folder/'source').rglob('*') if p.is_file()})
    assert manifest['binary_sha256'] == previous['parents']['brems']['binary_sha256']
    save(folder/'manifest.json', manifest)
    (folder/'run_native.py').write_text('''from pathlib import Path
import sys
folder = Path(__file__).resolve().parent
sys.path.insert(0, str(folder/'inference-source'))
from analysis import alpha_gamma_shared_tau_native as module
module.PROJECT = folder.parents[2]
if __name__ == '__main__':
    module.main()
''')
    (folder/'analyze_frozen.py').write_text('''from pathlib import Path
import sys
folder = Path(__file__).resolve().parent
sys.path.insert(0, str(folder/'inference-source'))
from analysis import alpha_gamma_shared_tau_study as module
module.PROJECT = folder.parents[2]
if __name__ == '__main__':
    module.main()
''')
    save(output/'production-plan.json', dict(native_experiment=str(folder.relative_to(PROJECT)),
        planned_windows=168, alpha=24, gamma=144, seed_base=262021000,
        first_four_checks_count_inside_test=True, all_production_after_freeze=True))


def run(folder, workers=4):
    folder = Path(folder).resolve()
    manifest = json.loads((folder/'manifest.json').read_text())
    output = PROJECT/manifest['analysis_output']
    frozen = json.loads((output/'freeze.json').read_text())
    assert 1 <= workers <= 4
    for relative, expected in frozen['frozen_files'].items():
        assert digest(PROJECT/relative) == expected, relative
    assert all(r['status'] == 'completed' for r in manifest['runs'])
    native.PROJECT = PROJECT
    complete = {r['label'] for r in manifest['runs']}
    plans = [r for r in manifest['planned_runs'] if r['label'] not in complete]
    plans.sort(key=lambda r: (not r.get('preflight', False), r['seed']))
    cache = Path('/tmp/gagg-shared-tau-native-20261020')/manifest['binary_sha256']
    cache.mkdir(parents=True, exist_ok=True)
    binary = cache/'gagg-brems'
    if not binary.exists():
        shutil.copy2(folder/'bin/gagg-brems', binary)
    assert digest(binary) == manifest['binary_sha256']
    env = dict(os.environ); env.update(manifest['parents']['brems']['environment']['g4_dataset_environment'])
    manifest['runner'] = dict(binary_cache=str(binary), no_symbolic_links=True, workers=workers,
                             stdout='Regular local file, then byte copy into retained experiment logs')
    active = []; last = -float('inf')
    def persist():
        save(folder/'manifest.json', manifest)
    try:
        while plans or active:
            ready = plans and (manifest['preflight_validated'] or plans[0].get('preflight', False))
            if ready and len(active) < workers and time.monotonic()-last >= 1.05:
                root = PROJECT/'g4/data'/datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss.root')
                if root.exists():
                    time.sleep(.2); continue
                planned = plans.pop(0)
                row = {**planned, 'root': str(root.relative_to(PROJECT)), 'status': 'running',
                       'sipm_seed': planned['seed']+100000, 'readout_seed': planned['seed']+200000,
                       'started': datetime.now().astimezone().isoformat()}
                macro = folder/'macros'/(row['label']+'.mac')
                log = folder/'logs'/(row['label']+'.log')
                local_log = cache/(root.stem+'.log')
                assert not macro.exists() and not log.exists() and not local_log.exists()
                text = native.macro(folder, row, root).replace('/gun/position 0 0 0 mm',
                    f"/gun/position {row['x_mm']} {row['y_mm']} 0 mm")
                macro.write_text(text)
                stream = local_log.open('x')
                proc = subprocess.Popen([str(binary), str(macro)], cwd=PROJECT/'g4/build', env=env,
                                        stdout=stream, stderr=subprocess.STDOUT)
                row.update(macro=str(macro.relative_to(PROJECT)), macro_sha256=digest(macro),
                           log=str(log.relative_to(PROJECT)), actual_stdout_file=str(local_log),
                           command=[str(binary), str(macro)], cwd='g4/build', pid=proc.pid,
                           binary_sha256=manifest['binary_sha256'])
                tick = time.monotonic(); last = tick
                manifest['runs'].append(row); active.append((proc, stream, row, tick)); persist()
                print('NATIVE_START', row['label'], row['root'], flush=True)
            for task in list(active):
                proc, stream, row, tick = task
                if proc.poll() is None:
                    continue
                stream.close(); active.remove(task)
                shutil.copy2(row['actual_stdout_file'], PROJECT/row['log'])
                text = Path(row['actual_stdout_file']).read_text(); root = PROJECT/row['root']
                ok = proc.returncode == 0 and root.exists() and 'ROOT file is saved' in text
                ok = ok and not any(s in text for s in ('COMMAND NOT FOUND', 'Error in <TFile', 'Error in <TTree', 'Error in <TBranch'))
                row.update(returncode=proc.returncode, elapsed_s=time.monotonic()-tick,
                           status='completed' if ok else 'failed', finished=datetime.now().astimezone().isoformat())
                if root.exists():
                    row.update(root_sha256=digest(root), root_bytes=root.stat().st_size)
                persist(); assert ok, row['log']
                meta = native.export_window(root, output/'native-inputs', allowed_particles=('alpha', 'gamma'))
                assert meta['seed'] == row['seed'] and meta['macro_contents'] == (PROJECT/row['macro']).read_text()
                row['validation'] = dict(macro_verified=True, charge_conserved=True,
                                         primaries=meta['primaries'], charge_pe=meta['charge_pe'])
                done_preflight = [r for r in manifest['runs'] if r.get('preflight') and r['status'] == 'completed' and r.get('validation')]
                if len(done_preflight) == 4:
                    manifest['preflight_validated'] = True
                persist(); print('NATIVE_DONE', row['label'], round(row['elapsed_s'], 2), 's', flush=True)
            if plans or active:
                time.sleep(.2)
    except BaseException:
        for proc, stream, row, tick in active:
            proc.terminate(); proc.wait(); stream.close()
            shutil.copy2(row['actual_stdout_file'], PROJECT/row['log'])
            row.update(status='interrupted', returncode=proc.returncode)
        persist(); raise
    rows = []
    for r in manifest['runs']:
        stem = Path(r['root']).stem
        rows.append(dict(root=r['root'], seed=r['seed'], particle=r['particle'],
            source_manifest=str((folder/'manifest.json').relative_to(PROJECT)), root_original_sha256=r['root_sha256'],
            **{name: str((output/'native-inputs'/sub/(stem+suffix)).relative_to(PROJECT))
               for name, sub, suffix in (('waveform', 'waveforms', '.npz'), ('truth', 'truth', '.npz'), ('metadata', 'metadata', '.json'))}))
    save(output/'native-records.json', {'test_single': rows})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=('prepare', 'run'))
    p.add_argument('--folder', required=True, type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--workers', default=4, type=int)
    a = p.parse_args()
    if a.stage == 'prepare':
        prepare(a.folder, a.output)
    else:
        run(a.folder, a.workers)


if __name__ == '__main__':
    main()
