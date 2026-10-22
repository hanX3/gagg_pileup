"""Controlled separation scans, with immutable models and independent transports.

Run from the project root. Truth and observations are saved in separate files.
The shared-cell domain is joint offline digitization of fresh G4 photons; the
linear domain explicitly interpolates separately digitized components.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import replace
from datetime import datetime
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import threading
import time

import numpy as np

PROJECT = Path(os.environ.get('GAGG_RESOLUTION_PROJECT', Path(__file__).resolve().parents[2])).resolve()
NAME = 'alpha_pair_resolution_20261022'
EXPERIMENT = PROJECT / 'pileup/experiments' / NAME
OUTPUT = PROJECT / 'jupyter/artifacts' / NAME
RECORD = PROJECT / 'jupyter/run_records' / NAME
PARENT = PROJECT / 'pileup/experiments/alpha_pair_curriculum_20261021'
PARENT_OUTPUT = PROJECT / 'jupyter/artifacts/alpha_pair_curriculum_20261021'
GAPS = (0., 2., 4., 6., 8., 10., 12., 16., 20., 30., 40., 60., 80., 120., 200., 400., 800.)
SAMPLES = 1024
GRID = np.arange(SAMPLES)*2.+.25
MODEL_NAMES = ('pair_ordered_physical_262050001', 'pair_ordered_physical_262050002')


def load(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def settings():
    return load(EXPERIMENT/'settings.json')


def protocol():
    return load(EXPERIMENT/'protocol.json')


def prepare():
    if EXPERIMENT.exists() or OUTPUT.exists() or RECORD.exists():
        raise FileExistsError('Use an independent destination; this experiment already exists.')
    for p in (EXPERIMENT, OUTPUT, RECORD):
        p.mkdir(parents=True)
    for name in ('macros', 'logs', 'bin', 'models'):
        (EXPERIMENT/name).mkdir()
    for name in ('assets', 'photons', 'figures'):
        (OUTPUT/name).mkdir()
    for name in ('codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb', 'Run_Catalog.ipynb'):
        shutil.copyfile(PROJECT/'jupyter'/name, RECORD/(name+'.before'))
    shutil.copyfile(PROJECT/'docs/literature-use.md', RECORD/'literature-use.before.md')
    shutil.copyfile(PARENT/'settings.json', EXPERIMENT/'settings.json')
    shutil.copytree(PARENT/'source', EXPERIMENT/'source')
    shutil.copyfile(PARENT/'bin/gagg-brems', EXPERIMENT/'bin/gagg-brems')
    (EXPERIMENT/'bin/gagg-brems').chmod(0o755)
    shutil.copytree(PARENT/'inference-source', EXPERIMENT/'replay-source')
    for p in (PROJECT/'pileup/analysis').glob('alpha_pair_resolution*.py'):
        shutil.copyfile(p, EXPERIMENT/'replay-source/pileup/analysis'/p.name)
    for name in MODEL_NAMES:
        folder = EXPERIMENT/'models'/name
        folder.mkdir()
        for file in ('model.pt', 'config.json', 'frozen.json', 'selection.json'):
            shutil.copyfile(PARENT_OUTPUT/'models'/name/file, folder/file)
        frozen = load(folder/'frozen.json')
        assert digest(folder/'model.pt') == frozen['weights_sha256']
        assert digest(folder/'config.json') == frozen['config_sha256']
    shutil.copyfile(PARENT_OUTPUT/'fitting-selection.json', EXPERIMENT/'fitting-selection.json')
    run = """from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for key in ('GAGG_RESOLUTION_PROJECT','GAGG_PAIR_PROJECT','GAGG_COMPONENT_PROJECT'):
    os.environ[key]=str(folder.parents[2])
sys.path.insert(0,str(folder/'replay-source'))
if __name__ == '__main__':
    if sys.argv[1] == 'fit':
        sys.argv.pop(1)
        from pileup.analysis.alpha_pair_resolution_fit import main
    else:
        from pileup.analysis.alpha_pair_resolution_data import main
    main()
"""
    (EXPERIMENT/'replay.py').write_text(run)
    conditions = [
        dict(name='reference', q1=4000., q2=4000., noise_rms=.1, ar1=0.),
        dict(name='low_equal', q1=1000., q2=1000., noise_rms=.1, ar1=0.),
        dict(name='high_equal', q1=16000., q2=16000., noise_rms=.1, ar1=0.),
        dict(name='weak_second', q1=4000., q2=1000., noise_rms=.1, ar1=0.),
        dict(name='weak_first', q1=1000., q2=4000., noise_rms=.1, ar1=0.),
        dict(name='white_1mV', q1=4000., q2=4000., noise_rms=1., ar1=0.),
        dict(name='ar1_0p98', q1=4000., q2=4000., noise_rms=.1, ar1=.98),
    ]
    p = dict(created=datetime.now().astimezone().isoformat(),
             purpose='Conditional two-event resolution; no universal deadtime inferred.',
             gaps_ns=list(GAPS), mean_conditions=conditions, repetitions_per_point=40,
             mean_single_controls_per_condition=80, mean_empty_controls_per_condition=80,
             native_energies_MeV=[2.5, 5.486, 8.0], native_roots_per_energy=80,
             shared_conditions=['low_low', 'mid_mid', 'high_high', 'low_high', 'high_low'],
             shared_replicates=40, linear_conditions=['mid_mid'],
             higher_noise_shared_conditions=['mid_mid'],
             models=list(MODEL_NAMES), models_retrained=False, thresholds_reselected=False,
             fitting_penalty=load(EXPERIMENT/'fitting-selection.json')['penalty'],
             samples=SAMPLES, sampling_ns=2., first_sample_ns=.25,
             first_photon_control_range_ns=[96., 256.],
             fixed_scintillation_taus_ns=[95., 351.], adequate_dynamic_range=True,
             truth_charge='Actual avalanche-equivalent pe; not deposited/source energy.',
             truth_time='First photon inside active plane, not pulse maximum.',
             q_threshold_pe=500., relative_charge_tolerance=.2,
             ordinary_time_tolerance_ns=10., strict_time_rule='min(10 ns, true separation/4), evaluation only',
             coincident_control='Delta0 excluded from individual-separation efficiency; sum/ambiguity diagnostic.',
             false_output_budget=.05,
             budget_rule='Report empty and single controls separately; inherited thresholds never tuned on these tests.',
             delta90='Smallest positive tested gap with point efficiency>=90% at this and every later tested gap; pointwise Wilson95 bounds separately, not simultaneous confidence.',
             independence='40 independent pairs within each native condition; same parents reused across gaps and some conditions. Never pool as independent exposure.',
             physical_groups='Source energies set conditions; realized Q and amplitude ratios must be reported separately.',
             source_rate='Controlled pairs and singles, not a Poisson rate scan.',
             mean_seed=262223000, native_seed_start=262220000, digitizer_seed_start=262230000,
             limits='Local pure-alpha unknown0/1/2; fixed tau; not mixed alpha/gamma, unknownK>=3, full20us or measured hardware.',
             development_cost_check=dict(mean_fit_seconds=.0327, joint_digitization_seconds=[.409,.430,.380],
                                         domain='Previously exposed data only; measured before this protocol.'))
    save(EXPERIMENT/'protocol.json', p)
    plans = []
    for energy in p['native_energies_MeV']:
        for i in range(p['native_roots_per_energy']):
            plans.append(dict(label=f'alpha_{energy:g}_{i:03d}', particle='alpha',
                              energy_MeV=energy, x_mm=0., y_mm=0., replicate=i,
                              role='resolution_test_donor', gamma_model='brems',
                              seed=p['native_seed_start']+len(plans)))
    # Check new G4 seeds against existing production manifests.
    seeds = {r['seed'] for r in plans}
    for manifest in (PROJECT/'pileup/experiments').glob('*/manifest.json'):
        for r in load(manifest).get('planned_runs', []):
            assert r.get('seed') not in seeds, str(manifest)
    parent = load(PARENT/'manifest.json')
    save(EXPERIMENT/'manifest.json', dict(created=p['created'], planned_runs=plans, runs=[],
                                         preflight_validated=False, max_workers=4,
                                         environment=parent['environment'], window_us=20.,
                                         binary_sha256=digest(EXPERIMENT/'bin/gagg-brems')))
    with (RECORD/'git-status-before.txt').open('w') as f:
        subprocess.run(['git', '-c', 'core.fileMode=false', 'status', '--short'], cwd=PROJECT,
                       stdout=f, check=True)
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=PROJECT, text=True).strip()
    files = {}
    for folder in ('source', 'replay-source', 'models'):
        for file in (EXPERIMENT/folder).rglob('*'):
            if file.is_file():
                files[str(file.relative_to(EXPERIMENT))] = digest(file)
    for file in ('settings.json', 'protocol.json', 'fitting-selection.json', 'replay.py', 'bin/gagg-brems'):
        files[file] = digest(EXPERIMENT/file)
    save(EXPERIMENT/'freeze.json', dict(frozen_at=datetime.now().astimezone().isoformat(),
                                      git_head=head, files=files, new_test_generation_started=False,
                                      inference_truth_inputs=[], inherited_model_training_unchanged=True))
    save(RECORD/'environment.json', dict(python=sys.version, executable=sys.executable,
                                        numpy=np.__version__, optional_torch_path=os.environ.get('PYTHONPATH'),
                                        cpu_threads=os.environ.get('OMP_NUM_THREADS')))
    print('PREPARED', len(plans), 'new ROOTs; frozen weights/source/protocol', flush=True)


def check_freeze():
    f = load(EXPERIMENT/'freeze.json')
    for path, sha in f['files'].items():
        assert digest(EXPERIMENT/path) == sha, path
    return f


def save_dataset(name, waves, charges, times, cases, components=None):
    folder = OUTPUT/'assets'
    paths = [folder/(name+suffix) for suffix in ('.waveforms.npz', '.truth.npz', '.json')]
    assert not any(p.exists() for p in paths), f'Dataset exists: {name}'
    np.savez_compressed(paths[0], voltage_mV=np.asarray(waves, np.float32))
    q = np.asarray(charges, np.float64)
    values = dict(charge_pe=q, times_ns=np.asarray(times, np.float64), active=q>=500.)
    if components is not None:
        values['components_mV'] = np.asarray(components, np.float32)
    np.savez_compressed(paths[1], **values)
    save(paths[2], dict(name=name, cases=cases, samples=SAMPLES, sampling_ns=2.,
                       waveform_sha256=digest(paths[0]), truth_sha256=digest(paths[1]),
                       inference_observation_file=paths[0].name,
                       truth_not_in_observation_file=True))
    print('DATASET', name, len(cases), flush=True)


def make_noise(rng, sigma, rho=0.):
    noise = rng.normal(0., sigma, SAMPLES)
    if rho:
        factor = np.sqrt(1-rho*rho)
        for i in range(1, SAMPLES):
            noise[i] = rho*noise[i-1]+factor*noise[i]
    return noise


def mean_controls():
    from .alpha_pair_fitting import PairFitter
    check_freeze()
    p = protocol()
    fitter = PairFitter(settings())
    y, qs, ts, cases, components = [], [], [], [], []
    for ci, condition in enumerate(p['mean_conditions']):
        for gap in GAPS:
            for rep in range(p['repetitions_per_point']):
                # The same onset/noise within a replicate isolates the gap change.
                seed = p['mean_seed']+ci*10000+rep
                rng = np.random.default_rng(seed)
                first = rng.uniform(96., 256.)
                times = np.array([first, first+gap])
                charge = np.array([condition['q1'], condition['q2']])
                curves = (fitter.columns(GRID, times)*charge[None]).T
                y.append(fitter.baseline+curves.sum(0)+make_noise(rng, condition['noise_rms'], condition['ar1']))
                qs.append(charge); ts.append(times); components.append(curves)
                cases.append(dict(kind='pair', condition=condition['name'], gap_ns=gap,
                                  replicate=rep, cluster=f"mean-{ci}-{rep}", seed=seed,
                                  q1=condition['q1'], q2=condition['q2'],
                                  noise_rms=condition['noise_rms'], ar1=condition['ar1']))
        for count, number in ((0, p['mean_empty_controls_per_condition']),
                              (1, p['mean_single_controls_per_condition'])):
            for rep in range(number):
                seed = p['mean_seed']+ci*10000+1000+count*1000+rep
                rng = np.random.default_rng(seed); first = rng.uniform(96., 256.)
                charge = np.zeros(2); times = np.array([first, first]); curves = np.zeros((2,SAMPLES))
                if count:
                    charge[0] = (condition['q1'], condition['q2'], condition['q1']+condition['q2'])[rep%3]
                    curves[0] = fitter.columns(GRID, [first])[:,0]*charge[0]
                y.append(fitter.baseline+curves.sum(0)+make_noise(rng, condition['noise_rms'], condition['ar1']))
                qs.append(charge); ts.append(times); components.append(curves)
                cases.append(dict(kind='single' if count else 'empty', condition=condition['name'],
                                  gap_ns=None, replicate=rep, cluster=f"mean-simple-{ci}-{count}-{rep}",
                                  seed=seed, noise_rms=condition['noise_rms'], ar1=condition['ar1']))
    save_dataset('mean', y, qs, ts, cases, components)


def native(workers=4):
    import awkward as ak
    from .alpha_gamma_regularized_campaign import macro
    from .alpha_rate_data import export_window, open_root
    check_freeze()
    assert 1 <= workers <= 4
    manifest = load(EXPERIMENT/'manifest.json')
    completed = {r['label'] for r in manifest['runs']}
    assert all(r['status'] == 'completed' for r in manifest['runs'])
    plans = [p for p in manifest['planned_runs'] if p['label'] not in completed]
    cache = Path(tempfile.mkdtemp(prefix='gagg-resolution-native-'))
    binary = cache/'gagg-brems'
    shutil.copyfile(EXPERIMENT/'bin/gagg-brems', binary); binary.chmod(0o755)
    env = os.environ.copy(); env.update(manifest['environment']['g4_dataset_environment'])
    lock = threading.Lock(); last_launch = [0.]

    def one(plan):
        with lock:
            while True:
                wait = 1.05-(time.monotonic()-last_launch[0])
                if wait > 0: time.sleep(wait)
                root = PROJECT/'g4/data'/(datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss')+'.root')
                if not root.exists(): break
                time.sleep(1.05)
            row = dict(plan, root=str(root.relative_to(PROJECT)))
            path = EXPERIMENT/'macros'/(row['label']+'.mac')
            log = EXPERIMENT/'logs'/(row['label']+'.log')
            assert not path.exists() and not log.exists()
            body = macro(EXPERIMENT, row, root)
            path.write_text(body)
            stream = log.open('w')
            started = time.monotonic()
            proc = subprocess.Popen([str(binary), str(path)], cwd=PROJECT/'g4/build',
                                    env=env, stdout=stream, stderr=subprocess.STDOUT)
            last_launch[0] = time.monotonic()
            row.update(started=datetime.now().astimezone().isoformat(),
                       sipm_seed=row['seed']+100000, readout_seed=row['seed']+200000,
                       macro=str(path.relative_to(PROJECT)), macro_sha256=digest(path),
                       log=str(log.relative_to(PROJECT)))
        code = proc.wait(); stream.close()
        if code:
            save(EXPERIMENT/'failed-run.json', dict(**row, returncode=code))
            raise RuntimeError('G4 failed: '+str(log))
        local = cache/root.name
        shutil.copyfile(root, local)
        sha = digest(root); assert sha == digest(local)
        meta = export_window(local, OUTPUT/'native-inputs', allowed_particles=('alpha',))
        with open_root(local) as f:
            info = ak.to_list(f['RunInfo'].arrays(library='ak'))[0]
            assert f['WaveformEvent'].num_entries == f['SiPMEvent'].num_entries == 1
            assert info['random_seed'] == row['seed'] and info['macro_contents'] == body
            assert info['sipm_avalanche_base_seed'] == row['sipm_seed']
            assert info['sipm_readout_base_seed'] == row['readout_seed']
        assert meta['primaries'] == 1 and meta['root_sha256'] == sha
        meta.update(root=row['root'], io_note='Verified regular local copy; no links')
        save(OUTPUT/'native-inputs/metadata'/(root.stem+'.json'), meta)
        row.update(status='completed', returncode=0, elapsed_s=time.monotonic()-started,
                   root_sha256=sha, root_bytes=root.stat().st_size,
                   finished=datetime.now().astimezone().isoformat(),
                   validation=dict(runinfo_macro_verified=True, seed_verified=True,
                                   charge_balance_error_C=meta['charge_balance_error_C']))
        return row

    if not manifest['preflight_validated']:
        row = one(plans.pop(0))
        manifest['runs'].append(row); manifest['preflight_validated'] = True
        save(EXPERIMENT/'manifest.json', manifest)
        print('NATIVE_PREFLIGHT_PASSED', round(row['elapsed_s'],2), flush=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(one, row):row for row in plans}
        failures = []
        for future in as_completed(futures):
            try:
                row = future.result()
            except Exception as error:
                failed = dict(plan=futures[future], error=repr(error))
                failures.append(failed)
                manifest.setdefault('failed_attempts', []).append(failed)
                save(EXPERIMENT/'manifest.json', manifest)
                print('NATIVE_FAILED', failed['plan']['label'], repr(error), flush=True)
                continue
            manifest['runs'].append(row)
            save(EXPERIMENT/'manifest.json', manifest)
            print('NATIVE', len(manifest['runs']), '/', len(manifest['planned_runs']),
                  row['label'], round(row['elapsed_s'],1), flush=True)
    if failures:
        raise RuntimeError(f'{len(failures)} native runs failed; successful runs are recorded. Preserve attempts and retry the same seeds separately.')
    records = []
    for r in sorted(manifest['runs'], key=lambda r:r['seed']):
        stem = Path(r['root']).stem
        records.append(dict(root=r['root'], seed=r['seed'], energy_MeV=r['energy_MeV'],
                            replicate=r['replicate'], particle='alpha', root_original_sha256=r['root_sha256'],
                            source_manifest=str((EXPERIMENT/'manifest.json').relative_to(PROJECT)),
                            waveform=str((OUTPUT/'native-inputs/waveforms'/(stem+'.npz')).relative_to(PROJECT)),
                            truth=str((OUTPUT/'native-inputs/truth'/(stem+'.npz')).relative_to(PROJECT)),
                            metadata=str((OUTPUT/'native-inputs/metadata'/(stem+'.json')).relative_to(PROJECT))))
    save(OUTPUT/'native-records.json', dict(test_single=records))
    print('NATIVE_COMPLETE', len(records), flush=True)


def donor_cache():
    from .sipm_response import SensorConfig, read_photon_hits, map_cells
    records = load(OUTPUT/'native-records.json')['test_single']
    for i, r in enumerate(records):
        path = OUTPUT/'photons'/f'{i:04}.npz'
        if path.exists(): continue
        with tempfile.TemporaryDirectory(prefix='gagg-resolution-photons-') as temp:
            local = Path(temp)/Path(r['root']).name
            shutil.copyfile(PROJECT/r['root'], local)
            assert digest(local) == r['root_original_sha256']
            hits = read_photon_hits(local)
        config = load(PROJECT/r['metadata'])['config']
        cells = map_cells(hits, SensorConfig(**config['sensor']))
        first = float(hits.time_ns[cells>=0].min())
        np.savez_compressed(path, time_ns=hits.time_ns-first, x_mm=hits.x_mm,
                            y_mm=hits.y_mm, energy_eV=hits.energy_eV)
        if (i+1)%20 == 0: print('DONORS', i+1, '/', len(records), flush=True)
    save(OUTPUT/'photons/manifest.json', dict(records=records,
         files={p.name:digest(p) for p in (OUTPUT/'photons').glob('*.npz')},
         alignment='First active-plane photon; no amplitude rescaling or photon thinning.'))
    return records


_DONORS = None
_CONFIG = None


def _init_digitizer():
    global _DONORS, _CONFIG
    records = load(OUTPUT/'native-records.json')['test_single']
    _DONORS = []
    for i in range(len(records)):
        with np.load(OUTPUT/'photons'/f'{i:04}.npz') as f:
            _DONORS.append({k:f[k].copy() for k in f.files})
    _CONFIG = load(PROJECT/records[0]['metadata'])['config']


def joint_case(case):
    from .sipm_response import PhotonHits, SensorConfig, AvalancheResult, simulate_avalanches, map_cells
    from .sipm_readout import ReadoutConfig, synthesize_readout
    first = case['first_ns']; gap = case.get('gap_ns') or 0.
    seed = case['seed']; window = 4096.
    arrays = {k:[] for k in ('time_ns','x_mm','y_mm','energy_eV')}; owners = []
    for j, index in enumerate(case['parents']):
        donor = _DONORS[index]
        times = donor['time_ns']+first+j*gap
        keep = (times>=0.) & (times<window)
        for key in arrays:
            arrays[key].append((times if key=='time_ns' else donor[key])[keep])
        owners.extend([j]*int(keep.sum()))
    arrays = {k:np.concatenate(v) if v else np.array([]) for k,v in arrays.items()}
    owners = np.asarray(owners, int)
    order = np.argsort(arrays['time_ns'],kind='stable')
    arrays = {k:v[order] for k,v in arrays.items()}; owners = owners[order]
    hits = PhotonHits(**arrays, window_ns=window)
    sensor = SensorConfig(**_CONFIG['sensor']); readout = ReadoutConfig(**_CONFIG['readout'])
    avalanches = simulate_avalanches(hits, sensor, seed+100000)
    output = synthesize_readout(avalanches, readout, seed+200000)
    rec = avalanches.records; assignment = np.full(len(rec), -1, int)
    direct = rec['kind']==0
    assignment[direct] = owners[rec['photon_index'][direct]]
    for i in np.flatnonzero(rec['parent_id']>=0):
        assert rec['parent_id'][i] < i
        assignment[i] = assignment[rec['parent_id'][i]]
    q = np.zeros(2); times = np.full(2, first); curves = np.zeros((2,SAMPLES))
    replay = replace(readout, baseline_mV=0., noise_rms_mV=0.)
    total = np.zeros_like(output.arrays['signal_mV'])
    cells = map_cells(hits, sensor)
    for owner in range(-1, len(case['parents'])):
        block = rec[assignment==owner]
        signal = synthesize_readout(AvalancheResult(block,window,avalanches.full_charge_C,np.array([]),{}), replay,0).arrays['signal_mV']
        total += signal
        if owner >= 0:
            q[owner] = block['charge_C'].sum()/avalanches.full_charge_C
            times[owner] = float(hits.time_ns[(owners==owner)&(cells>=0)].min())
            curves[owner] = signal[:SAMPLES]
    error = float(np.max(np.abs(total-output.arrays['signal_mV'])))
    assert error < 1e-8
    return (output.arrays['analog_mV'][:SAMPLES], q, times, curves,
            dict(case, actual_gap_ns=float(times[1]-times[0]) if len(case['parents'])==2 else None,
                 component_sum_max_error_mV=error, noise_rms=.1, ar1=0.))


def physical_controls(workers=4):
    check_freeze()
    records = donor_cache()
    p = protocol(); groups = {}
    for energy, name in zip(p['native_energies_MeV'], ('low','mid','high')):
        groups[name] = [i for i,r in enumerate(records) if r['energy_MeV']==energy]
        assert len(groups[name]) == 80
    pairs = {}
    for name in ('low','mid','high'):
        pairs[name+'_'+name] = list(zip(groups[name][::2], groups[name][1::2]))
    pairs['low_high'] = list(zip(groups['low'][:40], groups['high'][:40]))
    pairs['high_low'] = [(b,a) for a,b in pairs['low_high']]
    cases = []
    for ci, (condition, members) in enumerate(pairs.items()):
        for gap in GAPS:
            for rep, parents in enumerate(members):
                seed = p['digitizer_seed_start']+ci*10000+rep
                first = float(np.random.default_rng(seed).uniform(96.,256.))
                cases.append(dict(kind='pair', condition=condition, gap_ns=gap, replicate=rep,
                                  parents=list(parents), seed=seed, first_ns=first,
                                  cluster='parents-'+'-'.join(map(str,sorted(parents)))))
    for rep in range(80):
        cases.append(dict(kind='empty', condition='all', gap_ns=None, replicate=rep,
                          parents=[], seed=p['digitizer_seed_start']+100000+rep,
                          first_ns=128., cluster=f'empty-{rep}'))
    rows = [None]*len(cases)
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_digitizer,
                             mp_context=mp.get_context('spawn')) as pool:
        future_ids = {pool.submit(joint_case,c):i for i,c in enumerate(cases)}
        done = 0
        for future in as_completed(future_ids):
            rows[future_ids[future]] = future.result(); done += 1
            if done%100==0: print('SHARED',done,'/',len(cases),flush=True)
    # Original single-event voltage; crop relative to the known source trigger,
    # never by centering its true first-photon time.
    for i, r in enumerate(records):
        meta = load(PROJECT/r['metadata'])
        with np.load(PROJECT/r['waveform']) as f:
            t = f['time_ns']; y = f['analog_mV']
        with np.load(PROJECT/r['truth']) as f:
            q = float(f['charge_pe'][0]); first = float(f['first_active_ns'][0])
        pretrigger = float(np.random.default_rng(262280000+i).uniform(96.,256.))
        start = int(np.floor((5000.-pretrigger-t[0])/2.))
        origin = float(t[start]-.25)
        condition = ('low','mid','high')[p['native_energies_MeV'].index(r['energy_MeV'])]
        info = dict(kind='single', condition=condition, gap_ns=None, replicate=r['replicate'],
                    parents=[i], seed=r['seed'], cluster=f'parent-{i}',
                    crop='Nominal source trigger plus independent pretrigger jitter',
                    crop_first_index=start, origin_ns=origin, noise_rms=.1, ar1=0.)
        rows.append((y[start:start+SAMPLES], np.array([q,0.]), np.array([first-origin,0.]),
                     np.zeros((2,SAMPLES)), info))
    # The single-cell truth file omits component arrays: these rows intentionally
    # retain raw native voltage; individual pair components are saved separately.
    save_dataset('shared', [r[0] for r in rows], [r[1] for r in rows], [r[2] for r in rows],
                 [r[4] for r in rows])
    pair_only = [r for r in rows if r[4]['kind']=='pair']
    np.savez_compressed(OUTPUT/'assets/shared.pair-components.npz',
                        components_mV=np.asarray([r[3] for r in pair_only],np.float32))
    # Same microcell realizations plus independent electronics noise; no change
    # in input photons, event truth, detector model or ADC acceptance.
    chosen = [r for r in rows if (r[4]['kind']=='pair' and r[4]['condition']=='mid_mid')
              or r[4]['kind']=='empty' or (r[4]['kind']=='single' and r[4]['condition']=='mid')]
    noisy = []
    for i,r in enumerate(chosen):
        seed = 262290000+i
        noise = np.random.default_rng(seed).normal(0.,np.sqrt(1.-.1**2),SAMPLES)
        noisy.append((r[0]+noise,r[1],r[2],dict(r[4],noise_rms=1.,added_noise_seed=seed)))
    save_dataset('shared-noise', [r[0] for r in noisy], [r[1] for r in noisy],
                 [r[2] for r in noisy], [r[3] for r in noisy])
    save(OUTPUT/'shared-provenance.json', dict(native_parents=len(records), cases=len(rows),
         independent_pairs_per_condition=40, workers=workers, shared_pixels=True,
         actual_Q_not_rescaled=True, photon_time_shift_only=True,
         seeds_reused_across_gaps='Paired/common-random-number design; cluster uncertainty by parents.',
         scope='Joint offline digitization of fresh native G4 optical transports, not new native pair ROOTs.'))
    print('PHYSICAL_COMPLETE',len(rows),flush=True)


def linear_controls():
    from .alpha_pair_data import native_pool
    records = load(OUTPUT/'native-records.json')['test_single']
    destination = OUTPUT/'native-pool.npz'
    if not destination.exists():
        native_pool(records, destination)
    with np.load(destination) as f:
        pool = {k:f[k].copy() for k in f.files}
    original = load(OUTPUT/'assets/shared.json')['cases']
    chosen = [c for c in original if (c['kind']=='pair' and c['condition']=='mid_mid')
              or c['kind']=='empty' or (c['kind']=='single' and c['condition']=='mid')]
    with np.load(OUTPUT/'assets/shared.waveforms.npz') as f:
        observed = f['voltage_mV']
    with np.load(OUTPUT/'assets/shared.truth.npz') as f:
        trueq = f['charge_pe']; truet = f['times_ns']
    y,qs,ts,cases,components=[],[],[],[],[]
    pool_grid=np.arange(pool['component'].shape[1])*2.+.25
    for source_index,c in enumerate(original):
        if c not in chosen: continue
        if c['kind']!='pair':
            y.append(observed[source_index]);qs.append(trueq[source_index]);ts.append(truet[source_index])
            components.append(np.zeros((2,SAMPLES)));cases.append(dict(c,linear_component_available=False))
            continue
        charge=np.zeros(2); times=np.array([c['first_ns'],c['first_ns']+c['gap_ns']])
        curves=np.zeros((2,SAMPLES))
        for j,index in enumerate(c['parents']):
            charge[j]=pool['charge_pe'][index]
            curves[j]=np.interp(GRID-times[j]+pool['times_ns'][index],pool_grid,
                                pool['component'][index],left=0.,right=0.)
        # One original pre-event background, not the sum of two noise traces.
        nuisance=pool['background'][c['parents'][0],:SAMPLES]
        y.append(10.+curves.sum(0)+nuisance);qs.append(charge);ts.append(times);components.append(curves)
        cases.append(dict(c,linear_component_available=True,
                          scope='Unscaled separately digitized native components; linear interpolation after electronics.'))
    save_dataset('linear', y,qs,ts,cases,components)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('prepare','mean','native','physical','linear'))
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.stage=='prepare': prepare()
    elif args.stage=='mean': mean_controls()
    elif args.stage=='native': native(args.workers)
    elif args.stage=='physical': physical_controls(args.workers)
    else: linear_controls()


if __name__=='__main__':
    main()
