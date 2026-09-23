#!/usr/bin/env python3
"""Submit independent serial G4 runs: one waveform window per ROOT file."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime, timezone
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from threading import Lock
import time
import traceback
import uuid


def stamp():
    return datetime.now(timezone.utc).isoformat()


def write(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    temporary.replace(path)


def load_runner(case):
    spec = importlib.util.spec_from_file_location('frozen_runner', case / 'run.py')
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    runner.verify()
    build = runner.read(case / 'build/build_record.json')
    if build['bundle_sha256'] != runner.digest(case / 'bundle.json'):
        raise ValueError('Build belongs to another experiment snapshot.')
    if build['executable_sha256'] != runner.digest(case / 'build/gagg'):
        raise ValueError('Executable differs from its build record.')
    return runner, build


def prepare(args):
    project = Path(__file__).resolve().parents[1]
    case = (project / args.experiment).resolve()
    runner, build = load_runner(case)
    config = runner.read(case / 'config.json')
    if not 1 <= args.workers <= 8 or args.files_per_rate < 1:
        raise ValueError('Use 1..8 workers and at least one file per rate.')
    jobs = []
    if args.replay:
        original = runner.read(args.replay)
        if original['events'] != 1 or original['bundle_sha256'] != runner.digest(case / 'bundle.json'):
            raise ValueError('Replay needs a one-window manifest from the selected experiment snapshot.')
        point = next(point for point in config['points'] if point['id'] == original['point']['id'])
        jobs.append({'point': point, 'index': 0, 'events': 1, 'seed': original['seed'], 'status': 'queued'})
    else:
        for index in range(args.files_per_rate):
            for point in config['points']:
                seed = point['seed'] + args.seed_offset + 1000 * index
                if not 1 <= seed <= 2147483646:
                    raise ValueError('Seed outside the supported range.')
                jobs.append({'point': point, 'index': index, 'events': 1, 'seed': seed, 'status': 'queued'})
    if len({job['seed'] for job in jobs}) != len(jobs):
        raise ValueError('Each file needs a distinct seed.')
    name = 'alpha_one_window_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '_' + uuid.uuid4().hex[:6]
    directory = project / 'pileup/results/submissions' / name
    directory.mkdir(parents=True, exist_ok=False)
    shutil.copy2(__file__, directory / 'run_scan.py')
    record = {'submission_id': name, 'created_utc': stamp(), 'project_root': str(project),
              'experiment': str(case), 'output_directory': str(project / 'g4/data'),
              'files_per_rate': 1 if args.replay else args.files_per_rate, 'events_per_file': 1, 'max_workers': args.workers,
              'replay_of': str(args.replay.resolve()) if args.replay else None,
              'seed_offset': args.seed_offset, 'seed_rule': 'point seed + seed_offset + 1000 * file_index',
              'filename_format': 'gagg_waveform_YYYYMMDD_HHhMMmSSs.root', 'minimum_launch_spacing_seconds': 1.0,
              'bundle_sha256': runner.digest(case / 'bundle.json'), 'build_record': build,
              'driver_sha256': runner.digest(directory / 'run_scan.py'), 'config': config,
              'project_git_commit': runner.capture(['git', '-C', str(project), 'rev-parse', 'HEAD']),
              'project_worktree_status': runner.capture(['git', '-C', str(project), 'status', '--porcelain=v1']),
              'jobs': jobs}
    write(directory / 'submission.json', record)
    write(directory / 'status.json', {'status': 'prepared', 'submission_id': name, 'jobs': jobs})
    print(directory, flush=True)
    return directory


def execute(directory):
    submission = json.loads((directory / 'submission.json').read_text())
    case = Path(submission['experiment'])
    project = Path(submission['project_root'])
    output_directory = Path(submission['output_directory'])
    output_directory.mkdir(parents=True, exist_ok=True)
    for variable in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        os.environ[variable] = '1'
    runner, build = load_runner(case)
    if runner.digest(Path(__file__)) != submission['driver_sha256']:
        raise ValueError('Submission driver was modified.')
    if runner.digest(case / 'bundle.json') != submission['bundle_sha256']:
        raise ValueError('Experiment changed after submission preparation.')
    environment = runner.environment()
    lock = Lock()
    launch_lock = Lock()
    launch_state = {'last_time': 0.0, 'used_names': set()}
    jobs = submission['jobs']
    state = {'status': 'running', 'submission_id': submission['submission_id'], 'pid': os.getpid(),
             'started_utc': stamp(), 'max_workers': submission['max_workers'], 'events_per_file': 1,
             'files_per_rate': submission['files_per_rate'], 'jobs': jobs}
    write(directory / 'status.json', state)

    def update(job, **values):
        with lock:
            job.update(values)
            write(directory / 'status.json', state)

    def reserve_filename():
        # Serialize starts so second-resolution names stay short and unique.
        with launch_lock:
            delay = 1.0 - (time.monotonic() - launch_state['last_time'])
            if delay > 0:
                time.sleep(delay)
            while True:
                local_time = datetime.now().astimezone()
                path = output_directory / ('gagg_waveform_' + local_time.strftime('%Y%m%d_%Hh%Mm%Ss') + '.root')
                if path.name not in launch_state['used_names'] and not path.exists():
                    launch_state['used_names'].add(path.name)
                    launch_state['last_time'] = time.monotonic()
                    return path, local_time
                time.sleep(1.0)

    def run(job):
        if (directory / 'STOP_AFTER_CURRENT').exists():
            update(job, status='cancelled', reason='STOP_AFTER_CURRENT requested')
            return
        run_directory = directory / 'runs' / job['point']['id'] / f"{job['index']:04d}"
        run_directory.mkdir(parents=True, exist_ok=False)
        root_path, local_time = reserve_filename()
        if (directory / 'STOP_AFTER_CURRENT').exists():
            update(job, status='cancelled', reason='STOP_AFTER_CURRENT requested')
            return
        macro = (case / 'macros' / (job['point']['id'] + '.mac')).read_text().splitlines()
        replacements = {'/run/beamOn': '1', '/gagg/run/seed': str(job['seed']),
                        '/gagg/run/outputFile': str(root_path)}
        for command, value in replacements.items():
            positions = [i for i, line in enumerate(macro) if line.startswith(command + ' ')]
            if len(positions) != 1:
                raise ValueError(f'Expected exactly one {command} command.')
            macro[positions[0]] = command + ' ' + value
        macro_path = run_directory / 'run.mac'
        macro_path.write_text('\n'.join(macro) + '\n')
        command = [str(case / 'build/gagg'), str(macro_path)]
        record = {'status': 'running', 'started_utc': stamp(), 'experiment': submission['config']['experiment'],
                  'submission_id': submission['submission_id'], 'point': job['point'], 'events': 1,
                  'seed': job['seed'], 'config': submission['config'], 'command': command,
                  'output': str(root_path), 'root_relative_path': str(root_path.relative_to(project)),
                  'macro_relative_path': str(macro_path.relative_to(project)), 'filename_timezone': str(local_time.tzinfo),
                  'bundle_sha256': submission['bundle_sha256'], 'build_record': build,
                  'macro_sha256': runner.digest(macro_path), 'environment': environment,
                  'driver_sha256': submission['driver_sha256'], 'project_git_commit': submission['project_git_commit'],
                  'replay_of': submission.get('replay_of'),
                  'source_rate_definition': 'injected alpha particles per second'}
        manifest = run_directory / 'manifest.json'
        write(manifest, record)
        update(job, status='running', started_utc=record['started_utc'],
               root_file=record['root_relative_path'], manifest=str(manifest.relative_to(project)))
        started = time.monotonic()
        try:
            with (run_directory / 'simulation.log').open('w') as stream:
                process = subprocess.run(command, cwd=run_directory, stdout=stream, stderr=subprocess.STDOUT)
            record['exit_code'] = process.returncode
            if process.returncode or not root_path.is_file():
                raise RuntimeError(f'G4 failed: {run_directory / "simulation.log"}')
            update(job, status='validating')
            with (run_directory / 'validation.log').open('w') as stream:
                subprocess.run([sys.executable, str(case / 'analysis/check_run.py'), str(root_path),
                                '--output', str(run_directory / 'validation.json')],
                               stdout=stream, stderr=subprocess.STDOUT, check=True)
            report = json.loads((run_directory / 'validation.json').read_text())
            expected_gate = [submission['config']['pre_window_us'],
                             submission['config']['window_us'] - submission['config']['post_window_us']]
            if (report['n_windows'] != 1 or report['random_seed'] != job['seed'] or
                report['configured_injection_rate_hz'] != job['point']['rate_hz'] or
                report['window_us'] != submission['config']['window_us'] or report['analysis_gate_us'] != expected_gate):
                raise ValueError('ROOT settings do not match this single-window run manifest.')
            record.update(status='completed', root_sha256=runner.digest(root_path), root_bytes=root_path.stat().st_size)
            record['validation'] = report
            print(job['point']['id'], job['index'], root_path.name, 'validated', flush=True)
        except Exception as error:
            record.update(status='failed', error=f'{type(error).__name__}: {error}')
            traceback.print_exc()
        finally:
            record.update(finished_utc=stamp(), seconds=time.monotonic() - started)
            write(manifest, record)
            update(job, status=record['status'], finished_utc=record['finished_utc'],
                   seconds=record['seconds'], error=record.get('error'), validation=record.get('validation'))

    errors = []
    with ThreadPoolExecutor(max_workers=submission['max_workers']) as pool:
        futures = {pool.submit(run, job): job for job in jobs}
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as error:
                errors.append(str(error))
                update(futures[future], status='failed', error=str(error))
                traceback.print_exc()
    state['status'] = ('failed' if errors or any(job['status'] == 'failed' for job in jobs) else
                       'stopped' if any(job['status'] == 'cancelled' for job in jobs) else 'completed')
    state['finished_utc'] = stamp()
    state['errors'] = errors
    state['summary'] = []
    for point in submission['config']['points']:
        finished = [job for job in jobs if job['point']['id'] == point['id'] and job['status'] == 'completed']
        state['summary'].append({'point_id': point['id'], 'rate_hz': point['rate_hz'], 'files': len(finished),
                                 'windows': sum(job['validation']['n_windows'] for job in finished),
                                 'primary_particles': sum(job['validation']['primary_particles'] for job in finished)})
    write(directory / 'status.json', state)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment', default='pileup/experiments/alpha_p11b_v1')
    parser.add_argument('--files-per-rate', type=int, default=100)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--seed-offset', type=int, default=200)
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--replay', type=Path, help='Replay one saved single-window run with the same seed into a new ROOT file')
    parser.add_argument('--submission', type=Path, help='Execute an already prepared submission with its saved driver')
    args = parser.parse_args()
    directory = args.submission.resolve() if args.submission else prepare(args)
    if args.prepare_only:
        return
    with (directory / 'scan.log').open('a', buffering=1) as stream:
        with redirect_stdout(stream), redirect_stderr(stream):
            execute(directory)


if __name__ == '__main__':
    main()
