#!/usr/bin/env python3
"""Copied into each frozen experiment as run.py; needs no Git checkout."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time
import uuid

BASE = Path(__file__).resolve().parent


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def now():
    return datetime.now(timezone.utc).isoformat()


def capture(args):
    try:
        result = subprocess.run(args, capture_output=True, text=True, timeout=20)
        return (result.stdout + result.stderr).strip()
    except (OSError, subprocess.TimeoutExpired) as error:
        return str(error)


def environment():
    keys = ['G4NEUTRONHPDATA', 'G4LEDATA', 'G4LEVELGAMMADATA', 'G4RADIOACTIVEDATA',
            'G4PARTICLEXSDATA', 'G4PIIDATA', 'G4REALSURFACEDATA', 'G4SAIDXSDATA',
            'G4ABLADATA', 'G4INCLDATA', 'G4ENSDFSTATEDATA', 'G4CHANNELINGDATA',
            'CMAKE_PREFIX_PATH', 'LD_LIBRARY_PATH']
    return {'recorded_utc': now(), 'platform': platform.platform(),
            'python': sys.version, 'cmake': capture(['cmake', '--version']),
            'compiler': capture(['c++', '--version']),
            'geant4': capture(['geant4-config', '--version']),
            'root': capture(['root-config', '--version']),
            'dependency_environment': {key: os.environ[key] for key in keys if key in os.environ}}


def verify():
    manifest = read(BASE / 'bundle.json')
    for name, expected in manifest['files'].items():
        path = BASE / name
        if not path.is_file() or path.is_symlink() or digest(path) != expected:
            raise RuntimeError(f'冻结文件缺失或已修改：{name}；请从当前工程生成新的实验目录。')
    # CMake uses globbing: an added source file must not silently enter an old build.
    for directory in ['g4', 'macros', 'analysis']:
        actual = {str(path.relative_to(BASE)) for path in (BASE / directory).rglob('*')
                  if path.is_file() and '__pycache__' not in path.relative_to(BASE).parts}
        expected = {name for name in manifest['files'] if name.startswith(directory + '/')}
        if actual != expected:
            raise RuntimeError(f'冻结目录文件集合已改变：{directory}')
    return manifest, read(BASE / 'config.json')


def build(args):
    verify()
    directory = BASE / 'build'
    directory.mkdir(exist_ok=True)
    command = ['cmake', '-S', str(BASE / 'g4'), '-B', str(directory), '-DCMAKE_BUILD_TYPE=Release']
    for tool, variable, suffix in [('geant4-config', 'Geant4_DIR', 'lib/cmake/Geant4'),
                                   ('root-config', 'ROOT_DIR', 'cmake')]:
        if shutil.which(tool):
            candidate = Path(capture([tool, '--prefix'])) / suffix
            if candidate.is_dir():
                command.append(f'-D{variable}={candidate}')
    command.extend(args.cmake_arg)
    started = time.monotonic()
    with (directory / 'build.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run(['cmake', '--build', str(directory), '-j', str(args.jobs)],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    record = {'status': 'built', 'built_utc': now(), 'seconds': time.monotonic() - started,
              'bundle_sha256': digest(BASE / 'bundle.json'),
              'executable_sha256': digest(directory / 'gagg'), 'configure_command': command,
              'environment': environment(), 'linked_libraries': capture(['ldd', str(directory / 'gagg')])}
    write(directory / 'build_record.json', record)
    print('构建完成：', directory / 'gagg', flush=True)


def run_point(point_id, events=None, seed=None, replay_of=None):
    _, config = verify()
    points = {point['id']: point for point in config['points']}
    if point_id not in points:
        raise ValueError(f'未知条件 {point_id}，可用条件：{list(points)}')
    point = points[point_id]
    events = config['events'] if events is None else events
    seed = point['seed'] if seed is None else seed
    if not 1 <= events <= 2147483647 or not 1 <= seed <= 2147483646:
        raise ValueError('events 必须是正整数；seed 必须在 1..2147483646。')
    executable = BASE / 'build/gagg'
    if not executable.is_file() or not (BASE / 'build/build_record.json').is_file():
        raise RuntimeError('请先运行 python3 run.py build。')
    build_record = read(BASE / 'build/build_record.json')
    if build_record['bundle_sha256'] != digest(BASE / 'bundle.json') or build_record['executable_sha256'] != digest(executable):
        raise RuntimeError('可执行文件与冻结目录的构建记录不匹配，请重新 build。')
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '_' + uuid.uuid4().hex[:8]
    directory = BASE / 'runs' / point_id / run_id
    directory.mkdir(parents=True, exist_ok=False)
    macro = (BASE / 'macros' / f'{point_id}.mac').read_text()
    macro, event_count = re.subn(r'^/run/beamOn\s+\d+\s*$', f'/run/beamOn {events}', macro, flags=re.M)
    macro, seed_count = re.subn(r'^/gagg/run/seed\s+\d+\s*$', f'/gagg/run/seed {seed}', macro, flags=re.M)
    if event_count != 1 or seed_count != 1:
        raise RuntimeError('冻结宏必须恰好包含一次 beamOn 和一次 seed。')
    macro_path = directory / 'run.mac'
    macro_path.write_text(macro, encoding='utf-8')
    command = [str(executable), str(macro_path)]
    record = {'status': 'running', 'started_utc': now(), 'experiment': config['experiment'],
              'point': point, 'events': events, 'seed': seed, 'config': config,
              'bundle_sha256': digest(BASE / 'bundle.json'), 'build_record': build_record,
              'macro_sha256': digest(macro_path), 'command': command, 'environment': environment(),
              'replay_of': replay_of, 'output': 'waveform.root', 'source_rate_definition': 'injected alpha particles per second'}
    write(directory / 'manifest.json', record)
    started = time.monotonic()
    print('运行目录：', directory, flush=True)
    try:
        with (directory / 'simulation.log').open('w') as log:
            result = subprocess.run(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT)
        record['exit_code'] = result.returncode
        if result.returncode != 0 or not (directory / 'waveform.root').is_file():
            raise RuntimeError(f'仿真失败，请查看 {directory / "simulation.log"}')
        record['root_sha256'] = digest(directory / 'waveform.root')
        record['root_bytes'] = (directory / 'waveform.root').stat().st_size
        record['status'] = 'completed'
    except BaseException as error:
        record['status'] = 'interrupted' if isinstance(error, KeyboardInterrupt) else 'failed'
        record['error'] = str(error)
        raise
    finally:
        record['finished_utc'] = now()
        record['seconds'] = time.monotonic() - started
        write(directory / 'manifest.json', record)
    print('运行完成：', directory / 'waveform.root', flush=True)
    return directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    sub.add_parser('list', help='显示条件与预计注入粒子数，不运行仿真')
    sub.add_parser('verify', help='校验冻结文件 SHA256')
    builder = sub.add_parser('build', help='仅编译本实验保存的 g4 源码')
    builder.add_argument('--jobs', type=int, default=2)
    builder.add_argument('--cmake-arg', action='append', default=[])
    runner = sub.add_parser('run', help='在新的输出目录运行；不覆盖旧结果')
    choice = runner.add_mutually_exclusive_group(required=True)
    choice.add_argument('--point')
    choice.add_argument('--all', action='store_true')
    runner.add_argument('--events', type=int)
    runner.add_argument('--seed', type=int)
    replay = sub.add_parser('replay', help='按历史运行的宏参数与种子重跑，产生新目录')
    replay.add_argument('manifest', type=Path)
    args = parser.parse_args()
    _, config = verify()
    if args.action == 'verify':
        print('冻结文件校验通过。')
    elif args.action == 'list':
        gate = config['window_us'] - config['pre_window_us'] - config['post_window_us']
        print(f"{config['experiment']}: {config['events']} independent windows/point; full={config['window_us']} us; central={gate} us")
        for point in config['points']:
            print(f"{point['id']}: rate={point['rate_hz']:g} Hz, seed={point['seed']}, "
                  f"mean/window={point['rate_hz']*config['window_us']*1e-6:g}, "
                  f"mean/central={point['rate_hz']*gate*1e-6:g}")
    elif args.action == 'build':
        if args.jobs < 1:
            raise ValueError('--jobs 必须大于 0。')
        build(args)
    elif args.action == 'run':
        selected = [point['id'] for point in config['points']] if args.all else [args.point]
        if args.all and args.seed is not None:
            raise ValueError('--all 使用各条件原有的不同种子；覆盖种子时请使用 --point。')
        for point_id in selected:
            run_point(point_id, args.events, args.seed)
    else:
        record = read(args.manifest)
        if record['bundle_sha256'] != digest(BASE / 'bundle.json'):
            raise ValueError('历史运行属于另一个源码快照，请使用对应实验目录的 run.py。')
        run_point(record['point']['id'], record['events'], record['seed'], str(args.manifest.resolve()))


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(error, file=sys.stderr)
        sys.exit(1)
