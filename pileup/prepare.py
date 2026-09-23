#!/usr/bin/env python3
"""Create a new self-contained alpha pile-up experiment without running Geant4."""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from tools.run_experiment import digest, environment, write

ROOT = Path(__file__).resolve().parents[1]


def validate(config):
    if config['schema_version'] != 1 or config['source_model'] != 'p11b_effective_single_alpha':
        raise ValueError('This preparation step supports the existing p–11B alpha spectrum only.')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', config['experiment']):
        raise ValueError('Experiment name must contain only letters, digits, underscore and hyphen.')
    for name in ['window_us', 'pre_window_us', 'post_window_us', 'disk_radius_mm',
                 'gagg_max_step_um', 'mylar_front_max_step_um', 'mylar_side_max_step_um']:
        if not math.isfinite(config[name]) or config[name] < 0:
            raise ValueError(f'Invalid finite nonnegative value for {name}.')
    if config['window_us']*1e6 > 2**32-1 or config['window_us'] <= config['pre_window_us'] + config['post_window_us']:
        raise ValueError('Window must fit uint32 ps and contain a nonempty central gate.')
    window, pre, post = (math.floor(config[key]*1e6 + 0.5)
                         for key in ['window_us', 'pre_window_us', 'post_window_us'])
    if window <= pre + post:
        raise ValueError('Central gate must remain nonempty after rounding to picoseconds.')
    for name in ['gagg_max_step_um', 'mylar_front_max_step_um', 'mylar_side_max_step_um']:
        if config[name] <= 0:
            raise ValueError(f'{name} must be positive.')
    if not isinstance(config['events'], int) or isinstance(config['events'], bool) or not 1 <= config['events'] < 2**31:
        raise ValueError('events must be a positive 32-bit integer.')
    seen_ids, seen_seeds = set(), set()
    if not config['points']:
        raise ValueError('No rate points.')
    for point in config['points']:
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', point['id']) or point['id'] in seen_ids:
            raise ValueError('Point IDs must be valid and unique.')
        if not math.isfinite(point['rate_hz']) or point['rate_hz'] < 0:
            raise ValueError('rate_hz must be finite and nonnegative.')
        if not isinstance(point['seed'], int) or isinstance(point['seed'], bool) or not 1 <= point['seed'] <= 2147483646 or point['seed'] in seen_seeds:
            raise ValueError('Each point requires a distinct seed in 1..2147483646.')
        seen_ids.add(point['id'])
        seen_seeds.add(point['seed'])


def macro(config, point):
    return f'''# Frozen experiment: {config['experiment']}; condition: {point['id']}
# Rate is the injected alpha-particle rate, NOT reaction or detected-event rate.
# Primaries: independent Poisson(rate*T) and uniform emission times in [0,T).
# Central analysis gate: [{config['pre_window_us']:g}, {config['window_us']-config['post_window_us']:g}) us.
# Pre/post guards are an initial choice; reconstruction convergence is not yet tested.
/control/verbose 1
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/gagg/step/setGaggMaxStep {config['gagg_max_step_um']:g} um
/gagg/step/setMylarFrontMaxStep {config['mylar_front_max_step_um']:g} um
/gagg/step/setMylarSideMaxStep {config['mylar_side_max_step_um']:g} um
/run/initialize

/gagg/run/seed {point['seed']}
/gagg/run/windowLength {config['window_us']:g} us
/gagg/run/preWindow {config['pre_window_us']:g} us
/gagg/run/postWindow {config['post_window_us']:g} us
/gagg/run/outputFile waveform.root

/gagg/source/useDiskConeSource true
/gagg/source/aimAtGagg true
/gagg/source/diskRadius {config['disk_radius_mm']:g} mm
/gagg/pileup/enable true
/gagg/source/brems/enable false
/gagg/source/brems/rateHz 0
/gagg/source/c12Capture/enable false
/gagg/source/c12Capture/rateHz 0
/gagg/source/alphaP11B/enable true
/gagg/source/alphaP11B/rateHz {point['rate_hz']:.12g}

/run/printProgress 1
/run/beamOn {config['events']}
'''


def prepare(config, target):
    validate(config)
    if target.exists():
        raise FileExistsError(f'实验目录已存在，不覆盖：{target}；请使用新的名称。')
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.' + target.name + '_preparing_', dir=target.parent))
    simulation = staging/'g4'
    simulation.mkdir()
    for name in ['CMakeLists.txt', 'main.cc']:
        shutil.copy2(ROOT/'g4'/name, simulation/name)
    for name in ['include', 'src', 'macros']:
        shutil.copytree(ROOT/'g4'/name, simulation/name,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc', 'AGENTS.md', '*.ipynb'))
    (staging/'macros').mkdir()
    for point in config['points']:
        (staging/'macros'/f"{point['id']}.mac").write_text(macro(config, point), encoding='utf-8')
    (staging/'analysis').mkdir()
    shutil.copy2(ROOT/'pileup/tools/run_experiment.py', staging/'run.py')
    shutil.copy2(ROOT/'pileup/analysis/check_run.py', staging/'analysis/check_run.py')
    write(staging/'config.json', config)
    write(staging/'environment_prepare.json', environment())
    source_diff = subprocess.check_output(['git', 'diff', '--binary', 'HEAD', '--', 'g4'], cwd=ROOT)
    (staging/'source_changes.patch').write_bytes(source_diff)
    (staging/'README.md').write_text(f'''# {config['experiment']}

本目录是可独立编译和重跑的实验快照。`g4/`、宏、配置、运行器和检查脚本均有实际副本，不依赖项目根目录的后续修改，也不要求切换 Git 提交。复用旧结果时保留整个实验目录；Geant4、ROOT 与编译环境仍须可用，版本和依赖路径见环境及构建记录。

在本目录执行：

```bash
python3 run.py list
python3 run.py verify
python3 run.py build --jobs 2
python3 run.py run --point {config['points'][0]['id']} --events 2
# 确认成本后，运行完整一档：
python3 run.py run --point {config['points'][0]['id']}
# 或运行全部条件：
python3 run.py run --all
# 在新的输出目录重跑某次历史运行：
python3 run.py replay runs/<条件>/<运行目录>/manifest.json
python3 analysis/check_run.py runs/<条件>/<运行目录>/waveform.root
```

每次运行生成新的 `runs/<条件>/<UTC时间_随机后缀>/`，保存 ROOT、实际 `run.mac`、标准输出日志及 `manifest.json`。`--events` 与 `--seed` 的覆盖值会记录在宏和清单中。重放使用历史运行的实际事件数和种子；ROOT 文件本身有时间戳／UUID，比较重跑结果应比较物理数据分支，而不是要求文件字节一致。

默认每档 {config['events']} 个独立窗，完整窗 {config['window_us']:g} μs，中央区 [{config['pre_window_us']:g}, {config['window_us']-config['post_window_us']:g}) μs。只启用当前 p–¹¹B 的有效单 α 源；源注入率不等于融合反应率或形成可测信号的事件率。前后保护区不是已验证收敛的最终选择。

本阶段仅提供仿真数据入口和一致性检查；未知脉冲数的 α 堆积重建及精度评价将在阅读用户提供的文献结果后实现。生产扫描由用户提交。

冻结文件改变后 `verify` 和运行器会拒绝继续复用该快照。需要改模型／宏参数时，从主项目重新生成另一个实验目录。
''', encoding='utf-8')
    files = {str(path.relative_to(staging)): digest(path) for path in sorted(staging.rglob('*')) if path.is_file()}
    manifest = {'schema_version': 1, 'experiment': config['experiment'],
                'created_utc': datetime.now(timezone.utc).isoformat(),
                'source_git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                'source_worktree_status': subprocess.check_output(
                    ['git', '--no-optional-locks', 'status', '--porcelain=v1', '--', 'g4'], cwd=ROOT, text=True).splitlines(),
                'files': files, 'reconstruction_status': 'not_implemented_pending_user_literature'}
    write(staging/'bundle.json', manifest)
    if target.exists():
        raise FileExistsError(f'目标目录在准备期间出现：{target}')
    staging.rename(target)
    print('已准备实验（尚未运行）：', target)
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT/'pileup/configs/alpha_p11b_v1.json')
    parser.add_argument('--name', help='New immutable experiment name')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding='utf-8'))
    if args.name:
        config['experiment'] = args.name
    target = args.output.resolve() if args.output else ROOT/'pileup/experiments'/config['experiment']
    prepare(config, target)


if __name__ == '__main__':
    main()
