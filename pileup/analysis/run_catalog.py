"""Index central ROOT storage and exact saved macros without reading photon arrays."""
from pathlib import Path
import json
import uproot
from uproot.source.file import MemmapSource


def collect(project):
    project = Path(project).resolve()
    data = project / 'g4/data'
    records = {}
    paths = list((project / 'pileup/experiments').glob('*/runs/*/*/manifest.json'))
    paths += list((project / 'pileup/results/submissions').glob('*/runs/*/*/manifest.json'))
    for path in paths:
        record = json.loads(path.read_text())
        if record.get('status') == 'deleted':
            continue
        if 'root_relative_path' in record:
            root_path = project / record['root_relative_path']
        else:
            root_path = path.parent / record.get('output', 'waveform.root')
        root_path = root_path.resolve()
        if data not in root_path.parents:
            continue
        records[root_path] = (record, path)
    # Native single-particle campaigns keep multiple single-window runs in one
    # manifest. Use each run's saved macro, never today's same-named macro.
    for path in (project / 'jupyter/run_records').glob('sipm_native_*/manifest.json'):
        campaign = json.loads(path.read_text())
        for run in campaign.get('runs', []):
            if not run.get('root'):
                continue
            root_path = (project / run['root']).resolve()
            if data not in root_path.parents:
                continue
            record = {'status': 'completed' if run.get('returncode') == 0 else 'failed',
                      'events': 1, 'config': {'window_us': campaign['window_ns'] / 1000},
                      'seed': run['g4_seed'], 'macro_relative_path': run['macro'],
                      'campaign_id': path.parent.name, 'single_particle': True}
            records[root_path] = (record, path)
    # Native rate/ML campaigns save one manifest with explicit checks and runs.
    # Keep the saved macro and status, including a same-seed check's role.
    for path in (project / 'pileup/experiments').glob('*/manifest.json'):
        campaign = json.loads(path.read_text())
        for role in ('checks', 'runs'):
            for run in campaign.get(role, []):
                if not run.get('root'):
                    continue
                root_path = (project / run['root']).resolve()
                if data not in root_path.parents:
                    continue
                record = {'status': run.get('status', 'unknown'), 'events': 1,
                          'config': {'window_us': run.get('window_us', campaign.get('window_us'))},
                          'point': {'rate_hz': run.get('rate_hz')},
                          'seed': run.get('seed'), 'macro_relative_path': run.get('macro'),
                          'campaign_id': path.parent.name,
                          'native_waveform_label': run.get('label', ''),
                          'environment_check': role == 'checks'}
                records[root_path] = (record, path)
    files = {path.resolve() for path in data.rglob('*.root') if path.is_file()} | set(records)
    rows = []
    for path in sorted(files):
        record, manifest = records.get(path, ({}, None))
        macro = record.get('macro_relative_path')
        if not macro and manifest and (manifest.parent / 'run.mac').is_file():
            macro = str((manifest.parent / 'run.mac').relative_to(project))
        state = record.get('status', 'historical')
        row = {'root_file': str(path.relative_to(project)), 'status': state,
               'size_MiB': round(path.stat().st_size / 1024**2, 3) if path.is_file() else None,
               'waveform_events': record.get('events'), 'window_us': record.get('config', {}).get('window_us'),
               'alpha_injection_rate_Hz': record.get('point', {}).get('rate_hz'),
               'seed': record.get('seed'), 'macro_file': macro,
               'manifest': str(manifest.relative_to(project)) if manifest else None,
               'submission_id': record.get('submission_id'), 'macro_path_in_RunInfo': None,
               'campaign_id': record.get('campaign_id', record.get('submission_id')),
               'notes': ''}
        if state in ('running', 'validating'):
            row['notes'] = '运行或检查中；只读运行清单，不读取正在写入的 ROOT。'
        elif state in ('failed', 'interrupted'):
            row['notes'] = '失败或中断输出，未作为完整数据使用；保留供追溯。'
        else:
            try:
                with uproot.open(path, handler=MemmapSource) as root:
                    row['waveform_events'] = root['WaveformEvent'].num_entries if 'WaveformEvent' in root else None
                    if 'RunInfo' not in root or root['RunInfo'].num_entries != 1:
                        row['status'] = 'missing_RunInfo'
                        row['notes'] = 'ROOT 可打开，但没有唯一的 RunInfo 记录；运行条件及完整性尚未核实。'
                        rows.append(row)
                        continue
                    info = root['RunInfo']

                    def value(name):
                        if name not in info:
                            return None
                        result = info[name].array(library='np')[0]
                        return result.item() if hasattr(result, 'item') else result

                    row['waveform_events'] = root['WaveformEvent'].num_entries if 'WaveformEvent' in root else value('n_events')
                    length = value('t_length_ps')
                    row['window_us'] = length * 1e-6 if length is not None else None
                    row['seed'] = value('random_seed')
                    rate = value('alpha_p11b_rate_hz')
                    if rate is not None:
                        row['alpha_injection_rate_Hz'] = rate
                    row['macro_path_in_RunInfo'] = value('macro_file')
                    row['primary_particles'] = root['PrimaryPhoton'].num_entries if 'PrimaryPhoton' in root else None
                if record.get('single_particle'):
                    row['alpha_injection_rate_Hz'] = None
                    row['notes'] = '单粒子独立窗，无稳态注入率；实际宏、器件配置和种子见运行清单。'
                elif not macro:
                    row['notes'] = '历史文件：未找到完整原始宏快照；RunInfo 中的路径不等于当前同名宏。'
                elif row['waveform_events'] != 1:
                    row['notes'] = '调整前的多窗测试文件，保留为历史记录；当前生产要求一窗一文件。'
                else:
                    row['notes'] = '一窗一文件；宏和种子由本次实际运行清单记录。'
                    if record.get('native_waveform_label'):
                        row['notes'] += ' 标签：' + record['native_waveform_label'] + '。'
                    if record.get('environment_check'):
                        row['notes'] += ' 同种子环境对照，不计独立性能样本。'
            except Exception as error:
                row['status'] = 'unreadable'
                row['notes'] = f'{type(error).__name__}: {error}'
        rows.append(row)
    return rows
