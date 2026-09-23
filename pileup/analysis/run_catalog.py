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
                if not macro:
                    row['notes'] = '历史文件：未找到完整原始宏快照；RunInfo 中的路径不等于当前同名宏。'
                elif row['waveform_events'] != 1:
                    row['notes'] = '调整前的多窗测试文件，保留为历史记录；当前生产要求一窗一文件。'
                else:
                    row['notes'] = '一窗一文件；宏和种子由本次实际运行清单记录。'
            except Exception as error:
                row['status'] = 'unreadable'
                row['notes'] = f'{type(error).__name__}: {error}'
        rows.append(row)
    return rows
