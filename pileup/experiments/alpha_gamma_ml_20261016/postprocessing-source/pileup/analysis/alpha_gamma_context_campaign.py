"""Native position and controlled-history calibration; one ROOT per window.

The schedule interface is implemented only in this experiment's frozen G4
source. The transport and shared-cell digitizer physics are unchanged.
"""
import argparse
from pathlib import Path
from . import alpha_gamma_regularized_campaign as base
from .alpha_rate_campaign import digest, save

PROJECT = Path(__file__).resolve().parents[2]
_ordinary_macro = base.macro


def macro(folder, row, output):
    text = _ordinary_macro(folder, row, output)
    if 'particle' in row:
        pos = row.get('position_mm', [0., 0., 0.])
        text = text.replace('/gun/position 0 0 0 mm',
                            '/gun/position ' + ' '.join(map(str, pos)) + ' mm')
    if 'schedule' in row:
        path = folder/'schedules'/(row['label']+'.json')
        assert not path.exists()
        save(path, dict(primaries=row['schedule']))
        row['schedule_file'] = str(path.relative_to(PROJECT))
        row['schedule_sha256'] = digest(path)
        text = text.replace('/run/beamOn 1',
                            f'/gagg/source/scheduleFile {path}\n/run/beamOn 1')
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--stage', required=True)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    base.PROJECT = PROJECT
    base.macro = macro
    base.run(args.folder, args.stage, args.inputs, args.workers)


if __name__ == '__main__':
    main()
