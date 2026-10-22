"""Re-run saved frozen methods in a new directory without new G4 or training."""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys


experiment = Path(__file__).resolve().parent
project = experiment.parents[2]
os.environ['GAGG_ML_REFRESH_PROJECT'] = str(project)
sys.path.insert(0, str(experiment / 'inference-source'))
from pileup.analysis import alpha_ml_refresh as study
from pileup.analysis.alpha_rate_campaign import digest, save


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--phase', choices=['validation', 'common', 'confirmation'], default='confirmation')
    parser.add_argument('--one', help='Optional original ROOT filename stem')
    parser.add_argument('--workers', type=int, choices=[1, 2], default=1)
    args = parser.parse_args()
    study.check_frozen()
    source_output = study.OUTPUT
    original = study.load(source_output / args.phase / 'summary.json')
    rows = original['windows']
    if args.one:
        rows = [r for r in rows if Path(r['root']).stem == args.one]
        if len(rows) != 1:
            raise ValueError('Specify exactly one existing ROOT filename stem')
    output = (project / args.output).resolve() if not args.output.is_absolute() else args.output.resolve()
    if not output.is_relative_to(project):
        raise ValueError('Use a new output directory inside this project')
    output.mkdir(parents=True, exist_ok=False)
    for name in ('protocol.json', 'freeze.json'):
        shutil.copyfile(source_output / name, output / name)
    (output / 'pilots').mkdir()
    (output / args.phase).mkdir()
    save(output / 'replay-request.json', dict(
        source=str(source_output.relative_to(project)), phase=args.phase,
        roots=[r['root'] for r in rows], executable_sha256=digest(__file__),
        new_simulation=False, retraining=False, source_data_modified=False))
    study.OUTPUT = output
    study.records = lambda phase: rows
    study.evaluate(args.phase, args.workers)
    repeated = study.load(output / args.phase / 'summary.json')
    checks = []
    for row in repeated['windows']:
        old = next(r for r in original['windows'] if r['root'] == row['root'])
        for name, data in row['methods'].items():
            assert data['score'] == old['methods'][name]['score'], (row['root'], name)
            checks.append(dict(root=row['root'], method=name, original_score_exactly_reproduced=True))
    save(output / 'replay-validation.json', dict(status='passed', checks=checks,
         frozen_sources_only=True, new_simulation=False, retraining=False))
    print('REPLAY_PASSED', len(checks), flush=True)


if __name__ == '__main__':
    main()
