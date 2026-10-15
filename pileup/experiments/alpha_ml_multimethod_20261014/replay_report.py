"""Run post-test reports from frozen source, regardless of later project edits."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

FOLDER = Path(__file__).resolve().parent
PROJECT = FOLDER.parents[2]
for name in ('analysis-source.json', 'post-test-diagnostics-source.json'):
    manifest = json.loads((FOLDER / name).read_text())
    for relative, expected in manifest['files'].items():
        with (FOLDER / relative).open('rb') as stream:
            assert hashlib.file_digest(stream, 'sha256').hexdigest() == expected, relative
sys.path.insert(0, str(FOLDER / 'source'))
sys.path.append(str(PROJECT))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('report', 'diagnostics'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--benchmark', choices=('new_test', 'regression'), default='new_test')
    arguments = parser.parse_args()
    if arguments.stage == 'report':
        from pileup.analysis.alpha_ml_multimethod_report import timing_null, figures
        timing_null(arguments.output / arguments.benchmark)
        figures(arguments.output, arguments.benchmark)
    else:
        from pileup.analysis.alpha_ml_multimethod_diagnostics import diagnostics
        diagnostics(arguments.output, FOLDER, arguments.benchmark)
