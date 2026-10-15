"""Analyze completed windows with the frozen baseline, optionally during acquisition.

Run from the project root with tools/run_gagg_python.sh and BLAS threads = 1.
The acquisition manifest is read-only here. Same-seed checks are excluded.
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
import json
import sys
import time

FOLDER = Path(__file__).resolve().parent
PROJECT = FOLDER.parents[2]
sys.path.insert(0, str(FOLDER / 'source'))
sys.path.insert(0, str(FOLDER / 'analysis_source'))
sys.path.append(str(PROJECT))

from pileup.analysis import alpha_rate_study as study
from pileup.analysis.alpha_rate_campaign import digest, save


def available_memory():
    return int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines()
                    if x.startswith('MemAvailable:'))) * 1024


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wait', action='store_true', help='Wait for already planned acquisition windows')
    parser.add_argument('--output', type=Path, default=Path('jupyter/artifacts/alpha_intermediate_rates_20261014_114000'))
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    initial = json.loads((FOLDER / 'manifest.json').read_text())
    for relative, expected in initial['source_sha256'].items():
        assert digest(FOLDER / relative) == expected, relative
    analysis_source = json.loads((FOLDER / 'analysis-source.json').read_text())
    for relative, expected in analysis_source['files'].items():
        assert digest(FOLDER / relative) == expected, relative
    assert Path(study.__file__).is_relative_to(FOLDER / 'source')
    study.CALIBRATION = FOLDER / 'source/jupyter/artifacts/alpha_gamma_fit_20261013_165844/summary.json'
    h, cfg = study.setup(output)
    baseline = json.loads((PROJECT / 'jupyter/artifacts/alpha_rate_fit_20261013_211541/analysis_config.json').read_text())
    current = json.loads((output / 'analysis_config.json').read_text())
    for key in ('template', 'fit_config', 'scoring', 'calibration_sha256'):
        assert current[key] == baseline[key], key
    rows = []
    for plan in initial['planned_runs']:
        while True:
            manifest = json.loads((FOLDER / 'manifest.json').read_text())
            run = next((r for r in manifest['runs'] if r['label'] == plan['label']), None)
            if run and run['status'] == 'completed':
                break
            if run and run['status'] in ('failed', 'interrupted'):
                raise RuntimeError(f"Acquisition {run['label']} {run['status']}")
            if not args.wait:
                raise RuntimeError('Acquisition incomplete; use --wait or finish acquisition first')
            time.sleep(5)
        fit_path = output / 'fits' / (Path(run['root']).stem + '.json')
        if fit_path.exists():
            row = json.loads(fit_path.read_text())
            assert row['root_sha256'] == run['root_sha256']
        else:
            while available_memory() < 4 * 1024**3:
                print('Waiting for 4 GiB available memory before export', flush=True)
                time.sleep(10)
            single = output / 'run_inputs' / run['label']
            single.mkdir(parents=True, exist_ok=True)
            save(single / 'manifest.json', dict(runs=[run], original_manifest=str(FOLDER / 'manifest.json')))
            start = time.monotonic()
            row = study.analyze(single, output, h, cfg)['windows'][0]
            row['export_and_analysis_seconds'] = time.monotonic() - start
            save(fit_path, row)
        rows.append(row)
        progress = dict(state='in_progress', completed=len(rows), planned=len(initial['planned_runs']),
                        rates=study.summarize(rows), windows=rows)
        save(output / 'progress_summary.json', progress)
        print('ANALYSIS', len(rows), '/', len(initial['planned_runs']), run['label'],
              'ADC', row['adc_fit']['status'], 'analog',
              None if row['analog_fit'] is None else row['analog_fit']['status'], flush=True)
    manifest = json.loads((FOLDER / 'manifest.json').read_text())
    assert len(manifest['runs']) == len(rows) == 12
    assert all(r['status'] == 'completed' for r in manifest['runs'])
    analog_rows = [dict(r, adc_fit=r['analog_fit'], adc_score=r['analog_score'], central_clipped_samples=0)
                   for r in rows if r['analog_fit'] is not None]
    summary = dict(state='completed', campaign=str(FOLDER.relative_to(PROJECT)),
        campaign_manifest_sha256=digest(FOLDER / 'manifest.json'),
        environment_checks_excluded=[r['label'] for r in manifest['checks']],
        rates=study.summarize(rows), windows=rows,
        analog_diagnostic_rates=study.summarize(analog_rows),
        analog_status_counts=dict(Counter(r['analog_fit']['status'] for r in rows if r['analog_fit'])),
        source_sha256={str((FOLDER / p).relative_to(PROJECT)): v for p, v in manifest['source_sha256'].items()},
        analysis_driver_sha256=digest(Path(__file__)),
        analysis_source=analysis_source,
        caveats=[
            'Three independent 20 us windows per source rate; central [5,15) us evaluation only',
            'Source injection rate differs from charge-producing event rate and q>=500 pe target rate',
            'Any ADC rail sample rejects the full record under the frozen baseline policy',
            'Analog diagnostic uses pre-ADC simulation voltage, not information recoverable from clipped ADC',
            'Known pure-alpha template family; blind pulse count, times and charges; no truth or rate input',
            'Realized shared-cell charge, not deposited-energy reconstruction; no real-device or online throughput validation',
            'Controls reuse the previously validated unchanged baseline; no threshold tuning on these windows'])
    save(output / 'summary.json', summary)
    save(output / 'progress_summary.json', dict(state='completed', completed=12, planned=12))
    print('ANALYSIS COMPLETE', output / 'summary.json', flush=True)


if __name__ == '__main__':
    main()
