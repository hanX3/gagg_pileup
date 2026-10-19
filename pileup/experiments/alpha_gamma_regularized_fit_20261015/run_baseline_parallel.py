"""Run the frozen baseline on disjoint held-out windows, one process per window.

This launcher changes scheduling only. Model, score and checked NNLS solver are
imported from the pre-test source snapshot. Run from the project root.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import json
from pathlib import Path
import sys
import time

FOLDER = Path(__file__).resolve().parent
sys.path.insert(0, str(FOLDER / 'source'))
sys.path.insert(1, str(FOLDER.parents[2]))


def run_one(row, inputs, output):
    import numpy as np
    from pileup.analysis import alpha_gamma_pileup_fit as engine
    from pileup.analysis.alpha_gamma_nnls import reduced_nnls
    from pileup.analysis.alpha_gamma_study import calibration
    from pileup.analysis.alpha_rate_campaign import digest, save

    inputs, output = Path(inputs), Path(output)
    source = FOLDER / 'source/pileup/analysis'
    frozen = json.loads((FOLDER / 'frozen-fit-settings.json').read_text())
    for name in ('alpha_gamma_pileup_fit.py', 'alpha_gamma_nnls.py', 'alpha_gamma_study.py'):
        path = source / name
        assert digest(path) == frozen['source_sha256'][str(path.relative_to(FOLDER))]
    stem = Path(row['root']).stem
    destination = output / 'baseline' / f'{stem}.json'
    wavepath = inputs / 'waveforms' / f'{stem}.npz'
    if destination.exists():
        stored = json.loads(destination.read_text())
        assert stored['waveform_sha256'] == digest(wavepath)
        assert stored['baseline_fitter_sha256'] == digest(source / 'alpha_gamma_pileup_fit.py')
        assert stored['nnls_helper_sha256'] == digest(source / 'alpha_gamma_nnls.py')
        return row['label'], 'already completed'
    engine.nnls = reduced_nnls
    _, bank, config = calibration(FOLDER / 'baseline-settings.json')
    wave = dict(np.load(wavepath))
    started = time.perf_counter()
    fit = engine.fit_waveform(wave['time_ns'], wave['analog_mV'], bank, config)
    seconds = time.perf_counter() - started
    # No truth is loaded until blind reconstruction is finished.
    truthpath = inputs / 'truth' / f'{stem}.npz'
    truth = dict(np.load(truthpath))
    meta = json.loads((inputs / 'metadata' / f'{stem}.json').read_text())
    score = engine.score_mixed(fit, truth, meta['gate_ns'])
    model = fit.pop('model_mV')
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(destination.with_suffix('.npz'), time_ns=wave['time_ns'], model_mV=model)
    save(destination, dict(row=row, method='baseline', fit=fit, score=score,
        fit_seconds=seconds, fitter_sha256=frozen['fitter_sha256'],
        baseline_fitter_sha256=digest(source / 'alpha_gamma_pileup_fit.py'),
        baseline_linear_solver='Equivalent Gram/Cholesky or QR NNLS reduction with original-space KKT checks',
        nnls_helper_sha256=digest(source / 'alpha_gamma_nnls.py'),
        calibration_sha256=digest(FOLDER / 'calibration.json'),
        waveform_sha256=digest(wavepath), truth_sha256=digest(truthpath),
        fit_truth_inputs=[], gate_ns=meta['gate_ns'],
        launcher_sha256=digest(Path(__file__)), scheduling='Independent processes; BLAS/OMP threads=1'))
    return row['label'], dict(seconds=round(seconds, 2),
        recovered={p: (v['recovered_count'], v['target_count']) for p, v in score['by_particle'].items()},
        extra=score['unmatched_predictions'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4, choices=range(1, 5))
    args = parser.parse_args()
    rows = [r for r in json.loads((FOLDER / 'manifest.json').read_text())['runs'] if r['role'] == 'test_mixed']
    assert len(rows) == 8 and all(r['status'] == 'completed' for r in rows)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs = [pool.submit(run_one, r, str(args.inputs), str(args.output)) for r in rows]
        for job in as_completed(jobs):
            print('BASELINE', *job.result(), flush=True)


if __name__ == '__main__':
    main()
