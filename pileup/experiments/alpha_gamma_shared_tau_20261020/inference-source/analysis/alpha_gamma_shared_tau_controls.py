"""Matched Gaussian model and controlled-double diagnostics, not native G4."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path

import numpy as np

from .alpha_gamma_shared_tau import ProfileProblem, fit_joint_shared, fit_problem
from .alpha_gamma_shared_tau_study import METHODS, POSTS, load_results
from .alpha_gamma_short_window_study import save, digest


def synthetic_problem(post, q, fraction, rng, second=None):
    x = np.arange(-20., post+.1, 2.)
    readout = dict(rise_ns=2., decay_ns=20., lowpass_ns=2., aperture_ns=.5)
    p = ProfileProblem(x, np.zeros(len(x)), readout, 13.618501389, 10., .2, .1)
    y = 10.+p.columns((95., 351.), -6.3)@np.array([q*(1-fraction), q*fraction])
    if second is not None:
        gap, q2, fraction2 = second
        y += p.columns((95., 351.), -6.3+gap)@np.array([q2*(1-fraction2), q2*fraction2])
    return ProfileProblem(x, y+rng.normal(0., .1, len(x)), readout, 13.618501389, 10., .2, .1)


def synthetic_worker(task):
    i, particle, q, fraction, shared, inflation = task
    # One maximal waveform per independent seed; all prefixes slice this same
    # noise realization, rather than pretending each window is a new event.
    p = synthetic_problem(1000., q, fraction, np.random.default_rng(262021700+i))
    result = dict(seed=262021700+i, particle=particle, true_charge=q, prefixes=[])
    for post in POSTS:
        use = p.x <= post+1e-7
        prefix = ProfileProblem(p.x[use], p.y[use], p.readout, p.area_per_pe, p.baseline, .2, .1)
        fits = {}
        for method in ('fixed_tau', 'shared_tau', 'free_tau', 'constrained_tau'):
            initial = fits.get('free_tau') if method == 'constrained_tau' else fits.get('shared_tau')
            fit = fit_problem(prefix, method, shared, inflation, initial)
            fit['charge_error'] = fit['charge_pe']/q-1
            fit['tau_fast_error'] = fit['tau_fast_ns']/95.-1
            fit['tau_slow_error'] = fit['tau_slow_ns']/351.-1
            fits[method] = fit
        result['prefixes'].append(dict(post_ns=post, methods=fits))
    return result


def run(output, workers=4):
    output = Path(output).resolve()
    inflation = json.loads((output/'freeze.json').read_text())['prior_inflation']
    destination = output/'matched-model'
    destination.mkdir(exist_ok=True)
    if not (output/'matched-calibration.json').exists():
        problems = []
        index = 0
        for fraction in (.70, .21):
            for q in (500., 1500., 5000.):
                for _ in range(4):
                    problems.append(synthetic_problem(1000., q, fraction, np.random.default_rng(262021600+index)))
                    index += 1
        shared = fit_joint_shared(problems)
        shared['generator'] = 'Matched double exponentials plus white Gaussian electronics0.1mV; mean-area Q, not a realized avalanche process'
        save(output/'matched-calibration.json', shared)
    shared = json.loads((output/'matched-calibration.json').read_text())
    tasks = []
    for particle, fraction in (('alpha', .70), ('gamma', .21)):
        for q in (500., 1500., 5000.):
            for _ in range(8):
                index = len(tasks)
                tasks.append((index, particle, q, fraction, shared, inflation))
    missing = [t for t in tasks if not (destination/f'{t[0]:03}.json').exists()]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        pending = {pool.submit(synthetic_worker, t): t[0] for t in missing}
        for i, future in enumerate(as_completed(pending), 1):
            save(destination/f'{pending[future]:03}.json', future.result())
            if i % 8 == 0:
                print('MATCHED_MODEL', i, '/', len(missing), flush=True)
    rows = load_results(destination)
    groups = []
    for post in POSTS:
        for particle in ('alpha', 'gamma'):
            for method in METHODS[:-1]:
                selected = [r for r in rows if r['particle'] == particle]
                fits = [next(p for p in r['prefixes'] if p['post_ns'] == post)['methods'][method] for r in selected]
                errors = [f['charge_error'] for f in fits]
                groups.append(dict(post_ns=post, particle=particle, method=method, targets=len(fits),
                    charge_5pct_success=sum(abs(e) <= .05 for e in errors),
                    both_tau_10pct_success=sum(max(abs(f['tau_fast_error']), abs(f['tau_slow_error'])) <= .1 for f in fits),
                    charge_and_tau_success=sum(abs(f['charge_error']) <= .05 and max(abs(f['tau_fast_error']), abs(f['tau_slow_error'])) <= .1 for f in fits),
                    charge_bias=float(np.mean(errors)), charge_rmse=float(np.sqrt(np.mean(np.square(errors)))),
                    tau_bound_count=sum(f['tau_at_search_bound'] for f in fits)))
    save(output/'matched-model-summary.json', dict(independent_calibration=24, independent_test=48,
        groups=groups, known_class_input=False, class_accuracy_assessed=False,
        generator='Exactly matched double-exponential/readout mean plus0.1mV independent Gaussian noise, no photon/SiPM shot fluctuations',
        charge_scope='Known mean-response area parameter, different from native realized avalanche Q',
        inference='Unknown onset/amplitudes; lifetimes follow each method; observation origin is preset, not a measured trigger; no trigger-latency claim',
        scope='Conditional identifiability/optimizer control, not a native or hardware performance guarantee'))
    if not (output/'controlled-doubles.json').exists():
        doubles = []
        for particle, f1, f2 in (('alpha', .70, .21), ('gamma', .21, .70)):
            for gap in (50., 100., 200., 400.):
                for ratio in (.25, 1., 4.):
                    seed = 262021900+len(doubles)//2
                    maximum = synthetic_problem(300., 2000., f1, np.random.default_rng(seed),
                                                 second=(gap, 2000.*ratio, f2))
                    for post in (150., 300.):
                        use = maximum.x <= post+1e-7
                        p = ProfileProblem(maximum.x[use], maximum.y[use], maximum.readout, maximum.area_per_pe,
                                           maximum.baseline, .2, .1)
                        fits = {}
                        for method in METHODS[:-1]:
                            fit = fit_problem(p, method, shared, inflation, diagnostics=False)
                            fits[method] = dict(charge_pe=fit['charge_pe'], first_charge_relative_error=fit['charge_pe']/2000.-1,
                                                tau_fast_ns=fit['tau_fast_ns'], tau_slow_ns=fit['tau_slow_ns'])
                        doubles.append(dict(first_particle=particle, gap_ns=gap, charge_ratio=ratio,
                                             seed=seed, post_ns=post, fits=fits))
        save(output/'controlled-doubles.json', dict(independent_waveforms=24, rows=doubles,
            scope='Two deterministic mean pulses plus white Gaussian noise; SINGLE-pulse fitter scored for first Q only; no native microcell coupling or second-event recovery'))
    save(output/'controls-provenance.json', dict(source_sha256=digest(__file__),
        frozen_prior_inflation=inflation, sources='Independent synthetic seeds262021600–262021923'))


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--workers', default=4, type=int)
    a = p.parse_args(); run(a.output, a.workers)
