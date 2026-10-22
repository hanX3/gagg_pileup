"""Variable-projection double-decay fits with independently shared lifetimes.

The four physical fits share their data, readout and pilot-frozen noise weights.
Only the lifetime constraint changes. No truth is accepted by inference.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import minimize, minimize_scalar

from .alpha_gamma_short_window import find_trigger, fit_prefix
from .alpha_gamma_waveform_fit import sampled_component

PHYSICAL_METHODS = ('fixed_tau', 'shared_tau', 'constrained_tau', 'free_tau')
TAU_BOUNDS = ((25., 250.), (250.5, 2000.))


@dataclass
class ProfileProblem:
    x: np.ndarray
    y: np.ndarray
    readout: dict
    area_per_pe: float
    baseline: float
    baseline_sigma: float
    factor: object

    def __post_init__(self):
        self.x, self.y = np.asarray(self.x, float), np.asarray(self.y, float)
        if not np.isfinite(self.y).all() or len(self.x) != len(self.y):
            raise ValueError('Nonfinite or inconsistent observed prefix')
        self.z = self.whiten(self.y-self.baseline)
        self.one = self.whiten(np.ones(len(self.x)))
        self.c = float(self.one@self.one + 1/self.baseline_sigma**2)
        self.zy = float(self.one@self.z)

    def whiten(self, value):
        if np.isscalar(self.factor):
            return np.asarray(value)/self.factor
        return solve_triangular(self.factor, value, lower=True, check_finite=False)

    def columns(self, taus, onset):
        return self.area_per_pe*np.column_stack([
            sampled_component(self.x, onset, tau, self.readout) for tau in taus])

    def profile(self, taus, onset, details=False):
        raw = self.columns(taus, onset)
        X = self.whiten(raw)
        cross = X.T@self.one
        gram = X.T@X-np.outer(cross, cross)/self.c
        rhs = X.T@self.z-cross*self.zy/self.c
        candidates = []
        for active in ((0, 1), (0,), (1,), ()):
            a = np.zeros(2)
            if active:
                ii = list(active)
                try:
                    a[ii] = np.linalg.solve(gram[np.ix_(ii, ii)], rhs[ii])
                except np.linalg.LinAlgError:
                    continue
                if np.any(a < 0):
                    continue
            b = (self.zy-cross@a)/self.c
            residual = self.z-X@a-b*self.one
            objective = float(residual@residual+b*b/self.baseline_sigma**2)
            candidates.append((objective, a, float(b), active))
        objective, a, b, active = min(candidates, key=lambda r: r[0])
        if not details:
            return objective
        q = float(np.sum(a))
        return dict(objective=objective, data_objective=objective,
                    charge_pe=q, amplitudes_pe=a.tolist(),
                    slow_fraction=float(a[1]/q) if q else 0.,
                    baseline_mV=self.baseline+b, onset_relative_ns=float(onset),
                    tau_fast_ns=float(taus[0]), tau_slow_ns=float(taus[1]),
                    active_amplitudes=list(active), gram_condition=float(np.linalg.cond(gram)),
                    rms_mV=float(np.sqrt(np.mean((self.y-self.baseline-b-raw@a)**2))))


def search_time(problem, taus, guess=None):
    grid = np.arange(-20., 2.01, 2.) if guess is None else np.unique(
        np.clip(np.r_[guess+np.array([-2., 0., 2.]), -20., 2.], -20., 2.))
    scores = [problem.profile(taus, t) for t in grid]
    i = int(np.argmin(scores))
    onset = float(grid[i])
    lo, hi = grid[max(0, i-1)], grid[min(len(grid)-1, i+1)]
    if hi > lo:
        result = minimize_scalar(lambda t: problem.profile(taus, t), bounds=(lo, hi),
                                 method='bounded', options={'xatol': .01})
        if result.fun < scores[i]:
            onset = float(result.x)
    return problem.profile(taus, onset, True)


def prior_precision(calibration, inflation):
    covariance = np.asarray(calibration['log_tau_covariance'], float)*float(inflation)**2
    return np.linalg.inv(covariance)


def fit_problem(problem, method, shared, inflation=1., initial=None, diagnostics=True):
    if method not in PHYSICAL_METHODS:
        raise ValueError(method)
    center = np.log(shared['taus_ns'])
    if method in ('fixed_tau', 'shared_tau'):
        taus = (95., 351.) if method == 'fixed_tau' else np.exp(center)
        result = search_time(problem, taus)
        result.update(optimizer_converged=True, optimizer_message='Bounded onset profile')
    else:
        precision = prior_precision(shared, inflation) if method == 'constrained_tau' else np.zeros((2, 2))
        bounds = [(math.log(a), math.log(b)) for a, b in TAU_BOUNDS]+[(-4., .4)]
        def objective(theta):
            delta = theta[:2]-center
            return problem.profile(np.exp(theta[:2]), 5.*theta[2])+float(delta@precision@delta)
        onset = initial['onset_relative_ns'] if initial else search_time(problem, np.exp(center))['onset_relative_ns']
        starts = [np.r_[center, onset/5.]]
        if initial:
            starts.append(np.r_[np.log([initial['tau_fast_ns'], initial['tau_slow_ns']]), onset/5.])
        pairs = ((50., 280.), (160., 600.), (60., 1000.), (180., 1200.))
        if method == 'free_tau':
            starts.extend(np.r_[np.log(pair), onset/5.] for pair in pairs)
        starts = [np.clip(s, [b[0] for b in bounds], [b[1] for b in bounds]) for s in starts]
        starts.sort(key=objective)
        trials = []
        for start in starts[:(3 if method == 'free_tau' else 2)]:
            opt = minimize(objective, start, method='L-BFGS-B', bounds=bounds,
                           options={'maxiter': 160, 'ftol': 1e-11, 'gtol': 1e-5, 'maxls': 30})
            trials.append(opt)
        best = min(trials, key=lambda r: r.fun)
        if not best.success:
            opt = minimize(objective, best.x, method='Powell', bounds=bounds,
                           options={'maxiter': 80, 'maxfev': 1200, 'ftol': 1e-9, 'xtol': 1e-5})
            if opt.fun < best.fun:
                best = opt
        result = problem.profile(np.exp(best.x[:2]), 5.*best.x[2], True)
        result.update(objective=float(best.fun), lifetime_penalty=float(best.fun-result['data_objective']),
                      optimizer_converged=bool(best.success), optimizer_message=str(best.message))
    result.update(method=method, prior_inflation=float(inflation) if method == 'constrained_tau' else None,
                  onset_at_search_bound=bool(min(result['onset_relative_ns']+20., 2.-result['onset_relative_ns']) < .2),
                  tau_at_search_bound=bool(method in ('constrained_tau', 'free_tau') and any(
                      min(abs(math.log(v/a)), abs(math.log(v/b))) < .005
                      for v, (a, b) in zip((result['tau_fast_ns'], result['tau_slow_ns']), TAU_BOUNDS))),
                  particle=None)
    if diagnostics:
        result.update(local_information(problem, result, shared, inflation))
    return result


def local_information(problem, result, shared, inflation):
    taus = np.array([result['tau_fast_ns'], result['tau_slow_ns']])
    onset, a = result['onset_relative_ns'], np.array(result['amplitudes_pe'])
    X = problem.columns(taus, onset)
    dt = (problem.columns(taus, onset+.01)-problem.columns(taus, onset-.01))@a/.02
    derivatives = []
    for j in range(2):
        plus, minus = taus.copy(), taus.copy()
        plus[j] *= math.exp(.0005); minus[j] *= math.exp(-.0005)
        derivatives.append((problem.columns(plus, onset)-problem.columns(minus, onset))@a/.001)
    # Calibration influence always profiles amplitudes, baseline and onset.
    Jlocal = problem.whiten(np.c_[X, np.ones(len(X)), dt])
    Jtau = problem.whiten(np.column_stack(derivatives))
    Hlocal = Jlocal.T@Jlocal
    Hlocal[2, 2] += 1/problem.baseline_sigma**2
    cross = Jlocal.T@Jtau
    inverse = np.linalg.pinv(Hlocal, rcond=1e-12)
    Htau = Jtau.T@Jtau-cross.T@inverse@cross
    residual = problem.whiten(problem.y-result['baseline_mV']-X@a)
    residual_aug_score = Jlocal.T@residual
    residual_aug_score[2] -= (result['baseline_mV']-problem.baseline)/problem.baseline_sigma**2
    score = Jtau.T@residual-cross.T@inverse@residual_aug_score
    output = dict(profile_log_tau_information=Htau.tolist(), profile_log_tau_score=score.tolist(),
                  charge_interval95_pe=None, local_log_tau_covariance=None)
    free = result['method'] in ('free_tau', 'constrained_tau')
    if free:
        info = np.block([[Hlocal, cross], [cross.T, Jtau.T@Jtau]])
        if result['method'] == 'constrained_tau':
            info[-2:, -2:] += prior_precision(shared, inflation)
    else:
        info = Hlocal
    eig = np.linalg.eigvalsh(info)
    valid = bool(eig[0] > 0 and len(result['active_amplitudes']) == 2
                 and not result.get('tau_at_search_bound') and not result.get('onset_at_search_bound'))
    if valid:
        cov = np.linalg.inv(info)
        direction = np.r_[1., 1., np.zeros(len(info)-2)]
        variance = float(direction@cov@direction)
        sigma = math.sqrt(max(0., variance))
        q = result['charge_pe']
        output['charge_interval95_pe'] = [max(0., q-1.96*sigma), q+1.96*sigma]
        if free:
            output['local_log_tau_covariance'] = cov[-2:, -2:].tolist()
    return output


def fit_joint_shared(problems, initial=(95., 351.), iterations=8):
    """Block coordinate variable projection: global taus, event-local onset/A/B.

    ROOT bootstrap below uses per-event profiled scores, never sample-wise
    resampling of correlated voltage points. The uncertainty is linearized.
    """
    center = np.log(initial)
    fits = [search_time(p, np.exp(center)) for p in problems]
    trace = []
    bounds = [(math.log(a), math.log(b)) for a, b in TAU_BOUNDS]
    for iteration in range(iterations):
        onsets = [r['onset_relative_ns'] for r in fits]
        def objective(logtau):
            taus = np.exp(logtau)
            return sum(p.profile(taus, t) for p, t in zip(problems, onsets))
        opt = minimize(objective, center, bounds=bounds, method='L-BFGS-B',
                       options={'maxiter': 120, 'ftol': 1e-12, 'gtol': 1e-5, 'maxls': 30})
        old = center.copy(); center = opt.x
        fits = [search_time(p, np.exp(center), r['onset_relative_ns']) for p, r in zip(problems, fits)]
        trace.append(dict(iteration=iteration, taus_ns=np.exp(center).tolist(),
                          objective=float(sum(r['objective'] for r in fits)),
                          max_log_tau_change=float(np.max(abs(center-old))), optimizer_converged=bool(opt.success)))
        print('SHARED_TAU', trace[-1], flush=True)
        if iteration >= 2 and np.max(abs(center-old)) < 2e-5:
            break
    mock = {'taus_ns': np.exp(center).tolist(), 'log_tau_covariance': np.eye(2).tolist()}
    for p, r in zip(problems, fits):
        r.update(method='shared_tau', onset_at_search_bound=False, tau_at_search_bound=False)
        r.update(local_information(p, r, mock, 1.))
    information = sum(np.array(r['profile_log_tau_information']) for r in fits)
    scores = np.array([r['profile_log_tau_score'] for r in fits])
    inverse = np.linalg.inv(information)
    scores -= scores.mean(axis=0)
    rng = np.random.default_rng(262021501)
    weights = rng.multinomial(len(fits), np.full(len(fits), 1/len(fits)), size=2000)
    shifts = (weights-1)@scores@inverse
    bootstrap_cov = np.cov(shifts, rowvar=False, ddof=1)
    # Explicit engineering floor avoids treating finite/model-dependent
    # calibration precision as exact microscopic lifetime knowledge.
    covariance = bootstrap_cov + np.eye(2)*.02**2
    return dict(taus_ns=np.exp(center).tolist(), log_tau_covariance=covariance.tolist(),
                log_tau_bootstrap_covariance=bootstrap_cov.tolist(),
                local_information_covariance=inverse.tolist(), iterations=trace,
                uncertainty='2000 linearized ROOT bootstrap replicates, plus 0.02 log-tau SD floor',
                calibration_events=len(problems), calibration_fits=fits)


def build_fraction_classifier(rows):
    """Independent calibration only; three Q bins per class and physical method."""
    classifier = {}
    for method in rows[0]['fits']:
        classifier[method] = {}
        for particle in ('alpha', 'gamma'):
            selected = sorted([r for r in rows if r['particle'] == particle], key=lambda r: r['q'])
            bins = np.array_split(np.arange(len(selected)), 3)
            knots = []
            for group in bins:
                qs = [selected[i]['q'] for i in group]
                fs = [selected[i]['fits'][method]['slow_fraction'] for i in group]
                knots.append(dict(log_charge=float(np.mean(np.log(qs))), mean=float(np.mean(fs)),
                                  sd=max(.05, float(np.std(fs, ddof=1))), count=len(group)))
            classifier[method][particle] = knots
    return classifier


def classify(fit, classifier, method_key=None):
    model = classifier[method_key or fit['method']]
    q = math.log(max(1., fit['charge_pe']))
    scores = {}
    for particle, knots in model.items():
        centers = [k['log_charge'] for k in knots]
        mean = float(np.interp(q, centers, [k['mean'] for k in knots]))
        sd = float(np.interp(q, centers, [k['sd'] for k in knots]))
        scores[particle] = ((fit['slow_fraction']-mean)/sd)**2+2*math.log(sd)
    return min(scores, key=scores.get), scores


def fit_shared_prefix(time_ns, voltage_mV, calibration, shared, *, post_ns,
                      inflation=1., classifier=None, trigger=None, diagnostics=True,
                      candidate_inflations=None):
    t, y = np.asarray(time_ns), np.asarray(voltage_mV)
    trigger = trigger or find_trigger(t, y, 20., calibration.config.threshold_mV)
    if trigger is None:
        return dict(detected=False, pre_ns=20., post_ns=post_ns, methods={})
    base = fit_prefix(t, y, calibration, pre_ns=20., post_ns=post_ns, trigger=trigger,
                      methods=('physics_ols', 'template_gls'))
    mask = (t >= trigger['time_ns']-20.-1e-7) & (t <= trigger['time_ns']+post_ns+1e-7)
    x, values = t[mask]-trigger['time_ns'], y[mask]
    pilot = base['methods']['physics_ols']
    C = calibration.covariance(x, pilot['onset_relative_ns'], max(1., pilot['charge_pe']))
    problem = ProfileProblem(x, values, calibration.readout, calibration.meta['area_per_pe_mV_ns'],
                             calibration.meta['baseline_mV'], .2, cholesky(C, lower=True, check_finite=False))
    output = {'template_gls': base['methods']['template_gls']}
    for method in ('fixed_tau', 'shared_tau', 'free_tau', 'constrained_tau'):
        initial = output.get('free_tau') if method == 'constrained_tau' else output.get('shared_tau')
        fit = fit_problem(problem, method, shared, inflation, initial, diagnostics)
        if classifier:
            key = f'constrained_tau_{inflation:g}' if method == 'constrained_tau' else method
            fit['particle'], fit['class_scores'] = classify(fit, classifier, key)
        output[method] = fit
    for value in candidate_inflations or ():
        key = f'constrained_tau_{value:g}'
        if value == inflation:
            output[key] = dict(output['constrained_tau'])
            continue
        fit = fit_problem(problem, 'constrained_tau', shared, value, output['free_tau'], diagnostics)
        if classifier:
            fit['particle'], fit['class_scores'] = classify(fit, classifier, key)
        output[key] = fit
    return dict(detected=True, trigger=trigger, pre_ns=20., post_ns=float(post_ns),
                last_observed_ns=float(t[mask][-1]), sample_count=len(x), methods=output,
                fit_truth_inputs=[], covariance_pilot=pilot,
                hypothesis='At most one isolated pulse, unknown class/Q/onset; count assumed')
