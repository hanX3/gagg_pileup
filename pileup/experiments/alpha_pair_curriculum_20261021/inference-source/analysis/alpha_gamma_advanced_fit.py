"""Self-consistent conditional fits of isolated native GAGG/SiPM prefixes.

The emission physics remains two-component. Independent calibration supplies
the complete class/Q-dependent mean. All inference inputs are voltage-only;
calibration truth is never a per-event reconstruction input.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import minimize

from .alpha_gamma_short_window import (
    PrefixCalibration, PrefixConfig, build_calibration, find_trigger, fit_prefix,
)

PROJECT = Path(__file__).resolve().parents[1]
METHODS = ('legacy_gls', 'expanded_gls', 'conditional_mle', 'pca8_mle',
           'pca16_mle', 'bootstrap_mle', 'pca8_bootstrap_mle')


@dataclass(frozen=True)
class AdvancedConfig:
    bootstrap_repetitions: int = 32
    bootstrap_seed: int = 262020501
    max_iterations: int = 100
    log_charge_step: float = .003
    timing_step_ns: float = .04
    baseline_prior_sigma_mV: float = .2


def interpolation(grid, points):
    """Linear interpolation indices/weights with zero outside the support."""
    grid, points = np.asarray(grid), np.asarray(points)
    coordinate = (points - grid[0]) / (grid[1] - grid[0])
    left = np.floor(coordinate).astype(int)
    fraction = coordinate - left
    valid = (coordinate >= 0) & (coordinate <= len(grid)-1)
    left = np.clip(left, 0, len(grid)-2)
    fraction = np.where(coordinate >= len(grid)-1, 1., fraction)
    return left, np.where(valid, 1-fraction, 0.), np.where(valid, fraction, 0.)


def project_covariance(covariance, grid, points):
    """P K P.T without constructing dense interpolation matrix P.

    The covariance is restricted to the observed prefix before factorization.
    This preserves PSD and allows continuous onset updates.
    """
    i, a, b = interpolation(grid, points)
    return (covariance[np.ix_(i, i)] * np.outer(a, a)
            + covariance[np.ix_(i, i+1)] * np.outer(a, b)
            + covariance[np.ix_(i+1, i)] * np.outer(b, a)
            + covariance[np.ix_(i+1, i+1)] * np.outer(b, b))


def interpolate_rows(values, grid, points):
    i, a, b = interpolation(grid, points)
    return values[..., i] * a + values[..., i+1] * b


def lowrank_covariance(covariance, rank):
    """PCA low rank plus residual diagonal, preserving marginal variances."""
    eigenvalues, vectors = np.linalg.eigh(covariance)
    selected = np.argsort(eigenvalues)[-min(rank, len(eigenvalues)):]
    modes = vectors[:, selected] * np.sqrt(np.maximum(eigenvalues[selected], 0.))
    diagonal = np.maximum(np.diag(covariance) - np.sum(modes*modes, axis=1), 0.)
    result = modes @ modes.T
    result.flat[::len(result)+1] += diagonal
    return result, modes, diagonal


def gaussian_profile(values, mean, covariance, baseline, baseline_sigma):
    """Gaussian -2logL plus a common independent soft baseline constraint.

    Profile B exactly; keep logdet C when C depends on fitted Q/time. The
    omitted normalizing constants are identical across hypotheses/methods.
    """
    factor = cholesky(covariance, lower=True, check_finite=False)
    whitened = solve_triangular(factor, np.c_[np.asarray(values)-baseline-mean,
                                               np.ones(len(values))],
                                lower=True, check_finite=False)
    z, one = whitened.T
    denominator = float(one@one + 1/baseline_sigma**2)
    offset = float(one@z / denominator)
    residual = z-offset*one
    quadratic = float(residual@residual + offset**2/baseline_sigma**2)
    logdet = float(2*np.sum(np.log(np.diag(factor))))
    return {'objective': quadratic+logdet, 'quadratic': quadratic,
            'logdet_covariance': logdet, 'baseline_mV': baseline+offset}


def calibration_samples(records, project=PROJECT):
    """Independent calibration labels/Q/onsets are allowed here only."""
    grid = np.arange(-200., 3000.1, 2.)
    rows = []
    for record in records:
        with np.load(project / record['waveform']) as wave, np.load(project / record['truth']) as truth:
            q = float(truth['charge_pe'][0]); onset = float(truth['first_active_ns'][0])
            if q < 250. or not np.isfinite(onset):
                continue
            t, y = wave['time_ns'], wave['analog_mV']
            baseline = float(np.median(y[(t >= 500.) & (t < 4000.)]))
            curve = np.interp(grid+onset, t, y-baseline)
            rows.append({'particle': str(truth['particle_name'][0]), 'q': q,
                         'curve': curve, 'root': record['root']})
    return rows


def bootstrap_means(normalized_curves, weights, repetitions, rng):
    """Stratified nonparametric bootstrap; rows are independent ROOT events.

    Bins are fixed at the original calibration partition. This estimates mean
    uncertainty conditional on the partition, not all calibration systematics.
    """
    n = len(normalized_curves)
    indices = rng.integers(0, n, size=(repetitions, n))
    means = np.array([np.average(normalized_curves[s], axis=0, weights=weights[s]) for s in indices])
    # Center before Gaussian marginal approximation; use the finite-n correction
    # that is exact for equal-weight sample means, approximate for weighted means.
    means -= means.mean(axis=0)
    return means * math.sqrt(n / (n-1))


def build_advanced_calibration(records, readout, output, config=AdvancedConfig(), project=PROJECT):
    output = Path(output)
    absolute_records = [{**r, **{key: str(project/r[key]) for key in ('waveform', 'truth')}} for r in records]
    base = build_calibration(absolute_records, readout, 13.618501389, output, PrefixConfig())
    samples = calibration_samples(records, project)
    rng = np.random.default_rng(config.bootstrap_seed)
    extra = {}
    for particle in ('alpha', 'gamma'):
        group = sorted([r for r in samples if r['particle'] == particle], key=lambda r: r['q'])
        bins = np.array_split(np.arange(len(group)), len(base.meta['classes'][particle]['counts']))
        deviations = []
        for indices in bins:
            q = np.array([group[i]['q'] for i in indices])
            curves = np.array([group[i]['curve']/group[i]['q'] for i in indices])
            deviations.append(bootstrap_means(curves, np.minimum(q, 4000.),
                                             config.bootstrap_repetitions, rng))
        extra[particle+'_bootstrap_deviations'] = np.asarray(deviations)
        extra[particle+'_bootstrap_deviations'][..., base.grid < 0] = 0.
    for rank in (8, 16):
        covariance, modes, diagonal = lowrank_covariance(base.arrays['shot_covariance_per_pe'], rank)
        extra[f'pca{rank}_shot_covariance'] = covariance
        extra[f'pca{rank}_modes'] = modes
        extra[f'pca{rank}_diagonal'] = diagonal
    np.savez_compressed(output / 'advanced-arrays.npz', **extra)
    (output / 'advanced-calibration.json').write_text(json.dumps({
        'advanced_config': asdict(config),
        'bootstrap_repetitions': config.bootstrap_repetitions,
        'bootstrap_seed': config.bootstrap_seed,
        'bootstrap_unit': 'One independent ROOT; resample separately within fixed class/Q bins',
        'uncertainty_model': 'Gaussian approximation to centered bootstrap mean error; Q² scaling, correlated time samples',
        'limits': ['Not a literal Poisson Barlow-Beeston waveform fit',
                   'Partition, covariance fitting and detector systematics are not fully resampled',
                   'Calibration mean uncertainty is shared by events, not new independent electronics noise'],
        'calibration_roots': [r['root'] for r in samples],
    }, ensure_ascii=False, indent=2)+'\n')
    return AdvancedCalibration(base.meta, {**base.arrays, **extra}, config)


class AdvancedCalibration(PrefixCalibration):
    def __init__(self, meta, arrays, advanced_config=AdvancedConfig()):
        super().__init__(meta, arrays)
        self.advanced_config = advanced_config
        self.maximum_charge_pe = max(c['charge_range_pe'][-1] for c in meta['classes'].values())*2
        shrink = self.config.covariance_shrinkage
        self.diagonal_shot = shrink*np.diag(arrays['shot_covariance_per_pe'])
        self.correlated_shot = arrays['shot_covariance_per_pe']-np.diag(self.diagonal_shot)

    @classmethod
    def load(cls, folder):
        folder = Path(folder)
        base = PrefixCalibration.load(folder)
        with np.load(folder / 'advanced-arrays.npz') as extra:
            arrays = {**base.arrays, **dict(extra)}
        cfg = json.loads((folder / 'advanced-calibration.json').read_text())
        return cls(base.meta, arrays, AdvancedConfig(**cfg['advanced_config']))

    def mean_deviations(self, points, particle, q):
        centers = np.asarray(self.meta['classes'][particle]['charge_centers_pe'])
        charge = float(np.clip(q, centers[0], centers[-1]))
        right = int(np.clip(np.searchsorted(centers, charge, side='right'), 1, len(centers)-1))
        f = math.log(charge/centers[right-1]) / math.log(centers[right]/centers[right-1])
        deviations = self.arrays[particle+'_bootstrap_deviations']
        curve = deviations[right-1]*(1-f) + deviations[right]*f
        return interpolate_rows(curve, self.grid, points)

    def conditional_covariance(self, points, q, particle, method):
        grid = self.arrays['covariance_grid_ns']
        rank = 8 if method.startswith('pca8') else 16 if method.startswith('pca16') else None
        if rank is not None:
            modes = interpolate_rows(self.arrays[f'pca{rank}_modes'].T, grid, points)
            covariance = modes.T@modes
            diagonal = self.arrays[f'pca{rank}_diagonal']
        else:
            covariance = project_covariance(self.correlated_shot, grid, points)
            diagonal = self.diagonal_shot
        # The shrinkage diagonal represents independent innovations at observed
        # samples. Interpolate its variance, not a synthetic white-noise trace;
        # P diag(D) P.T would create artificial half-sample variance minima.
        covariance.flat[::len(covariance)+1] += np.interp(points, grid, diagonal, left=0., right=0.)
        covariance *= q
        covariance.flat[::len(covariance)+1] += self.meta['electronics_sigma_mV']**2
        if 'bootstrap' in method:
            deviations = self.mean_deviations(points, particle, q)
            covariance += q*q * (deviations.T@deviations) / (len(deviations)-1)
        return covariance


def curvature_interval(objective, optimum, bounds, config):
    """Observed local curvature; empirical coverage must be measured separately."""
    steps = np.array([config.log_charge_step, config.timing_step_ns])
    if any(optimum[i]-steps[i] <= bounds[i][0] or optimum[i]+steps[i] >= bounds[i][1] for i in (0, 1)):
        return None, None
    f0 = objective(optimum)
    H = np.empty((2, 2))
    for i in (0, 1):
        delta = np.eye(2)[i]*steps[i]
        H[i, i] = (objective(optimum+delta)-2*f0+objective(optimum-delta))/steps[i]**2
    d0, d1 = np.diag(steps)
    H[0, 1] = H[1, 0] = (objective(optimum+d0+d1)-objective(optimum+d0-d1)
                           -objective(optimum-d0+d1)+objective(optimum-d0-d1))/(4*steps[0]*steps[1])
    if not np.isfinite(H).all() or np.linalg.eigvalsh(H)[0] <= 0:
        return None, None
    variance = float(2*np.linalg.inv(H)[0, 0])
    sigma_log = math.sqrt(variance)
    if not math.isfinite(sigma_log) or sigma_log > 5:
        return None, None
    q = math.exp(optimum[0])
    interval = [math.exp(optimum[0]-1.95996398454*sigma_log),
                math.exp(optimum[0]+1.95996398454*sigma_log)]
    return q*sigma_log, interval


def fit_hypothesis(x, y, calibration, particle, method, initial, pre_ns):
    cfg = calibration.advanced_config
    bounds = [(0., math.log(max(10000., calibration.maximum_charge_pe))), (-min(pre_ns, 60.), 2.)]
    evaluations = 0

    def evaluate(theta, details=False):
        nonlocal evaluations
        evaluations += 1
        q, onset = math.exp(float(theta[0])), float(theta[1])
        points = x-onset
        mean = q*calibration.template(points, particle, q)
        covariance = calibration.conditional_covariance(points, q, particle, method)
        result = gaussian_profile(y, mean, covariance, calibration.meta['baseline_mV'],
                                  cfg.baseline_prior_sigma_mV)
        if details:
            result['rms_mV'] = float(np.sqrt(np.mean((y-result['baseline_mV']-mean)**2)))
            return result
        return result['objective']

    start = np.array([math.log(max(1., initial['charge_pe'])), initial['onset_relative_ns']])
    start = np.clip(start, [b[0] for b in bounds], [b[1] for b in bounds])
    trials = []
    for dq in (-.35, 0., .35):
        for dt in (-1., 0., 1.):
            candidate = np.clip(start+[dq, dt], [b[0] for b in bounds], [b[1] for b in bounds])
            trials.append((evaluate(candidate), candidate))
    trials.sort(key=lambda value: value[0])
    optimized = []
    for _, candidate in trials[:2]:
        result = minimize(evaluate, candidate, bounds=bounds, method='L-BFGS-B',
                          options={'maxiter': cfg.max_iterations, 'ftol': 1e-10, 'gtol': 1e-5,
                                   'maxls': 25, 'finite_diff_rel_step': 1e-5})
        optimized.append(result)
    best = min(optimized, key=lambda result: float(result.fun))
    if not best.success:
        refined = minimize(evaluate, best.x, bounds=bounds, method='Powell',
                           options={'maxiter': 60, 'maxfev': 600, 'xtol': 1e-4, 'ftol': 1e-8})
        if refined.fun <= best.fun:
            best = refined
    optimum = np.asarray(best.x, float)
    value = evaluate(optimum, details=True)
    sigma, interval = curvature_interval(evaluate, optimum, bounds, cfg)
    value.update({
        'particle': particle, 'charge_pe': math.exp(optimum[0]),
        'onset_relative_ns': float(optimum[1]),
        'onset_at_search_bound': bool(min(optimum[1]-bounds[1][0], bounds[1][1]-optimum[1]) < .2),
        'charge_at_search_bound': bool(min(optimum[0]-bounds[0][0], bounds[0][1]-optimum[0]) < .002),
        'local_charge_sigma_pe': sigma, 'charge_interval95_pe': interval,
        'interval_scope': 'Local observed curvature, conditional on class and calibrated model; coverage evaluated empirically',
        'optimizer_converged': bool(best.success), 'optimizer_message': str(best.message),
        'optimizer_iterations': int(best.nit), 'objective_evaluations': evaluations,
        'model': 'Self-consistent class/Q mean, continuous Q/onset covariance, profiled constrained baseline',
    })
    return value


def with_legacy_interval(result):
    result = dict(result)
    q, sigma = result['charge_pe'], result['local_charge_sigma_pe']
    result['charge_interval95_pe'] = ([max(0., q-1.95996398454*sigma), q+1.95996398454*sigma]
                                     if sigma is not None and math.isfinite(sigma) else None)
    result['optimizer_converged'] = None
    result['interval_scope'] = 'Original local linearized curvature; not coverage-calibrated'
    return result


def fit_advanced_prefix(time_ns, voltage_mV, calibration, legacy_calibration=None,
                        *, pre_ns=20., post_ns=150., methods=METHODS, trigger=None):
    if any(method not in METHODS for method in methods):
        raise ValueError('Unknown advanced method')
    t, y = np.asarray(time_ns, float), np.asarray(voltage_mV, float)
    trigger = trigger or find_trigger(t, y, pre_ns, calibration.config.threshold_mV)
    if trigger is None:
        return {'detected': False, 'pre_ns': pre_ns, 'post_ns': post_ns, 'methods': {}}
    keep = (t >= trigger['time_ns']-pre_ns-1e-7) & (t <= trigger['time_ns']+post_ns+1e-7)
    x, values = t[keep]-trigger['time_ns'], y[keep]
    if len(x) < 8 or not np.isfinite(values).all():
        raise ValueError('Insufficient finite prefix samples')
    expanded = fit_prefix(t, y, calibration, pre_ns=pre_ns, post_ns=post_ns,
                          trigger=trigger, methods=('template_gls',))['methods']['template_gls']
    output = {}
    if 'legacy_gls' in methods:
        if legacy_calibration is None:
            raise ValueError('Original calibration is required for the legacy control')
        output['legacy_gls'] = with_legacy_interval(fit_prefix(
            t, y, legacy_calibration, pre_ns=pre_ns, post_ns=post_ns, trigger=trigger,
            methods=('template_gls',))['methods']['template_gls'])
    if 'expanded_gls' in methods:
        output['expanded_gls'] = with_legacy_interval(expanded)
    for method in methods:
        if method in ('legacy_gls', 'expanded_gls'):
            continue
        hypotheses = [fit_hypothesis(x, values, calibration, particle, method,
                                     next(r for r in expanded['hypotheses'] if r['particle'] == particle), pre_ns)
                      for particle in ('alpha', 'gamma')]
        chosen = dict(min(hypotheses, key=lambda r: r['objective']))
        chosen['class_score_gap'] = abs(hypotheses[0]['objective']-hypotheses[1]['objective'])
        chosen['hypotheses'] = hypotheses
        output[method] = chosen
    return {'detected': True, 'trigger': trigger, 'pre_ns': pre_ns, 'post_ns': post_ns,
            'actual_post_ns': float(x[-1]), 'sample_count': len(x),
            'last_observed_ns': float(t[keep][-1]), 'methods': output}
