"""Causal isolated-pulse prefix fits, with separately trained response statistics.

Inference accepts voltage, sample times and a frozen calibration only. The
single-pulse hypothesis is explicit: this is not a blind multi-pulse counter.
Charge means realized photon-descendant avalanche charge in full-cell pe.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from pathlib import Path
import json

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import minimize_scalar
from numpy.lib.stride_tricks import sliding_window_view

from .alpha_gamma_waveform_fit import sampled_component


@dataclass(frozen=True)
class PrefixConfig:
    threshold_mV: float = 2.
    consecutive_samples: int = 3
    baseline_prior_sigma_mV: float = .2
    covariance_shrinkage: float = .5
    onset_step_ns: float = 2.
    onset_tolerance_ns: float = .12
    onset_lookback_ns: float = 60.
    minimum_calibration_charge_pe: float = 250.


METHODS = ('physics_ols', 'physics_gls', 'template_ols', 'template_gls')


def find_trigger(time_ns, voltage_mV, pre_ns, threshold_mV=2.):
    """First sustained crossing of a past-only rolling-median baseline.

    The decision time is the THIRD above-threshold sample. Its entire lookback,
    including both confirmation samples, fits inside pre_ns. Vectorization is
    equivalent to streaming; no full-record peak or future-derived threshold.
    """
    t, y = np.asarray(time_ns), np.asarray(voltage_mV)
    dt = float(t[1]-t[0])
    span = int(round(pre_ns/dt))
    width = span-2
    if width < 4 or len(y) < span+1:
        return None
    medians = np.median(sliding_window_view(y, width), axis=1)
    above = y[width:] - medians[:-1] > threshold_mV
    indices = np.flatnonzero(above[:-2] & above[1:-1] & above[2:]) + width+2
    indices = indices[indices >= span]
    if not len(indices):
        return None
    index = int(indices[0])
    return {'index': index, 'time_ns': float(t[index]),
            'threshold_mV_above_rolling_median': threshold_mV,
            'lookback_ns': float(span*dt), 'confirmation_latency_ns': 2*dt}


def linear_profile(columns, values, baseline, baseline_sigma, factor=None):
    """Exact 1/2 nonnegative amplitudes with an unconstrained soft baseline.

    The baseline prior is learned from independent calibration, not reused
    prefix samples. Eliminate the baseline with a Schur complement and enumerate
    amplitude active sets. Coefficients use pe units for both model families.
    """
    X, z = np.asarray(columns, float), np.asarray(values, float)-baseline
    one = np.ones(len(z))
    if factor is not None and np.isscalar(factor):
        X, z, one = X/factor, z/factor, one/factor
    elif factor is not None:
        X = solve_triangular(factor, X, lower=True, check_finite=False)
        z = solve_triangular(factor, z, lower=True, check_finite=False)
        one = solve_triangular(factor, one, lower=True, check_finite=False)
    c = float(one@one + 1/baseline_sigma**2)
    cross, zy = X.T@one, float(one@z)
    gram = X.T@X - np.outer(cross, cross)/c
    rhs = X.T@z - cross*zy/c
    candidates = []
    for mask in range(1 << X.shape[1]):
        active = np.flatnonzero([(mask >> i) & 1 for i in range(X.shape[1])])
        a = np.zeros(X.shape[1])
        if len(active):
            g = gram[np.ix_(active, active)]
            try:
                a[active] = np.linalg.solve(g, rhs[active])
            except np.linalg.LinAlgError:
                a[active] = np.linalg.lstsq(g, rhs[active], rcond=1e-12)[0]
            if np.any(a[active] < 0):
                continue
        b = (zy-cross@a)/c
        r = z-X@a-b*one
        objective = float(r@r + b*b/baseline_sigma**2)
        candidates.append((objective, a, float(b), active))
    objective, a, b, active = min(candidates, key=lambda x: x[0])
    variance = None
    if len(active):
        try:
            variance = float(np.ones(len(active)) @ np.linalg.inv(gram[np.ix_(active, active)]) @ np.ones(len(active)))
        except np.linalg.LinAlgError:
            pass
    return {'objective': objective, 'amplitudes_pe': a, 'baseline_mV': baseline+b,
            'active_amplitudes': active.tolist(), 'fixed_time_charge_variance': variance,
            'gram_condition': float(np.linalg.cond(gram))}


class PrefixCalibration:
    def __init__(self, meta, arrays):
        self.meta, self.arrays = meta, arrays
        self.config = PrefixConfig(**meta['prefix_config'])
        self.readout = meta['readout']
        self.grid = arrays['grid_ns']

    @classmethod
    def load(cls, folder):
        folder = Path(folder)
        return cls(json.loads((folder/'calibration.json').read_text()), dict(np.load(folder/'calibration-arrays.npz')))

    def template(self, x, particle, shape_q):
        centers = np.asarray(self.meta['classes'][particle]['charge_centers_pe'])
        means = self.arrays[particle+'_means']
        q = float(np.clip(shape_q, centers[0], centers[-1]))
        right = int(np.searchsorted(centers, q, side='right'))
        if right == 0:
            curve = means[0]
        elif right == len(centers):
            curve = means[-1]
        else:
            f = math.log(q/centers[right-1])/math.log(centers[right]/centers[right-1])
            curve = means[right-1]*(1-f) + means[right]*f
        return np.interp(x, self.grid, curve, left=0., right=0.)

    def covariance(self, x, onset_guess, charge_guess):
        # Restrict the covariance first, then factor. Never crop its inverse.
        grid = self.arrays['covariance_grid_ns']
        index = np.rint((x-onset_guess-grid[0])/(grid[1]-grid[0])).astype(int)
        index = np.clip(index, 0, len(grid)-1)
        cov = self.arrays['shot_covariance_per_pe'][np.ix_(index, index)] * max(0., charge_guess)
        cov = cov.copy()
        cov.flat[::len(cov)+1] += self.meta['electronics_sigma_mV']**2
        return cov


def _search(x, y, columns_at, baseline, sigma_b, pre_ns, cfg, factor=None):
    def solve(onset):
        columns = columns_at(onset)
        return linear_profile(columns, y, baseline, sigma_b, factor)
    grid = np.arange(-min(pre_ns,cfg.onset_lookback_ns), 2.+cfg.onset_step_ns/2, cfg.onset_step_ns)
    trials = [solve(t0) for t0 in grid]
    best_index = int(np.argmin([r['objective'] for r in trials]))
    low, high = max(0, best_index-1), min(len(grid)-1, best_index+1)
    best, t0 = trials[best_index], float(grid[best_index])
    if high > low:
        opt = minimize_scalar(lambda v: solve(v)['objective'], bounds=(grid[low], grid[high]),
                              method='bounded', options={'xatol': cfg.onset_tolerance_ns})
        candidate = solve(float(opt.x))
        if candidate['objective'] < best['objective']:
            best, t0 = candidate, float(opt.x)
    best['onset_relative_ns'] = t0
    best['onset_at_search_bound'] = bool(t0 <= grid[0]+.2 or t0 >= grid[-1]-.2)
    best['charge_pe'] = float(sum(best['amplitudes_pe']))
    # Include timing as a nuisance in the local curvature diagnostic. This is
    # not a calibrated confidence interval; development data set its multiplier.
    column = columns_at(t0)
    derivative = ((columns_at(t0+.1)-columns_at(t0-.1)) @ best['amplitudes_pe'])/.2
    J = np.c_[column, np.ones(len(x)), derivative]
    if factor is not None and np.isscalar(factor):
        J = J/factor
    elif factor is not None:
        J = solve_triangular(factor, J, lower=True, check_finite=False)
    prior = np.zeros(J.shape[1]); prior[column.shape[1]] = 1/sigma_b
    info = J.T@J + np.outer(prior, prior)
    contrast = np.r_[np.ones(column.shape[1]), 0., 0.]
    try:
        variance = float(contrast @ np.linalg.solve(info, contrast))
        best['local_charge_sigma_pe'] = float(math.sqrt(max(0., variance)))
    except np.linalg.LinAlgError:
        best['local_charge_sigma_pe'] = None
    prediction = best['baseline_mV'] + column @ best['amplitudes_pe']
    best['rms_mV'] = float(np.sqrt(np.mean((y-prediction)**2)))
    best['amplitudes_pe'] = best['amplitudes_pe'].tolist()
    return best


def fit_prefix(time_ns, voltage_mV, calibration, *, pre_ns=50., post_ns=100., trigger=None,
               methods=METHODS):
    """Fit one detected pulse without class/Q/time truth or future voltages."""
    if any(m not in METHODS for m in methods):
        raise ValueError('Unknown prefix method')
    t, y = np.asarray(time_ns, float), np.asarray(voltage_mV, float)
    cfg = calibration.config
    trigger = trigger or find_trigger(t, y, pre_ns, cfg.threshold_mV)
    if trigger is None:
        return {'detected': False, 'pre_ns': pre_ns, 'post_ns': post_ns, 'methods': {}}
    mask = (t >= trigger['time_ns']-pre_ns-1e-7) & (t <= trigger['time_ns']+post_ns+1e-7)
    x, values = t[mask]-trigger['time_ns'], y[mask]
    if len(x) < 8 or not np.isfinite(values).all():
        raise ValueError('Insufficient finite prefix samples')
    baseline = calibration.meta['baseline_mV']
    area = calibration.meta['area_per_pe_mV_ns']
    physics = lambda t0: area*np.column_stack([
        sampled_component(x, t0, tau, calibration.readout) for tau in (95., 351.)])
    white_sigma = calibration.meta['electronics_sigma_mV']
    base = _search(x, values, physics, baseline, cfg.baseline_prior_sigma_mV, pre_ns, cfg, white_sigma)
    # A single voltage-derived charge/time pilot fixes SAME weights for both
    # particle hypotheses. No parameter-dependent likelihood is claimed.
    q0 = max(base['charge_pe'], 1.)
    factor = None
    if any(m.endswith('_gls') for m in methods):
        cov = calibration.covariance(x, base['onset_relative_ns'], q0)
        factor = cholesky(cov, lower=True, check_finite=False)
    output = {}
    for method in methods:
        L = factor if method.endswith('_gls') else white_sigma
        if method == 'physics_ols':
            result = dict(base)
        elif method == 'physics_gls':
            result = _search(x, values, physics, baseline, cfg.baseline_prior_sigma_mV, pre_ns, cfg, L)
        else:
            hypotheses = []
            for particle in ('alpha', 'gamma'):
                columns = lambda t0, name=particle: calibration.template(x-t0, name, q0)[:,None]
                trial = _search(x, values, columns, baseline, cfg.baseline_prior_sigma_mV, pre_ns, cfg, L)
                trial['particle'] = particle
                hypotheses.append(trial)
            result = dict(min(hypotheses, key=lambda r:r['objective']))
            result['hypotheses'] = hypotheses
            result['class_score_gap'] = abs(hypotheses[0]['objective']-hypotheses[1]['objective'])
        if method.startswith('physics'):
            result['slow_fraction'] = result['amplitudes_pe'][1]/result['charge_pe'] if result['charge_pe'] else None
            # Two freely fitted coefficients do not provide an independently
            # validated classifier. Keep this result separate from template PSD.
            result['particle'] = None
        result['shape_charge_pilot_pe'] = q0
        result['weights'] = ('frozen shared empirical covariance' if method.endswith('_gls')
                             else 'unweighted voltage SSE / electronics variance')
        output[method] = result
    return {'detected':True, 'trigger':trigger, 'pre_ns':pre_ns, 'post_ns':post_ns,
            'actual_post_ns':float(x[-1]),
            'total_prefix_span_ns':float(x[-1]-x[0]), 'sample_count':len(x),
            'last_observed_ns':float(t[mask][-1]), 'methods':output,
            'fit_truth_inputs':[], 'hypothesis':'At most one isolated pulse; count is not fitted'}


def build_calibration(records, readout, area_per_pe, output, cfg=PrefixConfig()):
    """Training-only labels/Q/onsets construct native response means and noise."""
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    grid = np.arange(-200., 3000.1, 2.)
    cov_grid = grid[grid <= 1250.]
    samples, baseline_values, zero_sigma = [], [], []
    for record in records:
        with np.load(record['waveform']) as wave, np.load(record['truth']) as truth:
            t, y = wave['time_ns'], wave['analog_mV']
            quiet = y[(t >= 500.) & (t < 4000.)]
            baseline_values.append(float(np.median(quiet)))
            zero_sigma.append(float(np.median(abs(np.diff(quiet)))/.6744897501960817/np.sqrt(2)))
            q = float(truth['charge_pe'][0]); onset = float(truth['first_active_ns'][0])
            if q < cfg.minimum_calibration_charge_pe or not np.isfinite(onset):
                continue
            curve = np.interp(grid+onset, t, y-np.median(quiet))
            samples.append({'particle':str(truth['particle_name'][0]), 'q':q,
                            'curve':curve, 'root':record['root']})
    arrays = {'grid_ns':grid, 'covariance_grid_ns':cov_grid}
    meta = {'prefix_config':asdict(cfg), 'readout':readout, 'area_per_pe_mV_ns':area_per_pe,
            'baseline_mV':float(np.median(baseline_values)),
            'electronics_sigma_mV':float(np.median(zero_sigma)), 'classes':{},
            'training_rows':len(records), 'used_training_rows':len(samples),
            'calibration_truth_scope':'Class, realized charge and first active photon used for independent calibration only',
            'physics_taus_ns':[95.,351.],
            'covariance_model':'Population native residual covariance / Q, 50% diagonal shrinkage, electronics added; pilot Q/t fixed per prefix'}
    for particle in ('alpha','gamma'):
        group = sorted([r for r in samples if r['particle']==particle], key=lambda r:r['q'])
        if len(group)<9:
            raise ValueError(f'Insufficient independent {particle} calibration pulses')
        bins = np.array_split(np.arange(len(group)), min(3, len(group)//6))
        centers, means, counts = [], [], []
        for indices in bins:
            q = np.array([group[i]['q'] for i in indices])
            curves = np.array([group[i]['curve']/group[i]['q'] for i in indices])
            centers.append(float(np.exp(np.mean(np.log(q)))))
            # Bounded weights prevent one bright pulse dominating an entire bin.
            means.append(np.average(curves, axis=0, weights=np.minimum(q,4000.)))
            counts.append(len(indices))
        arrays[particle+'_means'] = np.array(means)
        arrays[particle+'_means'][:,grid<0] = 0.
        meta['classes'][particle] = {'charge_centers_pe':centers,'counts':counts,
                                     'charge_range_pe':[group[0]['q'],group[-1]['q']]}
    model = PrefixCalibration(meta, arrays)
    residuals = np.array([(r['curve'][:len(cov_grid)] - r['q']*model.template(cov_grid,r['particle'],r['q']))/np.sqrt(r['q'])
                          for r in samples])
    residuals -= np.mean(residuals,axis=0)
    covariance = residuals.T@residuals / max(1,len(samples)-6)
    electronic_training = meta['electronics_sigma_mV']**2*np.mean([1/r['q'] for r in samples])
    covariance.flat[::len(covariance)+1] -= electronic_training
    # Projection and shrinkage keep finite-training covariance positive definite.
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    covariance = (eigenvectors*np.maximum(eigenvalues,0.))@eigenvectors.T
    shrink = cfg.covariance_shrinkage
    covariance = (1-shrink)*covariance + shrink*np.diag(np.diag(covariance))
    arrays['shot_covariance_per_pe'] = covariance
    np.savez_compressed(output/'calibration-arrays.npz', **arrays)
    (output/'calibration.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2)+'\n')
    return PrefixCalibration(meta, arrays)
