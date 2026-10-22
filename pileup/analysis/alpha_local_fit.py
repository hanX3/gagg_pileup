"""Blind pure-alpha reconstruction with local pulse-number model searches.

Only the measured total voltage, its sample times and independent calibration
enter inference. The unchanged historical fitter supplies an optional warm
start. Fixed double-decay templates are retained; no truth or source rate is an
input. This is a project implementation, not a reproduction of a paper.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import time

import numpy as np
from scipy.optimize import LinearConstraint, minimize, nnls

from .alpha_pileup_fit import AlphaTemplate, FitConfig, candidate_peaks, fit_waveform
from .alpha_gamma_nnls import reduced_nnls


@dataclass(frozen=True)
class LocalConfig:
    report_charge_pe: float = 500.
    component_charge_pe: float = 100.
    min_separation_ns: float = 3.
    pair_distance_ns: float = 200.
    move_ns: float = 24.
    tail_ns: float = 1400.
    passes: int = 2
    max_iterations: int = 35
    split_gaps_ns: tuple = (6., 16., 40.)
    max_additions: int = 160
    residual_budget: int = 32
    base_addition_cost_mV2: float = .5
    addition_cost_mV2_per_pe: float = .01
    penalty_multiplier: float = 1.


def matrix_and_slope(template, relative):
    """The exact piecewise-linear interpolant and its within-segment slope."""
    relative = np.asarray(relative)
    h = template(relative)
    step = template.grid[1] - template.grid[0]
    ix = np.floor((relative-template.grid[0])/step).astype(np.int64)
    valid = (ix >= 0) & (ix < len(template.grid)-1)
    ix = np.clip(ix, 0, len(template.grid)-2)
    slope = np.where(valid, (template.shape[ix+1]-template.shape[ix])/step, 0.)
    return h, slope


def projected_value_gradient(t, target, times, template):
    """Unweighted variable projection over nonnegative peak amplitudes.

    The gradient is exact away from interpolation knots and active-set changes.
    It is also valid on an unchanged NNLS active set (envelope theorem).
    """
    h, slope = matrix_and_slope(template, t[:, None]-np.asarray(times)[None, :])
    amplitude = nnls(h, target, maxiter=max(100, 20*h.shape[1]))[0]
    residual = h@amplitude-target
    return float(residual@residual), -2*amplitude*(slope.T@residual), amplitude


def optimize_local(t, target, seeds, template, lower, upper, config):
    """Several initializations; jointly optimize ordered times, solve Q exactly."""
    best = None
    trials = []
    for seed in seeds:
        seed = np.sort(np.asarray(seed, float))
        k = len(seed)
        if not k or upper-lower < (k-1)*config.min_separation_ns:
            continue
        seed = np.clip(seed, lower, upper)
        if k > 1 and np.min(np.diff(seed)) < config.min_separation_ns-1e-7:
            continue
        origin = float(np.mean(seed))
        scale = max(projected_value_gradient(t, target, seed, template)[0], 1.)
        def objective(z):
            cost, gradient, _ = projected_value_gradient(t, target, origin+10*z, template)
            return cost/scale, gradient*10/scale
        constraints = ()
        if k > 1:
            diff = np.eye(k, k, 1)[:-1]-np.eye(k)[:-1]
            constraints = (LinearConstraint(diff, config.min_separation_ns/10, np.inf),)
        opt = minimize(objective, (seed-origin)/10, jac=True, method='SLSQP',
            bounds=[((lower-origin)/10, (upper-origin)/10)]*k,
            constraints=constraints,
            options={'maxiter': config.max_iterations, 'ftol': 1e-9})
        # A failed or budget-limited optimizer need not invalidate an improved
        # feasible solution. Record its status; also retain the original seed.
        for times in (seed, origin+10*opt.x):
            if np.any(times < lower-1e-6) or np.any(times > upper+1e-6):
                continue
            if k > 1 and np.min(np.diff(times)) < config.min_separation_ns-1e-5:
                continue
            cost, _, amplitude = projected_value_gradient(t, target, times, template)
            if best is None or cost < best['cost']:
                best = dict(cost=cost, times=times.copy(), amplitude=amplitude,
                            optimizer_success=bool(opt.success), iterations=int(opt.nit))
        trials.append(dict(success=bool(opt.success), iterations=int(opt.nit), status=int(opt.status)))
    return best, trials


def _fit_arrays(t, voltage, template, baseline, config, initial):
    start = time.monotonic()
    if initial is None:
        initial = fit_waveform(t, voltage, template, baseline)
    if initial['status'] != 'fitted':
        return {**initial, 'local_status': 'initial_fit_failed'}
    target = voltage-baseline.baseline_mV
    times = np.asarray(initial['times_ns'], float).copy()
    amplitudes = np.asarray(initial['peak_mV'], float).copy()
    h = template(t[:, None]-times[None, :])
    residual = target-h@amplitudes
    peak_per_pe = baseline.area_per_pe_mV_ns*template.unit_area_peak
    changes, diagnostics = [], []
    attempts = dict(same_k=0, split=0, residual_add=0, remove=0)
    additions = 0

    def update(indices, proposal, kind, reference_charge=None):
        nonlocal times, amplitudes, h, residual, additions
        if proposal is None:
            return False
        ix = np.asarray(indices, int)
        tt, aa = proposal['times'], proposal['amplitude']
        if len(aa) and np.min(aa/peak_per_pe) < config.component_charge_pe:
            return False
        old_times = times[ix].copy()
        local = residual + h[:, ix]@amplitudes[ix]
        new_h = template(t[:, None]-tt[None, :])
        new_residual = local-new_h@aa
        gain = float(residual@residual-new_residual@new_residual)
        delta_k = len(tt)-len(ix)
        qref = max(config.report_charge_pe, float(np.sum(amplitudes[ix])/peak_per_pe))
        if reference_charge is not None:
            qref = max(qref, reference_charge)
        penalty = config.penalty_multiplier*(config.base_addition_cost_mV2 +
                    config.addition_cost_mV2_per_pe*qref)
        entry=dict(kind=kind, gain_mV2=gain, reference_charge_pe=qref,
            delta_k=int(delta_k), required_gain_mV2=penalty*delta_k,
            accepted=False)
        diagnostics.append(entry)
        if gain <= penalty*delta_k+1e-7:
            return False
        keep = np.ones(len(times), bool); keep[ix] = False
        if len(tt) and len(times[keep]) and np.min(abs(times[keep, None]-tt)) < config.min_separation_ns:
            return False
        entry['accepted']=True
        if delta_k == 0:
            # The outside-neighbour bounds keep this ordered block in place.
            # Preserve all unrelated columns instead of re-interpolating NxK.
            times[ix]=tt;amplitudes[ix]=aa;h[:,ix]=new_h
        else:
            times = np.r_[times[keep], tt]; amplitudes = np.r_[amplitudes[keep], aa]
            h = np.column_stack((h[:,keep],new_h))
            order = np.argsort(times); times, amplitudes, h = times[order], amplitudes[order], h[:,order]
        residual = new_residual
        additions += max(0, delta_k)
        changes.append(dict(kind=kind, old_times_ns=old_times.tolist(), new_times_ns=tt.tolist(),
            gain_mV2=gain, penalty_mV2=penalty, delta_k=int(delta_k),
            optimizer_success=proposal.get('optimizer_success', True)))
        return True

    def local_target(ix, extra_times=()):
        anchors = np.r_[times[ix], extra_times]
        lower = max(t[0], float(anchors.min())-config.move_ns)
        upper = min(t[-1], float(anchors.max())+config.move_ns)
        outside = np.delete(times, ix)
        earlier, later = outside[outside < anchors.min()], outside[outside > anchors.max()]
        if len(earlier): lower = max(lower, earlier.max()+config.min_separation_ns)
        if len(later): upper = min(upper, later.min()-config.min_separation_ns)
        mask = (t >= lower-4) & (t <= upper+config.tail_ns)
        return lower, upper, t[mask], (residual+h[:, ix]@amplitudes[ix])[mask]

    for iteration in range(config.passes):
        # Adjacent pair updates apply at every multiplicity, including >64.
        for j in range(max(0, len(times)-1)):
            if times[j+1]-times[j] > config.pair_distance_ns:
                continue
            ix = np.array([j, j+1]); lower, upper, tx, values = local_target(ix)
            seeds = [times[ix]]
            for shift in (-4., 4.):
                seeds.append(times[ix]+np.array([shift, -shift]))
            proposal, trials = optimize_local(tx, values, seeds, template, lower, upper, config)
            attempts['same_k'] += len(trials)
            update(ix, proposal, 'same_k')

        if len(times):
            amplitudes = reduced_nnls(h, target, maxiter=max(100, 20*len(times)))[0]
            residual = target-h@amplitudes
        # Discover components with no initial candidate. This still uses only
        # observed residuals; all trials are checked in the original voltage.
        discovery = FitConfig(**{**asdict(baseline), 'min_charge_pe': config.report_charge_pe/2,
                                   'min_separation_ns': config.min_separation_ns})
        extra, _ = candidate_peaks(t, residual, template, discovery)
        ranked = []
        for onset in extra:
            if len(times) and np.min(abs(times-onset)) < config.min_separation_ns:
                continue
            col = template(t-onset)
            ranked.append((max(float(col@residual), 0.)**2/max(float(col@col), 1e-30), onset))
        for _, onset in sorted(ranked, reverse=True)[:config.residual_budget]:
            if additions >= config.max_additions:
                break
            if len(times) and np.min(abs(times-onset)) < config.min_separation_ns:
                continue
            nearest = int(np.argmin(abs(times-onset))) if len(times) else None
            ix = np.array([nearest], int) if nearest is not None else np.array([], int)
            lower, upper, tx, values = local_target(ix, [onset])
            seeds = [np.sort(np.r_[times[ix], onset])]
            proposal, trials = optimize_local(tx, values, seeds, template, lower, upper, config)
            attempts['residual_add'] += len(trials)
            update(ix, proposal, 'residual_add')

        # Actively split even when the inverse filter showed only one peak.
        # Process a snapshot in time order; find surviving anchors by proximity.
        anchors = times.copy()
        for anchor in anchors:
            if additions >= config.max_additions or not len(times):
                break
            j = int(np.argmin(abs(times-anchor)))
            if abs(times[j]-anchor) > config.move_ns:
                continue
            ix = np.array([j]); lower, upper, tx, values = local_target(ix)
            seeds = []
            for gap in config.split_gaps_ns:
                seeds.extend(([times[j]-gap/2, times[j]+gap/2],
                              [times[j], times[j]+gap], [times[j]-gap, times[j]]))
            proposal, trials = optimize_local(tx, values, seeds, template, lower, upper, config)
            attempts['split'] += len(trials)
            update(ix, proposal, 'split')

        # Remove/merge proposals compete with their current model. Complexity
        # penalties are local and charge-scaled, not a calibrated chi-square.
        j = 0
        while j < len(times):
            ix = np.array([j]); cost = float((residual+h[:, ix]@amplitudes[ix]) @
                                           (residual+h[:, ix]@amplitudes[ix]))
            remove = dict(times=np.array([]), amplitude=np.array([]), cost=cost)
            attempts['remove'] += 1
            if update(ix, remove, 'remove'):
                continue
            if j+1 < len(times) and times[j+1]-times[j] <= config.pair_distance_ns:
                ix = np.array([j, j+1]); lower, upper, tx, values = local_target(ix)
                center = float(np.average(times[ix], weights=np.maximum(amplitudes[ix], 1e-20)))
                proposal, trials = optimize_local(tx, values, [[center]], template, lower, upper, config)
                attempts['remove'] += len(trials)
                if update(ix, proposal, 'merge'):
                    continue
            j += 1
        if len(times):
            amplitudes = reduced_nnls(h, target, maxiter=max(100, 20*len(times)))[0]
            residual = target-h@amplitudes

    all_charge = amplitudes/peak_per_pe
    keep = all_charge >= config.report_charge_pe
    area = amplitudes[keep]/template.unit_area_peak
    return dict(status='fitted', times_ns=times[keep].tolist(), peak_mV=amplitudes[keep].tolist(),
        charge_pe=all_charge[keep].tolist(), area_mV_ns=area.tolist(),
        all_components=dict(times_ns=times.tolist(), peak_mV=amplitudes.tolist(), charge_pe=all_charge.tolist()),
        model_mV=target-residual+baseline.baseline_mV,
        rms_mV=float(np.sqrt(np.mean(residual**2))), local_changes=changes,
        selection_diagnostics=diagnostics, trials=attempts,
        refinement_seconds=time.monotonic()-start, initial_reported_components=len(initial['times_ns']),
        all_component_count=len(times), config=asdict(config),
        penalty_note='Local full-record SSE gain; charge scale frozen across each comparison; not a p-value')


def fit_local(time_ns, voltage_mV, template: AlphaTemplate, baseline=FitConfig(),
              config=LocalConfig(), *, initial=None):
    """Infer K/t/Q from voltage. Optional initial is an observed-waveform fit."""
    t, y = np.asarray(time_ns, float), np.asarray(voltage_mV, float)
    if t.ndim != 1 or y.shape != t.shape or len(t) < 10 or not np.isfinite(y).all():
        raise ValueError('Invalid voltage waveform')
    if not np.allclose(np.diff(t), t[1]-t[0], rtol=0, atol=1e-8):
        raise ValueError('Uniform samples required')
    if abs(t[0]) > t[1]-t[0]:
        raise ValueError('Full record beginning near zero required')
    if config.min_separation_ns <= 0 or config.penalty_multiplier <= 0:
        raise ValueError('Positive separation and selection penalty required')
    return _fit_arrays(t, y, template, baseline, config, initial)
