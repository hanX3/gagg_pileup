"""Controlled charge/class-prior ablation with fixed blind candidate onsets.

The candidate union, initial labels/charges, shape charges, whitening and time
grid are shared by every arm. Only two explicit penalty strengths vary. This
diagnostic does not reproduce the production hybrid's time search/additions.
No truth, source rate or known event count is accepted by these interfaces.
"""
from dataclasses import dataclass, asdict
import hashlib
import json
import numpy as np
from scipy.linalg import qr
from scipy.optimize import nnls
from scipy.signal import find_peaks
from .alpha_gamma_regularized_fit import CalibratedBank, FrozenWhitening


@dataclass(frozen=True)
class PriorConfig:
    charge_strength: float = 0.
    class_strength: float = 0.
    max_passes: int = 6
    pair_distance_ns: float = 100.
    relative_stop: float = 1e-7


ARMS = {
    'none': PriorConfig(),
    'charge': PriorConfig(charge_strength=10.),
    'class': PriorConfig(class_strength=10.),
    'both': PriorConfig(charge_strength=10., class_strength=10.),
}


def decode_cached(t, dense, thresholds, *, minimum_charge=250., spacing_ns=4.):
    """The existing decoder in NumPy; validate caches without importing Torch."""
    candidates = []
    for k in (0, 1):
        peaks = find_peaks(dense['event_score'][k], height=thresholds[k], distance=2)[0]
        for i in peaks:
            q = float(dense['charge_pe'][k, i])
            if q >= minimum_charge:
                candidates.append((float(dense['event_score'][k, i]),
                                   float(t[i]+dense['time_offset_ns'][k, i]), k, q))
    kept = []
    for entry in sorted(candidates, reverse=True):
        if all(abs(entry[1]-old[1]) >= spacing_ns for old in kept):
            kept.append(entry)
    kept.sort(key=lambda row: row[1])
    return dict(times_ns=[r[1] for r in kept],
                particle_name=[('alpha', 'gamma')[r[2]] for r in kept],
                charge_pe=[r[3] for r in kept], event_score=[r[0] for r in kept])


def prepare(time_ns, voltage_mV, calibration, neural, dense, baseline):
    """Freeze all reconstruction inputs shared by the four penalty arms."""
    t, voltage = np.asarray(time_ns, float), np.asarray(voltage_mV, float)
    if t.ndim != 1 or voltage.shape != t.shape or not np.isfinite(voltage).all():
        raise ValueError('Invalid observed waveform')
    if len(t) < 2 or not np.allclose(np.diff(t), 2., atol=1e-8, rtol=0):
        raise ValueError('The frozen network and prior mapping require 2 ns sampling')
    initial, neural = baseline.get('all_components', baseline), neural['all_components']
    entries = []
    for priority, fitted in ((2., neural), (1., initial)):
        for j, (onset, name, q) in enumerate(zip(fitted['times_ns'], fitted['particle_name'], fitted['charge_pe'])):
            entries.append((priority, onset, 0 if name == 'alpha' else 1, q, j if priority == 2 else -1))
    kept = []
    for row in sorted(entries, reverse=True):
        if row[3] >= 40. and all(abs(row[1]-other[1]) >= 4. for other in kept):
            kept.append(row)
    kept.sort(key=lambda row: row[1])
    if len(kept) > 400:
        raise ValueError('Candidate count exceeds the declared 400-component limit')
    times = np.array([r[1] for r in kept], float)
    types = np.array([r[2] for r in kept], int)
    initial_q = np.array([r[3] for r in kept], float)
    shapes = np.maximum(initial_q, 1.)
    prior = np.zeros(len(times)); unit_weight = np.zeros(len(times)); unit_class = np.zeros((len(times), 2))
    nt = np.asarray(neural['times_ns'], float)
    for j, onset in enumerate(times):
        if not len(nt):
            continue
        m = int(np.argmin(abs(nt-onset)))
        if abs(nt[m]-onset) > 8.:
            continue
        prior[j] = neural['charge_pe'][m]
        unit_weight[j] = 1./(.35*max(500., prior[j]))**2
        # Deliberately retain the production prior's onset-index convention;
        # this experiment changes neither alignment nor score calibration.
        index = int(np.clip(np.rint(nt[m]/2.), 0, len(t)-1))
        scores = np.maximum(dense['event_score'][:, index], 1e-4)
        unit_class[j] = -2.*np.log(scores/scores.sum())
    bank = CalibratedBank(calibration)
    grid = t[::2]
    pilot = np.asarray(baseline['model_mV'])[::2]-calibration['baseline_mV']
    whitening = FrozenWhitening(pilot, 4., calibration['noise'], smooth_ns=0.)
    target = whitening.apply(voltage[::2]-calibration['baseline_mV'])
    if len(times):
        design = np.column_stack([bank.column(grid, x, k, q)
                                  for x, q in zip(times, shapes) for k in (0, 1)])
        weighted = whitening.apply(design)
        orthogonal, upper = qr(weighted, mode='economic', check_finite=False)
        reduced_target = orthogonal.T @ target
        constant = float(np.linalg.norm(target-orthogonal @ reduced_target)**2)
    else:
        upper, reduced_target, constant = np.empty((0, 0)), np.empty(0), float(target @ target)
    metadata = dict(times_ns=times.tolist(), initial_particle_name=bank.names[types].tolist(),
                    initial_charge_pe=initial_q.tolist(), shape_charge_pe=shapes.tolist(),
                    neural_prior_charge_pe=prior.tolist(), charge_weight_per_strength=unit_weight.tolist(),
                    class_cost_per_strength=unit_class.tolist(), candidate_source='Blind neural/physical union',
                    fixed_times=True, fixed_shapes=True, fixed_whitening=True,
                    interpretation='Zero explicit priors still uses the common neural candidates/initial shapes')
    identity = hashlib.sha256(json.dumps(metadata, sort_keys=True).encode()).hexdigest()
    return dict(t=t, voltage=voltage, bank=bank, calibration=calibration, times=times,
                initial_types=types, shapes=shapes, initial_q=initial_q, prior=prior,
                unit_weight=unit_weight, unit_class=unit_class, R=upper,
                target=reduced_target, constant=constant, metadata=metadata, candidate_sha256=identity)


def _solve_charge(design, target, weights, prior):
    if np.any(weights):
        design = np.vstack([design, np.diag(np.sqrt(weights))])
        target = np.r_[target, np.sqrt(weights)*prior]
    return nnls(design, target, maxiter=max(1000, 30*design.shape[1]))[0]


def fit_fixed(common, config=PriorConfig()):
    """Fit amplitudes and discrete labels; every candidate time stays fixed."""
    if config.charge_strength < 0 or config.class_strength < 0:
        raise ValueError('Prior strengths must be nonnegative')
    times, R, target = common['times'], common['R'], common['target']
    n = len(times); kinds = common['initial_types'].copy()
    weights = config.charge_strength*common['unit_weight']
    prior = common['prior']; class_cost = config.class_strength*common['unit_class']
    columns = 2*np.arange(n)+kinds
    H = R[:, columns].copy()
    q = _solve_charge(H, target, weights, prior) if n else np.empty(0)
    residual = target-H @ q

    def cost():
        return float(residual @ residual + common['constant'] + np.sum(weights*(q-prior)**2)
                     + class_cost[np.arange(n), kinds].sum())

    history = [cost()]
    for _ in range(config.max_passes if n else 0):
        for j in range(n):
            local = residual+H[:, j]*q[j]
            best = None
            for kind in (int(kinds[j]), 1-int(kinds[j])):
                h = R[:, 2*j+kind]
                denominator = float(h @ h)+weights[j]
                amp = max(0., (float(h @ local)+weights[j]*prior[j])/max(denominator, 1e-30))
                error = local-h*amp
                value = float(error @ error+weights[j]*(amp-prior[j])**2+class_cost[j, kind])
                if best is None or value < best[0]:
                    best = value, kind, amp, h, error
            _, kinds[j], q[j], H[:, j], residual = best
        for a in range(n-1):
            b = a+1
            if times[b]-times[a] > config.pair_distance_ns:
                continue
            ix = np.array([a, b]); local = residual+H[:, ix] @ q[ix]
            current = float(residual @ residual+np.sum(weights[ix]*(q[ix]-prior[ix])**2)
                            + class_cost[ix, kinds[ix]].sum())
            best = None
            for pair in ((0, 0), (0, 1), (1, 0), (1, 1)):
                design = R[:, 2*ix+pair]
                amp = _solve_charge(design, local, weights[ix], prior[ix])
                error = local-design @ amp
                value = float(error @ error+np.sum(weights[ix]*(amp-prior[ix])**2)
                              + class_cost[ix, list(pair)].sum())
                if value < current:
                    current = value; best = pair, design, amp, error
            if best is not None:
                pair, H[:, ix], q[ix], residual = best
                kinds[ix] = pair
        q = _solve_charge(H, target, weights, prior)
        residual = target-H @ q; history.append(cost())
        if history[-1] > history[-2]+1e-7*max(1., history[-2]):
            raise RuntimeError('Fixed-candidate objective increased')
        if history[-2]-history[-1] <= config.relative_stop*max(1., history[-2]):
            break
    model = common['bank'].matrix(common['t'], times, kinds, common['shapes']) @ q + common['calibration']['baseline_mV']
    reported = q >= 500.
    return dict(status='fitted', times_ns=times[reported].tolist(),
                particle_name=common['bank'].names[kinds[reported]].tolist(), charge_pe=q[reported].tolist(),
                all_components=dict(times_ns=times.tolist(), particle_name=common['bank'].names[kinds].tolist(),
                                    charge_pe=q.tolist(), shape_charge_pe=common['shapes'].tolist()),
                model_mV=model, candidate_sha256=common['candidate_sha256'], config=asdict(config),
                objective_history=history, passes=len(history)-1,
                hit_pass_budget=len(history)-1 == config.max_passes,
                weighted_voltage_sse=float(residual @ residual+common['constant']),
                charge_penalty=float(np.sum(weights*(q-prior)**2)),
                class_penalty=float(class_cost[np.arange(n), kinds].sum()),
                rms_mV=float(np.sqrt(np.mean((common['voltage']-model)**2))), fit_truth_inputs=[],
                times_and_count_fixed=True, residual_additions=[], scores_are_calibrated_probabilities=False)
