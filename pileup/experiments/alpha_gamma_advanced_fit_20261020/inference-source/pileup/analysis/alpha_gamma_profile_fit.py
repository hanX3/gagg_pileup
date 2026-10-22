"""Whole-record variable projection for blind alpha/gamma pulse fitting.

Every proposed type or time change re-estimates all nonnegative amplitudes.
The covariance is frozen from a voltage-only pilot; no truth enters inference.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import nnls, minimize
from .alpha_gamma_regularized_fit import CalibratedBank, FrozenWhitening
from .alpha_gamma_model_selection_fit import fit_model_selection, ModelSelectionConfig
from .alpha_gamma_nnls import reduced_nnls


@dataclass(frozen=True)
class ProfileConfig:
    weighting: str = 'ar1'
    fit_stride: int = 2
    passes: int = 3
    pair_distance_ns: float = 250.
    time_radius_ns: float = 12.
    time_maxiter: int = 60


class ProfileSystem:
    """Cached equivalent NNLS objectives for discrete dictionary subsets."""
    def __init__(self, matrix, target):
        self.matrix = matrix
        self.target = target
        self.norms = np.linalg.norm(matrix, axis=0)
        self.design = matrix/self.norms
        self.gram = self.design.T@self.design
        self.rhs = self.design.T@target
        self.target_norm = float(target@target)
        self.calls = 0

    def solve(self, columns):
        self.calls += 1
        g = self.gram[np.ix_(columns, columns)]; b = self.rhs[columns]
        try:
            upper = cholesky(g, lower=False, check_finite=False)
            target = solve_triangular(upper.T, b, lower=True, check_finite=False)
            x = nnls(upper, target, maxiter=max(300, 20*len(columns)))[0]
        except (np.linalg.LinAlgError, RuntimeError):
            x = reduced_nnls(self.design[:, columns], self.target, maxiter=max(300, 20*len(columns)))[0]
        gradient = g@x-b
        relative_kkt = max(float(np.max(np.maximum(-gradient, 0.))),
            float(np.max(abs(gradient[x > 1e-7]))) if np.any(x > 1e-7) else 0.)/max(1.,float(np.max(abs(b))))
        if relative_kkt > 1e-7:
            x = reduced_nnls(self.design[:, columns], self.target, maxiter=max(300, 20*len(columns)))[0]
        cost = float(self.target_norm-2*b@x+x@(g@x))
        return x/self.norms[columns], cost


def _derivative(bank, t, times, kinds, shapes):
    columns = []
    for onset, kind, shape in zip(times, kinds, shapes):
        p = bank.calibration['classes'][bank.names[kind]]
        shape = np.clip(max(shape, 1.), *p['charge_range_pe'])
        slow = np.clip(p['slow_intercept']+p['slow_log_slope']*np.log(shape/1000.), .02, .98)
        x = t-onset-p['delay_ns']; fast, tail = bank.basis
        columns.append(-bank.area*((1-slow)*np.interp(x, fast.grid, fast.slope, left=0., right=0.)*fast.unit_area_peak+
            slow*np.interp(x, tail.grid, tail.slope, left=0., right=0.)*tail.unit_area_peak))
    return np.column_stack(columns)


def fit_profile(time_ns, voltage_mV, calibration, config=ProfileConfig(), *, initial=None):
    if config.weighting not in ('ar1', 'diagonal', 'uniform', 'shot'):
        raise ValueError('Unknown fixed covariance ablation')
    t_full = np.asarray(time_ns, float); y_full = np.asarray(voltage_mV, float)
    if initial is None:
        initial = fit_model_selection(t_full, y_full, calibration, ModelSelectionConfig(addition_cost=1., fit_stride=config.fit_stride))
    if initial['status'] != 'fitted' or not initial['all_components']['times_ns']:
        return initial
    bank = CalibratedBank(calibration); c = initial['all_components']
    times = np.array(c['times_ns']); kinds = np.array([0 if k == 'alpha' else 1 for k in c['particle_name']])
    charges = np.array(c['charge_pe']); shapes = np.array(c['shape_charge_pe'])
    t = t_full[::config.fit_stride]; y = y_full[::config.fit_stride]-calibration['baseline_mV']
    pilot = bank.matrix(t, times, kinds, shapes)@charges
    noise = dict(calibration['noise'])
    if config.weighting == 'diagonal': noise['shot_rho_at_2ns'] = 0.
    if config.weighting == 'shot':
        from .alpha_gamma_shot_covariance import ShotWhitening
        W = ShotWhitening(t, times, kinds, charges, shapes, calibration)
    else:
        W = FrozenWhitening(pilot, t[1]-t[0], noise, config.weighting != 'uniform', smooth_ns=0.)
    target = W.apply(y)
    history = []; changes = []; time_status = []; calls = 0
    for iteration in range(config.passes):
        design = W.apply(bank.matrix(t, np.repeat(times, 2), np.tile([0,1],len(times)), np.repeat(shapes, 2)))
        system = ProfileSystem(design, target)
        columns = 2*np.arange(len(times))+kinds
        charges, cost = system.solve(columns)
        history.append(cost)
        # A type switch profiles every other pulse amplitude, including weak ones.
        for j in np.argsort(charges)[::-1]:
            trial = columns.copy(); trial[j] ^= 1
            q, value = system.solve(trial)
            if value < cost-1e-8*max(1., cost):
                changes.append(dict(pass_index=iteration, indices=[int(j)], decrease=cost-value))
                columns, charges, cost = trial, q, value
        for j in range(len(times)-1):
            if times[j+1]-times[j] > config.pair_distance_ns or max(charges[j:j+2]) < 250.:
                continue
            best = None
            for a,b in ((0,0),(0,1),(1,0),(1,1)):
                trial = columns.copy(); trial[j:j+2] = [2*j+a,2*(j+1)+b]
                if np.array_equal(trial, columns): continue
                q, value = system.solve(trial)
                if value < cost-1e-8*max(1., cost): best = (trial, q, value); cost = value
            if best is not None:
                columns, charges, cost = best
                changes.append(dict(pass_index=iteration, indices=[j,j+1], pair=True))
        kinds = columns % 2; calls += system.calls; history.append(cost)
        starts = times.copy()
        lower = np.maximum(0., starts-config.time_radius_ns)
        upper = np.minimum(t[-1], starts+config.time_radius_ns)
        # Fixed half-gap bounds keep every candidate separated during joint search.
        middle = (starts[1:]+starts[:-1])/2
        lower[1:] = np.maximum(lower[1:], middle+4.)
        upper[:-1] = np.minimum(upper[:-1], middle-4.)
        upper = np.maximum(upper, lower+1e-8)
        def objective(onsets):
            matrix = W.apply(bank.matrix(t, onsets, kinds, shapes))
            q = reduced_nnls(matrix, target, maxiter=max(300, 20*len(times)))[0]
            residual = matrix@q-target
            derivative = W.apply(_derivative(bank,t,onsets,kinds,shapes))*q
            return float(residual@residual), 2*derivative.T@residual
        if config.time_maxiter:
            opt = minimize(objective, times, jac=True, method='L-BFGS-B', bounds=list(zip(lower,upper)),
                options=dict(maxiter=config.time_maxiter, ftol=1e-10, gtol=1e-5, maxls=25, maxcor=15))
            if opt.fun < cost: times = opt.x.copy()
            time_status.append(dict(success=bool(opt.success), iterations=int(opt.nit), message=str(opt.message)))
        else:
            time_status.append(dict(success=True, iterations=0, message='Retain voltage-derived pilot times'))
        matrix = W.apply(bank.matrix(t, times, kinds, shapes))
        charges = reduced_nnls(matrix, target, maxiter=max(300,20*len(times)))[0]
        cost = float(np.sum((matrix@charges-target)**2)); history.append(cost)
        trial_shapes = np.maximum(charges, 1.)
        trial = W.apply(bank.matrix(t, times, kinds, trial_shapes))
        q = reduced_nnls(trial, target, maxiter=max(300,20*len(times)))[0]
        value = float(np.sum((trial@q-target)**2))
        if value < cost: shapes, charges, cost = trial_shapes, q, value
        history.append(cost)
    keep = charges > 1e-6; times,kinds,charges,shapes = times[keep],kinds[keep],charges[keep],shapes[keep]
    model = calibration['baseline_mV']+bank.matrix(t_full,times,kinds,shapes)@charges
    report = charges >= 500.
    return dict(status='fitted',times_ns=times[report].tolist(),particle_name=bank.names[kinds[report]].tolist(),
        charge_pe=charges[report].tolist(),area_mV_ns=(charges[report]*bank.area).tolist(),model_mV=model,
        all_components=dict(times_ns=times.tolist(),particle_name=bank.names[kinds].tolist(),charge_pe=charges.tolist(),shape_charge_pe=shapes.tolist()),
        rms_mV=float(np.sqrt(np.mean((y_full-model)**2))),weighted_cost_history=history,
        selection=initial['selection'],config=asdict(config),fit_truth_inputs=[],type_changes=changes,
        time_optimization=time_status,profile_nnls_calls=calls,nuisance_component_count=int(np.sum(~report)),
        model_includes_nuisance_components=True,selected_components=len(times),
        objective='Fixed covariance; whole-record amplitude profile at every type/time trial; finite local discrete search')
