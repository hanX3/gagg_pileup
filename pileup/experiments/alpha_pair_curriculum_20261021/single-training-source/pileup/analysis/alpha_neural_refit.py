"""Waveform-derived neural proposals followed by unknown-K physical WLS.

The covariance pilot is kept separate from the enlarged candidate list. Neural
charges only gate proposals; they are never a physical amplitude prior. The
existing local WLS can add, split, merge and remove components using measured
residuals. No truth, source rate, or known multiplicity is accepted.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import numpy as np
from scipy.signal import find_peaks

from .alpha_correlated_fit import CorrelatedConfig, ResidualCalibration, fit_correlated
from .alpha_pileup_fit import FitConfig


@dataclass(frozen=True)
class NeuralProposalConfig:
    threshold: float = .65
    minimum_charge_pe: float = 250.
    edge_guard_samples: int = 262
    maximum_components: int = 1600


def propose_neural_times(time_ns, dense, config=NeuralProposalConfig()):
    """Keep scored peaks with a full training-context margin at record edges."""
    t = np.asarray(time_ns, float)
    arrays = {k: np.asarray(dense[k], float) for k in
              ('probability', 'charge_pe', 'time_offset_ns')}
    if t.ndim != 1 or len(t) < 2 or any(a.shape != t.shape for a in arrays.values()):
        raise ValueError('Neural outputs must match the measured time grid')
    if not np.isfinite(t).all() or not all(np.isfinite(a).all() for a in arrays.values()):
        raise ValueError('Finite measured grid and neural outputs required')
    if not 0 < config.threshold < 1 or config.edge_guard_samples < 0:
        raise ValueError('Invalid proposal threshold or context margin')
    peaks, _ = find_peaks(arrays['probability'], height=config.threshold, distance=4)
    times = t[peaks] + arrays['time_offset_ns'][peaks]
    margin = config.edge_guard_samples * float(t[1] - t[0])
    charge_pass = arrays['charge_pe'][peaks] >= config.minimum_charge_pe
    context_pass = (times >= t[0] + margin) & (times <= t[-1] - margin)
    keep = charge_pass & context_pass
    return times[keep], arrays['probability'][peaks][keep], dict(
        peak_count=len(peaks), charge_pass_count=int(charge_pass.sum()),
        excluded_edge_peaks=int(np.count_nonzero(charge_pass & ~context_pass)),
        context_margin_ns=margin, supported_peaks=int(keep.sum()))


def fit_neural_wls(time_ns, voltage_mV, dense, template, baseline: FitConfig,
                   calibration: ResidualCalibration, physical: CorrelatedConfig,
                   pilot, proposals=NeuralProposalConfig()):
    """Refine the union of observed-pilot and guarded neural event candidates."""
    if not physical.diagonal or physical.global_calls:
        raise ValueError('This declared hybrid uses diagonal WLS without DE')
    if pilot['status'] != 'fitted':
        return {**pilot, 'neural_refit_status': 'pilot_failed'}
    extra, probability, diagnostic = propose_neural_times(time_ns, dense, proposals)
    times = list(np.asarray(pilot['times_ns'], float))
    initial_count = len(times)
    if len(times) > proposals.maximum_components:
        raise ValueError('Observed pilot exceeds the explicit component budget')
    rejected_close = rejected_budget = 0
    for j in np.argsort(-probability, kind='stable'):
        onset = float(extra[j])
        if times and np.min(np.abs(np.asarray(times) - onset)) < physical.min_separation_ns:
            rejected_close += 1
            continue
        if len(times) >= proposals.maximum_components:
            rejected_budget += 1
            continue
        times.append(onset)
    times = np.sort(times)
    # fit_correlated uses model_mV only to define C, and times/peak_mV only
    # as initial candidates. It re-solves all amplitudes by whitened NNLS.
    # Keeping the original model here freezes C to the observed pilot even
    # though the proposal list has changed. No neural-Q penalty is introduced.
    seeded = {**pilot, 'times_ns': times.tolist(),
              'peak_mV': np.ones(len(times)).tolist()}
    fitted = fit_correlated(time_ns, voltage_mV, template, baseline, calibration,
                            physical, initial=seeded)
    fitted.update(neural_proposals=dict(
        **diagnostic, initial_physical_components=initial_count,
        added_neural_components=len(times) - initial_count,
        rejected_close=rejected_close, rejected_budget=rejected_budget,
        config=asdict(proposals)),
        neural_charge_prior=False, covariance_pilot='Unchanged observed-voltage pilot',
        fit_truth_inputs=[], rate_input_at_inference=False,
        method='Guarded CNN proposals + observed-pilot union + residual unknown-K WLS')
    return fitted
