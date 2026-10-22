"""Linear SiPM readout, finite ADC aperture, noise, quantization and clipping."""
from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve, lfilter


def _grid_ratio(length, step, name):
    if not math.isfinite(length) or not math.isfinite(step) or length <= 0 or step <= 0:
        raise ValueError(f"{name} and internal_dt_ns must be finite and positive")
    ratio = length / step
    integer = round(ratio)
    if integer < 1 or not math.isclose(ratio, integer, rel_tol=0, abs_tol=1e-7):
        raise ValueError(f"{name} must be a positive integer multiple of internal_dt_ns")
    return integer


@dataclass(frozen=True)
class ReadoutConfig:
    model: str
    internal_dt_ns: float
    sample_period_ns: float
    aperture_ns: float
    rise_ns: float
    decay_ns: float
    transimpedance_ohm: float
    lowpass_ns: float
    baseline_mV: float
    noise_rms_mV: float
    adc_bits: int
    adc_min_mV: float
    adc_max_mV: float
    template_csv: str | None = None
    max_internal_bins: int = 5_000_000

    def __post_init__(self):
        for key in ("internal_dt_ns", "sample_period_ns", "aperture_ns", "rise_ns",
                    "decay_ns", "transimpedance_ohm"):
            value = getattr(self, key)
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{key} must be finite and positive")
        if self.rise_ns >= self.decay_ns:
            raise ValueError("Require 0 < rise_ns < decay_ns")
        for key in ("lowpass_ns", "noise_rms_mV"):
            if not math.isfinite(getattr(self, key)) or getattr(self, key) < 0:
                raise ValueError(f"{key} must be finite and nonnegative")
        if (not np.isfinite([self.baseline_mV, self.adc_min_mV, self.adc_max_mV]).all() or
                self.adc_max_mV <= self.adc_min_mV):
            raise ValueError("Invalid baseline or ADC range")
        if isinstance(self.adc_bits, bool) or not isinstance(self.adc_bits, int) or not 1 <= self.adc_bits <= 32:
            raise ValueError("adc_bits must be an integer in [1,32]")
        if (isinstance(self.max_internal_bins, bool) or
                not isinstance(self.max_internal_bins, int) or self.max_internal_bins < 1):
            raise ValueError("max_internal_bins must be a positive integer")
        _grid_ratio(self.sample_period_ns, self.internal_dt_ns, "sample_period_ns")
        _grid_ratio(self.aperture_ns, self.internal_dt_ns, "aperture_ns")
        if self.aperture_ns > self.sample_period_ns:
            raise ValueError("ADC aperture must not exceed sample period")
        if self.model == "double_exponential":
            if self.template_csv is not None:
                raise ValueError("template_csv is unused for double_exponential")
        elif self.model == "measured_voltage":
            if not self.template_csv or self.lowpass_ns != 0:
                raise ValueError("Measured full-chain SPE requires template_csv and lowpass_ns=0")
        else:
            raise ValueError("model must be double_exponential or measured_voltage")


@dataclass
class WaveformResult:
    arrays: dict
    metadata: dict


def _current_bin_charge(records, dt, n, rise, decay):
    """Exact integral in each bin for continuous avalanche times, no time rounding."""
    indices = np.floor(records["time_ns"] / dt).astype(np.int64)
    residual = (indices + 1) * dt - records["time_ns"]
    charge = records["charge_C"]
    integrals = []
    for tau in (decay, rise):
        inject_end = np.bincount(indices, weights=charge * np.exp(-residual / tau), minlength=n)
        states_end = lfilter([1.0], [1.0, -math.exp(-dt / tau)], inject_end)
        states_start = np.concatenate(([0.0], states_end[:-1]))
        inject_integral = np.bincount(indices, weights=charge * -np.expm1(-residual / tau), minlength=n)
        integrals.append(tau * (states_start * -math.expm1(-dt / tau) + inject_integral))
    return (integrals[0] - integrals[1]) / (decay - rise)


def _template_integral(path):
    data = np.genfromtxt(Path(path), delimiter=",", names=True, dtype=float)
    if data.dtype.names != ("time_ns", "voltage_mV") or data.ndim != 1 or len(data) < 2:
        raise ValueError("SPE CSV needs time_ns,voltage_mV columns and at least two rows")
    t, v = data["time_ns"], data["voltage_mV"]
    if (not np.isfinite(t).all() or not np.isfinite(v).all() or
            t[0] != 0 or np.any(np.diff(t) <= 0)):
        raise ValueError("SPE must be finite, causal and strictly increasing from time_ns=0")
    slopes = np.diff(v) / np.diff(t)
    area = np.concatenate(([0.0], np.cumsum(np.diff(t) * (v[:-1] + v[1:]) / 2)))

    def integral(query):
        query = np.clip(query, 0, t[-1])
        i = np.clip(np.searchsorted(t, query, side="right") - 1, 0, len(t) - 2)
        dx = query - t[i]
        return area[i] + v[i] * dx + slopes[i] * dx**2 / 2

    return integral, float(t[-1]), float(v[-1])


def _lowpass_bin_means(signal, dt, tau):
    """One-pole RC with bin-constant input; return exact mean in each bin."""
    if tau == 0:
        return signal
    one_minus_a = -math.expm1(-dt / tau)
    end = lfilter([one_minus_a], [1.0, -math.exp(-dt / tau)], signal)
    start = np.concatenate(([0.0], end[:-1]))
    fraction = tau * one_minus_a / dt
    return fraction * start + (1 - fraction) * signal


def synthesize_readout(avalanches, config, seed):
    """Convert avalanches to a separately reproducible linear voltage/ADC record.

    ADC samples average the first aperture_ns of each sample_period_ns, followed
    by independent Gaussian voltage noise. Analytic current-bin integrals are
    exact. RC input is bin-constant; measured-template triggers round *up* to the
    internal grid (causal delay < dt), both requiring grid convergence checks.
    """
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("Readout seed must be a nonnegative integer")
    dt = config.internal_dt_ns
    n = _grid_ratio(avalanches.window_ns, dt, "window_ns")
    stride = _grid_ratio(config.sample_period_ns, dt, "sample_period_ns")
    aperture = _grid_ratio(config.aperture_ns, dt, "aperture_ns")
    if n > config.max_internal_bins:
        raise ValueError("Readout grid exceeds max_internal_bins")
    if n % stride:
        raise ValueError("Window length must be an integer multiple of ADC sample period")
    records = avalanches.records
    if (not np.isfinite(records["time_ns"]).all() or not np.isfinite(records["charge_C"]).all() or
            np.any(records["time_ns"] < 0) or np.any(records["time_ns"] >= avalanches.window_ns) or
            np.any(records["charge_C"] < 0) or not math.isfinite(avalanches.full_charge_C) or
            avalanches.full_charge_C <= 0):
        raise ValueError("Invalid avalanche times or charges")
    current, current_time = np.array([]), np.array([])
    metadata = {"model": config.model, "seed": int(seed), "internal_dt_ns": dt,
                "sample_period_ns": config.sample_period_ns, "aperture_ns": config.aperture_ns,
                "aperture_position": "start of each sample period",
                "noise_model": "independent Gaussian voltage noise after sample aperture",
                "avalanche_charge_C": float(records["charge_C"].sum()),
                "initial_electronics_state": "zero signal history"}
    if config.model == "double_exponential":
        charge_bins = _current_bin_charge(records, dt, n, config.rise_ns, config.decay_ns)
        current = charge_bins / (dt * 1e-9)
        current_time = (np.arange(n) + 0.5) * dt
        signal = _lowpass_bin_means(current * config.transimpedance_ohm * 1000,
                                    dt, config.lowpass_ns)
        remaining_ns = avalanches.window_ns - records["time_ns"]
        tail = float(np.sum(records["charge_C"] * (
            config.decay_ns * np.exp(-remaining_ns / config.decay_ns) -
            config.rise_ns * np.exp(-remaining_ns / config.rise_ns)) /
            (config.decay_ns - config.rise_ns)))
        metadata.update({"in_window_current_charge_C": float(charge_bins.sum()),
                         "current_tail_charge_C": tail,
                         "charge_balance_residual_C": float(charge_bins.sum() + tail - records["charge_C"].sum()),
                         "current_representation": "exact bin-averaged A before readout filter",
                         "current_time_rounding_ns": 0.0})
    else:
        integral, support_ns, end_mV = _template_integral(config.template_csv)
        length = min(n, math.ceil(support_ns / dt))
        kernel = np.diff(integral(np.arange(length + 1) * dt)) / dt
        index = np.ceil(records["time_ns"] / dt).astype(np.int64)
        inside = index < n
        pulses = np.bincount(index[inside], weights=records["charge_C"][inside] / avalanches.full_charge_C,
                             minlength=n)
        signal = fftconvolve(pulses, kernel, mode="full")[:n]
        metadata.update({"template_support_ns": support_ns, "template_last_voltage_mV": end_mV,
                         "template_after_support": "zero; longer tail is not inferred",
                         "maximum_alignment_delay_ns": dt,
                         "rounded_beyond_window": int((~inside).sum()),
                         "in_window_current_charge_C": None, "current_tail_charge_C": None})
    sampled_signal = signal.reshape(-1, stride)[:, :aperture].mean(axis=1)
    times = np.arange(len(sampled_signal)) * config.sample_period_ns + config.aperture_ns / 2
    rng = np.random.default_rng(seed)
    analog = sampled_signal + config.baseline_mV + rng.normal(0, config.noise_rms_mV, len(times))
    if not np.isfinite(analog).all():
        raise ValueError("Nonfinite analog output")
    levels = 2**config.adc_bits
    lsb = (config.adc_max_mV - config.adc_min_mV) / levels
    saturation = (analog < config.adc_min_mV) | (analog >= config.adc_max_mV)
    codes = np.clip(np.floor((analog - config.adc_min_mV) / lsb), 0, levels - 1).astype(
        np.uint16 if config.adc_bits <= 16 else np.uint32)
    metadata.update({"adc_lsb_mV": lsb, "adc_code_rule": "floor then clip; each code represents its voltage bin",
                     "samples": len(times), "saturated_samples": int(saturation.sum()),
                     "analog_min_mV": float(analog.min()), "analog_max_mV": float(analog.max())})
    return WaveformResult({"time_ns": times, "signal_mV": sampled_signal, "analog_mV": analog,
                           "adc_counts": codes, "saturated": saturation,
                           "current_time_ns": current_time, "current_A": current}, metadata)
