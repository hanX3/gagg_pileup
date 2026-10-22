"""Pulse duration, Poisson crowding, and conditional amplitude separation.

This is a diagnostic of a saved template and completed blind fits. It does not
launch Geant4, change the detector model, or pass truth to a reconstruction.
The two-column calculation assumes known times/count and white voltage noise;
its uncertainty multiplier is not a blind recovery efficiency.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

from analysis.alpha_gamma_waveform_fit import sampled_component


DEFAULT_INPUT = Path("jupyter/artifacts/alpha_rate_fit_20261013_211541")


def density(t_ns, shape):
    """Unit-area sampled response, including the recorded ADC aperture."""
    w = shape["slow_fraction"]
    return ((1 - w) * sampled_component(t_ns, 0., shape["tau_fast_ns"], shape["readout"])
            + w * sampled_component(t_ns, 0., shape["tau_slow_ns"], shape["readout"]))


def duration_metrics(shape, *, step_ns=.1, stop_ns=10000.):
    """Report explicit widths; exponentials have no finite full duration.

    Area durations start at -aperture/2, the support onset of the
    aperture-averaged template. Quantile times are relative to model t0=0.
    The optical lifetimes and the area-fraction fit are not re-estimated.
    """
    start = -shape["readout"]["aperture_ns"] / 2
    t = np.arange(start, stop_ns + step_ns / 2, step_ns)
    h = density(t, shape)
    cumulative = cumulative_trapezoid(h, t, initial=0.)
    area = float(cumulative[-1])
    peak = int(np.argmax(h))
    crossings = {}
    for fraction in (.01, .1, .5, .9):
        crossings[str(fraction)] = [
            float(np.interp(fraction * h[peak], h[:peak + 1], t[:peak + 1])),
            float(np.interp(fraction * h[peak], h[peak:][::-1], t[peak:][::-1]))]
    quantiles = {str(q): float(np.interp(q * area, cumulative, t))
                 for q in (.01, .05, .5, .9, .95, .99)}
    metrics = dict(
        support_start_ns=start, integrated_area=area,
        peak_time_ns=float(t[peak]), unit_area_peak_per_ns=float(h[peak]),
        fwhm_ns=crossings["0.5"][1] - crossings["0.5"][0],
        rise_10_90_ns=crossings["0.9"][0] - crossings["0.1"][0],
        area_90_duration_ns=quantiles["0.9"] - start,
        area_95_duration_ns=quantiles["0.95"] - start,
        area_99_duration_ns=quantiles["0.99"] - start,
        central_90_area_width_ns=quantiles["0.95"] - quantiles["0.05"],
        fall_to_one_percent_ns=crossings["0.01"][1] - start,
        area_over_peak_ns=area / float(h[peak]),
        area_squared_over_l2_ns=area**2 / float(np.trapezoid(h * h, t)),
        level_crossings_ns=crossings, area_quantile_times_ns=quantiles,
        integration_step_ns=step_ns, integration_stop_ns=stop_ns)
    return metrics, t, h


def poisson_overlap(rate_hz, width_ns):
    """Forward-neighbour probability for a stationary Poisson event stream."""
    occupancy = np.asarray(rate_hz) * np.asarray(width_ns) * 1e-9
    return occupancy, -np.expm1(-occupancy)


def pair_metrics(shape, gaps_ns):
    """Known-time, known-baseline, two-amplitude white-noise LS diagnostic.

    rho is the normalized inner product of sampled templates. The reported
    Gram condition number uses unit-norm columns. The standard-deviation
    factor compares each amplitude to fitting that same pulse in isolation.
    Actual SiPM residuals are neither white nor purely additive electronics.
    """
    dt = shape["readout"]["sample_period_ns"]
    t = np.arange(-100., 10000., dt) + .25
    a = density(t, shape)
    rows = []
    for gap in gaps_ns:
        if gap <= 0:
            raise ValueError("Positive gap required; coincident pulses are singular")
        b = density(t - gap, shape)
        rho = float(a @ b / np.sqrt((a @ a) * (b @ b)))
        rows.append(dict(gap_ns=float(gap), rho=rho,
                         normalized_gram_condition=(1 + rho) / (1 - rho),
                         amplitude_sd_factor=1 / np.sqrt(1 - rho**2)))
    return rows


def saved_recovery_by_gap(input_dir, summary, widths):
    """Evaluation only: nearest other signal among all q>0 alpha parents.

    Targets and time/charge success are exactly those of the saved pilot.
    Subthreshold neighbours count as neighbours; dark noise is not a parent.
    Five-microsecond guards make the width comparisons independent of record
    edges. Large last-bin gaps can be censored by the finite record.
    """
    rows, input_records = [], []
    threshold = 500.
    for window in summary["windows"]:
        path = input_dir / "truth" / (Path(window["root"]).stem + ".npz")
        with np.load(path) as data:
            arrays = {key: data[key] for key in ("first_active_ns", "emission_ns", "charge_pe")}
        tt, emission, charge = (arrays[k] for k in ("first_active_ns", "emission_ns", "charge_pe"))
        input_records.append(dict(path=str(path), selected_array_sha256={
            k: hashlib.sha256(np.ascontiguousarray(v).tobytes()).hexdigest() for k, v in arrays.items()}))
        lo, hi = window["gate_ns"]
        targets = np.flatnonzero((emission >= lo) & (emission < hi)
                                & (charge >= threshold) & np.isfinite(tt))
        assert targets.tolist() == window["adc_score"]["target_indices"]
        # Independently reconstruct success from the saved matching and charge.
        good = set()
        for i, j in window["adc_score"]["matched_pairs"]:
            assert abs(window["adc_fit"]["times_ns"][j] - tt[i]) <= 10.
            if abs(window["adc_fit"]["charge_pe"][j] / charge[i] - 1) <= .20:
                good.add(i)
        assert len(good) == window["adc_score"]["recovered_count"]
        order = np.flatnonzero((charge > 0) & np.isfinite(tt))
        order = order[np.argsort(tt[order], kind="stable")]
        previous, following = np.full(len(tt), np.inf), np.full(len(tt), np.inf)
        differences = np.diff(tt[order])
        previous[order[1:]] = differences
        following[order[:-1]] = differences
        for i in targets:
            # Actual acquisition is 0..20000 ns; summary supplies central gate.
            assert tt[i] > widths["area_99_duration_ns"]
            assert tt[i] < 20000. - widths["area_99_duration_ns"]
            gap = min(previous[i], following[i])
            rows.append(dict(source_rate_hz=window["rate_hz"], root=window["root"],
                             target_index=int(i), charge_pe=float(charge[i]),
                             nearest_gap_ns=float(gap) if np.isfinite(gap) else None,
                             next_gap_ns=float(following[i]) if np.isfinite(following[i]) else None,
                             recovered=i in good))
    rates, bins = [], []
    edges = [0., 12., 20., 50., 100., 200., 500., 1500., np.inf]
    for original in summary["rates"]:
        subset = [x for x in rows if x["source_rate_hz"] == original["rate_hz"]]
        gaps = np.array([x["nearest_gap_ns"] if x["nearest_gap_ns"] is not None else np.inf for x in subset])
        next_gaps = np.array([x["next_gap_ns"] if x["next_gap_ns"] is not None else np.inf for x in subset])
        recovered = np.array([x["recovered"] for x in subset])
        assert len(subset) == original["target_count"]
        assert recovered.sum() == original["recovered_count"]
        exposure = original["gate_time_ns"] * 1e-9
        row = {k: original[k] for k in (
            "rate_hz", "windows", "gate_time_ns", "central_injected", "central_nonzero_charge",
            "target_count", "recovered_count", "missed_count", "unmatched_predictions", "merge_flags", "split_flags",
            "recovery_fraction", "adc_rejected_windows", "clipped_fraction", "charge_relative_p16_p50_p84")}
        row.update(nonzero_charge_rate_hz=original["central_nonzero_charge"] / exposure,
                   threshold_target_rate_hz=original["target_count"] / exposure,
                   either_neighbour_within_fwhm=int(np.sum(gaps < widths["fwhm_ns"])),
                   either_neighbour_within_area99=int(np.sum(gaps < widths["area_99_duration_ns"])),
                   next_neighbour_within_area99=int(np.sum(next_gaps < widths["area_99_duration_ns"])),
                   recovered_with_neighbour_within_area99=int(np.sum(recovered & (gaps < widths["area_99_duration_ns"]))))
        rates.append(row)
        for lo, hi in zip(edges[:-1], edges[1:]):
            mask = (gaps >= lo) & (gaps < hi) if np.isfinite(hi) else gaps >= lo
            bins.append(dict(source_rate_hz=original["rate_hz"], lower_ns=lo,
                             upper_ns=float(hi) if np.isfinite(hi) else None,
                             targets=int(mask.sum()), recovered=int(recovered[mask].sum())))
    return dict(rates=rates, nearest_gap_bins=bins, targets=rows, inputs=input_records,
                neighbour_definition="Any other alpha with realized charge >0; first active-plane photon time",
                target_definition="Central source emission, realized charge >=500 pe, finite first photon time",
                success_definition="Saved one-to-one match within 10 ns and charge error <=20%",
                caveat="Post-hoc descriptive bins, small correlated windows; not an independent resolution calibration")


def validate_calculation(shape, widths, *, seed=202610141):
    """Convergence plus independent noise and Poisson simulation checks."""
    fine, _, _ = duration_metrics(shape, step_ns=.05, stop_ns=12000.)
    names = ("fwhm_ns", "rise_10_90_ns", "area_90_duration_ns", "area_99_duration_ns")
    errors = {k: abs(fine[k] - widths[k]) for k in names}
    assert max(errors.values()) < .01
    assert abs(widths["integrated_area"] - 1) < 1e-8
    rng = np.random.default_rng(seed)
    gaps = rng.exponential(size=(2, 500000))
    poisson = []
    for mu in (.1, 1., 3.):
        for name, observed, expected in (
            ("next", np.mean(gaps[0] < mu), -np.expm1(-mu)),
            ("neither_side", np.mean(np.min(gaps, axis=0) >= mu), np.exp(-2 * mu))):
            se = np.sqrt(expected * (1 - expected) / gaps.shape[1])
            assert abs(observed - expected) < 5 * se
            poisson.append(dict(mu=mu, statistic=name, observed=float(observed),
                                expected=float(expected), standard_errors=float((observed - expected) / se)))
    # Actual unconstrained linear least-squares fits, known K and arrival times.
    t = np.arange(-100., 10000., shape["readout"]["sample_period_ns"]) + .25
    h = density(t, shape)
    h /= h.max()
    noise_sd, trials = .1, 4000
    covariance = []
    for gap in (None, 4., 20., 100.):
        H = h[:, None] if gap is None else np.column_stack((h, density(t - gap, shape) / density(t, shape).max()))
        amps = np.array([100.] if gap is None else [100., 80.])
        inverse = np.linalg.pinv(H)
        noiseless = H @ amps
        estimates = []
        for _ in range(trials // 200):
            y = noiseless[:, None] + rng.normal(0., noise_sd, (len(t), 200))
            estimates.append(inverse @ y - amps[:, None])
        error = np.concatenate(estimates, axis=1)
        empirical = error.std(axis=1, ddof=1)
        predicted = noise_sd * np.sqrt(np.diag(np.linalg.inv(H.T @ H)))
        ratio = empirical / predicted
        assert np.all(abs(ratio - 1) < .08)
        assert np.all(abs(error.mean(axis=1) / predicted) < .08)
        covariance.append(dict(gap_ns=gap, empirical_sd_mV=empirical.tolist(),
                               predicted_sd_mV=predicted.tolist(), observed_over_predicted=ratio.tolist()))
    return dict(passed=True, seed=seed, width_convergence_absolute_ns=errors,
                poisson_trials=500000, poisson_checks=poisson,
                known_time_ls_trials_per_case=trials, known_time_ls_checks=covariance,
                noise_note="0.1 mV IID Gaussian only; no ADC quantization, photon noise, or SiPM fluctuation")


def figures(output, report, t, h):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 150})
    w = report["widths"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), layout="constrained")
    ax = axes[0]
    ax.plot(t / 1000, h / h.max(), color="#27669c", lw=2)
    left, right = w["level_crossings_ns"]["0.5"]
    ax.hlines(.5, left / 1000, right / 1000, color="#d87920", lw=4,
              label=f"FWHM = {w['fwhm_ns']:.1f} ns")
    ax.axvline(w["area_quantile_times_ns"]["0.9"] / 1000, color="#2d8a5d", ls="--",
               label=f"90% area duration = {w['area_90_duration_ns'] / 1000:.3f} us")
    ax.axvline(w["area_quantile_times_ns"]["0.99"] / 1000, color="#8b579c", ls="--",
               label=f"99% area duration = {w['area_99_duration_ns'] / 1000:.3f} us")
    ax.set(xlim=(-.025, 1.8), ylim=(0, 1.07), xlabel="Time relative to template t0 (us)",
           ylabel="Voltage / peak", title="Current mean GAGG + SiPM template")
    ax.legend(fontsize=9)
    rate = np.geomspace(1e4, 1e8, 301)
    ax = axes[1]
    for width, label, color in ((40., "40 ns reference window", "#888888"),
            (w["fwhm_ns"], "Current FWHM", "#d87920"),
            (w["area_90_duration_ns"], "Current 90% area duration", "#2d8a5d"),
            (w["area_99_duration_ns"], "Current 99% area duration", "#8b579c")):
        ax.semilogx(rate / 1e6, 100 * poisson_overlap(rate, width)[1], label=label, color=color)
    ax.set(xlim=(.01, 100), ylim=(0, 102), xlabel="Detector signal rate (Mcps)",
           ylabel="P(next signal within W) (%)", title="Poisson overlap; not fit failure probability")
    ax.legend(fontsize=9)
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(output / "duration_and_overlap.png")
    fig.savefig(output / "duration_and_overlap.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), layout="constrained")
    pairs = pair_metrics(report["template"], np.geomspace(2, 1000, 130))
    ax = axes[0]
    ax.loglog([x["gap_ns"] for x in pairs], [x["amplitude_sd_factor"] for x in pairs], color="#27669c", lw=2)
    ax.axhline(1, color="#888888", ls=":")
    ax.set(xlim=(2, 1000), xlabel="Two-pulse separation (ns)", ylabel="Amplitude SD / isolated-pulse SD",
           title="Known times and K=2; white noise only")
    ax = axes[1]
    for delta in (10., 20., 50., 100.):
        ax.semilogx(rate / 1e6, 100 * np.exp(-2 * rate * delta * 1e-9), label=f"delta = {delta:g} ns")
    ax.set(xlim=(.01, 100), ylim=(0, 102), xlabel="Detector signal rate (Mcps)",
           ylabel="P(no neighbour within +/- delta) (%)", title="Hypothetical separation cut; not efficiency")
    ax.legend(fontsize=9)
    for ax in axes:
        ax.grid(alpha=.2, which="both")
    fig.savefig(output / "separation_conditioning.png")
    fig.savefig(output / "separation_conditioning.pdf")
    plt.close(fig)

    bins = [x for x in report["saved_pilot"]["nearest_gap_bins"]
            if x["source_rate_hz"] == 1e7 and x["targets"]]
    fig, ax = plt.subplots(figsize=(8.4, 4.5), layout="constrained")
    values = [100 * x["recovered"] / x["targets"] for x in bins]
    labels = [f"{x['lower_ns']:g}-{x['upper_ns']:g}" if x["upper_ns"] is not None
              else f">={x['lower_ns']:g}" for x in bins]
    bars = ax.bar(labels, values, color="#3b86a7", width=.7)
    for bar, row in zip(bars, bins):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                f"{row['recovered']}/{row['targets']}", ha="center", va="bottom")
    ax.set(ylim=(0, 112), xlabel="Nearest other alpha with nonzero charge: gap (ns)",
           ylabel="Recovered targets (%)", title="Saved pilot at 10 Mcps injection: 3 windows, 204 targets")
    ax.text(.02, .96, "Time +/-10 ns AND charge +/-20%; post-hoc bins",
            transform=ax.transAxes, ha="left", va="top", fontsize=10)
    ax.grid(axis="y", alpha=.2)
    ax.set_axisbelow(True)
    fig.savefig(output / "saved_pilot_nearest_gap.png")
    fig.savefig(output / "saved_pilot_nearest_gap.pdf")
    plt.close(fig)


def analyze(input_dir=DEFAULT_INPUT):
    input_dir = Path(input_dir)
    settings = json.loads((input_dir / "analysis_config.json").read_text())
    previous = json.loads((input_dir / "summary.json").read_text())
    assert settings["fit_config"]["min_charge_pe"] == 500.
    assert settings["scoring"] == {"time_tolerance_ns": 10., "relative_charge_tolerance": .2}
    shape = settings["template"]
    widths, t, h = duration_metrics(shape)
    reference = []
    for rate in (1e5, 1e6, 1e7, 1e8):
        row = dict(detector_signal_rate_hz=rate, mean_gap_ns=1e9 / rate)
        for name in ("fwhm_ns", "area_90_duration_ns", "area_99_duration_ns"):
            mu, p = poisson_overlap(rate, widths[name])
            row[name] = dict(occupancy=float(mu), next_overlap_probability=float(p))
        reference.append(row)
    report = dict(template=shape, widths=widths, poisson_reference=reference,
                  conditional_pairs=pair_metrics(shape, (4., 8., 12., 20., 50., 100., 200., 500., 1000.)),
                  saved_pilot=saved_recovery_by_gap(input_dir, previous, widths),
                  validation=validate_calculation(shape, widths),
                  input_sha256={str(input_dir / name): hashlib.sha256((input_dir / name).read_bytes()).hexdigest()
                                for name in ("analysis_config.json", "summary.json")},
                  caveats=["Template calculation, not a detector pulse-width measurement",
                           "Poisson reference uses detector signal rate, not source injection or recovered rate",
                           "No scalar pulse width determines blind recovery; rising edge, noise, amplitudes, and nonlinearity matter",
                           "Known-time two-column LS is an optimistic conditional calculation, not an algorithm efficiency",
                           "Saved pilot nearest-neighbour bins are evaluation only and are not independent validation",
                           "No new Geant4 runs, blind refits, energy calibration, or material parameter changes"])
    return report, t, h


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    report, t, h = analyze(args.input)
    figures(args.output, report, t, h)
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(args.output / "template.npz", t_ns=t, unit_area_response_per_ns=h)
    print(json.dumps(dict(output=str(args.output), widths=report["widths"],
                          actual_rates=report["saved_pilot"]["rates"], validation_passed=report["validation"]["passed"]),
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
