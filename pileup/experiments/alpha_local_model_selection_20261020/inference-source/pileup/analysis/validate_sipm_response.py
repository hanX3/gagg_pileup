"""Reproducible synthetic physics checks and a three-model optical ROOT comparison."""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np
import uproot
from uproot.source.file import MemmapSource

from .sipm_response import (AvalancheResult, HC_EV_NM, PhotonHits, digitize_root,
                            load_config, read_photon_hits, sha256, simulate_avalanches, write_json)
from .sipm_readout import synthesize_readout


def check_recovery(base):
    """Two flashes hit each of 3000 independent cells once per flash."""
    s = replace(base, cells_x=60, cells_y=50, pde_mode="constant", constant_pde=1,
                pde_wavelength_nm=(), pde_overvoltage_V=(), pde_probabilities=(),
                recovery_enabled=True, dark_rate_hz=0, crosstalk_probability=0,
                afterpulse_probability=0, transit_delay_mean_ns=0, transit_delay_sigma_ns=0)
    ix, iy = np.meshgrid(np.arange(s.cells_x), np.arange(s.cells_y))
    x = s.x_min_mm+(ix.ravel()+.5)*(s.x_max_mm-s.x_min_mm)/s.cells_x
    y = s.y_min_mm+(iy.ravel()+.5)*(s.y_max_mm-s.y_min_mm)/s.cells_y
    rows = []
    for k, delay in enumerate(np.array([0, .01, .05, .1, .25, .5, 1, 2, 4, 8])*s.recovery_ns):
        n = s.n_cells
        h = PhotonHits(np.concatenate((np.zeros(n), np.full(n, delay))), np.tile(x, 2), np.tile(y, 2),
                       np.full(2*n, HC_EV_NM/500), max(1000, delay+1))
        a = simulate_avalanches(h, s, 261010200+k).records
        second = a[a["photon_index"] >= n]
        recovery = 1-np.exp(-delay/s.recovery_ns)
        effective_v = max(recovery*s.overvoltage_V-s.trigger_threshold_V, 0)
        probability = (1-np.exp(-effective_v/s.trigger_scale_V)) / (
            1-np.exp(-(s.overvoltage_V-s.trigger_threshold_V)/s.trigger_scale_V))
        error = abs(len(second)/n-probability)
        tolerance = 6*np.sqrt(probability*(1-probability)/n)+1/n
        if error > tolerance:
            raise AssertionError("Double-flash recovery differs from the stated model")
        measured_charge = float(second["charge_C"].sum()/s.full_charge_C/n)
        if len(second) and not np.allclose(second["charge_C"]/s.full_charge_C, recovery, rtol=1e-12, atol=1e-12):
            raise AssertionError("Recovery charge mismatch")
        rows.append({"delay_ns": float(delay), "recovery_fraction": float(recovery),
                     "trigger_expected": float(probability), "trigger_observed": len(second)/n,
                     "charge_per_second_photon_pe": measured_charge,
                     "charge_expected_pe": float(recovery*probability), "cells": n})
    return rows


def check_saturation(base):
    """Poisson flashes: exact independent-cell occupancy expectation and variance."""
    s = replace(base, cells_x=16, cells_y=16, pde_mode="constant", constant_pde=.35,
                pde_wavelength_nm=(), pde_overvoltage_V=(), pde_probabilities=(),
                recovery_enabled=True, dark_rate_hz=0, crosstalk_probability=0,
                afterpulse_probability=0, transit_delay_mean_ns=0, transit_delay_sigma_ns=0)
    rng, rows, trials = np.random.default_rng(261010300), [], 100
    for j, ratio in enumerate((.1, .5, 1, 2, 5, 10)):
        mean = s.n_cells*ratio
        counts = []
        for trial in range(trials):
            n = rng.poisson(mean)
            h = PhotonHits(np.zeros(n), rng.uniform(s.x_min_mm, s.x_max_mm, n),
                           rng.uniform(s.y_min_mm, s.y_max_mm, n), np.full(n, HC_EV_NM/500), 1000)
            counts.append(len(simulate_avalanches(h, s, 261010400+j*trials+trial).records))
        occupancy = 1-np.exp(-mean*s.constant_pde/s.n_cells)
        expected = s.n_cells*occupancy
        se = np.sqrt(s.n_cells*occupancy*(1-occupancy)/trials)
        if abs(np.mean(counts)-expected) > 6*se:
            raise AssertionError("Short-flash occupancy is outside statistical tolerance")
        rows.append({"mean_photons": mean, "cells": s.n_cells, "pde": s.constant_pde,
                     "trials": trials, "expected_avalanches": float(expected),
                     "observed_avalanches": float(np.mean(counts)), "standard_error": float(se)})
    return rows


def controlled_alpha_checks(root_path, base, readout, output):
    """Use labelled truth only to CONSTRUCT controlled inputs, not to reconstruct.

    A selected alpha's recorded photons are aligned to its first photon and
    copied twice. This is a paired optical-history benchmark, not two independent
    Geant4 alpha events and not blind separation of the input pileup window.
    """
    aggregate = read_photon_hits(root_path)
    names = ["primary_id", "primary_particle_name", "primary_time_ps", "primary_energy_MeV",
             "edep_total_MeV", "photon_arrival_t_ps", "photon_arrival_x_mm",
             "photon_arrival_y_mm", "photon_energy_eV"]
    with uproot.open(root_path, handler=MemmapSource, num_workers=1,
                     decompression_executor=uproot.TrivialExecutor(),
                     interpretation_executor=uproot.TrivialExecutor()) as root:
        primaries = root["PrimaryPhoton"].arrays(names, library="ak")
    selected = next((p for p in primaries if p["primary_particle_name"] == "alpha" and
                     len(p["photon_arrival_t_ps"]) >= 1000 and
                     p["primary_time_ps"]/1000 < aggregate.window_ns-5000), None)
    if selected is None:
        raise ValueError("No sufficiently long, nonempty alpha history for controlled checks")
    original_times = ak.to_numpy(selected["photon_arrival_t_ps"]).astype(float)/1000
    relative = original_times-original_times.min()
    x = ak.to_numpy(selected["photon_arrival_x_mm"])
    y = ak.to_numpy(selected["photon_arrival_y_mm"])
    energy = ak.to_numpy(selected["photon_energy_eV"])
    sensor = replace(base, dark_rate_hz=0, crosstalk_probability=0, afterpulse_probability=0)
    electronics = replace(readout, noise_rms_mV=0)
    directory = output/"controlled_alpha"
    directory.mkdir()
    provenance = {"construction": "truth-selected alpha; first-photon alignment; identical optical histories copied for pairs",
                  "input_root_sha256": sha256(root_path), "primary_id": int(selected["primary_id"]),
                  "primary_energy_MeV": float(selected["primary_energy_MeV"]),
                  "edep_MeV": float(selected["edep_total_MeV"]), "photons_per_copy": len(relative),
                  "first_photon_ns": float(original_times.min()), "relative_last_photon_ns": float(relative.max()),
                  "sensor": asdict(sensor), "readout": asdict(electronics),
                  "avalanche_seed": 261010900, "readout_seed": 261010901,
                  "initial_state": "fully recovered; no background or correlated noise"}
    write_json(directory/"construction.json", provenance)
    np.savez_compressed(directory/"source_photons.npz", relative_time_ns=relative,
                        x_mm=x, y_mm=y, photon_energy_eV=energy)
    rows = []
    for delay in (None, 0, 30, 100, 300, 1000):
        name = "single" if delay is None else f"pair_{delay}ns"
        case = directory/name
        case.mkdir()
        copies = 1 if delay is None else 2
        times = relative+5000 if delay is None else np.concatenate((relative+5000, relative+5000+delay))
        order = np.argsort(times, kind="stable")
        hits = PhotonHits(times[order], np.tile(x, copies)[order], np.tile(y, copies)[order],
                          np.tile(energy, copies)[order], aggregate.window_ns,
                          aggregate.gate_start_ns, aggregate.gate_end_ns)
        row = {"case": name, "separation_ns": delay, "copies": copies, "input_photons": len(times)}
        for mode, enabled in (("linear", False), ("recovery", True)):
            a = simulate_avalanches(hits, replace(sensor, recovery_enabled=enabled), 261010900)
            w = synthesize_readout(a, electronics, 261010901)
            row[f"{mode}_charge_pe"] = a.metadata["avalanche_charge_C"]/sensor.full_charge_C
            row[f"{mode}_avalanches"] = len(a.records)
            row[f"{mode}_saturated_samples"] = w.metadata["saturated_samples"]
            np.savez_compressed(case/f"{mode}_avalanches.npz", records=a.records)
            np.savez_compressed(case/f"{mode}_waveform.npz", **w.arrays)
            write_json(case/f"{mode}_metadata.json", {"avalanche": a.metadata, "readout": w.metadata,
                                                       "recovery_enabled": enabled})
        # Common photon/transit RNGs, monotone PDE and no secondary noise imply this inequality.
        if row["recovery_charge_pe"] > row["linear_charge_pe"]*(1+1e-12):
            raise AssertionError("Recovery-only charge exceeds coupled full-recovery baseline")
        row["charge_ratio"] = row["recovery_charge_pe"]/row["linear_charge_pe"]
        rows.append(row)
    write_json(directory/"summary.json", rows)
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    paired = rows[1:]
    ax.plot([r["separation_ns"] for r in paired], [r["charge_ratio"] for r in paired], "o-",
            label="Two copies of selected alpha optical history")
    ax.axhline(rows[0]["charge_ratio"], ls="--", color="gray", label="Single selected alpha")
    ax.set(xlabel="First-photon separation [ns]", ylabel="Charge: recovery / linear baseline",
           title=f"Controlled alpha inputs, {base.name}; one RNG realization")
    ax.grid(alpha=.25)
    ax.legend(fontsize=8)
    fig.savefig(output/"controlled_alpha.png", dpi=160)
    plt.close(fig)
    return {"construction": provenance, "cases": rows}


def validate(root, config_path, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    sensor, readout = load_config(config_path)
    (output/"validate_sipm_response.py").write_bytes(Path(__file__).read_bytes())
    recovery = check_recovery(sensor)
    saturation = check_saturation(sensor)
    controlled = controlled_alpha_checks(root, sensor, readout, output)
    variants = {
        "linear": replace(sensor, recovery_enabled=False, dark_rate_hz=0,
                          crosstalk_probability=0, afterpulse_probability=0),
        "recovery": replace(sensor, dark_rate_hz=0, crosstalk_probability=0, afterpulse_probability=0),
        "full": sensor,
    }
    results = {}
    for name, model in variants.items():
        config = output/f"{name}_config.json"
        write_json(config, {"schema_version": 1, "sensor": asdict(model), "readout": asdict(readout)})
        results[name] = digitize_root(root, config, output/name, 2610101011, 2610101012)
    full = results["full"]
    with np.load(output/"full/avalanches.npz", allow_pickle=False) as data:
        avalanches = AvalancheResult(data["records"], full["window_ns"], sensor.full_charge_C,
                                    data["last_fire_ns"], full["avalanche_result"])
    signals, grid_rows = [], []
    # Keep sampling/aperture fixed while refining the internal RC integration.
    for factor in (2, 1, .5, .25):
        dt = readout.internal_dt_ns*factor
        waveform = synthesize_readout(avalanches, replace(readout, internal_dt_ns=dt, noise_rms_mV=0), 2610101012)
        signals.append(waveform.arrays["signal_mV"])
        grid_rows.append({"internal_dt_ns": dt})
    norm = np.linalg.norm(signals[-1])
    for row, signal in zip(grid_rows, signals):
        row["relative_l2_to_finest"] = float(np.linalg.norm(signal-signals[-1])/norm) if norm else 0
    if not (grid_rows[0]["relative_l2_to_finest"] > grid_rows[1]["relative_l2_to_finest"] >
            grid_rows[2]["relative_l2_to_finest"]):
        raise AssertionError("ROOT waveform grid refinement did not converge")
    if grid_rows[1]["relative_l2_to_finest"] > .01:
        raise AssertionError("Configured readout time step differs by more than 1% in L2 from finer grid")
    summary = {"status": "passed", "sensor_name": sensor.name,
               "parameter_status": sensor.parameter_status,
               "synthetic_checks_note": "Analytical test limits only: double-flash uses 3000 cells and constant PDE=1; Poisson occupancy uses 256 cells and constant PDE=0.35; both disable dark/correlated noise and transit delay. They are not full-device predictions. The full ROOT response uses the supplied device configuration unchanged.",
               "input_root": full["input_root"], "input_sha256": full["input_sha256"],
               "double_flash": recovery, "poisson_flash": saturation, "controlled_alpha": controlled,
               "readout_grid": grid_rows,
               "models": {}, "elapsed_seconds": time.perf_counter()-started,
               "limitations": ["No real-device calibration or experimental validation",
                               "One optical ROOT window, not the 400-window production scan",
                               "Initial fully recovered cells; pre-window history/steady state not validated",
                               "No claim of reconstructed alpha energy resolution or pileup efficiency"]}
    for name, result in results.items():
        av = result["avalanche_result"]
        ro = result["readout_result"]
        summary["models"][name] = {"counters": av["counters"],
                                   "avalanche_charge_pe": av["avalanche_charge_C"]/sensor.full_charge_C,
                                   "gate_charge_pe": av["gate_avalanche_charge_C"]/sensor.full_charge_C,
                                   "saturated_samples": ro["saturated_samples"],
                                   "charge_balance_residual_C": ro.get("charge_balance_residual_C"),
                                   "elapsed_seconds": result["elapsed_seconds"]}
    write_json(output/"summary.json", summary)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    x = np.array([r["delay_ns"] for r in recovery])/sensor.recovery_ns
    axes[0].plot(x, [r["charge_expected_pe"] for r in recovery], label="Model expectation")
    axes[0].scatter(x, [r["charge_per_second_photon_pe"] for r in recovery], label="Simulated")
    axes[0].set(xlabel="Flash separation / recovery time", ylabel="Second-flash charge / photon [p.e.]",
                title="Controlled recovery limit (3000 cells, PDE=1)")
    x = np.array([r["mean_photons"] for r in saturation])
    axes[1].plot(x, [r["expected_avalanches"] for r in saturation], label="Occupancy expectation")
    axes[1].errorbar(x, [r["observed_avalanches"] for r in saturation],
                     yerr=[r["standard_error"] for r in saturation], fmt="o", label="Simulated mean +/- SE")
    axes[1].set(xlabel="Mean incident photons", ylabel="Fired cells", title="Short Poisson flash, 256 cells")
    for ax in axes:
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
    fig.savefig(output/"physics_checks.png", dpi=160)
    plt.close(fig)

    hits = read_photon_hits(root)
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), constrained_layout=True)
    bins = np.arange(0, hits.window_ns+2, 2)
    counts, _ = np.histogram(hits.time_ns, bins)
    axes[0].step((bins[:-1]+1)/1000, counts, where="mid", color="gray", lw=.8)
    axes[0].set(ylabel="Arriving photons / 2 ns", title=f"Optical ROOT to SiPM / ADC: {sensor.name}")
    for name in variants:
        with np.load(output/name/"waveform.npz", allow_pickle=False) as w:
            axes[1].plot(w["time_ns"]/1000, w["signal_mV"], lw=.8, label=name)
    with np.load(output/"full/waveform.npz", allow_pickle=False) as w:
        peak_time = w["time_ns"][int(np.argmax(w["signal_mV"]))]
        mask = (w["time_ns"] >= peak_time-100) & (w["time_ns"] <= peak_time+300)
        axes[2].plot(w["time_ns"][mask], w["analog_mV"][mask], label="Analog with baseline/noise", lw=1)
        voltage = readout.adc_min_mV+(w["adc_counts"].astype(float)+.5)*full["readout_result"]["adc_lsb_mV"]
        axes[2].step(w["time_ns"][mask], voltage[mask], where="mid", label="ADC bin center", lw=.8)
    axes[1].set(xlabel="Time [us]", ylabel="Signal [mV]")
    axes[2].set(xlabel="Time [ns]", ylabel="Voltage [mV]", title="Zoom around largest signal")
    for ax in axes:
        ax.grid(alpha=.2)
    axes[1].legend()
    axes[2].legend()
    fig.savefig(output/"root_waveforms.png", dpi=160)
    plt.close(fig)
    summary["artifacts_sha256"] = {name: sha256(output/name) for name in
                                   ("physics_checks.png", "controlled_alpha.png", "root_waveforms.png")}
    summary["elapsed_seconds"] = time.perf_counter()-started
    write_json(output/"summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = validate(args.root, args.config, args.output)
    print(json.dumps({k: result[k] for k in ("status", "parameter_status", "models", "readout_grid", "elapsed_seconds")},
                     indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
