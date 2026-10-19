"""Analyze S13360-6050PE native windows, PDE, occupancy and recovery sensitivity.

Run from the project root: python -m analysis.hamamatsu_s13360 <campaign>.
Only reads ROOT inputs; writes a separate analysis record and figures.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path
import json
import os
import shutil

import numpy as np

from .single_particle_sipm import PROJECT, array, digest, load_root, write_json
from pileup.analysis.sipm_response import (
    HC_EV_NM, PhotonHits, load_config, map_cells, read_photon_hits, simulate_avalanches,
)


def analyze(campaign, output=None):
    campaign = (PROJECT / campaign).resolve()
    manifest = json.loads((campaign / "manifest.json").read_text())
    config_path = campaign / "source" / manifest["config"]
    raw = json.loads(config_path.read_text())
    sensor, _ = load_config(config_path)
    assert sensor.name == "Hamamatsu_S13360-6050PE" and sensor.parameter_status == "datasheet_based"
    output = ((PROJECT / output).resolve() if output is not None else
              PROJECT / "jupyter/artifacts" / campaign.name / "device_characterization")
    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(__file__, output / "hamamatsu_s13360.py")
    shutil.copy2(config_path, output / "config.json")
    native_results, loaded = [], {}
    for run in manifest["runs"]:
        if run["mode"] != "enabled":
            continue
        data = load_root(PROJECT / run["root"])
        assert json.loads(str(data["RunInfo"]["sipm_config_json"])) == raw
        loaded[run["label"]] = data
        digi = data["SiPMEvent"]
        native_results.append({
            "label": run["label"], "particle": run["particle"], "root": run["root"],
            "root_sha256": digest(PROJECT / run["root"]),
            "edep_MeV": float(data["WaveformEvent"]["edep_total_MeV"]),
            "total_recorded_photons": int(digi["n_input_photons"]),
            "active_photons": int(digi["n_input_photons"]-digi["n_outside_sensor"]),
            "peak_mV": float(array(digi, "analog_mV").max()),
            "charge_C": float(digi["charge_total_C"]),
            "photon_avalanches": int(digi["n_photon_avalanches"]),
            "clipped_samples": int(digi["n_clipped_samples"]),
        })

    # Independent analytical occupancy reference for a fixed number of uniform photons.
    quiet = replace(sensor, dark_rate_hz=0., crosstalk_probability=0., afterpulse_probability=0.)
    rng = np.random.default_rng(261013101)
    occupancy = []
    p = float(np.interp(450., sensor.pde_wavelength_nm, sensor.pde_probabilities[-1]))
    for i, photons in enumerate((3600, 14400, 57600)):
        cells = rng.integers(sensor.n_cells, size=photons)
        x = sensor.x_min_mm + (cells % sensor.cells_x+.5) * (sensor.x_max_mm-sensor.x_min_mm)/sensor.cells_x
        y = sensor.y_min_mm + (cells // sensor.cells_x+.5) * (sensor.y_max_mm-sensor.y_min_mm)/sensor.cells_y
        hits = PhotonHits(np.full(photons, 10.), x, y, np.full(photons, HC_EV_NM/450.), 100.)
        observed = len(simulate_avalanches(hits, quiet, 261013201+i).records)
        n = sensor.n_cells
        a = (1-p/n)**photons
        b = (1-2*p/n)**photons
        expected = n*(1-a)
        variance = n*a*(1-a)+n*(n-1)*(b-a*a)
        assert abs(observed-expected) < 6*np.sqrt(variance)
        occupancy.append({"photons": photons, "expected_fired_pixels": expected,
                          "observed_fired_pixels": observed, "standard_deviation": float(np.sqrt(variance))})

    first = next(r for r in manifest["runs"] if r["label"] == "alpha_000")
    hits = read_photon_hits(PROJECT / first["root"])
    on_sensor = map_cells(hits, sensor) >= 0
    wavelengths = HC_EV_NM / hits.energy_eV[on_sensor]
    mean_pde = float(np.mean(np.interp(wavelengths, sensor.pde_wavelength_nm, sensor.pde_probabilities[-1])))
    sensitivity = []
    for recovery_ns in (20., 29., 40.):
        result = simulate_avalanches(hits, replace(sensor, recovery_ns=recovery_ns), 261013301)
        sensitivity.append({"recovery_ns": recovery_ns, "avalanche_count": len(result.records),
                            "charge_C": float(result.records["charge_C"].sum())})
    reference = sensitivity[1]["charge_C"]
    for row in sensitivity:
        row["relative_charge_change_from_29ns"] = row["charge_C"]/reference-1

    result = {"sensor_name": sensor.name, "parameter_status": sensor.parameter_status,
              "config_sha256": digest(config_path), "analysis_source_sha256": digest(__file__),
              "campaign": str(campaign.relative_to(PROJECT)), "native_device_results": native_results,
              "alpha_active_photon_fraction": float(on_sensor.mean()),
              "alpha_mean_fully_recovered_PDE": mean_pde,
              "short_flash_fixed_photon_occupancy": occupancy,
              "offline_recovery_sensitivity": sensitivity,
              "sensitivity_note": "Same saved alpha photon history and same NumPy seed; illustrative proxy sensitivity, not uncertainty bounds or a native/NumPy RNG equality test.",
              "manufacturer_evidence": raw["parameter_evidence"],
              "checks": {"named_device_loaded": True, "recorded_configs_identical": True,
                         "uniform_short_flash_matches_occupancy_reference": True},
              "limits": "Datasheet typical inputs with approximate PDE curves and 29 ns same-family recovery proxy. Ideal optical coupling; afterpulse/SPTR omitted; demonstrator readout. No device calibration or energy-resolution measurement."}
    write_json(output / "summary.json", result)

    os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / "jupyter/artifacts/.matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for voltage, row in zip(sensor.pde_overvoltage_V, sensor.pde_probabilities):
        axes[0].plot(sensor.pde_wavelength_nm, 100*np.asarray(row), marker=".", label=f"{voltage:g} V")
    axes[0].axvline(530, ls=":", color="gray", label="GAGG emission peak (~530 nm)")
    axes[0].set(xlabel="Wavelength (nm)", ylabel="Photon detection efficiency (%)",
                title="S13360-6050PE: approximate datasheet curve points")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=.2)
    for row in native_results:
        data = loaded[row["label"]]
        source_ns = float(data["PrimaryPhoton"]["primary_time_ps"])/1000
        axes[1].plot(array(data["SiPMEvent"], "sample_time_ns")-source_ns,
                     array(data["SiPMEvent"], "analog_mV"),
                     label=row["label"], lw=1)
    axes[1].set(xlim=(-20, 1500), xlabel="Time after source emission (ns)",
                ylabel="Readout voltage (mV)", title="S13360-6050PE: single-particle signals")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(output / "pde_and_waveforms.png", dpi=160)
    plt.close(fig)
    print(json.dumps({"output": str(output.relative_to(PROJECT)), "native_device_results": native_results,
                      "mean_PDE": mean_pde, "recovery_sensitivity": sensitivity}, ensure_ascii=False, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--output", type=Path, help="New output directory; existing results are never overwritten")
    args = parser.parse_args()
    analyze(args.campaign, args.output)
