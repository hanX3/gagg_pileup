"""Peak-normalized single-event shapes before and after the stored SiPM response.

Uses ROOT RunInfo parameters and saved waveforms; does not rerun digitization.
Time zero is the first photon on the active area, shared by all four stages of
each event. This truth-assisted display is not a blind particle classifier.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

from .single_particle_sipm import PROJECT, array, write_json


STAGES = (
    ("photons", "Before SiPM: active-area photons", "photons / bin"),
    ("current", "After SiPM: output current", "A"),
    ("analog", "After readout: analog - baseline", "mV"),
    ("adc", "After ADC: bin centre - baseline", "mV"),
)


def normalized_stages(data, photon_bin_ns=2., gate_stop_ns=2000., tail_start_ns=100.):
    """Return four sampled stages, independently scaled by their in-gate peaks.

    Voltage baselines are taken from the event's recorded readout configuration.
    Dark avalanches and electronic noise remain; negative voltage samples are
    retained. Tail fractions use sums on each uniform grid over [0, gate_stop).
    """
    if not (0 < photon_bin_ns <= tail_start_ns < gate_stop_ns):
        raise ValueError("Require 0 < photon bin <= tail start < gate stop")
    if not all(np.isclose(v / photon_bin_ns, round(v / photon_bin_ns))
               for v in (gate_stop_ns, tail_start_ns)):
        raise ValueError("Photon histogram gates must fall on bin boundaries")
    cfg = json.loads(str(data["RunInfo"]["sipm_config_json"]))
    sensor, readout = cfg["sensor"], cfg["readout"]
    optical, digi, primary = (data[k] for k in ("WaveformEvent", "SiPMEvent", "PrimaryPhoton"))
    if int(optical["n_primary"]) != 1:
        raise ValueError("This display requires one known primary per window")
    x, y = (array(optical, f"all_photon_arrival_{axis}_mm") for axis in ("x", "y"))
    active = ((x >= sensor["x_min_mm"]) & (x <= sensor["x_max_mm"]) &
              (y >= sensor["y_min_mm"]) & (y <= sensor["y_max_mm"]))
    if int(active.sum()) != int(digi["n_input_photons"] - digi["n_outside_sensor"]):
        raise ValueError("Active-area photon count disagrees with the native response")
    arrival = array(optical, "all_photon_arrival_t_ps")[active].astype(float) / 1000.
    if not len(arrival):
        raise ValueError("Cannot normalize an event with no active-area photons")
    origin = float(arrival.min())
    if origin + gate_stop_ns > float(data["RunInfo"]["t_length_ps"]) / 1000.:
        raise ValueError("The comparison gate extends beyond the saved window")
    relative = arrival - origin
    edges = np.arange(round(gate_stop_ns / photon_bin_ns) + 1) * photon_bin_ns
    counts, _ = np.histogram(relative[(relative >= 0) & (relative < gate_stop_ns)], edges)
    current_t = array(digi, "current_time_ns") - origin
    sample_t = array(digi, "sample_time_ns") - origin
    baseline = float(readout["baseline_mV"])
    lsb = (readout["adc_max_mV"] - readout["adc_min_mV"]) / 2**readout["adc_bits"]
    analog = array(digi, "analog_mV")
    adc_voltage = readout["adc_min_mV"] + (array(digi, "adc_counts") + .5) * lsb
    series = {
        "photons": ((edges[:-1] + edges[1:]) / 2, counts.astype(float)),
        "current": (current_t, array(digi, "current_A")),
        "analog": (sample_t, analog - baseline),
        "adc": (sample_t, adc_voltage - baseline),
    }
    stages, metrics = {}, {}
    for name, _, unit in STAGES:
        times, values = series[name]
        gate = (times >= 0) & (times < gate_stop_ns)
        tail = gate & (times >= tail_start_ns)
        if not np.any(gate) or not np.all(np.isfinite(values)):
            raise ValueError(f"Missing or nonfinite {name} waveform")
        peak = float(values[gate].max())
        total = float(values[gate].sum())
        if peak <= 0 or total <= 0:
            raise ValueError(f"Cannot normalize a nonpositive {name} signal")
        stages[name] = {"time_ns": times, "normalized": values / peak}
        metrics[name] = {
            "peak_before_normalization": peak, "peak_unit": unit,
            "peak_time_ns": float(times[gate][np.argmax(values[gate])]),
            "tail_fraction": float(values[tail].sum() / total),
        }
    stages["photons"]["bin_edges_ns"] = edges
    sample_gate = (sample_t >= 0) & (sample_t < gate_stop_ns)
    clipped = array(digi, "adc_clipped").astype(bool)
    valid = sample_gate & ~clipped
    quantization_error = (float(np.max(np.abs(adc_voltage[valid] - analog[valid])))
                          if valid.any() else None)
    report = {
        "particle": str(primary["primary_particle_name"]),
        "sensor_name": sensor["name"], "parameter_status": sensor["parameter_status"],
        "source_energy_MeV": float(primary["primary_energy_MeV"]),
        "deposited_energy_MeV": float(optical["edep_total_MeV"]),
        "active_area_mm": [sensor["x_max_mm"] - sensor["x_min_mm"],
                           sensor["y_max_mm"] - sensor["y_min_mm"]],
        "active_photons_full_window": int(active.sum()),
        "active_photons_in_gate": int(counts.sum()),
        "time_origin_ns": origin,
        "first_active_photon_after_source_ns": origin - float(primary["primary_time_ps"]) / 1000.,
        "baseline_subtracted_mV": baseline, "adc_lsb_mV": lsb,
        "clipped_samples_in_gate": int(clipped[sample_gate].sum()),
        "max_unclipped_adc_quantization_error_mV": quantization_error,
        "stages": metrics,
    }
    return stages, report


def plot_normalized_comparison(examples, output, photon_bin_ns=2., gate_stop_ns=2000.,
                              tail_start_ns=100.):
    """Save one four-panel PNG and numeric summary for the two loaded examples.

    The caller owns the analysis output directory; rerunning updates these two
    derived files. Original ROOT inputs and earlier figures are never changed.
    """
    if set(examples) != {"alpha", "gamma"}:
        raise ValueError("Provide exactly one alpha and one gamma example")
    configs = [json.loads(str(data["RunInfo"]["sipm_config_json"])) for data in examples.values()]
    if configs[0] != configs[1]:
        raise ValueError("Alpha and gamma must use the same recorded SiPM/readout configuration")
    comparison = {}
    for particle, data in examples.items():
        stages, report = normalized_stages(data, photon_bin_ns, gate_stop_ns, tail_start_ns)
        if report["particle"] != particle:
            raise ValueError(f"Expected {particle}, found {report['particle']}")
        comparison[particle] = (stages, report)

    os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / "jupyter/artifacts/.matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}):
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
        width, height = comparison["alpha"][1]["active_area_mm"]
        for ax, (stage, title, _) in zip(axes.flat, STAGES):
            for particle, color, label in (("alpha", "#C34B30", r"$\alpha$"),
                                            ("gamma", "#2267A0", r"$\gamma$")):
                values = comparison[particle][0][stage]
                if stage == "photons":
                    ax.stairs(values["normalized"], values["bin_edges_ns"],
                              color=color, label=label, lw=1.1, alpha=.9)
                else:
                    ax.plot(values["time_ns"], values["normalized"], color=color,
                            label=label, lw=1.25, alpha=.9)
            detail = (f"{photon_bin_ns:g} ns bins; {width:g} x {height:g} mm active area" if stage == "photons"
                      else "includes dark avalanches" if stage == "current"
                      else f"{configs[0]['readout']['sample_period_ns']:g} ns samples; pedestal removed")
            ax.set_title(f"{title}\n{detail}", fontsize=11)
            ax.set(xlim=(-20, gate_stop_ns), ylim=(-.025, 1.08), ylabel="Amplitude / own peak")
            ax.grid(alpha=.18)
            ax.legend(loc="upper right", frameon=False, fontsize=12)
        for ax in axes[-1]:
            ax.set_xlabel("Time after first active-area photon (ns)")
        fig.suptitle("S13360-6050PE: peak-normalized single-event shapes", fontsize=15)
        fig.text(.5, .012, "Same optical time zero across stages; one alpha and one selected gamma event.",
                 ha="center", fontsize=10, color="#444444")
        fig.tight_layout(rect=(0, .035, 1, .955))
        fig.savefig(output / "alpha_gamma_normalized.png", dpi=160)
        plt.close(fig)

    summary = {
        "normalization": "Each stage divided by its own positive peak in the comparison gate; voltage pedestal subtracted first",
        "time_alignment": "First active-area photon of each event, shared across all four stages; truth-assisted",
        "gate_ns": [0., gate_stop_ns], "tail_gate_ns": [tail_start_ns, gate_stop_ns],
        "photon_bin_ns": photon_bin_ns,
        "events": {particle: values[1] for particle, values in comparison.items()},
        "limits": "Two illustrative events with unequal deposited energies; peak normalization retains photon statistics and saturation effects. Dark avalanches and electronic noise retained. Tail gate is descriptive, not optimized or a PSD performance result.",
    }
    write_json(output / "summary.json", summary)
    return summary
