"""Reproducible alpha-response calibration from frozen, independent ROOT windows.

This module uses labelled calibration/validation data. It is not a blind pile-up
reconstructor. Test windows are assigned but never opened by ``prepare``.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import time

import awkward as ak
import numpy as np
import pandas as pd
from scipy.optimize import nnls
import uproot
from uproot.source.file import MemmapSource


@dataclass(frozen=True)
class ResponseConfig:
    split_tag: str = "alpha-response-20260930-v1"
    calibration_windows_per_rate: int = 30
    validation_windows_per_rate: int = 20
    processed_rates_hz: tuple = (100000, 1000000)
    bin_ns: int = 1
    response_ns: int = 5000
    extended_ns: int = 10000
    min_photons: int = 20  # Diagnostic selection, not an algorithm threshold.
    bootstrap_replicates: int = 300
    bootstrap_seed: int = 26093001


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    allow_nan=False) + "\n", encoding="utf-8")


def make_split(project, config):
    """Assign whole independent windows, using seed hashes only, before analysis."""
    catalog = json.loads((project / "jupyter/run_records/root_catalog.json").read_text())
    audit = json.loads((project / "jupyter/run_records/alpha_scan_completion_20260923.json").read_text())
    submissions = set(audit["submissions"])
    records = [r for r in catalog["files"] if r["status"] == "completed"
               and r["waveform_events"] == 1 and r["window_us"] == 20
               and r.get("submission_id") in submissions]
    if len(records) != 400 or len({r["seed"] for r in records}) != 400:
        raise ValueError("Expected the audited 400-window campaign with unique seeds.")
    result = []
    for rate in sorted({int(r["alpha_injection_rate_Hz"]) for r in records}):
        group = [dict(r) for r in records if r["alpha_injection_rate_Hz"] == rate]
        if len(group) != 100:
            raise ValueError("Expected 100 windows per rate.")
        for r in group:
            r["split_hash"] = hashlib.sha256(
                f"{config.split_tag}:{r['seed']}".encode()).hexdigest()
        group.sort(key=lambda r: r["split_hash"])
        for index, r in enumerate(group):
            a = config.calibration_windows_per_rate
            b = a + config.validation_windows_per_rate
            role = "calibration" if index < a else "validation" if index < b else "test"
            manifest = json.loads((project / r["manifest"]).read_text())
            if manifest["seed"] != r["seed"] or manifest["root_relative_path"] != r["root_file"]:
                raise ValueError("Catalog and saved manifest disagree.")
            result.append({"root_file": r["root_file"], "rate_hz": rate,
                           "seed": r["seed"], "role": role, "split_hash": r["split_hash"],
                           "manifest": r["manifest"], "macro_file": r["macro_file"],
                           "experiment": manifest["experiment"],
                           "production_root_sha256": manifest["root_sha256"],
                           "root_bytes": manifest["root_bytes"]})
    return result


def histogram_ps(times_ps, length_ns, bin_ns=1):
    """Exact left-closed/right-open bins on integer picoseconds; no rounding."""
    if length_ns <= 0 or bin_ns <= 0 or length_ns % bin_ns:
        raise ValueError("Length must be a positive integer multiple of bin width.")
    values = np.asarray(times_ps, dtype=np.int64)
    if np.any(values < 0) or np.any(values >= length_ns * 1000):
        raise ValueError("Photon outside the declared [0, T) window.")
    return np.bincount(values // (bin_ns * 1000),
                       minlength=length_ns // bin_ns).astype(np.uint32)


def rebin_counts(counts, factor):
    values = np.asarray(counts)
    if factor < 1 or values.shape[-1] % factor:
        raise ValueError("Rebinning factor must divide the final dimension.")
    return values.reshape(*values.shape[:-1], -1, factor).sum(axis=-1)


def extract_window(project, item, output, config):
    """Read only time and scalar branches; compare total/labelled histograms."""
    path = project / item["root_file"]
    if path.stat().st_size != item["root_bytes"]:
        raise ValueError(f"File size changed: {path}")
    digest = sha256(path)
    if digest != item["production_root_sha256"]:
        raise ValueError(f"ROOT hash differs from production: {path}")
    with uproot.open(path, handler=MemmapSource, num_workers=1,
                     decompression_executor=uproot.TrivialExecutor(),
                     interpretation_executor=uproot.TrivialExecutor()) as f:
        info = f["RunInfo"].arrays(library="np")
        scalar = lambda name: info[name][0].item()
        if f["WaveformEvent"].num_entries != 1 or scalar("random_seed") != item["seed"]:
            raise ValueError("Window count or seed mismatch.")
        if scalar("alpha_p11b_rate_hz") != item["rate_hz"]:
            raise ValueError("Rate mismatch.")
        if scalar("brems_enabled") or scalar("c12_capture_enabled"):
            raise ValueError("Expected pure alpha data.")
        length_ps = int(scalar("t_length_ps"))
        if length_ps % 1000:
            raise ValueError("Non-integer-ns record length.")
        length_ns = length_ps // 1000
        gate_start_ps = int(scalar("analysis_start_ps"))
        gate_end_ps = int(scalar("analysis_end_ps"))
        if (length_ns, gate_start_ps, gate_end_ps) != (20000, 5000000, 15000000):
            raise ValueError("Unexpected campaign window/gate.")
        branches = ["primary_id", "primary_particle_name", "primary_energy_MeV",
                    "primary_time_ps", "edep_total_MeV", "n_optical_photons_arrived",
                    "photon_arrival_t_ps"]
        particles = f["PrimaryPhoton"].arrays(branches, library="ak")
        arrivals = ak.to_numpy(f["WaveformEvent"]["all_photon_arrival_t_ps"].array()[0])
        waveform = histogram_ps(arrivals, length_ns, config.bin_ns)
        n_reported = int(f["WaveformEvent"]["n_optical_photons_arrived_total"].array(library="np")[0])
        n_primary = int(f["WaveformEvent"]["n_primary"].array(library="np")[0])
        edep_total = float(f["WaveformEvent"]["edep_total_MeV"].array(library="np")[0])
    if len(particles) != n_primary or int(waveform.sum()) != n_reported:
        raise ValueError("Total count or multiplicity mismatch.")
    summed = np.zeros_like(waveform, dtype=np.uint64)
    rows, samples, sample_local_indices = [], [], []
    for index, particle in enumerate(particles):
        if particle["primary_particle_name"] != "alpha":
            raise ValueError("Non-alpha primary.")
        t = ak.to_numpy(particle["photon_arrival_t_ps"]).astype(np.int64)
        q_recorded = int(particle["n_optical_photons_arrived"])
        if len(t) != q_recorded:
            raise ValueError("Labelled photon count mismatch.")
        summed += histogram_ps(t, length_ns, config.bin_ns)
        emitted = int(particle["primary_time_ps"])
        first = int(t.min()) if len(t) else None
        if first is not None and first < emitted:
            raise ValueError("Detected photon before source emission.")
        rel = t - first if first is not None else np.array([], dtype=np.int64)
        q1 = int(np.count_nonzero(rel < 1000000))
        q5 = int(np.count_nonzero(rel < config.response_ns * 1000))
        q10 = int(np.count_nonzero(rel < config.extended_ns * 1000))
        edep = float(particle["edep_total_MeV"])
        central = gate_start_ps <= emitted < gate_end_ps
        complete = first is not None and first + config.response_ns * 1000 <= length_ps
        complete10 = first is not None and first + config.extended_ns * 1000 <= length_ps
        selected = central and complete and edep > 0 and q5 >= config.min_photons
        row = {"root_file": item["root_file"], "seed": item["seed"],
               "role": item["role"], "rate_hz": item["rate_hz"],
               "primary_id": int(particle["primary_id"]),
               "source_energy_MeV": float(particle["primary_energy_MeV"]),
               "edep_MeV": edep, "emission_ns": emitted / 1000,
               "first_arrival_ns": first / 1000 if first is not None else np.nan,
               "first_delay_ns": (first-emitted) / 1000 if first is not None else np.nan,
               "central": central, "complete_5us": complete,
               "complete_10us": complete10, "response_selected": selected,
               "q_recorded": q_recorded, "q_1us": q1, "q_5us": q5, "q_10us": q10}
        rows.append(row)
        if selected:
            samples.append(histogram_ps(rel[rel < config.response_ns*1000],
                                        config.response_ns, config.bin_ns))
            sample_local_indices.append(index)
    if not np.array_equal(summed, waveform):
        raise ValueError("Total and per-primary contributions disagree by time bin.")
    if not np.isclose(sum(r["edep_MeV"] for r in rows), edep_total, rtol=3e-4, atol=2e-5):
        raise ValueError("Energy conservation between trees failed.")
    stem = path.stem
    # No labels, source rate, multiplicity, energy or first-arrival truth in inputs.
    np.savez_compressed(output / "waveforms" / f"{stem}.npz",
                        counts=waveform, bin_width_ns=np.array(config.bin_ns),
                        duration_ns=np.array(length_ns))
    # Scoring labels are a physically separate file, never part of blind inputs.
    np.savez_compressed(output / "truth" / f"{stem}.npz",
                        primary_id=np.array([r["primary_id"] for r in rows], dtype=np.int64),
                        emission_ns=np.array([r["emission_ns"] for r in rows]),
                        first_arrival_ns=np.array([r["first_arrival_ns"] for r in rows]),
                        edep_MeV=np.array([r["edep_MeV"] for r in rows]),
                        q_recorded=np.array([r["q_recorded"] for r in rows], dtype=np.int64),
                        q_5us=np.array([r["q_5us"] for r in rows], dtype=np.int64),
                        central=np.array([r["central"] for r in rows], dtype=bool))
    check = {"root_file": item["root_file"], "role": item["role"],
             "rate_hz": item["rate_hz"], "seed": item["seed"],
             "root_sha256_verified": digest, "n_primary": len(rows),
             "arrived_photons": n_reported, "selected_responses": len(samples),
             "histogram_equality_bin_ns": [1, 2, 5],
             "central_source_count": sum(r["central"] for r in rows),
             "central_signal_count": sum(r["central"] and r["edep_MeV"] > 0
                                         and r["q_recorded"] >= config.min_photons for r in rows)}
    for factor in (2, 5):
        if not np.array_equal(rebin_counts(summed, factor), rebin_counts(waveform, factor)):
            raise ValueError("Rebinned photon conservation failed.")
    return rows, samples, sample_local_indices, check


def fit_charge(energy, charge, quadratic=True):
    """Positive monotone Q=aE+bE²; charge weights are an approximation, not a likelihood claim."""
    e, q = np.asarray(energy), np.asarray(charge)
    if len(e) < 3 or np.any(e <= 0) or np.any(q <= 0):
        raise ValueError("Need at least three positive calibration observations.")
    design = np.column_stack([e, e*e]) if quadratic else e[:, None]
    weight = 1 / np.sqrt(q)
    coefficient, _ = nnls(design * weight[:, None], q * weight)
    return np.pad(coefficient, (0, 2-len(coefficient)))


def invert_charge(charge, coefficient, energy_range=None):
    """Analytic monotone inverse; optional range rejects extrapolation with NaN."""
    q = np.asarray(charge, dtype=float)
    a, b = coefficient
    if a < 0 or b < 0 or a+b <= 0 or np.any(q < 0):
        raise ValueError("Invalid monotone charge relation.")
    if b == 0:
        energy = q / a
    else:
        denominator = a + np.sqrt(a*a + 4*b*q)
        energy = np.divide(2*q, denominator, out=np.zeros_like(q), where=denominator > 0)
    if energy_range is not None:
        lo, hi = energy_range
        energy = np.where((energy >= lo) & (energy <= hi), energy, np.nan)
    return energy


def energy_metrics(energy, prediction):
    energy, prediction = np.asarray(energy), np.asarray(prediction)
    valid = np.isfinite(prediction) & (energy > 0)
    if not np.any(valid):
        return {"n": int(len(energy)), "n_in_calibration_range": 0}
    residual = prediction[valid] - energy[valid]
    relative = residual / energy[valid]
    return {"n": int(len(energy)), "n_in_calibration_range": int(valid.sum()),
            "bias_MeV": float(residual.mean()), "rmse_MeV": float(np.sqrt(np.mean(residual**2))),
            "median_relative_error": float(np.median(relative)),
            "relative_error_q16_q84": np.quantile(relative, [0.16, 0.84]).tolist(),
            "residual_q16_q84_MeV": np.quantile(residual, [0.16, 0.84]).tolist()}


def read_frozen_yield(project):
    """Read the actual frozen material table, without executing C++ or regenerating it."""
    paths = [project / f"pileup/experiments/alpha_p11b_{version}/g4/src/gagg_birks_yield_curves.inc"
             for version in ("v1", "v2")]
    if sha256(paths[0]) != sha256(paths[1]):
        raise ValueError("Production snapshots have different yield tables.")
    text = paths[1].read_text()

    def vector(name):
        match = re.search(r"\b" + re.escape(name) + r"\s*=\s*\{([^}]*)\}", text, re.S)
        if not match:
            raise ValueError(f"Missing frozen yield vector: {name}")
        return np.array([float(x.strip().replace("*MeV", ""))
                         for x in match[1].split(",") if x.strip()])

    energy, light = vector("gaggScintillationEnergy"), vector("gaggAlphaScintillationYield")
    if len(energy) != len(light) or not np.all(np.diff(energy)>0) or not np.all(np.diff(light)>0):
        raise ValueError("Frozen cumulative alpha yield must be strictly monotone.")
    return energy, light, {"path": str(paths[1].relative_to(project)), "sha256":sha256(paths[1])}


def invert_yield(charge, scale, energy_table, light_table, energy_range):
    """No extrapolation outside either the material table or calibration range."""
    if scale <= 0:
        raise ValueError("Positive effective collection scale required.")
    energy = np.interp(np.asarray(charge)/scale, light_table, energy_table, left=np.nan, right=np.nan)
    lo, hi = energy_range
    return np.where((energy>=lo)&(energy<=hi), energy, np.nan)


def calibrate(project, output, config):
    events = pd.read_csv(output / "events.csv")
    data = np.load(output / "response_samples.npz")
    selected = events.iloc[data["event_row_index"]].reset_index(drop=True)
    histograms = data["counts"].astype(float)
    fractions = histograms / histograms.sum(axis=1, keepdims=True)
    training = selected.role.eq("calibration").to_numpy()
    validation = ~training
    template = fractions[training].mean(axis=0)
    photon_weighted_template = histograms[training].sum(axis=0)
    photon_weighted_template /= photon_weighted_template.sum()
    e, q = selected.edep_MeV.to_numpy(), selected.q_5us.to_numpy()
    linear = fit_charge(e[training], q[training], quadratic=False)
    quadratic = fit_charge(e[training], q[training], quadratic=True)
    calibration_range = [float(e[training].min()), float(e[training].max())]
    table_energy, table_light, table_record = read_frozen_yield(project)
    expected_light = np.interp(e, table_energy, table_light)
    # A model-conditioned candidate, fitted ONLY on calibration windows. The scale
    # includes collection/window effects, and is not a physical SiPM PDE measurement.
    collection_scale = q[training].sum()/expected_light[training].sum()
    # All attempted inversions are retained, including unsupported extrapolations.
    selected["energy_linear_MeV"] = invert_charge(q, linear, calibration_range)
    selected["energy_quadratic_MeV"] = invert_charge(q, quadratic, calibration_range)
    selected["energy_yield_table_MeV"] = invert_yield(q, collection_scale, table_energy,
                                                      table_light, calibration_range)
    selected["tail_after_1us_fraction"] = 1 - selected.q_1us/selected.q_5us
    selected.to_csv(output / "calibration_events.csv", index=False)
    seeds = selected.seed.to_numpy()
    windows = pd.read_csv(output / "window_checks.csv")
    train_seeds = windows.loc[windows.role.eq("calibration"), "seed"].to_numpy()
    rng = np.random.default_rng(config.bootstrap_seed)
    coefficients, collection_scales = [], []
    # Resample independent windows, not photons or primaries from the same window.
    for _ in range(config.bootstrap_replicates):
        picked = rng.choice(train_seeds, size=len(train_seeds), replace=True)
        indices = np.concatenate([np.flatnonzero(training & (seeds == seed)) for seed in picked])
        coefficients.append(fit_charge(e[indices], q[indices], quadratic=True))
        collection_scales.append(q[indices].sum()/expected_light[indices].sum())
    coefficients = np.asarray(coefficients)
    np.savez_compressed(output / "alpha_response.npz", bin_edges_ns=np.arange(config.response_ns+1),
                        shape_per_bin=template, photon_weighted_shape_per_bin=photon_weighted_template,
                        charge_linear_coefficients=linear, charge_quadratic_coefficients=quadratic,
                        calibration_energy_range_MeV=np.array(calibration_range),
                        bootstrap_quadratic_coefficients=coefficients,
                        yield_energy_MeV=table_energy, yield_photons=table_light,
                        effective_collection_scale=np.array(collection_scale),
                        bootstrap_collection_scales=np.asarray(collection_scales))
    tail = selected.tail_after_1us_fraction.to_numpy()
    range_reports = []
    edges = np.array([0, .5, 1, 2, 3, 4, 5.5])
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = training & (e >= lo) & (e < hi)
        vm = validation & (e >= lo) & (e < hi)
        entry = {"edep_interval_MeV": [float(lo), float(hi)],
                 "n_calibration": int(m.sum()), "n_validation": int(vm.sum())}
        if np.any(m):
            shape = fractions[m].mean(axis=0)
            # TV at 5 ns is a descriptive distance, not a hypothesis-test p value.
            entry.update({"template_TV_at_5ns": float(.5*np.abs(rebin_counts(shape, 5)
                                               -rebin_counts(template, 5)).sum()),
                          "tail_after_1us_mean": float(tail[m].mean())})
        entry["quadratic_validation"] = energy_metrics(e[vm], selected.energy_quadratic_MeV.to_numpy()[vm])
        entry["yield_table_validation"] = energy_metrics(e[vm], selected.energy_yield_table_MeV.to_numpy()[vm])
        range_reports.append(entry)
    extended = selected.complete_10us.to_numpy() & training
    late = int((selected.q_10us-selected.q_5us)[extended].sum())
    extended_q = int(selected.q_10us[extended].sum())
    scalar_statistics = []
    for (role, rate), group in events.groupby(["role", "rate_hz"]):
        central = group[group.central]
        scalar_statistics.append({"role": role, "rate_hz": int(rate),
            "central_injected": len(central), "central_zero_edep": int((central.edep_MeV==0).sum()),
            "central_signal": int(((central.edep_MeV>0)&(central.q_recorded>=config.min_photons)).sum()),
            "usable_5us_responses": int(central.response_selected.sum())})
    predictions = {name: selected[f"energy_{name}_MeV"].to_numpy()
                   for name in ("linear", "quadratic", "yield_table")}
    common = validation & np.all(np.isfinite(np.column_stack(list(predictions.values()))),axis=1)
    result = {"calibration_n": int(training.sum()), "validation_n": int(validation.sum()),
              "calibration_energy_range_MeV": calibration_range,
              "charge_model": "Q_5us = a*E_dep + b*E_dep**2; a,b>=0; weighted least squares",
              "linear_coefficients": linear.tolist(), "quadratic_coefficients": quadratic.tolist(),
              "quadratic_coefficient_bootstrap_95pct": np.quantile(coefficients,[.025,.975],axis=0).tolist(),
              "linear_validation": energy_metrics(e[validation], selected.energy_linear_MeV.to_numpy()[validation]),
              "quadratic_validation": energy_metrics(e[validation], selected.energy_quadratic_MeV.to_numpy()[validation]),
              "yield_table_validation": energy_metrics(e[validation], predictions["yield_table"][validation]),
              "common_validation_comparison": {name: energy_metrics(e[common], values[common])
                                                for name, values in predictions.items()},
              "yield_table": table_record, "effective_collection_scale": float(collection_scale),
              "effective_collection_scale_bootstrap_95pct": np.quantile(collection_scales,[.025,.975]).tolist(),
              "calibration_tail_after_1us_mean": float(tail[training].mean()),
              "calibration_tail_after_1us_q16_q84": np.quantile(tail[training],[.16,.84]).tolist(),
              "calibration_delays_ns_q16_q50_q84": np.quantile(selected.first_delay_ns[training],[.16,.5,.84]).tolist(),
              "complete_10us_calibration_n": int(extended.sum()),
              "photons_5_to_10us": late, "photons_0_to_10us": extended_q,
              "observed_5_to_10us_fraction": late/extended_q if extended_q else None,
              "energy_groups": range_reports, "central_statistics": scalar_statistics,
              "limitations": ["Labelled response/calibration validation, not blind pile-up reconstruction.",
                "Only calibration/validation files at 1e5 and 1e6 Hz were opened; final test ROOTs remain unused.",
                "First-photon alignment is not source emission or crystal entrance timing.",
                "Q is counted in [first photon, first photon+5 us), not an ADC integral.",
                "Conditional charge-to-energy errors include response fluctuations and do not measure pile-up recovery.",
                "The observed 5-10 us tail is a finite-data diagnostic, not a prehistory convergence proof.",
                "The quadratic calibration is a project candidate, not a literature formula or universal alpha law.",
                "Frozen-yield calibration shares the simulation model; it does not independently validate real GAGG light yield.",
                "Using L_alpha(E_dep) neglects possible escaping energy/secondary-specific response; validation residuals remain."]}
    write_json(output / "summary.json", result)
    return result


def prepare(project, output, config=ResponseConfig()):
    project, output = Path(project).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / "waveforms").mkdir()
    (output / "truth").mkdir()
    start = time.monotonic()
    split = make_split(project, config)
    write_json(output / "split_manifest.json", {"config": asdict(config), "windows": split})
    chosen = [r for r in split if r["rate_hz"] in config.processed_rates_hz and r["role"] != "test"]
    rows, samples, indices, checks = [], [], [], []
    for i, item in enumerate(chosen):
        new_rows, new_samples, local_indices, check = extract_window(project, item, output, config)
        indices.extend(len(rows)+j for j in local_indices)
        rows.extend(new_rows)
        samples.extend(new_samples)
        checks.append(check)
        if (i+1) % 10 == 0:
            print(f"{i+1}/{len(chosen)} windows verified; {len(samples)} complete central responses", flush=True)
    pd.DataFrame(rows).to_csv(output / "events.csv", index=False)
    pd.DataFrame(checks).to_csv(output / "window_checks.csv", index=False)
    np.savez_compressed(output / "response_samples.npz", counts=np.asarray(samples, dtype=np.uint32),
                        event_row_index=np.asarray(indices, dtype=np.int64))
    result = calibrate(project, output, config)
    code_hash = sha256(__file__)
    (output / "alpha_response_source.py").write_bytes(Path(__file__).read_bytes())
    record = {"created_utc": datetime.now(timezone.utc).isoformat(), "config": asdict(config),
              "project_git_commit": subprocess.check_output(["git","rev-parse","HEAD"],cwd=project,text=True).strip(),
              "working_tree": subprocess.check_output(["git","-c","core.fileMode=false","status","--short"],cwd=project,text=True).splitlines(),
              "source": "pileup/analysis/alpha_response.py", "source_sha256": code_hash,
              "catalog_sha256": sha256(project / "jupyter/run_records/root_catalog.json"),
              "opened_windows": len(chosen), "reserved_test_windows": sum(r["role"]=="test" for r in split),
              "all_processed_root_hashes_verified": True, "max_workers": 1,
              "elapsed_seconds": time.monotonic()-start,
              "versions": {"python": platform.python_version(), "numpy":np.__version__,
                           "uproot":uproot.__version__, "awkward":ak.__version__, "pandas":pd.__version__}}
    write_json(output / "run_record.json", record)
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, required=True, help="New result directory; existing paths are refused.")
    args = parser.parse_args()
    prepare(args.project, args.output)
