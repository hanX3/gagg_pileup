"""Offline SiPM avalanches from one complete optical-photon ROOT window.

Configurations distinguish hypothetical, datasheet-based and calibrated devices.
Times are ns, positions mm, energies eV, and avalanche charges C. See
docs/sipm-response.md for probability conventions and physical limitations.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import hashlib
import heapq
import json
import math
from pathlib import Path
import platform
import subprocess
import time

import awkward as ak
import numpy as np
import scipy
import uproot
from uproot.source.file import MemmapSource

from .sipm_readout import ReadoutConfig, synthesize_readout


ELECTRON_CHARGE_C = 1.602176634e-19
HC_EV_NM = 1239.8419843320026
KINDS = {0: "photon", 1: "dark", 2: "crosstalk", 3: "afterpulse"}
AVALANCHE_DTYPE = np.dtype([
    ("time_ns", "f8"), ("cell_id", "i8"), ("charge_C", "f8"),
    ("kind", "i1"), ("parent_id", "i8"), ("photon_index", "i8"),
    ("recovery_fraction", "f8"),
])


def _positive(name, value, allow_zero=False):
    if not np.isfinite(value) or (value < 0 if allow_zero else value <= 0):
        raise ValueError(f"{name} must be finite and {'nonnegative' if allow_zero else 'positive'}")


@dataclass(frozen=True)
class SensorConfig:
    name: str
    parameter_status: str
    parameter_source: str
    temperature_C: float
    cells_x: int
    cells_y: int
    x_min_mm: float
    x_max_mm: float
    y_min_mm: float
    y_max_mm: float
    overvoltage_V: float
    gain_electrons: float
    recovery_ns: float
    # Relative trigger probability for constant PDE and noise candidates.
    trigger_scale_V: float
    trigger_threshold_V: float = 0.0
    recovery_enabled: bool = True
    pde_mode: str = "constant"
    constant_pde: float | None = None
    # Rows = overvoltages, columns = wavelengths; no extrapolation.
    pde_wavelength_nm: tuple = ()
    pde_overvoltage_V: tuple = ()
    pde_probabilities: tuple = ()
    dark_rate_hz: float = 0.0
    crosstalk_probability: float = 0.0
    crosstalk_delay_ns: float = 0.0
    afterpulse_probability: float = 0.0
    afterpulse_delay_ns: float = 1.0
    transit_delay_mean_ns: float = 0.0
    transit_delay_sigma_ns: float = 0.0
    max_candidates: int = 10_000_000

    def __post_init__(self):
        if not self.name or not self.parameter_source:
            raise ValueError("Device name and parameter_source are required")
        if self.parameter_status not in ("hypothetical", "datasheet_based", "calibrated"):
            raise ValueError("parameter_status must be hypothetical, datasheet_based or calibrated")
        for key in ("cells_x", "cells_y", "max_candidates"):
            v = getattr(self, key)
            if isinstance(v, bool) or not isinstance(v, int) or v <= 0:
                raise ValueError(f"{key} must be a positive integer")
        if not isinstance(self.recovery_enabled, bool):
            raise ValueError("recovery_enabled must be a boolean")
        for lo, hi in ((self.x_min_mm, self.x_max_mm), (self.y_min_mm, self.y_max_mm)):
            if not np.isfinite([lo, hi]).all() or hi <= lo:
                raise ValueError("Invalid sensor footprint")
        if not math.isfinite(self.temperature_C):
            raise ValueError("temperature_C must be finite")
        for key in ("overvoltage_V", "gain_electrons", "recovery_ns", "trigger_scale_V",
                    "afterpulse_delay_ns"):
            _positive(key, getattr(self, key))
        for key in ("trigger_threshold_V", "dark_rate_hz", "crosstalk_delay_ns",
                    "transit_delay_mean_ns", "transit_delay_sigma_ns"):
            _positive(key, getattr(self, key), allow_zero=True)
        if self.trigger_threshold_V >= self.overvoltage_V:
            raise ValueError("Trigger threshold must be below the operating overvoltage")
        for key in ("crosstalk_probability", "afterpulse_probability"):
            if not 0 <= getattr(self, key) < 1:
                raise ValueError(f"{key} must be in [0, 1)")
        if self.crosstalk_probability + self.afterpulse_probability >= 1:
            raise ValueError("The proposed secondary branching mean must be below 1")
        if self.transit_delay_sigma_ns > 0 and self.transit_delay_mean_ns <= 0:
            raise ValueError("A causal gamma transit distribution requires a positive mean")
        if self.pde_mode == "constant":
            if self.constant_pde is None or not 0 <= self.constant_pde <= 1:
                raise ValueError("constant_pde in [0,1] is required in constant mode")
            if any(len(v) for v in (self.pde_wavelength_nm, self.pde_overvoltage_V,
                                    self.pde_probabilities)):
                raise ValueError("Remove unused PDE tables in constant mode")
        elif self.pde_mode == "table":
            if self.constant_pde is not None:
                raise ValueError("constant_pde must be null in table mode")
            w, v, p = map(np.asarray, (self.pde_wavelength_nm, self.pde_overvoltage_V,
                                      self.pde_probabilities))
            if (w.ndim != 1 or v.ndim != 1 or len(w) < 2 or len(v) < 2 or
                    p.shape != (len(v), len(w))):
                raise ValueError("PDE table shape must be [overvoltage][wavelength]")
            if (not all(np.isfinite(a).all() for a in (w, v, p)) or
                    np.any(np.diff(w) <= 0) or np.any(np.diff(v) <= 0) or w[0] <= 0 or
                    v[0] != 0 or v[-1] < self.overvoltage_V or
                    np.any((p < 0) | (p > 1)) or np.any(p[0] != 0) or
                    np.any(np.diff(p, axis=0) < 0)):
                raise ValueError("PDE table must cover 0..operating voltage with physical probabilities")
        else:
            raise ValueError("pde_mode must be constant or table")

    @property
    def n_cells(self):
        return self.cells_x * self.cells_y

    @property
    def full_charge_C(self):
        return ELECTRON_CHARGE_C * self.gain_electrons


@dataclass
class PhotonHits:
    time_ns: np.ndarray
    x_mm: np.ndarray
    y_mm: np.ndarray
    energy_eV: np.ndarray | None
    window_ns: float
    gate_start_ns: float = 0.0
    gate_end_ns: float | None = None
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        _positive("window_ns", self.window_ns)
        if self.gate_end_ns is None:
            self.gate_end_ns = self.window_ns
        if not 0 <= self.gate_start_ns < self.gate_end_ns <= self.window_ns:
            raise ValueError("Invalid evaluation gate")
        for name in ("time_ns", "x_mm", "y_mm", "energy_eV"):
            value = getattr(self, name)
            if value is None and name == "energy_eV":
                continue
            value = np.asarray(value, dtype=np.float64)
            if value.ndim != 1 or not np.isfinite(value).all():
                raise ValueError(f"{name} must be a finite one-dimensional array")
            setattr(self, name, value)
        if any(len(v) != len(self.time_ns) for v in (self.x_mm, self.y_mm, self.energy_eV)
               if v is not None):
            raise ValueError("Photon arrays are not aligned")
        if (np.any(self.time_ns < 0) or np.any(self.time_ns >= self.window_ns) or
                np.any(np.diff(self.time_ns) < 0)):
            raise ValueError("Photon times must be sorted and in [0, window_ns)")
        if self.energy_eV is not None and np.any(self.energy_eV <= 0):
            raise ValueError("Photon energies must be positive")


@dataclass
class AvalancheResult:
    records: np.ndarray
    window_ns: float
    full_charge_C: float
    last_fire_ns: np.ndarray
    metadata: dict


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_photon_hits(path, *, allow_missing_energy=False):
    """Read only aggregate photon input, never use per-primary truth to digitize."""
    path = Path(path)
    with uproot.open(path, handler=MemmapSource, num_workers=1,
                     decompression_executor=uproot.TrivialExecutor(),
                     interpretation_executor=uproot.TrivialExecutor()) as root:
        if root["WaveformEvent"].num_entries != 1 or root["RunInfo"].num_entries != 1:
            raise ValueError("Expected one independent WaveformEvent and one RunInfo")
        info = ak.to_list(root["RunInfo"].arrays(library="ak"))[0]
        tree = root["WaveformEvent"]
        names = ["all_photon_arrival_t_ps", "all_photon_arrival_x_mm",
                 "all_photon_arrival_y_mm"]
        has_energy = "photon_energy_eV" in tree.keys()
        if not has_energy and not allow_missing_energy:
            raise ValueError("No photon_energy_eV: explicitly allow the constant-PDE approximation")
        if has_energy:
            names.append("photon_energy_eV")
        arrays = tree.arrays(names, library="ak")
        values = {key: ak.to_numpy(arrays[key][0]) for key in names}
        count = int(tree["n_optical_photons_arrived_total"].array(library="np")[0])
        event_id = int(tree["event_id"].array(library="np")[0])
    if count != len(values[names[0]]):
        raise ValueError("Reported photon count differs from the arrays")
    window_ns = float(info["t_length_ps"]) / 1000
    return PhotonHits(values[names[0]].astype(np.float64) / 1000,
                      values[names[1]], values[names[2]],
                      values.get("photon_energy_eV"), window_ns,
                      float(info.get("analysis_start_ps", 0)) / 1000,
                      float(info.get("analysis_end_ps", info["t_length_ps"])) / 1000,
                      {"run_info": info, "event_id": event_id,
                       "missing_energy_approximation": not has_energy})


def map_cells(hits, config):
    """Regular rectangular pitch cells; device PDE already includes fill factor."""
    inside = ((hits.x_mm >= config.x_min_mm) & (hits.x_mm <= config.x_max_mm) &
              (hits.y_mm >= config.y_min_mm) & (hits.y_mm <= config.y_max_mm))
    cells = np.full(len(hits.time_ns), -1, dtype=np.int64)
    x = ((hits.x_mm[inside] - config.x_min_mm) /
         (config.x_max_mm - config.x_min_mm) * config.cells_x).astype(np.int64)
    y = ((hits.y_mm[inside] - config.y_min_mm) /
         (config.y_max_mm - config.y_min_mm) * config.cells_y).astype(np.int64)
    cells[inside] = np.minimum(x, config.cells_x - 1) + config.cells_x * np.minimum(y, config.cells_y - 1)
    return cells


def trigger_fraction(recovery, config):
    """Normalized empirical voltage response; not the photon PDE itself."""
    v = recovery * config.overvoltage_V - config.trigger_threshold_V
    if v <= 0:
        return 0.0
    return -math.expm1(-v / config.trigger_scale_V) / -math.expm1(
        -(config.overvoltage_V - config.trigger_threshold_V) / config.trigger_scale_V)


def _neighbours(cell, config):
    x, y = cell % config.cells_x, cell // config.cells_x
    result = []
    if x > 0:
        result.append(cell - 1)
    if x + 1 < config.cells_x:
        result.append(cell + 1)
    if y > 0:
        result.append(cell - config.cells_x)
    if y + 1 < config.cells_y:
        result.append(cell + config.cells_x)
    return result


def simulate_avalanches(hits, config, seed):
    """Merge photons, primary dark proposals and recursive secondary proposals.

    At full recovery, dark_rate_hz is the primary dark avalanche rate. Each
    avalanche proposes at most one XT and one AP, with probabilities p*q/q0.
    Their target cell accepts with trigger_fraction(recovery). All mechanisms
    share cell state. Table PDE is sampled directly once, without pre-thinning.
    """
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("Seed must be a nonnegative integer")
    if hits.energy_eV is None and config.pde_mode != "constant":
        raise ValueError("Energy-resolved PDE requires photon_energy_eV")
    cells = map_cells(hits, config)
    rng_photon, rng_transit, rng_dark, rng_secondary, rng_trigger = [
        np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(5)]
    times = hits.time_ns.copy()
    if config.transit_delay_sigma_ns:
        shape = (config.transit_delay_mean_ns / config.transit_delay_sigma_ns) ** 2
        scale = config.transit_delay_sigma_ns ** 2 / config.transit_delay_mean_ns
        times += rng_transit.gamma(shape, scale, len(times))
    else:
        times += config.transit_delay_mean_ns
    if not np.isfinite(times).all():
        raise ValueError("Nonfinite transit delays")
    # Original photon indices remain usable after transit delays reorder photons.
    order = np.argsort(times, kind="stable")
    order = order[(cells[order] >= 0) & (times[order] < hits.window_ns)]
    if config.pde_mode == "table":
        wavelengths = HC_EV_NM / hits.energy_eV
        w = np.asarray(config.pde_wavelength_nm)
        volts = np.asarray(config.pde_overvoltage_V)
        table = np.asarray(config.pde_probabilities)
        if len(order) and (np.any(wavelengths[order] < w[0]) or np.any(wavelengths[order] > w[-1])):
            raise ValueError("PDE wavelength table does not cover all photons on the sensor")
        wavelength_index = np.clip(np.searchsorted(w, wavelengths, side="right") - 1, 0, len(w) - 2)
        wavelength_fraction = (wavelengths - w[wavelength_index]) / (w[wavelength_index + 1] - w[wavelength_index])
    last = np.full(config.n_cells, -np.inf)
    queue, records = [], []
    sequence = 0
    counters = {"input_photons": len(times), "outside_footprint": int(np.sum(cells < 0)),
                "photon_candidates_after_window": int(np.sum((cells >= 0) & (times >= hits.window_ns))),
                "processed_candidates": 0, "dark_proposals": 0,
                "crosstalk_proposals": 0, "afterpulse_proposals": 0,
                "crosstalk_after_window": 0, "afterpulse_after_window": 0}

    def enqueue(t, cell, kind, parent):
        nonlocal sequence
        if not math.isfinite(t):
            raise ValueError("Nonfinite candidate time")
        if t >= hits.window_ns:
            if kind in (2, 3):
                counters[KINDS[kind] + "_after_window"] += 1
            return
        sequence += 1
        heapq.heappush(queue, (t, sequence, cell, kind, parent))

    def next_dark(t):
        if config.dark_rate_hz:
            t += rng_dark.exponential(1e9 / config.dark_rate_hz)
            if t < hits.window_ns:
                enqueue(t, int(rng_dark.integers(config.n_cells)), 1, -1)

    next_dark(0.0)
    pointer = 0
    while pointer < len(order) or queue:
        next_time = times[order[pointer]] if pointer < len(order) else math.inf
        # Queued prompt secondaries precede remaining photons with equal times.
        if queue and queue[0][0] <= next_time:
            t, _, cell, kind, parent = heapq.heappop(queue)
            photon_index = -1
            if kind == 1:
                counters["dark_proposals"] += 1
                next_dark(t)
        else:
            photon_index = int(order[pointer])
            pointer += 1
            t, cell, kind, parent = float(times[photon_index]), int(cells[photon_index]), 0, -1
        counters["processed_candidates"] += 1
        if counters["processed_candidates"] > config.max_candidates:
            raise RuntimeError("Candidate limit reached; increase max_candidates explicitly after checking the model")
        r = (-math.expm1(-(t - last[cell]) / config.recovery_ns)
             if config.recovery_enabled else 1.0)
        if kind == 0:
            if config.pde_mode == "constant":
                probability = config.constant_pde * trigger_fraction(r, config)
            else:
                voltage = r * config.overvoltage_V
                vi = min(int(np.searchsorted(volts, voltage, side="right")) - 1, len(volts) - 2)
                vi = max(vi, 0)
                vf = (voltage - volts[vi]) / (volts[vi + 1] - volts[vi])
                wi, wf = wavelength_index[photon_index], wavelength_fraction[photon_index]
                low = table[vi, wi] * (1 - wf) + table[vi, wi + 1] * wf
                high = table[vi + 1, wi] * (1 - wf) + table[vi + 1, wi + 1] * wf
                probability = low * (1 - vf) + high * vf
            accepted = rng_photon.random() < probability
        else:
            accepted = rng_trigger.random() < trigger_fraction(r, config)
        if not accepted:
            continue
        aid = len(records)
        records.append((t, cell, config.full_charge_C * r, kind, parent, photon_index, r))
        last[cell] = t
        if config.crosstalk_probability and rng_secondary.random() < config.crosstalk_probability * r:
            neighbours = _neighbours(cell, config)
            if neighbours:
                target = neighbours[int(rng_secondary.integers(len(neighbours)))]
                delay = rng_secondary.exponential(config.crosstalk_delay_ns) if config.crosstalk_delay_ns else 0.0
                counters["crosstalk_proposals"] += 1
                enqueue(t + delay, target, 2, aid)
        if config.afterpulse_probability and rng_secondary.random() < config.afterpulse_probability * r:
            counters["afterpulse_proposals"] += 1
            enqueue(t + rng_secondary.exponential(config.afterpulse_delay_ns), cell, 3, aid)
    data = np.array(records, dtype=AVALANCHE_DTYPE)
    gate = (data["time_ns"] >= hits.gate_start_ns) & (data["time_ns"] < hits.gate_end_ns)
    counters.update({f"{name}_avalanches": int(np.sum(data["kind"] == k)) for k, name in KINDS.items()})
    return AvalancheResult(data, hits.window_ns, config.full_charge_C, last,
                           {"seed": int(seed), "counters": counters,
                            "initial_state": "fully recovered, no pre-window secondary history",
                            "gate_ns": [hits.gate_start_ns, hits.gate_end_ns],
                            "gate_avalanche_count": int(gate.sum()),
                            "gate_avalanche_charge_C": float(data["charge_C"][gate].sum()),
                            "avalanche_charge_C": float(data["charge_C"].sum()),
                            "kind_codes": KINDS})


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def load_config(path):
    path = Path(path)
    raw = json.loads(path.read_text(encoding="utf-8"))
    required = {"schema_version", "sensor", "readout"}
    allowed = required | {"parameter_evidence", "configuration_notes"}
    if not required <= set(raw) <= allowed or raw["schema_version"] != 1:
        raise ValueError("Expected schema_version=1 with sensor and readout sections")
    sensor = SensorConfig(**raw["sensor"])
    readout_data = dict(raw["readout"])
    if readout_data.get("template_csv"):
        readout_data["template_csv"] = str((path.parent / readout_data["template_csv"]).resolve())
    return sensor, ReadoutConfig(**readout_data)


def digitize_root(root_path, config_path, output, seed, readout_seed, *, allow_missing_energy=False):
    """Write a new NPZ/JSON result directory with input and source provenance."""
    root_path, config_path, output = map(Path, (root_path, config_path, output))
    sensor, readout = load_config(config_path)
    if allow_missing_energy and sensor.pde_mode != "constant":
        raise ValueError("--allow-missing-energy requires constant PDE mode")
    # Fail before creating a result directory when a required input is absent.
    if not root_path.is_file():
        raise FileNotFoundError(root_path)
    if readout.template_csv and not Path(readout.template_csv).is_file():
        raise FileNotFoundError(readout.template_csv)
    output.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    project = Path(__file__).resolve().parents[2]

    def display_path(p):
        try:
            return str(p.resolve().relative_to(project))
        except ValueError:
            return str(p.resolve())

    manifest = {"status": "running", "created_utc": datetime.now(timezone.utc).isoformat(),
                "input_root": display_path(root_path), "input_bytes": root_path.stat().st_size,
                "input_sha256": sha256(root_path), "parameter_status": sensor.parameter_status,
                "sensor": asdict(sensor), "readout": asdict(readout),
                "avalanche_seed": int(seed), "readout_seed": int(readout_seed),
                "python": platform.python_version(), "numpy": np.__version__,
                "scipy": scipy.__version__, "awkward": ak.__version__,
                "uproot": uproot.__version__, "platform": platform.platform(),
                "output": display_path(output)}
    (output / "config.json").write_bytes(config_path.read_bytes())
    snapshots = output / "source"
    snapshots.mkdir()
    source_paths = [Path(__file__), Path(__file__).with_name("sipm_readout.py"),
                    Path(__file__).with_name("test_sipm_response.py")]
    manifest["source_sha256"] = {}
    for path in source_paths:
        if path.exists():
            (snapshots / path.name).write_bytes(path.read_bytes())
            manifest["source_sha256"][path.name] = sha256(path)
    if readout.template_csv:
        template = Path(readout.template_csv)
        (output / "spe_template.csv").write_bytes(template.read_bytes())
        saved = json.loads((output / "config.json").read_text())
        saved["readout"]["template_csv"] = "spe_template.csv"
        write_json(output / "config.json", saved)
        manifest["template_sha256"] = sha256(template)
    try:
        manifest["git_commit"] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=project, text=True, stderr=subprocess.DEVNULL).strip()
        (output / "git_status.txt").write_text(subprocess.check_output(
            ["git", "status", "--short"], cwd=project, text=True), encoding="utf-8")
        (output / "source_changes.patch").write_text(subprocess.check_output(
            ["git", "diff", "--", "pileup/analysis/sipm_response.py", "pileup/analysis/sipm_readout.py",
             "pileup/analysis/test_sipm_response.py"], cwd=project, text=True), encoding="utf-8")
    except (OSError, subprocess.CalledProcessError):
        manifest["git_commit"] = None
    write_json(output / "manifest.json", manifest)
    try:
        hits = read_photon_hits(root_path, allow_missing_energy=allow_missing_energy)
        avalanches = simulate_avalanches(hits, sensor, seed)
        simulation_seconds = time.perf_counter() - start
        waveform = synthesize_readout(avalanches, readout, readout_seed)
        np.savez_compressed(output / "avalanches.npz", records=avalanches.records,
                            last_fire_ns=avalanches.last_fire_ns)
        np.savez_compressed(output / "waveform.npz", **waveform.arrays)
        if sha256(root_path) != manifest["input_sha256"]:
            raise RuntimeError("Input ROOT changed during digitization")
        manifest.update({"status": "completed", "input_metadata": hits.metadata,
                         "window_ns": hits.window_ns, "avalanche_result": avalanches.metadata,
                         "readout_result": waveform.metadata, "simulation_seconds": simulation_seconds,
                         "elapsed_seconds": time.perf_counter() - start})
        manifest["output_sha256"] = {name: sha256(output / name) for name in
                                     ("avalanches.npz", "waveform.npz", "config.json")}
        write_json(output / "manifest.json", manifest)
    except Exception as error:
        manifest.update({"status": "failed", "error": f"{type(error).__name__}: {error}",
                         "elapsed_seconds": time.perf_counter() - start})
        write_json(output / "manifest.json", manifest)
        raise
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path, help="New directory; existing outputs are refused")
    parser.add_argument("--seed", required=True, type=int, help="Avalanche RNG seed")
    parser.add_argument("--readout-seed", required=True, type=int)
    parser.add_argument("--allow-missing-energy", action="store_true")
    args = parser.parse_args()
    result = digitize_root(args.root, args.config, args.output, args.seed, args.readout_seed,
                           allow_missing_energy=args.allow_missing_energy)
    print(json.dumps({"status": result["status"], "output": result["output"],
                      "parameter_status": result["parameter_status"],
                      "counters": result["avalanche_result"]["counters"],
                      "adc_samples": result["readout_result"]["samples"],
                      "saturated_samples": result["readout_result"]["saturated_samples"],
                      "elapsed_seconds": result["elapsed_seconds"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
