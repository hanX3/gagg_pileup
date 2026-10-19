"""Run and audit native G4 SiPM single-particle windows. Work from project root.

python -m analysis.single_particle_sipm run --alpha 4 --gamma 32
python -m analysis.single_particle_sipm analyze jupyter/run_records/<campaign>
All ROOT files go to g4/data; configurations and source snapshots stay with runs.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import awkward as ak
import numpy as np
import uproot
from uproot.source.file import MemmapSource

from pileup.analysis.sipm_response import AvalancheResult, AVALANCHE_DTYPE, ELECTRON_CHARGE_C
from pileup.analysis.sipm_readout import ReadoutConfig, synthesize_readout

PROJECT = Path(__file__).resolve().parents[1]


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def command(args):
    return subprocess.run(args, cwd=PROJECT, check=True, text=True, capture_output=True).stdout.strip()


def json_library_version():
    """Record the standalone package or ROOT's bundled header version."""
    try:
        return command(["pkg-config", "--modversion", "nlohmann_json"])
    except (FileNotFoundError, subprocess.CalledProcessError):
        header = Path(command(["root-config", "--incdir"])) / "nlohmann/json.hpp"
        content = header.read_text()
        version = [re.search(rf"#define NLOHMANN_JSON_VERSION_{part}\s+(\d+)", content).group(1)
                   for part in ("MAJOR", "MINOR", "PATCH")]
        return ".".join(version) + " (ROOT-bundled header)"


def load_root(path):
    with uproot.open(path, handler=MemmapSource, num_workers=1,
                     decompression_executor=uproot.TrivialExecutor(),
                     interpretation_executor=uproot.TrivialExecutor()) as root:
        out = {}
        for name in ("RunInfo", "WaveformEvent", "PrimaryPhoton", "SiPMEvent"):
            if name not in root:
                continue
            if root[name].num_entries != 1:
                raise AssertionError(f"{path}: {name} must have one entry")
            out[name] = {key: value[0] for key, value in root[name].arrays(library="ak", how=dict).items()}
        return out


def array(row, name):
    return ak.to_numpy(row[name])


def same_tree(a, b):
    return a.keys() == b.keys() and all(ak.to_list(a[k]) == ak.to_list(b[k]) for k in a)


def run_campaign(args):
    if args.alpha < 1 or args.gamma < 1:
        raise ValueError("Require at least one alpha and one gamma")
    label = "sipm_native_" + datetime.now().strftime("%Y%m%d_%H%M%S")
    folder = PROJECT / "jupyter/run_records" / label
    folder.mkdir(parents=True, exist_ok=False)
    (folder / "macros").mkdir()
    (folder / "logs").mkdir()
    config_path = (PROJECT / args.config).resolve()
    selected_config = json.loads(config_path.read_text())
    config_relative = config_path.relative_to(PROJECT)
    sources = [PROJECT / "g4/main.cc", PROJECT / "g4/CMakeLists.txt", Path(__file__).resolve(), config_path,
               PROJECT / "g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json", PROJECT / "requirements-analysis.txt"]
    for directory, pattern in (("g4/include", "*"), ("g4/src", "*"), ("g4/macros", "*.mac"),
                               ("g4/tests", "*.cc"), ("pileup/analysis", "sipm_*.py")):
        sources.extend(p for p in (PROJECT / directory).glob(pattern) if p.is_file())
    hashes = {}
    for path in sorted(set(sources)):
        relative = path.relative_to(PROJECT)
        target = folder / "source" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        hashes[str(relative)] = digest(target)
    (folder / "working_changes.patch").write_text(command(["git", "diff", "--binary"]))
    (folder / "git_status.txt").write_text(command(["git", "status", "--short"]) + "\n")
    versions = {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "uproot", "awkward", "matplotlib")}
    versions.update({"python": sys.version, "geant4": command(["geant4-config", "--version"]),
                     "root": command(["root-config", "--version"]), "compiler": command(["c++", "--version"]),
                     "cmake": command(["cmake", "--version"]),
                     "nlohmann_json": json_library_version()})
    plan = []
    for particle, count, base in (("alpha", args.alpha, 261010200), ("gamma", args.gamma, args.gamma_seed_start-1)):
        for i in range(count):
            plan.append({"label": f"{particle}_{i:03d}", "particle": particle, "mode": "enabled",
                         "g4_seed": base+i+1, "sipm_seed": base+i+101, "readout_seed": base+i+201})
    first = plan[0]
    for mode in ("disabled", "repeat", "readout_changed"):
        item = dict(first, label="alpha_control_"+mode, mode=mode)
        if mode == "readout_changed":
            item["readout_seed"] += 100000
        plan.append(item)
    manifest = {"schema_version": 1, "created": datetime.now().astimezone().isoformat(),
                "purpose": "Fixed-count single-particle implementation validation; retain zero-deposition gammas",
                "sensor_name": selected_config["sensor"]["name"],
                "parameter_status": selected_config["sensor"]["parameter_status"],
                "git_head": command(["git", "rev-parse", "HEAD"]), "source_sha256": hashes,
                "binary_sha256": digest(PROJECT / "g4/build/gagg"), "environment": versions,
                "alpha_count": args.alpha, "gamma_count": args.gamma, "planned_runs": plan, "runs": [],
                "config": str(config_relative), "sample_note": args.sample_note,
                "initial_state": "Fully recovered cells, zero electronic history at t=0; DCR evolves for 5 us before source",
                "source": "One central +z particle from (0,0,0) mm at 5 us; alpha 5.486 MeV or gamma 661.657 keV",
                "window_ns": 20000, "gate_ns": [5000, 15000], "max_parallel_processes": 1}
    write_json(folder / "manifest.json", manifest)
    print("RUN_RECORD", folder.relative_to(PROJECT), flush=True)
    last_start = -float("inf")
    for item in plan:
        macro_text = (folder / "source/g4/macros" / f"run_single_{item['particle']}_sipm.mac").read_text()
        replacements = {"/gagg/run/seed": item["g4_seed"], "/gagg/sipm/seed": item["sipm_seed"],
                        "/gagg/sipm/readoutSeed": item["readout_seed"],
                        "/gagg/sipm/configFile": folder / "source" / config_relative,
                        "/gagg/sipm/enable": "false" if item["mode"] == "disabled" else "true"}
        macro_text = "\n".join(f"{line.split()[0]} {replacements[line.split()[0]]}" if line.split() and
                               line.split()[0] in replacements else line for line in macro_text.splitlines()) + "\n"
        macro = folder / "macros" / (item["label"] + ".mac")
        macro.write_text(macro_text)
        logfile = folder / "logs" / (item["label"] + ".log")
        time.sleep(max(0, 1.05 - (time.monotonic() - last_start)))
        last_start = time.monotonic()
        started = datetime.now().astimezone().isoformat()
        args_run = [str(PROJECT / "g4/build/gagg"), str(macro)]
        with logfile.open("w") as stream:
            process = subprocess.run(args_run, cwd=PROJECT / "g4/build", stdout=stream,
                                     stderr=subprocess.STDOUT, timeout=300)
        output = logfile.read_text()
        matches = re.findall(r"----> ROOT output: (.+)", output)
        root = Path(matches[0]).resolve() if len(matches) == 1 else None
        row = dict(item, started=started, elapsed_s=time.monotonic()-last_start, returncode=process.returncode,
                   macro=str(macro.relative_to(PROJECT)), log=str(logfile.relative_to(PROJECT)), command=args_run,
                   cwd="g4/build", macro_sha256=digest(macro),
                   root=str(root.relative_to(PROJECT)) if root else None,
                   root_sha256=digest(root) if root and root.exists() else None)
        manifest["runs"].append(row)
        write_json(folder / "manifest.json", manifest)
        if process.returncode or root is None or "ROOT file is saved" not in output or "COMMAND NOT FOUND" in output:
            raise RuntimeError(f"Run failed: {logfile}\n{output[-2500:]}")
        data = load_root(root)
        row["edep_MeV"] = float(data["WaveformEvent"]["edep_total_MeV"])
        row["n_photons"] = int(data["WaveformEvent"]["n_optical_photons_arrived_total"])
        write_json(folder / "manifest.json", manifest)
        print(item["label"], row["root"], f"edep={row['edep_MeV']:.6g} MeV, photons={row['n_photons']}", flush=True)
    return folder


def audit_response(data):
    info, optical, primary, d = (data[k] for k in ("RunInfo", "WaveformEvent", "PrimaryPhoton", "SiPMEvent"))
    config = json.loads(str(info["sipm_config_json"]))
    sensor, readout = config["sensor"], ReadoutConfig(**config["readout"])
    window = float(info["t_length_ps"]) / 1000
    t, q, kind, cell, parent, index, recovery = (array(d, k) for k in (
        "avalanche_time_ns", "avalanche_charge_C", "avalanche_kind", "avalanche_cell", "avalanche_parent",
        "avalanche_photon_index", "avalanche_recovery"))
    assert all(len(x) == len(t) for x in (q, kind, cell, parent, index, recovery))
    assert np.all(np.isfinite(t)) and np.all(np.diff(t) >= 0) and np.all((t >= 0) & (t < window))
    assert np.all((cell >= 0) & (cell < sensor["cells_x"]*sensor["cells_y"]))
    assert int(info["n_events"]) == int(optical["n_primary"]) == 1
    assert int(d["event_id"]) == int(optical["event_id"])
    assert int(d["n_input_photons"]) == len(optical["photon_energy_eV"])
    assert np.all(array(optical, "photon_energy_eV") > 0)
    assert np.all(np.diff(array(optical, "all_photon_arrival_t_ps").astype(np.int64)) >= 0)
    assert np.all((kind >= 0) & (kind <= 3))
    photon = kind == 0
    x, y = (array(optical, key) for key in ("all_photon_arrival_x_mm", "all_photon_arrival_y_mm"))
    on_sensor = ((x >= sensor["x_min_mm"]) & (x <= sensor["x_max_mm"]) &
                 (y >= sensor["y_min_mm"]) & (y <= sensor["y_max_mm"]))
    assert int(d["n_outside_sensor"]) == int((~on_sensor).sum())
    assert np.all((index[photon] >= 0) & (index[photon] < int(d["n_input_photons"])))
    assert np.all(on_sensor[index[photon]])
    assert np.all(index[~photon] == -1)
    assert len(set(index[photon])) == int(photon.sum())
    assert np.all(t[photon] >= array(optical, "all_photon_arrival_t_ps")[index[photon]]/1000)
    assert np.all(parent[kind < 2] == -1)
    assert np.all((parent[kind >= 2] >= 0) & (parent[kind >= 2] < np.flatnonzero(kind >= 2)))
    assert np.all(t[kind >= 2] >= t[parent[kind >= 2]])
    assert np.all(cell[kind == 3] == cell[parent[kind == 3]])
    xt, xt_parent = cell[kind == 2], cell[parent[kind == 2]]
    assert np.all(np.abs(xt % sensor["cells_x"] - xt_parent % sensor["cells_x"]) +
                  np.abs(xt // sensor["cells_x"] - xt_parent // sensor["cells_x"]) == 1)
    last = np.full(sensor["cells_x"]*sensor["cells_y"], -np.inf)
    for time_ns, cell_id, observed in zip(t, cell, recovery):
        expected = -np.expm1(-(time_ns-last[cell_id])/sensor["recovery_ns"]) if sensor["recovery_enabled"] else 1.
        assert np.isclose(observed, expected, atol=1e-14, rtol=1e-13)
        last[cell_id] = time_ns
    q0 = sensor["gain_electrons"]*ELECTRON_CHARGE_C
    assert np.allclose(q, q0*recovery, rtol=1e-13, atol=1e-28)
    for i, name in enumerate(("photon", "dark", "crosstalk", "afterpulse")):
        assert int(d[f"n_{name}_avalanches"]) == int((kind == i).sum())
    assert np.isclose(float(d["charge_total_C"]), q.sum(), rtol=1e-12, atol=1e-27)
    gate = (t >= float(info["analysis_start_ps"])/1000) & (t < float(info["analysis_end_ps"])/1000)
    assert np.isclose(float(d["charge_gate_C"]), q[gate].sum(), rtol=1e-12, atol=1e-27)
    balance = abs(float(d["charge_balance_error_C"])) / max(q.sum(), q0)
    assert balance < 5e-12
    records = np.zeros(len(t), dtype=AVALANCHE_DTYPE)
    records["time_ns"], records["charge_C"] = t, q
    avalanches = AvalancheResult(records, window, q0, last, {})
    offline = synthesize_readout(avalanches, replace(readout, noise_rms_mV=0), int(d["readout_seed"]))
    signal = array(d, "signal_mV")
    current_error = float(np.max(np.abs(array(d, "current_A")-offline.arrays["current_A"])))
    signal_error = float(np.max(np.abs(signal-offline.arrays["signal_mV"])))
    assert np.allclose(array(d, "current_A"), offline.arrays["current_A"], rtol=1e-10, atol=1e-17)
    assert np.allclose(signal, offline.arrays["signal_mV"], rtol=1e-10, atol=1e-12)
    assert np.array_equal(array(d, "sample_time_ns"), offline.arrays["time_ns"])
    analog = array(d, "analog_mV")
    levels = 2**readout.adc_bits
    codes = np.clip(np.floor((analog-readout.adc_min_mV)/(readout.adc_max_mV-readout.adc_min_mV)*levels), 0, levels-1)
    clipped = (analog < readout.adc_min_mV) | (analog >= readout.adc_max_mV)
    assert np.array_equal(array(d, "adc_counts"), codes)
    assert np.array_equal(array(d, "adc_clipped"), clipped.astype(int))
    assert int(d["n_clipped_samples"]) == int(clipped.sum())
    fine = synthesize_readout(avalanches, replace(readout, internal_dt_ns=readout.internal_dt_ns/2, noise_rms_mV=0), 0)
    convergence = float(np.max(np.abs(signal-fine.arrays["signal_mV"]))) / max(float(np.max(np.abs(signal))), 1e-12)
    assert convergence < .002
    return {"source_energy_MeV": float(primary["primary_energy_MeV"]), "edep_MeV": float(optical["edep_total_MeV"]),
            "n_photons": int(d["n_input_photons"]), "n_active_photons": int(on_sensor.sum()),
            "n_outside_sensor": int((~on_sensor).sum()),
            "n_avalanches": len(t), "charge_pe": float(q.sum()/q0),
            "analog_peak_mV": float(analog.max()), "clipped_samples": int(clipped.sum()),
            "charge_balance_relative": balance, "python_current_max_error_A": current_error,
            "python_signal_max_error_mV": signal_error, "half_dt_max_error_over_peak": convergence,
            "avalanche_counts": {name: int((kind == i).sum()) for i, name in enumerate(("photon", "dark", "crosstalk", "afterpulse"))}}


def analyze(folder):
    folder = Path(folder).resolve()
    manifest = json.loads((folder / "manifest.json").read_text())
    if len(manifest["runs"]) != len(manifest["planned_runs"]):
        raise ValueError("Incomplete campaign; inspect manifest and logs")
    artifacts = PROJECT / "jupyter/artifacts" / folder.name
    artifacts.mkdir(parents=True, exist_ok=True)
    results, loaded = [], {}
    for run in manifest["runs"]:
        root = PROJECT / run["root"]
        assert digest(root) == run["root_sha256"]
        data = load_root(root)
        loaded[run["label"]] = data
        if run["mode"] == "disabled":
            assert "SiPMEvent" not in data and not data["RunInfo"]["sipm_enabled"]
            continue
        results.append(dict(label=run["label"], particle=run["particle"], mode=run["mode"], root=run["root"], **audit_response(data)))
    base = loaded["alpha_000"]
    for mode in ("disabled", "repeat", "readout_changed"):
        control = loaded["alpha_control_"+mode]
        for name in ("WaveformEvent", "PrimaryPhoton"):
            assert same_tree(base[name], control[name]), f"Transport changed for {mode}"
    assert same_tree(base["SiPMEvent"], loaded["alpha_control_repeat"]["SiPMEvent"])
    changed = loaded["alpha_control_readout_changed"]["SiPMEvent"]
    for key in base["SiPMEvent"]:
        if key not in ("readout_seed", "analog_mV", "adc_counts", "adc_clipped", "n_clipped_samples"):
            assert ak.to_list(base["SiPMEvent"][key]) == ak.to_list(changed[key]), key
    assert not np.array_equal(array(base["SiPMEvent"], "adc_counts"), array(changed, "adc_counts"))
    gamma = [r for r in results if r["particle"] == "gamma" and r["mode"] == "enabled"]
    actual_sensor = json.loads(str(base["RunInfo"]["sipm_config_json"]))["sensor"]
    summary = {"parameter_status": actual_sensor["parameter_status"], "sensor_name": actual_sensor["name"],
               "run_record": str(folder.relative_to(PROJECT)), "results": results,
               "sample_note": manifest.get("sample_note", "Fixed-count incident sample"),
               "checks": {"all_native_responses_audited": True, "optical_transport_unchanged": True,
                          "native_rng_reproducible": True, "readout_seed_independent": True,
                          "python_readout_crosscheck": True, "half_dt_convergence_below_0_2_percent": True},
               "gamma_incident_count": len(gamma), "gamma_positive_edep_count": sum(r["edep_MeV"] > 0 for r in gamma),
               "gamma_with_arrival_photons": sum(r["n_photons"] > 0 for r in gamma)}
    selected = ["alpha_000"]
    if gamma and max(r["n_photons"] for r in gamma) > 0:
        selected.append(max(gamma, key=lambda r: r["n_photons"])["label"])
    summary["displayed_events"] = selected
    summary["display_selection"] = "First alpha; gamma with most recorded photons in this fixed incident sample. Illustrative conditional selection, not detection efficiency or PSD."
    os.environ.setdefault("MPLCONFIGDIR", str(PROJECT / "jupyter/artifacts/.matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for label in selected:
        data = loaded[label]
        d = data["SiPMEvent"]
        optical = data["WaveformEvent"]
        primary = data["PrimaryPhoton"]
        source_t = float(primary["primary_time_ps"])/1000
        edep = float(optical["edep_total_MeV"])
        name = f"{label}, deposited {edep:.4f} MeV"
        hist, edges = np.histogram(array(optical, "all_photon_arrival_t_ps")/1000-source_t, bins=np.arange(-50, 2002, 2))
        axes[0].stairs(hist, edges, label=name)
        axes[1].plot(array(d, "current_time_ns")-source_t, array(d, "current_A")*1e6, label=label, lw=1.2)
        cfg = json.loads(str(data["RunInfo"]["sipm_config_json"]))["readout"]
        axes[2].plot(array(d, "sample_time_ns")-source_t, array(d, "analog_mV"), label=label+" analog", lw=1)
        adc_mV = cfg["adc_min_mV"]+(array(d, "adc_counts")+.5)*(cfg["adc_max_mV"]-cfg["adc_min_mV"])/2**cfg["adc_bits"]
        axes[2].plot(array(d, "sample_time_ns")-source_t, adc_mV, label=label+" ADC bin centre", ls=":", lw=1)
    for ax, ylabel in zip(axes, ("Arriving photons / 2 ns", "SiPM current (uA)", "Readout voltage (mV)")):
        ax.set_ylabel(ylabel); ax.legend(fontsize=8); ax.grid(alpha=.2)
    axes[-1].set_xlabel("Time after source emission (ns)")
    axes[-1].set_xlim(-50, 2000)
    fig.suptitle(f"Native Geant4 SiPM digitization: {actual_sensor['name']} ({actual_sensor['parameter_status']})\nSingle incident alpha 5.486 MeV; gamma 661.657 keV (conditional display)")
    fig.tight_layout()
    fig.savefig(artifacts / "single_alpha_gamma_electric.png", dpi=160)
    plt.close(fig)
    write_json(artifacts / "summary.json", summary)
    print(json.dumps({k: v for k, v in summary.items() if k != "results"}, ensure_ascii=False, indent=2))
    print("ARTIFACTS", artifacts.relative_to(PROJECT))
    return summary


def compare_adc_ranges(before_folder, after_folder):
    """Pair stored runs by all three base seeds and compare every output branch."""
    before_folder, after_folder = Path(before_folder).resolve(), Path(after_folder).resolve()
    before = json.loads((before_folder / "manifest.json").read_text())
    after = json.loads((after_folder / "manifest.json").read_text())
    pairs = []
    for run in after["runs"]:
        if run["mode"] != "enabled":
            continue
        old = next(r for r in before["runs"] if r["mode"] == "enabled" and all(
            r[k] == run[k] for k in ("particle", "g4_seed", "sipm_seed", "readout_seed")))
        a, b = load_root(PROJECT / old["root"]), load_root(PROJECT / run["root"])
        for name in ("WaveformEvent", "PrimaryPhoton"):
            assert same_tree(a[name], b[name])
        ca = json.loads(str(a["RunInfo"]["sipm_config_json"]))
        cb = json.loads(str(b["RunInfo"]["sipm_config_json"]))
        assert ca["sensor"] == cb["sensor"]
        assert {k: v for k, v in ca["readout"].items() if k != "adc_max_mV"} == {
            k: v for k, v in cb["readout"].items() if k != "adc_max_mV"}
        changed = [k for k in a["SiPMEvent"] if ak.to_list(a["SiPMEvent"][k]) != ak.to_list(b["SiPMEvent"][k])]
        assert set(changed) <= {"adc_counts", "adc_clipped", "n_clipped_samples"}
        pairs.append({"old_root": old["root"], "new_root": run["root"],
                      "adc_max_mV": [ca["readout"]["adc_max_mV"], cb["readout"]["adc_max_mV"]],
                      "old_clipped": int(a["SiPMEvent"]["n_clipped_samples"]),
                      "new_clipped": int(b["SiPMEvent"]["n_clipped_samples"]),
                      "changed_branches": changed, "optical_and_analog_identical": True})
    if not pairs:
        raise ValueError("No paired enabled runs")
    target = PROJECT / "jupyter/artifacts" / after_folder.name / "adc_range_check.json"
    write_json(target, {"comparisons": pairs})
    return pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run")
    run.add_argument("--alpha", type=int, default=4)
    run.add_argument("--gamma", type=int, default=32)
    run.add_argument("--config", type=Path, default=Path("g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json"))
    run.add_argument("--gamma-seed-start", type=int, default=261010501)
    run.add_argument("--sample-note", default="Fixed-count incident sample; not conditioned on energy deposition")
    check = sub.add_parser("analyze")
    check.add_argument("folder", type=Path)
    ranges = sub.add_parser("compare-ranges")
    ranges.add_argument("before", type=Path)
    ranges.add_argument("after", type=Path)
    args = parser.parse_args()
    if args.command == "run":
        analyze(run_campaign(args))
    elif args.command == "compare-ranges":
        print(json.dumps(compare_adc_ranges(args.before, args.after), indent=2))
    else:
        analyze(args.folder)


if __name__ == "__main__":
    main()
