"""Complete-record pulse overlays from saved native signals and frozen fits.

All nonzero-charge alpha contributions, dark chains, reported fit components and
internal subthreshold fit components are retained. No inference or G4 is run.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime
import importlib.metadata
import json
from pathlib import Path
import platform
import time

import numpy as np

from .alpha_pileup_components import PROJECT, DEFAULT_SOURCE, LABELS, digest, relative, saved_model
from .alpha_full_signal_study import make_template
from .alpha_pileup_fit import score
from .sipm_readout import ReadoutConfig, synthesize_readout
from .sipm_response import AvalancheResult, ELECTRON_CHARGE_C


def create(output, *, rate_hz, figure_id, source=DEFAULT_SOURCE, method="diagonal"):
    started = time.monotonic()
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"Use a new output directory: {output}")
    if method not in ("local", "diagonal"):
        raise ValueError("Choose a previously saved local or diagonal physical fit")
    summary_path = source / "test-summary.json"
    summary = json.loads(summary_path.read_text())
    row = min((r for r in summary["windows"] if r["rate_hz"] == rate_hz), key=lambda r: r["seed"])
    waveform_path, truth_path = PROJECT / row["waveform"], PROJECT / row["truth"]
    for p, sha in ((waveform_path, row["waveform_sha256"]), (truth_path, row["truth_sha256"])):
        if digest(p) != sha:
            raise ValueError(f"Changed source input: {p}")
    stem = Path(row["root"]).stem
    metadata_path = source / "inputs/metadata" / f"{stem}.json"
    metadata = json.loads(metadata_path.read_text())
    if metadata["root_sha256"] != row["root_sha256"] or metadata["seed"] != row["seed"]:
        raise ValueError("The metadata and saved waveform refer to different ROOTs")
    settings_path = source / "baseline-settings.json"
    settings = json.loads(settings_path.read_text())
    template = make_template(settings)
    baseline = float(settings["fit_config"]["baseline_mV"])
    with np.load(waveform_path) as data:
        t, observed = data["time_ns"].copy(), data["analog_mV"].copy()
    with np.load(truth_path) as data:
        truth = {k: data[k].copy() for k in data.files}
    if set(truth["particle_name"]) - {"alpha"}:
        raise ValueError("This plotter is for the pure-alpha comparison")
    fitted = row["methods"][method]
    if json.loads(json.dumps(score(fitted["fit"], truth, row["gate_ns"]))) != fitted["score"]:
        raise ValueError("The original scoring did not reproduce")
    model, model_path = saved_model(source, stem, method)

    duration = float(metadata["duration_ns"])
    q0 = metadata["config"]["sensor"]["gain_electrons"] * ELECTRON_CHARGE_C
    readout = replace(ReadoutConfig(**metadata["config"]["readout"]), baseline_mV=0., noise_rms_mV=0.)
    owners = truth["avalanche_owner"]
    if np.any(owners < -1) or np.any(owners >= len(truth["charge_pe"])):
        raise ValueError("Invalid avalanche ancestry index")
    # Stable grouping avoids rescanning millions of avalanches for every alpha.
    # The readout only consumes these two fields, retaining original times/charges.
    order = np.argsort(owners, kind="stable")
    records = np.empty(len(order), dtype=[("time_ns", "f8"), ("charge_C", "f8")])
    records["time_ns"] = truth["avalanche_time_ns"][order]
    records["charge_C"] = truth["avalanche_charge_C"][order]
    counts = np.bincount(owners + 1, minlength=len(truth["charge_pe"]) + 1)
    bounds = np.concatenate(([0], np.cumsum(counts)))
    del order

    def voltage(block):
        result = synthesize_readout(AvalancheResult(block, duration, q0, np.array([]), {}), readout, 0)
        np.testing.assert_array_equal(result.arrays["time_ns"], t)
        return result.arrays["signal_mV"]

    indices = np.flatnonzero(truth["charge_pe"] > 0.)
    actual = np.empty((len(indices), len(t)), dtype=float)
    charge_error = 0.
    for k, i in enumerate(indices):
        block = records[bounds[i+1]:bounds[i+2]]
        actual[k] = voltage(block)
        charge_error = max(charge_error, abs(float(block["charge_C"].sum()/q0 - truth["charge_pe"][i])))
    dark = voltage(records[bounds[0]:bounds[1]])
    rebuilt = voltage(records)
    native_error = float(np.max(abs(rebuilt - truth["signal_mV"])))
    sum_error = float(np.max(abs(actual.sum(axis=0) + dark - truth["signal_mV"])))
    if native_error > 1e-8 or sum_error > 1e-8 or charge_error > 1e-7:
        raise ValueError("Full ancestry decomposition failed native voltage/charge checks")

    fit = fitted["fit"]
    full = fit.get("all_components") or fit
    fit_times = np.asarray(full.get("times_ns", []), float)
    fit_peaks = np.asarray(full.get("peak_mV", []), float)
    fit_charges = np.asarray(full.get("charge_pe", []), float)
    components = template(t[None, :] - fit_times[:, None]) * fit_peaks[:, None]
    model_error = float(np.max(abs(baseline + components.sum(axis=0) - model)))
    if model_error > 1e-8:
        raise ValueError("Full fit components do not reproduce the saved total fit")
    report_lookup = {(float(tt), float(aa)): j for j, (tt, aa) in enumerate(zip(fit["times_ns"], fit.get("peak_mV", [])))}
    report_indices = np.array([report_lookup.get((float(tt), float(aa)), -1)
                               for tt, aa in zip(fit_times, fit_peaks)], dtype=int)
    if set(report_indices[report_indices >= 0]) != set(range(len(fit["times_ns"]))):
        raise ValueError("Reported fit components are not an exact subset of the internal model")
    pairs = fitted["score"]["all_matched_pairs"]
    report_truth = {j: i for i, j in pairs}
    fit_truth_indices = np.array([report_truth.get(int(j), -1) for j in report_indices], dtype=int)
    threshold = 500.
    events = [dict(truth_index=int(i), primary_id=int(truth["primary_id"][i]),
                   first_active_ns=float(truth["first_active_ns"][i]),
                   emission_ns=float(truth["emission_ns"][i]), charge_pe=float(truth["charge_pe"][i]),
                   above_threshold=bool(truth["charge_pe"][i] >= threshold)) for i in indices]
    arrays = dict(time_ns=t, observed_analog_mV=observed, model_mV=model,
        native_signal_mV=truth["signal_mV"], actual_components_mV=actual, actual_truth_indices=indices,
        native_dark_mV=dark, fitted_components_mV=components, fitted_times_ns=fit_times,
        fitted_charge_pe=fit_charges, fitted_reported_indices=report_indices, fitted_truth_indices=fit_truth_indices)
    output.mkdir(parents=True)
    np.savez_compressed(output / "components.npz", **arrays)
    paths = (summary_path, waveform_path, truth_path, metadata_path, settings_path, model_path)
    manifest = dict(created_at=datetime.now().astimezone().isoformat(), source=relative(source),
        root=row["root"], root_sha256=row["root_sha256"], source_seed=row["seed"], rate_hz=rate_hz,
        figure_id=figure_id, method=method, display_interval_ns=[0., duration],
        fit_interval_ns=[0., duration], scoring_interval_ns=row["gate_ns"], baseline_mV=baseline,
        sample_count=len(t), source_freeze=summary["frozen_at"],
        selection="Lowest predeclared source seed at the rate; every recorded sample shown; no success-based selection",
        component_scope="Every alpha with realized Q>0, full dark-chain voltage and every internal fitted component, including subthreshold terms and guards",
        truth_use="Post-inference plotting and original scoring only; no inference rerun or retuning",
        native_nonzero_alpha_count=len(events), native_threshold_alpha_count=sum(e["above_threshold"] for e in events),
        fitted_reported_count=len(fit["times_ns"]), fitted_internal_count=len(fit_times),
        full_unmatched_reported_count=int(np.sum((report_indices >= 0) & (fit_truth_indices < 0))),
        events=events, score=fitted["score"],
        validation=dict(native_signal_max_error_mV=native_error, all_component_sum_max_error_mV=sum_error,
                        ancestry_charge_max_error_pe=charge_error, fit_model_max_error_mV=model_error,
                        original_score_reproduced=True, omitted_nonzero_alphas=0, omitted_fit_components=0),
        inputs={relative(p): digest(p) for p in paths},
        source_sha256={relative(Path(__file__)): digest(__file__),
                       "pileup/analysis/alpha_pileup_components.py": digest(PROJECT / "pileup/analysis/alpha_pileup_components.py"),
                       "pileup/analysis/sipm_readout.py": digest(PROJECT / "pileup/analysis/sipm_readout.py")},
        environment=dict(python=platform.python_version(), packages={p: importlib.metadata.version(p)
                         for p in ("numpy", "scipy", "matplotlib")}))
    render(output, arrays, manifest)
    manifest["elapsed_seconds"] = time.monotonic() - started
    manifest["outputs"] = {p.name: digest(p) for p in sorted(output.iterdir()) if p.is_file()}
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    return manifest


def render(output, arrays, manifest):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .14,
                         "savefig.dpi": 180, "path.simplify": True})
    t = arrays["time_ns"] / 1000.
    lo, hi = np.asarray(manifest["display_interval_ns"]) / 1000.
    gate_lo, gate_hi = np.asarray(manifest["scoring_interval_ns"]) / 1000.
    fig, axes = plt.subplots(3, 1, figsize=(17, 10), sharex=True, layout="constrained",
                             gridspec_kw={"height_ratios": [1., 1.3, .65]})
    axes[0].plot(t, arrays["observed_analog_mV"], color=".3", lw=.7, label="Observed voltage")
    axes[0].plot(t, arrays["model_mV"], color="#b33e52", lw=.95, label=f"Blind total fit ({LABELS[manifest['method']]})")
    axes[0].set(title="(a) Complete observed record and saved total fit", ylabel="Total (mV)")
    axes[0].legend(handles=[*axes[0].get_legend_handles_labels()[0], Patch(facecolor=".88", label="Guard regions")],
                   ncol=3, loc="upper right", fontsize=8)

    def color(index):
        return plt.get_cmap("hsv")((int(index) * .618033988749895) % 1.)

    # Every truth curve is drawn. Small signals and internal nuisance terms are
    # gray, distinct from threshold-qualified paired alpha components.
    for k, event in enumerate(manifest["events"]):
        c = color(event["truth_index"]) if event["above_threshold"] else ".72"
        axes[1].plot(t, arrays["actual_components_mV"][k], color=c, ls="--",
                     lw=.65 if event["above_threshold"] else .4, alpha=.85 if event["above_threshold"] else .6)
    axes[1].plot(t, arrays["native_dark_mV"], color=".45", ls="--", lw=.4, alpha=.6)
    for k, (report, truth_index) in enumerate(zip(arrays["fitted_reported_indices"], arrays["fitted_truth_indices"])):
        c, style, width = (".65", ":", .45) if report < 0 else (
            (color(truth_index), "-", .65) if truth_index >= 0 else ("#cf0084", ":", .8))
        axes[1].plot(t, arrays["fitted_components_mV"][k], color=c, ls=style, lw=width, alpha=.8)
    handles = [Line2D([], [], color=".2", ls="--", label="Actual alpha, Q >= 500 pe"),
               Line2D([], [], color=".2", label="Matched fit: same color"),
               Line2D([], [], color="#cf0084", ls=":", label="Unmatched reported fit"),
               Line2D([], [], color=".72", ls="--", label="Weak actual alpha / dark"),
               Line2D([], [], color=".65", ls=":", label="Internal subthreshold fit")]
    axes[1].legend(handles=handles, ncol=3, loc="upper right", fontsize=8)
    axes[1].set(title=f"(b) All actual and fitted components: {manifest['native_nonzero_alpha_count']} nonzero alphas, "
                      f"{manifest['fitted_internal_count']} fit terms; dashed truth / solid matched fit",
                ylabel="Components (mV)")
    low_y, high_y = axes[1].get_ylim()
    axes[1].set_ylim(low_y, high_y * 1.18 if high_y > 0 else 1.)
    axes[2].plot(t, arrays["observed_analog_mV"] - arrays["model_mV"], color=".3", lw=.7)
    axes[2].axhline(0., color=".45", ls=":", lw=.8)
    axes[2].set(title="(c) Observed minus saved total fit", ylabel="Residual (mV)", xlabel="Record time (us)")
    for ax in axes:
        ax.axvspan(lo, gate_lo, color=".9", alpha=.45, zorder=-10)
        ax.axvspan(gate_hi, hi, color=".9", alpha=.45, zorder=-10)
        ax.axvline(gate_lo, color=".6", ls=":", lw=.7)
        ax.axvline(gate_hi, color=".6", ls=":", lw=.7)
        ax.set_xlim(lo, hi)
    fig.suptitle(f"Figure 0003-{manifest['figure_id']}: {manifest['rate_hz']/1e6:g} MHz; full 0-{hi:g} us record\n"
                 "Every recorded sample and all alpha / fit components shown; central 5-15 us used for scoring", fontsize=13)
    fig.savefig(Path(output) / "pulse_recovery_full.png")
    fig.savefig(Path(output) / "pulse_recovery_full.svg")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--rate-mhz", type=float, required=True)
    parser.add_argument("--figure-id", required=True)
    parser.add_argument("--method", choices=("local", "diagonal"), default="diagonal")
    args = parser.parse_args()
    m = create(args.output, rate_hz=args.rate_mhz * 1e6, figure_id=args.figure_id, source=args.source, method=args.method)
    print(json.dumps(dict(rate_mhz=args.rate_mhz, root=m["root"], alphas=m["native_nonzero_alpha_count"],
                         fit_terms=m["fitted_internal_count"], seconds=m["elapsed_seconds"], validation=m["validation"])))


if __name__ == "__main__":
    main()
