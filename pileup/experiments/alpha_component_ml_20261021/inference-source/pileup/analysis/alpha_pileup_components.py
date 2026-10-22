"""Plot actual single-alpha contributions and previously frozen blind fits.

Run from the project root. This is postprocessing: it never fits a waveform,
selects a pulse count, runs G4, or modifies source data. Ancestry is used only to
display realized shared-cell contributions and to retain the original scoring.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform

import numpy as np

from .alpha_full_signal_study import make_template
from .alpha_pileup_fit import score
from .sipm_readout import ReadoutConfig, synthesize_readout
from .sipm_response import AvalancheResult, AVALANCHE_DTYPE, ELECTRON_CHARGE_C


PROJECT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = PROJECT / "jupyter/artifacts/alpha_correlated_global_20261020"
METHODS = ("local", "diagonal")
LABELS = {"local": "Local-K SSE", "diagonal": "Diagonal WLS"}
COLORS = {"local": "#2878b5", "diagonal": "#c43832"}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def relative(path):
    try:
        return str(Path(path).resolve().relative_to(PROJECT))
    except ValueError:
        return str(Path(path).resolve())


def saved_model(source, stem, method):
    path = (source / "test/windows" / f"{stem}_{method}.npz" if method == "local"
            else source / "test" / method / f"{stem}.npz")
    with np.load(path) as data:
        return data["model_mV"].copy(), path


def create(output, *, source=DEFAULT_SOURCE, rate_hz=1e7, start_ns=5000., stop_ns=7000.,
           overview_id="19-6", individual_id="19-7", layout="comparison", method="diagonal"):
    source, output = Path(source).resolve(), Path(output).resolve()
    if layout not in ("comparison", "triptych") or method not in METHODS:
        raise ValueError("Unknown component plot layout or physical method")
    if not np.isfinite([rate_hz, start_ns, stop_ns]).all() or rate_hz <= 0 or stop_ns <= start_ns:
        raise ValueError("Require a positive rate and finite increasing display bounds")
    if output.exists():
        raise FileExistsError(f"Use a new output directory: {output}")
    summary_path = source / "test-summary.json"
    summary = json.loads(summary_path.read_text())
    row = min((r for r in summary["windows"] if r["rate_hz"] == rate_hz), key=lambda r: r["seed"])
    if not row["gate_ns"][0] <= start_ns < stop_ns <= row["gate_ns"][1]:
        raise ValueError("The display interval must lie inside the existing scoring gate")
    stem = Path(row["root"]).stem
    waveform_path, truth_path = PROJECT / row["waveform"], PROJECT / row["truth"]
    for path, expected in ((waveform_path, row["waveform_sha256"]), (truth_path, row["truth_sha256"])):
        if digest(path) != expected:
            raise ValueError(f"Saved input identity mismatch: {path}")
    metadata_path = source / "inputs/metadata" / f"{stem}.json"
    metadata = json.loads(metadata_path.read_text())
    settings_path = source / "baseline-settings.json"
    settings = json.loads(settings_path.read_text())
    template = make_template(settings)
    baseline = settings["fit_config"]["baseline_mV"]
    with np.load(waveform_path) as data:
        t, observed = data["time_ns"].copy(), data["analog_mV"].copy()
    with np.load(truth_path) as data:
        truth = {key: data[key].copy() for key in data.files}
    if not t[0] <= start_ns < stop_ns < t[-1]:
        raise ValueError("Display interval outside the saved waveform")
    actual_time, actual_charge = truth["first_active_ns"], truth["charge_pe"]
    central = ((truth["emission_ns"] >= row["gate_ns"][0]) &
               (truth["emission_ns"] < row["gate_ns"][1]))
    target_indices = np.flatnonzero(central & np.isfinite(actual_time) & (actual_charge >= 500) &
                                   (actual_time >= start_ns) & (actual_time < stop_ns))
    target_indices = target_indices[np.argsort(actual_time[target_indices])]
    if not len(target_indices):
        raise ValueError("The fixed display interval has no threshold-qualified targets")
    names = {int(i): f"T{j+1}" for j, i in enumerate(target_indices)}

    # Readout is linear in the *realized* avalanche charge. No sensor simulation
    # is repeated, no microcells are reset, and no fresh electronic noise is added.
    readout = replace(ReadoutConfig(**metadata["config"]["readout"]), baseline_mV=0., noise_rms_mV=0.)
    q0 = metadata["config"]["sensor"]["gain_electrons"] * ELECTRON_CHARGE_C
    avalanche = np.zeros(len(truth["avalanche_time_ns"]), dtype=AVALANCHE_DTYPE)
    avalanche["time_ns"] = truth["avalanche_time_ns"]
    avalanche["charge_C"] = truth["avalanche_charge_C"]
    owner = truth["avalanche_owner"]

    def voltage(records):
        result = synthesize_readout(
            AvalancheResult(records, metadata["duration_ns"], q0, np.array([]), {}), readout, 0)
        np.testing.assert_array_equal(result.arrays["time_ns"], t)
        return result.arrays["signal_mV"]

    native_rebuilt = voltage(avalanche)
    actual_components = np.stack([voltage(avalanche[owner == i]) for i in target_indices])
    other_actual = voltage(avalanche[~np.isin(owner, target_indices)])
    native_error = float(np.max(abs(native_rebuilt - truth["signal_mV"])))
    additivity_error = float(np.max(abs(actual_components.sum(axis=0) + other_actual - native_rebuilt)))
    if native_error > 1e-8 or additivity_error > 1e-8:
        raise ValueError("Attributed voltage does not reproduce the saved native signal")
    charge_errors = [float(avalanche["charge_C"][owner == i].sum() / q0 - actual_charge[i])
                     for i in target_indices]
    if max(map(abs, charge_errors), default=0.) > 1e-7:
        raise ValueError("Attributed avalanche charge does not match the scoring truth")

    arrays = dict(time_ns=t, observed_analog_mV=observed, native_signal_mV=truth["signal_mV"],
                  target_indices=target_indices, actual_components_mV=actual_components,
                  other_actual_components_mV=other_actual)
    event_rows = [dict(label=names[int(i)], truth_index=int(i), primary_id=int(truth["primary_id"][i]),
                       first_active_ns=float(actual_time[i]), true_charge_pe=float(actual_charge[i]), methods={})
                  for i in target_indices]
    method_rows, input_files = {}, [summary_path, waveform_path, truth_path, metadata_path, settings_path]
    for method in METHODS:
        saved = row["methods"][method]
        fit, old_score = saved["fit"], saved["score"]
        rescored = score(fit, truth, row["gate_ns"])
        if json.loads(json.dumps(rescored)) != old_score:
            raise ValueError(f"Original scoring did not reproduce for {method}")
        model, model_path = saved_model(source, stem, method)
        input_files.append(model_path)
        full = fit["all_components"]
        full_model = baseline + template(t[:, None] - np.asarray(full["times_ns"])[None, :]) @ np.asarray(full["peak_mV"])
        model_error = float(np.max(abs(full_model - model)))
        if model_error > 1e-8:
            raise ValueError(f"All-component reconstruction does not match the saved {method} model")
        matched = dict(old_score["all_matched_pairs"])
        used = set(matched.values())
        times, peaks = np.asarray(fit["times_ns"]), np.asarray(fit["peak_mV"])
        displayed = np.flatnonzero((times >= start_ns) & (times < stop_ns))
        reported = template(t[:, None] - times[displayed][None, :]).T * peaks[displayed, None]
        matched_components = np.zeros_like(actual_components)
        fitted_times, fitted_charges = np.full(len(target_indices), np.nan), np.full(len(target_indices), np.nan)
        for k, (i, event) in enumerate(zip(target_indices, event_rows)):
            if int(i) not in matched:
                event["methods"][method] = dict(status="missed", prediction_index=None)
                continue
            j = matched[int(i)]
            dt, dq = float(times[j] - actual_time[i]), float(fit["charge_pe"][j] / actual_charge[i] - 1.)
            status = "complete" if abs(dq) <= .2 else "charge_failed"
            event["methods"][method] = dict(status=status, prediction_index=j,
                fitted_time_ns=float(times[j]), fitted_charge_pe=float(fit["charge_pe"][j]),
                time_error_ns=dt, relative_charge_error=dq)
            matched_components[k] = template(t - times[j]) * peaks[j]
            fitted_times[k], fitted_charges[k] = times[j], fit["charge_pe"][j]
        extra = [int(j) for j in displayed if j not in used]
        method_rows[method] = dict(
            target_count=len(target_indices),
            complete=sum(e["methods"][method]["status"] == "complete" for e in event_rows),
            charge_failed=sum(e["methods"][method]["status"] == "charge_failed" for e in event_rows),
            missed=sum(e["methods"][method]["status"] == "missed" for e in event_rows),
            extra_output_indices=extra, displayed_prediction_indices=displayed.tolist(),
            displayed_prediction_truth_labels={str(j): names.get(next((i for i, jj in matched.items() if jj == j), -1))
                                               for j in displayed},
            all_component_model_max_error_mV=model_error)
        arrays.update({f"{method}_model_mV": model,
            f"{method}_displayed_components_mV": reported,
            f"{method}_other_components_mV": model - baseline - reported.sum(axis=0),
            f"{method}_matched_components_mV": matched_components,
            f"{method}_matched_time_ns": fitted_times, f"{method}_matched_charge_pe": fitted_charges})

    output.mkdir(parents=True)
    np.savez_compressed(output / "components.npz", **arrays)
    manifest = dict(created_at=datetime.now().astimezone().isoformat(), source=relative(source),
        root=row["root"], root_sha256=row["root_sha256"], source_seed=row["seed"], rate_hz=rate_hz,
        display_interval_ns=[start_ns, stop_ns], original_fit_interval_ns=[0., metadata["duration_ns"]],
        original_scoring_gate_ns=row["gate_ns"], source_freeze=summary["frozen_at"],
        selection="Lowest predeclared source seed at the requested rate, fixed display crop; no selection by success",
        truth_use="Post-inference ancestry decomposition and unchanged scoring only; no inference rerun",
        component_definition="Actual photon-descendant avalanche charges; native shared-cell state retained; linear readout only, no added noise/baseline",
        figure_ids=[overview_id] if layout == "triptych" else [overview_id, individual_id],
        layout=layout, displayed_method=method if layout == "triptych" else None,
        methods=method_rows, events=event_rows,
        crop_is_not_a_refit=True, individual_display_relative_ns=[-30., 800.],
        caveats=["Fixed crop is an illustration, not a new efficiency sample",
                 "Full charge values do not come from integrating a truncated displayed curve",
                 "Other components include older tails, weak signals and dark chains; no tail is discarded",
                 "Matched colors and target labels use truth only after the original blind inference"],
        validation=dict(native_signal_max_error_mV=native_error, voltage_additivity_max_error_mV=additivity_error,
                        ancestry_charge_max_error_pe=max(map(abs, charge_errors), default=0.), original_scores_reproduced=True),
        inputs={relative(p): digest(p) for p in input_files},
        source_sha256={relative(Path(__file__)): digest(__file__),
                       "pileup/analysis/sipm_readout.py": digest(PROJECT / "pileup/analysis/sipm_readout.py"),
                       "pileup/analysis/alpha_pileup_fit.py": digest(PROJECT / "pileup/analysis/alpha_pileup_fit.py")},
        environment=dict(python=platform.python_version(), packages={p: importlib.metadata.version(p)
                         for p in ("numpy", "scipy", "matplotlib")}))
    if layout == "triptych":
        render_triptych(output, arrays, manifest, method)
    else:
        render(output, arrays, manifest, baseline)
    manifest["outputs"] = {p.name: digest(p) for p in sorted(output.iterdir()) if p.is_file()}
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    return manifest


def render(output, arrays, manifest, baseline):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    output = Path(output)
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .18, "savefig.dpi": 160})
    t, y = arrays["time_ns"], arrays["observed_analog_mV"]
    lo, hi = manifest["display_interval_ns"]
    show = (t >= lo) & (t < hi)
    events = manifest["events"]
    pulse_colors = {e["label"]: plt.get_cmap("tab10")(i % 10) for i, e in enumerate(events)}
    fig, axes = plt.subplots(5, 1, figsize=(14, 14), sharex=True, layout="constrained",
                             gridspec_kw={"height_ratios": [1.25, 1., 1., 1., .65]})
    axes[0].plot(t[show] / 1000, y[show], color=".25", lw=1., label="Observed pile-up voltage")
    for method in METHODS:
        axes[0].plot(t[show] / 1000, arrays[f"{method}_model_mV"][show],
                     color=COLORS[method], lw=1.2, ls="--" if method == "local" else "-", label=LABELS[method])
    axes[0].set_title("(a) One observed pile-up waveform and the saved complete fits", loc="left")
    axes[0].set_ylabel("Total voltage (mV)")
    axes[0].legend(ncol=3, loc="upper right", fontsize=9)
    axes[1].plot(t[show] / 1000, arrays["other_actual_components_mV"][show], color=".65", lw=1.,
                 label="Other alpha tails + weak signals + dark chains")
    for component, event in zip(arrays["actual_components_mV"], events):
        color = pulse_colors[event["label"]]
        axes[1].plot(t[show] / 1000, component[show], color=color, lw=1.2, label=event["label"])
        axes[1].axvline(event["first_active_ns"] / 1000, color=color, lw=.6, alpha=.5, ls=":")
    axes[1].set_title("(b) Actual single-alpha contributions (truth after inference; no baseline)", loc="left")
    axes[1].set_ylabel("Component (mV)")
    axes[1].legend(ncol=6, fontsize=8, loc="upper right")
    for ax, method, panel in zip(axes[2:4], METHODS, ("c", "d")):
        state = manifest["methods"][method]
        ax.plot(t[show] / 1000, arrays[f"{method}_other_components_mV"][show], color=".65", lw=1.,
                label="Earlier/later tails + subthreshold fit terms")
        for component, j in zip(arrays[f"{method}_displayed_components_mV"], state["displayed_prediction_indices"]):
            label = state["displayed_prediction_truth_labels"][str(j)]
            extra = j in state["extra_output_indices"]
            color = "#d60084" if extra else pulse_colors.get(label, ".45")
            label = f"Extra output {j}" if extra else (label or "Matched outside crop")
            ax.plot(t[show] / 1000, component[show], color=color, lw=1.15, label=label)
        ax.set_title(f"({panel}) {LABELS[method]}: fitted single-pulse curves; "
                     f"{state['complete']}/{state['target_count']} complete, "
                     f"{state['charge_failed']} charge failure, {state['missed']} missed, "
                     f"{len(state['extra_output_indices'])} extra", loc="left")
        ax.set_ylabel("Component (mV)")
        ax.legend(ncol=6, fontsize=8, loc="upper right")
    for method in METHODS:
        axes[4].plot(t[show] / 1000, (y - arrays[f"{method}_model_mV"])[show],
                     color=COLORS[method], lw=.9, label=LABELS[method])
    axes[4].axhline(0., color=".4", ls=":", lw=.8)
    axes[4].set(title="(e) Observed minus total fit", ylabel="Residual (mV)", xlabel="Record time (us)")
    axes[4].legend(ncol=2, fontsize=8)
    for ax in axes:
        ax.set_xlim(lo / 1000, hi / 1000)
    fig.suptitle(f"Figure 0003-{manifest['figure_ids'][0]}: {manifest['rate_hz']/1e6:g} MHz pile-up and pulse decomposition\n"
                 "Fixed crop of the first predeclared window; fitting used the complete 20 us record", fontsize=13)
    fig.savefig(output / "pileup_decomposition.png")
    plt.close(fig)

    columns, rows = 3, int(np.ceil(len(events) / 3))
    fig, axes = plt.subplots(rows, columns, figsize=(14, 3.5 * rows + .8), squeeze=False)
    for k, (ax, event) in enumerate(zip(axes.ravel(), events)):
        relative_t = t - event["first_active_ns"]
        show = (relative_t >= -30.) & (relative_t <= 800.)
        ax.plot(relative_t[show], arrays["actual_components_mV"][k, show], color=".2", lw=1.3)
        labels = []
        bad = False
        for method in METHODS:
            result = event["methods"][method]
            tag = "LK" if method == "local" else "WLS"
            if result["status"] == "missed":
                labels.append(f"{tag}: missed")
                bad = True
                continue
            ax.plot(relative_t[show], arrays[f"{method}_matched_components_mV"][k, show],
                    color=COLORS[method], lw=1.1, ls="--" if method == "local" else "-")
            passed = result["status"] == "complete"
            bad |= not passed
            labels.append(f"{tag}: dQ {100*result['relative_charge_error']:+.1f}%, "
                          f"dt {result['time_error_ns']:+.1f} ns " + ("[pass]" if passed else "[Q fail]"))
        ax.text(.98, .96, "\n".join(labels), transform=ax.transAxes, ha="right", va="top", fontsize=8,
                bbox=dict(facecolor="white", edgecolor=".8", alpha=.9))
        ax.set_title(f"({chr(97+k)}) {event['label']}: true Q={event['true_charge_pe']:.0f} pe\n"
                     f"First active photon at {event['first_active_ns']/1000:.6f} us",
                     fontsize=10, color="#a03020" if bad else ".15")
        ax.set(xlim=(-30, 800), xlabel="Time relative to this truth onset (ns)", ylabel="Single-pulse voltage (mV)")
        ax.axvline(0., color=".5", lw=.6, ls=":")
    for ax in axes.ravel()[len(events):]:
        ax.set_visible(False)
    fig.suptitle(f"Figure 0003-{manifest['figure_ids'][1]}: Each threshold-qualified alpha in the fixed crop\n"
                 "Actual shared-cell contribution vs matched blind-fit components; no amplitude normalization", fontsize=13, y=.995)
    handles = [Line2D([], [], color=".2", label="Actual attributed pulse"),
               Line2D([], [], color=COLORS["local"], ls="--", label=LABELS["local"]),
               Line2D([], [], color=COLORS["diagonal"], label=LABELS["diagonal"])]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .945), ncol=3, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, .905), h_pad=2.)
    fig.savefig(output / "individual_alpha_overlays.png")
    plt.close(fig)


def render_triptych(output, arrays, manifest, method="diagonal"):
    """Figure-6-2 style: total, paired component overlays, and residual."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .17, "savefig.dpi": 160})
    t, y = arrays["time_ns"], arrays["observed_analog_mV"]
    lo, hi = manifest["display_interval_ns"]
    shown = (t >= lo) & (t < hi)
    model = arrays[f"{method}_model_mV"]
    events = manifest["events"]
    state = manifest["methods"][method]
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True, layout="constrained")
    axes[0].plot(t[shown], y[shown], color=".3", lw=1., label="Observed voltage")
    axes[0].plot(t[shown], model[shown], color="#b33e52", lw=1.25,
                 label=f"Blind total fit ({LABELS[method]})")
    axes[0].set_title("(a) Observed pile-up waveform and blind total fit", loc="left", fontsize=10)
    axes[0].set_ylabel("Total (mV)")
    axes[0].legend(loc="upper right", fontsize=8)

    palette = plt.get_cmap("tab20")
    pulse_handles = []
    for k, event in enumerate(events):
        # Truth supplies plotting colors and matched pairs only after inference.
        color = palette((2 * k if k < 10 else 2 * (k - 10) + 1) % 20)
        actual = arrays["actual_components_mV"][k]
        axes[1].plot(t[shown], actual[shown], color=color, lw=1.15, ls="--")
        if event["methods"][method]["status"] != "missed":
            axes[1].plot(t[shown], arrays[f"{method}_matched_components_mV"][k, shown],
                         color=color, lw=1.1)
        pulse_handles.append(Line2D([], [], color=color, lw=1.1,
                                    label=f"alpha {event['truth_index']}"))
    extra = set(state["extra_output_indices"])
    for curve, j in zip(arrays[f"{method}_displayed_components_mV"], state["displayed_prediction_indices"]):
        if j in extra:
            axes[1].plot(t[shown], curve[shown], color="#cf0084", ls=":", lw=1.5)
    axes[1].set_title("(b) Same color: actual pulse (dashed) and fitted pulse (solid)",
                      loc="left", fontsize=10)
    axes[1].set_ylabel("Components (mV)")
    handles = [Line2D([], [], color=".25", ls="--", label="Actual attributed pulse"),
               Line2D([], [], color=".25", label="Matched fitted pulse")]
    if extra:
        handles.append(Line2D([], [], color="#cf0084", ls=":", label="Extra fitted output"))
    # Compact color labels identify individual pulses without splitting panels.
    handles += pulse_handles
    legend_rows = int(np.ceil(len(handles) / 5))
    axes[1].legend(handles=handles, ncol=5, loc="upper right", fontsize=6.7,
                   handlelength=1.8, columnspacing=.9)
    axes[1].set_ylim(top=axes[1].get_ylim()[1] * (1.12 + .05 * legend_rows))

    axes[2].plot(t[shown], (y - model)[shown], color=".3", lw=.9)
    axes[2].axhline(0., color=".5", ls=":", lw=.9)
    axes[2].set(title="(c) Observed minus total fit", ylabel="Residual (mV)", xlabel="Record time (ns)")
    for ax in axes:
        ax.set_xlim(lo, hi)
    fig.suptitle(f"Figure 0003-{manifest['figure_ids'][0]}: {manifest['rate_hz']/1e6:g} MHz; "
                 f"{LABELS[method]} pulse recovery", fontsize=12)
    fig.savefig(Path(output) / "pulse_recovery_triptych.png")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rate-mhz", type=float, default=10.)
    parser.add_argument("--start-us", type=float, default=5.)
    parser.add_argument("--stop-us", type=float, default=7.)
    parser.add_argument("--overview-id", default="19-6")
    parser.add_argument("--individual-id", default="19-7")
    parser.add_argument("--layout", choices=("comparison", "triptych"), default="comparison")
    parser.add_argument("--method", choices=METHODS, default="diagonal")
    args = parser.parse_args()
    result = create(args.output, source=args.source, rate_hz=args.rate_mhz * 1e6,
                    start_ns=args.start_us * 1000, stop_ns=args.stop_us * 1000,
                    overview_id=args.overview_id, individual_id=args.individual_id,
                    layout=args.layout, method=args.method)
    print(json.dumps(dict(root=result["root"], events=len(result["events"]),
                         methods=result["methods"], validation=result["validation"])))


if __name__ == "__main__":
    main()
