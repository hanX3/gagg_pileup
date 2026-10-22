"""Render the two dual-decay comparisons in Notebook 0002 from saved fits.

This is a presentation step: it reads frozen fit parameters and prediction
arrays, and never runs an optimizer or modifies a historical artifact.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DUAL_MODELS = ("two_fixed", "two_free")
PARTICLES = ("alpha", "gamma")


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def render_dual_decay_figures(summary_path, arrays_path, output_directory):
    """Redraw ADC fits/residuals and sensitivity using only the two saved models.

    The returned manifest identifies every input array used in the ADC figure.
    Area fractions and residual metrics in the second figure are read unchanged
    from the original report. Existing scientific results remain untouched.
    """
    summary_path, arrays_path = Path(summary_path), Path(arrays_path)
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    report = json.loads(summary_path.read_text())
    results = report["results"]
    colors = {"alpha": "#d7655a", "gamma": "#4d83bd"}
    styles = {
        "two_fixed": ("#202020", "-", "Two decays: fixed lifetimes"),
        "two_free": ("#d3922a", "--", "Two decays: free lifetimes"),
    }
    array_hashes = {}

    with np.load(arrays_path, allow_pickle=False) as arrays, plt.rc_context(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    ):
        fig, axes = plt.subplots(
            2, 2, figsize=(12, 7), sharex="col",
            gridspec_kw={"height_ratios": [3, 1]}, layout="constrained",
        )
        for column, particle in enumerate(PARTICLES):
            fits = results[particle]["voltage"]["adc"]
            reference = fits["two_fixed"]["reference"]
            names = [f"{particle}_time_ns", f"{particle}_adc_mV"]
            names += [f"{particle}_{model}_prediction_mV" for model in DUAL_MODELS]
            for name in names:
                values = arrays[name]
                if not np.all(np.isfinite(values)):
                    raise ValueError(f"Non-finite saved array: {name}")
                array_hashes[name] = hashlib.sha256(values.tobytes()).hexdigest()
            time = arrays[f"{particle}_time_ns"] - reference["anchor_ns"]
            observed = arrays[f"{particle}_adc_mV"]
            keep = (time >= -100) & (time <= 2000)
            axes[0, column].plot(
                time[keep], observed[keep] - reference["baseline_mV"],
                color=colors[particle], lw=.8, alpha=.8, label="Observed ADC",
            )
            for model in DUAL_MODELS:
                color, linestyle, label = styles[model]
                prediction = arrays[f"{particle}_{model}_prediction_mV"]
                if prediction.shape != observed.shape or time.shape != observed.shape:
                    raise ValueError(f"Incompatible saved arrays: {particle}/{model}")
                axes[0, column].plot(
                    time[keep], prediction[keep] - reference["baseline_mV"],
                    color=color, ls=linestyle, lw=1.4, label=label,
                )
                axes[1, column].plot(
                    time[keep], observed[keep] - prediction[keep],
                    color=color, ls=linestyle, lw=.8, alpha=.85,
                )
            axes[0, column].set(
                title=f"({chr(97 + column)}) {particle}: ADC and dual-decay fits",
                ylabel="Voltage above baseline (mV)",
            )
            axes[0, column].legend(fontsize=8)
            axes[1, column].axhline(0, color=".4", lw=.8)
            axes[1, column].set(
                title=f"({chr(99 + column)}) {particle}: residuals",
                xlabel="Time from observed threshold crossing (ns)",
                ylabel="Residual (mV)", xlim=(-100, 2000),
            )
        fig.suptitle("Two scintillation decays; full readout response and 0.5 ns aperture")
        fig.savefig(output / "adc_fits_residuals.png", dpi=170)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), layout="constrained")
        for particle in PARTICLES:
            result = results[particle]
            windows = sorted(
                result["adc_window_scan"] + [result["voltage"]["adc"]["two_fixed"]],
                key=lambda fit: fit["stop_ns"],
            )
            axes[0].plot(
                [fit["stop_ns"] for fit in windows],
                [fit["params"]["slow_fraction"] for fit in windows],
                "o-", color=colors[particle], label=f"{particle} ADC fit",
            )
            axes[0].axhline(
                result["optical_fixed"]["slow_fraction"], color=colors[particle],
                ls="--", alpha=.7, label=f"{particle} photon fit",
            )
            axes[1].plot(
                [0, 1],
                [result["voltage"]["adc"][model]["nrmse_percent_peak"]
                 for model in DUAL_MODELS],
                "o-", color=colors[particle], label=particle,
            )
        axes[0].set(
            title="(a) Long-window stability: fixed lifetimes",
            xlabel="ADC fit end from observed crossing (ns)",
            ylabel="Fitted slow area fraction", ylim=(0, 1),
        )
        axes[1].set(
            title="(b) Dual-decay fit residuals at 2000 ns",
            xticks=[0, 1], xticklabels=["Fixed lifetimes", "Free lifetimes"],
            xlim=(-.2, 1.2), ylabel="Residual RMS / peak (%)",
        )
        for axis in axes:
            axis.legend(fontsize=8)
            axis.grid(alpha=.2)
        fig.suptitle("Saved dual-decay results; descriptive residuals, no confidence interval")
        fig.savefig(output / "fit_sensitivity.png", dpi=170)
        plt.close(fig)

    manifest = {
        "operation": "Render frozen dual-decay fits only; no refitting",
        "models": list(DUAL_MODELS),
        "summary_sha256": _digest(summary_path),
        "arrays_sha256": _digest(arrays_path),
        "plotted_array_sha256": array_hashes,
        "renderer_sha256": _digest(Path(__file__)),
        "figures": {name: _digest(output / name) for name in (
            "adc_fits_residuals.png", "fit_sensitivity.png",
        )},
    }
    (output / "provenance.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    )
    return manifest
