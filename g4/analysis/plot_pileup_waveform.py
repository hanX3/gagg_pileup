#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_pileup_waveform.py

Read a pileup ROOT file and plot ONLY the first event in WaveformEvent.

The figure includes:
  - total waveform from WaveformEvent: black
  - individual primary waveforms from PrimaryPhoton:
      alpha: red
      x-ray / bremsstrahlung: olive
      12C capture / de-excitation gamma: blue
      other / unknown: gray

Command-line usage
------------------
python3 plot_pileup_waveform.py ../data/gagg_waveform_20260624_15h46m34s.root

Change histogram bin width:
python3 plot_pileup_waveform.py input.root --bin-width-ns 2.0

The PNG is saved directly in the current working directory. For example:

  input ROOT:
    ../data/gagg_waveform_20260624_15h46m34s.root

  output PNG:
    ./20260624_15h46m34s.png

Jupyter usage
-------------
Run the notebook from the project root directory, then use:

    from analysis.plot_pileup_waveform import plot_first_event

    fig = plot_first_event(
        "data/gagg_waveform_20260624_15h46m34s.root",
        bin_width_ns=1.0,
        save=False,
    )
    fig

Set save=True to save a PNG as well. Use output_file to specify its path.
"""

import argparse
from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import uproot


def to_str(x):
    if x is None:
        return ""
    return str(x)


def classify_source(primary_particle_name, primary_response_name):
    """
    Return the plotting color and legend label for one primary source.
    """
    pn = to_str(primary_particle_name).strip()
    rn = to_str(primary_response_name).strip()

    # Alpha particles.
    if pn and pn.lower() != "gamma":
        return "red", "alpha"

    # X-ray / bremsstrahlung component.
    if rn == "gamma_brems":
        return "olive", "x-ray"

    # 12C capture / de-excitation gamma component.
    if rn.startswith("gamma_c12_capture"):
        return "blue", "12C de-excitation gamma"

    # Fallback for other primary types.
    return "gray", rn if rn else (pn if pn else "other")


def make_output_filename(root_file_name):
    """
    Build the output PNG filename from the ROOT filename.

    Example:
      gagg_waveform_20260624_15h46m34s.root
        -> 20260624_15h46m34s.png

      test.root
        -> test.png
    """
    stem = Path(root_file_name).stem

    prefix = "gagg_waveform_"
    if stem.startswith(prefix):
        stem = stem[len(prefix):]

    return Path.cwd() / f"{stem}.png"


def make_histogram_from_ps(times_ps, edges_ns):
    """
    Convert photon-arrival times from ps to a binned waveform in ns.
    """
    if times_ps is None:
        return np.zeros(len(edges_ns) - 1, dtype=np.int64)

    times_ps = np.asarray(times_ps, dtype=np.float64)
    if times_ps.size == 0:
        return np.zeros(len(edges_ns) - 1, dtype=np.int64)

    times_ns = times_ps / 1000.0
    counts, _ = np.histogram(times_ns, bins=edges_ns)
    return counts


def get_runinfo_window_ns(root_file):
    """
    Try to read the waveform time-window length from RunInfo::t_length_ps.
    Return None if RunInfo or t_length_ps is not available.
    """
    try:
        if "RunInfo" not in root_file:
            return None

        run_info = root_file["RunInfo"]
        arrays = run_info.arrays(library="ak")

        if "t_length_ps" in arrays.fields and len(arrays["t_length_ps"]) > 0:
            return float(arrays["t_length_ps"][0]) / 1000.0
    except Exception:
        return None

    return None


def build_edges_ns(total_times_ps, bin_width_ns, forced_tmax_ns=None):
    """
    Build histogram bin edges in ns.

    If RunInfo::t_length_ps exists, use the full event window.
    Otherwise, use the maximum total photon-arrival time in the first event.
    """
    if forced_tmax_ns is not None:
        tmax_ns = forced_tmax_ns
    else:
        if total_times_ps is None or len(total_times_ps) == 0:
            tmax_ns = bin_width_ns
        else:
            tmax_ns = float(
                np.max(np.asarray(total_times_ps, dtype=np.float64)) / 1000.0
            )

        tmax_ns = max(tmax_ns, bin_width_ns)

    n_bins = max(int(np.ceil(tmax_ns / bin_width_ns)), 1)
    return np.linspace(0.0, n_bins * bin_width_ns, n_bins + 1)


def plot_first_event(
    root_file,
    bin_width_ns=1.0,
    output_file=None,
    save=False,
):
    """
    Plot the first WaveformEvent entry and return the Matplotlib figure.

    Parameters
    ----------
    root_file : str or pathlib.Path
        Input pileup ROOT file.
    bin_width_ns : float, optional
        Histogram bin width in ns. Default: 1.0.
    output_file : str or pathlib.Path or None, optional
        PNG output path. If omitted and save=True, the filename is generated
        from the ROOT filename in the current working directory.
    save : bool, optional
        Save the figure to a PNG file. Default: False, which is convenient
        for displaying the figure directly in Jupyter.

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure. In Jupyter, place ``fig`` on the last line of
        a cell to display it.
    """
    root_file = Path(root_file)

    if bin_width_ns <= 0:
        raise ValueError("bin_width_ns must be greater than zero.")
    if not root_file.exists():
        raise FileNotFoundError(f"Input ROOT file not found: {root_file}")

    with uproot.open(root_file) as root:
        if "WaveformEvent" not in root:
            raise RuntimeError("Tree 'WaveformEvent' not found.")
        if "PrimaryPhoton" not in root:
            raise RuntimeError("Tree 'PrimaryPhoton' not found.")

        waveform = root["WaveformEvent"].arrays(library="ak")
        primary = root["PrimaryPhoton"].arrays(library="ak")

        if len(waveform["event_id"]) == 0:
            raise RuntimeError("WaveformEvent tree is empty.")

        # Only the first entry of WaveformEvent is plotted.
        widx = 0
        event_id = int(waveform["event_id"][widx])

        total_times_ps = waveform["all_photon_arrival_t_ps"][widx]
        runinfo_tmax_ns = get_runinfo_window_ns(root)
        edges_ns = build_edges_ns(
            total_times_ps,
            bin_width_ns,
            forced_tmax_ns=runinfo_tmax_ns,
        )
        x_ns = edges_ns[:-1]
        total_counts = make_histogram_from_ps(total_times_ps, edges_ns)

        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot the individual primary waveforms that belong to the first event.
        p_event_ids = ak.to_numpy(primary["event_id"])
        pidxs = np.where(p_event_ids == event_id)[0]

        used_labels = set()

        for pidx in pidxs:
            times_ps = primary["photon_arrival_t_ps"][pidx]
            counts = make_histogram_from_ps(times_ps, edges_ns)

            particle_name = (
                primary["primary_particle_name"][pidx]
                if "primary_particle_name" in primary.fields
                else ""
            )
            response_name = (
                primary["primary_response_name"][pidx]
                if "primary_response_name" in primary.fields
                else ""
            )

            color, label = classify_source(particle_name, response_name)
            plot_label = label if label not in used_labels else None
            used_labels.add(label)

            ax.step(
                x_ns,
                counts,
                where="post",
                color=color,
                linewidth=1.2,
                alpha=0.75,
                label=plot_label,
            )

        # Plot the total waveform last so that it stays visible on top.
        ax.step(
            x_ns,
            total_counts,
            where="post",
            color="black",
            linewidth=2.0,
            label="total waveform",
            zorder=10,
        )

        ax.set_xlabel("Time [ns]")
        ax.set_ylabel("Photon counts / bin")
        ax.set_title(f"Pileup waveform, first event: event_id = {event_id}")
        ax.grid(True, alpha=0.25)
        ax.legend()
        fig.tight_layout()

    if save:
        outfile = (
            Path(output_file)
            if output_file is not None
            else make_output_filename(root_file)
        )
        outfile.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outfile, dpi=200)
        print(f"[saved] {outfile}")

    print(f"[info] plotted first WaveformEvent entry, event_id = {event_id}")
    print(f"[info] number of primary entries in this event = {len(pidxs)}")

    return fig


def main():
    parser = argparse.ArgumentParser(
        description="Plot only the first event from a pileup ROOT file."
    )
    parser.add_argument("root_file", help="Input ROOT file")
    parser.add_argument(
        "--bin-width-ns",
        type=float,
        default=1.0,
        help="Histogram bin width in ns. Default: 1.0",
    )
    args = parser.parse_args()

    fig = plot_first_event(
        args.root_file,
        bin_width_ns=args.bin_width_ns,
        save=True,
    )

    # The command-line program does not need to keep the figure open.
    plt.close(fig)


if __name__ == "__main__":
    main()
