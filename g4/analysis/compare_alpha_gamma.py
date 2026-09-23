#!/usr/bin/env python3
"""
Run example:

  cd <project>/g4/analysis
  python3 compare_alpha_gamma.py ../data/gamma.root ../data/alpha.root

Input ROOT tree:
  PrimaryPhoton

Required branches:
  edep_total_MeV
  n_optical_photons_arrived
  photon_arrival_t_ps

Output files:
  alpha_gamma_average_waveform.png
      Area-normalized average waveforms.

  alpha_gamma_psd_2d.png
      Literature-style charge-comparison PSD plot:
      Tail integral Q_tail vs total integral Q_total, with alpha and gamma
      events combined in one 2D histogram.

  alpha_gamma_psd_ratio_2d.png
      PSD-ratio plot: Q_tail / Q_total vs Q_total.

  alpha_gamma_psd_2d_separate.png
      Old-style separated gamma/alpha panels for diagnostic checking.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import uproot
import awkward as ak


TREE_NAME = "PrimaryPhoton"

EDEP_BRANCH = "edep_total_MeV"
NPH_BRANCH = "n_optical_photons_arrived"
TIME_BRANCH = "photon_arrival_t_ps"

# -------- waveform plot settings --------
PLOT_T_MIN_NS = -200.0
PLOT_T_MAX_NS = 1200.0
WAVEFORM_BIN_WIDTH_NS = 5.0

# -------- event selection --------
MIN_PHOTONS = 20

# -------- PSD definition --------
PSD_TOTAL_MIN_NS = 0.0
PSD_TOTAL_MAX_NS = 1000.0
PSD_TAIL_MIN_NS = 200.0
PSD_TAIL_MAX_NS = 1000.0


def set_root_like_style():
    plt.rcParams["figure.figsize"] = (8, 5)
    plt.rcParams["axes.linewidth"] = 1.2
    plt.rcParams["axes.grid"] = False
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.right"] = True
    plt.rcParams["font.size"] = 12
    plt.rcParams["legend.fontsize"] = 10


def load_tree(filename):
    f = uproot.open(filename)

    if TREE_NAME not in f:
        print(f"Available keys in {filename}:")
        for key in f.keys():
            print("  ", key)
        raise RuntimeError(f"Tree '{TREE_NAME}' not found in {filename}")

    tree = f[TREE_NAME]

    required = [EDEP_BRANCH, NPH_BRANCH, TIME_BRANCH]
    missing = [br for br in required if br not in tree.keys()]
    if missing:
        print(f"Available branches in {filename}:")
        for br in tree.keys():
            print("  ", br)
        raise RuntimeError(f"Missing branches: {missing}")

    return tree


def make_waveform_bins():
    bins = np.arange(PLOT_T_MIN_NS, PLOT_T_MAX_NS + WAVEFORM_BIN_WIDTH_NS, WAVEFORM_BIN_WIDTH_NS)
    return bins


def get_aligned_time_ns(t_ps):
    """
    Convert ps to ns, and align each event to its first-arriving photon time.
    This removes the alpha/gamma time-of-flight offset, so we compare pulse shapes.
    """
    if len(t_ps) == 0:
        return np.array([], dtype=float)

    t_ns = np.asarray(ak.to_numpy(t_ps), dtype=float) / 1000.0
    t0 = np.min(t_ns)
    t_ns = t_ns - t0
    return t_ns


def accumulate_average_waveform(tree):
    """
    For each valid event:
      1. align to first photon
      2. make a histogram in [-200, 1200] ns
      3. normalize each event histogram to area = 1
      4. average over all valid events
    """
    bins = make_waveform_bins()
    shape_sum = np.zeros(len(bins) - 1, dtype=float)

    valid_events = 0
    scanned_events = 0

    branches = [EDEP_BRANCH, NPH_BRANCH, TIME_BRANCH]

    for arrays in tree.iterate(
        branches,
        step_size="50 MB",
        library="ak",
    ):
        edep = ak.to_numpy(arrays[EDEP_BRANCH])
        nph = ak.to_numpy(arrays[NPH_BRANCH])
        times = arrays[TIME_BRANCH]

        scanned_events += len(edep)

        for i in range(len(edep)):
            if edep[i] <= 0.0:
                continue
            if nph[i] < MIN_PHOTONS:
                continue
            if len(times[i]) == 0:
                continue

            t_ns = get_aligned_time_ns(times[i])

            t_plot = t_ns[(t_ns >= PLOT_T_MIN_NS) & (t_ns < PLOT_T_MAX_NS)]
            hist, _ = np.histogram(t_plot, bins=bins)

            if hist.sum() <= 0:
                continue

            shape_sum += hist / hist.sum()
            valid_events += 1

    if valid_events == 0:
        return bins, np.zeros(len(bins) - 1), 0, scanned_events

    shape_avg = shape_sum / valid_events
    return bins, shape_avg, valid_events, scanned_events


def compute_psd_arrays(tree):
    """
    Compute PSD for each valid event.

      Q_total = number of arrived photons in [0, 1000) ns
      Q_tail  = number of arrived photons in [200, 1000) ns
      PSD     = Q_tail / Q_total

    Return arrays:
      q_total_arr
      q_tail_arr
      psd_arr
    """
    q_total_list = []
    q_tail_list = []
    psd_list = []

    branches = [EDEP_BRANCH, NPH_BRANCH, TIME_BRANCH]

    for arrays in tree.iterate(
        branches,
        step_size="50 MB",
        library="ak",
    ):
        edep = ak.to_numpy(arrays[EDEP_BRANCH])
        nph = ak.to_numpy(arrays[NPH_BRANCH])
        times = arrays[TIME_BRANCH]

        for i in range(len(edep)):
            if edep[i] <= 0.0:
                continue
            if nph[i] < MIN_PHOTONS:
                continue
            if len(times[i]) == 0:
                continue

            t_ns = get_aligned_time_ns(times[i])

            q_total = np.sum((t_ns >= PSD_TOTAL_MIN_NS) & (t_ns < PSD_TOTAL_MAX_NS))
            if q_total <= 0:
                continue

            q_tail = np.sum((t_ns >= PSD_TAIL_MIN_NS) & (t_ns < PSD_TAIL_MAX_NS))
            psd = q_tail / q_total

            q_total_list.append(q_total)
            q_tail_list.append(q_tail)
            psd_list.append(psd)

    return (
        np.asarray(q_total_list, dtype=float),
        np.asarray(q_tail_list, dtype=float),
        np.asarray(psd_list, dtype=float),
    )


def percentile_axis_max(*arrays, percentile=99.5, minimum=10.0, margin=1.08):
    """Return a robust axis maximum using a percentile of finite positive values."""
    valid = []
    for arr in arrays:
        arr = np.asarray(arr, dtype=float)
        arr = arr[np.isfinite(arr) & (arr > 0.0)]
        if len(arr) > 0:
            valid.append(arr)

    if not valid:
        return minimum

    all_values = np.concatenate(valid)
    return max(float(np.percentile(all_values, percentile)) * margin, minimum)


def plot_average_waveform(gamma_bins, gamma_avg, gamma_n, alpha_bins, alpha_avg, alpha_n):
    plt.figure()

    plt.stairs(
        gamma_avg,
        gamma_bins,
        linewidth=1.5,
        label=f"gamma, N={gamma_n}",
    )

    plt.stairs(
        alpha_avg,
        alpha_bins,
        linewidth=1.5,
        label=f"alpha, N={alpha_n}",
    )

    plt.xlabel("Aligned photon arrival time [ns]")
    plt.ylabel("Average area-normalized counts / bin")
    plt.title("Alpha vs gamma average waveform")
    plt.xlim(PLOT_T_MIN_NS, PLOT_T_MAX_NS)
    plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()

    outname = "alpha_gamma_average_waveform.png"
    plt.savefig(outname, dpi=300)
    plt.close()

    print("Saved:", outname)


def plot_psd_ratio_combined(gamma_q, gamma_psd, alpha_q, alpha_psd):
    """
    Standard PSD discrimination plot:
      x = Q_total
      y = Q_tail / Q_total

    Gamma and alpha events are drawn in the same axes.
    """
    if len(gamma_q) == 0 and len(alpha_q) == 0:
        print("No valid PSD events found.")
        return

    q_max = percentile_axis_max(gamma_q, alpha_q, percentile=99.5, minimum=10.0)

    plt.figure(figsize=(8, 5.5))

    if len(gamma_q) > 0:
        plt.scatter(
            gamma_q,
            gamma_psd,
            s=10,
            alpha=0.65,
            marker=".",
            label=f"gamma, N={len(gamma_q)}",
        )

    if len(alpha_q) > 0:
        plt.scatter(
            alpha_q,
            alpha_psd,
            s=10,
            alpha=0.65,
            marker=".",
            label=f"alpha, N={len(alpha_q)}",
        )

    plt.xlabel(f"Q_total = arrived photons in {PSD_TOTAL_MIN_NS:.0f}-{PSD_TOTAL_MAX_NS:.0f} ns")
    plt.ylabel("PSD = Q_tail / Q_total")
    plt.title(
        "Alpha and gamma PSD\n"
        f"Tail window: {PSD_TAIL_MIN_NS:.0f}-{PSD_TAIL_MAX_NS:.0f} ns, "
        f"Total window: {PSD_TOTAL_MIN_NS:.0f}-{PSD_TOTAL_MAX_NS:.0f} ns"
    )
    plt.xlim(0, q_max)
    plt.ylim(0, 1.0)
    plt.legend(loc="best")
    plt.tight_layout()

    outname = "alpha_gamma_psd_ratio_2d.png"
    plt.savefig(outname, dpi=300)
    plt.close()

    print("Saved:", outname)


def plot_tail_total_combined(gamma_q, gamma_tail, alpha_q, alpha_tail):
    """
    Literature-style charge-comparison plot:
      x = Q_total
      y = Q_tail

    This is the plot style commonly used in PSD papers. Alpha and gamma
    events are combined into one 2D histogram, so diagonal bands appear when
    Q_tail is approximately proportional to Q_total.
    """
    if len(gamma_q) == 0 and len(alpha_q) == 0:
        print("No valid PSD events found.")
        return

    q_max = percentile_axis_max(gamma_q, alpha_q, percentile=99.5, minimum=10.0)
    tail_max = percentile_axis_max(gamma_tail, alpha_tail, percentile=99.5, minimum=10.0)

    all_q = []
    all_tail = []

    if len(gamma_q) > 0:
        all_q.append(gamma_q)
        all_tail.append(gamma_tail)

    if len(alpha_q) > 0:
        all_q.append(alpha_q)
        all_tail.append(alpha_tail)

    all_q = np.concatenate(all_q)
    all_tail = np.concatenate(all_tail)

    xbins = np.linspace(0.0, q_max, 180)
    ybins = np.linspace(0.0, tail_max, 180)

    plt.figure(figsize=(8, 6))

    counts, xedges, yedges, image = plt.hist2d(
        all_q,
        all_tail,
        bins=[xbins, ybins],
        cmin=1,
        norm=LogNorm(),
    )

    cb = plt.colorbar(image)
    cb.set_label("Counts")

    # Add simple labels for the two known event classes.
    # These are labels only; the histogram itself is filled with the combined data.
    if len(gamma_q) > 0:
        gx = np.median(gamma_q)
        gy = np.median(gamma_tail)
        plt.text(gx, gy, "Gammas", fontsize=12, ha="left", va="bottom")

    if len(alpha_q) > 0:
        ax = np.median(alpha_q)
        ay = np.median(alpha_tail)
        plt.text(ax, ay, "Alphas", fontsize=12, ha="left", va="bottom")

    plt.xlabel(f"Total integral Q_total [{PSD_TOTAL_MIN_NS:.0f}-{PSD_TOTAL_MAX_NS:.0f} ns photons]")
    plt.ylabel(f"Tail integral Q_tail [{PSD_TAIL_MIN_NS:.0f}-{PSD_TAIL_MAX_NS:.0f} ns photons]")
    plt.title(
        "Alpha and gamma charge-comparison PSD\n"
        f"Tail window: {PSD_TAIL_MIN_NS:.0f}-{PSD_TAIL_MAX_NS:.0f} ns, "
        f"Total window: {PSD_TOTAL_MIN_NS:.0f}-{PSD_TOTAL_MAX_NS:.0f} ns"
    )
    plt.xlim(0, q_max)
    plt.ylim(0, tail_max)
    plt.tight_layout()

    outname = "alpha_gamma_psd_2d.png"
    plt.savefig(outname, dpi=300)
    plt.close()

    print("Saved:", outname)


def plot_psd_2d_separate(gamma_q, gamma_psd, alpha_q, alpha_psd):
    """
    Old diagnostic plot: gamma and alpha are shown in separated panels.
    This is kept only for checking individual distributions.
    """
    if len(gamma_q) == 0 and len(alpha_q) == 0:
        print("No valid PSD events found.")
        return

    q_max = percentile_axis_max(gamma_q, alpha_q, percentile=99.5, minimum=10.0)

    xbins = np.linspace(0, q_max, 120)
    ybins = np.linspace(0, 1.0, 120)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)

    h1 = axes[0].hist2d(gamma_q, gamma_psd, bins=[xbins, ybins], cmin=1)
    axes[0].set_title(f"gamma PSD 2D\nN={len(gamma_q)}")
    axes[0].set_xlabel("Q_total = arrived photons in 0-1000 ns")
    axes[0].set_ylabel("PSD = Q_tail / Q_total")

    h2 = axes[1].hist2d(alpha_q, alpha_psd, bins=[xbins, ybins], cmin=1)
    axes[1].set_title(f"alpha PSD 2D\nN={len(alpha_q)}")
    axes[1].set_xlabel("Q_total = arrived photons in 0-1000 ns")

    cb1 = fig.colorbar(h1[3], ax=axes[0])
    cb1.set_label("Counts")

    cb2 = fig.colorbar(h2[3], ax=axes[1])
    cb2.set_label("Counts")

    fig.suptitle("Separated PSD 2D diagnostic\nTail window: 200-1000 ns, Total window: 0-1000 ns", y=1.02)
    fig.tight_layout()

    outname = "alpha_gamma_psd_2d_separate.png"
    fig.savefig(outname, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Saved:", outname)


def print_psd_summary(label, q_total, q_tail, psd):
    if len(q_total) == 0:
        print(f"  {label}: no valid PSD events")
        return

    print(f"  {label}:")
    print(f"    N              = {len(q_total)}")
    print(f"    Q_total mean   = {np.mean(q_total):.6g}")
    print(f"    Q_total RMS    = {np.std(q_total):.6g}")
    print(f"    Q_tail mean    = {np.mean(q_tail):.6g}")
    print(f"    Q_tail RMS     = {np.std(q_tail):.6g}")
    print(f"    PSD mean       = {np.mean(psd):.6g}")
    print(f"    PSD RMS        = {np.std(psd):.6g}")


def main():
    if len(sys.argv) < 3:
        print("Usage:")
        print("  python3 compare_alpha_gamma.py ../data/gamma.root ../data/alpha.root")
        sys.exit(1)

    gamma_file = sys.argv[1]
    alpha_file = sys.argv[2]

    set_root_like_style()

    print("Gamma file:", gamma_file)
    print("Alpha file:", alpha_file)

    gamma_tree = load_tree(gamma_file)
    alpha_tree = load_tree(alpha_file)

    print("\nTree name:", TREE_NAME)
    print("Gamma entries:", gamma_tree.num_entries)
    print("Alpha entries:", alpha_tree.num_entries)
    print("Waveform plot window:", PLOT_T_MIN_NS, "to", PLOT_T_MAX_NS, "ns")
    print("Waveform bin width:", WAVEFORM_BIN_WIDTH_NS, "ns")
    print("Minimum photons per valid event:", MIN_PHOTONS)
    print("PSD definition:")
    print(f"  Q_total: {PSD_TOTAL_MIN_NS} to {PSD_TOTAL_MAX_NS} ns")
    print(f"  Q_tail : {PSD_TAIL_MIN_NS} to {PSD_TAIL_MAX_NS} ns")
    print("  PSD = Q_tail / Q_total")

    # ---- average waveform ----
    gamma_bins, gamma_avg, gamma_valid, gamma_scanned = accumulate_average_waveform(gamma_tree)
    alpha_bins, alpha_avg, alpha_valid, alpha_scanned = accumulate_average_waveform(alpha_tree)

    print("\nAverage waveform statistics:")
    print(f"  gamma scanned events: {gamma_scanned}")
    print(f"  gamma valid events  : {gamma_valid}")
    print(f"  alpha scanned events: {alpha_scanned}")
    print(f"  alpha valid events  : {alpha_valid}")

    plot_average_waveform(gamma_bins, gamma_avg, gamma_valid, alpha_bins, alpha_avg, alpha_valid)

    # ---- PSD ----
    gamma_q, gamma_tail, gamma_psd = compute_psd_arrays(gamma_tree)
    alpha_q, alpha_tail, alpha_psd = compute_psd_arrays(alpha_tree)

    print("\nPSD statistics:")
    print_psd_summary("gamma", gamma_q, gamma_tail, gamma_psd)
    print_psd_summary("alpha", alpha_q, alpha_tail, alpha_psd)

    # Main literature-style PSD figure.
    plot_tail_total_combined(gamma_q, gamma_tail, alpha_q, alpha_tail)

    # Additional diagnostic figures.
    plot_psd_ratio_combined(gamma_q, gamma_psd, alpha_q, alpha_psd)
    plot_psd_2d_separate(gamma_q, gamma_psd, alpha_q, alpha_psd)


if __name__ == "__main__":
    main()
