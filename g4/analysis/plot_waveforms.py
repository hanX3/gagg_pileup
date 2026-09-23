#!/usr/bin/env python3
"""
Run examples:

  cd <project>/g4/analysis

  python3 plot_waveforms.py ../data/gamma.root gamma
  python3 plot_waveforms.py ../data/alpha.root alpha

Input ROOT tree:

  PrimaryPhoton

Required branches:

  edep_total_MeV
  n_optical_photons_arrived
  photon_arrival_t_ps

Output files:

  gamma_raw_first10_waveforms.png
  gamma_normalized_first10_waveforms.png
  gamma_normalized_average_all_events.png

  alpha_raw_first10_waveforms.png
  alpha_normalized_first10_waveforms.png
  alpha_normalized_average_all_events.png
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import uproot
import awkward as ak


TREE_NAME = "PrimaryPhoton"

EDEP_BRANCH = "edep_total_MeV"
NPH_BRANCH = "n_optical_photons_arrived"
TIME_BRANCH = "photon_arrival_t_ps"

# Plot time window
T_MIN_NS = -200.0
T_MAX_NS = 1200.0
BIN_WIDTH_NS = 2.0

# Event selection
MIN_PHOTONS = 20
N_RAW_EVENTS = 10


def set_root_like_style():
    """Set a simple ROOT-like plotting style."""
    plt.rcParams["figure.figsize"] = (8, 5)
    plt.rcParams["axes.linewidth"] = 1.2
    plt.rcParams["axes.grid"] = False
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.right"] = True
    plt.rcParams["font.size"] = 12
    plt.rcParams["legend.fontsize"] = 8


def load_tree(filename):
    f = uproot.open(filename)

    if TREE_NAME not in f:
        print("Available keys:")
        for key in f.keys():
            print("  ", key)
        raise RuntimeError(f"Tree '{TREE_NAME}' not found in {filename}")

    tree = f[TREE_NAME]

    required = [EDEP_BRANCH, NPH_BRANCH, TIME_BRANCH]
    missing = [br for br in required if br not in tree.keys()]
    if missing:
        print("Available branches:")
        for br in tree.keys():
            print("  ", br)
        raise RuntimeError(f"Missing branches: {missing}")

    return tree


def make_time_bins():
    bins = np.arange(T_MIN_NS, T_MAX_NS + BIN_WIDTH_NS, BIN_WIDTH_NS)
    centers = 0.5 * (bins[:-1] + bins[1:])
    return bins, centers


def make_hist_from_ps(t_ps, bins):
    """
    Convert ps to ns and make a histogram without time alignment.

    Since bins start at -200 ns, the region [-200, 0) will naturally be zero
    for most events, which is what you want for plotting.
    """
    t_ns = np.asarray(ak.to_numpy(t_ps), dtype=float) / 1000.0
    t_ns = t_ns[(t_ns >= T_MIN_NS) & (t_ns < T_MAX_NS)]
    hist, _ = np.histogram(t_ns, bins=bins)
    return hist


def get_first_valid_events(tree, n_events):
    """
    Find the first n valid events.
    A valid event has:
      edep_total_MeV > 0
      n_optical_photons_arrived >= MIN_PHOTONS
      non-empty photon_arrival_t_ps
    """
    valid_events = []

    branches = [EDEP_BRANCH, NPH_BRANCH, TIME_BRANCH]
    global_entry = 0

    for arrays in tree.iterate(
        branches,
        step_size="50 MB",
        library="ak",
    ):
        edep = ak.to_numpy(arrays[EDEP_BRANCH])
        nph = ak.to_numpy(arrays[NPH_BRANCH])
        times = arrays[TIME_BRANCH]

        for i in range(len(edep)):
            entry = global_entry + i

            if edep[i] <= 0.0:
                continue
            if nph[i] < MIN_PHOTONS:
                continue
            if len(times[i]) == 0:
                continue

            valid_events.append(
                {
                    "entry": entry,
                    "edep": float(edep[i]),
                    "nph": int(nph[i]),
                    "times": times[i],
                }
            )

            if len(valid_events) >= n_events:
                return valid_events

        global_entry += len(edep)

    return valid_events


def plot_raw_first10(label, valid_events, bins):
    """
    Plot the first 10 valid raw waveforms.
    No time alignment.
    No normalization.
    """
    plt.figure()

    for ev in valid_events:
        hist = make_hist_from_ps(ev["times"], bins)

        plt.stairs(
            hist,
            bins,
            linewidth=1.0,
            label=f"Entry {ev['entry']}, Edep={ev['edep']:.4f} MeV, Nph={ev['nph']}",
        )

    plt.xlabel("Photon arrival time [ns]")
    plt.ylabel("Counts / bin")
    plt.title(f"{label}: raw waveforms of first {len(valid_events)} valid events")
    plt.xlim(T_MIN_NS, T_MAX_NS)
    plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()

    outname = f"{label}_raw_first10_waveforms.png"
    plt.savefig(outname, dpi=300)
    plt.close()

    print("Saved:", outname)


def plot_normalized_first10(label, valid_events, bins):
    """
    Plot the first 10 valid normalized waveforms.
    No time alignment.
    Each event is normalized by its own total counts in the time window.
    """
    plt.figure()

    for ev in valid_events:
        hist = make_hist_from_ps(ev["times"], bins)

        if hist.sum() <= 0:
            continue

        hist_norm = hist / hist.sum()

        plt.stairs(
            hist_norm,
            bins,
            linewidth=1.0,
            label=f"Entry {ev['entry']}, Edep={ev['edep']:.4f} MeV",
        )

    plt.xlabel("Photon arrival time [ns]")
    plt.ylabel("Area-normalized counts / bin")
    plt.title(f"{label}: normalized waveforms of first {len(valid_events)} valid events")
    plt.xlim(T_MIN_NS, T_MAX_NS)
    plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()

    outname = f"{label}_normalized_first10_waveforms.png"
    plt.savefig(outname, dpi=300)
    plt.close()

    print("Saved:", outname)


def plot_normalized_average_all(label, tree, bins):
    """
    Plot the normalized average waveform over all valid events.

    For each valid event:
      1. make histogram
      2. normalize histogram area to 1
      3. add to shape_sum

    Final waveform:
      shape_avg = shape_sum / number_of_valid_events
    """
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

            hist = make_hist_from_ps(times[i], bins)

            if hist.sum() <= 0:
                continue

            shape_sum += hist / hist.sum()
            valid_events += 1

    print("Scanned events:", scanned_events)
    print("Valid events used for normalized average:", valid_events)

    if valid_events == 0:
        print("No valid events found. No average waveform was produced.")
        return

    shape_avg = shape_sum / valid_events

    plt.figure()
    plt.stairs(shape_avg, bins, linewidth=1.5)

    plt.xlabel("Photon arrival time [ns]")
    plt.ylabel("Average area-normalized counts / bin")
    plt.title(f"{label}: normalized average waveform, all valid events\nN={valid_events}")
    plt.xlim(T_MIN_NS, T_MAX_NS)
    plt.ylim(bottom=0)
    plt.tight_layout()

    outname = f"{label}_normalized_average_all_events.png"
    plt.savefig(outname, dpi=300)
    plt.close()

    print("Saved:", outname)


def main():
    if len(sys.argv) < 3:
        print("Usage:")
        print("  python3 plot_waveform.py ../data/gamma.root gamma")
        print("  python3 plot_waveform.py ../data/alpha.root alpha")
        sys.exit(1)

    filename = sys.argv[1]
    label = sys.argv[2]

    set_root_like_style()

    print("Input ROOT file:", filename)
    print("Label:", label)

    tree = load_tree(filename)

    print("Tree:", TREE_NAME)
    print("Entries:", tree.num_entries)
    print("Time window:", T_MIN_NS, "to", T_MAX_NS, "ns")
    print("Bin width:", BIN_WIDTH_NS, "ns")
    print("Minimum photons per valid event:", MIN_PHOTONS)

    bins, _ = make_time_bins()

    valid_events = get_first_valid_events(tree, N_RAW_EVENTS)

    print(f"First {len(valid_events)} valid events:")
    for ev in valid_events:
        print(
            f"  Entry {ev['entry']}: "
            f"Edep={ev['edep']:.6f} MeV, "
            f"Nph={ev['nph']}"
        )

    if len(valid_events) == 0:
        print("No valid events found. Check edep or photon collection.")
        return

    plot_raw_first10(label, valid_events, bins)
    plot_normalized_first10(label, valid_events, bins)
    plot_normalized_average_all(label, tree, bins)


if __name__ == "__main__":
    main()
