#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_p11b_alpha_spectrum.py

Generate a first-version effective single-alpha energy spectrum for
p + 11B -> 3 alpha + 8.7 MeV.

This is a practical source-spectrum model for detector-response / pileup
studies, not a full three-body Dalitz generator.

Model:
  1. Primary alpha component:
       E_primary = 3.76 MeV

  2. Secondary alpha component from sequential 8Be* decay:
       E_secondary = E_boost + E_star
                     + 2*sqrt(E_boost*E_star)*cos(theta)

     cos(theta) is uniformly sampled in [-1, 1].

Default values:
  E_primary = 3.76 MeV
  E_star    = 1.515 MeV
  E_boost   = 0.94 MeV
  primary fraction = 1/3
  secondary fraction = 2/3

How to run:
  python3 make_p11b_alpha_spectrum.py

Recommended:
  python3 make_p11b_alpha_spectrum.py --n-alpha 1000000 --out-prefix p11b_alpha

With primary peak broadening:
  python3 make_p11b_alpha_spectrum.py --primary-sigma-mev 0.05

Outputs:
  <out-prefix>_histogram.csv
  <out-prefix>_spectrum_linear.png
  <out-prefix>_spectrum_logy.png

Optional:
  <out-prefix>_samples.csv, if --save-samples is used.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def sample_p11b_alpha_spectrum(
    n_alpha,
    primary_energy_mev=3.76,
    primary_fraction=1.0 / 3.0,
    secondary_estar_mev=1.515,
    secondary_eboost_mev=0.94,
    primary_sigma_mev=0.0,
    seed=12345,
):
    rng = np.random.default_rng(seed)

    is_primary = rng.random(n_alpha) < primary_fraction
    energies = np.empty(n_alpha, dtype=float)
    labels = np.empty(n_alpha, dtype=object)

    n_primary = np.count_nonzero(is_primary)
    if primary_sigma_mev > 0:
        primary_e = rng.normal(primary_energy_mev, primary_sigma_mev, n_primary)
        primary_e = np.clip(primary_e, 0.0, None)
    else:
        primary_e = np.full(n_primary, primary_energy_mev)

    energies[is_primary] = primary_e
    labels[is_primary] = "primary"

    n_secondary = n_alpha - n_primary
    cos_theta = rng.uniform(-1.0, 1.0, n_secondary)

    e_star = secondary_estar_mev
    e_boost = secondary_eboost_mev
    secondary_e = e_boost + e_star + 2.0 * np.sqrt(e_boost * e_star) * cos_theta

    energies[~is_primary] = secondary_e
    labels[~is_primary] = "secondary"

    return energies, labels


def build_histogram(energies, labels, bins, hist_min, hist_max):
    edges = np.linspace(hist_min, hist_max, bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    all_counts, _ = np.histogram(energies, bins=edges)
    primary_counts, _ = np.histogram(energies[labels == "primary"], bins=edges)
    secondary_counts, _ = np.histogram(energies[labels == "secondary"], bins=edges)

    df = pd.DataFrame({
        "E_low_MeV": edges[:-1],
        "E_high_MeV": edges[1:],
        "E_center_MeV": centers,
        "counts_total": all_counts,
        "counts_primary": primary_counts,
        "counts_secondary": secondary_counts,
    })

    total = all_counts.sum()
    df["prob_total_per_bin"] = all_counts / total if total > 0 else 0.0
    return df


def plot_spectrum(df, out_prefix):
    x = df["E_center_MeV"].to_numpy()
    y_total = df["counts_total"].to_numpy()
    y_primary = df["counts_primary"].to_numpy()
    y_secondary = df["counts_secondary"].to_numpy()

    plt.figure(figsize=(8, 5))
    plt.step(x, y_total, where="mid", label="total")
    plt.step(x, y_primary, where="mid", label="primary alpha")
    plt.step(x, y_secondary, where="mid", label="secondary alpha")
    plt.xlabel("Alpha energy [MeV]")
    plt.ylabel("Counts / bin")
    plt.title("Effective p-11B alpha energy spectrum")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_spectrum_linear.png", dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.step(x, y_total, where="mid", label="total")
    plt.step(x, y_primary, where="mid", label="primary alpha")
    plt.step(x, y_secondary, where="mid", label="secondary alpha")
    plt.xlabel("Alpha energy [MeV]")
    plt.ylabel("Counts / bin")
    plt.yscale("log")
    plt.title("Effective p-11B alpha energy spectrum")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_spectrum_logy.png", dpi=200)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Generate a first-version effective alpha energy spectrum for p + 11B -> 3 alpha."
    )
    parser.add_argument("--n-alpha", type=int, default=1_000_000)
    parser.add_argument("--primary-energy-mev", type=float, default=3.76)
    parser.add_argument("--primary-fraction", type=float, default=1.0 / 3.0)
    parser.add_argument("--secondary-estar-mev", type=float, default=1.515)
    parser.add_argument("--secondary-eboost-mev", type=float, default=0.94)
    parser.add_argument("--primary-sigma-mev", type=float, default=0.0)
    parser.add_argument("--bins", type=int, default=600)
    parser.add_argument("--emin-mev", type=float, default=0.0)
    parser.add_argument("--emax-mev", type=float, default=6.0)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--out-prefix", type=str, default="p11b_alpha")
    parser.add_argument("--save-samples", action="store_true")
    args = parser.parse_args()

    if not (0.0 <= args.primary_fraction <= 1.0):
        raise ValueError("--primary-fraction must be between 0 and 1.")

    e_star = args.secondary_estar_mev
    e_boost = args.secondary_eboost_mev
    sec_emin = (np.sqrt(e_star) - np.sqrt(e_boost)) ** 2
    sec_emax = (np.sqrt(e_star) + np.sqrt(e_boost)) ** 2

    print("=== p-11B effective alpha spectrum model ===")
    print(f"N alpha                  = {args.n_alpha}")
    print(f"Primary fraction         = {args.primary_fraction:.6f}")
    print(f"Primary alpha energy     = {args.primary_energy_mev:.6f} MeV")
    print(f"Primary sigma            = {args.primary_sigma_mev:.6f} MeV")
    print(f"Secondary E_star         = {e_star:.6f} MeV")
    print(f"Secondary E_boost        = {e_boost:.6f} MeV")
    print(f"Secondary energy range   = {sec_emin:.6f} -- {sec_emax:.6f} MeV")
    print(f"Secondary mean energy    = {e_star + e_boost:.6f} MeV")

    energies, labels = sample_p11b_alpha_spectrum(
        n_alpha=args.n_alpha,
        primary_energy_mev=args.primary_energy_mev,
        primary_fraction=args.primary_fraction,
        secondary_estar_mev=args.secondary_estar_mev,
        secondary_eboost_mev=args.secondary_eboost_mev,
        primary_sigma_mev=args.primary_sigma_mev,
        seed=args.seed,
    )

    print()
    print("=== sampled statistics ===")
    print(f"Mean energy              = {np.mean(energies):.6f} MeV")
    print(f"Min energy               = {np.min(energies):.6f} MeV")
    print(f"Max energy               = {np.max(energies):.6f} MeV")
    print(f"N primary                = {np.count_nonzero(labels == 'primary')}")
    print(f"N secondary              = {np.count_nonzero(labels == 'secondary')}")

    df = build_histogram(
        energies=energies,
        labels=labels,
        bins=args.bins,
        hist_min=args.emin_mev,
        hist_max=args.emax_mev,
    )

    hist_file = f"{args.out_prefix}_histogram.csv"
    df.to_csv(hist_file, index=False)
    plot_spectrum(df, args.out_prefix)

    print()
    print("=== output files ===")
    print(hist_file)
    print(f"{args.out_prefix}_spectrum_linear.png")
    print(f"{args.out_prefix}_spectrum_logy.png")

    if args.save_samples:
        sample_file = f"{args.out_prefix}_samples.csv"
        pd.DataFrame({
            "alpha_energy_MeV": energies,
            "component": labels,
        }).to_csv(sample_file, index=False)
        print(sample_file)


if __name__ == "__main__":
    main()
