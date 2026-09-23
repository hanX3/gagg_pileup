#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_brems_spectrum.py

Generate parameterized thermal bremsstrahlung photon spectra for different kTe.

The source photon number spectrum is approximated as

    dN/dE ∝ exp(-E/kTe) / E

where E and kTe are in the same energy unit. This is the photon-number
spectrum, suitable for drawing primary gamma/photon energies in Geant4.
It is not the detector-deposited-energy spectrum.

Example usage:

    python make_brems_spectrum.py \
        --kte-kev 40 100 160 300 \
        --emin-kev 10 \
        --emax-kev 5000 \
        --bins 1000 \
        --n-photons 1000000 \
        --out-prefix brems

Outputs:

    brems_spectra.csv
        Table of bin-center energy, normalized PDF, expected counts per bin,
        and CDF for each kTe.

    brems_spectra_linear.png
        Expected counts per bin versus photon energy, linear y-axis.

    brems_spectra_logy.png
        Expected counts per bin versus photon energy, log y-axis.

    brems_sampled_energies_kTeXXXkeV.csv  [optional]
        Random sampled photon energies, if --sample-size > 0 is set.

Notes:

1. The low-energy cutoff Emin is mandatory in practice because the photon
   number spectrum contains a 1/E factor and diverges as E -> 0.
2. For a 5 um Mylar front window, photons above about 10 keV are nearly
   transparent to the window, so Emin = 10 keV is a reasonable first setting.
3. The expected counts here are source photon counts. The GAGG deposited-energy
   spectrum will be different because of photoelectric absorption, Compton
   scattering, pair production, and escape.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def brems_pdf_unnormalized(energy_kev: np.ndarray, kte_kev: float) -> np.ndarray:
    """Unnormalized photon-number PDF: f(E) ∝ exp(-E/kTe) / E."""
    energy_kev = np.asarray(energy_kev, dtype=float)
    if kte_kev <= 0:
        raise ValueError("kte_kev must be positive.")
    if np.any(energy_kev <= 0):
        raise ValueError("All energies must be positive because the spectrum contains 1/E.")
    return np.exp(-energy_kev / kte_kev) / energy_kev


def make_energy_edges(emin_kev: float, emax_kev: float, bins: int, log_bins: bool) -> np.ndarray:
    if emin_kev <= 0:
        raise ValueError("emin_kev must be positive.")
    if emax_kev <= emin_kev:
        raise ValueError("emax_kev must be larger than emin_kev.")
    if bins <= 0:
        raise ValueError("bins must be positive.")

    if log_bins:
        return np.geomspace(emin_kev, emax_kev, bins + 1)
    return np.linspace(emin_kev, emax_kev, bins + 1)


def spectrum_for_kte(
    kte_kev: float,
    edges_kev: np.ndarray,
    n_photons: float,
) -> pd.DataFrame:
    """Return bin-center PDF, CDF, and expected counts for one kTe.

    The PDF is normalized over [Emin, Emax]. Counts are computed as

        counts_i = n_photons * pdf(E_i) * bin_width_i

    which is accurate enough for plotting and source-spectrum checks.
    For Geant4 sampling, use the CDF column or the sampling function below.
    """
    centers = 0.5 * (edges_kev[:-1] + edges_kev[1:])
    widths = edges_kev[1:] - edges_kev[:-1]

    f = brems_pdf_unnormalized(centers, kte_kev)
    norm = np.sum(f * widths)
    pdf = f / norm
    prob_per_bin = pdf * widths
    cdf = np.cumsum(prob_per_bin)
    cdf[-1] = 1.0
    counts = n_photons * prob_per_bin

    return pd.DataFrame(
        {
            "E_low_keV": edges_kev[:-1],
            "E_high_keV": edges_kev[1:],
            "E_center_keV": centers,
            f"pdf_kTe_{kte_kev:g}_keV": pdf,
            f"prob_per_bin_kTe_{kte_kev:g}_keV": prob_per_bin,
            f"expected_counts_kTe_{kte_kev:g}_keV": counts,
            f"cdf_kTe_{kte_kev:g}_keV": cdf,
        }
    )


def sample_energies_from_cdf(
    kte_kev: float,
    emin_kev: float,
    emax_kev: float,
    sample_size: int,
    cdf_grid_size: int = 200_000,
    seed: int | None = None,
) -> np.ndarray:
    """Draw random photon energies from f(E) ∝ exp(-E/kTe)/E.

    This uses numerical inverse-CDF sampling on a dense linear grid.
    """
    if sample_size <= 0:
        return np.array([], dtype=float)

    rng = np.random.default_rng(seed)
    grid = np.linspace(emin_kev, emax_kev, cdf_grid_size)
    f = brems_pdf_unnormalized(grid, kte_kev)

    # Build numerical CDF using trapezoidal cumulative integration.
    dE = np.diff(grid)
    area = 0.5 * (f[:-1] + f[1:]) * dE
    cdf = np.concatenate([[0.0], np.cumsum(area)])
    cdf /= cdf[-1]

    u = rng.random(sample_size)
    return np.interp(u, cdf, grid)


def merge_spectra(tables: list[pd.DataFrame]) -> pd.DataFrame:
    """Merge spectra tables that share the same energy binning."""
    if not tables:
        raise ValueError("No spectra tables to merge.")

    base_cols = ["E_low_keV", "E_high_keV", "E_center_keV"]
    merged = tables[0][base_cols].copy()

    for table in tables:
        for col in table.columns:
            if col not in base_cols:
                merged[col] = table[col].values
    return merged


def plot_spectra(
    spectra: pd.DataFrame,
    kte_values: list[float],
    out_prefix: str,
    log_y: bool,
) -> None:
    plt.figure(figsize=(8, 5))
    x = spectra["E_center_keV"].to_numpy()

    for kte in kte_values:
        y = spectra[f"expected_counts_kTe_{kte:g}_keV"].to_numpy()
        plt.plot(x, y, label=f"kTe = {kte:g} keV")

    plt.xlabel("Photon energy Eγ [keV]")
    plt.ylabel("Expected source counts / bin")
    plt.title("Thermal bremsstrahlung photon-number spectra")
    plt.legend()
    plt.grid(True, alpha=0.3)
    if log_y:
        plt.yscale("log")
        out_name = f"{out_prefix}_spectra_logy.png"
    else:
        out_name = f"{out_prefix}_spectra_linear.png"
    plt.tight_layout()
    plt.savefig(out_name, dpi=200)
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate Maxwellian thermal bremsstrahlung photon spectra for different kTe values."
    )
    parser.add_argument(
        "--kte-kev",
        nargs="+",
        type=float,
        default=[40.0, 100.0, 160.0, 300.0],
        help="kTe values in keV. Default: 40 100 160 300",
    )
    parser.add_argument("--emin-kev", type=float, default=10.0, help="Minimum photon energy in keV. Default: 10")
    parser.add_argument("--emax-kev", type=float, default=5000.0, help="Maximum photon energy in keV. Default: 5000")
    parser.add_argument("--bins", type=int, default=1000, help="Number of spectrum bins. Default: 1000")
    parser.add_argument(
        "--log-bins",
        action="store_true",
        help="Use logarithmic energy bins instead of linear bins.",
    )
    parser.add_argument(
        "--n-photons",
        type=float,
        default=1_000_000,
        help="Total source photon counts used to scale expected counts. Default: 1e6",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=0,
        help="If >0, also randomly sample this many photon energies for each kTe.",
    )
    parser.add_argument("--seed", type=int, default=12345, help="Random seed for sampled spectra. Default: 12345")
    parser.add_argument("--out-prefix", type=str, default="brems", help="Output filename prefix. Default: brems")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    for kte in args.kte_kev:
        if kte <= 0:
            raise ValueError("All kTe values must be positive.")

    edges = make_energy_edges(args.emin_kev, args.emax_kev, args.bins, args.log_bins)
    tables = [spectrum_for_kte(kte, edges, args.n_photons) for kte in args.kte_kev]
    spectra = merge_spectra(tables)

    csv_name = f"{args.out_prefix}_spectra.csv"
    spectra.to_csv(csv_name, index=False)

    plot_spectra(spectra, args.kte_kev, args.out_prefix, log_y=False)
    plot_spectra(spectra, args.kte_kev, args.out_prefix, log_y=True)

    if args.sample_size > 0:
        for idx, kte in enumerate(args.kte_kev):
            sampled = sample_energies_from_cdf(
                kte_kev=kte,
                emin_kev=args.emin_kev,
                emax_kev=args.emax_kev,
                sample_size=args.sample_size,
                seed=args.seed + idx,
            )
            sample_csv = f"{args.out_prefix}_sampled_energies_kTe{kte:g}keV.csv"
            pd.DataFrame({"primary_energy_keV": sampled, "primary_energy_MeV": sampled / 1000.0}).to_csv(
                sample_csv, index=False
            )

    print("Done.")
    print(f"Wrote: {csv_name}")
    print(f"Wrote: {args.out_prefix}_spectra_linear.png")
    print(f"Wrote: {args.out_prefix}_spectra_logy.png")
    if args.sample_size > 0:
        print("Wrote sampled energy CSV files for each kTe.")


if __name__ == "__main__":
    main()
