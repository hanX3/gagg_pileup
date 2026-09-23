#!/usr/bin/env python3
# =============================================================================
# make_birks_yield_curves.py
#
# Purpose:
#   Generate Geant4 particle-dependent total scintillation light-yield curves
#   for GAGG(Ce), using electron, alpha, and 16O mass stopping-power tables.
#
# Put this script in:
#   ./
#
# Required input files in the same directory:
#   electron_gagg_stopping.csv
#   alpha_gagg_stopping.csv
#   o16_gagg_stopping.csv
#
# Required input CSV format:
#   E_MeV,S_mass_MeV_cm2_g
#   0.01,13.09
#   ...
#
# Recommended command:
#   python3 make_birks_yield_curves.py \
#     --electron electron_gagg_stopping.csv \
#     --alpha alpha_gagg_stopping.csv \
#     --o16 o16_gagg_stopping.csv \
#     --out gagg_birks_yield_curves.inc \
#     --emax 10.0 \
#     --npoints 501 \
#     --plot
#
# Output files:
#   gagg_birks_yield_curves.inc
#       C++ include file containing total Geant4 yield vectors.
#
#   gagg_birks_yield_curves.csv
#       Numerical table for checking dL/dE and cumulative light yield.
#
#   gagg_birks_yield_curves.png
#       Diagnostic plot, generated when --plot is used.
#
# Physics model:
#   dL/dE = LY / (1 + a1 * S_mass)
#
# Default parameters:
#   LY = 46000 photons/MeV
#   a1 = 6.5e-3 g cm^-2 MeV^-1
#
# Important:
#   If you use the generated particle-dependent yield curves in Geant4,
#   do NOT also call:
#       G4Material::GetIonisation()->SetBirksConstant(...)
#   Otherwise quenching will be applied twice.
# =============================================================================

"""
Generate Geant4 particle-dependent total scintillation light-yield curves for
GAGG(Ce) from electron, alpha, and 16O mass stopping-power tables.

Typical usage
-------------
Put this script in the directory containing:

    electron_gagg_stopping.csv
    alpha_gagg_stopping.csv
    o16_gagg_stopping.csv

For example:

    cd ~/GAGG_pileup/stopping_power

    python3 make_birks_yield_curves.py \
      --electron electron_gagg_stopping.csv \
      --alpha alpha_gagg_stopping.csv \
      --o16 o16_gagg_stopping.csv \
      --out gagg_birks_yield_curves.inc \
      --emax 10.0 \
      --npoints 501 \
      --plot

Expected input CSV format
-------------------------
The electron, alpha, and 16O stopping-power files must have this header:

    E_MeV,S_mass_MeV_cm2_g
    0.01,13.09
    ...

Here:
    E_MeV              : kinetic energy [MeV]
                         For 16O, this is the total ion kinetic energy, not MeV/u.
    S_mass_MeV_cm2_g   : mass stopping power [MeV cm^2/g]

Birks model
-----------
    dL/dE = LY / (1 + a1 * S_mass)

where
    LY      : electron-equivalent light yield [photons/MeV]
    a1      : Birks parameter [g cm^-2 MeV^-1]
    S_mass  : mass stopping power [MeV cm^2/g]

Default parameters
------------------
    LY = 46000 photons/MeV
    a1 = 6.5e-3 g cm^-2 MeV^-1

Outputs
-------
    1) gagg_birks_yield_curves.inc
       C++ include file with total std::vector<G4double> arrays for Geant4.

    2) gagg_birks_yield_curves.csv
       Numerical table for checking the generated curves.

    3) gagg_birks_yield_curves.png
       Plot generated when --plot is used.

Important note
--------------
If these particle-dependent yield curves are used in Geant4, do not also apply
G4Material::GetIonisation()->SetBirksConstant(...), otherwise quenching will be
applied twice.

This script only generates total cumulative light-yield curves. Pulse-shape
component fractions should be configured separately in the Geant4 material
properties, for example through constants in Constants.hh.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Tuple

import numpy as np


def read_stopping_csv(path: str | Path) -> Tuple[np.ndarray, np.ndarray]:
    """Read stopping-power CSV with columns E_MeV and S_mass_MeV_cm2_g."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Cannot find file: {path}")

    data = np.genfromtxt(path, delimiter=",", names=True)

    if data.size == 0:
        raise ValueError(f"No data read from {path}")

    names = data.dtype.names
    if names is None:
        raise ValueError(f"{path} does not appear to have a CSV header")

    if "E_MeV" not in names or "S_mass_MeV_cm2_g" not in names:
        raise ValueError(
            f"{path} must contain columns: E_MeV,S_mass_MeV_cm2_g. "
            f"Found: {names}"
        )

    E = np.asarray(data["E_MeV"], dtype=float)
    S = np.asarray(data["S_mass_MeV_cm2_g"], dtype=float)

    mask = np.isfinite(E) & np.isfinite(S) & (E > 0.0) & (S > 0.0)
    E = E[mask]
    S = S[mask]

    order = np.argsort(E)
    E = E[order]
    S = S[order]

    # Remove duplicate energies if any.
    unique_E, unique_idx = np.unique(E, return_index=True)
    E = unique_E
    S = S[unique_idx]

    if len(E) < 2:
        raise ValueError(f"Not enough valid positive data points in {path}")

    return E, S


def loglog_interp(x: np.ndarray, xp: np.ndarray, fp: np.ndarray) -> np.ndarray:
    """Log-log interpolation with endpoint clamping."""
    x = np.asarray(x, dtype=float)
    xp = np.asarray(xp, dtype=float)
    fp = np.asarray(fp, dtype=float)

    # Clamp positive x into input range. This avoids extrapolation blow-ups.
    x_clamped = np.clip(x, xp[0], xp[-1])

    return np.exp(
        np.interp(
            np.log(x_clamped),
            np.log(xp),
            np.log(fp),
        )
    )


def make_energy_grid(emin: float, emax: float, npoints: int) -> np.ndarray:
    """Return energy grid in MeV. Includes 0 as first point."""
    if emin <= 0.0:
        raise ValueError("emin must be > 0 for logarithmic grid")
    if emax <= emin:
        raise ValueError("emax must be larger than emin")
    if npoints < 10:
        raise ValueError("npoints should be at least 10")

    positive_grid = np.logspace(math.log10(emin), math.log10(emax), npoints - 1)
    return np.concatenate(([0.0], positive_grid))


def birks_dlde(light_yield_per_mev: float,
               birks_a1_g_cm2_per_mev: float,
               stopping_mass: np.ndarray) -> np.ndarray:
    """Birks differential light output dL/dE [photons/MeV]."""
    return light_yield_per_mev / (1.0 + birks_a1_g_cm2_per_mev * stopping_mass)


def cumulative_trapezoid(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Cumulative trapezoidal integral with L(0)=0."""
    out = np.zeros_like(x, dtype=float)
    dx = np.diff(x)
    area = 0.5 * (y[1:] + y[:-1]) * dx
    out[1:] = np.cumsum(area)
    return out


def make_yield_curve(E_grid: np.ndarray,
                     E_stop: np.ndarray,
                     S_stop: np.ndarray,
                     light_yield_per_mev: float,
                     birks_a1_g_cm2_per_mev: float) -> Tuple[np.ndarray, np.ndarray]:
    """Return dL/dE on the grid and cumulative light-yield curve."""
    S_grid = np.empty_like(E_grid)

    # At exactly E=0, use the first available stopping point only to define
    # the first trapezoid. The cumulative yield itself remains L(0)=0.
    positive = E_grid > 0.0
    S_grid[0] = S_stop[0]
    S_grid[positive] = loglog_interp(E_grid[positive], E_stop, S_stop)

    dlde = birks_dlde(light_yield_per_mev, birks_a1_g_cm2_per_mev, S_grid)
    L = cumulative_trapezoid(E_grid, dlde)
    L[0] = 0.0
    return dlde, L


def format_cpp_vector(name: str, values: np.ndarray, unit: str | None = None) -> str:
    """Format a std::vector<G4double>. If unit is given, multiply each value by unit."""
    lines = [f"static const std::vector<G4double> {name} = {{"]
    chunks = []
    for v in values:
        if unit:
            chunks.append(f"{v:.12g}*{unit}")
        else:
            chunks.append(f"{v:.12g}")

    for i in range(0, len(chunks), 4):
        lines.append("  " + ", ".join(chunks[i:i + 4]) + ("," if i + 4 < len(chunks) else ""))
    lines.append("};")
    return "\n".join(lines)


def write_cpp_include(path: str | Path,
                      E_grid: np.ndarray,
                      L_e: np.ndarray,
                      L_a: np.ndarray,
                      L_o16: np.ndarray) -> None:
    """Write C++ include file containing total energy and yield vectors."""
    path = Path(path)

    text = []
    text.append("// Auto-generated by make_birks_yield_curves.py")
    text.append("// Energies are particle kinetic energies for Geant4 particle-dependent scintillation yield vectors.")
    text.append("// Yields are total cumulative photon yields L(E), not dL/dE.")
    text.append("// For 16O, the kinetic energy is the total ion kinetic energy, not MeV/u.")
    text.append("// Do not apply G4Material::GetIonisation()->SetBirksConstant() again if these curves already include Birks quenching.")
    text.append("")
    text.append("#ifndef GAGG_BIRKS_YIELD_CURVES_INC")
    text.append("#define GAGG_BIRKS_YIELD_CURVES_INC")
    text.append("")
    text.append("#include <vector>")
    text.append("")
    text.append(format_cpp_vector("gaggScintillationEnergy", E_grid, "MeV"))
    text.append("")
    text.append(format_cpp_vector("gaggGammaElectronScintillationYield", L_e))
    text.append("")
    text.append(format_cpp_vector("gaggAlphaScintillationYield", L_a))
    text.append("")
    text.append("// 16O total cumulative light-yield curve.")
    text.append("// Intended for O16 recoil/ion fallback, for example as IONSCINTILLATIONYIELD.")
    text.append(format_cpp_vector("gaggO16ScintillationYield", L_o16))
    text.append("")
    text.append("#endif  // GAGG_BIRKS_YIELD_CURVES_INC")
    text.append("")

    path.write_text("\n".join(text), encoding="utf-8")


def write_check_csv(path: str | Path,
                    E_grid: np.ndarray,
                    dlde_e: np.ndarray,
                    L_e: np.ndarray,
                    dlde_a: np.ndarray,
                    L_a: np.ndarray,
                    dlde_o16: np.ndarray,
                    L_o16: np.ndarray) -> None:
    """Write inspection CSV."""
    path = Path(path)
    arr = np.column_stack([
        E_grid,
        dlde_e,
        L_e,
        dlde_a,
        L_a,
        dlde_o16,
        L_o16,
        L_a / np.maximum(L_e, 1e-300),
        L_o16 / np.maximum(L_e, 1e-300),
    ])
    header = (
        "E_MeV,"
        "electron_dLdE_ph_per_MeV,electron_L_ph,"
        "alpha_dLdE_ph_per_MeV,alpha_L_ph,"
        "o16_dLdE_ph_per_MeV,o16_L_ph,"
        "alpha_over_electron_L,o16_over_electron_L"
    )
    np.savetxt(path, arr, delimiter=",", header=header, comments="", fmt="%.12g")


def make_plot(path: str | Path,
              E_grid: np.ndarray,
              dlde_e: np.ndarray,
              L_e: np.ndarray,
              dlde_a: np.ndarray,
              L_a: np.ndarray,
              dlde_o16: np.ndarray,
              L_o16: np.ndarray) -> None:
    """Make a quick diagnostic plot."""
    import matplotlib.pyplot as plt

    positive = E_grid > 0.0

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    axes[0].semilogx(E_grid[positive], dlde_e[positive], label="gamma/electron")
    axes[0].semilogx(E_grid[positive], dlde_a[positive], label="alpha")
    axes[0].semilogx(E_grid[positive], dlde_o16[positive], label="O16")
    axes[0].set_xlabel("Particle kinetic energy [MeV]")
    axes[0].set_ylabel("dL/dE [photons/MeV]")
    axes[0].set_title("Birks differential light yield")
    axes[0].legend()

    axes[1].plot(E_grid, L_e, label="gamma/electron")
    axes[1].plot(E_grid, L_a, label="alpha")
    axes[1].plot(E_grid, L_o16, label="O16")
    axes[1].set_xlabel("Particle kinetic energy [MeV]")
    axes[1].set_ylabel("Cumulative light yield L(E) [photons]")
    axes[1].set_title("Cumulative scintillation yield")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate Birks-quenched total GAGG scintillation yield curves for Geant4."
    )
    parser.add_argument("--electron", required=True, help="electron stopping CSV: E_MeV,S_mass_MeV_cm2_g")
    parser.add_argument("--alpha", required=True, help="alpha stopping CSV: E_MeV,S_mass_MeV_cm2_g")
    parser.add_argument("--o16", required=True, help="16O stopping CSV: E_MeV,S_mass_MeV_cm2_g. E_MeV is total ion kinetic energy, not MeV/u.")
    parser.add_argument("--out", default="gagg_birks_yield_curves.inc", help="output C++ include file")
    parser.add_argument("--csv", default="gagg_birks_yield_curves.csv", help="output inspection CSV")
    parser.add_argument("--plot-file", default="gagg_birks_yield_curves.png", help="output plot filename")
    parser.add_argument("--plot", action="store_true", help="also write diagnostic PNG plot")
    parser.add_argument("--emin", type=float, default=0.01, help="minimum positive energy in MeV")
    parser.add_argument("--emax", type=float, default=10.0, help="maximum energy in MeV")
    parser.add_argument("--npoints", type=int, default=501, help="number of output energy points including E=0")
    parser.add_argument("--light-yield", type=float, default=46000.0, help="electron-equivalent light yield [photons/MeV]")
    parser.add_argument("--birks-a1", type=float, default=6.5e-3, help="Birks a1 [g cm^-2 MeV^-1]")

    args = parser.parse_args()

    E_e, S_e = read_stopping_csv(args.electron)
    E_a, S_a = read_stopping_csv(args.alpha)
    E_o16, S_o16 = read_stopping_csv(args.o16)

    E_grid = make_energy_grid(args.emin, args.emax, args.npoints)

    dlde_e, L_e = make_yield_curve(E_grid, E_e, S_e, args.light_yield, args.birks_a1)
    dlde_a, L_a = make_yield_curve(E_grid, E_a, S_a, args.light_yield, args.birks_a1)
    dlde_o16, L_o16 = make_yield_curve(E_grid, E_o16, S_o16, args.light_yield, args.birks_a1)

    write_cpp_include(
        args.out,
        E_grid,
        L_e,
        L_a,
        L_o16,
    )

    write_check_csv(args.csv, E_grid, dlde_e, L_e, dlde_a, L_a, dlde_o16, L_o16)

    if args.plot:
        make_plot(args.plot_file, E_grid, dlde_e, L_e, dlde_a, L_a, dlde_o16, L_o16)

    print("Generated:")
    print(f"  {args.out}")
    print(f"  {args.csv}")
    if args.plot:
        print(f"  {args.plot_file}")

    # Print useful reference values.
    for x in [0.162, 0.662, 5.486, 10.0]:
        if 0.0 <= x <= args.emax:
            Le_x = np.interp(x, E_grid, L_e)
            La_x = np.interp(x, E_grid, L_a)
            Lo16_x = np.interp(x, E_grid, L_o16)
            print(
                f"E={x:.3f} MeV: "
                f"L_e={Le_x:.6g} ph, "
                f"L_alpha={La_x:.6g} ph, "
                f"L_O16={Lo16_x:.6g} ph, "
                f"alpha/e={La_x / Le_x if Le_x > 0 else float('nan'):.6g}, "
                f"O16/e={Lo16_x / Le_x if Le_x > 0 else float('nan'):.6g}"
            )


if __name__ == "__main__":
    main()
