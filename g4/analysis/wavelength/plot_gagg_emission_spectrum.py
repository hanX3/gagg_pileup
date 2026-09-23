#!/usr/bin/env python3
"""
复现 Geant4 G4Scintillation 对 GAGG 发射谱的抽样, 并画出光子能量与波长分布.

抽样算法严格照抄 Geant4 的两步:
  1) BuildThePhysicsTable(): 梯形法累积积分 -> CDF 节点
  2) PostStepDoIt():         对逆 CDF 做线性插值 -> 反变换抽样
第 2 步意味着每个区间内能量是均匀分布的, 所以等效 pdf 是阶梯函数,
高度 = 0.5*(I_i + I_{i+1}), 而不是原始强度的折线.

终端:
    python3 plot_gagg_emission_spectrum.py
    python3 plot_gagg_emission_spectrum.py --src ../src/DetectorConstruction.cc -n 2000000
    python3 plot_gagg_emission_spectrum.py --component 2 --out slow.png

Jupyter:
    %run plot_gagg_emission_spectrum.py
  或
    import plot_gagg_emission_spectrum as g
    r = g.run(n_samples=500_000)
    r["wavelength_nm"]        # 抽样得到的波长数组
"""

from __future__ import annotations

import os
import re
import sys

import numpy as np

# h*c, 单位 eV*nm. E[eV] = HC / lambda[nm]
HC_EV_NM = 1239.841984

# DetectorConstruction.cc 里的默认值, 作为找不到源码时的回退
FALLBACK_ENERGY_EV = [2.00, 2.15, 2.25, 2.34, 2.43, 2.55, 2.75, 3.00]
FALLBACK_INTENSITY = [0.05, 0.20, 0.55, 1.00, 0.55, 0.20, 0.05, 0.01]

# Jupyter 里直接改这里, 或调 run(**kwargs) 传参
CONFIG = dict(
    src=None,          # DetectorConstruction.cc 路径; None 表示自动搜索
    component=1,       # 1 -> scint_fast, 2 -> scint_slow
    n_samples=1_000_000,
    bins=200,
    seed=20260624,
    out="gagg_emission_spectrum.png",
    show=True,
)


# ---------------------------------------------------------------- 源码解析

def _find_source() -> str | None:
    """从常见相对位置找 DetectorConstruction.cc"""
    candidates = [
        "DetectorConstruction.cc",
        "src/DetectorConstruction.cc",
        "../src/DetectorConstruction.cc",
        "../../src/DetectorConstruction.cc",
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return None


def parse_spectrum(path: str | None, component: int = 1):
    """
    从 DetectorConstruction.cc 解析 photon_energy[] 和 scint_fast/slow[].

    返回 (energy_eV, intensity, source_label).
    解析失败时回退到硬编码值.
    """
    array_name = "scint_fast" if component == 1 else "scint_slow"

    if path is None:
        path = _find_source()

    if path and os.path.isfile(path):
        try:
            text = open(path, "r", errors="replace").read()

            def grab(name, unit=None):
                # 匹配  G4double name[...] = { ... };
                m = re.search(
                    r"G4double\s+" + re.escape(name) + r"\s*\[[^\]]*\]\s*=\s*\{(.*?)\}\s*;",
                    text, re.S)
                if not m:
                    return None
                body = m.group(1)
                if unit:
                    vals = re.findall(r"([-+0-9.eE]+)\s*\*\s*" + unit, body)
                else:
                    vals = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", body)
                return [float(v) for v in vals] or None

            E = grab("photon_energy", unit="eV")
            I = grab(array_name)

            if E and I and len(E) == len(I):
                return (np.asarray(E, float), np.asarray(I, float),
                        f"{os.path.basename(path)} :: {array_name}[{len(E)}]")
            print(f"[warn] {path} 里没解析出成对的 photon_energy/{array_name}, 用回退值",
                  file=sys.stderr)
        except OSError as exc:
            print(f"[warn] 读取 {path} 失败: {exc}, 用回退值", file=sys.stderr)
    else:
        print("[info] 未找到 DetectorConstruction.cc, 用脚本内硬编码值 "
              "(可用 --src 指定路径)", file=sys.stderr)

    return (np.asarray(FALLBACK_ENERGY_EV, float),
            np.asarray(FALLBACK_INTENSITY, float),
            f"内置回退值 :: {array_name}[{len(FALLBACK_ENERGY_EV)}]")


# ---------------------------------------------------------- Geant4 抽样复现

def build_cdf(E: np.ndarray, I: np.ndarray) -> np.ndarray:
    """梯形法累积积分, 对应 G4Scintillation::BuildThePhysicsTable()"""
    trap = 0.5 * np.diff(E) * (I[:-1] + I[1:])
    return np.concatenate([[0.0], np.cumsum(trap)])


def sample_energy(E: np.ndarray, cdf: np.ndarray, n: int, rng) -> np.ndarray:
    """
    反变换抽样. 对逆 CDF 线性插值, 对应
    G4PhysicsOrderedFreeVector::GetEnergy() 的行为.
    """
    u = rng.random(n) * cdf[-1]
    return np.interp(u, cdf, E)


def effective_pdf_E(E: np.ndarray, I: np.ndarray):
    """
    Geant4 实际抽样对应的等效 pdf: 阶梯函数.
    返回 (bin_edges, heights), heights 已按能量归一化.
    """
    h = 0.5 * (I[:-1] + I[1:])          # 每个区间的常数高度
    norm = np.sum(h * np.diff(E))
    return E, h / norm


def naive_pdf_E(E: np.ndarray, I: np.ndarray, n_grid: int = 2000):
    """把原始强度当折线 pdf (Geant4 并非如此抽样), 用于对比"""
    grid = np.linspace(E[0], E[-1], n_grid)
    pdf = np.interp(grid, E, I)
    pdf /= np.trapezoid(pdf, grid) if hasattr(np, "trapezoid") else np.trapz(pdf, grid)
    return grid, pdf


def step_pdf_at(E: np.ndarray, heights: np.ndarray, x: np.ndarray) -> np.ndarray:
    """在任意能量点上取阶梯 pdf 的值"""
    idx = np.clip(np.searchsorted(E, x, side="right") - 1, 0, len(heights) - 1)
    val = heights[idx]
    return np.where((x >= E[0]) & (x <= E[-1]), val, 0.0)


# -------------------------------------------------------------- 环境判定

def _in_notebook() -> bool:
    try:
        from IPython import get_ipython
        ip = get_ipython()
        return ip is not None and ip.__class__.__name__ == "ZMQInteractiveShell"
    except Exception:
        return False


def _setup_matplotlib(want_show: bool):
    """终端无显示时切 Agg, 避免报错"""
    import matplotlib
    if not _in_notebook():
        headless = (sys.platform.startswith("linux")
                    and not os.environ.get("DISPLAY")
                    and not os.environ.get("WAYLAND_DISPLAY"))
        if headless or not want_show:
            matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # 有中文字体就用上, 没有就保持默认 (图内文字本身是英文, 不影响)
    try:
        from matplotlib import font_manager
        have = {f.name for f in font_manager.fontManager.ttflist}
        for cand in ("Noto Sans CJK SC", "Source Han Sans SC", "WenQuanYi Zen Hei",
                     "Microsoft YaHei", "PingFang SC", "SimHei"):
            if cand in have:
                plt.rcParams["font.sans-serif"] = [cand] + plt.rcParams["font.sans-serif"]
                plt.rcParams["axes.unicode_minus"] = False
                break
    except Exception:
        pass

    return plt


# ------------------------------------------------------------------ 主流程

def run(**kwargs):
    cfg = {**CONFIG, **kwargs}

    E, I, label = parse_spectrum(cfg["src"], cfg["component"])
    if np.any(np.diff(E) <= 0):
        raise ValueError("photon_energy 必须严格单调递增")

    cdf = build_cdf(E, I)
    rng = np.random.default_rng(cfg["seed"])
    e_samp = sample_energy(E, cdf, int(cfg["n_samples"]), rng)
    lam_samp = HC_EV_NM / e_samp

    edges, heights = effective_pdf_E(E, I)

    # ---- 文字输出 ----
    lam_nodes = HC_EV_NM / E
    print(f"谱形来源      : {label}")
    print(f"能量范围      : {E[0]:.3f} - {E[-1]:.3f} eV")
    print(f"波长范围      : {lam_nodes[-1]:.1f} - {lam_nodes[0]:.1f} nm")
    print(f"抽样数        : {len(e_samp):,}")
    print()
    print(f"平均光子能量  : {e_samp.mean():.4f} eV  (-> {HC_EV_NM/e_samp.mean():.2f} nm)")
    print(f"平均波长      : {lam_samp.mean():.2f} nm   "
          f"(注意 != hc/<E>, 因为 1/E 非线性)")
    print(f"波长中位数    : {np.median(lam_samp):.2f} nm")
    print(f"波长标准差    : {lam_samp.std():.2f} nm")
    ipk = int(np.argmax(heights))
    print(f"最高阶梯区间  : {E[ipk]:.3f} - {E[ipk+1]:.3f} eV  "
          f"= {HC_EV_NM/E[ipk+1]:.1f} - {HC_EV_NM/E[ipk]:.1f} nm")
    print()

    # ---- 绘图 ----
    plt = _setup_matplotlib(cfg["show"])
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))

    # 左: 光子能量
    ax[0].hist(e_samp, bins=cfg["bins"], density=True, histtype="stepfilled",
               alpha=0.35, color="#377eb8", label=f"sampled (N={len(e_samp):,})")
    xs = np.linspace(E[0], E[-1], 4000)
    ax[0].plot(xs, step_pdf_at(edges, heights, xs), color="#e41a1c", lw=1.6,
               label="Geant4 effective pdf (staircase)")
    g, p = naive_pdf_E(E, I)
    ax[0].plot(g, p, color="#4daf4a", lw=1.2, ls="--",
               label="naive piecewise-linear (NOT Geant4)")
    ax[0].plot(E, I / np.sum(0.5*(I[:-1]+I[1:]) * np.diff(E)), "ko", ms=4,
               label=f"table nodes (n={len(E)})")
    ax[0].set_xlabel("photon energy [eV]")
    ax[0].set_ylabel("pdf [1/eV]")
    ax[0].set_title("emission spectrum vs photon energy")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.25)

    # 右: 波长, 需要雅可比 p(lam) = p(E) * hc / lam^2
    ax[1].hist(lam_samp, bins=cfg["bins"], density=True, histtype="stepfilled",
               alpha=0.35, color="#377eb8", label="sampled")
    lam_grid = np.linspace(lam_nodes[-1], lam_nodes[0], 4000)
    pdf_lam = step_pdf_at(edges, heights, HC_EV_NM / lam_grid) * HC_EV_NM / lam_grid**2
    ax[1].plot(lam_grid, pdf_lam, color="#e41a1c", lw=1.6,
               label="analytic pdf (with Jacobian)")
    ax[1].axvline(530.0, color="k", ls=":", lw=1,
                  label="GAGG(Ce) measured peak ~530 nm")
    ax[1].set_xlabel("wavelength [nm]")
    ax[1].set_ylabel("pdf [1/nm]")
    ax[1].set_title("emission spectrum vs wavelength")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.25)

    fig.suptitle(f"GAGG scintillation emission, component {cfg['component']}  |  {label}",
                 fontsize=10)
    fig.tight_layout()

    if cfg["out"]:
        fig.savefig(cfg["out"], dpi=150, bbox_inches="tight")
        print(f"图已保存: {cfg['out']}")

    if cfg["show"]:
        plt.show()

    return dict(energy_eV=E, intensity=I, cdf=cdf,
                sampled_energy_eV=e_samp, wavelength_nm=lam_samp,
                step_edges=edges, step_heights=heights, fig=fig, ax=ax)


def _cli():
    import argparse
    p = argparse.ArgumentParser(description="画 GAGG 闪烁发射谱 (复现 Geant4 抽样)")
    p.add_argument("--src", default=CONFIG["src"],
                   help="DetectorConstruction.cc 路径 (默认自动搜索)")
    p.add_argument("--component", type=int, choices=(1, 2), default=CONFIG["component"],
                   help="1 = scint_fast, 2 = scint_slow")
    p.add_argument("-n", "--n-samples", type=int, default=CONFIG["n_samples"])
    p.add_argument("--bins", type=int, default=CONFIG["bins"])
    p.add_argument("--seed", type=int, default=CONFIG["seed"])
    p.add_argument("--out", default=CONFIG["out"], help="输出 PNG; 空字符串表示不存")
    p.add_argument("--no-show", action="store_true", help="不弹窗, 只存文件")
    a = p.parse_args()
    return run(src=a.src, component=a.component, n_samples=a.n_samples,
               bins=a.bins, seed=a.seed, out=a.out or None, show=not a.no_show)


if __name__ == "__main__":
    if _in_notebook():
        _result = run()
    else:
        _result = _cli()
