#!/usr/bin/env python3
"""
画 GAGG:Ce 发射谱。

读入的是 C++ 数组文本文件（默认 gagg_emission_arrays.txt），从中解析
  G4double photon_energy[...] = { ... *eV, ... };   // 光子能量，eV
  G4double scint_fast[...]    = { ... };            // 相对强度
两个数组，然后画出：
  左图 - 能量域：解析出的谱形节点
  右图 - 波长域：按 Geant4 抽样逻辑（阶梯 pdf + 反变换抽样）得到的波长分布

这些数组就是将来喂进 DetectorConstruction.cc 的那份，所以画的就是模拟实际用的谱。

终端：
    python3 plot_wavelength_spectrum.py
    python3 plot_wavelength_spectrum.py --txt gagg_emission_arrays.txt
    python3 plot_wavelength_spectrum.py --intensity-array scint_slow --out s.png --no-show

Jupyter：
    %run plot_wavelength_spectrum.py
  或
    import plot_wavelength_spectrum as g
    r = g.run(txt="gagg_emission_arrays.txt")
    r["energy_eV"], r["intensity"], r["wavelength_nm"]   # 解析出的节点
    r["sampled_wavelength_nm"]                            # 抽样得到的波长数组
"""

from __future__ import annotations

import os
import re
import sys

import numpy as np

HC_EV_NM = 1239.841984  # h*c，E[eV] = HC / lambda[nm]

CONFIG = dict(
    txt=None,                       # C++ 数组 txt 路径；None 时自动搜索
    energy_array="photon_energy",   # 能量数组名
    intensity_array="scint_fast",   # 强度数组名（可换 scint_slow）
    n_samples=1_000_000,
    bins=160,
    seed=20260626,
    out="wavelength_spectrum.png",
    show=True,
)


# ------------------------------------------------------------ 解析 txt

def _find_txt():
    for c in ("gagg_emission_arrays.txt", "gagg_emission_arrays_64.txt",
              "gagg_emission_arrays_40.txt"):
        if os.path.isfile(c):
            return c
    return None


def parse_cpp_array(text, name):
    """
    从 C++ 源文本里提取名为 name 的数组，返回浮点列表。
    支持 G4double name[...] = { 1.55*eV, ... }; 或不带单位的 { 0.02, ... };
    单位记号（*eV, *nm 等）与注释都会被忽略，只取数值。
    """
    m = re.search(
        r"\b" + re.escape(name) + r"\s*\[[^\]]*\]\s*=\s*\{(.*?)\}\s*;",
        text, re.S)
    if not m:
        return None
    body = m.group(1)
    body = re.sub(r"//[^\n]*", "", body)
    body = re.sub(r"/\*.*?\*/", "", body, flags=re.S)
    nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", body)
    return [float(x) for x in nums] if nums else None


def load_arrays(path, energy_name, intensity_name):
    if path is None:
        path = _find_txt()
    if path is None or not os.path.isfile(path):
        raise FileNotFoundError(
            "找不到 txt。用 --txt 指定路径，或把文件放到当前目录。")
    text = open(path, "r", errors="replace").read()

    E = parse_cpp_array(text, energy_name)
    I = parse_cpp_array(text, intensity_name)
    if E is None:
        raise ValueError("txt 里没找到数组 '%s'" % energy_name)
    if I is None:
        raise ValueError("txt 里没找到数组 '%s'" % intensity_name)
    if len(E) != len(I):
        raise ValueError("两个数组长度不一致：%s=%d, %s=%d"
                         % (energy_name, len(E), intensity_name, len(I)))

    E = np.asarray(E, float)
    I = np.asarray(I, float)
    if np.any(np.diff(E) < 0):
        o = np.argsort(E)
        E, I = E[o], I[o]
    return E, I, path


# ------------------------------------------------ Geant4 抽样复现（阶梯）

def build_cdf(E, I):
    trap = 0.5 * np.diff(E) * (I[:-1] + I[1:])
    return np.concatenate([[0.0], np.cumsum(trap)])


def sample_energy(E, cdf, n, rng):
    return np.interp(rng.random(n) * cdf[-1], cdf, E)


def stair_pdf_lambda(E, I, lam_x):
    h = 0.5 * (I[:-1] + I[1:])
    norm = np.sum(h * np.diff(E))
    h = h / norm
    Ex = HC_EV_NM / lam_x
    idx = np.clip(np.searchsorted(E, Ex, side="right") - 1, 0, len(h) - 1)
    pdf_E = np.where((Ex >= E[0]) & (Ex <= E[-1]), h[idx], 0.0)
    return pdf_E * HC_EV_NM / lam_x**2


# -------------------------------------------------------------- 环境处理

def _in_notebook():
    try:
        from IPython import get_ipython
        ip = get_ipython()
        return ip is not None and ip.__class__.__name__ == "ZMQInteractiveShell"
    except Exception:
        return False


def _get_plt(want_show):
    import matplotlib
    if not _in_notebook():
        headless = (sys.platform.startswith("linux")
                    and not os.environ.get("DISPLAY")
                    and not os.environ.get("WAYLAND_DISPLAY"))
        if headless or not want_show:
            matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


# ------------------------------------------------------------------ 主流程

def run(**kwargs):
    cfg = dict(CONFIG)
    cfg.update(kwargs)

    E, I, path = load_arrays(cfg["txt"], cfg["energy_array"],
                             cfg["intensity_array"])
    I_norm = I / I.max()
    lam = HC_EV_NM / E

    cdf = build_cdf(E, I_norm)
    rng = np.random.default_rng(cfg["seed"])
    e_samp = sample_energy(E, cdf, int(cfg["n_samples"]), rng)
    lam_samp = HC_EV_NM / e_samp

    ipk = int(np.argmax(I_norm))
    print("txt 文件          : %s" % path)
    print("能量数组          : %s  (%d 点)" % (cfg["energy_array"], len(E)))
    print("强度数组          : %s" % cfg["intensity_array"])
    print("能量范围          : %.4f - %.4f eV" % (E[0], E[-1]))
    print("波长范围          : %.1f - %.1f nm" % (lam.min(), lam.max()))
    print("峰值波长          : %.1f nm" % lam[ipk])
    print("抽样平均波长      : %.2f nm  (红端拖尾使其 > 峰值)" % lam_samp.mean())
    print("抽样波长中位数    : %.2f nm" % np.median(lam_samp))

    plt = _get_plt(cfg["show"])
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))

    ax[0].plot(E, I_norm, "-o", ms=3.5, color="#1a7f37", lw=1.2,
               label="%s (%d pts)" % (cfg["intensity_array"], len(E)))
    ax[0].set_xlabel("photon energy [eV]")
    ax[0].set_ylabel("relative intensity")
    ax[0].set_title("emission spectrum (energy domain)")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.25)

    ax[1].hist(lam_samp, bins=cfg["bins"], density=True, histtype="stepfilled",
               alpha=0.35, color="#377eb8", label="Geant4-style sampled")
    lg = np.linspace(lam.min(), lam.max(), 4000)
    ax[1].plot(lg, stair_pdf_lambda(E, I_norm, lg), color="#e41a1c", lw=1.4,
               label="analytic staircase pdf")
    ax[1].axvline(540, ls=":", color="k", lw=1, label="Kobayashi RL peak 540 nm")
    ax[1].set_xlabel("wavelength [nm]")
    ax[1].set_ylabel("pdf [1/nm]")
    ax[1].set_title("resulting sampled spectrum (wavelength)")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.25)

    fig.suptitle("GAGG:Ce emission spectrum (from C++ array file)", fontsize=11)
    fig.tight_layout()

    if cfg["out"]:
        fig.savefig(cfg["out"], dpi=150, bbox_inches="tight")
        print("图已保存          : %s" % cfg["out"])
    if cfg["show"]:
        plt.show()

    return dict(energy_eV=E, intensity=I_norm, wavelength_nm=lam, cdf=cdf,
                sampled_energy_eV=e_samp, sampled_wavelength_nm=lam_samp,
                fig=fig, ax=ax)


def _cli():
    import argparse
    p = argparse.ArgumentParser(description="从 C++ 数组 txt 画 GAGG 发射谱")
    p.add_argument("--txt", default=CONFIG["txt"],
                   help="C++ 数组 txt 路径（默认自动搜索当前目录）")
    p.add_argument("--energy-array", default=CONFIG["energy_array"],
                   help="能量数组名（默认 photon_energy）")
    p.add_argument("--intensity-array", default=CONFIG["intensity_array"],
                   help="强度数组名（默认 scint_fast，可换 scint_slow）")
    p.add_argument("--n-samples", type=int, default=CONFIG["n_samples"])
    p.add_argument("--bins", type=int, default=CONFIG["bins"])
    p.add_argument("--seed", type=int, default=CONFIG["seed"])
    p.add_argument("--out", default=CONFIG["out"],
                   help="输出 PNG；传空串则不保存")
    p.add_argument("--no-show", action="store_true", help="不弹窗，只存文件")
    a = p.parse_args()
    return run(txt=a.txt, energy_array=a.energy_array,
               intensity_array=a.intensity_array, n_samples=a.n_samples,
               bins=a.bins, seed=a.seed, out=a.out or None, show=not a.no_show)


if __name__ == "__main__":
    if _in_notebook():
        _result = run()
    else:
        _result = _cli()
