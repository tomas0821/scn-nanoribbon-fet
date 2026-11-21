"""T1 figures and tables from fits/t1/ (rows/, sk/, dft.json, cg/, corr/, pitch/, cg_array/). Reads stored JSON only.

    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/t1_figs.py figs    # manuscript figures only; writes nothing in fits/
    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/t1_figs.py         # legacy: tables + Fig4, Fig5, Fig7 (+ parallel plate)
    ... t1_figs.py corr                                                    # legacy: fits/t1/corr/summary.md + Fig7b
    ... t1_figs.py pitch                                                   # legacy: Fig8
Figures (sized for the sn-jnl / Supplementary text widths, scripts/figstyle.py; re-rendered 2026-11-19):
  figures/Fig4_Eg_W.png (main, \\linewidth), Fig5_valley_mass_W.png (main, \\linewidth),
  Fig7b_Ion_W_corrected.png (main Fig. 7, \\linewidth), Fig7_Ion_W.png (SM, \\linewidth), Fig8_Ion_pitch.png (SM).
The per-width array curves (fits/t1/cg_array/array.json, periodic Laplace C_G at s = 2 and 7.7 nm) are drawn in
Fig7b (× B, 2D), Fig7 (ballistic) and Fig8 (× B, 2D), next to the parallel-plate and isolated-gate curves.
Tables (legacy modes only): fits/t1/t1_table.md, fits/t1/cg/summary.md, fits/t1/corr/summary.md.
Labels: "ballistic, ideal-electrostatics upper bound" (PREREG_t1).
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as FS  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
T1 = ROOT / "fits" / "t1"
FIG = ROOT / "figures"
IRDS = 1.979
# Okabe–Ito (CVD-safe), fixed per family
COL = {"100": "#0072B2", "110": "#D55E00", "110h": "#009E73", "100h": "#CC79A7"}
LAB = {"100": "[100] bare", "110": "[110] bare", "110h": "[110] H (in-sample only)", "100h": "[100] H"}   # tables
FLAB = {"100": "[100] bare", "110": "[110] bare", "110h": "[110] H (in-sample)", "100h": "[100] H"}     # figures
SKC = "#56B4E9"
INK = "#333333"
IRDS_KW = dict(color="#8c8c8c", lw=0.8, ls=(0, (4, 2.5)), zorder=1)
ARR = {"s2": dict(ls="-", lw=2.6, alpha=0.4, solid_capstyle="round"),            # array under one gate, s = 2 nm
       "sGAA": dict(ls=(0, (0.05, 1.6)), lw=2.6, alpha=0.5, dash_capstyle="round")}   # s = 2 t_HfO2 = 7.7 nm


def rows(fam):
    out = [json.loads(p.read_text()) for p in sorted((T1 / "rows").glob(f"{fam}_n*.json"))]
    return sorted(out, key=lambda r: r["n"])


def sk_rows():
    return sorted((json.loads(p.read_text()) for p in (T1 / "sk").glob("100_n*.json")), key=lambda r: r["n"])


def series(rs, f):
    return np.array([r["W_A"] / 10 for r in rs]), np.array([f(r) for r in rs])


def cg_rows(geom, fam):
    return sorted((json.loads(p.read_text()) for p in (T1 / "cg" / geom).glob(f"{fam}_n*.json")), key=lambda r: r["n"])


def corr_rows(fam):
    return sorted((json.loads(p.read_text()) for p in (T1 / "corr").glob(f"{fam}_n*.json")), key=lambda r: r["n"])


def array_data():
    f = T1 / "cg_array" / "array.json"
    return json.loads(f.read_text()) if f.exists() else None


def style(ax):
    ax.grid(color="#e5e5e5", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def irds_handle():
    return Line2D([], [], **{k: v for k, v in IRDS_KW.items() if k != "zorder"})


def arr_handle(s, color=INK):
    return Line2D([], [], color=color, **ARR[s])


def square(ax, x, y, c, mk="s", ms=3.8, clip=True):
    ax.plot(x, y, mk, mfc="white", mec=c, mew=0.9, ms=ms, zorder=6, clip_on=clip)


# ==================================================================================================== Fig. 4
def fig4(fams, sk, dft):
    """Band gap against width; (a) HSE06, (b) PBE. Error bands (PREREG_t1b §2): constant, symmetric, the larger
    hold-out gap error per direction; sign not used."""
    FS.use(7.0)
    band = {"100": 0.025, "110": 0.032, "110h": 0.032}
    ref = {"HSE": {"sheet": 1.415, "slab": 1.106}, "PBE": {"sheet": 0.404, "slab": 0.120}}   # A2 / A3 DFT
    XINF = 6.85
    fig, axes = plt.subplots(1, 2, figsize=(FS.width(1.0), 3.45))
    for i, (ax, func) in enumerate(zip(axes, ("HSE", "PBE"))):
        for f in fams:
            x, y = series(rows(f), lambda r: r[func]["gap"])
            ax.fill_between(x, y - band[f], y + band[f], color=COL[f], alpha=0.2, lw=0)
            ax.plot(x, y, "-o", color=COL[f], ms=2.0, lw=1.0)
            if f == "110":
                x2, y2 = series(rows(f), lambda r: r[func]["no_edge_s"]["gap"])
                ax.plot(x2, y2, ":", color=COL[f], lw=1.1)
        if func == "HSE" and sk:
            x, y = series(sk, lambda r: r["SK_HSE"]["gap"])
            ax.plot(x, y, "--", color=SKC, lw=1.0)
        if dft:
            for f, c in dft["calibration"].items():
                square(ax, c["W_A"] / 10, c["HSE_grid_gap"] if func == "HSE" else c["PBE_W2"]["gap"], COL[f], "s", 4.2)
            for f, h in dft["holdouts"].items():
                square(ax, h["W_A"] / 10, h["HSE_grid_gap"] if func == "HSE" else h["PBE_W2"]["gap"], COL[f], "D", 4.0)
        ax.plot(XINF, ref[func]["sheet"], "h", mfc=INK, mec=INK, ms=4.6, zorder=5)
        ax.plot(XINF, ref[func]["slab"], "p", mfc="white", mec=INK, mew=0.9, ms=4.6, zorder=5)
        ax.axvline(6.45, color="#bdbdbd", lw=0.6, ls=(0, (1, 1.5)), zorder=0)
        ax.set_xlim(0.8, 7.1)
        ax.set_xticks([1, 2, 3, 4, 5, 6, XINF], ["1", "2", "3", "4", "5", "6", "∞"])
        ax.set_xlabel("$W$ (nm)")
        ax.set_ylabel(f"{'HSE06' if func == 'HSE' else 'PBE'} gap (eV)")
        FS.letter(ax, "ab"[i])
        FS.style(ax)
    h = [Line2D([], [], color=COL[f], marker="o", ms=2.0, lw=1.0) for f in fams]
    l = [FLAB[f] for f in fams]
    h += [Line2D([], [], color=COL["110"], ls=":", lw=1.1), Line2D([], [], color=SKC, ls="--", lw=1.0),
          Patch(color=INK, alpha=0.2, lw=0),
          Line2D([], [], ls="", marker="s", mfc="white", mec=INK, mew=0.9, ms=4.2),
          Line2D([], [], ls="", marker="D", mfc="white", mec=INK, mew=0.9, ms=4.0),
          Line2D([], [], ls="", marker="s", mfc="white", mec=COL["100h"], mew=0.9, ms=4.2),
          Line2D([], [], ls="", marker="h", mfc=INK, mec=INK, ms=4.6),
          Line2D([], [], ls="", marker="p", mfc="white", mec=INK, mew=0.9, ms=4.6)]
    l += ["[110] bare, edge $s$ left out", "[100] bare, SK-HSE06", "hold-out gap error",
          "DFT calibration ribbon", "DFT hold-out", "[100] H (DFT only)", "sheet ($W \\to \\infty$)",
          "three-plane slab ($W \\to \\infty$)"]
    fig.legend(h, l, loc="outside lower center", ncol=3, columnspacing=0.9, handlelength=1.7)  # column-major: 4 rows
    FS.save(fig, FIG / "Fig4_Eg_W.png")


# ==================================================================================================== Fig. 5
def fig5(fams, sk):
    """(a) conduction-valley offset E_c(X) − E_c(Γ), (b) electron and (c) hole band-edge masses, HSE06."""
    FS.use(7.0)
    fig, axes = plt.subplots(1, 3, figsize=(FS.width(1.0), 2.55))
    ax = axes[0]
    for f in fams:
        x, y = series(rows(f), lambda r: r["HSE"]["valley_XG"])
        ax.plot(x, y, "-o", color=COL[f], ms=2.0, lw=1.0)
        if f == "110":
            x2, y2 = series(rows(f), lambda r: r["HSE"]["no_edge_s"]["valley_XG"])
            ax.plot(x2, y2, ":", color=COL[f], lw=1.1)
    ax.axhline(0, color="#999999", lw=0.6)
    for a in axes:
        a.axvspan(1.0, 2.5, color="#efefef", lw=0, zorder=0)
    axes[1].text(1.75, 0.5, "DFT\nbracket", transform=axes[1].get_xaxis_transform(), fontsize=8.0, color="#777777",
                 va="center", ha="center", linespacing=1.0)
    ax.set_ylabel("$E_\\mathrm{c}(\\mathrm{X}) - E_\\mathrm{c}(\\Gamma)$ (eV)")
    for ax, key, lab in ((axes[1], "m_CB", "electron mass $|m^*|$ ($m_\\mathrm{e}$)"),
                         (axes[2], "m_VB", "hole mass $|m^*|$ ($m_\\mathrm{e}$)")):
        for f in fams:
            x, y = series(rows(f), lambda r: abs(r["HSE"][key]))
            ax.plot(x, y, "-o", color=COL[f], ms=2.0, lw=1.0)
        if sk:
            x, y = series(sk, lambda r: abs(r["SK_HSE"][key]))
            ax.plot(x, y, "--", color=SKC, lw=1.0)
        ax.set_ylabel(lab)
        if key == "m_VB":
            ax.set_yscale("log")              # [110] bare n = 4 has a nearly flat VB top (|m| ≈ 38)
    for i, ax in enumerate(axes):
        ax.set_xlabel("$W$ (nm)")
        ax.set_xticks([1, 2, 3, 4, 5, 6])
        FS.letter(ax, "abc"[i])
        FS.style(ax)
    h = [Line2D([], [], color=COL[f], marker="o", ms=2.0, lw=1.0) for f in fams]
    l = [FLAB[f] for f in fams]
    h += [Line2D([], [], color=COL["110"], ls=":", lw=1.1), Line2D([], [], color=SKC, ls="--", lw=1.0)]
    l += ["[110] bare, edge $s$ left out", "[100] bare, SK-HSE06"]
    fig.legend(h, l, loc="outside lower center", ncol=3)
    FS.save(fig, FIG / "Fig5_valley_mass_W.png")


# ==================================================================================================== Fig. 7 (SM)
def fig7(fams):
    """Ballistic I_ON(W) at EOT 0.6 nm under every C_G (decision 2026-11-16; array added 2026-11-19).
    Solid: parallel plate; shaded to the isolated ribbon with a single HfO2 gate (dashed); dotted: isolated GAA;
    thick translucent: an array under one gate (periodic Laplace C_G), s = 2 nm (solid) and 7.7 nm (dotted).
    Calibration DFT points: square = parallel plate, ▲ single, ▼ GAA."""
    FS.use(8.0)
    dft = json.loads((T1 / "dft.json").read_text())["calibration"]
    dcg = {g: json.loads(p.read_text())["points"] for g in ("single", "gaa")
           if (p := T1 / "cg" / g / "dft.json").exists()}
    arr = array_data()
    fig, axes = plt.subplots(2, len(fams), figsize=(FS.width(1.0, sm=True), 4.75), sharex=True, sharey=True,
                             squeeze=False)
    for j, f in enumerate(fams):
        pts = [f] + (["100h"] if f == "100" else [])
        for i, c in enumerate(("n", "p")):
            ax = axes[i, j]
            k = f"{c} EOT 0.6"
            x, ypp = series(rows(f), lambda r: r["HSE"]["tob"][k]["I_on"])
            ax.plot(x, ypp, "-", color=COL[f], lw=1.3, zorder=3)
            s_ = {r["n"]: r["tob"][k]["I_on"] for r in cg_rows("single", f)}
            g_ = {r["n"]: r["tob"][k]["I_on"] for r in cg_rows("gaa", f)}
            ns = [r["n"] for r in rows(f)]
            if s_:
                m = np.array([n in s_ for n in ns])
                ys = np.array([s_[n] for n in ns if n in s_])
                ax.fill_between(x[m], ypp[m], ys, color=COL[f], alpha=0.13, lw=0)
                ax.plot(x[m], ys, "--", color=COL[f], lw=1.0)
            if g_:
                m = np.array([n in g_ for n in ns])
                ax.plot(x[m], [g_[n] for n in ns if n in g_], ":", color=COL[f], lw=1.1)
            if arr:
                ra = arr["families"][f]
                for s in ("s2", "sGAA"):
                    ax.plot([r["W_nm"] for r in ra], [r[s][c]["per_W"]["ballistic"] for r in ra], color=COL[f],
                            zorder=2, **ARR[s])
            for p in pts:
                if p not in dft:
                    continue
                W = dft[p]["W_A"] / 10
                square(ax, W, dft[p]["HSE"]["tob"][k]["I_on"], COL[p], "s", 4.2)
                for g, mk in (("single", "^"), ("gaa", "v")):
                    if g in dcg and p in dcg[g]:
                        square(ax, W, dcg[g][p]["tob"][k]["I_on"], COL[p], mk, 4.2)
            ax.axhline(IRDS, **IRDS_KW)
            ax.set_yscale("log")
            FS.letter(ax, "abcdef"[3 * i + j], header=FLAB[f] if i == 0 else None, header_color=COL[f])
            if i == 1:
                ax.set_xlabel("$W$ (nm)")
            if j == 0:
                ax.set_ylabel(f"{c}-FET $I_\\mathrm{{ON}}$ (mA/μm)")
            ax.set_xticks([1, 2, 3, 4, 5, 6])
            FS.style(ax)
    h = [Line2D([], [], color=INK, lw=1.3), Patch(color=INK, alpha=0.13, lw=0),
         Line2D([], [], color=INK, ls="--", lw=1.0), Line2D([], [], color=INK, ls=":", lw=1.1),
         arr_handle("s2"), arr_handle("sGAA"),
         Line2D([], [], ls="", marker="s", mfc="white", mec=INK, mew=0.9, ms=4.2),
         Line2D([], [], ls="", marker="^", mfc="white", mec=INK, mew=0.9, ms=4.2),
         Line2D([], [], ls="", marker="v", mfc="white", mec=INK, mew=0.9, ms=4.2),
         Line2D([], [], ls="", marker="s", mfc="white", mec=COL["100h"], mew=0.9, ms=4.2), irds_handle()]
    l = ["parallel-plate $C_\\mathrm{G}$", "range to the isolated single gate",
         "isolated ribbon, single gate", "isolated ribbon, gate all around",
         "array under one gate, $s$ = 2 nm", "array under one gate, $s$ = 7.7 nm",
         "DFT ribbon, parallel plate", "DFT ribbon, single gate", "DFT ribbon, gate all around",
         "[100] H, DFT point only", "IRDS, 1.979 mA/μm"]
    fig.legend(h, l, loc="outside lower center", ncol=3)   # column-major, 4 rows
    FS.save(fig, FIG / "Fig7_Ion_W.png")


# ==================================================================================================== Fig. 7 (main)
def fig7b(fams):
    """Main-text Fig. 7: per family (columns) and carrier (rows): ballistic thin, B-corrected bold, B range (2D → 3D)
    shaded, for parallel plate (solid) and the isolated single gate (dashed); the array under one gate (× B, 2D;
    thick translucent: s = 2 nm solid, 7.7 nm dotted); flagged widths ticked along the bottom axis."""
    FS.use(7.0)
    dft = json.loads((T1 / "corr" / "dft.json").read_text())["points"] if (T1 / "corr" / "dft.json").exists() else {}
    arr = array_data()
    fams = [f for f in fams if corr_rows(f)]
    YLO, YHI = 0.04, 25.0
    fig, axes = plt.subplots(2, len(fams), figsize=(FS.width(1.0), 4.3), sharex=True, sharey=True, squeeze=False)
    for j, f in enumerate(fams):
        rs = corr_rows(f)
        x = np.array([r["W_nm"] for r in rs])
        for i, c in enumerate(("n", "p")):
            ax = axes[i, j]
            for g, ls in (("parallel_plate", "-"), ("single", "--")):
                k = f"{g} EOT 0.6"
                ball = np.array([r["carriers"][c]["I_on"][k]["ballistic"] for r in rs])
                b2 = np.array([r["carriers"][c]["I_on"][k]["B2d"] for r in rs])
                b3 = np.array([r["carriers"][c]["I_on"][k]["B3d"] for r in rs])
                ax.plot(x, ball, ls, color=COL[f], lw=0.6, alpha=0.55)
                ax.fill_between(x, b2, b3, color=COL[f], alpha=0.15, lw=0)
                ax.plot(x, b2, ls, color=COL[f], lw=1.25, zorder=3)
            if arr:
                ra = arr["families"][f]
                for s in ("s2", "sGAA"):
                    ax.plot([r["W_nm"] for r in ra], [r[s][c]["per_W"]["B2d"] for r in ra], color=COL[f], zorder=2,
                            **ARR[s])
            for p in ([f] + (["100h"] if f == "100" else [])):
                if p in dft:
                    d = dft[p]["carriers"][c]["I_on"]
                    square(ax, dft[p]["W_nm"], d["parallel_plate EOT 0.6"]["B2d"], COL[p], "s")
                    square(ax, dft[p]["W_nm"], d["single EOT 0.6"]["B2d"], COL[p], "^")
            flag = np.array([r["carriers"][c]["flags"]["btbt"] for r in rs])
            if flag.any():
                ax.plot(x[flag], np.full(flag.sum(), 0.035), "|", color="#7a7a7a", ms=4.5, mew=0.8,
                        transform=ax.get_xaxis_transform(), clip_on=False)
            ax.axhline(IRDS, **IRDS_KW)
            for r in rs:                       # flat building-block band (bare [110] n = 4 VB): say so
                v = r["carriers"][c]["I_on"]["parallel_plate EOT 0.6"]["B2d"]
                if v < 0.25:
                    ax.annotate(f"flat band, α {r['carriers'][c]['alpha_3D']:.1f}", (r["W_nm"], v),
                                (r["W_nm"] + 0.6, 0.12), fontsize=8.0, color="#555555", va="center",
                                arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5, shrinkA=1, shrinkB=2))
            ax.set_yscale("log")
            ax.set_ylim(YLO, YHI)
            ax.set_xticks([1, 2, 3, 4, 5, 6])
            FS.letter(ax, "abcdef"[3 * i + j], header=FLAB[f] if i == 0 else None, header_color=COL[f])
            if i == 1:
                ax.set_xlabel("$W$ (nm)")
            if j == 0:
                ax.set_ylabel(f"{c}-FET $I_\\mathrm{{ON}}$ (mA/μm)")
            FS.style(ax)
    # legend columns (filled column-major, 4 rows): line weight | C_G | markers
    h = [Line2D([], [], color=INK, lw=0.6, alpha=0.55), Line2D([], [], color=INK, lw=1.25),
         Patch(color=INK, alpha=0.15, lw=0), irds_handle(),
         Line2D([], [], color=INK, ls="-", lw=1.0), Line2D([], [], color=INK, ls="--", lw=1.0),
         arr_handle("s2"), arr_handle("sGAA"),
         Line2D([], [], ls="", marker="s", mfc="white", mec=INK, mew=0.9, ms=3.8),
         Line2D([], [], ls="", marker="^", mfc="white", mec=INK, mew=0.9, ms=3.8),
         Line2D([], [], ls="", marker="s", mfc="white", mec=COL["100h"], mew=0.9, ms=3.8),
         Line2D([], [], ls="", marker="|", color="#7a7a7a", ms=4.5, mew=0.8)]
    l = ["ballistic", "$\\times B$ (2D, $L$ = 15 nm)", "range to bulk 3D $B$", "IRDS, 1.979 mA/μm",
         "parallel-plate $C_\\mathrm{G}$", "isolated, single gate",
         "array, $s$ = 2 nm ($\\times B$)", "array, $s$ = 7.7 nm ($\\times B$)",
         "DFT, parallel plate", "DFT, single gate", "[100] H (DFT only)",
         "band-to-band flag"]
    fig.legend(h, l, loc="outside lower center", ncol=3, columnspacing=0.9, handlelength=1.7)
    FS.save(fig, FIG / "Fig7b_Ion_W_corrected.png")


# ==================================================================================================== Fig. 8 (SM)
def fig8():
    """figures/Fig8_Ion_pitch.png from fits/t1/pitch/pitch.json (PREREG_t1b §1) and the array C_G
    (fits/t1/cg_array/array.json; PREREG_array). × B (2D) throughout."""
    FS.use(8.0)
    d = json.loads((T1 / "pitch" / "pitch.json").read_text())
    arr = array_data()
    sg = d["spacings_nm"]["sGAA"]
    cols = [("per_W", None, "per unit width"), ("fp", "s2", "per footprint, $s$ = 2 nm"),
            ("fp", "sGAA", f"per footprint, $s$ = {sg:.1f} nm")]
    gls = {"parallel_plate": "-", "single": "--", "gaa": ":"}
    fams = list(d["families"])
    fig, axes = plt.subplots(2, 3, figsize=(FS.width(1.0, sm=True), 5.0), sharex=True, sharey=True)
    lo = 0.04
    for j, (mode, sk, title) in enumerate(cols):
        geoms = ("parallel_plate", "single", "gaa") if sk != "s2" else ("parallel_plate", "single")
        for i, c in enumerate(("n", "p")):
            ax = axes[i, j]
            for f in fams:
                rs = d["families"][f]
                x = np.array([r["W_nm"] for r in rs])
                get = (lambda r, g: r["per_W"][c][g]["B2d"]) if mode == "per_W" else \
                      (lambda r, g: r["fp"][c][sk][g]["B2d"])
                y = {g: np.array([get(r, g) for r in rs]) for g in geoms}
                ax.fill_between(x, y["parallel_plate"], y["single"], color=COL[f], alpha=0.12, lw=0)
                for g in geoms:
                    ax.plot(x, y[g], gls[g], color=COL[f], lw=1.2 if g == "parallel_plate" else 1.0)
                if arr:
                    ra = arr["families"][f]
                    for s in (("s2", "sGAA") if mode == "per_W" else (sk,)):
                        ax.plot([r["W_nm"] for r in ra], [r[s][c]["per_W" if mode == "per_W" else "fp"]["B2d"]
                                                          for r in ra], color=COL[f], zorder=2, **ARR[s])
                for r in rs:                         # flat-band point below the axis: say so
                    v = get(r, "parallel_plate")
                    if v < lo:
                        ax.annotate(f"{v:.3f} (flat band)", (r["W_nm"], lo), (r["W_nm"] + 0.9, lo * 1.2),
                                    fontsize=8.0, color="#555555", va="bottom",
                                    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5))
            for p, q in d["dft"].items():
                v = q["per_W"][c] if mode == "per_W" else q["fp"][c][sk]
                square(ax, q["W_nm"], v["parallel_plate"]["B2d"], COL[p], "s", 3.8)
                square(ax, q["W_nm"], v["single"]["B2d"], COL[p], "^", 3.8)
            ax.axhline(IRDS, **IRDS_KW)
            ax.set_yscale("log")
            ax.set_ylim(lo, 15)
            ax.set_xticks([1, 2, 3, 4, 5, 6])
            FS.letter(ax, "abcdef"[3 * i + j], header=title if i == 0 else None)
            if i == 1:
                ax.set_xlabel("$W$ (nm)")
            if j == 0:
                ax.set_ylabel(f"{c}-FET $I_\\mathrm{{ON}} \\times B$ (mA/μm)")
            FS.style(ax)
    h = [Line2D([], [], color=COL[f], lw=1.2) for f in fams]
    l = [FLAB[f] for f in fams]
    h += [Line2D([], [], color=INK, ls="-", lw=1.2), Line2D([], [], color=INK, ls="--", lw=1.0),
          Line2D([], [], color=INK, ls=":", lw=1.0), arr_handle("s2"), arr_handle("sGAA"),
          Line2D([], [], ls="", marker="s", mfc="white", mec=INK, mew=0.9, ms=3.8),
          Line2D([], [], ls="", marker="^", mfc="white", mec=INK, mew=0.9, ms=3.8),
          Line2D([], [], ls="", marker="s", mfc="white", mec=COL["100h"], mew=0.9, ms=3.8), irds_handle()]
    l += ["parallel-plate $C_\\mathrm{G}$", "isolated ribbon, single gate", "isolated ribbon, GAA (not at $s$ = 2 nm)",
          "array under one gate, $s$ = 2 nm", "array under one gate, $s$ = 7.7 nm",
          "DFT ribbon, parallel plate", "DFT ribbon, single gate", "[100] H, DFT point only", "IRDS, 1.979 mA/μm"]
    fig.legend(h, l, loc="outside lower center", ncol=3)   # column-major, 4 rows
    FS.save(fig, FIG / "Fig8_Ion_pitch.png")


# ==================================================================================================== legacy outputs
# The tables and the parallel-plate-only Fig. 7 below are unchanged from the version of 2026-11-17 (backup
# scripts/t1_figs.py.pre_sn_2026-11-19.bak); they run only in the legacy modes.
def fig7_pp(fams, sk, dft):
    """figures/Fig7_Ion_W_parallel_plate.png (not in the manuscript), matplotlib default style."""
    with plt.style.context("default"):
        # ---------------------------------------------------------------- Fig 7 (partial): I_ON(W), parallel-plate C_G
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), sharey=True)
        for ax, c in zip(axes, ("n", "p")):
            for f in fams:
                for eot, ls in (("0.6", "-"), ("1.0", "--")):
                    x, y = series(rows(f), lambda r: r["HSE"]["tob"][f"{c} EOT {eot}"]["I_on"])
                    ax.plot(x, y, ls, marker="o" if eot == "0.6" else None, ms=3, color=COL[f], lw=1.3)
            if sk:
                x, y = series(sk, lambda r: r["SK_HSE"]["tob"][f"{c} EOT 0.6"]["I_on"])
                ax.plot(x, y, "-.", color=SKC, lw=1.2)
            if dft:
                for f, d in dft["calibration"].items():
                    ax.plot(d["W_A"] / 10, d["HSE"]["tob"][f"{c} EOT 0.6"]["I_on"], "s", mfc="white", mec=COL[f],
                            mew=1.4, ms=7, zorder=5)
            ax.axhline(IRDS, color="#999999", lw=0.9, ls=(0, (4, 3)))
            ax.text(ax.get_xlim()[1] if False else 5.9, IRDS * 1.02, "IRDS 1.979 mA/μm", fontsize=7.5, color="#666666",
                    ha="right", va="bottom")
            ax.set_xlabel("W (nm)")
            ax.set_title(f"{c}-FET", fontsize=10)
            style(ax)
        axes[0].set_ylabel("I_ON (mA/μm), HSE06")
        h = [plt.Line2D([], [], color=COL[f], lw=1.4) for f in fams]
        l = [f"{LAB[f]}" for f in fams]
        h += [plt.Line2D([], [], color=INK, ls="-", marker="o", ms=3), plt.Line2D([], [], color=INK, ls="--"),
              plt.Line2D([], [], color=SKC, ls="-."),
              plt.Line2D([], [], ls="", marker="s", mfc="white", mec=INK, ms=7)]
        l += ["EOT 0.6 nm", "EOT 1 nm", "[100] bare, SK-HSE (EOT 0.6)", "DFT-Wannier calibration ribbon (EOT 0.6)"]
        fig.legend(h, l, loc="lower center", ncol=4, fontsize=8, frameon=False)
        fig.suptitle("Ballistic, ideal-electrostatics upper bound; parallel-plate C_G", fontsize=9, color="#555555", y=0.995)
        fig.tight_layout(rect=(0, 0.12, 1, 0.97))
        fig.savefig(FIG / "Fig7_Ion_W_parallel_plate.png", dpi=200)
        plt.close(fig)



def tables(fams, sk, dft):
    # ---------------------------------------------------------------- table
    lines = ["# T1 width sweep: HSE06 building-block device Hamiltonians (`scripts/t1_bb_widths.py`)", "",
             "Ballistic, ideal-electrostatics upper bound; parallel-plate C_G. I_ON in mA/μm; IRDS: ≥ 1.979 at EOT 0.6.", ""]
    for f in fams:
        lines += [f"## {LAB[f]}", "",
                  "| n | W (nm) | gap HSE (eV) | direct | E_c(X) − E_c(Γ) | m_VB | m_CB | I_ON n 1.0 / 0.6 | I_ON p 1.0 / 0.6 | IRDS n / p | gap PBE |",
                  "|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in rows(f):
            m, t = r["HSE"], r["HSE"]["tob"]
            lines.append(f"| {r['n']} | {r['W_A'] / 10:.2f} | {m['gap']:.3f} | {'yes' if m['direct'] else 'no'} | "
                         f"{m['valley_XG']:+.3f} | {m['m_VB']:+.3f} | {m['m_CB']:+.3f} | "
                         f"{t['n EOT 1.0']['I_on']:.2f} / {t['n EOT 0.6']['I_on']:.2f} | "
                         f"{t['p EOT 1.0']['I_on']:.2f} / {t['p EOT 0.6']['I_on']:.2f} | "
                         f"{'yes' if t['irds']['n'] else 'no'} / {'yes' if t['irds']['p'] else 'no'} | {r['PBE']['gap']:.3f} |")
        lines.append("")
    if sk:
        lines += ["## [100] bare, SK-HSE cross-check", "", "| n | W (nm) | gap | m_VB | m_CB | I_ON n / p (EOT 0.6) |",
                  "|---|---|---|---|---|---|"]
        for r in sk:
            m = r["SK_HSE"]
            lines.append(f"| {r['n']} | {r['W_A'] / 10:.2f} | {m['gap']:.3f} | {m['m_VB']:+.3f} | {m['m_CB']:+.3f} | "
                         f"{m['tob']['n EOT 0.6']['I_on']:.2f} / {m['tob']['p EOT 0.6']['I_on']:.2f} |")
        lines.append("")
    if dft:
        lines += ["## DFT points", "", "| ribbon | variant | W (nm) | HSE grid gap | Wannier gap | I_ON n / p (EOT 0.6) |",
                  "|---|---|---|---|---|---|"]
        for f, c in dft["calibration"].items():
            lines.append(f"| {c['key']} | {c['variant']} | {c['W_A'] / 10:.2f} | {c['HSE_grid_gap']:.3f} | {c['HSE']['gap']:.3f} | "
                         f"{c['HSE']['tob']['n EOT 0.6']['I_on']:.2f} / {c['HSE']['tob']['p EOT 0.6']['I_on']:.2f} |")
        for f, h in dft["holdouts"].items():
            lines.append(f"| {h['key']} | HSE gap point | {h['W_A'] / 10:.2f} | {h['HSE_grid_gap']:.3f} | | |")
    (T1 / "t1_table.md").write_text("\n".join(lines) + "\n")
    print("wrote fits/t1/t1_table.md")


def cg_summary(fams):
    dft = json.loads((T1 / "dft.json").read_text())["calibration"]
    dcg = {g: json.loads(p.read_text())["points"] for g in ("single", "gaa")
           if (p := T1 / "cg" / g / "dft.json").exists()}
    lines = ["# T1 I_ON under the three C_G conventions (EOT 0.6 / 1 nm; mA/μm)", "",
             "Settings: `fits/PREREG_cg.md`. IRDS ≥ 1.979 at EOT 0.6. C/C_pp at EOT 0.6 in brackets.", ""]
    for f in fams:
        lines += [f"## {LAB[f]}", "", "| n | W (nm) | n: pp / single / GAA (EOT 0.6) | p: pp / single / GAA (EOT 0.6) | "
                  "n EOT 1: pp / single / GAA | p EOT 1: pp / single / GAA | IRDS n, p (pp / single / GAA) |",
                  "|---|---|---|---|---|---|---|"]
        s_, g_ = {r["n"]: r for r in cg_rows("single", f)}, {r["n"]: r for r in cg_rows("gaa", f)}
        for r in rows(f):
            if r["n"] not in s_ or r["n"] not in g_:
                continue
            t, ts, tg = r["HSE"]["tob"], s_[r["n"]]["tob"], g_[r["n"]]["tob"]
            yn = lambda d, c: "y" if d["irds"][c] else "n"   # noqa: E731
            lines.append(
                f"| {r['n']} | {r['W_A'] / 10:.2f} | {t['n EOT 0.6']['I_on']:.2f} / {ts['n EOT 0.6']['I_on']:.2f} "
                f"({ts['n EOT 0.6']['C_ratio']:.1f}×) / {tg['n EOT 0.6']['I_on']:.2f} ({tg['n EOT 0.6']['C_ratio']:.1f}×) | "
                f"{t['p EOT 0.6']['I_on']:.2f} / {ts['p EOT 0.6']['I_on']:.2f} / {tg['p EOT 0.6']['I_on']:.2f} | "
                f"{t['n EOT 1.0']['I_on']:.2f} / {ts['n EOT 1.0']['I_on']:.2f} / {tg['n EOT 1.0']['I_on']:.2f} | "
                f"{t['p EOT 1.0']['I_on']:.2f} / {ts['p EOT 1.0']['I_on']:.2f} / {tg['p EOT 1.0']['I_on']:.2f} | "
                f"{yn(t, 'n')}{yn(t, 'p')} / {yn(ts, 'n')}{yn(ts, 'p')} / {yn(tg, 'n')}{yn(tg, 'p')} |")
        lines.append("")
    if len(dcg) == 2:
        lines += ["## Calibration DFT-Wannier ribbons (relaxed W)", "",
                  "| ribbon | variant | W (nm) | n: pp / single / GAA (EOT 0.6) | p: pp / single / GAA (EOT 0.6) | "
                  "IRDS n, p (pp / single / GAA) |", "|---|---|---|---|---|---|"]
        yn = lambda d, c: "y" if d["irds"][c] else "n"   # noqa: E731
        for p, c in dft.items():
            t, ts, tg = c["HSE"]["tob"], dcg["single"][p]["tob"], dcg["gaa"][p]["tob"]
            lines.append(f"| {c['key']} | {c['variant']} | {c['W_A'] / 10:.2f} | "
                         f"{t['n EOT 0.6']['I_on']:.2f} / {ts['n EOT 0.6']['I_on']:.2f} / {tg['n EOT 0.6']['I_on']:.2f} | "
                         f"{t['p EOT 0.6']['I_on']:.2f} / {ts['p EOT 0.6']['I_on']:.2f} / {tg['p EOT 0.6']['I_on']:.2f} | "
                         f"{yn(t, 'n')}{yn(t, 'p')} / {yn(ts, 'n')}{yn(ts, 'p')} / {yn(tg, 'n')}{yn(tg, 'p')} |")
        lines.append("")
    (T1 / "cg" / "summary.md").write_text("\n".join(lines) + "\n")
    print("wrote fits/t1/cg/summary.md")


def corr_summary(fams):
    """fits/t1/corr/summary.md (PREREG_corr)."""
    dft = json.loads((T1 / "corr" / "dft.json").read_text())["points"] if (T1 / "corr" / "dft.json").exists() else {}
    L = ["# T1 corrections: Fröhlich ballisticity and tunnelling guard (`fits/PREREG_corr.md`)", "",
         "B from Sohier's 2D Fröhlich model in vacuum (main) and the bulk 3D rate (other end); L = 15 nm. "
         "R = OFF-state source–drain tunnelling ratio (transfer matrix) at λ_s = 1 / 2 / 3 nm. BTBT = drain-side "
         "interband upper bound ÷ I_OFF at λ_s = 2 nm (window = E_bF + V_DD − E_g). I_ON in mA/μm at EOT 0.6: "
         "ballistic → B-corrected → B × re-pinned (λ_s = 2 nm). IRDS ≥ 1.979.", ""]
    for f in fams:
        rs = corr_rows(f)
        if not rs:
            continue
        L += [f"## {LAB[f]}", "", "| n | W (nm) | E_g | carrier | m | α_3D | λ0 2D / 3D (nm) | B 2D / 3D | R (1/2/3 nm) | "
              "ΔV_T (mV) | window (eV) | BTBT ÷ I_OFF | I_ON pp | I_ON single | I_ON GAA | IRDS pp / single / GAA (ball, B, full) |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in rs:
            for c, x in r["carriers"].items():
                I = x["I_on"]
                tri = lambda g: (f"{I[f'{g} EOT 0.6']['ballistic']:.2f} → {I[f'{g} EOT 0.6']['B2d']:.2f} → "   # noqa: E731
                                 f"{I[f'{g} EOT 0.6']['full']:.2f}")
                yn = lambda g: "".join("y" if x["irds"][g][k] else "n" for k in ("ballistic", "B2d", "full"))  # noqa: E731
                L.append(f"| {r['n']} | {r['W_nm']:.2f} | {r['Eg']:.3f} | {c} | {x['m']:.2f} | {x['alpha_3D']:.2f} | "
                         f"{x['frohlich']['2d']['lam0_nm']:.1f} / {x['frohlich']['3d']['lam0_nm']:.1f} | "
                         f"{x['frohlich']['2d']['B']:.2f} / {x['frohlich']['3d']['B']:.2f} | "
                         f"{x['sdt']['1']['R_tm']:.2f} / {x['sdt']['2']['R_tm']:.2f} / {x['sdt']['3']['R_tm']:.2f} | "
                         f"{x['dVT_mV']:.1f} | {x['btbt']['2']['window_eV']:+.3f} | {x['btbt_bound_over_IOFF_lam2']:.1e} | "
                         f"{tri('parallel_plate')} | {tri('single')} | {tri('gaa')} | "
                         f"{yn('parallel_plate')} / {yn('single')} / {yn('gaa')} |")
        L.append("")
    if dft:
        L += ["## Calibration DFT-Wannier ribbons", "", "| ribbon | carrier | m | B 2D / 3D | R (2 nm) | window (eV) | "
              "BTBT ÷ I_OFF | I_ON pp | I_ON single | I_ON GAA |", "|---|---|---|---|---|---|---|---|---|---|"]
        for p, d in dft.items():
            for c, x in d["carriers"].items():
                I = x["I_on"]
                tri = lambda g: (f"{I[f'{g} EOT 0.6']['ballistic']:.2f} → {I[f'{g} EOT 0.6']['B2d']:.2f} → "   # noqa: E731
                                 f"{I[f'{g} EOT 0.6']['full']:.2f}")
                L.append(f"| {d['key']} ({d['variant']}) | {c} | {x['m']:.2f} | {x['frohlich']['2d']['B']:.2f} / "
                         f"{x['frohlich']['3d']['B']:.2f} | {x['sdt']['2']['R_tm']:.2f} | {x['btbt']['2']['window_eV']:+.3f} | "
                         f"{x['btbt_bound_over_IOFF_lam2']:.1e} | {tri('parallel_plate')} | {tri('single')} | {tri('gaa')} |")
        L.append("")
    (T1 / "corr" / "summary.md").write_text("\n".join(L) + "\n")

    print("wrote fits/t1/corr/summary.md")


def main():
    """Legacy default mode: tables and the T1 figures (as on 2026-11-17, plus the re-rendered figures)."""
    fams = [f for f in ("100", "110", "110h") if rows(f)]
    sk = sk_rows()
    dft = json.loads((T1 / "dft.json").read_text()) if (T1 / "dft.json").exists() else None
    fig7_pp(fams, sk, dft)
    fig4(fams, sk, dft)
    fig5(fams, sk)
    tables(fams, sk, dft)
    if all((T1 / "cg" / g).exists() for g in ("single", "gaa")):
        fig7(fams)
        cg_summary(fams)


def figs():
    """Manuscript figures only (Fig. 4, 5, 7 main, S3, S5). Reads fits/t1; writes only figures/."""
    fams = [f for f in ("100", "110", "110h") if rows(f)]
    sk = sk_rows()
    dft = json.loads((T1 / "dft.json").read_text()) if (T1 / "dft.json").exists() else None
    fig4(fams, sk, dft)
    fig5(fams, sk)
    fig7(fams)
    fig7b(fams)
    fig8()


if __name__ == "__main__":
    if sys.argv[1:] == ["figs"]:
        figs()
    elif sys.argv[1:] == ["corr"]:
        corr_summary(["100", "110", "110h"])
        fig7b(["100", "110", "110h"])
    elif sys.argv[1:] == ["pitch"]:
        fig8()
    else:
        main()
