"""A4b-1 figure: NC PBE bands of the four calibration ribbons against the W2 Wannier bands (the chosen recipe).

    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/plot_a4b1.py      → figures/A4b1_W90_vs_DFT.png

Reads dft/A4b/outputs/bands_<r>.out, a4b_windows.json, a4b_criteria.json and w90/<r>/W2/<r>_W2_hr.dat.
"""
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt      # noqa: E402
import numpy as np                   # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
A4B = ROOT / "dft" / "A4b"
spec = importlib.util.spec_from_file_location("parse_a4b", A4B / "parse_a4b.py")
PA = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PA)
from scnfet import w90 as W          # noqa: E402

RIBBONS = [("r100_3", "[100] n = 3, bare (15 Å)"), ("r100_3h", "[100] n = 3, H (20 Å)"),
           ("r110_4", "[110] n = 4, bare (20 Å)"), ("r110_4h", "[110] n = 4, H (20 Å)")]
WIN = json.loads((A4B / "a4b_windows.json").read_text())
CRIT = json.loads((A4B / "a4b_criteria.json").read_text())
C_DFT, C_W90, C_DEV = "#8a8f98", "#d0661f", "#2f6db5"

fig, axes = plt.subplots(1, 4, figsize=(12.5, 4.6), sharey=True)
for ax, (key, title) in zip(axes, RIBBONS):
    win = WIN[key]
    st = PA.MB.geometry(key)
    kp, ep, _ = PA.P.bands(PA.out(f"bands_{key}"))
    Lx = st.cell[0, 0]
    kx = kp[:, 0]
    t = kx / kx.max()                                             # Γ = 0, X̄ = 1
    hr = W.read_hr(A4B / "w90" / key / "W2" / f"{key}_W2_hr.dat")
    tf = np.linspace(0, 1, 301)
    Ew = np.linalg.eigvalsh(W.hk_1d(hr, tf * kx.max() * Lx / (2 * np.pi), axis=0))
    vbm, cbm = win["vbm"], win["cbm"]
    v2 = win["variants"]["W2"]
    ax.axhspan(-0.5, cbm - vbm + 0.5, color=C_DEV, alpha=0.07, lw=0)
    ax.axhline(v2["froz"] - vbm, color=C_W90, lw=1.0, ls="--")
    nm = win.get("first_nonmodel", {})
    ax.plot(t, ep - vbm, "o", ms=2.6, mfc="none", mec=C_DFT, mew=0.8)
    ax.plot(tf, Ew - vbm, "-", color=C_W90, lw=1.1)
    c = CRIT[key]["W2"]
    ax.set_title(title, fontsize=10.5)
    ax.text(0.97, 0.40, f"device max {c['P2']['max_meV']:.1f} meV\nfrozen RMS {c['P3']['rms_meV']:.2f} meV\n"
            f"gap {c['P4']['d_meV']:+.2f} meV", transform=ax.transAxes, fontsize=8.5, va="bottom", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#d0d0d0", lw=0.6))
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 1], ["Γ", "X̄"])
    ax.grid(axis="y", color="#ebebec", lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
axes[0].set_ylim(-1.5, 2.6)
axes[0].set_ylabel("E − VBM (eV)")
handles = [plt.Line2D([], [], ls="none", marker="o", mfc="none", mec=C_DFT, label="NC PBE (Quantum ESPRESSO, 31 k)"),
           plt.Line2D([], [], color=C_W90, lw=1.1, label="Wannier90, recipe W2"),
           plt.Line2D([], [], color=C_W90, lw=1.0, ls="--", label="W2 frozen-window top"),
           plt.Rectangle((0, 0), 1, 1, color=C_DEV, alpha=0.12, label="device window (VBM − 0.5 to CBM + 0.5 eV)")]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=(0, 0.07, 1, 1))
outp = ROOT / "figures" / "A4b1_W90_vs_DFT.png"
fig.savefig(outp, dpi=200)
print(outp)
