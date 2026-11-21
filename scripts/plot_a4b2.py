"""A4b-2 figure: the HSE06 device Hamiltonians of the four calibration ribbons (the chosen Wannier variant per ribbon)
on a fine Γ–X̄ path, with the HSE06 eigenvalues on the IBZ grid as points.

    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/plot_a4b2.py      → figures/A4b2_HSE_device_Hamiltonians.png
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

RIBBONS = [("r100_3", "[100] n = 3, bare", "W3"), ("r100_3h", "[100] n = 3, H", "W2"),
           ("r110_4", "[110] n = 4, bare", "W1"), ("r110_4h", "[110] n = 4, H", "W2")]
WIN = json.loads((A4B / "a4b2_windows.json").read_text())
INFO = json.loads((A4B / "a4b_inputs.json").read_text())
C_DFT, C_W90, C_DEV, C_NM = "#8a8f98", "#2f6db5", "#2f6db5", "#c2410c"

fig, axes = plt.subplots(1, 4, figsize=(12.5, 4.6), sharey=True)
for ax, (key, title, v) in zip(axes, RIBBONS):
    w = WIN[key]
    st = PA.MB.geometry(key)
    Lx = st.cell[0, 0]
    kh, eh, _ = PA.P.bands(PA.hout(f"hse_{key}"))
    th = kh[:, 0] / (np.pi / Lx)
    hr = W.read_hr(PA.hdir(key, v) / f"{key}_H{v}_hr.dat")
    tf = np.linspace(0, 1, 401)
    Ew = np.linalg.eigvalsh(W.hk_1d(hr, tf * 0.5, axis=0))
    vbm, cbm = w["vbm"], w["cbm"]
    ne = INFO[key]["nexcl"]
    ax.axhspan(-0.5, cbm - vbm + 0.5, color=C_DEV, alpha=0.07, lw=0)
    ax.plot(tf, Ew - vbm, "-", color=C_W90, lw=1.1)
    ax.plot(np.repeat(th[:, None], eh.shape[1] - ne, 1), eh[:, ne:] - vbm, "o", ms=3.0, mfc="none", mec=C_DFT, mew=0.8)
    fn = w["nonmodel"]["0.25"]
    if fn is not None and fn["E_eV"] - vbm < 3.0:
        ax.axhline(fn["E_eV"] - vbm, color=C_NM, lw=0.9, ls=":")
    ax.set_title(f"{title}: {v}", fontsize=10.5)
    ax.text(0.97, 0.40, f"HSE gap {w['gap_grid']:.2f} eV", transform=ax.transAxes, fontsize=8.5, va="bottom", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#d0d0d0", lw=0.6))
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 1], ["Γ", "X̄"])
    ax.grid(axis="y", color="#ebebec", lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
axes[0].set_ylim(-1.5, 3.4)
axes[0].set_ylabel("E − VBM (eV, HSE06)")
handles = [plt.Line2D([], [], ls="none", marker="o", mfc="none", mec=C_DFT, label="HSE06 eigenvalues (IBZ grid)"),
           plt.Line2D([], [], color=C_W90, lw=1.1, label="Wannier90 device Hamiltonian"),
           plt.Line2D([], [], color=C_NM, lw=0.9, ls=":", label="first state flagged by H3 (s_out ≥ 0.25)"),
           plt.Rectangle((0, 0), 1, 1, color=C_DEV, alpha=0.12, label="device window")]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=(0, 0.07, 1, 1))
outp = ROOT / "figures" / "A4b2_HSE_device_Hamiltonians.png"
fig.savefig(outp, dpi=200)
print(outp)
