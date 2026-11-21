"""Main-text Fig. 2: building-block ribbons against the PBE-Wannier hold-outs (figures/Fig3b_holdouts_building_blocks.png).

    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/fig_holdouts_bb.py

The figure used to be drawn inside scripts/gate_bb_holdouts.py (`run`), which is checksum-frozen pre-registration code
(fits/PREREG_bb_holdouts.code.md5), so it is redrawn here instead. This script imports that file read-only and calls
only its band-structure helpers (`build`, `bands_and_metrics`), exactly as `run` does for the plot: the PBE-Wannier W2
hold-out Hamiltonian and the registered building-block recipe ("BB PBE": W1 edges, centre-row alignment, 20 Å
cutoff) on the 601-point Γ–X path. The band structures are not stored, so they are recomputed, as
`gate_bcd_sk_holdouts.py figs` does for its figures; no ToB is run and nothing is written in fits/. As a check, the
recomputed gap, valley offset and masses must match the stored fits/bb/gate_bb.json to 1e-9.
Sized for 0.95 \\linewidth of the sn-jnl text width (scripts/figstyle.py).
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as FS  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, [str(path)]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


GB = _load("gate_bb_holdouts", ROOT / "scripts" / "gate_bb_holdouts.py")      # frozen; read-only use
W = GB.W
NAME = {"r110_8": "[110] $n$ = 8", "r100_6": "[100] $n$ = 6"}
KEYS = ("gap", "valley_XG", "m_VB", "m_CB")


def bands():
    res = json.loads((ROOT / "fits" / "bb" / "gate_bb.json").read_text())
    kf = np.c_[np.linspace(0, 0.5, GB.NK_PATH), np.zeros((GB.NK_PATH, 2))]
    plot = {}
    for kh, kc in GB.HOLD.items():
        st = GB.MBH.geometry(kh)
        nv = 3 * st.species.count("N")
        Lx = st.cell[0, 0]
        wfmap = json.loads((GB.M5.TWIN / "w90" / kh / "wf_map.json").read_text())["B"]
        hr_w = W.read_hr(GB.M5.TWIN / "w90" / kh / "W2" / f"{kh}_W2_hr.dat")
        Ew, ref = GB.bands_and_metrics(hr_w, [(m["atom"], m["orbital"]) for m in wfmap], Lx, nv, kf)
        hr_b, wf_b, _, _ = GB.build(GB.MBH, kh, kc, "PBE", var=None, align="centre")
        E, m = GB.bands_and_metrics(hr_b, wf_b, Lx, nv, kf)
        stored = res["ribbons"][kh]
        for k in KEYS:                                   # same data as the registered test
            assert abs(ref[k] - stored["wannier"][k]) < 1e-9, (kh, "W2", k, ref[k], stored["wannier"][k])
            assert abs(m[k] - stored["BB PBE"]["metrics"][k]) < 1e-9, (kh, "BB", k, m[k], stored["BB PBE"]["metrics"][k])
        print(f"{kh}: W2 gap {ref['gap']:.4f}, BB gap {m['gap']:.4f} eV (match gate_bb.json)")
        plot[kh] = (Ew - Ew[:, nv - 1].max(), E - E[:, nv - 1].max(), ref["gap"])
    return plot


def draw(plot):
    FS.use(7.5)
    REF, BBC = "#8a8a8a", "#c2410c"
    fig, axes = plt.subplots(1, 2, figsize=(FS.width(0.95), 2.75))
    t = np.linspace(0, 0.5, GB.NK_PATH)
    for i, (ax, kh) in enumerate(zip(axes, ("r110_8", "r100_6"))):
        Ew, E, g = plot[kh]
        ax.plot(t, Ew, color=REF, lw=1.3)
        ax.plot(t, E, color=BBC, lw=0.75, ls=(0, (3.5, 1.8)))
        ax.set_ylim(-1.5, g + 1.5)
        ax.set_xlim(0, 0.5)
        ax.set_xticks([0, 0.5], ["Γ", "X"])
        ax.set_ylabel("$E - E_\\mathrm{VBM}$ (eV)")
        ax.axhline(0, color="#bbbbbb", lw=0.5, zorder=0)
        FS.letter(ax, "ab"[i], NAME[kh])
    h = [Line2D([], [], color=REF, lw=1.3), Line2D([], [], color=BBC, lw=0.75, ls=(0, (3.5, 1.8)))]
    fig.legend(h, ["PBE-Wannier hold-out (W2)", "building blocks (W1 edges + sheet interior)"],
               loc="outside lower center", ncol=2)
    FS.save(fig, ROOT / "figures" / "Fig3b_holdouts_building_blocks.png")


if __name__ == "__main__":
    draw(bands())
