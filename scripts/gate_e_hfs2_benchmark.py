"""Claims gate (e): the ToB code against the HfS2 zigzag benchmark (Matić & Poljak, Electronics 13, 1048, 2024).

Conventions and targets: 00_Master_Notebook.md, "Pre-registered device conventions" (signed off 2026-11-11, frozen
before this run). Report, never tune: a mismatch beyond ±10 % is flagged LR10 with its suspected cause, and no
convention changes.

Inputs: benchmark/quden/materials/MATERIALS_QUDEN/HfS2XNR/<N>HfS2XNR_hr.dat (QUDEN_MATERIALS, CC BY-NC-SA 4.0,
Poljak & Matić) and QUDEN's own parameters for HfS2XNR (src/QUDEN25_06/material_selection.m, transcribed below).
QUDEN_MATERIALS is not redistributed with this repository: obtain it from its authors' distribution and place the
HfS2XNR files at the path above before running. The --plot-only mode needs only the stored gate_e.json.

    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/gate_e_hfs2_benchmark.py
Writes fits/gate_e/gate_e.json and figures/Fig9_HfS2_benchmark.png.
    ... gate_e_hfs2_benchmark.py --plot-only   # redraw the figure from the stored gate_e.json (no ToB, no write in fits/)
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

from scnfet import tob, w90 as W

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as FS  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
HR = ROOT / "benchmark/quden/materials/MATERIALS_QUDEN/HfS2XNR"
OUT = ROOT / "fits" / "gate_e"
# QUDEN material_selection.m, case 'HfS2XNR' (index = N; N = 7, 9 not provided)
A_TRANSPORT = 3.62964875e-10                                       # m, period along the ribbon
B_QUDEN = {2: 13.60962775, 3: 19.97545273, 4: 26.25095712, 5: 32.53698781, 6: 38.826392317, 8: 51.39819749,
           10: 63.973507214}                                       # Å, QUDEN's ribbon width = Matić's quoted W
EG_QUDEN = {2: 1.0119454761, 3: 1.0008983932, 4: 1.0076835848, 5: 0.9979979553, 6: 0.9839042448, 8: 0.9838726618,
            10: 0.9719015604}
CB_QUDEN = {2: 31, 3: 43, 4: 55, 5: 67, 6: 79, 8: 103, 10: 127}  # 1-based index of the lowest conduction band
NCELL, NADD = 6, 2                                                 # atoms per N, extra atoms (edge)
# Registered targets (read off the paper 2026-11-11), EOT 1 nm, mA/μm
TARGETS = {("n", 2): 1.66, ("n", 3): 1.74, ("n", 10): 1.14, ("p", 2): 1.18, ("p", 10): 1.5}
P_NARROW = 1.12                                                    # "pFET … 1.12 at narrow widths" (descriptive)
BAND = 0.10
NK, NK2, CHUNK = 2000, 4000, 200
SPAN, NSWEEP = 4.0, 381                                            # search range only (as T1)


def bands(hr, nk):
    k = np.arange(nk) / nk
    kf = np.c_[np.zeros(nk), k, np.zeros(nk)]                     # QUDEN ribbons run along y (R = (0, R, 0))
    return np.vstack([np.linalg.eigvalsh(W.hk_frac(hr, kf[i:i + CHUNK])) for i in range(0, nk, CHUNK)])


def device(E, nv, Wm, nat, carrier, eot=1e-9):
    ch = tob.Channel(E=np.sort(E, axis=1), a=A_TRANSPORT, n_val=nv)
    f = tob.ToBFET(ch, W=Wm, EOT=eot, carrier=carrier, n_atoms_per_cell=nat)
    m = f.metrics(span=SPAN, n_sweep=NSWEEP)
    return dict(I_on=float(m["I_on"]), I_off=float(m["I_off"]), SS=float(m["SS_mV_dec"]), v_inj=float(m["v_inj"]),
                n_SD_cm2=float(f.n_SD / Wm * 1e-4), Eg=float(ch.Eg))


def main():
    if not HR.is_dir():
        sys.exit(f"{HR.relative_to(ROOT)} not found. QUDEN_MATERIALS (Poljak & Matić, CC BY-NC-SA 4.0) is not "
                 "redistributed here: download it from its authors' distribution and place the HfS2XNR "
                 "Hamiltonians there (see README.md, 'HfS2 benchmark').")
    OUT.mkdir(parents=True, exist_ok=True)
    res = dict(source="QUDEN_MATERIALS 2025_01 HfS2XNR (Poljak & Matić, CC BY-NC-SA 4.0)",
               w_precheck=("No structure files in QUDEN or QUDEN_MATERIALS: our-rule W cannot be computed. By the "
                           "pre-registered branch the benchmark is judged on Matić's W (= QUDEN b), which validates "
                           "the code, not the width definition."),
               conventions="notebook 'Pre-registered device conventions' (2026-11-11): alpha_G 1, alpha_D 0, "
                           "C_G = eps0*3.9*W/EOT, EOT 1 nm, VDD 0.7 V, I_OFF 10 nA/um, 300 K, doping 0.01 per atom",
               targets={f"{c} N={n}": v for (c, n), v in TARGETS.items()}, band=BAND, widths={})
    for N in sorted(B_QUDEN):
        t0 = time.time()
        hr = W.read_hr(HR / f"{N}HfS2XNR_hr.dat")
        assert hr.num_wann == 22 * N + 6, f"N={N}: {hr.num_wann} WFs"     # Hf d5 + 2 S p3 per f.u.; 2 edge S (p3)
        assert not hr.R[:, [0, 2]].any(), "expected R along y only"
        nv = CB_QUDEN[N] - 1
        assert nv == 3 * (4 * N + NADD), "valence count ≠ S p bands"
        nat = NCELL * N + NADD                                       # channel atoms (the 2 extra are S: 3 WFs each)
        Wm = B_QUDEN[N] * 1e-10
        E = bands(hr, NK)
        Eg = float(E[:, nv:].min() - E[:, :nv].max())
        row = dict(N=N, W_nm=B_QUDEN[N] / 10, num_wann=hr.num_wann, n_val=nv, n_atoms=nat, Eg=Eg,
                   Eg_quden=EG_QUDEN[N], dEg_meV=1e3 * (Eg - EG_QUDEN[N]))
        E2 = bands(hr, NK2)
        for c in ("n", "p"):
            d = device(E, nv, Wm, nat, c)
            d2 = device(E2, nv, Wm, nat, c)
            d["k4000_rel"] = d2["I_on"] / d["I_on"] - 1
            d["eot0.6"] = device(E, nv, Wm, nat, c, eot=0.6e-9)["I_on"]       # descriptive
            tgt = TARGETS.get((c, N))
            if tgt is not None:
                d["target"] = tgt
                d["rel"] = d["I_on"] / tgt - 1
                d["ok"] = bool(abs(d["rel"]) <= BAND)
            row[c] = d
        row["seconds"] = time.time() - t0
        res["widths"][N] = row
        print(f"N {N:2d} W {B_QUDEN[N] / 10:.2f} nm: Eg {Eg:.4f} (QUDEN {EG_QUDEN[N]:.4f}, {row['dEg_meV']:+.1f} meV) | "
              + " | ".join(f"{c}: I_ON {row[c]['I_on']:.3f}" + (f" vs {row[c]['target']:.2f} ({row[c]['rel']:+.1%}) "
                                                                 f"{'ok' if row[c]['ok'] else 'OUT'}" if 'target' in row[c] else "")
                           + f" [k4000 {row[c]['k4000_rel']:+.1e}; SS {row[c]['SS']:.1f}; n_SD {row[c]['n_SD_cm2']:.2e} cm⁻²]"
                           for c in ("n", "p")) + f" | {row['seconds']:.0f} s", flush=True)
    pn = [res["widths"][N]["p"]["I_on"] for N in (2, 3, 4)]
    res["p_narrow_descriptive"] = dict(min_N2_4=min(pn), paper=P_NARROW, rel=min(pn) / P_NARROW - 1)
    checks = [res["widths"][N][c] for (c, N) in TARGETS]
    res["all_within_band"] = bool(all(x["ok"] for x in checks))
    res["max_abs_rel"] = float(max(abs(x["rel"]) for x in checks))
    (OUT / "gate_e.json").write_text(json.dumps(res, indent=1, default=float))
    print(f"gate (e): {'all registered points within ±10 %' if res['all_within_band'] else 'LR10: outside ±10 %'}; "
          f"max |rel| {res['max_abs_rel']:.1%}; p narrow (descriptive) min N=2–4 {min(pn):.3f} vs {P_NARROW}")
    figure(res)


def figure(res):
    """Main-text Fig. 3, at 0.6 \\linewidth (re-rendered 2026-11-19 for the sn-jnl width; no in-image title)."""
    FS.use(7.5)
    COL = {"n": "#0072B2", "p": "#D55E00"}
    fig, ax = plt_subplots(figsize=(FS.width(0.6), 2.35))
    Ns = sorted(res["widths"])
    x = [res["widths"][N]["W_nm"] for N in Ns]
    for c in ("n", "p"):
        y = [res["widths"][N][c]["I_on"] for N in Ns]
        ax.plot(x, y, "-o", color=COL[c], ms=2.6, lw=1.0, label=f"{c}-FET, this work")
        for (cc, N), v in TARGETS.items():
            if cc == c:
                xw = res["widths"][N]["W_nm"]
                ax.errorbar(xw, v, yerr=BAND * v, fmt="s", mfc="white", mec=COL[c], ecolor=COL[c], ms=4.2, capsize=2,
                            mew=0.9, elinewidth=0.8, capthick=0.8, zorder=5)
    ax.plot([], [], "s", mfc="white", mec="#333333", mew=0.9, ms=4.2, label="Matić and Poljak, ±10 %")
    ax.set_xlabel("$W$ (nm)")
    ax.set_ylabel("$I_\\mathrm{ON}$ (mA/μm)")
    FS.style(ax)
    ax.legend(loc="upper right", handlelength=1.6)
    FS.save(fig, ROOT / "figures" / "Fig9_HfS2_benchmark.png")


def plt_subplots(**kw):
    import matplotlib.pyplot as plt
    return plt.subplots(**kw)


def plot_only():
    res = json.loads((OUT / "gate_e.json").read_text())
    res["widths"] = {int(k): v for k, v in res["widths"].items()}       # JSON keys are strings
    figure(res)


if __name__ == "__main__":
    plot_only() if "--plot-only" in sys.argv else main()
