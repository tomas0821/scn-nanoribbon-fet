"""Plan Fig. 6: transfer characteristics I_D(V_GS), n- and p-FET, one ribbon per family at the hold-out widths.

    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/t1_fig6.py
Registered conventions only (notebook "Pre-registered device conventions"): ballistic ToB, ideal electrostatics,
parallel-plate C_G, EOT 0.6 nm, I_OFF pinned at V_GS = 0, V_DD 0.7 V. Bands from the T1 cache (fits/t1/cache/).
Check: the curve at V_GS = V_DD reproduces the stored I_ON (fits/t1/rows) within 0.5 %.
Writes figures/Fig6_transfer.png and fits/t1/fig6.json.

    ... t1_fig6.py --plot-only    # redraw figures/Fig6_transfer.png from the stored fits/t1/fig6.json (no ToB, no write
                                  # in fits/); print size for the Supplementary Material (scripts/figstyle.py)
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import sys  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as FS  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
T1 = ROOT / "fits" / "t1"
PICK = {"100": 6, "110": 8, "110h": 7}          # W ≈ 2.3–2.4 nm, the hold-out widths
COL = {"100": "#0072B2", "110": "#D55E00", "110h": "#009E73"}
LAB = {"100": "[100] bare", "110": "[110] bare", "110h": "[110] H (in-sample)"}
EOT, VDD, IOFF, IRDS = 0.6, 0.7, 10e-9 / 1e-6, 1.979
A_TO_MA_UM = 1e3 / 1e6                            # A/m → mA/μm


def main():
    from scnfet import tob
    out = {}
    vrel = np.linspace(-0.3, VDD, 106)
    for f, n in PICK.items():
        z = np.load(T1 / "cache" / f"{f}_n{n}.npz")
        row = json.loads((T1 / "rows" / f"{f}_n{n}.json").read_text())
        ch = tob.Channel(E=np.sort(z["E"], axis=1), a=float(z["Lx"]) * 1e-10, n_val=int(z["nv"]))
        out[f] = dict(n=n, W_nm=float(z["W_A"]) / 10)
        for c in ("n", "p"):
            s = 1.0 if c == "n" else -1.0
            fet = tob.ToBFET(ch, W=float(z["W_A"]) * 1e-10, EOT=EOT * 1e-9, carrier=c, n_atoms_per_cell=int(z["n_host"]))
            m = fet.metrics(VDD=VDD, I_OFF=IOFF, span=4.0, n_sweep=381)
            ref = row["HSE"]["tob"][f"{c} EOT 0.6"]["I_on"]
            assert abs(m["I_on"] / ref - 1) < 5e-3, (f, c, m["I_on"], ref)
            curves = {}
            for vds in (VDD, 0.05):
                I = np.abs(fet.transfer(s * vrel, s * vds, shift=m["shift"])) / fet.W * A_TO_MA_UM
                curves[f"{vds}"] = I.tolist()
            out[f][c] = dict(I_on=m["I_on"], I_on_rows=ref, SS=m["SS_mV_dec"], V_GS=vrel.tolist(), I=curves)
    (T1 / "fig6.json").write_text(json.dumps(out, indent=1))
    for f in PICK:
        print(f, out[f]["W_nm"], {c: (round(out[f][c]["I_on"], 3), round(out[f][c]["SS"], 1)) for c in "np"})
    draw(out)


def draw(out):
    """(a, b) log and (c, d) linear transfer characteristics, n-FET left, p-FET right (re-rendered 2026-11-19)."""
    FS.use(8.0)
    fig, AX = plt.subplots(2, 2, figsize=(FS.width(1.0, sm=True), 4.9), sharex=True)
    for f in PICK:
        for j, c in enumerate(("n", "p")):
            d = out[f][c]
            v = np.array(d["V_GS"])
            for vds, ls, lw in (("0.7", "-", 1.3), ("0.05", "--", 0.9)):
                I = np.array(d["I"][vds])
                AX[0, j].semilogy(v, I, ls, color=COL[f], lw=lw)
                AX[1, j].plot(v, I, ls, color=COL[f], lw=lw)
    for j, c in enumerate(("n", "p")):
        ax = AX[0, j]
        ax.axhline(IOFF * A_TO_MA_UM, color="#999999", lw=0.7, ls=":")
        ax.axhline(IRDS, color="#8c8c8c", lw=0.8, ls=(0, (4, 2.5)))
        ax.text(-0.29, IRDS * 1.6, "IRDS $I_\\mathrm{ON}$, 1.979 mA/μm", fontsize=8.0, color=FS.MUTED, va="bottom")
        ax.text(0.68, IOFF * A_TO_MA_UM * 1.8, "$I_\\mathrm{OFF}$ = 10 nA/μm (pinned)", fontsize=8.0, color=FS.MUTED,
                va="bottom", ha="right")
        ax.set_ylim(3e-11, 60)
        AX[1, j].set_xlabel("$V_\\mathrm{GS}$ (V)" if c == "n" else "$|V_\\mathrm{GS}|$ (V)")
        FS.letter(ax, "ab"[j], header=f"{c}-FET")
        FS.letter(AX[1, j], "cd"[j])
        AX[1, j].axhline(IRDS, color="#8c8c8c", lw=0.8, ls=(0, (4, 2.5)))
        AX[1, j].set_ylim(0, 4.2)
        for a in AX[:, j]:
            a.axvline(VDD, color="#cccccc", lw=0.6, zorder=0)
            FS.style(a)
    AX[0, 0].set_ylabel("$|I_\\mathrm{D}|$ (mA/μm)")
    AX[1, 0].set_ylabel("$|I_\\mathrm{D}|$ (mA/μm)")
    h = [Line2D([], [], color=COL[f], lw=1.3) for f in PICK]
    l = [f"{LAB[f]}, $W$ = {out[f]['W_nm']:.2f} nm" for f in PICK]
    h += [Line2D([], [], color="#333333", lw=1.3), Line2D([], [], color="#333333", lw=0.9, ls="--")]
    l += ["$|V_\\mathrm{DS}|$ = 0.7 V", "$|V_\\mathrm{DS}|$ = 0.05 V"]
    fig.legend(h, l, loc="outside lower center", ncol=3)
    FS.save(fig, ROOT / "figures" / "Fig6_transfer.png")


def plot_only():
    draw(json.loads((T1 / "fig6.json").read_text()))


if __name__ == "__main__":
    plot_only() if "--plot-only" in sys.argv else main()
