"""Gate (a) trade-off of every bulk HSE06 SK variant tried on 2026-11-11 (LR1 tests): ±1 eV RMS against the worst
transport-mass error. Reads fits/A1_bulk/sk_fits.json, sk_fit_extensions.json, sk_fit_s.json, sk_fit_t3.json.
Writes figures/A1_SK_gate_tradeoff.png and prints the table.

Run from the project root: python3 scripts/a1_sk_tradeoff_plot.py   (needs matplotlib: ~/.venvs/scnfet/bin/python)
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
F = ROOT / "fits" / "A1_bulk"
FIG = ROOT / "figures" / "A1_SK_gate_tradeoff.png"

FAMILY = {  # variant: (family, weighting, free parameters)
    "H": ("13-parameter (plan)", "balanced", 12), "Hn": ("13-parameter (plan)", "near-gap", 12),
    "Hm": ("13-parameter (plan)", "mass ×3", 12), "H2": ("13-parameter (plan)", "balanced", 9),
    "pp2": ("longer two-centre shells", "balanced", 14), "pd2": ("longer two-centre shells", "balanced", 14),
    "dd3": ("longer two-centre shells", "balanced", 15), "all": ("longer two-centre shells", "balanced", 19),
    "Hn_all": ("longer two-centre shells", "near-gap", 19),
    "S": ("+ Sc 4s", "balanced", 16), "Sn": ("+ Sc 4s", "near-gap", 16),
    "T": ("+ three-centre", "balanced", 15), "Tn": ("+ three-centre", "near-gap", 15), "Tm": ("+ three-centre", "mass ×10", 15),
    "TS": ("+ three-centre + Sc 4s", "balanced", 19), "TSn": ("+ three-centre + Sc 4s", "near-gap", 19),
    "TSm5": ("+ three-centre + Sc 4s", "mass ×5", 19), "TSm10": ("+ three-centre + Sc 4s", "mass ×10", 19),
}


OFFSET = {"pp2": (1.0, 1.0, "left"), "H": (1.2, -0.2, "left"), "pd2": (-1.2, -0.9, "right"), "Hm": (-1.2, -1.6, "right"),
          "dd3": (-1.2, 0.9, "right"), "Sn": (-1.2, 0.9, "right"), "TS": (1.2, 0.9, "left"), "Tn": (-1.2, -1.4, "right"),
          "Hn_all": (0.6, 1.3, "left"), "TSn": (-1.2, 0.0, "right"), "Hn": (1.2, 0.0, "left"), "T": (1.2, 0.0, "left"),
          "TSm5": (1.2, 0.0, "left"), "TSm10": (-1.2, 0.0, "right"), "Tm": (1.2, 0.0, "left"), "all": (1.2, 0.0, "left"),
          "S": (1.2, 0.0, "left"), "H2": (-1.2, 0.0, "right")}


def load():
    res = {}
    for f in ("sk_fits.json", "sk_fit_extensions.json", "sk_fit_s.json", "sk_fit_t3.json"):
        d = json.loads((F / f).read_text())
        res.update(d.get("results", d))
    rows = []
    for v, (fam, wt, npar) in FAMILY.items():
        if v not in res:
            continue
        s = res[v]["score"]
        rms = max(s["path"]["rms_pm1eV_meV"], s["grid"]["rms_pm1eV_meV"])
        m = {k: 100 * x["rel"] for k, x in s["masses"].items()}
        worst = max(m, key=lambda k: abs(m[k]))
        rows.append(dict(v=v, fam=fam, wt=wt, npar=npar, rms=rms, worst=abs(m[worst]), worst_name=worst, m=m,
                         gx=s["gaps"]["GX"]["sk"], g12=s["gaps"]["G12"]["sk"] - s["gaps"]["G12"]["ref"]))
    return rows


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = load()
    for r in rows:
        print(f"{r['v']:7s} {r['fam']:26s} {r['wt']:9s} {r['npar']:2d}  RMS {r['rms']:5.1f}  worst {r['worst_name']} {r['worst']:5.1f} %  "
              f"m_t {r['m']['m_t']:+5.1f} hh {r['m']['hh100']:+5.1f} lh {r['m']['lh100']:+5.1f} m_l {r['m']['m_l']:+5.1f}  "
              f"Γ–X {r['gx']:.3f}  Γ12 err {r['g12']:+.2f}")
    ink, ink2, grid = "#0b0b0b", "#52514e", "#e4e3df"
    col = {"13-parameter (plan)": "#8a8984", "longer two-centre shells": "#c9c8c3", "+ Sc 4s": "#2a78d6",
           "+ three-centre": "#eb6834", "+ three-centre + Sc 4s": "#1baf7a"}
    mk = {"balanced": "o", "near-gap": "s", "mass ×3": "D", "mass ×5": "^", "mass ×10": "v"}
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": ink2, "axes.labelcolor": ink, "xtick.color": ink2,
                         "ytick.color": ink2, "axes.linewidth": 0.8})
    fig, ax = plt.subplots(figsize=(7.2, 5.4), constrained_layout=True)
    ax.axvspan(0, 30, ymin=0, ymax=10 / 48, color="#1baf7a", alpha=0.10, lw=0, zorder=0)
    ax.text(1.5, 9.0, "gate (a) passes", color=ink2, fontsize=8, va="top")
    ax.axvline(30, color=ink2, lw=0.8, ls=(0, (3, 3)), zorder=1)
    ax.axhline(10, color=ink2, lw=0.8, ls=(0, (3, 3)), zorder=1)
    for r in rows:
        ax.scatter(r["rms"], r["worst"], s=60, marker=mk[r["wt"]], color=col[r["fam"]], edgecolor="white", lw=1.5, zorder=3)
        dx, dy, ha = OFFSET.get(r["v"], (1.0, 0.6, "left"))
        ax.annotate(r["v"], (r["rms"], r["worst"]), xytext=(r["rms"] + dx, r["worst"] + dy), fontsize=7.5, color=ink2,
                    ha=ha, va="center")
    for fam, c in col.items():
        ax.scatter([], [], s=50, marker="o", color=c, label=fam)
    for wt, m in mk.items():
        ax.scatter([], [], s=40, marker=m, facecolor="none", edgecolor=ink2, label=f"weighting: {wt}")
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 48)
    ax.set_xlabel("RMS within ±1 eV of the band edges (meV; worse of path and grid)")
    ax.set_ylabel("worst transport-mass error (%)")
    ax.set_title("Bulk HSE06: every SK variant against gate (a)", loc="left", color=ink, fontsize=10)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(color=grid, lw=0.6, zorder=0)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=7.5, labelcolor=ink)
    fig.savefig(FIG, dpi=200, bbox_inches="tight", facecolor="white")
    print("wrote", FIG)


if __name__ == "__main__":
    main()
