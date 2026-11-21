"""T0: exercise the SK builder and fitter on the literature HSE06 bulk targets (Deng 2015, Mu 2021).

Debugging run for the scaffold (plan §4 T0), not a production fit: the production targets are our own NC-HSE06
(A1/A2). Every variant has the band-character targets (X3 = d_xy, Γ25′ = t2g, strict), the window floors (no
extra model band inside the target window) and free dd signs (pd and pp keep Harrison signs; see sk.CANONICAL_SIGN):
  A  Deng path, uniform σ (50 meV)
  C  Deng path, near-gap weighted (σ 20 meV within 1.5 eV below the VBM / 2.5 eV above it, 150 meV elsewhere)
  D  as C, with |dd2_*| ≤ 0.1 eV (the A_skmodel "physical 3NN" question)
  E  as C, 2NN only (dd2 = 0): the m_l lock
  Mu Mu high-symmetry table + Mu masses
Energies are VBM-referenced, so Ep is fixed at 0 (energy gauge).
History (2026-11-11): without character targets the fit landed with the CBM at X as d_z² (swapped characters);
without floors it put a spurious second CB at X 0.5–0.9 eV above the CBM. Both are kept out by construction now.

Run from the project root:  PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/t0_sk_literature_fit.py
Writes fits/T0_literature/sk_lit_fits.json and figures/T0_SK_literature_fit.png.
"""
import json
import time
from pathlib import Path

import numpy as np

from scnfet import sk, skfit, targets_lit as tl

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "fits" / "T0_literature"
FIG = ROOT / "figures" / "T0_SK_literature_fit.png"
N_STARTS = 24

# toy fit of A_skmodel.md (Γ25′ variant), shifted so that Ep = 0, plus small 3NN
X0 = dict(Ep=0.0, Et2g=4.73, Eeg=6.95, pd_s=-2.40, pd_p=1.38, dd1_s=-0.12, dd1_p=0.58, dd1_d=-0.22,
          pp_s=0.80, pp_p=-0.18, dd2_s=-0.02, dd2_p=0.02, dd2_d=0.02)
DD2 = ("dd2_s", "dd2_p", "dd2_d")


def near_gap(T):
    return ~T.floor & (((T.band < 3) & (T.E > -1.5)) | ((T.band >= 3) & (T.E < 2.5)))


def bounds(names, fixed, dd2_cap=None):
    free = [n for n in names if n not in fixed]
    lo, hi = skfit.sign_bounds(names, free)
    if dd2_cap is not None:
        for i, n in enumerate(free):
            if n in DD2:
                lo[i], hi[i] = -dd2_cap, dd2_cap
    return lo, hi


def key_quantities(m, x, a):
    P = sk.fcc_points(a)
    eG = m.eigvals(x, P["G"][None])[0]
    eX = m.eigvals(x, P["X"][None])[0]
    vbm = eG[2]
    q = dict(gap_GX=eX[3] - vbm, gap_XX=eX[3] - eX[2], gap_GG=eG[3] - vbm, G12=eG[6] - vbm,
             m_l=m.mass(x, P["X"], [0, 0, 1], 3), m_t=m.mass(x, P["X"], [1, 0, 0], 3),
             X_cb2=eX[4] - vbm)                                    # second CB at X (Deng 6.14, Mu 5.95)
    for lab, u in (("100", [1, 0, 0]), ("110", [1, 1, 0])):
        for b in (2, 1, 0):
            q[f"m_h{b}_{lab}"] = -m.mass(x, P["G"], u, b)       # hole masses, positive
    return {k: float(v) for k, v in q.items()}


def run_variant(name, m, T, fixed, dd2_cap, a, deng_eval):
    x0 = m.x_from(X0)
    if "dd2_s" in fixed or dd2_cap is not None:
        x0[[m.names.index(n) for n in DD2]] = 0.0
    t = time.time()
    r = skfit.fit(m, T, x0, fixed=fixed, bounds=bounds(m.names, fixed, dd2_cap), n_starts=N_STARTS, seed=0)
    rec = dict(variant=name, targets=T.label, fixed=list(fixed), dd2_cap=dd2_cap, seconds=round(time.time() - t, 1),
               cost=float(r.cost), n_res=int(r.n_res), params={k: float(v) for k, v in r.params().items()},
               sigma={k: (None if not np.isfinite(v) else float(v)) for k, v in r.sigma().items()},
               at_bound=r.at_bound, start_costs=[float(c) for c, _ in r.starts],
               near_best_spread={k: [float(lo), float(hi)] for k, (lo, hi) in r.near_best_spread.items()},
               key=key_quantities(m, r.x, a), report=skfit.report(m, r.x, T))
    if deng_eval is not None:                       # every variant scored on the same Deng path, uniform
        vals, _ = skfit.model_values(m, r.x, deng_eval)
        eq = ~deng_eval.floor
        err = (vals - deng_eval.E)[eq]
        ng = near_gap(deng_eval)[eq]
        rec["deng_floor_violations"] = int(np.sum(deng_eval.floor & (vals < deng_eval.E - 1e-6)))
        rec["deng_path_rms_meV"] = dict(all=1e3 * float(np.sqrt(np.mean(err ** 2))),
                                         near_gap=1e3 * float(np.sqrt(np.mean(err[ng] ** 2))),
                                         near_gap_max=1e3 * float(np.max(np.abs(err[ng]))))
    print(f"{name:3s} cost {r.cost:9.1f}  {rec['seconds']:5.1f} s  at_bound {r.at_bound}  "
          + "  ".join(f"{k} {v:.3f}" for k, v in rec["key"].items() if k in ("gap_GX", "gap_GG", "m_l", "m_t", "X_cb2"))
          + f"  X3 dxy {rec['report']['chars']['X3 = d_xy'][0]:.2f}"
          + (f"  Deng RMS all/near {rec['deng_path_rms_meV']['all']:.0f}/{rec['deng_path_rms_meV']['near_gap']:.0f} meV"
             if deng_eval is not None else ""), flush=True)
    return rec, r.x


def plot(models, results, xs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, ink2, grid, ref = "#0b0b0b", "#52514e", "#e4e3df", "#9b9a96"
    col = {"A": "#2a78d6", "C": "#eb6834", "Mu": "#1baf7a"}
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": ink2, "axes.labelcolor": ink, "xtick.color": ink2,
                         "ytick.color": ink2, "axes.linewidth": 0.8})
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True, constrained_layout=True)

    def draw(ax, a, label_sets, title):
        P = sk.fcc_points(a)
        kp, s, ticks = sk.path(P, ["L", "G", "X", "U", None, "K", "G"], n=120)
        for name in label_sets:
            E = models[name].eigvals(xs[name], kp)
            E = E - models[name].eigvals(xs[name], P["G"][None])[0, 2]
            for b in range(E.shape[1]):
                ax.plot(s, E[:, b], color=col[name], lw=2.0, zorder=3, solid_capstyle="round",
                        label=LABEL[name] if b == 0 else None)
        tick_s, tick_l = [], []
        for t, l in ticks:
            l = "Γ" if l == "G" else l
            if tick_s and abs(t - tick_s[-1]) < 1e-9:
                tick_l[-1] = f"{tick_l[-1]}|{l}"
            else:
                tick_s.append(t)
                tick_l.append(l)
        for t in tick_s:
            ax.axvline(t, color=grid, lw=0.8, zorder=0)
        ax.axhline(0, color=grid, lw=0.8, zorder=0)
        ax.set_xticks(tick_s, tick_l)
        ax.set_xlim(tick_s[0], tick_s[-1])
        ax.set_ylim(-6.5, 7.0)
        ax.set_title(title, loc="left", color=ink, fontsize=10)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        return kp, s, tick_s

    # (a) Deng: digitised curves in grey on the same path coordinate
    ax = axes[0]
    kp, s, tick_s = draw(ax, 4.50, ["A", "C"], "(a) Deng 2015 HSE06, a = 4.50 Å")
    xpt = {"L": 95.8, "G": 141.5, "X": 194.3, "UK": 212.9, "G2": 268.9}
    seg_pt = [(xpt["L"], xpt["G"]), (xpt["G"], xpt["X"]), (xpt["X"], xpt["UK"]), (xpt["UK"], xpt["G2"])]
    first = True
    for px, pe in tl.deng_paths():
        sx = np.interp(px, [p for seg in seg_pt for p in seg],
                       [tick_s[0], tick_s[1], tick_s[1], tick_s[2], tick_s[2], tick_s[3], tick_s[3], tick_s[4]])
        ax.plot(sx, pe, color=ref, lw=4.0, alpha=0.55, zorder=2, solid_capstyle="round",
                label="Deng HSE06 (digitised)" if first else None)
        first = False
    ax.axhline(5.5, color=ink2, lw=0.8, ls=(0, (3, 3)), zorder=1, label="window top (5.5 eV)")
    ax.set_ylabel("E − E_VBM (eV)")
    ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.62), frameon=False, fontsize=8, labelcolor=ink)

    # (b) Mu: table points
    ax = axes[1]
    draw(ax, 4.48, ["Mu"], "(b) Mu 2021 HSE06 (semicore frozen), a = 4.48 Å")
    P = sk.fcc_points(4.48)
    kp2, s2, ticks2 = sk.path(P, ["L", "G", "X", "U", None, "K", "G"], n=120)
    spos = {"L": ticks2[0][0], "G": ticks2[1][0], "X": ticks2[2][0], "K": ticks2[3][0]}
    first = True
    for lab, pts in tl.MU_TABLE.items():
        for b, e in pts:
            ax.plot(spos[lab], e, "o", ms=6, mfc="none", mec=ink, mew=1.2, zorder=4,
                    label="Mu HSE06 (high-symmetry points)" if first else None)
            first = False
    ax.legend(loc="lower right", frameon=False, fontsize=8, labelcolor=ink)

    txt = []
    for name in ("A", "C", "Mu"):
        k = results[name]["key"]
        txt.append(f"{LABEL[name].split(':')[0]}: Γ–X {k['gap_GX']:.2f}, Γ–Γ {k['gap_GG']:.2f}, X 2nd CB {k['X_cb2']:.2f} eV; "
                   f"m_l {k['m_l']:.2f}, m_t {k['m_t']:.2f}")
    fig.text(0.01, -0.02, "   ·   ".join(txt) + "\nTargets: Deng Γ–X 0.94 (dig.), Γ–Γ 3.73, X 2nd CB 6.14, m_l 1.44, m_t 0.24;  "
             "Mu Γ–X 0.79, Γ–Γ 3.58, X 2nd CB 5.95, m_l 1.43, m_t 0.17.  Debug fits of the T0 scaffold, not production.",
             fontsize=7.5, color=ink2, va="top")
    FIG.parent.mkdir(exist_ok=True)
    fig.savefig(FIG, dpi=200, bbox_inches="tight", facecolor="white")
    print("wrote", FIG)


LABEL = {"A": "A: all-window", "C": "C: near-gap weighted", "Mu": "Mu: high-symmetry fit"}


def main():
    import sys
    if "--plot-only" in sys.argv:
        res = json.loads((OUT / "sk_lit_fits.json").read_text())["results"]
        models = {"A": sk.build(sk.rocksalt(4.50)), "C": sk.build(sk.rocksalt(4.50)), "Mu": sk.build(sk.rocksalt(4.48))}
        plot(models, res, {n: models[n].x_from(res[n]["params"]) for n in models})
        return
    m45 = sk.build(sk.rocksalt(4.50))
    m448 = sk.build(sk.rocksalt(4.48))
    deng_eval = tl.deng_path_targets()
    results, xs = {}, {}

    def deng_T(weighted):
        T = tl.deng_path_targets()
        if weighted:
            T.sigma = np.where(T.floor, 0.05, np.where(near_gap(T), 0.02, 0.15))
        return T

    for name, weighted, fixed, cap in (("A", False, ("Ep",), None), ("C", True, ("Ep",), None),
                                       ("D", True, ("Ep",), 0.1), ("E", True, ("Ep",) + DD2, None)):
        results[name], xs[name] = run_variant(name, m45, deng_T(weighted), fixed, cap, 4.50, deng_eval)
    results["Mu"], xs["Mu"] = run_variant("Mu", m448, tl.mu_targets(), ("Ep",), None, 4.48, None)

    OUT.mkdir(parents=True, exist_ok=True)
    meta = dict(note=__doc__, n_starts=N_STARTS, x0=X0, model="scnfet.sk rocksalt_shells(a0), Harrison η "
                + json.dumps(sk.DEFAULT_ETA), deng_targets=dict(n_k=int(len(deng_eval.k)), n=int(len(deng_eval.E))))
    (OUT / "sk_lit_fits.json").write_text(json.dumps(dict(meta=meta, results=results), indent=1))
    print("wrote", OUT / "sk_lit_fits.json")
    plot({"A": m45, "C": m45, "Mu": m448}, results, xs)


if __name__ == "__main__":
    main()
