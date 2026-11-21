"""A1 → SK, step 2: fit the bulk SK model to our own bands (plan D3), HSE06 (production) and PBE (the twin, D7).

Bulk only, so preliminary: the production fit is bulk + sheet together once A2 exists.

Targets (VBM-referenced; the model's VBM is pinned by `ref`, so Ep = 0 fixes the energy gauge):
  * path L–Γ–X–W–K–Γ: Wannier-interpolated bands (`hse_k16` / `pbe_k16`) up to VBM + 4 eV. Below that the
    16³ interpolation is good to ≤ 6.8 meV (A1 round 2);
  * grid: the DFT eigenvalues of the 8 model bands (QE bands 6–13) at the HSE 16³ IBZ points (145 k) / the PBE
    16³ grid reduced by O_h (sorted |k|), at every energy up to the frozen-window top (6.5 / 5.4 eV). These are
    exact and off-path. Above the window, the band gets a floor (model ≥ window top): no spurious band inside;
  * σ: 20 meV near the gap (VB ≥ −1.5 eV, CB ≤ CBM + 2.5 eV), 100 meV elsewhere;
  * masses from the Wannier model (band-edge curvature): m_l, m_t at X (rel σ 0.03); heavy/light hole along [100]
    and [110] at Γ (rel σ 0.05);
  * characters: X3 = d_xy, Γ25′ = t2g (strict), as in T0.
Variants: H (HSE, 13 parameters), H2 (HSE, 2NN only: the m_l lock on our data), P (PBE twin, 13 parameters).
Two seeds per variant, each with multi-start: the T0 Deng fit C and the Wannier projection (scripts/a1_sk_projection.py).
Scored on the gate (a) metric (plan §4): RMS ≤ 30 meV within ±1 eV of the band edges; transport masses ±10 %.

Run from the project root: PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/a1_sk_fit.py [--plot-only]
Writes fits/A1_bulk/sk_fits.json and figures/A1_SK_bulk_fit.png.
"""
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

from scnfet import sk, skfit, w90, targets_lit as tl

ROOT = Path(__file__).resolve().parents[1]
A1 = ROOT / "dft" / "A1"
OUT = ROOT / "fits" / "A1_bulk"
FIG = ROOT / "figures" / "A1_SK_bulk_fit.png"
A = 4.501
N_STARTS = 16
E_SPLIT = 4.0
SIG_NEAR, SIG_FAR = 0.02, 0.10
CASES = {"HSE06": dict(seed="hse_k16", grid="hse_ibz_k16.out", froz=6.5, outer=13.5),
         "PBE": dict(seed="pbe_k16", grid="pbe_nscf_k16.out", froz=5.4, outer=11.8)}
# outer: below the first band outside the 9-band Sc d + N p + Sc 4s manifold (QE band 15: ≥ 13.9 / 12.2 eV), the A1
# outer windows. Used as the grid-target cap when the model has the Sc 4s pseudo-site (nb = 9).
DD2 = ("dd2_s", "dd2_p", "dd2_d")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


PROJ = _load("a1_sk_projection", ROOT / "scripts" / "a1_sk_projection.py")
PARSE = _load("parse_a1_w90", A1 / "parse_a1_w90.py")


def grid_eigs(fname, nb=8):
    k, e = PARSE.qe_bands(A1 / "outputs" / fname)
    u = 2 * np.pi / A
    key = np.round(np.sort(np.abs(k / u), axis=1), 4)
    _, first = np.unique(key, axis=0, return_index=True)
    first = np.sort(first)
    return k[first], e[first][:, 5:5 + nb]


class Case:
    def __init__(self, label, model):
        c = CASES[label]
        self.label, self.froz = label, c["froz"]
        self.nb = model.norb                       # 8 (Sc d + N p) or 9 (+ Sc 4s)
        self.cap = self.froz if self.nb == 8 else c["outer"]
        _, _, _, self.hm, _ = PROJ.to_model(c["seed"], model)
        self.st = model.structure
        self.vbm = self.E_w90(np.zeros((1, 3)))[0, 2]
        kg, eg = grid_eigs(c["grid"], self.nb)
        self.kg, self.eg = kg, eg - eg[:, 2].max()
        P = sk.fcc_points(A)
        self.masses = [("m_l", P["X"], [0, 0, 1], 3, 0.03), ("m_t", P["X"], [1, 0, 0], 3, 0.03),
                       ("hh100", P["G"], [1, 0, 0], 2, 0.05), ("lh100", P["G"], [1, 0, 0], 0, 0.05),
                       ("hh110", P["G"], [1, 1, 0], 2, 0.05), ("lh110", P["G"], [1, 1, 0], 0, 0.05)]
        self.mass_ref = {m[0]: self.mass_w90(*m[1:4]) for m in self.masses}

    def E_w90(self, k):
        return np.linalg.eigvalsh(w90.hk_frac(self.hm, self.st.kfrac(k)))

    def mass_w90(self, k0, u, band, h=0.004):
        u = np.asarray(u, float) / np.linalg.norm(u)
        E = self.E_w90(np.array([k0 - h * u, k0, k0 + h * u]))[:, band]
        return sk.HB2M / ((E[0] + E[2] - 2 * E[1]) / h ** 2)

    def path(self, n):
        return sk.path(sk.fcc_points(A), ["L", "G", "X", "W", "K", "G"], n=n)

    def targets(self):
        kp, _, _ = self.path(30)
        ep = self.E_w90(kp) - self.vbm
        cbm = min(ep[:, 3].min(), self.eg[:, 3].min())
        ks, ik, band, E, sig, grp, flo = [], [], [], [], [], [], []

        def add(kk, ee, src):
            for i in range(len(kk)):
                n0 = len(ks)
                ks.append(kk[i])
                for b in range(ee.shape[1]):
                    e = ee[i, b]
                    if src == "path" and e > E_SPLIT:
                        continue
                    near = (b < 3 and e >= -1.5) or (b >= 3 and e <= cbm + 2.5)
                    ik.append(n0)
                    band.append(b)
                    if src == "grid" and e > self.cap:
                        E.append(self.cap)
                        sig.append(0.05)
                        grp.append("floor")
                        flo.append(True)
                    else:
                        E.append(e)
                        sig.append(SIG_NEAR if near else SIG_FAR)
                        grp.append(("VB" if b < 3 else "CB") + ("_near" if near else "_far"))
                        flo.append(False)
        add(kp, ep, "path")
        add(self.kg, self.eg, "grid")
        P = sk.fcc_points(A)
        masses = [skfit.MassTarget(n, np.asarray(k0), np.asarray(u, float), b, self.mass_ref[n], rs)
                  for n, k0, u, b, rs in self.masses]
        return skfit.Targets(f"{self.label} bulk (W90 path ≤ {E_SPLIT} eV + DFT grid ≤ {self.cap} eV, {self.nb} bands)", np.array(ks), np.array(ik),
                             np.array(band), np.array(E), np.array(sig), masses=masses, ref=(P["G"], 2),
                             group=np.array(grp), floor=np.array(flo), chars=tl.character_targets(A))

    def score(self, model, x):
        """Gate (a): RMS within ±1 eV of the band edges, on the dense path and the grid; masses; key gaps."""
        kp, _, _ = self.path(80)
        out = {}
        for nm, kk, ref in (("path", kp, self.E_w90(kp) - self.vbm), ("grid", self.kg, self.eg)):
            e = (model.eigvals(x, kk) - model.eigvals(x, np.zeros((1, 3)))[0, 2])[:, :ref.shape[1]]
            cbm = ref[:, 3].min()
            m = np.zeros_like(ref, bool)
            m[:, :3] = ref[:, :3] >= -1.0
            m[:, 3:] = ref[:, 3:] <= cbm + 1.0
            err = (e - ref)[m]
            inwin = ref <= (E_SPLIT if nm == "path" else self.cap)
            out[nm] = dict(rms_pm1eV_meV=1e3 * float(np.sqrt(np.mean(err ** 2))),
                           max_pm1eV_meV=1e3 * float(np.max(np.abs(err))),
                           rms_window_meV=1e3 * float(np.sqrt(np.mean((e - ref)[inwin] ** 2))))
        P = sk.fcc_points(A)
        out["masses"] = {}
        for n, k0, u, b, _ in self.masses:
            mm = float(model.mass(x, k0, u, b))
            ref = self.mass_ref[n]
            out["masses"][n] = dict(sk=mm, ref=float(ref), rel=float(mm / ref - 1))
        eG, eX = (model.eigvals(x, P[q][None])[0] for q in ("G", "X"))
        rG, rX = (self.E_w90(P[q][None])[0] - self.vbm for q in ("G", "X"))
        out["gaps"] = {q: dict(sk=float(s), ref=float(r)) for q, s, r in (
            ("GX", eX[3] - eG[2], rX[3]), ("XX", eX[3] - eX[2], rX[3] - rX[2]), ("GG", eG[3] - eG[2], rG[3]),
            ("G12", eG[6] - eG[2], rG[6]))}
        out["gate_a"] = dict(rms_ok=bool(max(out["path"]["rms_pm1eV_meV"], out["grid"]["rms_pm1eV_meV"]) <= 30),
                             masses_ok=bool(all(abs(v["rel"]) <= 0.10 for v in out["masses"].values())))
        return out


def seeds(model):
    t0 = json.loads((ROOT / "fits" / "T0_literature" / "sk_lit_fits.json").read_text())["results"]["C"]["params"]
    pj = json.loads((OUT / "sk_projection.json").read_text())
    out = {"T0-C": model.x_from(t0)}
    for lab in ("HSE06", "PBE"):
        p = dict(pj[lab]["projected_params"])
        ep = p["Ep"]
        for k in ("Ep", "Et2g", "Eeg"):
            p[k] -= ep                                    # Ep = 0 gauge
        out[f"proj-{lab}"] = model.x_from(p)
    return out


def run(model, case, T, fixed, x0s, name):
    best, rec_starts = None, {}
    lo, hi = skfit.sign_bounds(model.names, [n for n in model.names if n not in fixed])
    t = time.time()
    for sname, x0 in x0s.items():
        x0 = x0.copy()
        for n in fixed:
            if n in DD2:
                x0[model.names.index(n)] = 0.0
        r = skfit.fit(model, T, x0, fixed=fixed, bounds=(lo, hi), n_starts=N_STARTS, seed=1)
        rec_starts[sname] = float(r.cost)
        if best is None or r.cost < best.cost:
            best, bseed = r, sname
    sc = case.score(model, best.x)
    rep = skfit.report(model, best.x, T)
    rec = dict(variant=name, case=case.label, fixed=list(fixed), seconds=round(time.time() - t, 1),
               cost=float(best.cost), n_res=int(best.n_res), seed_costs=rec_starts, best_seed=bseed,
               params={k: float(v) for k, v in best.params().items()},
               sigma={k: (None if not np.isfinite(v) else float(v)) for k, v in best.sigma().items()},
               at_bound=best.at_bound,
               near_best_spread={k: [float(a), float(b)] for k, (a, b) in best.near_best_spread.items()},
               score=sc, fit_report=dict(rms_meV=rep["rms_meV"], max_meV=rep["max_meV"],
                                         floor_violations=rep["floor_violations"], chars=rep["chars"]))
    m = sc["masses"]
    print(f"{name:3s} {case.label:5s} cost {best.cost:9.1f} ({bseed}) {rec['seconds']:6.1f} s | ±1 eV RMS path/grid "
          f"{sc['path']['rms_pm1eV_meV']:.0f}/{sc['grid']['rms_pm1eV_meV']:.0f} meV (max {sc['path']['max_pm1eV_meV']:.0f}) | "
          f"Γ–X {sc['gaps']['GX']['sk']:.3f}/{sc['gaps']['GX']['ref']:.3f} | "
          + " ".join(f"{n} {v['sk']:.3f}/{v['ref']:.3f}" for n, v in m.items())
          + f" | gate a {sc['gate_a']} | at_bound {best.at_bound}", flush=True)
    return rec, best.x


def plot(model, cases, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, ink2, grid, ref = "#0b0b0b", "#52514e", "#e4e3df", "#9b9a96"
    col = {"H": "#eb6834", "P": "#2a78d6", "H2": "#1baf7a", "Hn": "#1baf7a"}
    name = {"H": "H balanced", "Hn": "Hn near-gap", "H2": "H2 2NN only", "P": "P balanced"}
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": ink2, "axes.labelcolor": ink, "xtick.color": ink2,
                         "ytick.color": ink2, "axes.linewidth": 0.8})
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.6), sharey=True, constrained_layout=True)
    ext = OUT / "sk_fit_extensions.json"
    if ext.exists() and "Hn" in json.loads(ext.read_text()):
        res = {**res, "Hn": json.loads(ext.read_text())["Hn"]}
    hse_vars = ("H", "Hn") if "Hn" in res else ("H", "H2")
    for ax, (lab, vars_, title) in zip(axes, (("HSE06", hse_vars, "(a) HSE06: SK fits against our Wannier bands"),
                                              ("PBE", ("P",), "(b) PBE twin: SK fit against the PBE Wannier bands"))):
        c = cases[lab]
        kp, s, ticks = c.path(120)
        ew = c.E_w90(kp) - c.vbm
        for b in range(8):
            ax.plot(s, ew[:, b], color=ref, lw=4.0, alpha=0.55, zorder=2, solid_capstyle="round",
                    label=f"{lab} Wannier90 (16³)" if b == 0 else None)
        for v in vars_:
            r = res[v]
            x = model.x_from(r["params"])
            e = model.eigvals(x, kp)
            e -= model.eigvals(x, np.zeros((1, 3)))[0, 2]
            sc = r["score"]
            lw, ls = (1.8, "-") if v in ("H", "P") else (1.3, (0, (4, 2)))
            mm = sc["masses"]
            for b in range(8):
                ax.plot(s, e[:, b], color=col[v], lw=lw, ls=ls, zorder=3,
                        label=(f"SK {name[v]}: ±1 eV RMS {sc['path']['rms_pm1eV_meV']:.0f} meV; "
                               f"m_t {100 * mm['m_t']['rel']:+.0f} %, hh {100 * mm['hh100']['rel']:+.0f} %, "
                               f"lh {100 * mm['lh100']['rel']:+.0f} %" if b == 0 else None))
        tick_s, tick_l = [], []
        for t, l in ticks:
            tick_s.append(t)
            tick_l.append("Γ" if l == "G" else l)
        for t in tick_s:
            ax.axvline(t, color=grid, lw=0.8, zorder=0)
        ax.axhline(0, color=grid, lw=0.8, zorder=0)
        ax.axhline(E_SPLIT, color=ink2, lw=0.8, ls=(0, (3, 3)), zorder=1,
                   label=f"path targets end (+{E_SPLIT:.0f} eV); grid targets to +{c.froz} eV")
        ax.set_xticks(tick_s, tick_l)
        ax.set_xlim(tick_s[0], tick_s[-1])
        ax.set_ylim(-6.5, 8.0)
        ax.set_title(title, loc="left", color=ink, fontsize=10)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.07), frameon=False, fontsize=7.5, labelcolor=ink)
    axes[0].set_ylabel("E − E_VBM (eV)")
    FIG.parent.mkdir(exist_ok=True)
    fig.savefig(FIG, dpi=200, bbox_inches="tight", facecolor="white")
    print("wrote", FIG)


def main():
    model = sk.build(sk.rocksalt(A), shells=sk.rocksalt_shells(A))
    cases = {lab: Case(lab, model) for lab in CASES}
    if "--plot-only" in sys.argv:
        plot(model, cases, json.loads((OUT / "sk_fits.json").read_text())["results"])
        return
    x0s = seeds(model)
    res = {}
    for name, lab, fixed in (("H", "HSE06", ("Ep",)), ("H2", "HSE06", ("Ep",) + DD2), ("P", "PBE", ("Ep",))):
        T = cases[lab].targets()
        res[name], _ = run(model, cases[lab], T, fixed, x0s, name)
        res[name]["n_targets"] = dict(n=int(len(T.E)), n_k=int(len(T.k)), floors=int(T.floor.sum()))
    meta = dict(note=__doc__, a=A, n_starts_per_seed=N_STARTS, e_split=E_SPLIT, sigma_near=SIG_NEAR, sigma_far=SIG_FAR,
                model="scnfet.sk rocksalt_shells(4.501), Harrison η " + json.dumps(sk.DEFAULT_ETA) + " (scale 1 at bulk)",
                mass_refs={lab: c.mass_ref for lab, c in cases.items()})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "sk_fits.json").write_text(json.dumps(dict(meta=meta, results=res), indent=1, default=float))
    print("wrote", OUT / "sk_fits.json")
    plot(model, cases, res)


if __name__ == "__main__":
    main()
