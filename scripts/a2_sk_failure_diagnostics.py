"""Why the pre-registered joint fit failed gate (a): post-hoc diagnostics only (no gate is re-judged here).

  1. Transfer, bond by bond: each NN Sc–N, 2NN N–N and 2NN Sc–Sc block of the Wannier H(R) is projected onto the
     two-centre SK form (least squares in the bond frame; the residual is the non-two-centre part). Bulk (A1) and
     sheet (A2) values are compared with the Harrison scaling the model uses (default η, not yet fitted).
  2. A sheet-only fit (same recipe and targets, no bulk): can this SK form fit the sheet at all?
  3. Orbital character of the sheet states at Γ: Wannier vs the joint fit.

Run from the project root: OMP_NUM_THREADS=1 PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/a2_sk_failure_diagnostics.py
Writes fits/A2_joint/diagnostics.json.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from scnfet import sk, skfit

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("J", ROOT / "scripts" / "a2_sk_fit_bulk_sheet.py")
J = importlib.util.module_from_spec(spec)
spec.loader.exec_module(J)
OUT = ROOT / "fits" / "A2_joint"
L = {"Sc": 2, "N": 1}
NPAR = {(2, 1): 2, (1, 2): 2, (1, 1): 2, (2, 2): 3}


def bonds(hm, st, orb_site, pair, dwin):
    """All (i, j, R, vector) with species pair and |d| in dwin, i < j or i == j with R ≠ 0, one direction only."""
    out = []
    starts = {s: int(np.where(orb_site == s)[0][0]) for s in range(len(st.species))}
    for r, R in enumerate(hm.R):
        for i, si in enumerate(st.species):
            for j, sj in enumerate(st.species):
                if (si, sj) != pair:
                    continue
                v = R @ st.cell + st.pos[j] - st.pos[i]
                d = np.linalg.norm(v)
                if dwin[0] <= d <= dwin[1]:
                    ni, nj = len(sk.ORBITALS[si]), len(sk.ORBITALS[sj])
                    B = hm.H[r, starts[i]:starts[i] + ni, starts[j]:starts[j] + nj].real
                    out.append((i, j, r, v, B))
    return out


def project(B, li, lj, v):
    """Least-squares two-centre (σ, π[, δ]) for one block; returns V and the residual fraction ‖B − B_SK‖/‖B‖."""
    c = v / np.linalg.norm(v)
    n = NPAR[(li, lj)]
    basis = [sk.sk_block(li, lj, c, tuple(np.eye(n)[k])).ravel() for k in range(n)]
    A = np.array(basis).T
    V, *_ = np.linalg.lstsq(A, B.ravel(), rcond=None)
    res = np.linalg.norm(B.ravel() - A @ V) / max(np.linalg.norm(B), 1e-12)
    return V, res


def shell_table(hm, st, orb_site, label):
    rows = {}
    for name, pair, win, li, lj in (("pd", ("Sc", "N"), (1.9, 2.5), 2, 1), ("pp", ("N", "N"), (2.9, 3.4), 1, 1),
                                    ("dd1", ("Sc", "Sc"), (2.9, 3.4), 2, 2)):
        for i, j, r, v, B in bonds(hm, st, orb_site, pair, win):
            d = np.linalg.norm(v)
            kind = "in-plane" if abs(v[2]) < 0.3 else ("along z" if np.hypot(v[0], v[1]) < 0.3 else "inter-plane")
            if label == "bulk":
                kind = "bulk"
            V, res = project(B, li, lj, v)
            rows.setdefault((name, kind), []).append((d, V, res))
    out = {}
    for (name, kind), rr in rows.items():
        d = np.array([x[0] for x in rr])
        V = np.array([x[1] for x in rr])
        res = np.array([x[2] for x in rr])
        out[f"{name} {kind}"] = dict(n=len(rr), d=float(d.mean()), V=V.mean(0).round(4).tolist(),
                                     V_spread=(V.max(0) - V.min(0)).round(4).tolist(), resid=float(res.mean()))
    return out


def character(model_or_none, H_k, lab, nstates):
    e, U = np.linalg.eigh(H_k)
    groups = {"Sc t2g-like (x2-y2, xz, yz)": ["dx2-y2", "dyz", "dzx"], "Sc eg in-plane (xy)": ["dxy"],
              "Sc z2": ["dz2"], "N pxy": ["px", "py"], "N pz": ["pz"]}
    rows = []
    for n in range(nstates):
        w = np.abs(U[:, n]) ** 2
        rows.append((float(e[n]), {g: round(float(sum(w[i] for i, l in enumerate(lab) if l.split(":")[1] in o)), 3)
                                   for g, o in groups.items()}))
    return rows


def main():
    lab = sys.argv[1] if len(sys.argv) > 1 else "HSE06"
    mb, ms = J.T3.make_model(False), J.sheet_model()
    names = ms.names
    cb, cs = J.F.Case(lab, mb), J.SheetCase(lab, ms)
    rec = json.loads((OUT / f"sk_joint_{lab}.json").read_text())
    x = np.array([rec["params"][n] for n in names])
    out = {"label": lab}

    # 1. bond-resolved transfer
    tb = shell_table(cb.hm, mb.structure, mb.orb_site, "bulk")
    ts = shell_table(cs.hm, ms.structure, ms.orb_site, "sheet")
    eta = sk.DEFAULT_ETA
    print(f"== 1. Wannier hoppings projected on two-centre SK ({lab}); Harrison η = {eta} ==")
    print(f"{'bond':22s} {'d (Å)':>6s}  {'Wannier V (σ, π[, δ]) eV':34s} {'non-2c':>6s} | {'bulk × (d_b/d)^η':28s} | ratio σ")
    trans = {}
    for k, v in {**tb, **ts}.items():
        name = k.split()[0]
        b = tb[f"{name} bulk"]
        pred = (np.array(b["V"]) * (b["d"] / v["d"]) ** eta[name]).round(3)
        ratio = v["V"][0] / pred[0] if abs(pred[0]) > 1e-6 else float("nan")
        trans[k] = dict(**v, harrison_pred=pred.tolist(), ratio_sigma=round(float(ratio), 3))
        print(f"{k:22s} {v['d']:6.3f}  {str([round(t, 3) for t in v['V']]):34s} {v['resid']:6.2f} | "
              f"{str(pred.tolist()):28s} | {ratio:5.2f}")
    out["transfer"] = trans

    # 2. sheet-only fit (diagnostic)
    Ts = cs.targets()
    base = json.loads((ROOT / "fits" / "A1_bulk" / "sk_fit_t3.json").read_text())["Tn" if lab == "HSE06" else "PT"]["params"]
    proj = json.loads((ROOT / "fits" / "A1_bulk" / "sk_projection.json").read_text())[lab]["projected_params"]
    shifts = J.seed_shifts(cs, proj)
    free = [n for n in names if n not in J.FIXED]
    lo, hi = skfit.sign_bounds(names, free)
    x0s = [np.array([{**base, **shifts}.get(n, 0.0) for n in names]), x.copy()]
    best = None
    for x0 in x0s:
        for n in J.FIXED:
            x0[names.index(n)] = 0.0
        r = skfit.fit_joint([(ms, Ts)], names, x0, fixed=J.FIXED, bounds=(lo, hi), n_starts=J.N_STARTS, seed=2)
        if best is None or r.cost < best.cost:
            best = r
    ss = cs.score(ms, best.x)
    out["sheet_only"] = dict(cost=float(best.cost), rms=[ss["path"]["rms_pm1eV_meV"], ss["grid"]["rms_pm1eV_meV"]],
                             masses={k: round(100 * v["rel"], 1) for k, v in ss["masses"].items()},
                             levels={k: [round(v["sk"], 3), round(v["ref"], 3)] for k, v in ss["levels"].items()},
                             params={k: round(float(v), 3) for k, v in best.params().items()}, gate_a=ss["gate_a"])
    xb_sheet = best.x[[names.index(n) for n in mb.names]]
    sb = cb.score(mb, xb_sheet)
    out["sheet_only"]["bulk_with_these_params"] = dict(
        rms=[sb["path"]["rms_pm1eV_meV"], sb["grid"]["rms_pm1eV_meV"]],
        masses={k: round(100 * v["rel"], 1) for k, v in sb["masses"].items()})
    print("\n== 2. Sheet-only fit (diagnostic, same recipe, no bulk) ==")
    print(json.dumps({k: v for k, v in out["sheet_only"].items() if k != "params"}, indent=1))
    print("params", out["sheet_only"]["params"])

    # 3. character at Γ
    G = np.zeros((1, 3))
    Hw = J.w90.hk_frac(cs.hm, ms.structure.kfrac(G))[0]
    Hj = ms.hk(x, G)[0]
    Hs = ms.hk(best.x, G)[0]
    out["gamma_character"] = {}
    print("\n== 3. Sheet states at Γ (E − VBM, eV; orbital weights) ==")
    for nm, H in (("Wannier", Hw), ("joint fit", Hj), ("sheet-only fit", Hs)):
        rows = character(None, H, ms.orb_label, 10)
        vbm = rows[J.NV_SHEET - 1][0]
        out["gamma_character"][nm] = [(round(e - vbm, 3), w) for e, w in rows]
        print(nm)
        for e, w in rows:
            print(f"  {e - vbm:+7.3f}  " + "  ".join(f"{g.split(' (')[0]} {v:.2f}" for g, v in w.items() if v >= 0.05))
    (OUT / f"diagnostics_{lab}.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
