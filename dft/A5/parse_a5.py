#!/usr/bin/env python3
"""A5 analysis. stage2: the A5-1b runs at the relaxed hold-outs (README.md, "A5-1b staged").

- Gate M criterion 1 as registered (A5 README item 2): magnetic if ΔE_mag = E_NM − min(E_FM, E_AFM-edges) > 5 meV per
  edge atom at k 15; edge atoms as Gate M (outermost host row of each edge in the ideal cut, both planes).
- PBE gap on the 31-k Γ–X̄ path and on the k 15 grid; band-edge character (edge weight, |ψ|², orbital shares), as
  parse_a4a criterion 4 and edge_bands.py. Descriptive.
- Band survey for A5-2's NSCF band count, as A4a's for A4b-1: smallest N whose band minimum on the path is ≥ VBM + L,
  L = the PBE outer-window level of cost test T2 (parse_a4a.L_PBE), also at L + 0.5 and L + 1.0 eV.

The A4a helpers are reused unchanged; only their output directory and structure source are pointed at A5.

    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P dft/A5/parse_a5.py        (from the project root, after copying back)
Writes dft/A5/a5_stage2.json.
"""
import importlib.util
import json
import re
import sys
import types
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_argv, sys.argv = sys.argv, ["x"]
M5 = _load("make_a5", HERE / "make_a5.py")
EB = _load("edge_bands", ROOT / "dft" / "A4a" / "edge_bands.py")     # loads its own parse_a4a as EB.pa
sys.argv = _argv
PA = EB.pa
PA.HERE = HERE                                       # PA.out(tag) → dft/A5/outputs/<tag>.out
PA.MK = types.SimpleNamespace(structure=M5.structure, ALL=M5.HOLDOUTS)   # PA.edge_masks → the A5 ideal cuts
RY_EV = PA.RY_EV
THRESH = dict(mag_meV_per_edge=5.0)
COMPARE = {"r110_8": "r110_4", "r100_6": "r100_3"}   # the A4a calibration ribbon with the same edge


def character(key, nocc, kp, ep, left, right, hmask):
    """Band-by-band character at Γ and X̄, bands nocc−4 … nocc+3 (1-based numbers in the output)."""
    cw = EB.comps(f"projwfc_{key}")
    assert len(cw) == len(kp), "projwfc k count differs from the bands run"
    vbm = ep[:, nocc - 1].max()
    res = {}
    for ik, lab in ((0, "G"), (len(kp) - 1, "X")):
        rows = []
        for b in range(nocc - 4, nocc + 4):
            E, wo, wa, p2 = cw[ik][b]
            assert abs(E - ep[ik, b]) < 2e-3
            tot = sum(wo.values())
            sh = {f"{s} {l}": round(wo.get((s, l), 0.0) / tot, 3) for s, l in
                  (("Sc", "s"), ("Sc", "p"), ("Sc", "d"), ("N", "s"), ("N", "p"))}
            rows.append(dict(band=b + 1, E_rel_vbm=round(E - vbm, 4), edge=round(PA.edge_weight(wa, left, right, hmask)["edge"], 3),
                             psi2=p2, **sh))
        res[lab] = rows
    return res


def stage2():
    res = {}
    for key in M5.HOLDOUTS:
        if not PA.done(f"nm_{key}"):
            continue
        left, right, hmask = PA.edge_masks(key)
        n_edge = int((left | right).sum())
        nat = len(left)
        nocc = PA.nocc_of(f"nm_{key}")
        r = dict(n_edge_atoms=n_edge, nocc=nocc)
        runs = {t: PA.scf(t) for t in (f"nm_{key}", f"fm_{key}", f"afm_edges_{key}") if PA.done(t)}
        r["runs"] = runs
        # 1. magnetism, k 15 (registered)
        mags = {t: v for t, v in runs.items() if v.get("abs_mag") is not None}
        if len(mags) == 2:
            nm = runs[f"nm_{key}"]
            tmin = min(mags, key=lambda t: mags[t]["E_Ry"])
            dE = 1e3 * RY_EV * (nm["E_Ry"] - mags[tmin]["E_Ry"]) / n_edge
            r["c1_k15"] = dict(dE_meV_per_edge=dE, lowest=tmin, lowest_abs_mag=mags[tmin]["abs_mag"],
                               max_site_moment=max(v["max_site_moment"] for v in mags.values()),
                               all_collapsed=all(v["abs_mag"] < 0.01 for v in mags.values()),
                               per_run_meV_per_edge={t: 1e3 * RY_EV * (nm["E_Ry"] - v["E_Ry"]) / n_edge
                                                     for t, v in mags.items()})
            r["c1_magnetic"] = bool(dE > THRESH["mag_meV_per_edge"])
        # gap, survey, character
        if PA.done(f"bands_{key}"):
            kp, ep, _ = PA.P.bands(PA.out(f"bands_{key}"))
            x = kp[:, 0] / kp[-1, 0]
            vb, cb = ep[:, nocc - 1], ep[:, nocc]
            iv, ic = int(np.argmax(vb)), int(np.argmin(cb))
            ef = runs[f"nm_{key}"]["E_F_eV"]
            gap = float(cb[ic] - vb[iv])
            r["gap"] = dict(gap_path_eV=gap, k_vbm=round(float(x[iv]), 4), k_cbm=round(float(x[ic]), 4),
                            direct=bool(iv == ic), direct_gap_G_eV=float(cb[0] - vb[0]),
                            direct_gap_X_eV=float(cb[-1] - vb[-1]), gap_grid_k15_eV=PA.grid_gap(f"nm_{key}", nocc)[0],
                            E_F_eV=ef, vbm_eV=float(vb[iv]), metallic=bool(gap <= 0 or vb.max() >= ef or cb.min() <= ef),
                            nbnd=int(ep.shape[1]))
            emin = ep.min(axis=0) - vb[iv]
            n_sc = M5.structure(key).species.count("Sc")
            sv = {}
            for lab, L in (("L", PA.L_PBE), ("L+0.5", PA.L_PBE + 0.5), ("L+1.0", PA.L_PBE + 1.0)):
                ok = np.where(emin >= L)[0]
                N = int(ok[0]) + 1 if len(ok) else None
                sv[lab] = dict(level_rel_vbm=round(L, 3), nbnd=N, per_fu=round(N / n_sc, 2) if N else None)
            r["survey"] = dict(levels=sv, top_band_min_rel_vbm=float(emin[-1]), n_sc=n_sc)
            if PA.done(f"projwfc_{key}"):
                pw = PA.projwfc(f"projwfc_{key}")
                assert len(pw) == len(kp)
                assert abs(pw[iv][1][nocc - 1][0] - vb[iv]) < 2e-3 and abs(pw[ic][1][nocc][0] - cb[ic]) < 2e-3
                ew_v = PA.edge_weight(pw[iv][1][nocc - 1][1], left, right, hmask)
                ew_c = PA.edge_weight(pw[ic][1][nocc][1], left, right, hmask)
                along = {}
                for lab, b in (("VB-1", nocc - 2), ("VB", nocc - 1), ("CB", nocc), ("CB+1", nocc + 1)):
                    along[lab] = dict(G=round(PA.edge_weight(pw[0][1][b][1], left, right, hmask)["edge"], 3),
                                      X=round(PA.edge_weight(pw[-1][1][b][1], left, right, hmask)["edge"], 3),
                                      max=round(max(PA.edge_weight(s[b][1], left, right, hmask)["edge"] for _, s in pw), 3),
                                      E_G_rel_vbm=round(pw[0][1][b][0] - vb[iv], 3),
                                      E_X_rel_vbm=round(pw[-1][1][b][0] - vb[iv], 3))
                r["edges"] = dict(vbm=ew_v, cbm=ew_c, psi2_vbm=pw[iv][1][nocc - 1][2], psi2_cbm=pw[ic][1][nocc][2],
                                  uniform_share=n_edge / nat, along_path=along,
                                  character=character(key, nocc, kp, ep, left, right, hmask))
        res[key] = r
    res["thresholds"] = THRESH
    res["L_pbe_rel_vbm"] = PA.L_PBE
    (HERE / "a5_stage2.json").write_text(json.dumps(res, indent=1))
    report(res)
    return res


def report(res):
    a4a = json.loads((ROOT / "dft" / "A4a" / "a4a_stage2.json").read_text())
    for key in [k for k in M5.HOLDOUTS if k in res]:
        r = res[key]
        print(f"\n=== {key}  (edge atoms {r['n_edge_atoms']}, nocc {r['nocc']})")
        for t, v in r["runs"].items():
            print(f"   {t:20s} E {v['E_Ry']:.8f} Ry  conv {v['converged']}"
                  + (f"  |m| {v['abs_mag']:.3f} (iter1 {v['abs_mag_first_iter']:.2f}) site max {v['max_site_moment']}"
                     if "abs_mag" in v else ""))
        if "c1_k15" in r:
            c = r["c1_k15"]
            print(f"   C1 k15: ΔE_mag {c['dE_meV_per_edge']:+.4f} meV/edge atom (lowest {c['lowest']}, |m| "
                  f"{c['lowest_abs_mag']}); per run " + ", ".join(f"{t} {v:+.4f}" for t, v in c["per_run_meV_per_edge"].items())
                  + f"; all collapsed {c['all_collapsed']} → magnetic {r['c1_magnetic']}")
        if "gap" in r:
            c = r["gap"]
            ref = a4a.get(COMPARE[key], {}).get("c3", {})
            print(f"   PBE gap {c['gap_path_eV']:.4f} eV ({'direct' if c['direct'] else 'indirect'}; VBM k {c['k_vbm']}, "
                  f"CBM k {c['k_cbm']} of Γ–X̄); Γ {c['direct_gap_G_eV']:.4f}, X̄ {c['direct_gap_X_eV']:.4f}; grid "
                  f"{c['gap_grid_k15_eV']:.4f}; metallic {c['metallic']}   [{COMPARE[key]}: {ref.get('gap_path_eV', float('nan')):.4f}]")
            s = r["survey"]
            print(f"   band survey (nbnd {c['nbnd']}, top band min VBM+{s['top_band_min_rel_vbm']:.2f} eV): "
                  + "; ".join(f"{k} (VBM+{v['level_rel_vbm']:.2f}) N {v['nbnd']} = {v['per_fu']}/f.u." for k, v in s["levels"].items()))
        if "edges" in r:
            e = r["edges"]
            print(f"   edge weight VBM {e['vbm']['edge']:.3f}, CBM {e['cbm']['edge']:.3f} (uniform {e['uniform_share']:.3f}); "
                  f"|ψ|² VBM {e['psi2_vbm']:.3f}, CBM {e['psi2_cbm']:.3f}")
            for lab, v in e["along_path"].items():
                print(f"     {lab:5s} Γ {v['E_G_rel_vbm']:+.3f} eV edge {v['G']:.2f} | X̄ {v['E_X_rel_vbm']:+.3f} eV edge {v['X']:.2f} | max {v['max']:.2f}")
            for lab in ("G", "X"):
                print(f"   -- {'Γ' if lab == 'G' else 'X̄'}: band E−VBM edge |ψ|²  Sc s  Sc p  Sc d   N s   N p")
                for row in e["character"][lab]:
                    print(f"      {row['band']:4d} {row['E_rel_vbm']:+7.3f} {row['edge']:.2f} {row['psi2']:.2f}  "
                          f"{row['Sc s']:.2f}  {row['Sc p']:.2f}  {row['Sc d']:.2f}  {row['N s']:.2f}  {row['N p']:.2f}")


# ------------------------------------------------------------------------------------- A5-2 job 1: the twin windows
def _pb():
    """parse_a4b, unchanged, pointed at dft/A5/twin: its make_a4b is make_a5._make_a4b(TWIN) (the A5 relaxed geometry and
    ideal cuts) with r110_8's edge Sc s in set B, as make_a5.twin set it when it wrote the inputs."""
    _a, sys.argv = sys.argv, ["x"]
    PB = _load("parse_a4b", ROOT / "dft" / "A4b" / "parse_a4b.py")
    sys.argv = _a
    MB = M5._make_a4b(M5.TWIN)
    MB.EDGE_S = set(MB.EDGE_S) | {"r110_8"}
    PB.MB, PB.HERE = MB, M5.TWIN
    return PB


def twin_windows(keys):
    """parse_a4b.py windows on the A5 twin → twin/a4b_windows.json and twin/w90/<r>/<v>/<r>_<v>.win. parse_a4b reads
    its run facts from <dir>/a4b_inputs.json in A4b's format; that file is written here from a5_inputs.json."""
    info = json.loads((HERE / "a5_inputs.json").read_text())
    facts = {k: dict(nocc=t["nocc_nc"], nexcl=t["nexcl"], nk_w90=t["nk_w90"], nbnd_nscf=t["nscf"],
                     num_wann_B=t["num_wann_B"], num_wann_S=t["num_wann_S"])
             for k, t in ((k, info[k]["twin"]) for k in M5.HOLDOUTS)}
    (M5.TWIN / "a4b_inputs.json").write_text(json.dumps(facts, indent=1))
    PB = _pb()
    PB.windows(keys)
    wins = json.loads((M5.TWIN / "a4b_windows.json").read_text())
    for k in keys:
        if k not in wins:
            continue
        w = wins[k]
        fn = w["nonmodel"]["B"]["0.25"]
        inside = fn is not None and fn["E_eV"] <= w["cbm"] + PB.DEVICE[1]
        print(f"   {k}: amendment 1 (PBE) — first non-model state inside the device window (≤ CBM + 0.5): {inside}")


# ------------------------------------------------------------------------------------ A5-2 job 2: the twin criteria
def h4_pbe(PB, key, win, variants, nk=601):
    """H4 (amendment 2, registered for the hold-outs) on the PBE variants: every pair agrees within 5 meV over the part of
    the device window [VBM − 0.5, CBM + 0.5] both freeze, on a 601-point Γ–X̄ path. parse_a4b.h4_pairs with the PBE
    variant directories (w90/<r>/<v>/<r>_<v>_hr.dat) in place of the HSE ones."""
    from scnfet import w90 as W
    t = np.linspace(0, 0.5, nk)
    E = {v: np.linalg.eigvalsh(W.hk_1d(W.read_hr(PB.MB.W90 / key / v / f"{key}_{v}_hr.dat"), t, axis=0)) for v in variants}
    vbm, cbm = win["vbm"], win["cbm"]
    res = {}
    for i, a in enumerate(variants):
        for b in variants[i + 1:]:
            n = min(E[a].shape[1], E[b].shape[1])
            top = min(cbm + PB.DEVICE[1], win["variants"][a]["froz"], win["variants"][b]["froz"])
            m = (E[a][:, :n] >= vbm - PB.DEVICE[0]) & (E[a][:, :n] <= top)
            d = np.abs(E[a][:, :n] - E[b][:, :n])
            dmax = 1e3 * float(d[m].max()) if m.any() else float("nan")
            k_at = float(t[np.unravel_index(np.argmax(np.where(m, d, -1)), d.shape)[0]] / 0.5) if m.any() else None
            res[f"{a}-{b}"] = dict(max_meV=dmax, k_frac_at_max=k_at, top_rel_cbm=top - cbm, n=int(m.sum()),
                                   ok=bool(dmax <= PB.H4_MEV))
    return res


def twin_criteria(keys):
    """parse_a4b.py criteria, unchanged, on the A5 twin (P1–P6, decision rule W2 → W1 → back for review) → twin/a4b_criteria.json;
    then H4 between the variants run → twin/a5_h4.json. A hold-out whose H4 fails goes back for review (README A5-2)."""
    PB = _pb()
    PB.criteria(keys)
    wins = json.loads((M5.TWIN / "a4b_windows.json").read_text())
    crit = json.loads((M5.TWIN / "a4b_criteria.json").read_text())
    p = M5.TWIN / "a5_h4.json"
    old = json.loads(p.read_text()) if p.exists() else {}
    for k in keys:
        if k not in crit:
            continue
        vs = [v for v in ("W1", "W2", "W3") if v in crit[k]]
        h4 = old[k] = h4_pbe(PB, k, wins[k], vs)
        for pair, c in h4.items():
            print(f"   {k} H4 {pair}: max |ΔE| {c['max_meV']:.2f} meV (at k {c['k_frac_at_max']:.3f} of Γ–X̄) over "
                  f"[VBM − 0.5, CBM + {c['top_rel_cbm']:.2f}] ({c['n']} values) → {'ok' if c['ok'] else 'FAIL: back for review'}")
    p.write_text(json.dumps(old, indent=1, default=float))


# ---------------------------------------------------------------------------------------- A5-3: HSE gap single points
def hse_points(keys):
    """A5-3 as registered (README): the gap from the IBZ grid eigenvalues, only with QE's converged '!!' energy (the
    parse_a4b check); E_vac from the HSE potential (parse_a4b.vacuum_level); the first non-model state from projwfc on
    the HSE SCF (parse_a4b.first_nonmodel, set B with r110_8's edge Sc s) against the device window
    [VBM − 0.5, CBM + 0.5] — reported, not a stop; band-edge character (edge weight, |ψ|², orbital shares), descriptive.
    → hse/a5_hse.json."""
    PB = _pb()
    info = json.loads((HERE / "a5_inputs.json").read_text())
    twin = json.loads((M5.TWIN / "a4b_windows.json").read_text())
    H = HERE / "hse"
    p = H / "a5_hse.json"
    res = json.loads(p.read_text()) if p.exists() else {}
    for key in keys:
        f = H / "outputs" / f"hse_{key}.out"
        txt = f.read_text() if f.exists() else ""
        if "!!   total energy" not in txt or re.search(r"Maximum CPU time exceeded", txt):
            print(f"{key}: no converged HSE run ('!!' missing or a max_seconds stop); skipped")
            continue
        st = PB.MB.geometry(key)
        nocc = info[key]["twin"]["nocc_nc"]
        dexx = [float(x) for x in re.findall(r"est\. exchange err \(dexx\)\s*=\s*([-\d.E+]+)", txt)]
        kh, eh, _ = PB.P.bands(f)
        kx = kh[:, 0] * st.cell[0, 0] / (2 * np.pi)
        iv, ic = int(np.argmax(eh[:, nocc - 1])), int(np.argmin(eh[:, nocc]))
        vbm, cbm = float(eh[iv, nocc - 1]), float(eh[ic, nocc])
        vac = PB.vacuum_level(H / f"hpot_{key}.cube.gz")
        hproj = H / "outputs" / f"hproj_{key}.out"
        nm = PB.first_nonmodel(key, st, nocc, vbm, "B", hproj)
        fn = nm[f"{PB.S_OUT:.2f}"]
        # band-edge character from the same projwfc run
        left, right, hmask = PA.edge_masks(key)
        _h, PA.HERE = PA.HERE, H
        pw = PA.projwfc(f"hproj_{key}")
        PA.HERE = _h
        sts = PB.projwfc_states(None, hproj)
        assert len(pw) == len(sts) == len(kh)

        def edge(ik, b):
            E, wa, p2 = pw[ik][1][b]
            E2, w, _ = sts[ik][b]
            assert abs(E - eh[ik, b]) < 2e-3 and abs(E2 - E) < 1e-6
            tot = sum(w.values())
            sh = {}
            for (a, sp, l), v in w.items():
                sh[f"{sp} {l}"] = sh.get(f"{sp} {l}", 0.0) + v / tot
            return dict(E_rel_vbm=round(E - vbm, 4), edge=round(PA.edge_weight(wa, left, right, hmask)["edge"], 3),
                        psi2=p2, shares={k_: round(v, 3) for k_, v in sorted(sh.items(), key=lambda t: -t[1])[:4]})
        q = twin.get(key, {})
        r = res[key] = dict(
            nbnd=int(eh.shape[1]), nk_ibz=int(eh.shape[0]), exx_outer_steps=len(dexx), dexx_final_Ry=dexx[-1],
            E_Ry=float(re.findall(r"^!!\s+total energy\s+=\s+(-?[\d.]+)", txt, re.M)[-1]),
            gap_grid=cbm - vbm, k_vbm_frac=round(float(abs(kx[iv])), 4), k_cbm_frac=round(float(abs(kx[ic])), 4),
            direct=bool(iv == ic), vbm=vbm, cbm=cbm, vacuum=vac, E_vac_rel_vbm=vac["E_vac_eV"] - vbm, nonmodel=nm,
            device_window=dict(top_rel_cbm=PB.DEVICE[1],
                               first_nonmodel_rel_cbm=None if fn is None else fn["E_eV"] - cbm,
                               inside=bool(fn is not None and fn["E_eV"] <= cbm + PB.DEVICE[1]),
                               inside_at_0p20=bool(nm["0.20"] is not None and nm["0.20"]["E_eV"] <= cbm + PB.DEVICE[1]),
                               inside_at_0p30=bool(nm["0.30"] is not None and nm["0.30"]["E_eV"] <= cbm + PB.DEVICE[1])),
            edges=dict(VB=edge(iv, nocc - 1), CB=edge(ic, nocc), VB_1=edge(iv, nocc - 2), CB_1=edge(ic, nocc + 1),
                       uniform_share=float((left | right).sum() / len(left))),
            pbe_twin=dict(gap=q.get("gap_path"), E_vac_rel_vbm=q.get("E_vac_rel_vbm"),
                          first_nonmodel_rel_cbm=None if not q else (q["nonmodel"]["B"]["0.25"]["E_eV"] - q["cbm"])))
        fnr = r["device_window"]["first_nonmodel_rel_cbm"]
        print(f"=== {key} HSE: {len(dexx)} EXX steps, dexx {dexx[-1]:.1e} Ry; gap (IBZ grid) {r['gap_grid']:.4f} eV "
              f"({'direct' if r['direct'] else 'indirect'}; VBM k {r['k_vbm_frac']}, CBM k {r['k_cbm_frac']} of 2π/a) — "
              f"PBE twin {r['pbe_twin']['gap']:.4f}; E_vac − VBM {r['E_vac_rel_vbm']:.3f} (PBE {r['pbe_twin']['E_vac_rel_vbm']:.3f})")
        for thr, x in nm.items():
            print(f"   first non-model (s_out ≥ {thr}): " + ("none" if x is None else
                  f"CBM + {x['E_eV'] - cbm:.3f} eV, band {x['band']}, k {x['k_index']}, |ψ|² {x['psi2']:.2f}, s_out {x['s_out']}"))
        print(f"   device window (≤ CBM + 0.5): first non-model state {'INSIDE (registered limitation)' if r['device_window']['inside'] else 'outside'}"
              f" (CBM + {fnr:.3f}; PBE twin CBM + {r['pbe_twin']['first_nonmodel_rel_cbm']:.3f})")
        for lab, e in r["edges"].items() if isinstance(r["edges"], dict) else []:
            if lab == "uniform_share":
                continue
            print(f"   {lab:4s} E−VBM {e['E_rel_vbm']:+.3f}, edge weight {e['edge']:.2f} (uniform {r['edges']['uniform_share']:.2f}),"
                  f" |ψ|² {e['psi2']:.2f}, {e['shares']}")
    p.write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    modes = dict(twin_windows=twin_windows, twin_criteria=twin_criteria, hse_points=hse_points)
    if len(sys.argv) > 1 and sys.argv[1] in modes:
        modes[sys.argv[1]](sys.argv[2:] or list(M5.HOLDOUTS))
    else:
        stage2()
