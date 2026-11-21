"""T1 corrections: Fröhlich ballisticity B and the tunnelling guard. Settings: fits/PREREG_corr.md (frozen first).

    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/t1_corr.py run FAM … | dft
    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/t1_corr.py summary

One JSON per ribbon in fits/t1/corr/ (resumable). No band computation for the families: the cached HSE bands of
fits/t1/cache are used.
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from scnfet import frohlich as FR, guard as GD, tob, w90 as W

ROOT = Path(__file__).resolve().parents[1]
T1 = ROOT / "fits" / "t1"
OUT = T1 / "corr"
PREREG, PREREG_MD5 = ROOT / "fits" / "PREREG_corr.md", ROOT / "fits" / "PREREG_corr.md5"
L_NM, VDD, IOFF = 15.0, 0.7, 10e-9 / 1e-6          # nm, V, A/m
LAMS = (1.0, 2.0, 3.0)
GEOMS = ("parallel_plate", "single", "gaa")
IRDS = 1.979
SPAN = 4.0
KT = tob.KB * 300.0


def _frozen():
    md5 = hashlib.md5(PREREG.read_bytes()).hexdigest()
    if not PREREG_MD5.exists() or PREREG_MD5.read_text().split()[0] != md5:
        sys.exit("fits/PREREG_corr.md is not frozen")
    return md5


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, [str(path)]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


def pin(fet, target, VD):
    f = lambda v: np.log(abs(fet.solve(v, VD)[2])) - np.log(target)
    lo, hi = -SPAN, SPAN
    while f(lo) * f(hi) > 0:
        lo, hi = 1.5 * lo, 1.5 * hi
    return brentq(f, lo, hi, xtol=1e-7)


def ribbon(E, nv, Lx_A, W_A, nat, masses, cgrow, stored, label):
    """masses: {'n': m_CB, 'p': m_VB}; cgrow: cg_table entry; stored: {geom: {'n EOT 0.6': I_on, …}}."""
    ch = tob.Channel(E=np.sort(E, axis=1), a=Lx_A * 1e-10, n_val=nv)
    Wm = W_A * 1e-10
    target = IOFF * Wm
    m_r = abs(masses["n"]) * abs(masses["p"]) / (abs(masses["n"]) + abs(masses["p"]))
    out = dict(Eg=float(ch.Eg), W_nm=W_A / 10, m_r=m_r, carriers={})
    for c in ("n", "p"):
        s = 1.0 if c == "n" else -1.0
        VD = s * VDD
        m = abs(masses[c])
        r = dict(m=m)
        # --- OFF state from the ToB (parallel plate, EOT 0.6) and the stored-I_ON check for every geometry
        on, v_off = {}, {}
        for g in GEOMS:
            Cg = None if g == "parallel_plate" else cgrow["EOT 0.6"][g]["C_F_per_m"]
            fet = tob.ToBFET(ch, W=Wm, EOT=0.6e-9, carrier=c, n_atoms_per_cell=nat, C_G=Cg)
            v_off[g] = pin(fet, target, VD)
            on[g] = abs(fet.solve(v_off[g] + s * VDD, VD)[2]) / Wm * 1e-3
            ref = stored[g][f"{c} EOT 0.6"]
            assert abs(on[g] / ref - 1) < 5e-3, f"{label} {c} {g}: pin I_ON {on[g]:.4f} vs stored {ref:.4f}"
            if g == "parallel_plate":
                U_off = fet.solve(v_off[g], VD)[0]
                edge = ch.Ec if c == "n" else ch.Ev
                E_bF = (edge + U_off - fet.muS) if c == "n" else (fet.muS - edge - U_off)
                dS = (fet.muS - edge) if c == "n" else (edge - fet.muS)
            r.setdefault("_fets", {})[g] = fet
        r.update(E_bF=float(E_bF), delta_S=float(dS), I_on_pin=on)
        # --- Fröhlich
        fr = FR.Frohlich(m)
        r["alpha_3D"] = float(fr.alpha)
        r["frohlich"] = {mod: dict(lam0_nm=float(l0), B=float(FR.ballisticity(l0, L_NM)))
                         for mod, l0 in ((mod, fr.lam0(mod)) for mod in ("2d", "3d"))}
        # --- source–drain tunnelling
        r["sdt"] = {f"{lam:g}": dict(zip(("R_tm", "R_wkb"), map(float, GD.sdt(m, lam, E_bF, dS, VDD, L_NM))))
                    for lam in LAMS}
        R2 = r["sdt"]["2"]["R_tm"]
        r["dVT_mV"] = float(1e3 * KT * np.log(R2))
        # re-pin: thermionic I_OFF target / R; subthreshold is exponential in V_GS (ideal electrostatics), so the
        # shift is kT ln R; checked against the target on the parallel-plate device
        rep = {}
        for g in GEOMS:
            fet = r["_fets"][g]
            v2 = v_off[g] - s * KT * np.log(R2)
            if g == "parallel_plate":
                chk = abs(fet.solve(v2, VD)[2]) / (target / R2)
                assert abs(chk - 1) < 0.01, f"{label} {c}: re-pin {chk}"
            rep[g] = abs(fet.solve(v2 + s * VDD, VD)[2]) / Wm * 1e-3
        r["I_on_repinned_lam2"] = rep
        # --- band-to-band
        r["btbt"] = {f"{lam:g}": dict(zip(("I_bound_A", "window_eV"), map(float, GD.gidl_bound(
            m_r, ch.Eg, lam, E_bF, dS, VDD, L_NM)))) for lam in LAMS}
        r["btbt_bound_over_IOFF_lam2"] = r["btbt"]["2"]["I_bound_A"] / target
        # --- corrected I_ON
        B2, B3 = r["frohlich"]["2d"]["B"], r["frohlich"]["3d"]["B"]
        corr = {}
        for g in GEOMS:
            for eot in ("1.0", "0.6"):
                I = stored[g][f"{c} EOT {eot}"]
                corr[f"{g} EOT {eot}"] = dict(ballistic=I, B2d=B2 * I, B3d=B3 * I)
            corr[f"{g} EOT 0.6"]["full"] = B2 * rep[g]
        r["I_on"] = corr
        r["irds"] = {g: dict(ballistic=bool(corr[f"{g} EOT 0.6"]["ballistic"] >= IRDS),
                             B2d=bool(corr[f"{g} EOT 0.6"]["B2d"] >= IRDS),
                             full=bool(corr[f"{g} EOT 0.6"]["full"] >= IRDS)) for g in GEOMS}
        r["flags"] = dict(alpha_gt_1=bool(fr.alpha > 1), W_lt_2nm=bool(W_A < 20), sdt=bool(R2 > 1.5),
                          btbt=bool(r["btbt_bound_over_IOFF_lam2"] > 0.5), gap_lt_0p9=bool(ch.Eg < 0.9))
        del r["_fets"]
        out["carriers"][c] = r
    return out


def stored_rows(fam, n):
    r = json.loads((T1 / "rows" / f"{fam}_n{n}.json").read_text())
    st = {"parallel_plate": {k: v["I_on"] for k, v in r["HSE"]["tob"].items() if k != "irds"}}
    for g in ("single", "gaa"):
        t = json.loads((T1 / "cg" / g / f"{fam}_n{n}.json").read_text())["tob"]
        st[g] = {k: v["I_on"] for k, v in t.items() if k != "irds"}
    return r, st


def run(fams):
    md5 = _frozen()
    OUT.mkdir(parents=True, exist_ok=True)
    tab = json.loads((T1 / "cg" / "cg_table.json").read_text())["widths"]
    for fam in fams:
        if fam == "dft":
            run_dft(md5, tab)
            continue
        for p in sorted((T1 / "rows").glob(f"{fam}_n*.json"), key=lambda p: int(p.stem.split("_n")[1])):
            n = int(p.stem.split("_n")[1])
            f = OUT / f"{fam}_n{n}.json"
            if f.exists():
                continue
            r, st = stored_rows(fam, n)
            z = np.load(T1 / "cache" / f"{fam}_n{n}.npz")
            res = ribbon(z["E"], int(z["nv"]), float(z["Lx"]), float(z["W_A"]), int(z["n_host"]),
                         {"n": r["HSE"]["m_CB"], "p": r["HSE"]["m_VB"]}, tab[f"{fam}_n{n}"], st, f"{fam} n {n}")
            f.write_text(json.dumps(dict(prereg_md5=md5, fam=fam, n=n, **res), indent=1))
            log(f"{fam} n {n:2d}", res)


def run_dft(md5, tab):
    f = OUT / "dft.json"
    if f.exists():
        return
    T = _load("t1", ROOT / "scripts" / "t1_bb_widths.py")
    dft = json.loads((T1 / "dft.json").read_text())["calibration"]
    sg = {g: json.loads((T1 / "cg" / g / "dft.json").read_text())["points"] for g in ("single", "gaa")}
    out = {}
    for fam, c in dft.items():
        st = T.MB.geometry(c["key"])
        E = T.G.w_eigs(W.read_hr(T.G.hr_path(T.MB, c["key"], c["variant"])), T.KG)
        stored = {"parallel_plate": {k: v["I_on"] for k, v in c["HSE"]["tob"].items() if k != "irds"}}
        for g in ("single", "gaa"):
            stored[g] = {k: v["I_on"] for k, v in sg[g][fam]["tob"].items() if k != "irds"}
        res = ribbon(E, T.nval(st), st.cell[0, 0], c["W_A"], T.n_host(st),
                     {"n": c["HSE"]["m_CB"], "p": c["HSE"]["m_VB"]}, tab[f"dft_{fam}"], stored, c["key"])
        out[fam] = dict(key=c["key"], variant=c["variant"], **res)
        log(f"DFT {c['key']}", res)
    f.write_text(json.dumps(dict(prereg_md5=md5, points=out), indent=1))


def log(tag, res):
    for c, r in res["carriers"].items():
        fl = ",".join(k for k, v in r["flags"].items() if v) or "-"
        print(f"  {tag} {c}: m {r['m']:.2f} α {r['alpha_3D']:.2f} λ0 2D/3D {r['frohlich']['2d']['lam0_nm']:.1f}/"
              f"{r['frohlich']['3d']['lam0_nm']:.1f} nm B {r['frohlich']['2d']['B']:.2f}/{r['frohlich']['3d']['B']:.2f} | "
              f"E_bF {r['E_bF']:.3f} R(1,2,3) {r['sdt']['1']['R_tm']:.2f}/{r['sdt']['2']['R_tm']:.2f}/{r['sdt']['3']['R_tm']:.2f} "
              f"BTBT/I_OFF {r['btbt_bound_over_IOFF_lam2']:.1e} | I_ON pp {r['I_on']['parallel_plate EOT 0.6']['ballistic']:.2f}"
              f"→B {r['I_on']['parallel_plate EOT 0.6']['B2d']:.2f}→full {r['I_on']['parallel_plate EOT 0.6']['full']:.2f}; "
              f"single {r['I_on']['single EOT 0.6']['ballistic']:.2f}→{r['I_on']['single EOT 0.6']['B2d']:.2f} | flags {fl}",
              flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "run":
        run(sys.argv[2:] or ["100", "110", "110h", "dft"])
    else:
        sys.exit(f"unknown mode {mode}")
