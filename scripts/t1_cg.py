"""T1 with the geometric gate capacitance (Laplace C_G). Settings: fits/PREREG_cg.md (frozen before any device number).

Modes (PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/t1_cg.py <mode>):
    table           C_G(W) for every T1 width and the calibration DFT points: main geometries + brackets
    cache FAM …     HSE building-block eigenvalues on the 2000-k ToB grid (same code as T1), with a check that the
                    parallel-plate I_ON of the first width reproduces fits/t1/rows exactly
    tob GEOM        ToB with the Laplace C_G (GEOM = single | gaa) from the cache, plus the calibration DFT points
    figs            figures/Fig7_Ion_W.png (parallel plate, single gate, GAA) and fits/t1/cg/summary.md
Outputs in fits/t1/cg/.
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from scnfet import cg, tob, w90 as W

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "fits" / "t1" / "cg"
CACHE = ROOT / "fits" / "t1" / "cache"
PREREG, PREREG_MD5 = ROOT / "fits" / "PREREG_cg.md", ROOT / "fits" / "PREREG_cg.md5"
EPS_HFO2, D_RIB = 25.0, 0.217
EOTS = (1.0, 0.6)
FAMS = ("100", "110", "110h")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, [str(path)]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


def _frozen():
    md5 = hashlib.md5(PREREG.read_bytes()).hexdigest()
    if not PREREG_MD5.exists() or PREREG_MD5.read_text().split()[0] != md5:
        sys.exit("fits/PREREG_cg.md is not frozen")
    return md5


def rows(fam):
    return sorted((json.loads(p.read_text()) for p in (ROOT / "fits/t1/rows").glob(f"{fam}_n*.json")),
                  key=lambda r: r["n"])


def geometries(W_nm, eot):
    t = eot * EPS_HFO2 / 3.9
    pp = cg.parallel_plate(W_nm, eot)
    out = dict(parallel_plate=pp,
               single=cg.cg_single(W_nm, t, d=D_RIB, eps_below=EPS_HFO2, eps_above=EPS_HFO2),
               gaa=cg.cg_gaa(W_nm, t, d=D_RIB, eps=EPS_HFO2),
               single_vacuum_above=cg.cg_single(W_nm, t, d=D_RIB, eps_below=EPS_HFO2, eps_above=1.0),
               single_sio2=cg.cg_single(W_nm, eot, d=D_RIB, eps_below=3.9, eps_above=3.9),
               gaa_sio2=cg.cg_gaa(W_nm, eot, d=D_RIB, eps=3.9))
    return {k: dict(C_F_per_m=v, ratio_pp=v / pp) for k, v in out.items()}


def table():
    md5 = _frozen()
    res = dict(prereg_md5=md5, eps_hfo2=EPS_HFO2, d_nm=D_RIB, widths={})
    for fam in FAMS:
        for r in rows(fam):
            W_nm = r["W_A"] / 10
            res["widths"][f"{fam}_n{r['n']}"] = dict(W_nm=W_nm, **{f"EOT {e}": geometries(W_nm, e) for e in EOTS})
    dft = json.loads((ROOT / "fits/t1/dft.json").read_text())
    for fam, c in dft["calibration"].items():
        W_nm = c["W_A"] / 10
        res["widths"][f"dft_{fam}"] = dict(W_nm=W_nm, key=c["key"], variant=c["variant"],
                                           **{f"EOT {e}": geometries(W_nm, e) for e in EOTS})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "cg_table.json").write_text(json.dumps(res, indent=1))
    for k, v in res["widths"].items():
        print(f"  {k:10s} W {v['W_nm']:.2f} nm | " + " | ".join(
            f"EOT {e}: single {v[f'EOT {e}']['single']['ratio_pp']:.2f} gaa {v[f'EOT {e}']['gaa']['ratio_pp']:.2f} "
            f"(vac-above {v[f'EOT {e}']['single_vacuum_above']['ratio_pp']:.2f}, SiO2 {v[f'EOT {e}']['single_sio2']['ratio_pp']:.2f}"
            f"/{v[f'EOT {e}']['gaa_sio2']['ratio_pp']:.2f})" for e in EOTS), flush=True)


def cache(fams):
    _frozen()
    T = _load("t1", ROOT / "scripts" / "t1_bb_widths.py")
    CACHE.mkdir(parents=True, exist_ok=True)
    for fam in fams:
        for i, r in enumerate(rows(fam)):
            f = CACHE / f"{fam}_n{r['n']}.npz"
            if f.exists():
                continue
            st0, Tr, hr, wf, info = T.bb_hr(fam, r["n"], "HSE06")
            E = T.G.w_eigs(hr, T.KG)
            np.savez_compressed(f, E=E, nv=r["nv"], Lx=r["Lx"], W_A=r["W_A"], n_host=r["n_host"])
            if i == 0:          # reproduce the stored parallel-plate I_ON exactly
                chk = T.G.ion(E, r["nv"], r["Lx"], r["W_A"], r["n_host"])["n EOT 0.6"]["I_on"]
                ref = r["HSE"]["tob"]["n EOT 0.6"]["I_on"]
                assert abs(chk / ref - 1) < 1e-9, f"{fam} n {r['n']}: cache I_ON {chk} vs rows {ref}"
                print(f"  {fam} n {r['n']}: cache reproduces the parallel-plate I_ON ({chk:.4f})", flush=True)
            print(f"  cached {fam} n {r['n']} ({E.shape[1]} bands)", flush=True)


def tob_geom(geom):
    md5 = _frozen()
    T = _load("t1", ROOT / "scripts" / "t1_bb_widths.py")
    tab = json.loads((OUT / "cg_table.json").read_text())["widths"]
    d = OUT / geom
    d.mkdir(parents=True, exist_ok=True)

    def run(E, nv, Lx_A, W_A, nat, key):
        ch = tob.Channel(E=np.sort(E, axis=1), a=Lx_A * 1e-10, n_val=nv)
        res = {}
        for c in ("n", "p"):
            for eot in EOTS:
                Cg = tab[key][f"EOT {eot}"][geom]["C_F_per_m"]
                f = tob.ToBFET(ch, W=W_A * 1e-10, EOT=eot * 1e-9, carrier=c, n_atoms_per_cell=nat, C_G=Cg)
                m = f.metrics(span=T.G.TOB_SPAN, n_sweep=T.G.TOB_NSWEEP)
                res[f"{c} EOT {eot}"] = dict(I_on=float(m["I_on"]), SS=float(m["SS_mV_dec"]), v_inj=float(m["v_inj"]),
                                             C_ratio=tab[key][f"EOT {eot}"][geom]["ratio_pp"])
        res["irds"] = {c: bool(res[f"{c} EOT 0.6"]["I_on"] >= T.G.IRDS_ION) for c in ("n", "p")}
        return res

    for fam in FAMS:
        for r in rows(fam):
            f = d / f"{fam}_n{r['n']}.json"
            if f.exists():
                continue
            z = np.load(CACHE / f"{fam}_n{r['n']}.npz")
            res = run(z["E"], int(z["nv"]), float(z["Lx"]), float(z["W_A"]), int(z["n_host"]), f"{fam}_n{r['n']}")
            f.write_text(json.dumps(dict(prereg_md5=md5, fam=fam, n=r["n"], W_A=r["W_A"], geom=geom, tob=res), indent=1))
            print(f"  {geom} {fam} n {r['n']:2d}: I_ON n {res['n EOT 1.0']['I_on']:.2f}/{res['n EOT 0.6']['I_on']:.2f} "
                  f"p {res['p EOT 1.0']['I_on']:.2f}/{res['p EOT 0.6']['I_on']:.2f} (C/C_pp {res['n EOT 0.6']['C_ratio']:.2f})",
                  flush=True)
    dft = json.loads((ROOT / "fits/t1/dft.json").read_text())
    f = d / "dft.json"
    if not f.exists():
        out = {}
        for fam, c in dft["calibration"].items():
            kc, var = c["key"], c["variant"]
            st = T.MB.geometry(kc)
            hr = W.read_hr(T.G.hr_path(T.MB, kc, var))
            E = T.G.w_eigs(hr, T.KG)
            out[fam] = dict(key=kc, variant=var, W_A=c["W_A"],
                            tob=run(E, T.nval(st), st.cell[0, 0], c["W_A"], T.n_host(st), f"dft_{fam}"))
            print(f"  {geom} DFT {kc}: I_ON n/p (0.6) {out[fam]['tob']['n EOT 0.6']['I_on']:.2f}/"
                  f"{out[fam]['tob']['p EOT 0.6']['I_on']:.2f}", flush=True)
        f.write_text(json.dumps(dict(prereg_md5=md5, points=out), indent=1))


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "table":
        table()
    elif mode == "cache":
        cache(sys.argv[2:] or FAMS)
    elif mode == "tob":
        tob_geom(sys.argv[2])
    else:
        sys.exit(f"unknown mode {mode}")
