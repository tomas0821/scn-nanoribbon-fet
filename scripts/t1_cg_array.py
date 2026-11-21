"""T1 with the C_G of a periodic ribbon array (single planar gate). Settings: fits/PREREG_array.md (frozen first).

Modes (PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/t1_cg_array.py <mode>):
    table        C_G per ribbon at s = 2.0 nm and s = 2 t_HfO2 (EOT 0.6), every T1 width + the calibration DFT points
    tob S        ToB at EOT 0.6 with the array C_G (S = s2 | sGAA), from fits/t1/cache, plus the DFT points
Outputs in fits/t1/cg_array/.
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from scnfet import cg, tob, w90 as W
from scnfet.cg_array import cg_single_array

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "fits" / "t1" / "cg_array"
CACHE = ROOT / "fits" / "t1" / "cache"
PREREG, PREREG_MD5 = ROOT / "fits" / "PREREG_array.md", ROOT / "fits" / "PREREG_array.md5"
EPS_HFO2, D_RIB, EOT = 25.0, 0.217, 0.6
T_HFO2 = EOT * EPS_HFO2 / 3.9
SPACINGS = {"s2": 2.0, "sGAA": 2 * T_HFO2}
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
        sys.exit("fits/PREREG_array.md is not frozen")
    return md5


def rows(fam):
    return sorted((json.loads(p.read_text()) for p in (ROOT / "fits/t1/rows").glob(f"{fam}_n*.json")),
                  key=lambda r: r["n"])


def entry(W_nm):
    pp = cg.parallel_plate(W_nm, EOT)
    iso = cg.cg_single(W_nm, T_HFO2, d=D_RIB, eps_below=EPS_HFO2, eps_above=EPS_HFO2)
    out = dict(W_nm=W_nm, C_pp=pp, C_iso=iso)
    for k, s in SPACINGS.items():
        c = cg_single_array(W_nm, s, T_HFO2, d=D_RIB, eps_below=EPS_HFO2, eps_above=EPS_HFO2)
        out[k] = dict(s_nm=s, C_F_per_m=c, ratio_pp=c / pp, ratio_iso=c / iso)
    return out


def table():
    md5 = _frozen()
    res = dict(prereg_md5=md5, eot=EOT, t_hfo2_nm=T_HFO2, spacings_nm=SPACINGS, widths={})
    for fam in FAMS:
        for r in rows(fam):
            res["widths"][f"{fam}_n{r['n']}"] = entry(r["W_A"] / 10)
    dft = json.loads((ROOT / "fits/t1/dft.json").read_text())
    for fam, c in dft["calibration"].items():
        res["widths"][f"dft_{fam}"] = dict(key=c["key"], variant=c["variant"], **entry(c["W_A"] / 10))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "cg_array_table.json").write_text(json.dumps(res, indent=1))
    for k, v in res["widths"].items():
        print(f"  {k:10s} W {v['W_nm']:.2f} | iso/pp {v['C_iso'] / v['C_pp']:.2f} | " + " | ".join(
            f"{s}: /pp {v[s]['ratio_pp']:.2f} /iso {v[s]['ratio_iso']:.2f}" for s in SPACINGS), flush=True)


def tob_s(skey):
    md5 = _frozen()
    T = _load("t1", ROOT / "scripts" / "t1_bb_widths.py")
    tab = json.loads((OUT / "cg_array_table.json").read_text())["widths"]
    d = OUT / skey
    d.mkdir(parents=True, exist_ok=True)

    def run(E, nv, Lx_A, W_A, nat, key):
        ch = tob.Channel(E=np.sort(E, axis=1), a=Lx_A * 1e-10, n_val=nv)
        res = {}
        Cg = tab[key][skey]["C_F_per_m"]
        for c in ("n", "p"):
            f = tob.ToBFET(ch, W=W_A * 1e-10, EOT=EOT * 1e-9, carrier=c, n_atoms_per_cell=nat, C_G=Cg)
            m = f.metrics(span=T.G.TOB_SPAN, n_sweep=T.G.TOB_NSWEEP)
            res[f"{c} EOT {EOT}"] = dict(I_on=float(m["I_on"]), SS=float(m["SS_mV_dec"]), v_inj=float(m["v_inj"]),
                                         C_ratio_pp=tab[key][skey]["ratio_pp"], C_ratio_iso=tab[key][skey]["ratio_iso"])
        res["irds"] = {c: bool(res[f"{c} EOT {EOT}"]["I_on"] >= T.G.IRDS_ION) for c in ("n", "p")}
        return res

    for fam in FAMS:
        for r in rows(fam):
            f = d / f"{fam}_n{r['n']}.json"
            if f.exists():
                continue
            z = np.load(CACHE / f"{fam}_n{r['n']}.npz")
            res = run(z["E"], int(z["nv"]), float(z["Lx"]), float(z["W_A"]), int(z["n_host"]), f"{fam}_n{r['n']}")
            f.write_text(json.dumps(dict(prereg_md5=md5, fam=fam, n=r["n"], W_A=r["W_A"], spacing=skey,
                                         s_nm=SPACINGS[skey], tob=res), indent=1))
            print(f"  {skey} {fam} n {r['n']:2d}: I_ON n {res['n EOT 0.6']['I_on']:.2f} p {res['p EOT 0.6']['I_on']:.2f} "
                  f"(C/C_pp {res['n EOT 0.6']['C_ratio_pp']:.2f})", flush=True)
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
            print(f"  {skey} DFT {kc}: I_ON n/p {out[fam]['tob']['n EOT 0.6']['I_on']:.2f}/"
                  f"{out[fam]['tob']['p EOT 0.6']['I_on']:.2f}", flush=True)
        f.write_text(json.dumps(dict(prereg_md5=md5, points=out), indent=1))


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "table":
        table()
    elif mode == "tob":
        tob_s(sys.argv[2])
    else:
        sys.exit(f"unknown mode {mode}")
