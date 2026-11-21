"""HSE cost tests: the decisions pre-registered in README.md (written before any output).

    cd dft/A2/costtest && PYTHONPATH=../../.. ~/.venvs/scnfet/bin/python -P parse_costtest.py
(after copying back outputs/ and w90/*/{*.wout,*_hr.dat}). Writes costtest.json.
"""
import importlib.util
import json
import re
from pathlib import Path

import numpy as np

from scnfet import w90

HERE = Path(__file__).resolve().parent
A2 = HERE.parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


P = _load("parse_a2", A2 / "parse_a2.py")
M = _load("make_a2", A2 / "make_a2.py")
T = _load("qe_timing", A2.parent / "A4a" / "qe_timing.py")
MK = _load("make_costtest", HERE / "make_costtest.py")
CELL = np.diag([M.A, M.A, M.C])
WIN = json.loads((A2 / "a2_windows.json").read_text())
S1 = json.loads((A2 / "a2_stage1.json").read_text())
QUAD, LIN = (40 / 24) ** 2, 40 / 24
BOUND = float(np.sqrt(QUAD * LIN))                        # 2.15


def t1():
    rows = {}
    runs = {tag: (HERE / "outputs" / f"{tag}.out", HERE / "inputs" / f"{tag}.in") for tag in MK.T1}
    runs["A2_hse_k8_q4"] = (A2 / "outputs" / "hse_k8_q4.out", A2 / "inputs" / "hse_k8_q4.in")
    for tag, (out, inp) in runs.items():
        d = T.parse(out, inp)
        txt = out.read_text()
        node = re.search(r"Parallel version.*?\n", txt)
        rows[tag] = dict(nbnd=d["nbnd"], nks=d["nks"], nq=d.get("nq"), outer=d.get("outer"), vexx_s=d.get("vexx"),
                         vexx_calls=d.get("vexx_calls"), t_call_s=d["vexx"] / d["vexx_calls"], wall_s=d["wall"],
                         scf_iters=d["iters"], E_tot_Ry=float(re.findall(r"!!?\s+total energy\s+=\s+(-?[\d.]+)", txt)[-1]))
    r_fix = rows["ct_fix40"]["t_call_s"] / rows["ct_fix24"]["t_call_s"]
    r_sm = rows["ct_sm40"]["t_call_s"] / rows["ct_sm24"]["t_call_s"]
    # eigenvalues: smearing vs fixed at the same nbnd, all bands and k, VBM aligned
    eq = {}
    for nb in (24, 40):
        _, ef, _ = P.bands(HERE / "outputs" / f"ct_fix{nb}.out")
        _, es, _ = P.bands(HERE / "outputs" / f"ct_sm{nb}.out")
        sh = ef[:, P.NV - 1].max() - es[:, P.NV - 1].max()
        eq[f"nbnd{nb}"] = dict(max_dE_meV=1e3 * float(np.abs(ef - es - sh).max()), vbm_shift_meV=1e3 * float(sh),
                               gap_fix=float(ef[:, P.NV].min() - ef[:, P.NV - 1].max()),
                               gap_sm=float(es[:, P.NV].min() - es[:, P.NV - 1].max()))
    same = all(v["max_dE_meV"] <= 1.0 for v in eq.values())
    return dict(rows=rows, r_fix=r_fix, r_sm=r_sm, quadratic=QUAD, linear=LIN, boundary=BOUND,
                fixed_is_quadratic=bool(r_fix >= BOUND), smearing_is_linear=bool(r_sm < BOUND),
                eigen_fixed_vs_smearing=eq, eigen_same=bool(same),
                smearing_removes_empty_band_cost=bool(r_sm < BOUND and same),
                control_vs_A2=rows["ct_fix24"]["t_call_s"] / rows["A2_hse_k8_q4"]["t_call_s"])


def centres(wout):
    txt = Path(wout).read_text()
    blk = txt[txt.rfind("Final State"):]
    c = np.array([[float(v) for v in m] for m in re.findall(
        r"WF centre and spread\s+\d+\s+\(\s*(-?[\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+)\s*\)\s+([\d.]+)", blk)])
    om = float(re.search(r"Omega Total\s*=\s*([\d.]+)", blk).group(1))
    z = M.zfrac()
    atoms = np.array([(0, 0, z["Sc1"]), (0.5, 0.5, z["N1"]), (0, 0, z["N2"]), (0.5, 0.5, z["Sc2"])]) @ CELL
    inv = np.linalg.inv(CELL)
    off, dz = [], []
    for r in c[:, :3]:
        d = []
        for a in atoms:
            f = (r - a) @ inv
            f -= np.round(f)
            d.append(f @ CELL)
        j = int(np.argmin([np.linalg.norm(v) for v in d]))
        off.append(float(np.linalg.norm(d[j])))
        dz.append(float(d[j][2]))
    return dict(max_offset_A=float(max(off)), sc_max_abs_dz_A=float(np.abs(dz[:10]).max()),
                dz_A=[round(v, 4) for v in dz], spreads_A2=c[:, 3].round(4).tolist(),
                max_spread_A2=float(c[:, 3].max()), omega_total=om)


def interp(hr_path, kcart):
    hr = w90.read_hr(hr_path)
    return np.linalg.eigvalsh(w90.hk_frac(hr, kcart @ CELL.T / (2 * np.pi)))


def t2():
    win = json.loads((HERE / "windows.json").read_text())
    kp, ep, _ = P.bands(A2 / "outputs" / "pbe_bands.out")
    refp = ep[:, P.NEXCL:P.NEXCL + P.NMODEL]
    vbm_p = S1["pbe_bands"]["vbm"]
    kh, eh, _ = P.bands(A2 / "outputs" / "hse_ibz_k16.out")
    refh = eh[:, P.NEXCL:P.NEXCL + P.NMODEL]
    vbm_h, gap_h = S1["hse_ibz_k16"]["vbm"], S1["hse_ibz_k16"]["gap"]
    base_h = interp(A2 / "w90" / "hse_k16" / "hse_k16_hr.dat", kp)        # baseline HSE W90 on the PBE path
    # VB / CB-to-CBM+1 masks on the path, from the baseline HSE W90 bands
    edge_h = (base_h - vbm_h <= 0.01) | ((base_h - vbm_h >= gap_h - 0.01) & (base_h - vbm_h <= gap_h + 1.0))

    def pbe_metrics(hr_path):
        E = interp(hr_path, kp)
        inw = refp <= WIN["pbe_k16"]["froz"]
        err = np.abs(E - refp)
        gap_p = S1["pbe_bands"]["gap"]
        vb = inw & (refp - vbm_p <= 0.01)
        cb = inw & (refp - vbm_p >= gap_p - 0.01) & (refp - vbm_p <= gap_p + 1.0)
        return dict(rms_frozen_meV=1e3 * float(np.sqrt(np.mean(err[inw] ** 2))), max_frozen_meV=1e3 * float(err[inw].max()),
                    max_vb_meV=1e3 * float(err[vb].max()), max_cb_le_cbm_plus_1_meV=1e3 * float(err[cb].max()))

    def hse_metrics(hr_path):
        Eg = interp(hr_path, kh)
        inw = refh <= WIN["hse_k16"]["froz"]
        Ep = interp(hr_path, kp)
        return dict(grid_max_meV=1e3 * float(np.abs(Eg - refh)[inw].max()),
                    vs_baseline_edges_max_meV=1e3 * float(np.abs(Ep - base_h)[edge_h].max()))

    out = {"baseline_pbe_k16": dict(outer_rel_vbm=WIN["pbe_k16"]["outer_rel"], nbnd_equiv=44,
                                    **pbe_metrics(A2 / "w90" / "pbe_k16" / "pbe_k16_hr.dat"),
                                    centres=centres(A2 / "w90" / "pbe_k16" / "pbe_k16.wout")),
           "baseline_hse_k16": dict(outer_rel_vbm=WIN["hse_k16"]["outer_rel"], nbnd_equiv=48,
                                    **hse_metrics(A2 / "w90" / "hse_k16" / "hse_k16_hr.dat"),
                                    centres=centres(A2 / "w90" / "hse_k16" / "hse_k16.wout"))}
    for v, r in win.items():
        d = HERE / "w90" / v
        hr = d / f"{r['seed']}_hr.dat"
        row = dict(nbnd=r["nbnd"], outer_rel_vbm=r["outer_rel_vbm"])
        if not hr.exists():
            row["missing"] = True
        else:
            row.update(pbe_metrics(hr) if r["seed"] == "pbe_k16" else hse_metrics(hr))
            row["centres"] = centres(d / f"{r['seed']}.wout")
        out[v] = row
    # acceptance (README): PBE twin (1) + HSE (2), (3) at the same N
    acc = {}
    for N in MK.NBND:
        p = out.get(f"pbe_k16_nb{N}", out["baseline_pbe_k16"] if N == 44 else None)
        h = out.get(f"hse_k16_nb{N}")
        if not p or not h or p.get("missing") or h.get("missing"):
            continue
        ok1 = p["max_vb_meV"] <= 2 and p["max_cb_le_cbm_plus_1_meV"] <= 2 and p["rms_frozen_meV"] <= 3
        ok2 = h["grid_max_meV"] <= 1
        ok3 = h["vs_baseline_edges_max_meV"] <= 5
        acc[N] = dict(ok=bool(ok1 and ok2 and ok3), pbe_twin=bool(ok1), hse_grid=bool(ok2), hse_vs_baseline=bool(ok3),
                      hse_outer_rel_vbm=h["outer_rel_vbm"])
    good = [N for N, a in acc.items() if a["ok"]]
    out["acceptance"] = acc
    out["smallest_acceptable"] = dict(nbnd=min(good), hse_outer_rel_vbm=acc[min(good)]["hse_outer_rel_vbm"]) if good else None
    return out


def main():
    res = dict(T1=t1(), T2=t2())
    (HERE / "costtest.json").write_text(json.dumps(res, indent=1))
    a = res["T1"]
    print(f"T1: r_fix {a['r_fix']:.2f}  r_sm {a['r_sm']:.2f}  (quadratic {QUAD:.2f}, linear {LIN:.2f}, boundary "
          f"{BOUND:.2f}); eigen same: {a['eigen_same']}; smearing removes empty-band cost: "
          f"{a['smearing_removes_empty_band_cost']}; control/A2 {a['control_vs_A2']:.2f}")
    for t, r in a["rows"].items():
        print(f"   {t:13s} nbnd {r['nbnd']} outer {r['outer']} t_call {r['t_call_s']:.3f}s vexx {r['vexx_s']:.0f}s wall "
              f"{r['wall_s']:.0f}s")
    for k, v in res["T2"].items():
        if k in ("acceptance", "smallest_acceptable"):
            print(k, v)
        else:
            c = v.get("centres", {})
            print(f"   {k:17s} " + " ".join(f"{kk} {vv:.2f}" for kk, vv in v.items() if isinstance(vv, float))
                  + (f"  Sc|dz|max {c['sc_max_abs_dz_A']:.2f} Å, spread max {c['max_spread_A2']:.2f}, Ω {c['omega_total']:.2f}"
                     if c else ""))


if __name__ == "__main__":
    main()
