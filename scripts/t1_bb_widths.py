"""T1: building-block device Hamiltonians across 1–6 nm. Settings: fits/PREREG_t1.md (frozen before any width result).

Modes (run from the project root; PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/t1_bb_widths.py <mode>):
    hcalib              the H rule: H calibration ribbons rebuilt from their own halves + the sheet (PBE W1, HSE HW1)
    geomcheck           the D4-copy geometry against the two relaxed hold-outs
    sweep FAM [FAM …]   widths of a family (100, 110, 100h, 110h): HSE bands + ToB, PBE bands. Resumable.
    dft                 overlay points: HSE device Hamiltonians of the calibration ribbons (ToB), DFT gaps
    sk                  SK-HSE overlay on bare [100] (gate (d) recipe), bands + ToB
    figs                figures/Fig4_Eg_W.png, Fig5_valley_mass_W.png, Fig7_Ion_W_parallel_plate.png; fits/t1/t1_table.md
Outputs in fits/t1/ (rows/<fam>_n<n>.json, sk/, dft.json, hcalib.json, geomcheck.json).
"""
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

from scnfet import bb, sk, structures as S, w90 as W

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "fits" / "t1"
PREREG, PREREG_MD5 = ROOT / "fits" / "PREREG_t1.md", ROOT / "fits" / "PREREG_t1.md5"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, [str(path)]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


B = _load("gate_bb", ROOT / "scripts" / "gate_bb_holdouts.py")       # the frozen building-block recipe
G, MB, MBH, M5 = B.G, B.MB, B.MBH, B.M5
MK = MB.MK
FAM = {"100": ("100", False, "r100_3", range(3, 15)), "110": ("110", False, "r110_4", range(4, 21)),
       "100h": ("100", True, "r100_3h", range(3, 15)), "110h": ("110", True, "r110_4h", range(4, 21))}
FUNC_VAR = {"HSE06": "HW1", "PBE": "W1"}                                 # the recipe's edge variants
DEVICE_VAR = {"r100_3": "HW3", "r110_4": "HW1", "r100_3h": "HW2", "r110_4h": "HW2"}   # A4b-2 device Hamiltonians
LABEL = {"100": "[100] bare", "110": "[110] bare", "100h": "[100] H", "110h": "[110] H"}
NK_PATH, NK_TOB = G.NK_PATH, G.NK_TOB
KF = np.c_[np.linspace(0, 0.5, NK_PATH), np.zeros((NK_PATH, 2))]
KG = np.c_[np.arange(NK_TOB) / NK_TOB, np.zeros((NK_TOB, 2))]
H_LIMITS = dict(gap_meV=50.0, valley_meV=50.0, mass_rel=0.10)


def _frozen():
    md5 = hashlib.md5(PREREG.read_bytes()).hexdigest()
    if not PREREG_MD5.exists() or PREREG_MD5.read_text().split()[0] != md5:
        sys.exit("fits/PREREG_t1.md is not frozen (md5 missing or changed)")
    return md5


# ------------------------------------------------------------------------------------------------- geometry
def target(fam, n):
    d, h, _, _ = FAM[fam]
    st0 = S.ribbon(d, n, a_p=MK.A_P, d=MK.D_PLANE, vacuum=20.0, passivate=h)
    return st0, bb.Ribbon(pos=st0.pos, species=list(st0.species), Lx=float(st0.cell[0, 0]))


def cal_disp(kc):
    """The calibration ribbon's relaxation displacements (relaxed − ideal cut), rigid shift (host mean) removed."""
    st_i, st_r = MK.structure(kc), MB.geometry(kc)
    d = st_r.pos - st_i.pos
    L = st_i.cell[0, 0]
    d[:, 0] -= L * np.round(d[:, 0] / L)
    host = np.array([s != "H" for s in st_i.species])
    return d - d[host].mean(axis=0)


def bb_geometry(st0, T, kc):
    """D4's copy step: the ideal cut with the matched calibration atoms' displacements on the edge-region atoms."""
    C, _ = B.ideal(MK, kc)
    disp = cal_disp(kc)
    pos = st0.pos.copy()
    for a, (side, c) in enumerate(bb.edge_map(T, C)):
        if c is not None:
            pos[a] += disp[c]
    return sk.Structure(st0.cell, list(st0.species), pos, pbc=st0.pbc, frame=st0.frame)


def n_host(st):
    return sum(1 for s in st.species if s != "H")


def nval(st):
    return MB.nocc(st) - MB.nexcl(st)


# ------------------------------------------------------------------------------------------------- building blocks
_SRC = {}


def sources(kc, func):
    key = (kc, func)
    if key not in _SRC:
        sheet = B.sheet_source(func)
        csrc, var = B.cal_source(kc, func, FUNC_VAR[func])
        C, _ = B.ideal(MK, kc)
        _SRC[key] = (sheet, csrc, C, bb.centre_shift(C, csrc, sheet))
    return _SRC[key]


def bb_hr(fam, n, func):
    st0, T = target(fam, n)
    kc = FAM[fam][2]
    sheet, csrc, C, shift = sources(kc, func)
    emap = bb.edge_map(T, C)
    s_atoms = [a for a, (_, c) in enumerate(emap) if c is not None and "s" in csrc.idx.get(c, {})
               and T.species[a] == "Sc"]
    hr, wf, used = bb.assemble(T, C, csrc, sheet, B.q_sheet(st0), B.RCUT, s_atoms=s_atoms, sheet_shift=shift)
    return st0, T, hr, wf, dict(sheet_shift_eV=shift, n_s=len(s_atoms), n_sheet_pairs=sum(v == "S" for v in used.values()),
                                n_pairs=len(used))


def tob_all(E_grid, nv, Lx, Wd, nat):
    res = G.ion(E_grid, nv, Lx, Wd, nat)
    res["irds"] = {c: bool(res[f"{c} EOT 0.6"]["I_on"] >= G.IRDS_ION) for c in ("n", "p")}
    return res


# ------------------------------------------------------------------------------------------------- modes
def hcalib():
    """The H rule (PREREG_t1): rebuild each H calibration ribbon from its own halves plus the sheet, against itself."""
    _frozen()
    res = {}
    for fam in ("100h", "110h"):
        kc = FAM[fam][2]
        ok_all = True
        for func in ("PBE", "HSE06"):
            sheet, csrc, C, shift = sources(kc, func)
            st = MB.geometry(kc)
            nv = nval(st)
            Ew = B.G.w_eigs(csrc.hr, KF) - csrc.E_vac
            ref = G.ribbon_metrics(Ew, C.Lx, nv)
            hr, wf, used = bb.assemble(C, C, csrc, sheet, B.q_sheet(MK.structure(kc)), B.RCUT, sheet_shift=shift)
            Eb = G.w_eigs(hr, KF)
            m = G.ribbon_metrics(Eb, C.Lx, nv)
            d = G.compare(m, ref)
            ok = bool(abs(d["gap_meV"]) <= H_LIMITS["gap_meV"] and abs(d["valley_meV"]) <= H_LIMITS["valley_meV"]
                      and max(abs(d["m_VB_rel"]), abs(d["m_CB_rel"])) <= H_LIMITS["mass_rel"])
            ok_all &= ok
            res[f"{kc} {func} {FUNC_VAR[func]}"] = dict(own=ref, rebuild=m, d=d, ok=ok, sheet_shift_eV=shift,
                                                         window_dev_meV=B.window_dev(Eb, Ew, nv),
                                                         n_sheet_pairs=sum(v == "S" for v in used.values()))
            print(f"  {kc} {func} ({FUNC_VAR[func]}, shift {shift:+.3f} eV): gap {m['gap']:.3f} vs {ref['gap']:.3f} "
                  f"({d['gap_meV']:+.0f} meV) valley {d['valley_meV']:+.0f} m_VB {d['m_VB_rel']:+.1%} "
                  f"m_CB {d['m_CB_rel']:+.1%} → {'ok' if ok else 'FAIL'}")
        res[f"{fam} curve"] = ok_all
        print(f"  {LABEL[fam]}: {'labelled building-block curve' if ok_all else 'NO curve (DFT point only)'}")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "hcalib.json").write_text(json.dumps(res, indent=1, default=float))


def geomcheck():
    _frozen()
    res = {}
    for fam, kh, n in (("110", "r110_8", 8), ("100", "r100_6", 6)):
        st0, T = target(fam, n)
        g = bb_geometry(st0, T, FAM[fam][2])
        rel = M5.relaxed(kh)
        ideal = M5.structure(kh)
        d1 = g.pos - st0.pos
        d2 = rel.pos - ideal.pos
        L = st0.cell[0, 0]
        d2[:, 0] -= L * np.round(d2[:, 0] / L)
        d2 -= d2.mean(axis=0)
        dev = np.linalg.norm(d1 - d2, axis=1)
        res[kh] = dict(W_bb_A=S.width(g), W_relaxed_A=S.width(rel), W_ideal_A=S.width(st0),
                       max_atom_dev_A=float(dev.max()), edge_region_atoms=int(sum(c is not None for _, c in
                                                                                  bb.edge_map(T, B.ideal(MK, FAM[fam][2])[0]))))
        print(f"  {kh}: W copy-step {S.width(g):.3f} Å vs relaxed {S.width(rel):.3f} Å (ideal {S.width(st0):.3f}); "
              f"max atom deviation {dev.max():.3f} Å")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "geomcheck.json").write_text(json.dumps(res, indent=1))


def sweep(fams):
    md5 = _frozen()
    hc = json.loads((OUT / "hcalib.json").read_text())
    rows = OUT / "rows"
    rows.mkdir(parents=True, exist_ok=True)
    for fam in fams:
        if fam.endswith("h") and not hc[f"{fam} curve"]:
            print(f"  {fam}: skipped by the H rule")
            continue
        for n in FAM[fam][3]:
            f = rows / f"{fam}_n{n}.json"
            if f.exists():
                continue
            t0 = time.time()
            st0, T, hr, wf, info = bb_hr(fam, n, "HSE06")
            geo = bb_geometry(st0, T, FAM[fam][2])
            Wd, Lx, nv, nat = S.width(geo), T.Lx, nval(st0), n_host(st0)
            E, sw = B.eig_s(hr, wf, KF)
            m = G.ribbon_metrics(E, Lx, nv, sw)
            Eg = G.w_eigs(hr, KG)
            row = dict(prereg_md5=md5, fam=fam, n=n, W_A=Wd, W_ideal_A=S.width(st0), Lx=Lx, nv=nv, n_host=nat,
                       nw=len(wf), build_HSE=info, HSE=m)
            row["HSE"]["tob"] = tob_all(Eg, nv, Lx, Wd, nat)
            del Eg
            st0, T, hr, wf, info = bb_hr(fam, n, "PBE")
            E, sw = B.eig_s(hr, wf, KF)
            row["build_PBE"], row["PBE"] = info, G.ribbon_metrics(E, Lx, nv, sw)
            row["seconds"] = time.time() - t0
            f.write_text(json.dumps(row, indent=1, default=float))
            t = row["HSE"]["tob"]
            print(f"  {fam} n {n:2d} W {Wd / 10:.2f} nm: HSE gap {m['gap']:.3f} valley {m['valley_XG']:+.3f} "
                  f"m_VB {m['m_VB']:+.3f} m_CB {m['m_CB']:+.3f} | I_ON n {t['n EOT 1.0']['I_on']:.2f}/{t['n EOT 0.6']['I_on']:.2f} "
                  f"p {t['p EOT 1.0']['I_on']:.2f}/{t['p EOT 0.6']['I_on']:.2f} | PBE gap {row['PBE']['gap']:.3f} "
                  f"| {row['seconds']:.0f} s", flush=True)


def dft():
    """Overlay points: ToB on the calibration ribbons' HSE device Hamiltonians (A4b-2 decisions); DFT gaps."""
    md5 = _frozen()
    res = dict(prereg_md5=md5, calibration={}, holdouts={})
    w2 = json.loads((ROOT / "dft/A4b/a4b2_windows.json").read_text())
    for fam, (_, _, kc, _) in FAM.items():
        st = MB.geometry(kc)
        var = DEVICE_VAR[kc]
        wfmap = json.loads((MB.W90 / kc / "wf_map.json").read_text())["S" if var.endswith("3") else "B"]
        hr = W.read_hr(G.hr_path(MB, kc, var))
        assert hr.num_wann == len(wfmap)
        nv = nval(st)
        Lx, Wd, nat = st.cell[0, 0], S.width(st), n_host(st)
        E = G.w_eigs(hr, KF)
        m = G.ribbon_metrics(E, Lx, nv)
        m["tob"] = tob_all(G.w_eigs(hr, KG), nv, Lx, Wd, nat)
        Ep = G.w_eigs(W.read_hr(G.hr_path(MB, kc, "W2")), KF)
        res["calibration"][fam] = dict(key=kc, variant=var, W_A=Wd, HSE_grid_gap=w2[kc]["gap_grid"], HSE=m,
                                       PBE_W2=G.ribbon_metrics(Ep, Lx, nv))
        print(f"  {kc} ({var}): W {Wd / 10:.2f} nm, gap {m['gap']:.3f} (grid {w2[kc]['gap_grid']:.3f}); I_ON n "
              f"{m['tob']['n EOT 0.6']['I_on']:.2f} p {m['tob']['p EOT 0.6']['I_on']:.2f} (EOT 0.6)", flush=True)
    hse = json.loads((ROOT / "dft/A5/hse/a5_hse.json").read_text())
    bd = json.loads((ROOT / "fits/holdouts/gate_bd.json").read_text())
    for kh, fam in (("r110_8", "110"), ("r100_6", "100")):
        rel = M5.relaxed(kh)
        res["holdouts"][fam] = dict(key=kh, W_A=S.width(rel), HSE_grid_gap=hse[kh]["gap_grid"],
                                    PBE_W2=bd["ribbons"][kh]["wannier"])
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "dft.json").write_text(json.dumps(res, indent=1, default=float))


def sk_overlay():
    """SK-HSE (joint fit, gate (d) recipe: sheet surface shifts on all atoms, HW2 edge blocks of r100_3) on bare [100],
    on the same copy-step geometry."""
    md5 = _frozen()
    corr, _, _ = G.edge_corrections(MB, "r100_3", G.EDGE_SRC["HSE06"]["r100_3"])
    x_all = G.params("HSE06")
    d = OUT / "sk"
    d.mkdir(parents=True, exist_ok=True)
    for n in FAM["100"][3]:
        f = d / f"100_n{n}.json"
        if f.exists():
            continue
        st0, T = target("100", n)
        geo = bb_geometry(st0, T, "r100_3")
        edge = {a: corr[("L" if T.row[a] == 0 else "R", T.species[a], int(T.plane[a]))]
                for a in range(len(T.species)) if T.row[a] in (0, T.nrow - 1)}
        md = G.model(geo, range(len(geo.species)), edge)
        x = md.x_from(x_all)
        nv, Lx, Wd, nat = nval(st0), T.Lx, S.width(geo), n_host(st0)
        m = G.ribbon_metrics(G.sk_eigs(md, x, KF), Lx, nv)
        m["tob"] = tob_all(G.sk_eigs(md, x, KG), nv, Lx, Wd, nat)
        f.write_text(json.dumps(dict(prereg_md5=md5, n=n, W_A=Wd, SK_HSE=m), indent=1, default=float))
        print(f"  SK-HSE [100] n {n}: gap {m['gap']:.3f} I_ON n/p (0.6) {m['tob']['n EOT 0.6']['I_on']:.2f}/"
              f"{m['tob']['p EOT 0.6']['I_on']:.2f}", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "sweep":
        sweep(sys.argv[2:])
    else:
        dict(hcalib=hcalib, geomcheck=geomcheck, dft=dft, sk=sk_overlay)[mode]()
