"""Building-block Wannier ribbons (plan LR2; user decision 2026-11-16): the pre-registered test on the PBE hold-outs.

Recipe and criteria: fits/PREREG_bb_holdouts.md, frozen before `run` reads any hold-out Wannier Hamiltonian.

Modes (run from the project root):
    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/gate_bb_holdouts.py selftest   # sheet + calibration only
    ... calib   # in-sample: building blocks rebuild the calibration ribbons; on-site and hopping transfer
    ... run     # the hold-outs: refuses unless fits/PREREG_bb_holdouts.md matches its frozen md5
Writes fits/bb/{selftest,calib,gate_bb}.json; `run` also figures/Fig3b_holdouts_building_blocks.png.
"""
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np

from scnfet import bb, w90 as W
from scnfet.structures import R45, slab as slab_structure, width

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "fits" / "bb"
PREREG = ROOT / "fits" / "PREREG_bb_holdouts.md"
PREREG_MD5 = ROOT / "fits" / "PREREG_bb_holdouts.md5"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, [str(path)]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


G = _load("gate_bcd", ROOT / "scripts" / "gate_bcd_sk_holdouts.py")   # metrics, ToB, geometry helpers (unchanged)
MB, MBH, M5 = G.MB, G.MBH, G.M5
FUNCS = G.FUNCS
HOLD = G.HOLD                                     # hold-out → calibration ribbon with the same (bare) edge
# The recipe's edge source: W1 / HW1, A2's window rule, the rule the sheet's Wannier functions were made with, so the
# two sources share one Wannier gauge (calib: with W2 the Sc z² elements disagree by 0.2–0.9 eV)
CAL_VAR = {"PBE": {"r100_3": "W1", "r110_4": "W1"}, "HSE06": {"r100_3": "HW1", "r110_4": "HW1"}}
SHEET = {"PBE": ROOT / "dft/A2/w90/pbe_k16/pbe_k16", "HSE06": ROOT / "dft/A2/w90/hse_k16/hse_k16"}
RCUT = 20.0                                       # Å; in-sample: 12 Å moved the r110_4 gap by 21 meV, 20 Å by ≤ 1.2 meV
LIMITS = G.LIMITS
# in-sample variants: the A4b SK-block variants (W2 / HW2 / HW1) and W1 / HW1, A2's window rule, as the sheet ran
CALIB_VARS = {"PBE": {"r100_3": ("W2", "W1"), "r110_4": ("W2", "W1")},
              "HSE06": {"r100_3": ("HW2", "HW1"), "r110_4": ("HW1",)}}
NK_PATH, NK_TOB = G.NK_PATH, G.NK_TOB


# ------------------------------------------------------------------------------------------------- sources
def win_atoms(path):
    """(cell rows, species, positions Å) from a Wannier90 .win (unit_cell_cart in Å; atoms_cart or atoms_frac)."""
    txt = Path(path).read_text()
    cb = re.search(r"begin unit_cell_cart(.*?)end unit_cell_cart", txt, re.S | re.I).group(1).split("\n")
    cell = np.array([[float(v) for v in l.split()] for l in cb if len(l.split()) == 3])
    m = re.search(r"begin atoms_(cart|frac)(.*?)end atoms_\1", txt, re.S | re.I)
    rows = [l.split() for l in m.group(2).split("\n") if len(l.split()) == 4]
    sp = [r[0] for r in rows]
    pos = np.array([[float(v) for v in r[1:]] for r in rows])
    if m.group(1).lower() == "frac":
        pos = pos @ cell
    return cell, sp, pos


def sheet_source(func):
    seed = SHEET[func]
    cell, sp, pos = win_atoms(seed.with_suffix(".win"))
    if func == "PBE":
        ev = json.loads((ROOT / "dft/A2/a2_stage1.json").read_text())["vacuum"]["E_vac_eV"]
    else:
        ev = json.loads((ROOT / "dft/A2/a2_w90.json").read_text())["criterion3_hse_vacuum"]["E_vac_eV"]
    src = bb.Source(hr=W.read_hr(Path(f"{seed}_hr.dat")), cell=cell, pos=pos, species=sp, wf=bb.wf_list(sp), E_vac=ev)
    return bb.Sheet(src, a_p=float(cell[0, 0]), plane_z=float(pos[:, 2].mean()))


def cal_source(key, func, var=None):
    var = var or CAL_VAR[func][key]
    hp = G.hr_path(MB, key, var)
    cell, sp, pos = win_atoms(hp.parent / f"{key}_{var}.win")
    st = MB.geometry(key)
    assert sp == list(st.species), f"{key}: .win species order differs"
    wfmap = json.loads((MB.W90 / key / "wf_map.json").read_text())["B"]
    if func == "PBE":
        w = json.loads((ROOT / "dft/A4b/a4b_windows.json").read_text())[key]
        ev = w["vbm"] + w["E_vac_rel_vbm"]
    else:
        ev = json.loads((ROOT / "dft/A4b/a4b2_windows.json").read_text())[key]["vacuum"]["E_vac_eV"]
    src = bb.Source(hr=W.read_hr(hp), cell=cell, pos=pos, species=sp,
                    wf=[(m["atom"], m["orbital"]) for m in wfmap], E_vac=ev)
    return src, var


def ideal(mk, key):
    st0 = mk.structure(key)
    return bb.Ribbon(pos=st0.pos, species=list(st0.species), Lx=float(st0.cell[0, 0])), st0


def q_sheet(st):
    """Target-frame vector → sheet-frame vector: d_s = F_S F_Tᵀ d (frames hold the cubic axes as columns)."""
    F_t = np.eye(3) if st.frame is None else np.asarray(st.frame, float)
    return R45 @ F_t.T


def s_atoms_of(mb, key, src=None):
    """Atoms with an Sc s WF: the edge Sc of bare [110] (make_a4b EDGE_S)."""
    st = mb.geometry(key)
    left, right = mb.edge_rows(key, st)
    if src is not None:                            # calibration: whatever the WF map has
        return sorted(a for a, o in src.wf if o == "s")
    return [i for i, s in enumerate(st.species) if s == "Sc" and (left[i] or right[i])]


def build(target_mb, target_key, cal_key, func, var=None, align="centre", h=None):
    """Building-block H(R) for a target ribbon. align: "centre" (the recipe: the sheet's on-site energies matched to
    the calibration ribbon's centre rows) or "vacuum" (descriptive: the two vacuum levels aligned)."""
    sheet = sheet_source(func)
    csrc, var = cal_source(cal_key, func, var)
    C, _ = ideal(MB.MK, cal_key)
    T, st0 = ideal(target_mb.MK, target_key)
    s_atoms = s_atoms_of(target_mb, target_key) if any(o == "s" for _, o in csrc.wf) else []
    shift = bb.centre_shift(C, csrc, sheet) if align == "centre" else 0.0
    hr, wf, used = bb.assemble(T, C, csrc, sheet, q_sheet(st0), RCUT, s_atoms=s_atoms, h=h, sheet_shift=shift)
    return hr, wf, used, dict(variant=var, sheet_shift_eV=shift, align=align,
                              n_sheet_pairs=sum(1 for v in used.values() if v == "S"), n_pairs=len(used))


def eig_s(hr, wf, kf):
    s_idx = [j for j, (_, o) in enumerate(wf) if o == "s"]
    if not s_idx:
        return G.w_eigs(hr, kf), None
    E, V = G.w_eigs(hr, kf, vecs=True)
    return E, (np.abs(V[:, s_idx, :]) ** 2).sum(axis=1)


def window_dev(E, Eref, nv, lo=1.0, hi=1.0):
    """max |ΔE| (meV) over bands whose reference lies within [VBM − lo, CBM + hi] (sorted, same band index; both on
    the same energy zero)."""
    vbm, cbm = Eref[:, nv - 1].max(), Eref[:, nv:].min()
    m = (Eref >= vbm - lo) & (Eref <= cbm + hi)
    return float(1e3 * np.abs(E - Eref)[m].max())


# ------------------------------------------------------------------------------------------------- modes
def selftest():
    res = {}
    # 1. rotation + lookup: a [100]-aligned conventional sheet cell assembled from rotated sheet elements reproduces
    #    the primitive sheet at the folded k (Γ_c ← Γ_p ∪ M̄_p; generic k_c ← R45 k_c ∪ R45 (k_c − G_c))
    for func in FUNCS:
        sh = sheet_source(func)
        a_p = sh.a_p
        conv = slab_structure(2, a_p=a_p, cell="conventional")
        Q = q_sheet(conv)
        sp = list(conv.species)
        cls = [(s, int(z > conv.pos[:, 2].mean())) for s, z in zip(sp, conv.pos[:, 2])]
        wf = bb.wf_list(sp)
        orbs = {}
        for j, (a, o) in enumerate(wf):
            orbs.setdefault(a, []).append((j, o))
        a_c = conv.cell[0, 0]
        devs = []
        for kc in (np.zeros(3), np.array([0.31, 0.17, 0.0])):
            Hk = np.zeros((len(wf), len(wf)), complex)
            for a in range(len(sp)):
                ia, oa = zip(*orbs[a])
                Ca = bb.orb_rot(Q, list(oa))
                for b in range(len(sp)):
                    ib, ob = zip(*orbs[b])
                    Cb = bb.orb_rot(Q, list(ob))
                    for i in range(-6, 7):
                        for j in range(-6, 7):
                            d = conv.pos[b] + np.array([i * a_c, j * a_c, 0]) - conv.pos[a]
                            if np.linalg.norm(d[:2]) > RCUT:
                                continue
                            blk = sh.block(cls[a], cls[b], Q @ d, list(oa), list(ob))
                            if blk is not None:
                                Hk[np.ix_(ia, ib)] += np.exp(1j * kc @ d) * (Ca.T @ blk @ Cb)
            e_c = np.linalg.eigvalsh(0.5 * (Hk + Hk.conj().T))
            Gc = np.array([2 * np.pi / a_c, 0, 0])
            e_p = []
            for kk in (kc, kc - Gc):
                kp = Q @ kk
                Hp = np.zeros((16, 16), complex)
                s = sh.src
                for iR, R in enumerate(s.hr.R):
                    for m in range(4):
                        for n in range(4):
                            d = s.pos[n] + R @ s.cell - s.pos[m]
                            if np.linalg.norm(d[:2]) > RCUT:
                                continue
                            im = [s.idx[m][o] for o in bb.ORB[s.species[m]]]
                            jn = [s.idx[n][o] for o in bb.ORB[s.species[n]]]
                            Hp[np.ix_(im, jn)] += np.exp(1j * kp @ d) * s.hr.H[iR][np.ix_(im, jn)]
                e_p.append(np.linalg.eigvalsh(0.5 * (Hp + Hp.conj().T)))
            e_p = np.sort(np.concatenate(e_p)) - sh.src.E_vac
            devs.append(float(np.abs(e_c - e_p).max()))
        res[f"rotation {func}"] = dict(max_dev_eV=max(devs), ok=bool(max(devs) < 1e-8))
    # 2. orb_rot is orthogonal and a 90° rotation maps d_xz → d_yz (sign as the real harmonics)
    Q90 = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1.0]])
    Cd = bb.orb_rot(Q90, list(bb.ORB["Sc"]))
    res["orb_rot"] = dict(orthogonal=bool(np.allclose(Cd.T @ Cd, np.eye(5))),
                          ok=bool(np.allclose(Cd.T @ Cd, np.eye(5)) and np.allclose(bb.orb_rot(np.eye(3), ["px", "py", "pz"]), np.eye(3))))
    # 3. the calibration ribbon rebuilt from itself alone (h = N_C/2, the sheet never consulted for in-half pairs)
    #    reproduces its own Wannier bands where only in-half pairs exist: checked in calib; here the atom maps
    for kh, kc in HOLD.items():
        T, _ = ideal(MBH.MK, kh)
        C, _ = ideal(MB.MK, kc)
        res[f"rows {kh} / {kc}"] = dict(T_rows=int(T.nrow), C_rows=int(C.nrow), h=int(C.nrow // 2),
                                       ok=bool((T.nrow - C.nrow) % 2 == 0))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "selftest.json").write_text(json.dumps(res, indent=1))
    for k, v in res.items():
        print(f"  {k}: {v}")
    assert all(v["ok"] for v in res.values()), "selftest failed"


def calib():
    """In-sample (sheet + calibration ribbons only):
    * range: the largest |H| beyond RCUT in the sheet and the calibration ribbons;
    * on-site transfer: per row, the vacuum-referred mean Sc d and N p on-site energies against the sheet's;
    * hopping transfer: centre-row pairs of the calibration ribbon against the sheet's element for the same bond;
    * rebuild: building blocks for the calibration ribbon itself, with h = N_C/2 (the recipe: only cross-half pairs
      from the sheet) and h = 1 (only the outermost row from the ribbon), against its own Wannier bands."""
    res = {}
    t = np.linspace(0, 0.5, NK_PATH)
    kf = np.c_[t, np.zeros((NK_PATH, 2))]
    for func in FUNCS:
        sh = sheet_source(func)
        s = sh.src
        far = []
        for iR, R in enumerate(s.hr.R):
            for m in range(4):
                for n in range(4):
                    d = s.pos[n] + R @ s.cell - s.pos[m]
                    if np.linalg.norm(d[:2]) > RCUT:
                        far.append(np.abs(s.hr.H[iR][np.ix_([s.idx[m][o] for o in bb.ORB[s.species[m]]],
                                                            [s.idx[n][o] for o in bb.ORB[s.species[n]]])]).max())
        res[f"range sheet {func}"] = dict(max_beyond_rcut_meV=1e3 * float(max(far)))
        for kc, var_ in [(kc, v) for kc in HOLD.values() for v in CALIB_VARS[func][kc]]:
            csrc, var = cal_source(kc, func, var_)
            C, st0 = ideal(MB.MK, kc)
            Q = q_sheet(st0)
            # range
            far = 0.0
            for iR, R in enumerate(csrc.hr.R):
                for a in range(len(C.species)):
                    for b in range(len(C.species)):
                        d = csrc.pos[b] + R @ csrc.cell - csrc.pos[a]
                        if np.linalg.norm(d) > RCUT:
                            ia = list(csrc.idx[a].values())
                            ib = list(csrc.idx[b].values())
                            far = max(far, float(np.abs(csrc.hr.H[iR][np.ix_(ia, ib)]).max()))
            # on-site transfer by row
            cls = [(sp, int(p)) for sp, p in zip(C.species, C.plane)]
            onsite = {}
            for a in range(len(C.species)):
                o = list(bb.ORB[C.species[a]])
                mine = np.diag(csrc.block(a, a, (0, 0, 0), o, o)).mean()
                ref = np.diag(sh.block(cls[a], cls[a], np.zeros(3), o, o)).mean()
                onsite.setdefault(f"row {int(C.row[a])} {C.species[a]}", []).append(1e3 * float(mine - ref))
            onsite = {k: round(float(np.mean(v)), 1) for k, v in sorted(onsite.items(), key=lambda kv: int(kv[0].split()[1]))}
            # hopping transfer: pairs with both atoms in the two centre rows, |d| ≤ 4.5 Å (NN, 2NN, 3NN)
            hop = {}
            cen = [a for a in range(len(C.species)) if C.row[a] in (C.nrow // 2 - 1, C.nrow // 2)]
            for a in cen:
                for b in cen:
                    for Rx in (-1, 0, 1):
                        d = C.pos[b] + np.array([Rx * C.Lx, 0, 0]) - C.pos[a]
                        r = np.linalg.norm(d)
                        if r < 1e-6 or r > 4.5:
                            continue
                        oa, ob = list(bb.ORB[C.species[a]]), list(bb.ORB[C.species[b]])
                        dc = csrc.pos[b] - csrc.pos[a]
                        Rc = int(np.rint((d[0] - dc[0]) / C.Lx))
                        mine = csrc.block(a, b, (Rc, 0, 0), oa, ob)
                        bs = sh.block(cls[a], cls[b], Q @ d, oa, ob)
                        ref = bb.orb_rot(Q, oa).T @ bs @ bb.orb_rot(Q, ob)
                        shell = f"{C.species[a]}-{C.species[b]} {r:.2f}"
                        hop[shell] = max(hop.get(shell, 0.0), 1e3 * float(np.abs(mine - ref).max()))
            # rebuild
            Ew_hr = csrc.hr
            nv = 3 * C.species.count("N")
            wf_c = csrc.wf
            Ew, sw = eig_s(Ew_hr, wf_c, kf)
            Ew = Ew - csrc.E_vac                       # the building blocks are vacuum-referred
            ref_m = G.ribbon_metrics(Ew, C.Lx, nv, sw)
            rows = {}
            for h in sorted({C.nrow // 2, 1}, reverse=True):
                s_atoms = sorted(a for a, o in csrc.wf if o == "s")
                hr_b, wf_b, used = bb.assemble(C, C, csrc, sh, Q, RCUT, s_atoms=s_atoms, h=h)
                Eb, swb = eig_s(hr_b, wf_b, kf)
                m = G.ribbon_metrics(Eb, C.Lx, nv, swb)
                rows[f"h = {h}"] = dict(metrics=m, vs_own=G.compare(m, ref_m), window_dev_meV=window_dev(Eb, Ew, nv),
                                        n_sheet_pairs=sum(1 for v in used.values() if v == "S"),
                                        n_pairs=len(used))
            res[f"{kc} {func} {var}"] = dict(variant=var, E_vac=csrc.E_vac, sheet_E_vac=s.E_vac,
                                       centre_shift_eV=bb.centre_shift(C, csrc, sh),
                                       max_beyond_rcut_meV=1e3 * far, onsite_minus_sheet_meV=onsite,
                                       centre_hop_dev_meV={k: round(v, 1) for k, v in hop.items()},
                                       own=ref_m, rebuild=rows)
            print(f"  {kc} {func} ({var}): beyond rcut {1e3 * far:.2f} meV | on-site − sheet (meV) {onsite}")
            print(f"     centre-row hoppings − sheet (max meV): {res[f'{kc} {func} {var}']['centre_hop_dev_meV']}")
            for lab, r in rows.items():
                d = r["vs_own"]
                print(f"     rebuild {lab}: gap {r['metrics']['gap']:.3f} vs {ref_m['gap']:.3f} ({d['gap_meV']:+.0f} meV) "
                      f"valley {d['valley_meV']:+.0f} m_VB {d['m_VB_rel']:+.1%} m_CB {d['m_CB_rel']:+.1%} "
                      f"window dev {r['window_dev_meV']:.0f} meV ({r['n_sheet_pairs']} of {r['n_pairs']} pairs from the sheet)")
        print(f"  sheet {func}: beyond rcut {res[f'range sheet {func}']['max_beyond_rcut_meV']:.2f} meV")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "calib.json").write_text(json.dumps(res, indent=1, default=float))


def _frozen_md5():
    md5 = hashlib.md5(PREREG.read_bytes()).hexdigest()
    if not PREREG_MD5.exists() or PREREG_MD5.read_text().split()[0] != md5:
        sys.exit("fits/PREREG_bb_holdouts.md is not frozen (md5 missing or changed): refusing to run the hold-outs")
    return md5


def bands_and_metrics(hr, wf, Lx, nv, kf):
    E, sw = eig_s(hr, wf, kf)
    return E, G.ribbon_metrics(E, Lx, nv, sw)


def tob_check(hr, ref, nv, Lx, Wd, nat, kg):
    I = G.ion(G.w_eigs(hr, kg), nv, Lx, Wd, nat)
    rel = {k: I[k]["I_on"] / ref["I_on"][k]["I_on"] - 1 for k in I}
    verdict = {c: (I[f"{c} EOT 0.6"]["I_on"] >= G.IRDS_ION, ref["I_on"][f"{c} EOT 0.6"]["I_on"] >= G.IRDS_ION)
               for c in ("n", "p")}
    return I, rel, verdict


def run():
    md5 = _frozen_md5()
    res = dict(prereg_md5=md5, limits=LIMITS, rcut_A=RCUT, recipe_variants=CAL_VAR, ribbons={})
    hse = json.loads((ROOT / "dft" / "A5" / "hse" / "a5_hse.json").read_text())
    kf = np.c_[np.linspace(0, 0.5, NK_PATH), np.zeros((NK_PATH, 2))]
    kg = np.c_[np.arange(NK_TOB) / NK_TOB, np.zeros((NK_TOB, 2))]
    plot = {}
    for kh, kc in HOLD.items():
        st = MBH.geometry(kh)
        nv = 3 * st.species.count("N")
        assert nv == MBH.nocc(st) - MBH.nexcl(st)
        Lx, Wd, nat = st.cell[0, 0], width(st), len(st.species)
        # reference: the PBE-Wannier hold-out W2 (A5-2 decision), exactly as gate (d) read it
        wfmap = json.loads((M5.TWIN / "w90" / kh / "wf_map.json").read_text())["B"]
        hr_w = W.read_hr(M5.TWIN / "w90" / kh / "W2" / f"{kh}_W2_hr.dat")
        Ew, ref = bands_and_metrics(hr_w, [(m["atom"], m["orbital"]) for m in wfmap], Lx, nv, kf)
        ref["I_on"] = G.ion(G.w_eigs(hr_w, kg), nv, Lx, Wd, nat)
        r = res["ribbons"][kh] = dict(calibration=kc, W_A=Wd, reference="PBE-Wannier W2", wannier=ref)
        n_s_ref = sum(1 for m in wfmap if m["orbital"] == "s")
        for lab, func, var, align in (("BB PBE", "PBE", None, "centre"),                 # registered
                                      ("BB HSE06", "HSE06", None, "centre"),             # registered (gap)
                                      ("BB PBE, vacuum-aligned", "PBE", None, "vacuum"), # descriptive
                                      ("BB PBE, W2 edges", "PBE", "W2", "centre")):      # descriptive
            hr_b, wf_b, used, info = build(MBH, kh, kc, func, var=var, align=align)
            assert sum(1 for _, o in wf_b if o == "s") == n_s_ref, "edge-s count differs from the reference"
            E, m = bands_and_metrics(hr_b, wf_b, Lx, nv, kf)
            row = dict(build=info, metrics=m)
            if func == "PBE":
                d = G.compare(m, ref)
                I, rel, verdict = tob_check(hr_b, ref, nv, Lx, Wd, nat, kg)
                d["ok"]["I_on"] = bool(max(abs(v) for v in rel.values()) <= LIMITS["ion_rel"])
                d["ok"]["irds"] = bool(all(a == b for a, b in verdict.values()))
                d["pass"] = bool(all(d["ok"].values()))
                row.update(vs_W2=d, I_on=I, I_on_rel=rel, irds_verdict_bb_vs_w=verdict)
                if lab == "BB PBE":
                    plot[kh] = (Ew - Ew[:, nv - 1].max(), E - E[:, nv - 1].max(), nv)
            else:
                g = hse[kh]["gap_grid"]
                row["vs_HSE_gap"] = dict(hse_gap=g, d_meV=1e3 * (m["gap"] - g),
                                         ok=bool(abs(1e3 * (m["gap"] - g)) <= LIMITS["hse_gap_meV"]))
            r[lab] = row
        r["pass"] = bool(r["BB PBE"]["vs_W2"]["pass"] and r["BB HSE06"]["vs_HSE_gap"]["ok"])
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "gate_bb.json").write_text(json.dumps(res, indent=1, default=float))
    report(res)
    figures(plot, res)


def report(res):
    for kh, r in res["ribbons"].items():
        w = r["wannier"]
        print(f"=== {kh} (edges from {r['calibration']}; W {r['W_A'] / 10:.2f} nm): W2 gap {w['gap']:.3f}, valley X−Γ "
              f"{w['valley_XG']:.3f}, m_VB {w['m_VB']:+.3f}, m_CB {w['m_CB']:+.3f}"
              + (f" | W2 without edge-s bands: gap {w['no_edge_s']['gap']:.3f}, valley {w['no_edge_s']['valley_XG']:.3f}"
                 if "no_edge_s" in w else ""))
        for lab in ("BB PBE", "BB PBE, vacuum-aligned", "BB PBE, W2 edges"):
            m, d, b = r[lab]["metrics"], r[lab]["vs_W2"], r[lab]["build"]
            print(f"   {lab} ({b['variant']}, shift {b['sheet_shift_eV']:+.3f} eV): gap {m['gap']:.3f} ({d['gap_meV']:+.0f} meV) "
                  f"valley {m['valley_XG']:.3f} ({d['valley_meV']:+.0f}) m_VB {m['m_VB']:+.3f} ({d['m_VB_rel']:+.1%}) "
                  f"m_CB {m['m_CB']:+.3f} ({d['m_CB_rel']:+.1%}) | I_ON rel "
                  + " ".join(f"{k} {v:+.1%}" for k, v in r[lab]["I_on_rel"].items())
                  + f" | IRDS (BB, W) {r[lab]['irds_verdict_bb_vs_w']} | ok {d['ok']} → {'PASS' if d['pass'] else 'FAIL'}"
                  + (f" | no edge-s: gap {m['no_edge_s']['gap']:.3f} valley {m['no_edge_s']['valley_XG']:.3f}"
                     if "no_edge_s" in m else ""))
        h = r["BB HSE06"]
        print(f"   BB HSE06 ({h['build']['variant']}, shift {h['build']['sheet_shift_eV']:+.3f} eV): gap "
              f"{h['metrics']['gap']:.3f} vs HSE {h['vs_HSE_gap']['hse_gap']:.3f} ({h['vs_HSE_gap']['d_meV']:+.0f} meV) "
              f"{'ok' if h['vs_HSE_gap']['ok'] else 'FAIL'}")
        print("   W2 I_ON (mA/μm): " + " ".join(f"{k} {v['I_on']:.3f}" for k, v in w["I_on"].items())
              + f" | hold-out {'PASSES' if r['pass'] else 'FAILS'}")


def figures(plot, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    REF, BBC = "#3b3b3b", "#c2410c"
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.4))
    t = np.linspace(0, 0.5, NK_PATH)
    for ax, kh in zip(axes, ("r110_8", "r100_6")):
        Ew, E, nv = plot[kh]
        ax.plot(t, Ew, color=REF, lw=1.4)
        ax.plot(t, E, color=BBC, lw=1.0, ls="--")
        g = res["ribbons"][kh]["wannier"]["gap"]
        ax.set_ylim(-1.5, g + 1.5)
        ax.set_xlim(0, 0.5)
        ax.set_xticks([0, 0.5], ["Γ", "X"])
        ax.set_title({"r110_8": "[110] n = 8", "r100_6": "[100] n = 6"}[kh], fontsize=10)
        ax.set_ylabel("E − E_VBM (eV)")
        ax.axhline(0, color="#bbbbbb", lw=0.6, zorder=0)
    h = [plt.Line2D([], [], color=REF, lw=1.4), plt.Line2D([], [], color=BBC, lw=1.0, ls="--")]
    fig.legend(h, ["PBE-Wannier hold-out (W2)", "building blocks: calibration edges (W1) + sheet interior"],
               loc="lower center", ncol=2, fontsize=8.5, frameon=False)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(ROOT / "figures" / "Fig3b_holdouts_building_blocks.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    dict(selftest=selftest, calib=calib, run=run)[sys.argv[1]]()
