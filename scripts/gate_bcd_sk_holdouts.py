"""Claims gates (b), (c) and (d) on the laptop: the joint-fit SK models (fits/A2_joint) on the 3-plane slab and the
ribbons. Recipe and criteria: fits/PREREG_sk_holdouts.md, frozen before `run` computed any SK hold-out number.

Modes (run from the project root):
    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/gate_bcd_sk_holdouts.py selftest   # no hold-out data used
    ... calib   # gate (c), in-sample: SK + edge blocks on the calibration ribbons the blocks come from
    ... run     # gates (b) and (d): refuses unless fits/PREREG_sk_holdouts.md matches its frozen md5
Writes fits/holdouts/{selftest,gate_c,gate_bd}.json; `run` also figures/Fig2c_slab_holdout.png and
figures/Fig3_holdouts_PBE_twin.png.
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from scnfet import sk, tob, w90 as W
from scnfet.structures import R45, slab as slab_structure, width

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "fits" / "holdouts"
PREREG = ROOT / "fits" / "PREREG_sk_holdouts.md"
PREREG_MD5 = ROOT / "fits" / "PREREG_sk_holdouts.md5"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, [str(path)]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


A2S = _load("a2_sk_fit_bulk_sheet", ROOT / "scripts" / "a2_sk_fit_bulk_sheet.py")   # SHIFTS, A_BULK
PB = _load("parse_a4b", ROOT / "dft" / "A4b" / "parse_a4b.py")                      # edge_mass, HB2M
MB = PB.MB                                                                          # make_a4b: calibration ribbons
M5 = _load("make_a5", ROOT / "dft" / "A5" / "make_a5.py")
MBH = M5._make_a4b(M5.TWIN)                                                         # make_a4b on the hold-outs
P3 = _load("parse_a3", ROOT / "dft" / "A3" / "parse_a3.py")                         # frac_path, M3.relaxed

A_BULK = A2S.A_BULK
NAMES = sk.PARAM_NAMES + sk.T3_PARAM_NAMES
FUNCS = ("PBE", "HSE06")
HOLD = {"r110_8": "r110_4", "r100_6": "r100_3"}          # hold-out → calibration ribbon with the same (bare) edge
# Edge-block sources: the variants the A4b decisions named for the SK blocks (A4b-1: W2; A4b-2: [100] W2, [110] W1)
EDGE_SRC = {"PBE": {"r100_3": "W2", "r110_4": "W2"}, "HSE06": {"r100_3": "HW2", "r110_4": "HW1"}}
NK_PATH, NK_TOB, CHUNK = 601, 2000, 200
LIMITS = dict(gap_meV=50.0, valley_meV=50.0, mass_rel=0.10, ion_rel=0.15, hse_gap_meV=80.0)
TOB_SPAN, TOB_NSWEEP = 4.0, 381                          # gate sweep ±4 V at the default ≈ 21 mV step (q·n_SD/C_G ≈ 2 V)
IRDS_ION = 1.979                                          # mA/μm at EOT 0.6 nm (notebook, device conventions)
SK_ORDER = {"Sc": ("dxy", "dyz", "dxz", "dx2-y2", "dz2"), "N": ("px", "py", "pz")}   # sk.ORBITALS, W90 names


def params(func):
    return json.loads((ROOT / "fits" / "A2_joint" / f"sk_joint_{func}.json").read_text())["params"]


# ------------------------------------------------------------------------------------------ the SK structures
DS = sk.d_rep(R45.T)                                      # sheet frame (A2: primitive cell, frame R45)


def shift_terms(st, sites):
    """The sheet's fitted C4v on-site shifts (A2S.SHIFTS, given in the sheet-cell frame) on `sites`, rotated into the
    structure's frame: S_t = T S Tᵀ with T = D_t D_sᵀ, D = sk.d_rep(frameᵀ) as in sk.crystal_field_projectors. The N
    patterns (p_z; p_x + p_y) are invariant under rotations about z, which is all that separates the frames here."""
    F = np.eye(3) if st.frame is None else np.asarray(st.frame, float)
    assert np.allclose(F[:, 2], [0, 0, 1]) and np.allclose(R45[:, 2], [0, 0, 1]), "frames must differ by a z rotation"
    T = (np.eye(5) if st.frame is None else sk.d_rep(F.T)) @ DS.T
    terms = {}
    for name, (sp, M) in A2S.SHIFTS.items():
        Mt = T @ M @ T.T if sp == "Sc" else M
        terms[name] = {i: Mt for i in sites if st.species[i] == sp}
    return terms


def model(st, surface_sites, edge=None):
    return sk.build(st, shells=sk.rocksalt_shells(A_BULK), names=NAMES, onsite_terms=shift_terms(st, surface_sites),
                    onsite_const=edge or {}, three_centre=sk.three_centre_rocksalt(A_BULK))


def kcart_of(st, kfrac):
    return 2 * np.pi * np.atleast_2d(kfrac) @ np.linalg.inv(st.cell).T


def sk_eigs(md, x, kfrac):
    kc = kcart_of(md.structure, kfrac)
    return np.vstack([md.eigvals(x, kc[i:i + CHUNK]) for i in range(0, len(kc), CHUNK)])


def w_eigs(hr, kfrac, vecs=False):
    out_e, out_v = [], []
    for i in range(0, len(kfrac), CHUNK):
        Hk = W.hk_frac(hr, kfrac[i:i + CHUNK])
        if vecs:
            e, v = np.linalg.eigh(Hk)
            out_e.append(e)
            out_v.append(v)
        else:
            out_e.append(np.linalg.eigvalsh(Hk))
    return (np.vstack(out_e), np.concatenate(out_v)) if vecs else np.vstack(out_e)


def ideal_rows(mb, key):
    """Row (by ideal-cut y) and plane (by ideal-cut z) of every atom; atom order is that of the relaxed geometry."""
    st0 = mb.MK.structure(key)
    y, z = np.round(st0.pos[:, 1], 3), np.round(st0.pos[:, 2], 3)
    uy, uz = np.unique(y), np.unique(z)
    assert len(uz) == 2, f"{key}: {len(uz)} planes"
    return np.searchsorted(uy, y), len(uy), np.searchsorted(uz, z)


def onsite_blocks(hr_path, wfmap):
    """R = 0 on-site blocks of every host atom in the SK orbital order (W90 real harmonics, signs +1 as A1 showed)."""
    hr = W.read_hr(hr_path)
    assert hr.num_wann == len(wfmap), f"{hr_path}: {hr.num_wann} WFs, map {len(wfmap)}"
    i0 = int(np.where((hr.R == 0).all(axis=1))[0][0])
    H0 = hr.H[i0] / hr.deg[i0]
    out = {}
    for a in sorted({m["atom"] for m in wfmap if m["species"] in SK_ORDER}):
        sp = next(m["species"] for m in wfmap if m["atom"] == a)
        idx = [next(j for j, m in enumerate(wfmap) if m["atom"] == a and m["orbital"] == o) for o in SK_ORDER[sp]]
        B = H0[np.ix_(idx, idx)]
        assert np.abs(B.imag).max() < 1e-5
        out[a] = B.real
    return out


def hr_path(mb, key, var):
    sub = mb.W90 / key / ("H" if var.startswith("H") else "") / var.lstrip("H")
    return sub / f"{key}_{var}_hr.dat"


def edge_corrections(mb, key, var):
    """{(side, species, plane): B} with B = H_W(edge atom) − mean H_W(same species, the two centre rows)."""
    st = mb.geometry(key)
    left, right = mb.edge_rows(key, st)
    row, nrow, plane = ideal_rows(mb, key)
    wfmap = json.loads((mb.W90 / key / "wf_map.json").read_text())["B"]
    blk = onsite_blocks(hr_path(mb, key, var), wfmap)
    centre = [a for a in blk if row[a] in (nrow // 2 - 1, nrow // 2)]
    ref = {sp: np.mean([blk[a] for a in centre if st.species[a] == sp], axis=0) for sp in SK_ORDER}
    corr, atoms = {}, {}
    for a in blk:
        if left[a] or right[a]:
            k = ("L" if left[a] else "R", st.species[a], int(plane[a]))
            assert k not in corr, f"{key}: two edge atoms keyed {k}"
            corr[k], atoms[k] = blk[a] - ref[st.species[a]], a
    return corr, atoms, dict(centre_atoms=centre, n_rows=nrow)


def neighbours(st, a, rcut=3.4):
    T = np.array([[i, 0, 0] for i in (-1, 0, 1)]) @ st.cell
    out = []
    for j, s in enumerate(st.species):
        for t in T:
            d = st.pos[j] + t - st.pos[a]
            if 1e-6 < np.linalg.norm(d) < rcut:
                out.append((s, d))
    return out


def env_match(st1, a1, st2, a2):
    """Max distance (Å) between matched same-species neighbour vectors of two atoms (inf if the counts differ)."""
    n1, n2 = neighbours(st1, a1), neighbours(st2, a2)
    if sorted(s for s, _ in n1) != sorted(s for s, _ in n2):
        return float("inf")
    return float(max(min(np.linalg.norm(d - e) for t, e in n2 if t == s) for s, d in n1))


def ribbon_edge_blocks(mb_t, key_t, corr, atoms_c=None, st_c=None):
    """Edge blocks on the target ribbon's edge atoms, matched to the calibration edge atom by (side, species, plane)."""
    st = mb_t.geometry(key_t)
    left, right = mb_t.edge_rows(key_t, st)
    _, _, plane = ideal_rows(mb_t, key_t)
    B, env = {}, {}
    for a in range(len(st.species)):
        if left[a] or right[a]:
            k = ("L" if left[a] else "R", st.species[a], int(plane[a]))
            B[a] = corr[k]
            if atoms_c is not None:
                env[a] = env_match(st, a, st_c, atoms_c[k])
    return B, env


# ------------------------------------------------------------------------------------------------- metrics
def ribbon_metrics(E, Lx, nv, s_weight=None):
    """E on the 601-point path Γ → X (k_frac 0 … ½). Gap, valley offset E_c(X) − E_c(Γ), band-edge masses (edge_mass)."""
    k = np.linspace(0, 0.5, NK_PATH) * 2 * np.pi / Lx
    iv, ic = nv - 1, nv
    m = dict(gap=float(E[:, ic].min() - E[:, iv].max()), direct=bool(np.argmin(E[:, ic]) == np.argmax(E[:, iv])),
             k_cbm=float(np.argmin(E[:, ic]) / (NK_PATH - 1) / 2), k_vbm=float(np.argmax(E[:, iv]) / (NK_PATH - 1) / 2),
             valley_XG=float(E[-1, ic] - E[0, ic]), vb_GX=float(E[0, iv] - E[-1, iv]))
    for lab, b, kind in (("VB", iv, "v"), ("CB", ic, "c")):
        mm, i0, n, _ = PB.edge_mass(k, E, b, kind)
        m[f"m_{lab}"] = float(mm)
    if s_weight is not None:                         # descriptive: lowest CB with edge-s weight < 0.5 (Wannier only)
        # (no mass: the selected band can switch where the edge s hybridises, so its curvature is not meaningful)
        Ec = np.array([E[i, nv:][s_weight[i, nv:] < 0.5].min() for i in range(len(E))])
        m["no_edge_s"] = dict(gap=float(Ec.min() - E[:, iv].max()), valley_XG=float(Ec[-1] - Ec[0]))
    return m


def compare(sk_m, ref_m):
    d = dict(gap_meV=1e3 * (sk_m["gap"] - ref_m["gap"]), valley_meV=1e3 * (sk_m["valley_XG"] - ref_m["valley_XG"]),
             m_VB_rel=sk_m["m_VB"] / ref_m["m_VB"] - 1, m_CB_rel=sk_m["m_CB"] / ref_m["m_CB"] - 1)
    d["ok"] = dict(gap=abs(d["gap_meV"]) <= LIMITS["gap_meV"], valley=abs(d["valley_meV"]) <= LIMITS["valley_meV"],
                   masses=max(abs(d["m_VB_rel"]), abs(d["m_CB_rel"])) <= LIMITS["mass_rel"])
    return d


def ion(E_grid, nv, Lx, Wd, natoms):
    ch = tob.Channel(E=np.sort(E_grid, axis=1), a=Lx * 1e-10, n_val=nv)
    res = {}
    for carrier in ("n", "p"):
        for eot in (1.0, 0.6):
            f = tob.ToBFET(ch, W=Wd * 1e-10, EOT=eot * 1e-9, carrier=carrier, n_atoms_per_cell=natoms)
            m = f.metrics(span=TOB_SPAN, n_sweep=TOB_NSWEEP)
            res[f"{carrier} EOT {eot}"] = dict(I_on=float(m["I_on"]), SS=float(m["SS_mV_dec"]), v_inj=float(m["v_inj"]))
    return res


def slab_metrics(E, st, nv):
    """Γ–X̄–M̄–Γ, 200 points per segment: gap, E_c(M̄) − E_c(Γ), E_c(X̄) − E_c(Γ), VB/CB (VB−1) masses along Γ–X̄, Γ–M̄."""
    t = P3.frac_path()
    kc = kcart_of(st, t)
    iv, ic = nv - 1, nv
    m = dict(gap=float(E[:, ic].min() - E[:, iv].max()), McG=float(E[400, ic] - E[0, ic]),
             XcG=float(E[200, ic] - E[0, ic]))
    for sn, ix in {"G-X": np.arange(0, 201), "G-M": np.arange(600, 399, -1)}.items():
        s = np.linalg.norm(kc[ix] - kc[ix[0]], axis=1)
        for lab, b, kind in (("VB", iv, "v"), ("CB", ic, "c"), ("VB-1", iv - 1, "v")):
            mm, *_ = PB.edge_mass(s, E[ix], b, kind)
            m[f"m_{lab} {sn}"] = float(mm)
    return m


def slab_sk(func):
    st = P3.M3.relaxed()
    z = st.pos[:, 2]
    zs = np.sort(z)
    plane = np.searchsorted(zs[np.r_[False, np.diff(zs) > 1.0]], z, side="right")     # planes by z gaps > 1 Å
    assert sorted(set(plane)) == [0, 1, 2], plane
    surf = [i for i in range(len(z)) if plane[i] in (0, 2)]
    md = model(st, surf)
    return st, md, md.x_from(params(func)), surf


# ------------------------------------------------------------------------------------------------- modes
def selftest():
    res = {}
    # 1. the shift rotation: a [100]-aligned sheet supercell (frame None) reproduces the primitive sheet (frame R45)
    #    at the folded k: Γ_c ← Γ_p ∪ M̄_p, and a generic k_c ← R45 k_c ∪ R45 (k_c − G_c)
    prim, conv = slab_structure(2, cell="primitive"), slab_structure(2, cell="conventional")
    for func in FUNCS:
        x = params(func)
        mp, mc = model(prim, range(len(prim.species))), model(conv, range(len(conv.species)))
        xp, xc = mp.x_from(x), mc.x_from(x)
        a_c = conv.cell[0, 0]
        dev = []
        for kc in (np.zeros(3), np.array([0.31, 0.17, 0.0])):
            Gc = np.array([2 * np.pi / a_c, 0, 0])
            e_c = mc.eigvals(xc, kc[None])[0]
            e_p = np.sort(np.concatenate([mp.eigvals(xp, (R45 @ kk)[None])[0] for kk in (kc, kc - Gc)]))
            dev.append(float(np.abs(e_c - e_p).max()))
        # the same check without rotating the shifts must fail (so the test can see a wrong frame)
        bad = sk.build(conv, shells=sk.rocksalt_shells(A_BULK), names=NAMES, three_centre=sk.three_centre_rocksalt(A_BULK),
                       onsite_terms={n: {i: M for i in range(len(conv.species)) if conv.species[i] == sp}
                                     for n, (sp, M) in A2S.SHIFTS.items()})
        e_bad = bad.eigvals(bad.x_from(x), np.zeros((1, 3)))[0]
        e_p0 = np.sort(np.concatenate([mp.eigvals(xp, k[None])[0] for k in (np.zeros(3), np.array([np.pi / prim.cell[0, 0]] * 2 + [0]))]))
        res[f"rotation {func}"] = dict(max_dev_eV=max(dev), unrotated_dev_eV=float(np.abs(e_bad - e_p0).max()),
                                       ok=bool(max(dev) < 1e-8 and np.abs(e_bad - e_p0).max() > 1e-3))
    # 2. onsite_const lands on the right block
    st = MB.geometry("r100_3")
    a = st.species.index("Sc")
    M = np.arange(25.0).reshape(5, 5) / 100
    M = M + M.T
    m1, m0 = model(st, range(len(st.species)), {a: M}), model(st, range(len(st.species)))
    o = np.where(m1.orb_site == a)[0]
    D = m1.C[0] - m0.C[0]
    rest = np.ones(D.shape, bool)
    rest[np.ix_(o, o)] = False
    res["onsite_const"] = dict(site=a, ok=bool(np.allclose(D[np.ix_(o, o)], M) and np.allclose(D[rest], 0)))
    # 3. edge atoms of each hold-out map to a calibration edge atom with the same local environment (geometry only)
    for kh, kc in HOLD.items():
        corr, atoms, _ = edge_corrections(MB, kc, EDGE_SRC["PBE"][kc])
        _, env = ribbon_edge_blocks(MBH, kh, corr, atoms, MB.geometry(kc))
        res[f"edge map {kh} ← {kc}"] = dict(n_edge_atoms=len(env), max_env_dev_A=max(env.values()),
                                            ok=bool(max(env.values()) < 0.25))
    # 4. the slab's surface planes
    st, md, x, surf = slab_sk("PBE")
    res["slab surface sites"] = dict(n=len(surf), of=len(st.species), ok=len(surf) == 4)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "selftest.json").write_text(json.dumps(res, indent=1))
    for k, v in res.items():
        print(f"  {k}: {v}")
    assert all(v["ok"] for v in res.values()), "selftest failed"


def ribbon_sk(mb, key, func, corr):
    st = mb.geometry(key)
    B = ribbon_edge_blocks(mb, key, corr)[0] if corr is not None else {}
    md = model(st, range(len(st.species)), B)
    return st, md, md.x_from(params(func))


def calib():
    """Gate (c), in-sample and descriptive: SK (with / without its own edge blocks) vs the Wannier Hamiltonian the
    blocks were read from, on the calibration ribbons."""
    res = {}
    t = np.linspace(0, 0.5, NK_PATH)
    for kc in HOLD.values():
        st = MB.geometry(kc)
        nv = 3 * st.species.count("N")
        assert nv == MB.nocc(st) - MB.nexcl(st)
        Lx = st.cell[0, 0]
        kf = np.c_[t, np.zeros((NK_PATH, 2))]
        for func in FUNCS:
            var = EDGE_SRC[func][kc]
            corr, _, info = edge_corrections(MB, kc, var)
            hr = W.read_hr(hr_path(MB, kc, var))
            wfmap = json.loads((MB.W90 / kc / "wf_map.json").read_text())["B"]
            s_idx = [j for j, m in enumerate(wfmap) if m["orbital"] == "s"]
            if s_idx:
                Ew, Vw = w_eigs(hr, kf, vecs=True)
                sw = (np.abs(Vw[:, s_idx, :]) ** 2).sum(axis=1)
                del Vw
            else:
                Ew, sw = w_eigs(hr, kf), None
            ref = ribbon_metrics(Ew, Lx, nv, sw)
            row = dict(ref=var, wannier=ref, edge_blocks={f"{k[0]} {k[1]} plane {k[2]}": dict(
                diag=np.round(np.diag(B), 3).tolist(), max_offdiag=float(np.abs(B - np.diag(np.diag(B))).max()))
                for k, B in corr.items()}, centre_atoms=info["centre_atoms"])
            for lab, c in (("SK + edge blocks", corr), ("SK, no edge blocks", None)):
                _, md, x = ribbon_sk(MB, kc, func, c)
                m = ribbon_metrics(sk_eigs(md, x, kf), Lx, nv)
                row[lab] = dict(metrics=m, vs_wannier=compare(m, ref))
            res[f"{kc} {func}"] = row
            a, b = row["SK + edge blocks"], row["SK, no edge blocks"]
            print(f"  {kc} {func} ({var}): gap W {ref['gap']:.3f} | SK+B {a['metrics']['gap']:.3f} "
                  f"({a['vs_wannier']['gap_meV']:+.0f} meV) | SK {b['metrics']['gap']:.3f} ({b['vs_wannier']['gap_meV']:+.0f}); "
                  f"valley W {ref['valley_XG']:.3f} SK+B {a['metrics']['valley_XG']:.3f}; m_CB W {ref['m_CB']:+.3f} "
                  f"SK+B {a['metrics']['m_CB']:+.3f}; m_VB W {ref['m_VB']:+.3f} SK+B {a['metrics']['m_VB']:+.3f}"
                  + (f" | W without edge-s bands: gap {ref['no_edge_s']['gap']:.3f} valley {ref['no_edge_s']['valley_XG']:.3f}"
                     if "no_edge_s" in ref else ""))
    (OUT / "gate_c.json").write_text(json.dumps(res, indent=1, default=float))


def _frozen_md5():
    md5 = hashlib.md5(PREREG.read_bytes()).hexdigest()
    if not PREREG_MD5.exists() or PREREG_MD5.read_text().split()[0] != md5:
        sys.exit("fits/PREREG_sk_holdouts.md is not frozen (md5 missing or changed): refusing to run the hold-outs")
    return md5


def run():
    md5 = _frozen_md5()
    res = dict(prereg_md5=md5, limits=LIMITS, tob_sweep=dict(span_V=TOB_SPAN, n=TOB_NSWEEP), ribbons={}, slab={})
    hse = json.loads((ROOT / "dft" / "A5" / "hse" / "a5_hse.json").read_text())
    t = np.linspace(0, 0.5, NK_PATH)
    kf = np.c_[t, np.zeros((NK_PATH, 2))]
    kg = np.c_[np.arange(NK_TOB) / NK_TOB, np.zeros((NK_TOB, 2))]
    plot = {}
    for kh, kc in HOLD.items():
        st = MBH.geometry(kh)
        nv = 3 * st.species.count("N")
        assert nv == MBH.nocc(st) - MBH.nexcl(st)
        Lx, Wd, nat = st.cell[0, 0], width(st), len(st.species)
        wfmap = json.loads((M5.TWIN / "w90" / kh / "wf_map.json").read_text())["B"]
        hr = W.read_hr(M5.TWIN / "w90" / kh / "W2" / f"{kh}_W2_hr.dat")
        s_idx = [j for j, m in enumerate(wfmap) if m["orbital"] == "s"]
        if s_idx:
            Ew, Vw = w_eigs(hr, kf, vecs=True)
            sw = (np.abs(Vw[:, s_idx, :]) ** 2).sum(axis=1)
            del Vw
        else:
            Ew, sw = w_eigs(hr, kf), None
        ref = ribbon_metrics(Ew, Lx, nv, sw)
        r = res["ribbons"][kh] = dict(calibration=kc, W_A=Wd, reference="PBE-Wannier W2", wannier=ref)
        r["wannier"]["I_on"] = ion(w_eigs(hr, kg), nv, Lx, Wd, nat)
        out = {}
        for func in FUNCS:
            corr, atoms, _ = edge_corrections(MB, kc, EDGE_SRC[func][kc])
            for lab, c in (("SK", corr), ("SK, no edge blocks", None)):
                _, md, x = ribbon_sk(MBH, kh, func, c)
                E = sk_eigs(md, x, kf)
                m = ribbon_metrics(E, Lx, nv)
                row = dict(metrics=m)
                if func == "PBE":
                    row["vs_W2"] = compare(m, ref)
                    row["I_on"] = ion(sk_eigs(md, x, kg), nv, Lx, Wd, nat)
                    rel = {k: row["I_on"][k]["I_on"] / ref["I_on"][k]["I_on"] - 1 for k in row["I_on"]}
                    verdict = {c_: (row["I_on"][f"{c_} EOT 0.6"]["I_on"] >= IRDS_ION,
                                    ref["I_on"][f"{c_} EOT 0.6"]["I_on"] >= IRDS_ION) for c_ in ("n", "p")}
                    row["I_on_rel"] = rel
                    row["irds_verdict_sk_vs_w"] = verdict
                    row["vs_W2"]["ok"]["I_on"] = bool(max(abs(v) for v in rel.values()) <= LIMITS["ion_rel"])
                    row["vs_W2"]["ok"]["irds"] = bool(all(a == b for a, b in verdict.values()))
                    row["vs_W2"]["pass"] = bool(all(row["vs_W2"]["ok"].values()))
                    if lab == "SK":
                        plot[kh] = (Ew - Ew[:, nv - 1].max(), E - E[:, nv - 1].max(), nv, Lx)
                else:
                    g = hse[kh]["gap_grid"]
                    row["vs_HSE_gap"] = dict(hse_gap=g, d_meV=1e3 * (m["gap"] - g),
                                             ok=bool(abs(1e3 * (m["gap"] - g)) <= LIMITS["hse_gap_meV"]))
                out[f"{lab} {func}"] = row
        # descriptive: edge-block transfer (D4 drift check), the hold-out's own W2 on-site blocks vs the calibration's
        own, own_atoms, _ = edge_corrections(MBH, kh, "W2")
        cal, _, _ = edge_corrections(MB, kc, EDGE_SRC["PBE"][kc])
        r["edge_block_transfer_PBE"] = {f"{k[0]} {k[1]} plane {k[2]}": dict(
            max_abs_diff_eV=float(np.abs(own[k] - cal[k]).max()),
            diag_own=np.round(np.diag(own[k]), 3).tolist(), diag_cal=np.round(np.diag(cal[k]), 3).tolist()) for k in own}
        r.update(out)
        # ToB k-sampling check (conventions: ≥ 2000 k, convergence-checked): W2 at 4000 k
        kg2 = np.c_[np.arange(2 * NK_TOB) / (2 * NK_TOB), np.zeros((2 * NK_TOB, 2))]
        i2 = ion(w_eigs(hr, kg2), nv, Lx, Wd, nat)
        r["tob_k_check"] = {k: i2[k]["I_on"] / ref["I_on"][k]["I_on"] - 1 for k in i2}
    # slab, gate (b)
    st, _, _, surf = slab_sk("PBE")
    nv = 3 * st.species.count("N")
    t3 = P3.frac_path()
    Ew = w_eigs(W.read_hr(ROOT / "dft/A3/w90/slab3/W2/slab3_W2_hr.dat"), t3)
    Ehw = w_eigs(W.read_hr(ROOT / "dft/A3/w90/slab3/H/W2/slab3_HW2_hr.dat"), t3)
    ref = slab_metrics(Ew, st, nv)
    href = slab_metrics(Ehw, st, nv)
    hgap = json.loads((ROOT / "dft/A3/a3_hwindows.json").read_text())["gap_grid"]
    dft = json.loads((ROOT / "dft/A3/a3_criteria.json").read_text())["W2"]["valleys"]
    res["slab"] = dict(surface_sites=surf, wannier_W2=ref, wannier_HW2=href, hse_grid_gap=hgap,
                       dft_valleys=dft)
    for func in FUNCS:
        _, md, x, _ = slab_sk(func)
        E = sk_eigs(md, x, t3)
        m = slab_metrics(E, st, nv)
        row = dict(metrics=m)
        if func == "PBE":
            reg = ("m_VB G-X", "m_CB G-X", "m_VB G-M", "m_CB G-M")
            d = dict(gap_meV=1e3 * (m["gap"] - ref["gap"]), McG_meV=1e3 * (m["McG"] - ref["McG"]),
                     XcG_meV=1e3 * (m["XcG"] - ref["XcG"]), XcG_vs_dft_meV=1e3 * (m["XcG"] - dft["Ec(X)-Ec(G)"]["dft"]),
                     masses_rel={k: m[k] / ref[k] - 1 for k in reg + ("m_VB-1 G-X", "m_VB-1 G-M")})
            d["ok"] = dict(gap=abs(d["gap_meV"]) <= LIMITS["gap_meV"],
                           valleys=max(abs(d["McG_meV"]), abs(d["XcG_meV"])) <= LIMITS["valley_meV"],
                           masses=max(abs(d["masses_rel"][k]) for k in reg) <= LIMITS["mass_rel"])
            d["pass"] = bool(all(d["ok"].values()))
            row["vs_W2"] = d
            plot["slab3"] = (Ew - Ew[:, nv - 1].max(), E - E[:, nv - 1].max(), nv, None)
        else:
            row["vs_HSE_gap"] = dict(hse_gap=hgap, d_meV=1e3 * (m["gap"] - hgap),
                                     ok=bool(abs(1e3 * (m["gap"] - hgap)) <= LIMITS["hse_gap_meV"]))
            plot["slab3 HSE"] = (Ehw - Ehw[:, nv - 1].max(), E - E[:, nv - 1].max(), nv, None)
        res["slab"][f"SK {func}"] = row
    (OUT / "gate_bd.json").write_text(json.dumps(res, indent=1, default=float))
    report(res)
    figures(plot, res)


def report(res):
    for kh, r in res["ribbons"].items():
        w = r["wannier"]
        print(f"=== {kh} (edge blocks from {r['calibration']}; W {r['W_A'] / 10:.2f} nm): W2 gap {w['gap']:.3f}, "
              f"valley X−Γ {w['valley_XG']:.3f}, m_VB {w['m_VB']:+.3f}, m_CB {w['m_CB']:+.3f}"
              + (f" | W2 without edge-s bands: gap {w['no_edge_s']['gap']:.3f}, valley {w['no_edge_s']['valley_XG']:.3f}"
                 if "no_edge_s" in w else ""))
        for lab in ("SK PBE", "SK, no edge blocks PBE"):
            m, d = r[lab]["metrics"], r[lab]["vs_W2"]
            print(f"   {lab}: gap {m['gap']:.3f} ({d['gap_meV']:+.0f} meV) valley {m['valley_XG']:.3f} "
                  f"({d['valley_meV']:+.0f}) m_VB {m['m_VB']:+.3f} ({d['m_VB_rel']:+.0%}) m_CB {m['m_CB']:+.3f} "
                  f"({d['m_CB_rel']:+.0%}) | I_ON rel " + " ".join(f"{k} {v:+.0%}" for k, v in r[lab]["I_on_rel"].items())
                  + f" | IRDS (SK, W) {r[lab]['irds_verdict_sk_vs_w']} | ok {d['ok']} → {'PASS' if d['pass'] else 'FAIL'}")
        for lab in ("SK HSE06", "SK, no edge blocks HSE06"):
            h = r[lab]["vs_HSE_gap"]
            print(f"   {lab}: gap {r[lab]['metrics']['gap']:.3f} vs HSE {h['hse_gap']:.3f} ({h['d_meV']:+.0f} meV) "
                  f"{'ok' if h['ok'] else 'FAIL'}")
        print("   W2 I_ON (mA/μm): " + " ".join(f"{k} {v['I_on']:.3f}" for k, v in w["I_on"].items())
              + " | k 4000 vs 2000: " + " ".join(f"{v:+.1e}" for v in r["tob_k_check"].values()))
    s = res["slab"]
    p, h = s["SK PBE"], s["SK HSE06"]
    d = p["vs_W2"]
    print(f"=== slab3: W2 gap {s['wannier_W2']['gap']:.3f}, M̄ {s['wannier_W2']['McG']:.3f}, X̄ {s['wannier_W2']['XcG']:.3f} | "
          f"SK PBE gap {p['metrics']['gap']:.3f} ({d['gap_meV']:+.0f} meV), M̄ {p['metrics']['McG']:.3f} ({d['McG_meV']:+.0f}), "
          f"X̄ {p['metrics']['XcG']:.3f} ({d['XcG_meV']:+.0f}; vs DFT {d['XcG_vs_dft_meV']:+.0f}); masses "
          + " ".join(f"{k[2:]} {v:+.0%}" for k, v in d["masses_rel"].items()) + f" | ok {d['ok']} → {'PASS' if d['pass'] else 'FAIL'}")
    print(f"   SK HSE06 gap {h['metrics']['gap']:.3f} vs HSE {h['vs_HSE_gap']['hse_gap']:.4f} "
          f"({h['vs_HSE_gap']['d_meV']:+.0f} meV) {'ok' if h['vs_HSE_gap']['ok'] else 'FAIL'}")


def figures(plot, res):
    """SM figures (re-rendered 2026-11-19 at the Supplementary text width, scripts/figstyle.py; no in-image titles):
    Fig3_holdouts_PBE_twin.png at 0.95 \\linewidth: (a) [110] n = 8, (b) [100] n = 6 (SK-PBE dashed blue, PBE-Wannier
    grey), (c) SK-HSE06 against the HSE06 gap points; Fig2c_slab_holdout.png at 0.95 \\linewidth: the three-plane slab,
    (a) PBE, (b) HSE06 (SK dashed blue, Wannier reference black)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import figstyle as FS
    FS.use(8.0)
    GREY, BLACK, SKC = "#8a8a8a", "#222222", "#2a6fdb"
    NAME = {"r110_8": "[110] $n$ = 8", "r100_6": "[100] $n$ = 6"}
    fig, axes = plt.subplots(1, 3, figsize=(FS.width(0.95, sm=True), 2.75), gridspec_kw=dict(width_ratios=[1, 1, 0.85]))
    for i, (ax, kh) in enumerate(zip(axes, ("r110_8", "r100_6"))):
        Ew, E, nv, Lx = plot[kh]
        t = np.linspace(0, 0.5, NK_PATH)
        ax.plot(t, Ew, color=GREY, lw=1.3)
        ax.plot(t, E, color=SKC, lw=0.8, ls=(0, (3.5, 1.8)))
        g = res["ribbons"][kh]["wannier"]["gap"]
        ax.set_ylim(-1.5, g + 1.5)
        ax.set_xlim(0, 0.5)
        ax.set_xticks([0, 0.5], ["Γ", "X"])
        ax.set_ylabel("$E - E_\\mathrm{VBM}$ (eV)")
        ax.axhline(0, color="#bbbbbb", lw=0.5, zorder=0)
        FS.letter(ax, "ab"[i], NAME[kh])
    ax = axes[2]
    labs = ["[110] $n$ = 8", "[100] $n$ = 6", "three-plane\nslab"]
    hs = [res["ribbons"]["r110_8"]["SK HSE06"]["vs_HSE_gap"]["hse_gap"], res["ribbons"]["r100_6"]["SK HSE06"]["vs_HSE_gap"]["hse_gap"],
          res["slab"]["hse_grid_gap"]]
    ss = [res["ribbons"]["r110_8"]["SK HSE06"]["metrics"]["gap"], res["ribbons"]["r100_6"]["SK HSE06"]["metrics"]["gap"],
          res["slab"]["SK HSE06"]["metrics"]["gap"]]
    y = np.arange(3)
    ax.barh(y + 0.18, hs, 0.34, color=BLACK, label="HSE06 (DFT)")
    ax.barh(y - 0.18, ss, 0.34, color=SKC, label="SK-HSE06")
    for yy, a, b in zip(y, hs, ss):
        ax.text(max(a, b) + 0.04, yy, f"{1e3 * (b - a):+.0f} meV", va="center", fontsize=8.0, color="#333333")
    ax.set_yticks(y, labs)
    ax.invert_yaxis()
    ax.set_xlabel("HSE06 gap (eV)")
    ax.set_xlim(0, max(hs + ss) * 1.62)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    FS.letter(ax, "c")
    from matplotlib.patches import Patch
    h = [Line2D([], [], color=GREY, lw=1.3), Line2D([], [], color=SKC, lw=0.8, ls=(0, (3.5, 1.8))),
         Patch(color=BLACK, lw=0), Patch(color=SKC, lw=0)]
    fig.legend(h, ["PBE-Wannier hold-out (W2)", "SK-PBE with edge blocks", "HSE06 gap point (DFT)", "SK-HSE06"],
               loc="outside lower center", ncol=2)
    FS.save(fig, ROOT / "figures" / "Fig3_holdouts_PBE_twin.png")

    fig, axes = plt.subplots(1, 2, figsize=(FS.width(0.95, sm=True), 2.85))
    for i, (ax, key, title) in enumerate(((axes[0], "slab3", "PBE"), (axes[1], "slab3 HSE", "HSE06"))):
        Ew, E, nv, _ = plot[key]
        s = np.arange(len(Ew))
        ax.plot(s, Ew, color=BLACK, lw=1.1)
        ax.plot(s, E, color=SKC, lw=0.8, ls=(0, (3.5, 1.8)))
        ax.set_xticks([0, 200, 400, 600], ["Γ", "X̄", "M̄", "Γ"])
        for xv in (200, 400):
            ax.axvline(xv, color="#dddddd", lw=0.5, zorder=0)
        ax.set_xlim(0, 600)
        ax.set_ylim(-2.0, 4.5)
        ax.set_ylabel("$E - E_\\mathrm{VBM}$ (eV)")
        ax.axhline(0, color="#bbbbbb", lw=0.5, zorder=0)
        FS.letter(ax, "ab"[i], title)
    h = [Line2D([], [], color=BLACK, lw=1.1), Line2D([], [], color=SKC, lw=0.8, ls=(0, (3.5, 1.8)))]
    fig.legend(h, ["Wannier reference (W2 / HW2)", "SK-PBE / SK-HSE06, surface on-site shifts from the sheet"],
               loc="outside lower center", ncol=2)
    FS.save(fig, ROOT / "figures" / "Fig2c_slab_holdout.png")


def figs():
    """Re-plot from gate_bd.json, recomputing only the band structures (no ToB)."""
    _frozen_md5()
    res = json.loads((OUT / "gate_bd.json").read_text())
    kf = np.c_[np.linspace(0, 0.5, NK_PATH), np.zeros((NK_PATH, 2))]
    plot = {}
    for kh, kc in HOLD.items():
        st = MBH.geometry(kh)
        nv = 3 * st.species.count("N")
        Ew = w_eigs(W.read_hr(M5.TWIN / "w90" / kh / "W2" / f"{kh}_W2_hr.dat"), kf)
        corr, _, _ = edge_corrections(MB, kc, EDGE_SRC["PBE"][kc])
        _, md, x = ribbon_sk(MBH, kh, "PBE", corr)
        E = sk_eigs(md, x, kf)
        plot[kh] = (Ew - Ew[:, nv - 1].max(), E - E[:, nv - 1].max(), nv, st.cell[0, 0])
    t3 = P3.frac_path()
    for key, func, path in (("slab3", "PBE", "dft/A3/w90/slab3/W2/slab3_W2_hr.dat"),
                            ("slab3 HSE", "HSE06", "dft/A3/w90/slab3/H/W2/slab3_HW2_hr.dat")):
        st, md, x, _ = slab_sk(func)
        nv = 3 * st.species.count("N")
        Ew, E = w_eigs(W.read_hr(ROOT / path), t3), sk_eigs(md, x, t3)
        plot[key] = (Ew - Ew[:, nv - 1].max(), E - E[:, nv - 1].max(), nv, None)
    figures(plot, res)


if __name__ == "__main__":
    dict(selftest=selftest, calib=calib, run=run, figs=figs)[sys.argv[1]]()
