"""A4a analysis. stage1: the three relaxations (convergence, magnetisation history, edge geometry = Gate M criterion 6,
a preliminary gap on the SCF grid). Gate M criteria 1–5 need A4a-2 (README.md).

    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P dft/A4a/parse_a4a.py        (after copying outputs/ back)
Writes dft/A4a/a4a_stage1.json.
"""
import importlib.util
import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


MK = _load("make_a4a", HERE / "make_a4a.py")
P = _load("parse_a2", HERE.parent / "A2" / "parse_a2.py")
BOND = 2.6          # Å: Sc–N neighbour cut (ideal in-plane 2.199, interplane 2.1685)
WALL_H = {"r100_3": 2 + 58 / 60 + 16 / 3600, "r110_4": 1 + 42 / 60 + 54 / 3600, "r100_3h": 4 + 22 / 60 + 13 / 3600}


def bonds(st, i):
    d = st.pos - st.pos[i]
    d[:, 0] -= st.cell[0, 0] * np.round(d[:, 0] / st.cell[0, 0])          # periodic along the ribbon (x)
    r = np.linalg.norm(d, axis=1)
    return [(j, float(r[j])) for j in np.argsort(r) if 0 < r[j] < BOND and st.species[j] != "H"]


def stage1():
    res = {}
    for key, (_, _, h) in MK.RIBBONS.items():
        txt = (HERE / "outputs" / f"relax_{key}.out").read_text()
        r = dict(converged="End final coordinates" in txt,
                 bfgs_steps=[int(v) for v in re.findall(r"bfgs converged in\s+\d+ scf cycles and\s+(\d+) bfgs steps", txt)],
                 E_final_Ry=float(re.findall(r"Final energy\s+=\s+(-?[\d.]+)", txt)[-1]),
                 core_h=round(64 * WALL_H[key], 1))
        if not h:
            ab = [float(v) for v in re.findall(r"absolute magnetization\s+=\s+(-?[\d.]+)", txt)]
            r.update(abs_mag_first_iter=ab[0], abs_mag_max=max(ab), abs_mag_last=ab[-1],
                     zero_moment_restart=bool(re.search(r"final configuration with zero\s+absolute magnetization", txt)))
        st0, st = MK.structure(key), MK.relaxed(key)
        sp = np.array(st.species)
        disp = np.linalg.norm(st.pos - st0.pos, axis=1)
        host = sp != "H"
        y = st0.pos[:, 1]
        ylo, yhi = y[host].min(), y[host].max()
        edge = host & ((np.abs(y - ylo) < 0.1) | (np.abs(y - yhi) < 0.1))
        inner = host & ~edge
        eb = sorted({round(b, 3) for i in np.where(edge)[0] for _, b in bonds(st, i)})
        ib = sorted({round(b, 3) for i in np.where(inner)[0] for _, b in bonds(st, i)})
        # plane separation at the edges and in the middle row (Sc/N directly above each other)
        zsep = {}
        for lab, mask in (("edge", edge), ("inner", inner)):
            gaps = []
            for i in np.where(mask)[0]:
                for j, b in bonds(st, i):
                    d = st.pos[j] - st.pos[i]
                    if abs(d[2]) > 1.5 and abs(d[1]) < 0.5:
                        gaps.append(abs(float(d[2])))
            zsep[lab] = [round(min(gaps), 3), round(max(gaps), 3)] if gaps else None
        r.update(n_atoms=len(sp), max_disp_A=float(disp.max()), max_disp_edge_A=float(disp[edge].max()),
                 max_disp_inner_A=float(disp[inner].max()) if inner.any() else None,
                 edge_bonds_A=[eb[0], eb[-1]], inner_bonds_A=[ib[0], ib[-1]], plane_sep_A=zsep,
                 width_A=float(st.pos[host, 1].max() - st.pos[host, 1].min()), width0_A=float(yhi - ylo))
        if h:
            hb = {"N": [], "Sc": []}
            for i in np.where(sp == "H")[0]:
                d = st.pos - st.pos[i]
                d[:, 0] -= st.cell[0, 0] * np.round(d[:, 0] / st.cell[0, 0])
                rr = np.linalg.norm(d, axis=1)
                rr[i] = 9e9
                j = int(np.argmin(rr))
                hb[sp[j]].append(round(float(rr[j]), 3))
            r["H_bond_A"] = {k: [min(v), max(v)] for k, v in hb.items() if v}
            r["max_disp_H_A"] = float(disp[sp == "H"].max())
        # preliminary gap: last SCF, 8 irreducible k of the 15-point grid (X̄ is not on it)
        k, e, _ = P.bands(HERE / "outputs" / f"relax_{key}.out")
        nel = float(re.findall(r"number of electrons\s*=\s*([\d.]+)", txt)[-1])
        nocc = int(round(nel)) // 2
        e = e[-8:] if len(e) > 8 else e                      # last k block (nspin 1 at the end in every case)
        vb, cb = e[:, nocc - 1], e[:, nocc]
        r.update(nocc=nocc, gap_grid_eV=float(cb.min() - vb.max()), k_vbm=int(np.argmax(vb)), k_cbm=int(np.argmin(cb)),
                 direct_gap_G_eV=float(cb[0] - vb[0]))
        res[key] = r
    (HERE / "a4a_stage1.json").write_text(json.dumps(res, indent=1))
    for key, r in res.items():
        print(key, json.dumps(r))
    return res


# ---------------------------------------------------------------------------------------------------------------
# stage 2: Gate M criteria 1–5 (README.md, pre-registered 2026-11-11). Conventions the README leaves open, fixed
# 2026-11-12 before criteria 2–5 were computed (the criterion-1 energies had been seen in a job status check):
#   - edge atoms = host atoms of the outermost row of each edge in the ideal cut (|y − y_edge| < 0.1 Å), both planes:
#     8 for r100_3 and r100_3h, 4 for r110_4. For r100_3h the H weight is reported separately and also added.
#   - k-test gap (criterion 2) = NM gap on the k 15 points, from the k 15 and the k 30 densities (the k 15 grid is a
#     subset of the k 30 one, so this isolates density convergence); both full-grid gaps are reported too.
#   - criterion 3 metallic flag: band nocc or nocc+1 reaches E_F (of the NM k 15 SCF) anywhere on the path.
#   - criterion 4 weight = Σ edge-atom Löwdin components / Σ all printed components, for band nocc at the VBM k and
#     band nocc+1 at the CBM k on the 31-k path (projwfc ran on the bands run); the uniform share is reported beside it.
#   - criterion 5 = |grid gap(vac 20) − grid gap(vac 15)| for r110_4 NM at k 15 (identical k sets).
# Decisions after the results (2026-11-12): (a) for H ribbons the H atoms count as the edge row in criterion 4
# (the flag uses edge_plus_H; the host-only share is kept); (b) criterion 5 is also run on the [100] ribbons
# (vac20_r100_3, vac20_r100_3h), same 5 meV limit; [110] ribbons use 20 Å; (c) the [110] H ribbon r110_4h (20 Å) goes
# through the same analysis once its runs exist.
#   - A4b band survey (descriptive): smallest N whose band minimum on the path is ≥ VBM + L, L = the sheet's PBE level
#     for the accepted N 40 window (cost test T2: 10.23 + 0.05 eV), also at L + 0.5 and L + 1.0 eV.
RY_EV = P.RY
WIN = json.loads((HERE.parent / "A2" / "costtest" / "windows.json").read_text())
L_PBE = WIN["pbe_k16_nb40"]["outer_rel_vbm"] + 0.05
SHEET_PBE_GAP = 0.404
THRESH = dict(mag_meV_per_edge=5.0, ktest_mag_meV=1.0, ktest_gap_meV=10.0, edge_weight=0.5, vac_gap_meV=5.0)


def out(tag):
    return HERE / "outputs" / f"{tag}.out"


def done(tag):
    p = out(tag)
    return p.exists() and "JOB DONE" in p.read_text()


def scf(tag):
    txt = out(tag).read_text()
    r = dict(E_Ry=float(re.findall(r"^!+\s+total energy\s+=\s+(-?[\d.]+)", txt, re.M)[-1]),
             converged="convergence has been achieved" in txt,
             E_F_eV=float(re.findall(r"the Fermi energy is\s+(-?[\d.]+)", txt)[-1]))
    ab = re.findall(r"absolute magnetization\s+=\s+(-?[\d.]+)", txt)
    if ab:
        blk = txt[txt.rfind("Magnetic moment per site"):]
        site = [float(v) for v in re.findall(r"magn=\s*(-?[\d.]+)", blk.split("\n\n")[0])]
        r.update(abs_mag=float(ab[-1]), tot_mag=float(re.findall(r"total magnetization\s+=\s+(-?[\d.]+)", txt)[-1]),
                 abs_mag_first_iter=float(ab[0]), max_site_moment=float(np.abs(site).max()) if site else None)
    return r


def nocc_of(tag):
    return int(round(float(re.search(r"number of electrons\s*=\s*([\d.]+)", out(tag).read_text()).group(1)))) // 2


def edge_masks(key):
    st0 = MK.structure(key)
    sp = np.array(st0.species)
    host = sp != "H"
    y = st0.pos[:, 1]
    ylo, yhi = y[host].min(), y[host].max()
    left = host & (np.abs(y - ylo) < 0.1)
    right = host & (np.abs(y - yhi) < 0.1)
    return left, right, sp == "H"


def grid_gap(tag, nocc, kref=None):
    k, e, _ = P.bands(out(tag))
    if kref is not None:                                   # restrict to the points of another grid
        sel = [int(np.argmin(np.linalg.norm(k - q, axis=1))) for q in kref]
        assert all(np.linalg.norm(k[s] - q) < 1e-6 for s, q in zip(sel, kref)), "k 15 point missing from the k 30 grid"
        k, e = k[sel], e[sel]
    return float(e[:, nocc].min() - e[:, nocc - 1].max()), k


def projwfc(tag):
    """[(k, [(E, {atom: weight}, psi2), ...]), ...] from a projwfc.x stdout (atoms 1-based → 0-based)."""
    txt = out(tag).read_text()
    atom = {int(s): int(a) - 1 for s, a in re.findall(r"state #\s*(\d+): atom\s+(\d+)", txt)}
    res = []
    for kb in re.split(r"\n k =", txt.split("Lowdin Charges")[0])[1:]:
        kv = [float(v) for v in kb.split("\n")[0].split()]
        st = []
        for m in re.finditer(r"==== e\(\s*\d+\) =\s*(-?[\d.]+) eV ====(.*?)\|psi\|\^2 = ([\d.]+)", kb, re.S):
            w = {}
            for c, s in re.findall(r"([\d.]+)\*\[#\s*(\d+)\]", m.group(2)):
                w[atom[int(s)]] = w.get(atom[int(s)], 0.0) + float(c)
            st.append((float(m.group(1)), w, float(m.group(3))))
        res.append((kv, st))
    return res


def edge_weight(w, left, right, hmask):
    tot = sum(w.values())
    f = lambda m: sum(v for a, v in w.items() if m[a]) / tot
    return dict(edge=f(left | right), left=f(left), right=f(right), H=f(hmask), edge_plus_H=f(left | right | hmask))


def stage2():
    res = {}
    for key, (_, _, h) in MK.ALL.items():
        if not done(f"nm_{key}"):
            continue
        r = {}
        left, right, hmask = edge_masks(key)
        n_edge = int((left | right).sum())
        nat = len(left)
        r["n_edge_atoms"] = n_edge
        nocc = nocc_of(f"nm_{key}")
        r["nocc"] = nocc
        # 1. magnetism (k 15: FM, AFM-edges, AFM-planes; k 30: FM, AFM-edges; H ribbon: FM at k 15 only)
        runs = {t: scf(t) for t in ([f"nm_{key}", f"fm_{key}"] + ([] if h else [
                f"afm_edges_{key}", f"afm_planes_{key}", f"nm_k30_{key}", f"fm_k30_{key}", f"afm_edges_k30_{key}"]))
                if done(t)}
        r["runs"] = runs

        def dmag(kk):
            nm = runs.get(f"nm{kk}_{key}")
            mags = {t: v for t, v in runs.items() if v.get("abs_mag") is not None and (("_k30_" in t) == bool(kk))}
            if not nm or not mags:
                return None
            tmin = min(mags, key=lambda t: mags[t]["E_Ry"])
            return dict(dE_meV_per_edge=1e3 * RY_EV * (nm["E_Ry"] - mags[tmin]["E_Ry"]) / n_edge, lowest=tmin,
                        lowest_abs_mag=mags[tmin]["abs_mag"], max_site_moment=max(v["max_site_moment"] for v in mags.values()),
                        all_collapsed=all(v["abs_mag"] < 0.01 for v in mags.values()))
        r["c1_k15"] = dmag("")
        r["c1_k30"] = None if h else dmag("_k30")
        c1 = r["c1_k30"] or r["c1_k15"]
        r["c1_magnetic"] = bool(c1["dE_meV_per_edge"] > THRESH["mag_meV_per_edge"])
        r["c1_basis"] = "k30" if r["c1_k30"] else "k15 (H ribbon: no k 30 magnetic runs, as pre-registered)"
        # 2. k-test (bare ribbons)
        if not h and r["c1_k30"]:
            g15, k15 = grid_gap(f"nm_{key}", nocc)
            g30_common, _ = grid_gap(f"nm_k30_{key}", nocc, kref=k15)
            g30_full, _ = grid_gap(f"nm_k30_{key}", nocc)
            ddm = abs(r["c1_k30"]["dE_meV_per_edge"] - r["c1_k15"]["dE_meV_per_edge"])
            r["c2"] = dict(dE_mag_change_meV_per_edge=ddm, gap_k15_eV=g15, gap_k30_on_k15_points_eV=g30_common,
                           gap_k30_full_grid_eV=g30_full, gap_change_meV=1e3 * abs(g30_common - g15),
                           ok=bool(ddm <= THRESH["ktest_mag_meV"] and 1e3 * abs(g30_common - g15) <= THRESH["ktest_gap_meV"]))
        # 3. PBE gap from the 31-k path + the SCF grid
        if done(f"bands_{key}"):
            kp, ep, _ = P.bands(out(f"bands_{key}"))
            x = kp[:, 0] / kp[-1, 0]
            vb, cb = ep[:, nocc - 1], ep[:, nocc]
            iv, ic = int(np.argmax(vb)), int(np.argmin(cb))
            ef = runs[f"nm_{key}"]["E_F_eV"]
            gap = float(cb[ic] - vb[iv])
            r["c3"] = dict(gap_path_eV=gap, k_vbm=round(float(x[iv]), 4), k_cbm=round(float(x[ic]), 4),
                           direct=bool(iv == ic), direct_gap_G_eV=float(cb[0] - vb[0]),
                           direct_gap_X_eV=float(cb[-1] - vb[-1]), gap_grid_k15_eV=grid_gap(f"nm_{key}", nocc)[0],
                           E_F_eV=ef, vbm_eV=float(vb[iv]), metallic=bool(gap <= 0 or vb.max() >= ef or cb.min() <= ef),
                           vs_sheet_eV=gap - SHEET_PBE_GAP, nbnd=int(ep.shape[1]))
            # band survey for A4b's nbnd (descriptive)
            emin = ep.min(axis=0) - vb[iv]
            n_sc = int((np.array(MK.structure(key).species) == "Sc").sum())
            sv = {}
            for lab, L in (("L", L_PBE), ("L+0.5", L_PBE + 0.5), ("L+1.0", L_PBE + 1.0)):
                ok = np.where(emin >= L)[0]
                N = int(ok[0]) + 1 if len(ok) else None
                sv[lab] = dict(level_rel_vbm=round(L, 3), nbnd=N, per_fu=round(N / n_sc, 2) if N else None)
            r["survey"] = dict(levels=sv, top_band_min_rel_vbm=float(emin[-1]), n_sc=n_sc)
            # 4. edge states
            if done(f"projwfc_{key}"):
                pw = projwfc(f"projwfc_{key}")
                assert len(pw) == len(kp), "projwfc k count differs from the bands run"
                assert abs(pw[iv][1][nocc - 1][0] - vb[iv]) < 2e-3 and abs(pw[ic][1][nocc][0] - cb[ic]) < 2e-3
                ew_v = edge_weight(pw[iv][1][nocc - 1][1], left, right, hmask)
                ew_c = edge_weight(pw[ic][1][nocc][1], left, right, hmask)
                along = {}
                for lab, b in (("VB-1", nocc - 2), ("VB", nocc - 1), ("CB", nocc), ("CB+1", nocc + 1)):
                    along[lab] = dict(G=round(edge_weight(pw[0][1][b][1], left, right, hmask)["edge"], 3),
                                      X=round(edge_weight(pw[-1][1][b][1], left, right, hmask)["edge"], 3),
                                      max=round(max(edge_weight(s[b][1], left, right, hmask)["edge"] for _, s in pw), 3),
                                      E_G_rel_vbm=round(pw[0][1][b][0] - vb[iv], 3),
                                      E_X_rel_vbm=round(pw[-1][1][b][0] - vb[iv], 3))
                ek = "edge_plus_H" if h else "edge"
                r["c4"] = dict(vbm=ew_v, cbm=ew_c, psi2_vbm=pw[iv][1][nocc - 1][2], psi2_cbm=pw[ic][1][nocc][2],
                               uniform_share=n_edge / nat, uniform_share_host=n_edge / int((~hmask).sum()),
                               along_path=along, edge_definition=ek,
                               edge_state_sets_gap=bool(ew_v[ek] > THRESH["edge_weight"] or ew_c[ek] > THRESH["edge_weight"]),
                               host_row_only_flag=bool(ew_v["edge"] > THRESH["edge_weight"] or ew_c["edge"] > THRESH["edge_weight"]))
        res[key] = r
    # 5. vacuum: r110_4 (pre-registered); r100_3, r100_3h (decision after it failed)
    res["c5_vacuum"] = {}
    for key in MK.RIBBONS:
        if done(f"vac20_{key}"):
            n = nocc_of(f"nm_{key}")
            g15, g20 = grid_gap(f"nm_{key}", n)[0], grid_gap(f"vac20_{key}", n)[0]
            res["c5_vacuum"][key] = dict(gap_vac15_eV=g15, gap_vac20_eV=g20, change_meV=1e3 * abs(g20 - g15),
                                         dE_tot_meV=1e3 * RY_EV * (scf(f"vac20_{key}")["E_Ry"] - scf(f"nm_{key}")["E_Ry"]),
                                         ok=bool(1e3 * abs(g20 - g15) <= THRESH["vac_gap_meV"]))
    # 5b. 20 → 25 Å on the two ribbons that failed 15 → 20 Å (2026-11-13, with A4b-1; same 5 meV limit). Γ eigenvalue
    # shifts, aligned on the median shift of the occupied bands (as vac_bands.py), say which bands move.
    res["c5_vacuum25"] = {}
    for key in ("r100_3h", "r110_4"):
        if done(f"vac25_{key}") and done(f"vac20_{key}"):
            n = nocc_of(f"nm_{key}")
            g20, g25 = grid_gap(f"vac20_{key}", n)[0], grid_gap(f"vac25_{key}", n)[0]
            E = {}
            for t in (f"vac20_{key}", f"vac25_{key}"):
                k, e, _ = P.bands(out(t))
                E[t] = e[int(np.argmin(np.linalg.norm(k, axis=1)))]
            d = E[f"vac25_{key}"] - E[f"vac20_{key}"]
            d -= np.median(d[:n])
            res["c5_vacuum25"][key] = dict(gap_vac20_eV=g20, gap_vac25_eV=g25, change_meV=1e3 * abs(g25 - g20),
                                           dE_tot_meV=1e3 * RY_EV * (scf(f"vac25_{key}")["E_Ry"] - scf(f"vac20_{key}")["E_Ry"]),
                                           gamma_shift_meV_bands_nocc_m2_to_p1={int(b + 1): round(1e3 * float(d[b]), 2)
                                                                                for b in range(n - 2, n + 1)},
                                           change_15_20_meV=res["c5_vacuum"].get(key, {}).get("change_meV"),
                                           ok=bool(1e3 * abs(g25 - g20) <= THRESH["vac_gap_meV"]))
    res["thresholds"] = THRESH
    res["L_pbe_rel_vbm"] = L_PBE
    (HERE / "a4a_stage2.json").write_text(json.dumps(res, indent=1))
    report(res)
    return res


def report(res):
    for key in [k for k in MK.ALL if k in res]:
        r = res[key]
        print(f"\n=== {key}  (edge atoms {r['n_edge_atoms']}, nocc {r['nocc']})")
        for t, v in r["runs"].items():
            print(f"   {t:22s} E {v['E_Ry']:.8f} Ry  conv {v['converged']}"
                  + (f"  |m| {v['abs_mag']:.2f} (iter1 {v['abs_mag_first_iter']:.2f}) site max {v['max_site_moment']}"
                     if "abs_mag" in v else ""))
        for kk in ("c1_k15", "c1_k30"):
            if r.get(kk):
                c = r[kk]
                print(f"   {kk}: ΔE_mag {c['dE_meV_per_edge']:+.3f} meV/edge atom (lowest {c['lowest']}, |m| "
                      f"{c['lowest_abs_mag']}), all collapsed {c['all_collapsed']}")
        print(f"   C1 magnetic: {r['c1_magnetic']} ({r['c1_basis']})")
        if "c2" in r:
            c = r["c2"]
            print(f"   C2 k-test: ΔΔE_mag {c['dE_mag_change_meV_per_edge']:.3f} meV/edge, gap k15 {c['gap_k15_eV']:.4f} → "
                  f"k30 (same k) {c['gap_k30_on_k15_points_eV']:.4f} [full k30 grid {c['gap_k30_full_grid_eV']:.4f}] "
                  f"Δ {c['gap_change_meV']:.1f} meV → ok {c['ok']}")
        if "c3" in r:
            c = r["c3"]
            print(f"   C3 PBE gap {c['gap_path_eV']:.4f} eV ({'direct' if c['direct'] else 'indirect'}; VBM k {c['k_vbm']}, "
                  f"CBM k {c['k_cbm']} of Γ–X̄); Γ {c['direct_gap_G_eV']:.4f}, X̄ {c['direct_gap_X_eV']:.4f}; grid "
                  f"{c['gap_grid_k15_eV']:.4f}; metallic {c['metallic']}; vs sheet {c['vs_sheet_eV']:+.3f} eV")
            s = r["survey"]
            print(f"   band survey (nbnd {c['nbnd']}, top band min VBM+{s['top_band_min_rel_vbm']:.2f} eV): "
                  + "; ".join(f"{k} (VBM+{v['level_rel_vbm']:.2f}) N {v['nbnd']} = {v['per_fu']}/f.u." for k, v in s["levels"].items()))
        if "c4" in r:
            c = r["c4"]
            print(f"   C4 edge weight VBM {c['vbm']['edge']:.3f} (L {c['vbm']['left']:.3f} R {c['vbm']['right']:.3f} H "
                  f"{c['vbm']['H']:.3f}), CBM {c['cbm']['edge']:.3f} (L {c['cbm']['left']:.3f} R {c['cbm']['right']:.3f} H "
                  f"{c['cbm']['H']:.3f}); uniform {c['uniform_share']:.3f}; |psi|² {c['psi2_vbm']:.3f}/{c['psi2_cbm']:.3f}"
                  f" → edge state sets gap ({c['edge_definition']}): {c['edge_state_sets_gap']}")
            for b, v in c["along_path"].items():
                print(f"      {b:5s} edge Γ {v['G']:.3f} X̄ {v['X']:.3f} max {v['max']:.3f}   E−VBM Γ {v['E_G_rel_vbm']:+.3f} "
                      f"X̄ {v['E_X_rel_vbm']:+.3f}")
    for key, c in res["c5_vacuum"].items():
        print(f"C5 vacuum {key}: gap {c['gap_vac15_eV']:.4f} → {c['gap_vac20_eV']:.4f} eV, Δ {c['change_meV']:.1f} meV "
              f"→ ok {c['ok']}  (ΔE_tot {c['dE_tot_meV']:+.1f} meV)")
    for key, c in res.get("c5_vacuum25", {}).items():
        print(f"C5 vacuum 20→25 Å {key}: gap {c['gap_vac20_eV']:.4f} → {c['gap_vac25_eV']:.4f} eV, Δ {c['change_meV']:.1f} meV "
              f"(15→20: {c['change_15_20_meV']:.1f}) → ok {c['ok']}  (ΔE_tot {c['dE_tot_meV']:+.1f} meV); Γ shifts "
              f"{c['gamma_shift_meV_bands_nocc_m2_to_p1']}")


if __name__ == "__main__":
    import sys
    stage2() if "stage2" in sys.argv else stage1()
