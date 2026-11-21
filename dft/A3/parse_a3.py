#!/usr/bin/env python3
"""A3 analysis: the 3-plane slab hold-out (README.md, frozen 2026-11-15).

  vacuum    the vacuum check (A2 criterion 2 on the slab): PBE gap and valley offsets from the SCF grids of scf_slab3
            (20 Å) and vac_slab3 (27.5 Å) → a3_vacuum.json. Γ, X̄ and M̄ are matched modulo the reciprocal lattice
            (QE lists X̄ as (0, −0.5)).
  windows   after the PBE twin (outputs/, pot_slab3.cube.gz, w90/slab3/*.nnkp copied back): vacuum level, first
            non-model state, the W1 / W2 / W3 frozen and outer windows and their band-count checks, and the variant
            .win files (w90/slab3/<v>/slab3_<v>.win) → a3_windows.json.

The window rules are A4b's as frozen on 13 Nov (dft/A4b/README.md; parse_a4b.py ribbon_windows), in 2D:
  W1  frozen top = min(E_vac − 0.5 eV, minimum of the first band beyond the model count − 0.1 eV)
  W2  frozen top = min(W1, first non-model state − 0.1 eV); non-model: s_out = |ψ|² × (Löwdin share outside the
      projection set) + (1 − |ψ|²) ≥ 0.25, empty states on the Γ–X̄–M̄–Γ path (also reported at 0.20 and 0.30)
  W3  as W1, projection set S (B + Sc s on every Sc); run only under amendment 1's rule (README A3-3)
  outer = min(VBM + 10.23 eV, lowest energy of the top NSCF band − 0.05 eV); a cap below VBM + 10.23 is flagged
  E_vac: A2's rule for the sheet (parse_a2.vacuum_level): the bare + Hartree potential averaged over the plane, then
  over the outer 5 % of the cell at each end in z (the slab is centred); the spread is reported.

    ~/.venvs/scnfet/bin/python -P dft/A3/parse_a3.py {vacuum|windows}          (from the project root)
"""
import gzip
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, [str(path)]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


M3 = _load("make_a3", HERE / "make_a3.py")
PB = _load("parse_a4b", ROOT / "dft" / "A4b" / "parse_a4b.py")      # frozen constants and the projwfc parser
P = PB.P                                                             # parse_a2: bands(), alat_of()
KEY = M3.KEY
OUTER_PBE, VAC_MARGIN, NONMODEL_MARGIN, S_OUT, CHECK_THR = (PB.OUTER_PBE, PB.VAC_MARGIN, PB.NONMODEL_MARGIN,
                                                            PB.S_OUT, PB.CHECK_THR)
DEVICE = PB.DEVICE                                                   # [VBM − 0.5, CBM + 0.5] eV


def out(tag):
    return HERE / "outputs" / f"{tag}.out"


def nocc():
    st = M3.relaxed()
    return st, M3._mb().nocc(st)


# ---------------------------------------------------------------------------------------------------------- vacuum
def _grid_point(k, a, frac):
    """Index of the k point equal to `frac` (units of 2π/a, in-plane) modulo the reciprocal lattice, or None."""
    f = k[:, :2] * a / (2 * np.pi)
    d = f - np.array(frac)
    d -= np.round(d)
    n = np.linalg.norm(d, axis=1)
    return int(np.argmin(n)) if n.min() < 1e-4 else None


def levels(tag, nv):
    k, e, a = P.bands(out(tag))
    vbm, cbm = e[:, nv - 1].max(), e[:, nv].min()
    iG, iX, iM = (_grid_point(k, a, f) for f in ((0, 0), (0.5, 0), (0.5, 0.5)))
    if iX is None:
        iX = _grid_point(k, a, (0, 0.5))
    assert None not in (iG, iX, iM), f"{tag}: Γ, X̄ or M̄ missing from the grid"
    r = dict(nocc=nv, gap_grid=float(cbm - vbm), gap_G=float(e[iG, nv] - e[iG, nv - 1]),
             CB_X_minus_CB_G=float(e[iX, nv] - e[iG, nv]), Mbar_offset=float(e[iM, nv] - e[iG, nv]),
             VB_M_minus_VB_G=float(e[iM, nv - 1] - e[iG, nv - 1]),
             k_vbm=(k[np.argmax(e[:, nv - 1]), :2] * a / (2 * np.pi)).round(4).tolist(),
             k_cbm=(k[np.argmin(e[:, nv]), :2] * a / (2 * np.pi)).round(4).tolist())
    return r


def vacuum():
    _, nv = nocc()
    runs = {t: levels(t, nv) for t in (f"scf_{KEY}", f"vac_{KEY}")}
    a, b = runs[f"scf_{KEY}"], runs[f"vac_{KEY}"]
    ch = {q: round(1e3 * (b[q] - a[q]), 1) for q in ("gap_grid", "gap_G", "Mbar_offset", "CB_X_minus_CB_G")}
    res = dict(runs=runs, change_meV=ch, ok=bool(abs(ch["gap_grid"]) <= 5 and abs(ch["Mbar_offset"]) <= 5),
               note="A2 criterion 2 on the slab (README run-level criteria)")
    p = HERE / "a3_vacuum.json"
    if p.exists():
        old = json.loads(p.read_text())
        same = all(abs(old["runs"][t][q] - runs[t][q]) < 1e-9 for t in runs for q in ("gap_grid", "Mbar_offset"))
        print(f"a3_vacuum.json exists; this run {'reproduces it' if same else 'DIFFERS from it'}")
    p.write_text(json.dumps(res, indent=1))
    for t, r in runs.items():
        print(f"{t}: gap {r['gap_grid']:.4f} (Γ {r['gap_G']:.4f}) eV, M̄ offset {r['Mbar_offset']:.4f}, "
              f"X̄ − Γ {r['CB_X_minus_CB_G']:.4f}")
    print(f"change (meV): {ch} → {'pass' if res['ok'] else 'FAIL'}")


# --------------------------------------------------------------------------------------------------------- windows
def vacuum_level(cube_gz):
    """parse_a2.vacuum_level (A2's rule), reading the gzipped cube."""
    with gzip.open(cube_gz, "rt") as f:
        lines = f.read().splitlines()
    nat = abs(int(lines[2].split()[0]))
    n = [int(lines[3 + i].split()[0]) for i in range(3)]
    zat = np.array([float(lines[6 + j].split()[4]) for j in range(nat)]) * PB.BOHR
    dz = float(lines[5].split()[3]) * PB.BOHR
    vals = np.array(" ".join(lines[6 + nat:]).split(), float).reshape(n)        # x slowest, z fastest
    vz = vals.mean(axis=(0, 1)) * PB.RY
    nz = n[2]
    edge = np.r_[vz[: nz // 20], vz[-nz // 20:]]
    zs = np.r_[np.arange(nz // 20), np.arange(nz - nz // 20, nz)] * dz
    dd = np.abs(zs[:, None] - zat[None, :])
    gap_to_slab = float(np.minimum(dd, nz * dz - dd).min())               # periodic in z
    return dict(E_vac_eV=float(edge.mean()), spread_meV=1e3 * float(edge.max() - edge.min()), n_points=int(len(edge)),
                min_dist_to_atom_A=gap_to_slab)


def in_model(seed):
    def f(atom, sp, l):
        return (sp, l) in (("Sc", "d"), ("N", "p")) or (seed == "S" and (sp, l) == ("Sc", "s"))
    return f


def first_nonmodel(nv, vbm, cbm, seed, nk_path, path=None):
    inm = in_model(seed)
    states = PB.projwfc_states(None, path or out(f"projwfc_{KEY}"))
    assert len(states) == nk_path, f"projwfc has {len(states)} k, the run {nk_path}"
    rows = []
    for ik, sts in enumerate(states):
        for b in range(nv, len(sts)):
            E, w, p2 = sts[b]
            tot = sum(w.values())
            outside = sum(v for (a, sp, l), v in w.items() if not inm(a, sp, l)) / tot if tot > 0 else 1.0
            top = sorted(((f"{sp} {l}", v / tot) for (a, sp, l), v in w.items()), key=lambda t: -t[1]) if tot > 0 else []
            agg = {}
            for lab, v in top:
                agg[lab] = agg.get(lab, 0.0) + v
            rows.append((E, p2 * outside + (1 - p2), ik, b + 1, p2,
                         {k_: round(v, 3) for k_, v in sorted(agg.items(), key=lambda t: -t[1])[:3]}))
    rows.sort(key=lambda r: r[0])
    res = {}
    for thr in (S_OUT, *CHECK_THR):
        hit = next((r for r in rows if r[1] >= thr), None)
        res[f"{thr:.2f}"] = None if hit is None else dict(
            E_eV=hit[0], E_rel_vbm=hit[0] - vbm, E_rel_cbm=hit[0] - cbm, s_out=round(hit[1], 3), k_index=hit[2],
            band=hit[3], psi2=hit[4], top=hit[5])
    return res


def windows():
    info = json.loads((HERE / "a3_inputs.json").read_text())["stage2"]
    st, nv = nocc()
    ne = info["nexcl"]
    kp, ep, _ = P.bands(out(f"bands_{KEY}"))
    kg, eg, _ = P.bands(out(f"nscf_{KEY}"))
    assert eg.shape == (M3.KW ** 2, info["nbnd"]), f"NSCF grid {eg.shape}"
    vbm = float(max(ep[:, nv - 1].max(), eg[:, nv - 1].max()))
    cbm = float(min(ep[:, nv].min(), eg[:, nv].min()))
    vac = vacuum_level(HERE / f"pot_{KEY}.cube.gz")
    nm = {s: first_nonmodel(nv, vbm, cbm, s, len(kp)) for s in ("B", "S")}
    r = dict(nocc=nv, vbm=vbm, cbm=cbm, gap=cbm - vbm, gap_path=float(ep[:, nv].min() - ep[:, nv - 1].max()),
             vacuum=vac, E_vac_rel_vbm=vac["E_vac_eV"] - vbm, path_top_band_min_rel_vbm=float(ep[:, -1].min() - vbm),
             nonmodel=nm)
    top_min = float(eg[:, -1].min())
    outer = min(vbm + OUTER_PBE, top_min - 0.05)
    r.update(outer=outer, outer_rel_vbm=outer - vbm, outer_capped=bool(outer < vbm + OUTER_PBE - 1e-9),
             nscf_top_band_min_rel_vbm=top_min - vbm)
    fn = nm["B"][f"{S_OUT:.2f}"]
    r["amendment1_device_window"] = dict(
        device_top=cbm + DEVICE[1], first_nonmodel_B=None if fn is None else fn["E_eV"],
        inside=bool(fn is not None and fn["E_eV"] <= cbm + DEVICE[1]))
    var = {}
    for v, seed in (("W1", "B"), ("W2", "B"), ("W3", "S")):
        nw = info[f"num_wann_{seed}"]
        first_idx = float(eg[:, ne + nw].min()) if ne + nw < eg.shape[1] else np.inf
        froz = min(vac["E_vac_eV"] - VAC_MARGIN, first_idx - NONMODEL_MARGIN)
        binding = "vacuum" if froz == vac["E_vac_eV"] - VAC_MARGIN else "band count"
        if v == "W2" and fn is not None and fn["E_eV"] - NONMODEL_MARGIN < froz:
            froz, binding = fn["E_eV"] - NONMODEL_MARGIN, "non-model state"
        win = eg[:, ne:]
        n_froz, n_outer = (win <= froz).sum(axis=1), (win <= outer).sum(axis=1)
        var[v] = dict(seed=seed, num_wann=nw, froz=froz, froz_rel_vbm=froz - vbm, froz_rel_cbm=froz - cbm,
                      binding=binding, max_frozen_bands=int(n_froz.max()), min_outer_bands=int(n_outer.min()),
                      ok_counts=bool(n_froz.max() <= nw <= n_outer.min()), vac_margin=vac["E_vac_eV"] - froz,
                      covers_device=bool(froz >= cbm + DEVICE[1]))
        d = M3.W90 / KEY / v
        d.mkdir(parents=True, exist_ok=True)
        base = (M3.W90 / KEY / f"{KEY}_{seed}.win").read_text()
        text = re.sub(r"dis_froz_max = .*", f"dis_froz_max = {froz:.4f}", base)
        text = re.sub(r"dis_win_max = .*", f"dis_win_max = {outer:.4f}", text)
        (d / f"{KEY}_{v}.win").write_text(text)
    r["variants"] = var
    r["path_covers_frozen"] = bool(ep[:, -1].min() >= max(w["froz"] for w in var.values()) + 0.1)
    (HERE / "a3_windows.json").write_text(json.dumps(r, indent=1, default=float))
    print(f"=== {KEY}: VBM {vbm:.3f}, gap {r['gap']:.3f} eV (path {r['gap_path']:.3f}); E_vac − VBM "
          f"{r['E_vac_rel_vbm']:.3f} eV (spread {vac['spread_meV']:.2f} meV, {vac['n_points']} z points, ≥ "
          f"{vac['min_dist_to_atom_A']:.1f} Å from the atoms); outer VBM + {r['outer_rel_vbm']:.2f}"
          f"{' CAPPED' if r['outer_capped'] else ''}; path covers frozen + 0.1: {r['path_covers_frozen']}")
    for s in ("B", "S"):
        for thr, f in nm[s].items():
            print(f"   first non-model ({s}, s_out ≥ {thr}): " + ("none" if f is None else
                  f"VBM + {f['E_rel_vbm']:.3f} (CBM + {f['E_rel_cbm']:.3f}), band {f['band']}, k {f['k_index']}, "
                  f"|ψ|² {f['psi2']:.2f}, s_out {f['s_out']}, {f['top']}"))
    a1 = r["amendment1_device_window"]
    print(f"   amendment 1 (PBE): non-model state inside the device window (≤ CBM + 0.5): {a1['inside']}")
    for v, w in var.items():
        print(f"   {v} ({w['seed']}, {w['num_wann']} WFs): frozen VBM + {w['froz_rel_vbm']:.3f} = CBM + "
              f"{w['froz_rel_cbm']:.3f} eV [{w['binding']}], vacuum margin {w['vac_margin']:.2f}; bands frozen ≤ "
              f"{w['max_frozen_bands']}, in outer ≥ {w['min_outer_bands']} → {'ok' if w['ok_counts'] else 'FAIL'}")


# ----------------------------------------------------------------------------------------- criteria (PBE twin, 2D)
# README "Run-level criteria" (frozen 15 Nov): P1–P5 as frozen, the reference the NC PBE bands on Γ–X̄–M̄–Γ; P5 at both
# band edges along Γ–X̄ and Γ–M̄; P6 descriptive; H4 between variants that freeze the same window, 200 points per
# segment; decision rule W2 → W1 → back for review. Limits and helpers are parse_a4b's (LIM, edge_mass, wout_centres).
ORB = {"s": ["s"], "p": ["pz", "px", "py"], "d": ["dz2", "dxz", "dyz", "dx2-y2", "dxy"]}   # Wannier90's m order
NSEG = 200


def frac_path(n=NSEG):
    """Γ–X̄–M̄–Γ in fractional coordinates, n points per segment (3n + 1 in all)."""
    nodes = np.array([[0, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0], [0, 0, 0]], float)
    return np.vstack([nodes[i] + np.outer(np.linspace(0, 1, n, endpoint=False), nodes[i + 1] - nodes[i]) for i in range(3)]
                     + [nodes[-1:]])


def wf_map(win_text, st):
    """Atom and orbital of each WF from the .win projections block (Wannier90 order: per projection line, the atoms of
    that species in atoms_cart order, then each l and m of the line)."""
    blk = re.search(r"begin projections\s*\n(.*?)end projections", win_text, re.S | re.I).group(1)
    m = []
    for line in blk.splitlines():
        line = line.split("!")[0].strip()
        if not line:
            continue
        sp, ls = (t.strip() for t in line.split(":"))
        for a, s in enumerate(st.species):
            if s == sp:
                for l in ls.split(";"):
                    m += [dict(atom=a, species=sp, orbital=o) for o in ORB[l.strip()]]
    return m


def centres(wout, win_text, st):
    """P6 (descriptive): WF offsets from their atoms. With one k along z, the Berry-phase centre of a WF on a plane at
    exactly z = c/2 sits on the branch cut: Wannier90 then reports z ≈ 0 and a spread ≈ (c/2)². Those WFs (the middle
    plane of the centred slab) are flagged and only their in-plane offset is kept. The Hamiltonian does not use the
    centres (use_ws_distance = .false.; kz = 0 only)."""
    c, spr, om = PB.wout_centres(wout)
    mp = wf_map(win_text, st)
    assert len(mp) == len(c), f"{len(mp)} projections, {len(c)} WFs"
    cell = np.diag(st.cell)
    offs = []
    for (r_, s_), m in zip(zip(c, spr), mp):
        d = r_ - st.pos[m["atom"]]
        d -= cell * np.round(d / cell)
        cut = bool(abs(st.pos[m["atom"]][2] - cell[2] / 2) < 1e-3)
        assert np.linalg.norm(d[:2]) < 0.5, f"WF on atom {m['atom']} {m['orbital']}: in-plane offset {d[:2]}"
        offs.append(dict(**m, branch_cut=cut, offset_A=float(np.linalg.norm(d[:2]) if cut else np.linalg.norm(d)),
                         d_A=d.round(4).tolist(), spread_A2=None if cut else float(s_), spread_reported_A2=float(s_)))
    ok = [o for o in offs if not o["branch_cut"]]
    host = [o for o in ok if not (o["species"] == "Sc" and o["orbital"] in ("dz2", "s"))]
    return dict(omega_total_reported=om, n_branch_cut=len(offs) - len(ok),
                omega_off_cut=float(sum(o["spread_A2"] for o in ok)),
                host_nonz2_max_offset_A=max(o["offset_A"] for o in host),
                z2_offsets_A=[round(o["offset_A"], 3) for o in ok if o["orbital"] == "dz2"],
                s_offsets_A=[round(o["offset_A"], 3) for o in ok if o["orbital"] == "s"],
                max_spread_off_cut_A2=max(o["spread_A2"] for o in ok), wf=offs)


def _fit_curv(s, e, idx, i0):
    return np.linalg.lstsq(np.vstack([np.ones(len(idx)), (s[idx] - s[i0]) ** 2]).T, e[idx], rcond=None)[0][1]


def criteria():
    from scnfet import w90 as W
    info = json.loads((HERE / "a3_inputs.json").read_text())["stage2"]
    win = json.loads((HERE / "a3_windows.json").read_text())
    st, nv = nocc()
    ne = info["nexcl"]
    kp, ep, _ = P.bands(out(f"bands_{KEY}"))
    kf = kp @ st.cell.T / (2 * np.pi)
    iG, iX, iM = 0, 60, 120
    assert np.allclose(kf[[iG, iX, iM, 180]], frac_path(1), atol=1e-4), "path nodes"
    seg = {"G-X": np.arange(0, 61), "G-M": np.arange(180, 119, -1)}       # Γ first in both
    vbm, cbm = win["vbm"], win["cbm"]
    iv, ic = nv - 1 - ne, nv - ne
    res = {}
    for v in ("W1", "W2", "W3"):
        d = M3.W90 / KEY / v
        seed = f"{KEY}_{v}"
        if not ((d / f"{seed}_hr.dat").exists() and (d / f"{seed}.wout").exists()):
            continue
        var = win["variants"][v]
        nw = var["num_wann"]
        hr = W.read_hr(d / f"{seed}_hr.dat")
        assert hr.num_wann == nw
        Ew = np.linalg.eigvalsh(W.hk_frac(hr, kf))
        ref = ep[:, ne:ne + nw]
        nb = min(ref.shape[1], nw)
        ref, Ew = ref[:, :nb], Ew[:, :nb]
        err = np.abs(Ew - ref)
        dev = (ref >= vbm - DEVICE[0]) & (ref <= cbm + DEVICE[1])
        froz = ref <= var["froz"]
        gap_dft, gap_w = ref[:, ic].min() - ref[:, iv].max(), Ew[:, ic].min() - Ew[:, iv].max()
        masses = {}
        for sn, ix in seg.items():
            s = np.linalg.norm(kp[ix] - kp[ix[0]], axis=1)
            for lab, b, kind in (("VB", iv, "v"), ("CB", ic, "c"), ("VB-1", iv - 1, "v")):
                m_d, i0, n, idx = PB.edge_mass(s, ref[ix], b, kind)
                m_w = PB.HB2M / _fit_curv(s, Ew[ix][:, b], idx, i0)
                masses[f"{lab} {sn}"] = dict(m_dft=float(m_d), m_w90=float(m_w), rel=float(m_w / m_d - 1),
                                             k_index=int(ix[i0]), at_gamma=bool(i0 == 0), n_points=n)
        reg = [k for k in masses if not k.startswith("VB-1")]                 # registered: the two band edges
        val = {lab: dict(dft=float(ref[i, ic] - ref[iG, ic]), w90=float(Ew[i, ic] - Ew[iG, ic]))
               for lab, i in (("Ec(M)-Ec(G)", iM), ("Ec(X)-Ec(G)", iX))}
        r = res[v] = dict(
            P1=dict(vac_margin_eV=var["vac_margin"], counts_ok=var["ok_counts"], outer_capped=win["outer_capped"],
                    ok=bool(var["vac_margin"] >= VAC_MARGIN - 1e-9 and var["ok_counts"] and not win["outer_capped"])),
            P2=dict(max_meV=1e3 * float(err[dev].max()), n=int(dev.sum()), ok=bool(1e3 * err[dev].max() <= PB.LIM["device_max_meV"])),
            P3=dict(rms_meV=1e3 * float(np.sqrt(np.mean(err[froz] ** 2))), max_meV=1e3 * float(err[froz].max()),
                    n=int(froz.sum()), ok=bool(1e3 * np.sqrt(np.mean(err[froz] ** 2)) <= PB.LIM["frozen_rms_meV"])),
            P4=dict(gap_dft=float(gap_dft), gap_w90=float(gap_w), d_meV=1e3 * float(gap_w - gap_dft),
                    ok=bool(1e3 * abs(gap_w - gap_dft) <= PB.LIM["gap_meV"])),
            P5=dict(masses=masses, registered=reg, vb_degeneracy_G_meV=1e3 * float(ref[iG, iv] - ref[iG, iv - 1]),
                    ok=bool(all(abs(masses[k]["rel"]) <= PB.LIM["mass_rel"] for k in reg))),
            P6=centres(d / f"{seed}.wout", (d / f"{seed}.win").read_text(), st),
            valleys=val)
        r["pass_P1_P5"] = bool(all(r[k]["ok"] for k in ("P1", "P2", "P3", "P4", "P5")))
    if res.get("W2", {}).get("pass_P1_P5"):
        res["decision"] = "W2"
    elif res.get("W1", {}).get("pass_P1_P5"):
        res["decision"] = "W1"
    else:
        res["decision"] = "neither: LR6, back for review" if "W1" in res and "W2" in res else "incomplete"
    # H4: variants that freeze the same window (here: every pair that both freeze the whole device window)
    t = frac_path()
    E = {v: np.linalg.eigvalsh(W.hk_frac(W.read_hr(M3.W90 / KEY / v / f"{KEY}_{v}_hr.dat"), t))
         for v in ("W1", "W2", "W3") if v in res}
    h4 = {}
    vs = list(E)
    for i, a in enumerate(vs):
        for b in vs[i + 1:]:
            n = min(E[a].shape[1], E[b].shape[1])
            top = min(cbm + DEVICE[1], win["variants"][a]["froz"], win["variants"][b]["froz"])
            m = (E[a][:, :n] >= vbm - DEVICE[0]) & (E[a][:, :n] <= top)
            dd = np.abs(E[a][:, :n] - E[b][:, :n])
            h4[f"{a}-{b}"] = dict(max_meV=1e3 * float(dd[m].max()), top_rel_cbm=top - cbm, n=int(m.sum()),
                                  k_index_at_max=int(np.unravel_index(np.argmax(np.where(m, dd, -1)), dd.shape)[0]),
                                  ok=bool(1e3 * dd[m].max() <= PB.H4_MEV))
    res["H4"] = dict(points_per_segment=NSEG, pairs=h4)
    (HERE / "a3_criteria.json").write_text(json.dumps(res, indent=1, default=float))
    print(f"=== {KEY}: decision {res['decision']}")
    for v in ("W1", "W2", "W3"):
        if v not in res:
            continue
        c = res[v]
        ms = " ".join(f"{k} {c['P5']['masses'][k]['rel']:+.1%}" for k in c["P5"]["registered"])
        print(f"   {v}: P1 {c['P1']['ok']} | P2 device max {c['P2']['max_meV']:.2f} meV ({c['P2']['n']}) {c['P2']['ok']} | "
              f"P3 frozen RMS {c['P3']['rms_meV']:.2f} (max {c['P3']['max_meV']:.1f}) {c['P3']['ok']} | P4 gap "
              f"{c['P4']['d_meV']:+.2f} meV {c['P4']['ok']} | P5 {ms} {c['P5']['ok']}")
        print(f"       masses (m_e, DFT / W90): " + "; ".join(
            f"{k} {m['m_dft']:+.3f}/{m['m_w90']:+.3f}{'' if m['at_gamma'] else ' (extremum off Γ, k ' + str(m['k_index']) + ')'}"
            for k, m in c["P5"]["masses"].items()) + f"; VB degeneracy at Γ {c['P5']['vb_degeneracy_G_meV']:.2f} meV")
        p6 = c["P6"]
        print(f"       valleys (DFT / W90): " + "; ".join(f"{k} {x['dft']:.4f}/{x['w90']:.4f}" for k, x in c["valleys"].items())
              + f" | P6 host max {p6['host_nonz2_max_offset_A']:.3f} Å, z² {p6['z2_offsets_A']}, max spread "
              f"{p6['max_spread_off_cut_A2']:.2f} Å² off the cut ({p6['n_branch_cut']} WFs on the z = c/2 branch cut; "
              f"reported Ω {p6['omega_total_reported']:.1f}, off-cut sum {p6['omega_off_cut']:.1f} Å²)")
    for pair, c in h4.items():
        print(f"   H4 {pair}: max |ΔE| {c['max_meV']:.2f} meV over [VBM − 0.5, CBM + {c['top_rel_cbm']:.2f}] ({c['n']} values,"
              f" {NSEG}/segment) → {'ok' if c['ok'] else 'FAIL: back for review'}")


# ---------------------------------------------------------------------------------------- HSE windows (H1, H3, a1)
# README A3-3 (frozen 15 Nov): windows by the same rules on the HSE eigenvalues (E_vac from the HSE potential,
# non-model states from projwfc on the HSE SCF), outer VBM + 11.05 eV not capped (H1); H3 / amendment 1: an HSE
# non-model state inside [VBM − 0.5, CBM + 0.5] → W3 (set S) as well, used if it passes H1, H2 and covers the window,
# otherwise W2 with a stated limitation (reported, not a stop). Writes w90/slab3/H/<v>/slab3_H<v>.win (the layout of
# dft/A4b/run_a4b2_w90.sh) and a3_hwindows.json.
def hwindows():
    info = json.loads((HERE / "a3_inputs.json").read_text())["stage2"]
    pbe = json.loads((HERE / "a3_windows.json").read_text())
    st, nv = nocc()
    ne = info["nexcl"]
    txt = out(f"hse_{KEY}").read_text()
    if "!!   total energy" not in txt or re.search(r"Maximum CPU time exceeded", txt):
        raise SystemExit(f"{KEY}: the HSE SCF did not converge (no '!!' energy, or a max_seconds stop); not used")
    dexx = [float(x) for x in re.findall(r"est\. exchange err \(dexx\)\s*=\s*([-\d.E+]+)", txt)]
    kh, eh, _ = P.bands(out(f"hse_{KEY}"))
    vbm, cbm = float(eh[:, nv - 1].max()), float(eh[:, nv].min())
    lv = levels(f"hse_{KEY}", nv)
    vac = vacuum_level(HERE / f"hpot_{KEY}.cube.gz")
    nm = {s: first_nonmodel(nv, vbm, cbm, s, len(kh), out(f"hproj_{KEY}")) for s in ("B", "S")}
    top_min = float(eh[:, -1].min())
    outer = min(vbm + PB.HSE_OUTER, top_min - 0.05)
    r = dict(nbnd=int(eh.shape[1]), nk_ibz=int(eh.shape[0]), exx_outer_steps=len(dexx), dexx_final_Ry=dexx[-1],
             vbm=vbm, cbm=cbm, gap_grid=cbm - vbm, levels=lv, vacuum=vac, E_vac_rel_vbm=vac["E_vac_eV"] - vbm,
             nonmodel=nm, outer=outer, outer_rel_vbm=outer - vbm,
             outer_capped=bool(outer < vbm + PB.HSE_OUTER - 1e-9), top_band_min_rel_vbm=top_min - vbm)
    fn = nm["B"][f"{S_OUT:.2f}"]
    dev_top = cbm + DEVICE[1]
    r["H3"] = dict(device_top=dev_top, first_nonmodel_rel_cbm=None if fn is None else fn["E_eV"] - cbm,
                   inside=bool(fn is not None and fn["E_eV"] <= dev_top),
                   inside_at_0p20=bool(nm["B"]["0.20"] is not None and nm["B"]["0.20"]["E_eV"] <= dev_top),
                   inside_at_0p30=bool(nm["B"]["0.30"] is not None and nm["B"]["0.30"]["E_eV"] <= dev_top))
    H = M3.W90 / KEY / "H"
    bases = {s: (H / f"{KEY}_{n}.win").read_text() for s, n in (("B", "H"), ("S", "HS"))}
    for s, b in bases.items():
        if "@KLIST@" in b:
            raise SystemExit(f"{KEY}: the set-{s} HSE .win still has @KLIST@; copy back the cluster copy")
    variants = [("W1", "B"), ("W2", "B")] + ([("W3", "S")] if r["H3"]["inside"] else [])
    stale = H / "W3" / f"{KEY}_HW3.win"
    if not r["H3"]["inside"] and stale.exists():
        raise SystemExit(f"{stale} exists but amendment 1 does not apply; remove it by hand before the job")
    var = {}
    for v, seed in variants:
        nw = info[f"num_wann_{seed}"]
        first_idx = float(eh[:, ne + nw].min()) if ne + nw < eh.shape[1] else np.inf
        froz = min(vac["E_vac_eV"] - VAC_MARGIN, first_idx - NONMODEL_MARGIN)
        binding = "vacuum" if froz == vac["E_vac_eV"] - VAC_MARGIN else "band count"
        if v == "W2" and fn is not None and fn["E_eV"] - NONMODEL_MARGIN < froz:
            froz, binding = fn["E_eV"] - NONMODEL_MARGIN, "non-model state"
        w = eh[:, ne:]
        n_froz, n_outer = (w <= froz).sum(axis=1), (w <= outer).sum(axis=1)
        var[v] = dict(seed=seed, num_wann=nw, froz=froz, froz_rel_vbm=froz - vbm, froz_rel_cbm=froz - cbm,
                      binding=binding, covers_device=bool(froz >= dev_top), max_frozen_bands=int(n_froz.max()),
                      min_outer_bands=int(n_outer.min()), ok_counts=bool(n_froz.max() <= nw <= n_outer.min()),
                      vac_margin=vac["E_vac_eV"] - froz)
        var[v]["H1_ok"] = bool(var[v]["vac_margin"] >= VAC_MARGIN - 1e-9 and var[v]["ok_counts"] and not r["outer_capped"])
        d = H / v
        d.mkdir(parents=True, exist_ok=True)
        text = re.sub(r"dis_froz_max = .*", f"dis_froz_max = {froz:.4f}", bases[seed])
        text = re.sub(r"dis_win_max = .*", f"dis_win_max = {outer:.4f}", text)
        (d / f"{KEY}_H{v}.win").write_text(text)
    r["variants"] = var
    (HERE / "a3_hwindows.json").write_text(json.dumps(r, indent=1, default=float))
    pl = json.loads((HERE / "a3_vacuum.json").read_text())["runs"][f"scf_{KEY}"]
    print(f"=== {KEY} HSE: {len(dexx)} EXX steps, dexx {dexx[-1]:.1e} Ry; gap (IBZ grid) {r['gap_grid']:.4f} eV (PBE "
          f"{pbe['gap']:.4f}); Γ gap {lv['gap_G']:.4f}; M̄ offset {lv['Mbar_offset']:.4f} (PBE {pl['Mbar_offset']:.4f}); "
          f"Ec(X̄) − Ec(Γ) {lv['CB_X_minus_CB_G']:.4f} (PBE {pl['CB_X_minus_CB_G']:.4f}); VBM k {lv['k_vbm']}, CBM k {lv['k_cbm']}")
    print(f"   E_vac − VBM {r['E_vac_rel_vbm']:.3f} eV (PBE {pbe['E_vac_rel_vbm']:.3f}; spread {vac['spread_meV']:.2f} meV); "
          f"outer VBM + {r['outer_rel_vbm']:.2f}{' CAPPED' if r['outer_capped'] else ''} (top band min VBM + "
          f"{r['top_band_min_rel_vbm']:.2f})")
    for s in ("B", "S"):
        for thr, f in nm[s].items():
            print(f"   first non-model ({s}, s_out ≥ {thr}): " + ("none" if f is None else
                  f"CBM + {f['E_rel_cbm']:.3f} eV, band {f['band']}, k {f['k_index']}, |ψ|² {f['psi2']:.2f}, "
                  f"s_out {f['s_out']}, {f['top']}"))
    h3 = r["H3"]
    print(f"   H3 / amendment 1: non-model state inside the device window: {h3['inside']} (at 0.20 / 0.30: "
          f"{h3['inside_at_0p20']} / {h3['inside_at_0p30']}) → W3 {'written' if h3['inside'] else 'not run'}")
    for v, w in var.items():
        print(f"   H{v} ({w['seed']}, {w['num_wann']} WFs): frozen CBM + {w['froz_rel_cbm']:.3f} eV [{w['binding']}], vacuum "
              f"margin {w['vac_margin']:.2f}; bands frozen ≤ {w['max_frozen_bands']}, in outer ≥ {w['min_outer_bands']}; "
              f"covers device {w['covers_device']} → H1 {'ok' if w['H1_ok'] else 'FAIL'}")


# --------------------------------------------------------------------------------------- HSE criteria (H2, H4)
# README A3-3 (frozen 15 Nov): H2 Wannier90 reproduces the HSE grid eigenvalues inside the frozen window to ≤ 1 meV;
# H4 between variants that freeze the same window, ≤ 5 meV over the device window, 200 points per segment; decision
# W2 → W1 → back for review. Gap, valleys, masses (both edges, Γ–X̄ and Γ–M̄) and centres are descriptive, on the fine
# path (no HSE band path exists; the PBE twin tested them).
def hcriteria():
    from scnfet import w90 as W
    info = json.loads((HERE / "a3_inputs.json").read_text())["stage2"]
    win = json.loads((HERE / "a3_hwindows.json").read_text())
    st, nv = nocc()
    ne = info["nexcl"]
    kh, eh, _ = P.bands(out(f"hse_{KEY}"))
    kf = kh @ st.cell.T / (2 * np.pi)
    t = frac_path()
    tc = t @ np.linalg.inv(st.cell.T) * 2 * np.pi                        # Å⁻¹, for the masses
    iG, iX, iM = 0, NSEG, 2 * NSEG
    seg = {"G-X": np.arange(0, NSEG + 1), "G-M": np.arange(3 * NSEG, 2 * NSEG - 1, -1)}
    iv, ic = nv - 1 - ne, nv - ne
    vbm, cbm = win["vbm"], win["cbm"]
    res, E = {}, {}
    for v in ("W1", "W2", "W3"):
        if v not in win["variants"]:
            continue
        d = M3.W90 / KEY / "H" / v
        seed = f"{KEY}_H{v}"
        if not ((d / f"{seed}_hr.dat").exists() and (d / f"{seed}.wout").exists()):
            print(f"H{v}: Wannier output missing")
            continue
        var = win["variants"][v]
        nw = var["num_wann"]
        hr = W.read_hr(d / f"{seed}_hr.dat")
        assert hr.num_wann == nw
        # H2 on the full open_grid (the .win k list) with the .eig eigenvalues (12 digits) and the frozen top as written
        # in the .win: the set Wannier90 froze. QE's stdout has 4 decimals, so an IBZ state within 0.1 meV of the top
        # can be counted on the wrong side; the IBZ check is kept as a cross-check (H2_ibz).
        wt = (d / f"{seed}.win").read_text()
        fw = float(re.search(r"dis_froz_max = (\S+)", wt).group(1))
        kw = np.array([[float(x) for x in l.split()[:3]] for l in
                       re.search(r"begin kpoints\s*\n(.*?)end kpoints", wt, re.S | re.I).group(1).strip().splitlines()])
        eig = np.loadtxt(M3.W90 / KEY / "H" / f"{KEY}_H.eig")
        nbw = int(eig[:, 0].max())
        Ee = eig[:, 2].reshape(len(kw), nbw)[:, :nw]                      # bands after exclude_bands, ascending
        assert np.allclose(eig[:, 1].reshape(len(kw), nbw)[:, 0], np.arange(1, len(kw) + 1))
        Ewg = np.linalg.eigvalsh(W.hk_frac(hr, kw))
        fzg = Ee <= fw
        errg = np.abs(Ewg - Ee)[fzg]
        ref = eh[:, ne:ne + nw]
        Eg = np.linalg.eigvalsh(W.hk_frac(hr, kf))
        froz = ref <= var["froz"]
        err = np.abs(Eg - ref)[froz]
        edge = np.abs(ref - fw)[froz] < 1e-4                               # IBZ states at the top, within stdout precision
        Ef = E[v] = np.linalg.eigvalsh(W.hk_frac(hr, t))
        masses = {}
        for sn, ix in seg.items():
            s = np.linalg.norm(tc[ix] - tc[ix[0]], axis=1)
            for lab, b, kind in (("VB", iv, "v"), ("CB", ic, "c"), ("VB-1", iv - 1, "v")):
                m_, i0, n, _ = PB.edge_mass(s, Ef[ix], b, kind)
                masses[f"{lab} {sn}"] = dict(m=float(m_), at_gamma=bool(i0 == 0), n_points=n)
        res[v] = dict(H1_ok=var["H1_ok"],
                      H2=dict(max_meV=1e3 * float(errg.max()), n=int(errg.size), nk=len(kw), froz_win=fw,
                              ok=bool(1e3 * errg.max() <= PB.H2_MEV)),
                      H2_ibz=dict(max_meV=1e3 * float(err.max()), n=int(err.size),
                                  max_meV_without_top_states=1e3 * float(err[~edge].max()), n_top_states=int(edge.sum())),
                      gap_path=float(Ef[:, ic].min() - Ef[:, iv].max()),
                      vbm_at_G=bool(np.argmax(Ef[:, iv]) == iG), cbm_at_G=bool(np.argmin(Ef[:, ic]) == iG),
                      valleys={"Ec(M)-Ec(G)": float(Ef[iM, ic] - Ef[iG, ic]), "Ec(X)-Ec(G)": float(Ef[iX, ic] - Ef[iG, ic])},
                      masses=masses, P6=centres(d / f"{seed}.wout", (d / f"{seed}.win").read_text(), st))
    ok = {v: r["H1_ok"] and r["H2"]["ok"] for v, r in res.items()}
    h4 = {}
    vs = list(E)
    for i, a in enumerate(vs):
        for b in vs[i + 1:]:
            n = min(E[a].shape[1], E[b].shape[1])
            top = min(cbm + DEVICE[1], win["variants"][a]["froz"], win["variants"][b]["froz"])
            m = (E[a][:, :n] >= vbm - DEVICE[0]) & (E[a][:, :n] <= top)
            dd = np.abs(E[a][:, :n] - E[b][:, :n])
            h4[f"{a}-{b}"] = dict(max_meV=1e3 * float(dd[m].max()), top_rel_cbm=top - cbm, n=int(m.sum()),
                                  ok=bool(1e3 * dd[m].max() <= PB.H4_MEV))
    if win["H3"]["inside"]:
        dec = "amendment 1 applies: W3 if it passes H1, H2 and covers the device window (decide by hand)"
    elif ok.get("W2"):
        dec = "W2"
    elif ok.get("W1"):
        dec = "W1"
    else:
        dec = "neither passes H1/H2: back for review" if "W1" in res and "W2" in res else "incomplete"
    if any(not c["ok"] for c in h4.values()):
        dec += "; H4 FAILS: back for review"
    res.update(decision=dec, H4=dict(points_per_segment=NSEG, pairs=h4))
    (HERE / "a3_hcriteria.json").write_text(json.dumps(res, indent=1, default=float))
    print(f"=== {KEY} HSE: decision {dec}")
    for v in ("W1", "W2", "W3"):
        if v not in res:
            continue
        c = res[v]
        print(f"   H{v}: H1 {c['H1_ok']} | H2 full grid max {c['H2']['max_meV']:.3f} meV ({c['H2']['n']} values, {c['H2']['nk']} k) "
              f"{c['H2']['ok']} [IBZ stdout {c['H2_ibz']['max_meV']:.3f}; without the {c['H2_ibz']['n_top_states']} state(s) at the "
              f"top {c['H2_ibz']['max_meV_without_top_states']:.3f}] | "
              f"gap {c['gap_path']:.4f} eV (VBM at Γ {c['vbm_at_G']}, CBM at Γ {c['cbm_at_G']}) | valleys "
              + ", ".join(f"{k} {x:.4f}" for k, x in c["valleys"].items()))
        print("       m* " + "; ".join(f"{k} {m['m']:+.3f}{'' if m['at_gamma'] else ' (off Γ)'}" for k, m in c["masses"].items())
              + f" | z² {c['P6']['z2_offsets_A']} Å, off-cut Ω {c['P6']['omega_off_cut']:.1f} Å² ({c['P6']['n_branch_cut']} on the cut)")
    for pair, c in h4.items():
        print(f"   H4 {pair}: max |ΔE| {c['max_meV']:.2f} meV over [VBM − 0.5, CBM + {c['top_rel_cbm']:.2f}] → {'ok' if c['ok'] else 'FAIL'}")


if __name__ == "__main__":
    modes = dict(vacuum=vacuum, windows=windows, criteria=criteria, hwindows=hwindows, hcriteria=hcriteria)
    if len(sys.argv) > 1 and sys.argv[1] in modes:
        modes[sys.argv[1]]()
    else:
        print(__doc__)
