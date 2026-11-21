"""Post-hoc diagnostic (not registered; changes no verdict): where is the H self-rebuild window deviation?

`t1_bb_widths.py hcalib` stores window_dev_meV = max |ΔE| over bands within [VBM − 1, CBM + 1] eV, compared band index
by band index after sorting: 312 meV for [110] H and 285 meV for [100] H in HSE06 (53 / 62 meV in PBE). This script
rebuilds the same two ribbons with the same code path and reports, per functional:
  * the (band, k) of the largest sorted-index deviation, its energy from the VBM / CBM, and the band's orbital character;
  * the same maximum restricted to the device window [VBM − 0.5, CBM + 0.5] eV;
  * an overlap-matched deviation (bands paired by eigenvector overlap; both Hamiltonians share one WF basis), which
    separates a real band shift from an index offset caused by a misplaced band crossing in or out of the window.
Run: PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/diag_hcalib_dev.py   → fits/t1/hcalib_dev.json
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, [str(path)]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


T1 = _load("t1", ROOT / "scripts" / "t1_bb_widths.py")
B, G, MB, MK = T1.B, T1.G, T1.MB, T1.MK
from scnfet import bb  # noqa: E402

KF = T1.KF
kx = KF[:, 0]


def character(V, wf, C, atoms_edge):
    """Weights of one eigenvector (column) on: H s, Sc d, N p, Sc s; edge-row atoms vs centre rows."""
    w = np.abs(V) ** 2
    out = {"H s": 0.0, "Sc d": 0.0, "N p": 0.0, "Sc s": 0.0, "edge rows": 0.0, "centre rows": 0.0}
    for j, (a, o) in enumerate(wf):
        sp = C.species[a]
        lab = "H s" if sp == "H" else ("Sc s" if o == "s" else ("Sc d" if sp == "Sc" else "N p"))
        out[lab] += w[j]
        out["edge rows" if a in atoms_edge else "centre rows"] += w[j]
    return {k: round(float(v), 3) for k, v in out.items()}


def one(kc, func):
    sheet, csrc, C, shift = T1.sources(kc, func)
    st = MB.geometry(kc)
    nv = T1.nval(st)
    Er, Vr = G.w_eigs(csrc.hr, KF, vecs=True)
    Er = Er - csrc.E_vac
    hr, wf, used = bb.assemble(C, C, csrc, sheet, B.q_sheet(MK.structure(kc)), B.RCUT, sheet_shift=shift)
    Eb, Vb = G.w_eigs(hr, KF, vecs=True)
    # same WF basis? (assemble uses wf_list order; the calibration source uses its wf_map)
    same_basis = list(map(tuple, csrc.wf)) == list(map(tuple, wf))
    perm = None
    if not same_basis:
        pos = {tuple(x): i for i, x in enumerate(csrc.wf)}
        perm = [pos[tuple(x)] for x in wf]               # rebuild index j ↔ reference index perm[j]
    vbm, cbm = Er[:, nv - 1].max(), Er[:, nv:].min()
    # edge rows: the outermost host row on each side plus the H row (row 0 / last are H for the H ribbons)
    edge_atoms = {a for a in range(len(C.species)) if C.row[a] in (0, 1, C.nrow - 2, C.nrow - 1)}
    res = dict(nv=int(nv), vbm_eV=float(vbm), cbm_eV=float(cbm), gap_eV=float(cbm - vbm), nrow=int(C.nrow),
               same_wf_order=bool(same_basis), sheet_shift_eV=float(shift))
    dE = Eb - Er
    for lab, (lo, hi) in (("window_1eV", (1.0, 1.0)), ("device_window_0p5eV", (0.5, 0.5))):
        m = (Er >= vbm - lo) & (Er <= cbm + hi)
        A = np.where(m, np.abs(dE), -1)
        ik, ib = np.unravel_index(np.argmax(A), A.shape)
        vec = Vr[ik][:, ib]
        res[lab] = dict(max_dev_meV=float(1e3 * abs(dE[ik, ib])), sign=float(np.sign(dE[ik, ib])), k_frac=float(kx[ik]),
                        band=int(ib), band_rel=(f"VB-{nv - 1 - ib}" if ib < nv else f"CB+{ib - nv}"),
                        E_ref_minus_VBM_eV=float(Er[ik, ib] - vbm), E_ref_minus_CBM_eV=float(Er[ik, ib] - cbm),
                        char_ref=character(vec, csrc.wf, C, edge_atoms),
                        max_dev_at_band_edges_meV=dict(
                            VB_top=float(1e3 * np.abs(dE[:, nv - 1]).max()), CB_bottom=float(1e3 * np.abs(dE[:, nv]).max())))
    # overlap matching: pair reference band i with rebuild band j maximising Σ|<ψ_i|φ_j>|² at each k
    m1 = (Er >= vbm - 1.0) & (Er <= cbm + 1.0)
    best = dict(dev=-1.0)
    per_band = {}
    for ik in range(len(kx)):
        Vr_k = Vr[ik] if perm is None else Vr[ik][perm, :]
        O = np.abs(Vr_k.conj().T @ Vb[ik]) ** 2           # (ref band, rebuild band)
        ri, cj = linear_sum_assignment(-O)
        for i, j in zip(ri, cj):
            if not m1[ik, i]:
                continue
            d = abs(Eb[ik, j] - Er[ik, i])
            per_band.setdefault(int(i), []).append(d)
            if d > best["dev"]:
                best = dict(dev=d, ik=ik, i=int(i), j=int(j), ov=float(O[i, j]))
    ik, i = best["ik"], best["i"]
    res["overlap_matched_1eV"] = dict(
        max_dev_meV=float(1e3 * best["dev"]), k_frac=float(kx[ik]), band=i,
        band_rel=(f"VB-{nv - 1 - i}" if i < nv else f"CB+{i - nv}"), rebuild_band=best["j"], overlap=round(best["ov"], 3),
        E_ref_minus_VBM_eV=float(Er[ik, i] - vbm), E_ref_minus_CBM_eV=float(Er[ik, i] - cbm),
        char_ref=character(Vr[ik][:, i], csrc.wf, C, edge_atoms),
        per_band_max_meV={(f"VB-{nv - 1 - b}" if b < nv else f"CB+{b - nv}"): round(1e3 * max(v), 1)
                          for b, v in sorted(per_band.items())})
    # energy profile: max |ΔE| (sorted index) by distance from the nearer band edge, 0.1 eV bins out to 1 eV
    prof = {}
    for side, sel, dist in (("CB", np.s_[:, nv:], lambda E: E - cbm), ("VB", np.s_[:, :nv], lambda E: vbm - E)):
        D, A = dist(Er[sel]), np.abs(dE[sel])
        prof[side] = {f"{b / 10:.1f}-{(b + 1) / 10:.1f}": round(float(1e3 * A[(D >= b / 10) & (D < (b + 1) / 10)].max()), 1)
                      if ((D >= b / 10) & (D < (b + 1) / 10)).any() else None for b in range(10)}
    res["profile_meV_by_eV_from_edge"] = prof
    # ToB on-current (registered settings, parallel plate): rebuild vs the ribbon's own HW1 / W1 Hamiltonian
    Lx, Wd, nat = st.cell[0, 0], T1.S.width(st), T1.n_host(st)
    t_own = T1.tob_all(G.w_eigs(csrc.hr, T1.KG), nv, Lx, Wd, nat)
    t_bb = T1.tob_all(G.w_eigs(hr, T1.KG), nv, Lx, Wd, nat)
    res["tob_I_on"] = {c: dict(own=t_own[c]["I_on"], rebuild=t_bb[c]["I_on"], rel=t_bb[c]["I_on"] / t_own[c]["I_on"] - 1)
                       for c in ("n EOT 0.6", "p EOT 0.6", "n EOT 1.0", "p EOT 1.0")}
    return res


if __name__ == "__main__":
    out = {}
    for kc in ("r110_4h", "r100_3h"):
        for func in ("HSE06", "PBE"):
            r = one(kc, func)
            out[f"{kc} {func}"] = r
            w, d, o = r["window_1eV"], r["device_window_0p5eV"], r["overlap_matched_1eV"]
            print(f"{kc} {func}: ±1 eV max {w['max_dev_meV']:.0f} meV at {w['band_rel']} k={w['k_frac']:.3f} "
                  f"(E−VBM {w['E_ref_minus_VBM_eV']:+.2f}, E−CBM {w['E_ref_minus_CBM_eV']:+.2f}); "
                  f"±0.5 eV max {d['max_dev_meV']:.0f} meV at {d['band_rel']} k={d['k_frac']:.3f}; "
                  f"overlap-matched max {o['max_dev_meV']:.0f} meV at {o['band_rel']} k={o['k_frac']:.3f} "
                  f"(ov {o['overlap']})", flush=True)
            print("   char (±1 eV argmax):", w["char_ref"], " band edges:", w["max_dev_at_band_edges_meV"], flush=True)
            print("   profile CB:", r["profile_meV_by_eV_from_edge"]["CB"], flush=True)
            print("   profile VB:", r["profile_meV_by_eV_from_edge"]["VB"], flush=True)
            print("   I_ON:", {c: f"{v['own']:.3f}→{v['rebuild']:.3f} ({v['rel']:+.1%})" for c, v in r["tob_I_on"].items()},
                  flush=True)
    (ROOT / "fits" / "t1" / "hcalib_dev.json").write_text(json.dumps(out, indent=1, default=float))
