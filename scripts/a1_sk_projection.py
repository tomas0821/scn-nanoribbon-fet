"""A1 → SK, step 1: map the bulk Wannier H(R) onto the SK model's orbital order and lattice gauge, check the
orbital signs, and project it onto the SK forms (plan D3: the fit's initial guess).

  * W90 order for `Sc:d; N:p`: dz2, dxz, dyz, dx2-y2, dxy | pz, px, py (w90.W90_ORDER); the model's
    (sk.ORBITALS): dxy, dyz, dzx, dx2-y2, dz2 | px, py, pz.
  * W90 cell: rows (−a/2,0,a/2), (0,a/2,a/2), (−a/2,a/2,0), N at (½,½,½) → (−a/2,a/2,a/2); the model's: fcc rows
    a/2(0,1,1)…, N at (a/2,0,0). `w90.regauge` keeps every element's Cartesian bond vector.
  * Sign check: projection-only WFs inherit the phases of the trial orbitals, so the W90 real harmonics should map
    with all signs +1. Tested by flipping each orbital in turn: the in-model residual must be smallest with none.
  * Projection: `skfit.project_hr` on the shells NN pd, 2NN dd/pp, 3NN dd (Harrison scale 1 at the bulk a).
    Reported: in-model RMS (three-centre / crystal-field content the two-centre form cannot carry), out-of-model
    norm by shell distance (longer-range hops dropped by the truncation), and the projected model's bands.

Run from the project root: PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/a1_sk_projection.py
Writes fits/A1_bulk/sk_projection.json.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

from scnfet import sk, skfit, w90

ROOT = Path(__file__).resolve().parents[1]
A1 = ROOT / "dft" / "A1"
OUT = ROOT / "fits" / "A1_bulk"
A = 4.501
SEEDS = {"HSE06": "hse_k16", "PBE": "pbe_k16"}
MODEL_ORDER = [o for o in sk.ORBITALS["Sc"]] + [o for o in sk.ORBITALS["N"]]
W90_LIST = list(w90.W90_ORDER["d"]) + list(w90.W90_ORDER["p"])
PERM = [W90_LIST.index(o.replace("dzx", "dxz")) for o in MODEL_ORDER]


def win_geometry(win):
    txt = win.read_text()
    cell = np.array([[float(v) for v in l.split()] for l in
                     re.search(r"begin unit_cell_cart\s*\n\s*ang\s*\n(.*?)end unit_cell_cart", txt, re.S | re.I)
                     .group(1).strip().splitlines()])
    atoms = re.search(r"begin atoms_frac\s*\n(.*?)end atoms_frac", txt, re.S | re.I).group(1).strip().splitlines()
    frac = {l.split()[0]: np.array([float(v) for v in l.split()[1:4]]) for l in atoms}
    return cell, {s: f @ cell for s, f in frac.items()}


def to_model(seed, model):
    d = A1 / "w90" / seed
    hr = w90.read_hr(d / f"{seed}_hr.dat")
    cell_w, pos_w = win_geometry(d / f"{seed}.win")
    cen_w = np.array([pos_w["Sc"]] * 5 + [pos_w["N"]] * 3)
    st = model.structure
    cen_m = np.array([st.pos[0]] * 5 + [st.pos[1]] * 3)
    return hr, cell_w, cen_w, w90.regauge(hr, cell_w, cen_w, st.cell, cen_m, perm=PERM), cen_m


def kpath(model, n=80):
    P = sk.fcc_points(A)
    return sk.path(P, ["L", "G", "X", "W", "K", "G"], n=n)


def main():
    model = sk.build(sk.rocksalt(A), shells=sk.rocksalt_shells(A))
    st = model.structure
    out = {"a": A, "perm_model_from_w90": PERM, "model_order": MODEL_ORDER, "w90_order": W90_LIST}
    kp, s, ticks = kpath(model)
    for lab, seed in SEEDS.items():
        hr, cell_w, cen_w, hm, cen_m = to_model(seed, model)
        # 1. the regauged H gives the same bands (same Cartesian k)
        e_w = np.linalg.eigvalsh(w90.hk_frac(hr, kp @ cell_w.T / (2 * np.pi)))
        e_m = np.linalg.eigvalsh(w90.hk_frac(hm, st.kfrac(kp)))
        vbm = np.linalg.eigvalsh(w90.hk_frac(hm, np.zeros((1, 3))))[0, 2]
        rec = {"regauge_max_dE_meV": 1e3 * float(np.max(np.abs(e_w - e_m))), "vbm_eV": float(vbm)}
        # 2. on-site block: cubic symmetry → diagonal, t2g and e_g each degenerate
        H0 = hm.H[np.where((hm.R == 0).all(1))[0][0]].real
        rec["onsite_diag_rel_vbm"] = dict(zip(MODEL_ORDER, np.round(np.diag(H0) - vbm, 4).tolist()))
        same = np.equal.outer(np.r_[[0] * 5, [1] * 3], np.r_[[0] * 5, [1] * 3])      # same-site blocks only
        rec["onsite_max_offdiag_meV"] = 1e3 * float(np.max(np.abs((H0 - np.diag(np.diag(H0)))[same])))
        # 2b. two-centre test of the 2NN Sc–Sc block along (1,1,0)/√2: yz,yz = (π+δ)/2 and yz,zx = (π−δ)/2 give π, δ;
        #     two-centre then predicts x²−y²,x²−y² = π and xy,xy = 3σ/4 + δ/4
        R2 = np.rint((np.array([A / 2, A / 2, 0]) - st.pos[0] + st.pos[0]) @ np.linalg.inv(st.cell)).astype(int)
        B = hm.H[np.where((hm.R == R2).all(1))[0][0]].real[:5, :5]
        pi_, de_ = B[1, 1] + B[1, 2], B[1, 1] - B[1, 2]
        rec["dd1_block_110"] = dict(xy_xy=B[0, 0], yz_yz=B[1, 1], yz_zx=B[1, 2], x2y2=B[3, 3], z2_z2=B[4, 4], xy_z2=B[0, 4],
                                    pi_from_t2g=pi_, delta_from_t2g=de_, x2y2_two_centre_pred=pi_)
        # 3. sign check: residual of the projection with each orbital flipped
        x, info = skfit.project_hr(model, hm)
        flips = {}
        for i, o in enumerate(MODEL_ORDER):
            sg = np.ones(8)
            sg[i] = -1
            hf = w90.regauge(hr, cell_w, cen_w, st.cell, cen_m, perm=PERM, signs=sg)
            flips[o] = 1e3 * skfit.project_hr(model, hf)[1]["rms_in_model"]
        rec["sign_check_rms_in_model_meV"] = {"none": 1e3 * info["rms_in_model"], **flips}
        rec["signs_ok"] = bool(1e3 * info["rms_in_model"] < min(flips.values()))
        # 4. projected parameters (on-site relative to the W90 VBM)
        p = dict(zip(model.names, x))
        for k in ("Ep", "Et2g", "Eeg"):
            p[k] -= vbm
        rec["projected_params"] = {k: round(float(v), 4) for k, v in p.items()}
        rec["projection"] = {"rms_in_model_meV": 1e3 * info["rms_in_model"],
                             "out_of_model_norm_eV": info["out_of_model"], "n_elements": info["n_elements"]}
        # 5. out-of-model content by bond length (Cartesian), Frobenius norm per shell
        R_cart = hm.R @ st.cell
        dist = np.linalg.norm(R_cart[:, None, None, :] + cen_m[None, None, :, :] - cen_m[None, :, None, :], axis=-1)
        used = np.zeros(hm.H.shape, bool)
        rmap = {tuple(r): i for i, r in enumerate(hm.R)}
        for r_i, R in enumerate(model.R):
            j = rmap.get(tuple(R))
            if j is not None:
                used[j] |= (np.abs(model.G[:, r_i]).sum(0) > 0) | (np.abs(model.C[r_i]) > 0)
        Hr = np.abs(hm.H.real)
        shells = {}
        for dval in np.unique(np.round(dist[dist < 2.01 * A], 3)):
            m = np.isclose(dist, dval, atol=2e-3)
            shells[f"{dval:.3f}"] = dict(in_model=float(np.sqrt(np.sum(Hr[m & used] ** 2))),
                                         dropped=float(np.sqrt(np.sum(Hr[m & ~used] ** 2))),
                                         max_dropped_meV=1e3 * float(np.max(np.r_[0, Hr[m & ~used]])))
        rec["by_distance_A"] = shells
        # 6. the projected SK model against the Wannier bands (VBM-referenced each)
        e_sk = model.eigvals(x, kp)
        e_sk -= model.eigvals(x, np.zeros((1, 3)))[0, 2]
        e_ref = e_m - vbm
        err = e_sk - e_ref
        cbm = e_ref[:, 3].min()
        ng = np.zeros_like(err, bool)
        ng[:, :3] = e_ref[:, :3] >= -1.0
        ng[:, 3:] = e_ref[:, 3:] <= cbm + 1.0
        P = sk.fcc_points(A)
        eG, eX = (model.eigvals(x, P[q][None])[0] for q in ("G", "X"))
        rec["projected_model_vs_w90"] = dict(
            rms_all_meV=1e3 * float(np.sqrt(np.mean(err ** 2))), max_all_meV=1e3 * float(np.max(np.abs(err))),
            rms_pm1eV_meV=1e3 * float(np.sqrt(np.mean(err[ng] ** 2))),
            gap_GX=float(eX[3] - eG[2]), gap_XX=float(eX[3] - eX[2]), gap_GG=float(eG[3] - eG[2]),
            m_l=float(model.mass(x, P["X"], [0, 0, 1], 3)), m_t=float(model.mass(x, P["X"], [1, 0, 0], 3)))
        out[lab] = rec
        print(lab, json.dumps({k: v for k, v in rec.items() if k != "by_distance_A"}, indent=1))
        print(lab, "by distance (Å): in-model / dropped norm (eV), max dropped (meV)")
        for dv, v in shells.items():
            print(f"   {dv:>7s}  {v['in_model']:.3f} / {v['dropped']:.3f}   {v['max_dropped_meV']:.1f}")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "sk_projection.json").write_text(json.dumps(out, indent=1))
    print("wrote", OUT / "sk_projection.json")


if __name__ == "__main__":
    main()
