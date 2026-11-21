"""A1 stage-1 analysis: gaps, the pre-registered criteria 1–3 (dft/A1/README.md) and the Wannier window survey.

Run from the project root:  ~/.venvs/scnfet/bin/python -P dft/A1/parse_a1.py
Writes dft/A1/a1_stage1.json.
"""
import itertools
import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
A = np.array([[-0.5, 0, 0.5], [0, 0.5, 0.5], [-0.5, 0.5, 0]])        # QE ibrav = 2 primitive vectors / a
NOCC = 8                                                              # 16 valence electrons


def parse(tag):
    txt = (OUT / f"{tag}.out").read_text()
    assert "JOB DONE" in txt and "convergence NOT achieved" not in txt, tag
    i = max(txt.rfind("End of self-consistent calculation"), txt.rfind("End of band structure calculation"))
    blk = txt[i:]
    ks, eig = [], []
    for m in re.finditer(r"k =\s*(-?[\d.]+)\s*(-?[\d.]+)\s*(-?[\d.]+) \(\s*\d+ PWs\)\s+bands \(ev\):\s+(.*?)(?=\n\s*\n\s*(?:occupation|k =|the Fermi|highest|Writing|$))", blk, re.S):
        ks.append([float(m.group(j)) for j in (1, 2, 3)])
        eig.append([float(v) for v in re.findall(r"-?\d+\.\d+", m.group(4))])
    r = dict(k=np.array(ks), eig=np.array(eig))
    m = re.search(r"highest occupied, lowest unoccupied level \(ev\):\s+(-?[\d.]+)\s+(-?[\d.]+)", txt)
    if m:
        r["homo"], r["lumo"] = float(m.group(1)), float(m.group(2))
    m = re.findall(r"the Fermi energy is\s+(-?[\d.]+)", txt)
    if m:
        r["ef"] = float(m[-1])
    t = re.search(r"PWSCF\s*:.*?([\d.hms ]+)WALL", txt)
    r["wall"] = t.group(1).strip() if t else None
    return r


def special(k):
    """Indices of Γ and of the X points (|k| = 1 along a Cartesian axis, 2π/a units)."""
    iG = np.where(np.linalg.norm(k, axis=1) < 1e-6)[0]
    iX = np.where(np.isclose(np.sort(np.abs(k), axis=1), [0, 0, 1], atol=1e-5).all(axis=1))[0]
    return iG, iX


def gaps(r):
    e, k = r["eig"], r["k"]
    iG, iX = special(k)
    G, X = e[iG[0]], e[iX[0]]
    return dict(GX=X[NOCC] - G[NOCC - 1], XX=X[NOCC] - X[NOCC - 1], GG=G[NOCC] - G[NOCC - 1],
                vbm=e[:, NOCC - 1].max(), cbm=e[:, NOCC].min())


def equivalent(k1, k2, ops):
    """True if k2 = R k1 + G for some cubic operation R (48) and reciprocal-lattice vector G."""
    for R in ops:
        d = R @ k1 - k2
        f = A @ d                     # crystal coordinates of d: integers ⇔ d is a reciprocal-lattice vector
        if np.allclose(f, np.round(f), atol=1e-5):
            return True
    return False


def cubic_ops():
    ops = []
    for perm in itertools.permutations(range(3)):
        for sg in itertools.product([1, -1], repeat=3):
            M = np.zeros((3, 3))
            for i, p in enumerate(perm):
                M[i, p] = sg[i]
            ops.append(M)
    return ops


def main():
    R = {t: parse(t) for t in ("pbe_scf", "pbe_scf_bands", "pbe_nscf", "hse_ibz", "hse_nosym", "hse_q8", "hse_fock320")}
    g = {t: gaps(R[t]) for t in ("pbe_nscf", "hse_ibz", "hse_nosym", "hse_q8", "hse_fock320")}
    # criterion 1: route N (full grid) against route O (= IBZ eigenvalues, symmetry-unfolded)
    ops = cubic_ops()
    ibz, full = R["hse_ibz"], R["hse_nosym"]
    dmax = 0.0
    nomatch = 0
    for kf, ef in zip(full["k"], full["eig"]):
        j = next((j for j, ki in enumerate(ibz["k"]) if equivalent(ki, kf, ops)), None)
        if j is None:
            nomatch += 1
            continue
        dmax = max(dmax, float(np.max(np.abs(ef[:20] - ibz["eig"][j][:20]))))
    c1 = dict(max_meV=1e3 * dmax, unmatched=nomatch, n_full=len(full["k"]), n_ibz=len(ibz["k"]), ok=dmax <= 0.002 and nomatch == 0)
    ref = g["hse_ibz"]
    dq = {q: 1e3 * (g["hse_q8"][q] - ref[q]) for q in ("GX", "XX", "GG")}
    df = {q: 1e3 * (g["hse_fock320"][q] - ref[q]) for q in ("GX", "XX", "GG")}
    c2 = dict(delta_meV=dq, ok=all(abs(v) <= 10 for v in dq.values()))
    c3 = dict(delta_meV=df, ok=all(abs(v) <= 10 for v in df.values()))
    # Wannier window survey on the full grids (band index 1-based in the text)
    survey = {}
    for t in ("pbe_nscf", "hse_nosym"):
        e = R[t]["eig"]
        ref_e = g[t]["vbm"]
        survey[t] = dict(vbm=ref_e, band6_min=e[:, 5].min() - ref_e, band13_max=e[:, 12].max() - ref_e,
                         band14_min=e[:, 13].min() - ref_e, band14_max=e[:, 13].max() - ref_e,
                         band30_min=e[:, 29].min() - ref_e, nbnd=e.shape[1],
                         per_band_min_max=[[float(e[:, b].min() - ref_e), float(e[:, b].max() - ref_e)] for b in range(5, e.shape[1])])
    out = dict(gaps=g, criterion1_routes=c1, criterion2_q=c2, criterion3_ecutfock=c3, survey=survey,
               wall={t: R[t]["wall"] for t in R}, pbe_ef=R["pbe_scf"].get("ef"))
    (HERE / "a1_stage1.json").write_text(json.dumps(out, indent=1, default=float))
    print("gaps (eV):")
    for t, v in g.items():
        print(f"  {t:12s} Γ–X {v['GX']:+.4f}  X–X {v['XX']:.4f}  Γ–Γ {v['GG']:.4f}")
    print("criterion 1 (routes O vs N, bands 1–20, ≤ 2 meV):", c1)
    print("criterion 2 (q 4³ → 8³, ≤ 10 meV):", {k: round(v, 1) for k, v in dq.items()}, "PASS" if c2["ok"] else "FAIL")
    print("criterion 3 (ecutfock 160 → 320, ≤ 10 meV):", {k: round(v, 2) for k, v in df.items()}, "PASS" if c3["ok"] else "FAIL")
    for t, s in survey.items():
        print(f"window survey {t} (eV rel. VBM): band6 min {s['band6_min']:.2f}, band13 max {s['band13_max']:.2f}, "
              f"band14 min {s['band14_min']:.2f} / max {s['band14_max']:.2f}, band30 min {s['band30_min']:.2f}")
        print("   bands 6–20 [min, max]:", [[round(a, 2), round(b, 2)] for a, b in s["per_band_min_max"][:15]])
    print("wall:", out["wall"])


if __name__ == "__main__":
    main()
