"""Literature HSE06 targets for bulk rocksalt ScN: debugging targets for the SK fitter scaffold (plan §4, T0).

These are cross-checks, not production fit targets. The production model is fitted to our own semicore-consistent
NC-HSE06 (plan D3; LM7: the literature HSE sets differ by a rigid ~0.13 eV gap shift from the semicore treatment,
and their m_t differ by 40 %). Energies in eV with the bulk VBM (Γ15) = 0, masses in m_e.

Sources:
  * Deng et al., PRB 91, 045104 (2015), HSE06 at a = 4.50 Å. Fig. 1(a) is vector graphics in the author PDF; its 11
    band paths were extracted exactly into `data/deng2015_fig1a_paths.json` (PDF points). Calibration: path
    L–Γ–X–U|K–Γ at x = 95.8, 141.5, 194.3, 212.9, 268.9 pt; E = (122.7 − y)/4.192 eV. This reproduces the
    tabulated digitised values (Γ25′ 3.73, Γ12 4.91, X3 0.94, X5′ −1.09, L3′ −2.93 …) to 0.01 eV.
    Masses from the text: m_l 1.44 ± 0.02, m_t 0.24 ± 0.04.
  * Mu et al., PRB 104, 075118 (2021), HSE06 with the Sc 3s3p semicore frozen, a = 4.48 Å. High-symmetry energies
    digitised from the Fig. 1 raster (±0.02 eV); masses from Table I: m_l 1.43, m_t 0.17.
"""
import json
from pathlib import Path

import numpy as np

from scnfet import sk
from scnfet.skfit import CharTarget, MassTarget, Targets

DATA = Path(__file__).parent / "data"

# Deng Fig. 1(a) frame calibration (PDF points)
_DENG_X = {"L": 95.8, "G": 141.5, "X": 194.3, "U": 212.9, "K": 212.9, "G2": 268.9}
_DENG_SEGMENTS = [("L", "G"), ("G", "X"), ("X", "U"), ("K", "G2")]


def character_targets(a):
    """Rocksalt-ScN band characters that energies alone cannot fix: the CBM at X_z is X3 = d_xy,
    and the lowest CB at Γ is Γ25′ = t2g."""
    P = sk.fcc_points(a)
    return [CharTarget("X3 = d_xy", P["X"], 3, (":dxy",), 0.9),
            CharTarget("Γ25′ = t2g", P["G"], 3, (":dxy", ":dyz", ":dzx"), 0.9)]


def _deng_E(y):
    return (122.7 - np.asarray(y)) / ((185.5 - 80.7) / 25)


def deng_paths():
    """The 11 digitised curves as (x_pt, E_eV) arrays, sorted by x, duplicate points removed."""
    raw = json.loads((DATA / "deng2015_fig1a_paths.json").read_text())
    out = []
    for p in raw:
        p = np.unique(np.asarray(p, float), axis=0)
        p = p[np.argsort(p[:, 0], kind="stable")]
        out.append((p[:, 0], _deng_E(p[:, 1])))
    return out


def deng_k_of_x(xpt, a=4.50):
    """Map a Deng-figure x coordinate to a Cartesian k on L–Γ–X–U|K–Γ (X = X_z)."""
    P = sk.fcc_points(a)
    P = {**P, "G2": P["G"]}
    for s0, s1 in _DENG_SEGMENTS:
        x0, x1 = _DENG_X[s0], _DENG_X[s1]
        if x0 <= xpt <= x1 and x1 > x0:
            t = (xpt - x0) / (x1 - x0)
            return P[s0] + t * (P[s1] - P[s0])
    raise ValueError(f"x = {xpt} outside the Deng path")


def deng_path_targets(a=4.50, window=(-8.0, 5.5), step=1.5, sigma_vb=0.05, sigma_cb=0.05, sigma_hi=0.10,
                      e_hi=4.0, with_masses=True):
    """Dense band targets from Deng Fig. 1(a): at each sampled x, the curve energies inside `window`, sorted, are
    bands 0, 1, 2 (VB: E ≤ 0.2) and 3, 4, … (CB). The window keeps only the 8 model bands: it drops N 2s
    (< −13 eV) and the Sc 4s-mixed bands (Γ1 7.07, X 6.14+, L 6.37+). Points where the VB count is not 3
    (a gap in the digitisation) are skipped. CB energies above `e_hi` get `sigma_hi`."""
    paths = deng_paths()
    xs = []
    for s0, s1 in _DENG_SEGMENTS:
        x0, x1 = _DENG_X[s0], _DENG_X[s1]
        n = max(int((x1 - x0) / step), 2)
        xs += list(np.linspace(x0 + 0.05, x1 - 0.05, n))
    ks, ik, band, E, sig, grp, flo = [], [], [], [], [], [], []
    skipped = 0
    for xpt in xs:
        vals = []
        for px, pe in paths:
            if px[0] <= xpt <= px[-1]:
                vals.append(float(np.interp(xpt, px, pe)))
        vals = np.sort([v for v in vals if window[0] <= v <= window[1]])
        vb, cb = vals[vals <= 0.2], vals[vals > 0.2]
        if len(vb) != 3 or len(cb) > 5:
            skipped += 1
            continue
        i = len(ks)
        ks.append(deng_k_of_x(xpt, a))
        for b, v in enumerate(np.concatenate([vb, cb])):
            ik.append(i)
            band.append(b)
            E.append(v)
            sig.append(sigma_vb if b < 3 else (sigma_cb if v <= e_hi else sigma_hi))
            grp.append("VB" if b < 3 else "CB")
            flo.append(False)
        # No extra bands inside the window: the next model band must lie above the window top. Without this the
        # fit is free to drop an unconstrained band into the gap region (it put a second CB at X at 1.4–1.8 eV).
        if 3 + len(cb) <= 7:
            ik.append(i)
            band.append(3 + len(cb))
            E.append(window[1])
            sig.append(sigma_cb)
            grp.append("floor")
            flo.append(True)
    P = sk.fcc_points(a)
    masses = []
    if with_masses:
        masses = [MassTarget("m_l (X→Γ)", P["X"], np.array([0, 0, 1.0]), 3, 1.44, 0.03),
                  MassTarget("m_t (X→W)", P["X"], np.array([1.0, 0, 0]), 3, 0.24, 0.10)]
    T = Targets("Deng 2015 HSE06 (Fig. 1a, vector)", np.array(ks), np.array(ik), np.array(band), np.array(E),
                np.array(sig), masses=masses, ref=(P["G"], 2), group=np.array(grp), floor=np.array(flo),
                chars=character_targets(a))
    T.skipped = skipped
    return T


# Mu 2021 Fig. 1 high-symmetry energies (B_targets.md table A), and Deng's for the same points
MU_TABLE = {
    "G": [(0, 0.0), (1, 0.0), (2, 0.0), (3, 3.58), (4, 3.58), (5, 3.58), (6, 4.74), (7, 4.74)],
    "X": [(0, -3.61), (1, -1.12), (2, -1.12), (3, 0.79)],
    "L": [(0, -5.70), (1, -3.02), (2, -3.02)],
    "K": [(0, -3.44), (1, -3.41), (2, -2.13), (3, 3.38)],
}
DENG_TABLE = {
    "G": [(0, 0.0), (1, 0.0), (2, 0.0), (3, 3.73), (4, 3.73), (5, 3.73), (6, 4.91), (7, 4.91)],
    "X": [(0, -3.52), (1, -1.09), (2, -1.09), (3, 0.94)],
    "L": [(0, -5.61), (1, -2.93), (2, -2.93)],
    "K": [(0, -3.35), (1, -3.33), (2, -2.06), (3, 3.43)],
}


def table_targets(table, label, a, ml, mt, sigma=0.03, ml_rel=0.03, mt_rel=0.10, floor_E=5.5):
    """High-symmetry energies + masses. At every listed point the next band above the listed ones gets a floor at
    `floor_E` (both HSE sets have no band between their listed CBs and ≈ 5.9 eV at X, L, K). Point tables alone
    leave the bands free between the points; use them as cross-checks, not as a fit's only data."""
    P = sk.fcc_points(a)
    names = list(table)
    k = np.array([P[n] for n in names])
    ik, band, E, grp, flo = [], [], [], [], []
    for i, n in enumerate(names):
        for b, v in table[n]:
            ik.append(i)
            band.append(b)
            E.append(v)
            grp.append("VB" if b < 3 else "CB")
            flo.append(False)
        nxt = max(b for b, _ in table[n]) + 1
        if nxt <= 7:
            ik.append(i)
            band.append(nxt)
            E.append(floor_E)
            grp.append("floor")
            flo.append(True)
    masses = [MassTarget("m_l (X→Γ)", P["X"], np.array([0, 0, 1.0]), 3, ml, ml_rel),
              MassTarget("m_t (X→W)", P["X"], np.array([1.0, 0, 0]), 3, mt, mt_rel)]
    return Targets(label, k, np.array(ik), np.array(band), np.array(E), np.full(len(E), sigma), masses=masses,
                   ref=(P["G"], 2), group=np.array(grp), floor=np.array(flo), chars=character_targets(a))


def mu_targets():
    return table_targets(MU_TABLE, "Mu 2021 HSE06 (semicore frozen)", 4.48, ml=1.43, mt=0.17)
