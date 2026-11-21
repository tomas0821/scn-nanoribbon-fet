"""A1-2 analysis: Wannier90 quality (pre-registered criteria 4 and 5), HSE06 gaps and masses from the production
_hr.dat, and the bulk band figure.

Run from the project root, after copying w90/*/ back from the cluster:
    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P dft/A1/parse_a1_w90.py
Writes dft/A1/a1_w90.json and figures/A1_bulk_HSE_W90.png.
"""
import json
import re
from pathlib import Path

import numpy as np

from scnfet import w90
from scnfet.sk import HB2M

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
A_EXP = 4.501
A_RED = np.array([[-0.5, 0, 0.5], [0, 0.5, 0.5], [-0.5, 0.5, 0]])        # QE ibrav = 2, rows = a_i / a
A_ANG = A_RED * A_EXP
FROZ = {"pbe": 5.4, "hse": 6.5}
TWIN, PROD = "pbe", "hse_q8"                                         # frozen window tops, eV above the VBM
PATH = [("L", (0.5, 0.5, 0.5)), ("Γ", (0, 0, 0)), ("X", (0, 0, 1)), ("W", (0.5, 0, 1)), ("K", (0.75, 0, 0.75)), ("Γ", (0, 0, 0))]


def hk(hr, kcart):
    """H(k) for Cartesian k in 1/Å; R in lattice units, so k·R = 2π (A k / 2π)·R."""
    f = np.atleast_2d(kcart) @ A_ANG.T / (2 * np.pi)
    ph = np.exp(2j * np.pi * f @ hr.R.T) / hr.deg
    H = np.einsum("kr,rmn->kmn", ph, hr.H)
    return 0.5 * (H + np.conj(np.transpose(H, (0, 2, 1))))


def eig(hr, kcart):
    return np.linalg.eigvalsh(hk(hr, kcart))


def kpath(n=60):
    u = 2 * np.pi / A_EXP
    ks, ss, ticks, s0 = [], [], [(0.0, PATH[0][0])], 0.0
    for (l1, c1), (l2, c2) in zip(PATH[:-1], PATH[1:]):
        p1, p2 = u * np.array(c1, float), u * np.array(c2, float)
        t = np.linspace(0, 1, n, endpoint=False)         # each segment: its start included, its end left to the next
        L = np.linalg.norm(p2 - p1)
        ks.append(p1 + np.outer(t, p2 - p1)); ss.append(s0 + t * L)
        s0 += L
        ticks.append((s0, l2))
    ks.append((u * np.array(PATH[-1][1], float))[None]); ss.append([s0])
    return np.vstack(ks), np.concatenate(ss), ticks


def qe_bands(path):
    """Eigenvalues and Cartesian k (1/Å) from a QE bands/scf output (verbosity high)."""
    txt = Path(path).read_text()
    i = max(txt.rfind("End of self-consistent calculation"), txt.rfind("End of band structure calculation"))
    ks, es = [], []
    for m in re.finditer(r"k =\s*(-?[\d.]+)\s*(-?[\d.]+)\s*(-?[\d.]+) \(\s*\d+ PWs\)\s+bands \(ev\):\s+(.*?)(?=\n\s*\n\s*(?:occupation|k =|the Fermi|highest|Writing|$))", txt[i:], re.S):
        ks.append([float(m.group(j)) for j in (1, 2, 3)])
        es.append([float(v) for v in re.findall(r"-?\d+\.\d+", m.group(4))])
    return np.array(ks) * 2 * np.pi / A_EXP, np.array(es)


def centres(seed):
    """Final WF centres (Å) and spreads (Å²) from the .wout."""
    txt = (HERE / "w90" / seed / f"{seed}.wout").read_text()
    blk = txt[txt.rfind("Final State"):]
    c = [(float(a), float(b), float(d), float(sp)) for a, b, d, sp in
         re.findall(r"WF centre and spread\s+\d+\s+\(\s*(-?[\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+)\s*\)\s+([\d.]+)", blk)]
    om = re.search(r"Omega Total\s*=\s*([\d.]+)", blk)
    return np.array(c), (float(om.group(1)) if om else None)


def nearest_atom_dist(r):
    """Distance from r to the nearest Sc (fcc sites) or N (fcc + (a/2, 0, 0) equivalent) site."""
    best = 1e9
    rng = range(-2, 3)
    for n1 in rng:
        for n2 in rng:
            for n3 in rng:
                T = np.array([n1, n2, n3]) @ A_ANG
                for tau in (np.zeros(3), np.array([0.5, 0.5, 0.5]) @ A_ANG):
                    best = min(best, np.linalg.norm(r - T - tau))
    return best


def main():
    H = {s: w90.read_hr(HERE / "w90" / s / f"{s}_hr.dat") for s in ("pbe", "hse_n", "hse_o", "hse_q8")}
    # production: the 16³ Wannier grid when it exists (round-2 diagnostic: the 8³ grid limits the interpolation)
    twin = "pbe_k16" if (HERE / "w90" / "pbe_k16" / "pbe_k16_hr.dat").exists() else "pbe"
    prod = "hse_k16" if (HERE / "w90" / "hse_k16" / "hse_k16_hr.dat").exists() else "hse_q8"
    for s_ in {twin, prod}:
        H.setdefault(s_, w90.read_hr(HERE / "w90" / s_ / f"{s_}_hr.dat"))
    H["twin"], H["prod"] = H[twin], H[prod]
    out = {"twin_seed": twin, "production_seed": prod}
    # criterion 4: PBE twin against the QE PBE bands, inside the frozen window
    kq, eq = qe_bands(HERE / "outputs" / "pbe_scf_bands.out")
    _, enscf = qe_bands(HERE / "outputs" / "pbe_nscf.out")
    vbm_pbe = enscf[:, 7].max()
    ew = eig(H["twin"], kq)
    ref = eq[:, 5:13]                                         # QE bands 6–13 ↔ the 8 WFs
    inwin = ref <= vbm_pbe + FROZ["pbe"]
    err = np.abs(ew - ref)[inwin]
    out["criterion4"] = dict(max_meV=1e3 * float(err.max()), rms_meV=1e3 * float(np.sqrt(np.mean(err ** 2))), n=int(err.size),
                             ok=bool(err.max() <= 0.010 and np.sqrt(np.mean(err ** 2)) <= 0.003))
    # criterion 5: routes O and N give the same H(R)
    same_R = H["hse_n"].R.shape == H["hse_o"].R.shape and np.array_equal(H["hse_n"].R, H["hse_o"].R)
    dH = float(np.max(np.abs(H["hse_n"].H - H["hse_o"].H))) if same_R else float("nan")
    dq = float(np.max(np.abs(H["hse_q8"].H - H["hse_o"].H))) if np.array_equal(H["hse_q8"].R, H["hse_o"].R) else float("nan")
    out["criterion5"] = dict(same_R=bool(same_R), max_dH_meV=1e3 * dH, ok=bool(same_R and dH <= 0.001),
                             q8_vs_q4_max_dH_meV=1e3 * dq)
    # centres and spreads
    cs = {}
    for s in [x for x in H if x not in ("twin", "prod")]:
        c, om = centres(s)
        d = [nearest_atom_dist(r[:3]) for r in c]
        cs[s] = dict(max_offset_A=float(max(d)), spreads_A2=[round(float(x), 4) for x in c[:, 3]], omega_total=om)
    out["centres"] = cs
    out["criterion4"]["centres_ok"] = bool(all(v["max_offset_A"] <= 0.05 for v in cs.values()))
    # HSE06 production: gaps and masses from hse_q8
    u = 2 * np.pi / A_EXP
    G, X = np.zeros(3), u * np.array([0, 0, 1.0])
    hq = H["prod"]
    eG, eX = eig(hq, G)[0], eig(hq, X)[0]
    vbm = eG[2]

    def mass(k0, d, band, h=0.004):
        d = np.asarray(d, float) / np.linalg.norm(d)
        e = eig(hq, np.array([k0 - h * d, k0, k0 + h * d]))[:, band]
        return HB2M / ((e[0] + e[2] - 2 * e[1]) / h ** 2)
    hole = {f"{lab}_band{b}": -mass(G, dirn, b) for lab, dirn in (("100", [1, 0, 0]), ("110", [1, 1, 0]), ("111", [1, 1, 1])) for b in (2, 1, 0)}
    out["production"] = dict(gap_GX=float(eX[3] - vbm), gap_XX=float(eX[3] - eX[2]), gap_GG=float(eG[3] - vbm),
                         G12=float(eG[6] - vbm), m_l=float(mass(X, [0, 0, 1], 3)), m_t=float(mass(X, [1, 0, 0], 3)),
                         hole_masses={k: float(v) for k, v in hole.items()})
    for s in ("hse_o", "pbe"):
        hh = H[s]
        eG, eX = eig(hh, G)[0], eig(hh, X)[0]
        out[s + "_gaps"] = dict(gap_GX=float(eX[3] - eG[2]), gap_XX=float(eX[3] - eX[2]), gap_GG=float(eG[3] - eG[2]))
    (HERE / "a1_w90.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
    global TWIN, PROD
    TWIN, PROD = twin, prod
    plot(H, kq, eq, vbm_pbe)


def plot(H, kq, eq, vbm_pbe):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scnfet import targets_lit as tl
    ink, ink2, grid, ref = "#17181a", "#4d4f53", "#e6e6e4", "#b8b7b2"
    c1, c2 = "#2a78d6", "#eb6834"
    plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": ink2, "xtick.color": ink2, "ytick.color": ink2, "axes.linewidth": 0.8})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), sharey=True, constrained_layout=True)
    kp, s, ticks = kpath(90)

    def frame(ax, title):
        for t, _ in ticks:
            ax.axvline(t, color=grid, lw=0.8, zorder=0)
        ax.axhline(0, color=grid, lw=0.8, zorder=0)
        ax.set_xticks([t for t, _ in ticks], [l for _, l in ticks])
        ax.set_xlim(s[0], s[-1]); ax.set_ylim(-6.5, 10.5)
        ax.set_title(title, loc="left", fontsize=9.5, color=ink)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    # (a) PBE twin: Wannier vs QE bands on the same path (QE path in the same segment order)
    ax = axes[0]
    e = eig(H["twin"], kp) - vbm_pbe
    for b in range(8):
        ax.plot(s, e[:, b], color=c1, lw=2, zorder=3, label="Wannier90 (8 WFs)" if b == 0 else None)
    # QE path coordinate: cumulative distance along kq
    sq = np.r_[0, np.cumsum(np.linalg.norm(np.diff(kq, axis=0), axis=1))]
    for b in range(5, 16):
        ax.plot(sq, eq[:, b] - vbm_pbe, ".", color=ink, ms=2.2, zorder=4, label="QE PBE bands" if b == 5 else None)
    ax.axhline(FROZ["pbe"], color=ink2, lw=0.8, ls=(0, (3, 3)), label=f"frozen window top (+{FROZ['pbe']} eV)")
    frame(ax, f"(a) PBE twin ({'16³' if 'k16' in TWIN else '8³'} Wannier grid), a = 4.501 Å")
    ax.set_ylabel("E − E_VBM (eV)")
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")

    # (b) HSE06 production (q 8³) with Deng's HSE06 bands (a = 4.50 Å) on the same path coordinate
    ax = axes[1]
    G = np.zeros(3)
    vb = eig(H["prod"], G[None])[0, 2]
    e = eig(H["prod"], kp) - vb
    xpt = [95.8, 141.5, 141.5, 194.3, 194.3, 212.9, 212.9, 268.9]
    # Deng's path is L–Γ–X–U|K–Γ; ours is L–Γ–X–W–K–Γ, so overlay Deng on the shared L–Γ–X segments only
    tk = [t for t, _ in ticks]
    first = True
    for px, pe in tl.deng_paths():
        m = px <= 194.3
        if m.sum() < 2:
            continue
        sx = np.interp(px[m], [95.8, 141.5, 194.3], [tk[0], tk[1], tk[2]])
        ax.plot(sx, pe[m], color=ref, lw=5, alpha=0.8, zorder=2, solid_capstyle="round", label="Deng 2015 HSE06 (L–Γ–X)" if first else None)
        first = False
    for b in range(8):
        ax.plot(s, e[:, b], color=c2, lw=2, zorder=3, label=f"ours: HSE06 q 8³, {'16³' if 'k16' in PROD else '8³'} Wannier grid" if b == 0 else None)
    ax.axhline(FROZ["hse"], color=ink2, lw=0.8, ls=(0, (3, 3)), label=f"frozen window top (+{FROZ['hse']} eV)")
    frame(ax, "(b) HSE06 (production, q 8³), a = 4.501 Å")
    ax.legend(frameon=False, fontsize=7.5, loc="lower center", bbox_to_anchor=(0.4, 0.0))
    out = ROOT / "figures" / "A1_bulk_HSE_W90.png"
    fig.savefig(out, dpi=200, facecolor="white")
    print("wrote", out)


if __name__ == "__main__":
    main()
