"""Cost model for A4a (Gate M), A4b and A5, and the LR1 DFT-Wannier ribbon options (2026-11-11; updated 2026-11-12
with the A4a-1 actuals and the HSE cost tests, dft/A2/costtest/).

PBE: the heavy-plan cost law (ln c = −18.270 + 1.288 ln N_PW + 1.436 ln nbnd, core-s
per (k, spin, SCF iteration)), fitted on in-house runs of ribbon size. Our A0/A2 sheet runs are too small to refit it
(seconds of wall, overhead-dominated); they sit ×1.3–3.2 above it, so PBE is quoted as law × [1, 2].

HSE: fitted on our 13 A1/A2 HSE06 runs (qe_timing.py). The QE 7.5 EXX cost scales as nbnd² (all bands, empty ones
included), not nbnd·nocc, **with fixed occupations and with gaussian smearing alike** (cost test T1: 24 → 40 bands
costs ×2.90 fixed, ×2.91 smearing; quadratic predicts 2.78):
    core-s per (k_IBZ, q, outer step) = C · nbnd² · N_r ln N_r,   C = 2.8–5.9 ns (mean 4.5)
The upper end comes from the cost tests: the same input ran its EXX part ×1.68 slower on one compute node (12 Nov) than on
another (11 Nov), with the non-EXX part identical; nodes are the same hardware. N_r is the EXX FFT grid at ecutfock 160.
Outer EXX steps 6–7; + 10 % for the non-EXX part.

W90 band count (cost test T2): the sheet's Wannier bands pass with the outer window down to VBM + 11.05 eV (HSE),
i.e. 20 bands per f.u.; 17–18 per f.u. fail badly. Ribbons carry ≈ 1.5× more vacuum per f.u., so more free-electron
bands fall below the same level: 20–24 per f.u. here, fixed per ribbon by the A4a-2 PBE band survey.

    ~/.venvs/scnfet/bin/python -I dft/A4a/cost_a4.py        (from the project root)
Writes dft/A4a/cost_a4.json.
"""
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from scnfet import structures as S          # noqa: E402

BOHR = 0.529177
H = 3600.0
A_P = 3.117899 * 4.501 / 4.5134            # A2 convention: A0 PBE sheet scaled by bulk a_exp / a_PBE
D_PLANE = 2.1685                             # A2 scaled plane spacing (Sc/N 2.168/2.169 Å)
VAC = 15.0                                   # structures.ribbon default, in y and z
ZV = dict(Sc=11, N=5, H=1)                   # valence: Sc spn PAW / SG15 both 11
CO, PBE_FACTOR = (-18.270, 1.288, 1.436), (1.0, 2.0)
C_EXX = (2.8e-9, 5.9e-9)                     # s per (k, q, outer) per nbnd² per N_r ln N_r (13 runs, 2 nodes)
A4A1_ACTUAL = {"[100] n=3": 190.2, "[110] n=4": 109.8, "[100] n=3 H": 279.7}   # jobs 142068–142070, 64 × wall
# jobs 142490–142492 (12 Nov), 64 × sacct Elapsed; the vacuum check ran inside the [110] n=4 job
A4A2_ACTUAL = {"[100] n=3": 64 * (3 + 7 / 60 + 15 / 3600), "[110] n=4": 64 * (1 + 30 / 60 + 39 / 3600),
               "[100] n=3 H": 64 * 57 / 60}
# job 142557 (12 Nov), 64 × 37 min 04 s: NM SCF at 20 Å on [100] n=3 bare and H (decision after Gate M)
VAC20_ACTUAL = 64 * (37 / 60 + 4 / 3600)
# job 142699 (13 Nov), 64 × 2 h 20 min 47 s: relax_r110_4h ([110] n=4 H, 20 Å; 31 BFGS steps)
H110_RELAX_ACTUAL = 64 * (2 + 20 / 60 + 47 / 3600)
# job 142745 (13 Nov), 64 × 53 min 35 s: r110_4h NM, FM start, bands (nbnd 316), projwfc
H110_GATEM_ACTUAL = 64 * (53 / 60 + 35 / 3600)
# r110_4h band survey at 20 Å (a4a_stage2.json): smallest N with band minimum ≥ VBM + 10.28 / + 11.28 eV (PBE)
H110_NBND_MEASURED = (201, 217)
# job 142949 (13 Nov), 64 × 44 min 16 s: NM SCF at 25 Å on r100_3h and r110_4
VAC25_ACTUAL = 64 * (44 / 60 + 16 / 3600)
# A4b-1 job 1 (13 Nov, jobs 142945–142948), 64 × sacct Elapsed: NC PBE twin scf → pot → bands → projwfc → nscf →
# w90 -pp → pw2wannier90 (B, S). It replaces the PBE-W90 (nscf + p2w) the budget carried for these four ribbons.
A4B1_ACTUAL = {"r100_3": 64 * (31 / 60 + 7 / 3600), "r100_3h": 64 * (54 / 60 + 59 / 3600),
               "r110_4": 64 * (30 / 60 + 40 / 3600), "r110_4h": 64 * (30 / 60 + 1 / 3600)}
A4B1_W90_ACTUAL = 64 * (17 / 60 + 8 / 3600)   # job 143002 (13 Nov), the 12 Wannier90 variants, 17 min 08 s
# A4b-2 (13–14 Nov), 64 × sacct Elapsed: HSE06 job 1 per ribbon (jobs 143006–143009), then the Wannier90 jobs
# (143231 r100_3 + r110_4h, 23 min 19 s; 143234 r110_4, 18 min 53 s; 143311 r100_3h, 16 min 02 s) and the sheet projwfc
# (143232, 1 min)
A4B2_ACTUAL = {"r100_3": 64 * (7 + 6 / 60 + 5 / 3600), "r110_4h": 64 * (13 + 27 / 60 + 12 / 3600),
               "r110_4": 64 * (16 + 50 / 60 + 3 / 3600), "r100_3h": 64 * (21 + 15 / 60 + 57 / 3600)}
A4B2_W90_ACTUAL = 64 * ((23 / 60 + 19 / 3600) + (18 / 60 + 53 / 3600) + (16 / 60 + 2 / 3600) + 1 / 60)
SPENT = dict(A0=9.7, A1=145, A2=88, A4a1=sum(A4A1_ACTUAL.values()), cost_tests=18.2, A4a2=sum(A4A2_ACTUAL.values()),
             vac20=VAC20_ACTUAL, h110_relax=H110_RELAX_ACTUAL, h110_gateM=H110_GATEM_ACTUAL,
             vac25=VAC25_ACTUAL, A4b1=sum(A4B1_ACTUAL.values()) + A4B1_W90_ACTUAL,
             A4b2=sum(A4B2_ACTUAL.values()) + A4B2_W90_ACTUAL)
# A5-1 relaxations (jobs 143316 r110_8 3 h 37 min 35 s, 143317 r100_6 4 h 08 min 08 s; 14–15 Nov), 64 × wall
A5_RELAX_ACTUAL = {"r110_8": 64 * (3 + 37 / 60 + 35 / 3600), "r100_6": 64 * (4 + 8 / 60 + 8 / 3600)}
SPENT["A5_relax"] = sum(A5_RELAX_ACTUAL.values())
# A5-1b estimate from measured runs: each relax's first (from-scratch) NM SCF, CPU s to convergence, × 1.07 wall/CPU
# (A4a-2), times the A4a-2 bare/H ribbons' wall ratios to their NM SCF: FM 1.73–1.81, AFM edges 1.67–1.70,
# bands (28 bands per f.u.) 2.13–2.74, projwfc 0.07–0.09; × 1.3 on top for SCF iteration counts in the wider ribbons
A5_NM_SCF_CPU_S = {"r110_8": 1648.7, "r100_6": 2451.5}
A5_STAGE2_RATIO = (1 + 1.73 + 1.67 + 2.13 + 0.07, (1 + 1.81 + 1.70 + 2.74 + 0.09) * 1.3)
# A5-1b actual (jobs 143514 r110_8 3 h 47 min 22 s, 143515 r100_6 5 h 26 min 02 s; 15 Nov), 64 × wall
A5_1B_ACTUAL = {"r110_8": 64 * (3 + 47 / 60 + 22 / 3600), "r100_6": 64 * (5 + 26 / 60 + 2 / 3600)}
SPENT["A5_1b"] = sum(A5_1B_ACTUAL.values())
# A3-1, the slab relaxation (job 143638, 15 Nov, 1 min 54 s), 64 × wall
A3_RELAX_ACTUAL = 64 * (1 / 60 + 54 / 3600)
SPENT["A3_relax"] = A3_RELAX_ACTUAL
# A5-2 job 1, the NC PBE twins (jobs 143633 r110_8 1 h 36 min 35 s, 143634 r100_6 1 h 59 min 48 s; 15 Nov), and the
# A3-2 PBE twin + vacuum check (143662, 8 min 34 s; 143664, 13 s), 64 × wall
A5_2_JOB1_ACTUAL = {"r110_8": 64 * (1 + 36 / 60 + 35 / 3600), "r100_6": 64 * (1 + 59 / 60 + 48 / 3600)}
SPENT["A5_2_job1"] = sum(A5_2_JOB1_ACTUAL.values())
A3_2_PBE_ACTUAL = 64 * (8 / 60 + 34 / 3600) + 64 * 13 / 3600
SPENT["A3_2_pbe"] = A3_2_PBE_ACTUAL
# A5-2 job 2 + the slab's PBE Wannier90 (143673, 1 h 09 min 30 s), the A3-2 HSE06 run (143663, 3 h 17 min 18 s),
# A5-3 r110_8 HSE gap point (143636, two nodes, 8 h 59 min 59 s); 64 × wall per node. A5-3 r100_6 (143637) still running.
A5_2_JOB2_ACTUAL = 64 * (1 + 9 / 60 + 30 / 3600)
SPENT["A5_2_job2"] = A5_2_JOB2_ACTUAL
A3_2_HSE_ACTUAL = 64 * (3 + 17 / 60 + 18 / 3600)
SPENT["A3_2_hse"] = A3_2_HSE_ACTUAL
A5_3_ACTUAL = {"r110_8": 128 * (8 + 59 / 60 + 59 / 3600)}
SPENT["A5_3_r110_8"] = A5_3_ACTUAL["r110_8"]
# A3-3 HSE Wannier90 of the slab (143764, 16 Nov, 2 min 14 s), 64 × wall
A3_3_HW90_ACTUAL = 64 * (2 / 60 + 14 / 3600)
SPENT["A3_3_hw90"] = A3_3_HW90_ACTUAL
# A5-3 r100_6 HSE gap point (143637, 15–16 Nov, one node, 20 h 14 min 20 s), 64 × wall
A5_3_ACTUAL["r100_6"] = 64 * (20 + 14 / 60 + 20 / 3600)
SPENT["A5_3_r100_6"] = A5_3_ACTUAL["r100_6"]


def spent_upto(key):
    """Core-h spent through SPENT[key], inclusive (SPENT is in run order): for blocks that record a past state."""
    keys = list(SPENT)
    return sum(SPENT[k] for k in keys[:keys.index(key) + 1])


# EXX constant recalibrated on the four A4b-2 ribbons (core-s per outer step per IBZ k per q per nbnd² N_r ln N_r):
# 15 Å [100] bare 6.0e-9; the three 20 Å runs 1.05–1.06e-8 on three nodes (the 13-run law: 2.8–5.9e-9)
C_EXX_RIBBON = (6.0e-9, 1.1e-8)
OUTER = (6, 7)
NONEXX = 1.10


def c_kit(pw, nb):
    return math.exp(CO[0] + CO[1] * math.log(pw) + CO[2] * math.log(nb))


def good_fft(x):
    n = int(math.ceil(x))
    while True:
        m = n
        for p in (2, 3, 5):
            while m % p == 0:
                m //= p
        if m == 1:
            return n
        n += 1


def ribbon(direction, n, h=False, vacuum=None):
    st = S.ribbon(direction, n, a_p=A_P, d=D_PLANE, vacuum=vacuum or VAC, passivate=h)
    sp = list(st.species)
    L, Wy, Wz = st.cell[0, 0], st.cell[1, 1], st.cell[2, 2]
    nSc, nN, nH = sp.count("Sc"), sp.count("N"), sp.count("H")
    nel = sum(ZV[s] for s in sp)
    return dict(label=f"[{direction}] n={n}{' H' if h else ''}", direction=direction, n=n, h=h, nat=len(sp), nSc=nSc,
                nN=nN, nH=nH, nel=nel, nocc=nel // 2, L=L, Wy=Wy, Wz=Wz, V=L * Wy * Wz,
                W_nm=S.width(st) / 10)


def npw(V, L, ecut, nk_line):
    kmax = (nk_line // 2) / nk_line * 2 * math.pi / (L / BOHR)
    return (V / BOHR ** 3) * (math.sqrt(ecut) + kmax) ** 3 / (6 * math.pi ** 2)


def pbe(r, kind, nk, iters, nspin=1, nbnd=None, ecut=90.0):
    """core-h range for one PBE run: nk irreducible (or explicit) k, iters SCF iterations (NSCF: 18 equivalents)."""
    pw = npw(r["V"], r["L"], ecut, 15)
    nb = nbnd or int(math.ceil(1.2 * r["nel"] / 2))
    pen = 1.20 if (nspin == 2 and kind == "relax") else (1.08 if nspin == 2 else 1.0)
    c = c_kit(pw, nb) * nk * nspin * iters * pen / H
    return [c * f for f in PBE_FACTOR]


def hse(r, nk_full, nq, nbnd):
    g = [good_fft(x / BOHR * math.sqrt(160.0) / math.pi) for x in (r["L"], r["Wy"], r["Wz"])]
    Nr = g[0] * g[1] * g[2]
    nk_ibz = nk_full // 2 + 1
    per = nbnd ** 2 * Nr * math.log(Nr)
    lo = OUTER[0] * nk_ibz * nq * per * C_EXX[0] * NONEXX / H
    hi = OUTER[1] * nk_ibz * nq * per * C_EXX[1] * NONEXX / H
    mem = Nr * nbnd * nk_full * 16 / 1e9          # exxbuff if it holds all nbnd bands (one copy per pool), GB
    return dict(grid=g, Nr=Nr, nk_ibz=nk_ibz, nq=nq, nbnd=nbnd, core_h=[lo, hi], exxbuff_GB=mem,
                wall_h_64=[lo / 64, hi / 64])


def kgrid(r):
    """HSE/W90 grid along x, the sheet's density (16 × 3.11 Å): [110] 16 k / q 8, [100] 12 k / q 6."""
    return (16, 8) if r["direction"] == "110" else (12, 6)


def nbnd_w90(r):
    """Bands for W90 with the outer window at VBM + 11.05 eV (HSE; cost test T2): 20 per f.u. on the sheet, up to 24
    with the ribbons' extra vacuum; + 2 per H. Fixed per ribbon by the A4a-2 PBE band survey."""
    return [20 * r["nSc"] + 2 * r["nH"], 24 * r["nSc"] + 2 * r["nH"]]


def nbnd_survey(r):
    """A4a-2 bands run: enough bands to see where the W90 outer window falls (28 per f.u. + 1 per H); above 15 Å vacuum
    the empty bands scale with the cell cross-section (make_a4a.py, 2026-11-13)."""
    nb = 28 * r["nSc"] + r["nH"]
    ref = ribbon(r["direction"], r["n"], r["h"])
    f = r["Wy"] * r["Wz"] / (ref["Wy"] * ref["Wz"])
    return nb if f < 1 + 1e-6 else r["nocc"] + int(round((nb - r["nocc"]) * f))


def add(a, b):
    return [a[0] + b[0], a[1] + b[1]]


def a4a_block(r, spin_full=True):
    """A4a-1 relax + A4a-2 spin single points, k-test, bands. Bare: relax FM start (nspin 2); H: relax nspin 1."""
    rows = {}
    if r["h"]:
        rows["relax (nspin 1)"] = pbe(r, "relax", 8, 300)
        rows["SCF FM check"] = pbe(r, "scf", 8, 87, nspin=2)
    else:
        rows["relax (FM start, nspin 2)"] = pbe(r, "relax", 8, 480, nspin=2)
        rows["SCF AFM edges + AFM planes"] = [2 * x for x in pbe(r, "scf", 8, 87, nspin=2)]
        rows["k-test k30: NM + best magnetic"] = add(pbe(r, "scf", 16, 40), pbe(r, "scf", 16, 87, nspin=2))
    rows["SCF NM k15"] = pbe(r, "scf", 8, 40)
    rows["bands Γ–X 31 k + projwfc (band survey)"] = pbe(r, "nscf", 31, 18, nbnd=nbnd_survey(r))
    return rows


def a5_relax(r, magnetic=False):
    return pbe(r, "relax", 8, 480 if magnetic else 300, nspin=2 if magnetic else 1)


def pbe_w90(r):
    nk, _ = kgrid(r)
    nb = nbnd_w90(r)[1]
    x = pbe(r, "nscf", nk, 18, nbnd=nb)
    return [v * 1.15 for v in x]                    # + pw2wannier90 (15 %), as F2


def main():
    R = {k: ribbon(*k) for k in [("100", 3, False), ("110", 4, False), ("100", 3, True), ("110", 8, False),
                                 ("100", 6, False), ("110", 11, False), ("100", 8, False)]}
    out = {"assumptions": dict(A_P=A_P, D=D_PLANE, VAC=VAC, PBE_FACTOR=PBE_FACTOR, C_EXX=C_EXX, OUTER=OUTER)}
    print("Ribbons (scaled A2 geometry, 15 Å vacuum):")
    for r in R.values():
        print(f"  {r['label']:12s} W {r['W_nm']:.2f} nm  nat {r['nat']:3d} (Sc {r['nSc']}, H {r['nH']})  nel {r['nel']}  "
              f"cell {r['L']:.2f} × {r['Wy']:.1f} × {r['Wz']:.1f} Å")
    # A4a
    tot = [0, 0]
    print("\nA4a (PBE, PAW 90/720 Ry, k 15 along x; 8 irreducible):")
    out["A4a"] = {}
    for k in [("100", 3, False), ("110", 4, False), ("100", 3, True)]:
        rows = a4a_block(R[k])
        lab = R[k]["label"]
        est1 = [v for n, v in rows.items() if n.startswith("relax")][0]
        rows = {n: v for n, v in rows.items() if not n.startswith("relax")}
        s = [sum(v[0] for v in rows.values()), sum(v[1] for v in rows.values())]
        out["A4a"][lab] = dict(A4a1_estimate=est1, A4a1_actual=A4A1_ACTUAL[lab], A4a2_rows=rows, A4a2_total=s)
        out["A4a"][lab]["A4a2_actual"] = A4A2_ACTUAL[lab]
        tot = add(tot, [A4A1_ACTUAL[lab] + A4A2_ACTUAL[lab]] * 2)
        print(f"  {lab:12s} A4a-1 actual {A4A1_ACTUAL[lab]:5.0f} (estimate {est1[0]:.0f}–{est1[1]:.0f});  A4a-2 actual "
              f"{A4A2_ACTUAL[lab]:4.0f} (estimate {s[0]:4.0f}–{s[1]:4.0f}):  " + "; ".join(f"{n} {v[0]:.0f}–{v[1]:.0f}" for n, v in rows.items()))
    vac = pbe(R[("110", 4, False)], "scf", 8, 40)
    vac = [v * (VAC + 5 + R[("110", 4, False)]["Wy"] - VAC) / R[("110", 4, False)]["Wy"] * 1.4 for v in vac]
    out["A4a"]["vacuum check [110] n=4, 20 Å"] = vac
    out["A4a_total"] = tot
    print(f"  vacuum check (NM SCF at 20 Å, inside the [110] n=4 A4a-2 job): estimate {vac[0]:.0f}–{vac[1]:.0f}\n  A4a total (actual): "
          f"{tot[0]:.0f}–{tot[1]:.0f} core-h (plan 600–1,200)")
    # HSE-W90 per ribbon + PBE twin W90
    print("\nHSE06-W90 per ribbon (NC 80 Ry, ecutfock 160, IBZ + open_grid, fixed occupations):")
    out["hse_w90"] = {}
    for k, r in R.items():
        nk, nq = kgrid(r)
        lo_b, hi_b = nbnd_w90(r)
        a, b = hse(r, nk, nq, lo_b), hse(r, nk, nq, hi_b)
        rng = [a["core_h"][0], b["core_h"][1]]
        out["hse_w90"][r["label"]] = dict(nk=nk, nq=nq, nbnd=[lo_b, hi_b], grid=a["grid"], core_h=rng,
                                         wall_h_64=[rng[0] / 64, rng[1] / 64], exxbuff_GB=[a["exxbuff_GB"], b["exxbuff_GB"]],
                                         pbe_w90=pbe_w90(r))
        print(f"  {r['label']:12s} k {nk} (IBZ {a['nk_ibz']}) q {nq} nbnd {lo_b}–{hi_b} grid {a['grid']}: "
              f"{rng[0]:6.0f}–{rng[1]:6.0f} core-h (wall @64 {rng[0]/64:5.1f}–{rng[1]/64:5.1f} h; exxbuff "
              f"{a['exxbuff_GB']:.0f}–{b['exxbuff_GB']:.0f} GB)   PBE-W90 {out['hse_w90'][r['label']]['pbe_w90'][0]:.0f}–"
              f"{out['hse_w90'][r['label']]['pbe_w90'][1]:.0f}")
    # HSE gap single point (plan A5): SCF with nbnd = nocc + 10 %, k IBZ only (no W90)
    gap = {}
    for k in [("110", 8, False), ("100", 6, False)]:
        r = R[k]
        nk, nq = kgrid(r)
        g = hse(r, nk, nq, int(r["nocc"] * 1.1) + 8)
        gap[r["label"]] = g["core_h"]
    out["hse_gap_single_points"] = gap
    # scenarios
    H_ = out["hse_w90"]
    plan = dict()
    plan["A4b (plan): HSE-W90 on the 3 calibration ribbons"] = [sum(H_[l]["core_h"][i] + H_[l]["pbe_w90"][i]
                                                                 for l in ("[100] n=3", "[110] n=4", "[100] n=3 H"))
                                                             for i in (0, 1)]
    a5_pbe = [sum(a5_relax(R[k])[i] + pbe(R[k], "scf", 8, 40)[i] + H_[R[k]["label"]]["pbe_w90"][i]
                  for k in [("110", 8, False), ("100", 6, False)]) for i in (0, 1)]
    plan["A5 (plan): hold-outs PBE relax + SCF + PBE-W90 (NM)"] = a5_pbe
    plan["A5 (plan): HSE gap single points"] = [sum(v[i] for v in gap.values()) for i in (0, 1)]
    lr1 = [sum(H_[l]["core_h"][i] for l in ("[110] n=8", "[100] n=6")) - plan["A5 (plan): HSE gap single points"][i]
           for i in (0, 1)]
    plan["LR1 option A: upgrade the hold-outs' HSE gap points to HSE-W90 (+)"] = lr1
    wide = [sum(a5_relax(R[k])[i] + pbe(R[k], "scf", 8, 40)[i] + H_[R[k]["label"]]["pbe_w90"][i] + H_[R[k]["label"]]["core_h"][i]
                for k in [("110", 11, False), ("100", 8, False)]) for i in (0, 1)]
    plan["LR1 option B: + [110] n=11 and [100] n=8 (3.3 nm), PBE relax + PBE-W90 + HSE-W90 (+)"] = wide
    out["scenarios"] = plan
    print("\nStages (core-h):")
    for k, v in plan.items():
        print(f"  {k:92s} {v[0]:6.0f}–{v[1]:6.0f}")
    spent = sum(SPENT.values())
    base = add(add(tot, plan["A4b (plan): HSE-W90 on the 3 calibration ribbons"]),
               add(plan["A5 (plan): hold-outs PBE relax + SCF + PBE-W90 (NM)"], plan["A5 (plan): HSE gap single points"]))
    optA = add(base, lr1)
    optB = add(optA, wide)
    pre = spent - SPENT["A4a1"] - SPENT["A4a2"]      # A0–A2 + cost tests, outside the A4/A5 blocks
    out["totals"] = dict(spent=spent, plan_A4_A5=base, with_LR1_A=optA, with_LR1_A_B=optB,
                         project_with_LR1_A=[pre + v for v in optA], project_with_LR1_A_B=[pre + v for v in optB])
    print(f"\n  spent so far ({', '.join(f'{k} {v:.0f}' for k, v in SPENT.items())}): {spent:.0f}")
    print(f"  remaining A4b + A5 as planned: {base[0] - tot[0]:.0f}–{base[1] - tot[1]:.0f}")
    print(f"  A4a + A4b + A5 as planned (A4a at its actual cost): {base[0]:.0f}–{base[1]:.0f}")
    print(f"  ... + LR1 option A: {optA[0]:.0f}–{optA[1]:.0f}   (+ A3 3-plane slab, not costed here)")
    print(f"  ... + LR1 options A + B: {optB[0]:.0f}–{optB[1]:.0f}")
    print(f"  PROJECT TOTAL incl. A0–A2, the cost tests and the post-Gate-M runs ({pre:.0f}): option A {pre + optA[0]:.0f}–{pre + optA[1]:.0f};"
          f" A + B {pre + optB[0]:.0f}–{pre + optB[1]:.0f}   (+ A3, not costed)")
    # Gate M criterion 5 failed (r110_4 gap 15 → 20 Å moved 12.3 meV; README): option A with 20 Å vacuum
    five = [("100", 3, False), ("110", 4, False), ("100", 3, True), ("110", 8, False), ("100", 6, False)]
    hold = {("110", 8, False), ("100", 6, False)}

    def option_a(vac):
        res = {}
        for k in five:
            r = ribbon(*k, vacuum=vac)
            nk, nq = kgrid(r)
            lo_b, hi_b = nbnd_w90(r)
            c = add([hse(r, nk, nq, lo_b)["core_h"][0], hse(r, nk, nq, hi_b)["core_h"][1]], pbe_w90(r))
            if k in hold:
                c = add(c, add(a5_relax(r), pbe(r, "scf", 8, 40)))
            res[k] = c
        return res
    a15, a20 = option_a(VAC), option_a(VAC + 5)
    s15 = [sum(v[i] for v in a15.values()) for i in (0, 1)]
    assert all(abs(s15[i] - (optA[i] - tot[i])) < 1e-6 * optA[i] for i in (0, 1)), "option-A per-ribbon sum mismatch"
    s20 = [sum(v[i] for v in a20.values()) for i in (0, 1)]
    # vacuum per ribbon after the 20 Å checks (job 142557, 12 Nov; README Results): [100] bare passes (0.0 meV) → 15 Å,
    # the [100] n=6 hold-out shares its edges → 15 Å; [100] H fails (17.7 meV) → 20 Å; [110] → 20 Å (12.3 meV)
    VAC_SET = {("100", 3, False): VAC, ("110", 4, False): VAC + 5, ("100", 3, True): VAC + 5,
               ("110", 8, False): VAC + 5, ("100", 6, False): VAC}
    by = {VAC: a15, VAC + 5: a20, VAC + 10: option_a(VAC + 10)}
    # A4b-2 as staged (13 Nov): HSE nbnd per calibration ribbon from the A4b-1 NC NSCF survey + 0.5 eV margin
    # (make_a4b.NB_HSE), in place of the 20–24 per f.u. model; the hold-outs keep the model.
    NB_A4B2 = {("100", 3, False): 267, ("110", 4, False): 214, ("100", 3, True): 273}

    def a4b2_row(k):
        r_ = ribbon(*k, vacuum=VAC_SET[k])
        nk_, nq_ = kgrid(r_)
        return add(hse(r_, nk_, nq_, NB_A4B2[k])["core_h"], pbe_w90(r_))
    dec = [sum((a4b2_row(k) if k in NB_A4B2 else by[VAC_SET[k]][k])[i] for k in five) for i in (0, 1)]
    a4b2_hse = [sum(hse(ribbon(*k, vacuum=VAC_SET[k]), *kgrid(ribbon(*k, vacuum=VAC_SET[k])), NB_A4B2[k])["core_h"][i]
                    for k in NB_A4B2) for i in (0, 1)]
    # LR3 (decision 12 Nov): [110] n=4 H at 20 Å, A4a (relax nspin 1 + FM check + NM + bands) and HSE-W90 + PBE-W90
    rh = ribbon("110", 4, True, vacuum=VAC + 5)
    nk, nq = kgrid(rh)
    lo_b, hi_b = nbnd_w90(rh)
    # the law misses relaxations by a ribbon-dependent factor: actual / law-low = ×1.47 for the [110] n=4 bare relax,
    # ×5.1 for the [100] n=3 H relax (more BFGS steps); r110_4h is both, so its relax is bracketed by the two.
    # The rest scales like the [100] H A4a-2 rows (actual / law-low 0.62).
    rows_h = a4a_block(rh)
    law = {k: a4a_block(R[k]) for k in (("110", 4, False), ("100", 3, True))}
    f_rel = [A4A1_ACTUAL["[110] n=4"] / law[("110", 4, False)]["relax (FM start, nspin 2)"][0],
             A4A1_ACTUAL["[100] n=3 H"] / law[("100", 3, True)]["relax (nspin 1)"][0]]
    f_2 = A4A2_ACTUAL["[100] n=3 H"] / sum(v[0] for n_, v in law[("100", 3, True)].items() if not n_.startswith("relax"))
    rest = f_2 * sum(v[0] for n_, v in rows_h.items() if not n_.startswith("relax"))
    h110_a4a = [0.0, 0.0]                        # remaining: none; relax and Gate M runs are spent (estimate rest above)
    model_nb = (lo_b, hi_b)
    lo_b, hi_b = 205, 205                        # A4b-2 as staged (make_a4b.NB_HSE); survey 201–217, model 168–200
    h110_w90 = add([hse(rh, nk, nq, lo_b)["core_h"][0], hse(rh, nk, nq, hi_b)["core_h"][1]], pbe_w90(rh))
    # optional 25 Å check (NM SCF k 15) on the two ribbons that failed at 20 Å (scaled like the 20 Å checks:
    # estimate 25–50 for the [100] pair, actual 40)
    chk25 = [sum(pbe(ribbon(*k, vacuum=VAC + 10), "scf", 8, 40)[i] for k in (("110", 4, False), ("100", 3, True)))
             for i in (0, 1)]
    vac25 = {k: (VAC + 10 if v > VAC else VAC) for k, v in VAC_SET.items()}
    dec25 = [sum(by[vac25[k]][k][i] for k in five) for i in (0, 1)]
    rh25 = ribbon("110", 4, True, vacuum=VAC + 10)
    nk25, nq25 = kgrid(rh25)
    lo25, hi25 = nbnd_w90(rh25)
    h110_w90_25 = add([hse(rh25, nk25, nq25, lo25)["core_h"][0], hse(rh25, nk25, nq25, hi25)["core_h"][1]], pbe_w90(rh25))
    proj = lambda rem: [pre + tot[i] + rem[i] for i in (0, 1)]          # noqa: E731
    rem_dec = add(dec, add(h110_a4a, h110_w90))
    rem_25 = add(add(dec25, add(h110_a4a, h110_w90_25)), chk25)
    # spent since (13 Nov): the 25 Å check (chk25) ran; A4b-1 did the PBE-W90 of the four calibration ribbons at the
    # decided vacuum. In the 25 Å scenario only [100] bare (15 Å) keeps its A4b-1 run; the others would redo it.
    a4b_ribbons = {("100", 3, False): ("100", 3, False), ("110", 4, False): ("110", 4, False),
                   ("100", 3, True): ("100", 3, True)}
    pw_done = [sum(pbe_w90(ribbon(*k, vacuum=VAC_SET[k]))[i] for k in a4b_ribbons) + pbe_w90(rh)[i] for i in (0, 1)]
    pw_keep25 = pbe_w90(ribbon("100", 3, False, vacuum=VAC))
    rem_dec = [rem_dec[i] - pw_done[i] for i in (0, 1)]
    rem_25 = [rem_25[i] - chk25[i] - pw_keep25[i] for i in (0, 1)]
    out["vacuum"] = dict(per_ribbon_15={str(k): v for k, v in a15.items()}, per_ribbon_20={str(k): v for k, v in a20.items()},
                         per_ribbon_25={str(k): v for k, v in by[VAC + 10].items()},
                         vacuum_set={str(k): v for k, v in VAC_SET.items()}, vac20_checks_actual=VAC20_ACTUAL,
                         remaining_optA_15=s15, remaining_optA_20_all=s20, remaining_optA_decided=dec,
                         h110_A4a=h110_a4a, h110_hse_pbe_w90_20=h110_w90, h110_hse_pbe_w90_25=h110_w90_25,
                         check_25A=chk25, remaining_decided_plus_h110=rem_dec, remaining_25A_scenario=rem_25,
                         project_optA_20_all=proj(s20), project_decided_plus_h110=proj(rem_dec),
                         project_25A_scenario=proj(rem_25), a4b1_replaced_pbe_w90=pw_done,
                         vac25_check_actual=VAC25_ACTUAL)
    print("\nVacuum (Gate M criterion 5): remaining A4b + A5 + LR1 option A (5 ribbons)")
    print(f"  15 Å everywhere: {s15[0]:.0f}–{s15[1]:.0f};  20 Å everywhere: {s20[0]:.0f}–{s20[1]:.0f};  as decided "
          f"([100] bare 15 Å, [100] H and [110] 20 Å): {dec[0]:.0f}–{dec[1]:.0f}")
    print(f"  + [110] n=4 H (LR3) at 20 Å: relax actual {H110_RELAX_ACTUAL:.0f} (bracket {rows_h['relax (nspin 1)'][0] * f_rel[0]:.0f}–"
          f"{rows_h['relax (nspin 1)'][0] * f_rel[1]:.0f}), Gate M runs actual {H110_GATEM_ACTUAL:.0f} (estimate {rest:.0f}; nbnd "
          f"{nbnd_survey(rh)}); both spent. HSE-W90 + PBE-W90 at the measured nbnd {lo_b}–{hi_b} (model {model_nb[0]}–{model_nb[1]}) "
          f"{h110_w90[0]:.0f}–{h110_w90[1]:.0f}  →  remaining {rem_dec[0]:.0f}–{rem_dec[1]:.0f}")
    print(f"  PROJECT TOTAL (spent {spent:.0f} incl. the 20 Å checks {VAC20_ACTUAL:.0f} and the [110] H relax): 15 Å option A "
          f"{proj(s15)[0]:.0f}–{proj(s15)[1]:.0f}; as decided + [110] n=4 H {proj(rem_dec)[0]:.0f}–{proj(rem_dec)[1]:.0f}")
    print(f"  if 20 Å is not enough either: 25 Å check on [110] n=4 + [100] n=3 H {chk25[0]:.0f}–{chk25[1]:.0f}; then 25 Å "
          f"for those (and [110] n=8, [110] n=4 H; HSE-W90 [110] n=4 H {h110_w90_25[0]:.0f}–{h110_w90_25[1]:.0f}): "
          f"remaining {rem_25[0]:.0f}–{rem_25[1]:.0f}, project {proj(rem_25)[0]:.0f}–{proj(rem_25)[1]:.0f}")
    print(f"  (A4b-1 spent {SPENT['A4b1']:.0f} and the 25 Å check {VAC25_ACTUAL:.0f}; removed from 'remaining': the PBE-W90 "
          f"law rows A4b-1 replaced, {pw_done[0]:.0f}–{pw_done[1]:.0f}, and the 25 Å check. The 25 Å scenario's redo of "
          f"A4b-1 is counted as PBE-W90 only.)")
    h4 = hse(rh, nk, nq, 205)["core_h"]
    a4b2 = [a4b2_hse[i] + h4[i] for i in (0, 1)]
    out["vacuum"]["A4b2_hse_scf_as_staged"] = a4b2
    print(f"  A4b-2 as staged (nbnd 267 / 273 / 214 / 205): HSE06 SCFs {a4b2[0]:.0f}–{a4b2[1]:.0f} core-h "
          f"(+ the pp.x / projwfc / open_grid / pw2wannier90 steps)")
    print("  (hold-outs: band counts from the per-f.u. model; more vacuum adds free-electron bands, not counted)")

    # ---- after A4b-2 (14 Nov): what is left is A5 (the two hold-outs) and A3, at the recalibrated EXX constant
    def hse_rib(r_, nk_, nq_, nb_):
        g = [good_fft(x / BOHR * math.sqrt(160.0) / math.pi) for x in (r_["L"], r_["Wy"], r_["Wz"])]
        Nr = g[0] * g[1] * g[2]
        per = nb_ ** 2 * Nr * math.log(Nr) * (nk_ // 2 + 1) * nq_ * NONEXX / H
        return [OUTER[0] * per * C_EXX_RIBBON[0], OUTER[1] * per * C_EXX_RIBBON[1]]
    # HSE nbnd per Sc from A4b-2: [110] at 20 Å 25.6–26.8 (205, 214 for 8 Sc), [100] bare at 15 Å 22.25 (267 for 12 Sc);
    # wider ribbons carry relatively fewer vacuum bands, so the range runs down to 24 / 20 per Sc
    A5 = {("110", 8, False, 20.0): (24.0, 26.75), ("100", 6, False, 15.0): (20.0, 22.25)}
    a5 = {}
    for (d_, n_, h_, v_), (lo_s, hi_s) in A5.items():
        r_ = ribbon(d_, n_, h_, vacuum=v_)
        nk_, nq_ = kgrid(r_)
        nbs = (int(round(lo_s * r_["nSc"])), int(round(hi_s * r_["nSc"])))
        hx = hse_rib(r_, nk_, nq_, nbs[0])[0], hse_rib(r_, nk_, nq_, nbs[1])[1]
        pre_ = add(a5_relax(r_), add(pbe(r_, "scf", 8, 40), pbe_w90(r_)))
        a5[f"[{d_}] n={n_} ({v_:.0f} Å)"] = dict(nbnd=nbs, hse=list(hx), pbe_side=pre_, total=add(list(hx), pre_),
                                                  wall_h_64=[hx[0] / 64, hx[1] / 64])
    rem_a5 = [sum(v_["total"][i] for v_ in a5.values()) for i in (0, 1)]
    A3_EST = (100.0, 200.0)
    spent_b2 = spent_upto("A4b2")                  # as written on 14 Nov, before A5/A3 ran
    proj_now = [spent_b2 + rem_a5[i] + A3_EST[i] for i in (0, 1)]
    out["after_A4b2"] = dict(C_EXX_ribbon=C_EXX_RIBBON, A4b2_actual=dict(A4B2_ACTUAL, w90=A4B2_W90_ACTUAL), A5=a5,
                             remaining_A5=rem_a5, A3_estimate=A3_EST, spent=spent_b2, project_total=proj_now)
    print(f"\nAfter A4b-2 (spent {spent_b2:.0f}; A4b-2 itself {SPENT['A4b2']:.0f}), recalibrated EXX constant "
          f"{C_EXX_RIBBON[0]:.1e}–{C_EXX_RIBBON[1]:.1e}:")
    for k_, v_ in a5.items():
        print(f"  A5 {k_}: nbnd {v_['nbnd'][0]}–{v_['nbnd'][1]}, HSE {v_['hse'][0]:.0f}–{v_['hse'][1]:.0f} core-h "
              f"(one node {v_['wall_h_64'][0]:.0f}–{v_['wall_h_64'][1]:.0f} h), + relax/SCF/PBE-W90 "
              f"{v_['pbe_side'][0]:.0f}–{v_['pbe_side'][1]:.0f}")
    print(f"  remaining A5 {rem_a5[0]:.0f}–{rem_a5[1]:.0f} + A3 ≈ {A3_EST[0]:.0f}–{A3_EST[1]:.0f}  →  PROJECT TOTAL "
          f"{proj_now[0]:.0f}–{proj_now[1]:.0f} core-h")

    # ---- after the A5-1 relaxations (15 Nov): A5 is option 2 (decision 14 Nov; dft/A5/README.md)
    st2 = {k_: [64 * A5_NM_SCF_CPU_S[k_] * 1.07 / 3600 * f for f in A5_STAGE2_RATIO] for k_ in A5_NM_SCF_CPU_S}
    st2_tot = [sum(v_[i] for v_ in st2.values()) for i in (0, 1)]
    spent_r = spent_upto("A5_relax")               # this block as written on 15 Nov, before A5-1b ran
    out["after_A5_relax"] = dict(relax_actual=A5_RELAX_ACTUAL, stage2_estimate=st2, stage2_total=st2_tot, spent=spent_r,
                                 stage2_actual=A5_1B_ACTUAL, spent_after_stage2=spent_upto("A5_1b"))
    print(f"\nAfter the A5-1 relaxations (spent {spent_r:.0f}; relaxations " +
          ", ".join(f"{k_} {v_:.0f}" for k_, v_ in A5_RELAX_ACTUAL.items()) + "):")
    for k_, v_ in st2.items():
        print(f"  A5-1b {k_} (NM/FM/AFM-edges, bands survey, projwfc): {v_[0]:.0f}–{v_[1]:.0f} core-h "
              f"(one node {v_[0] / 64:.1f}–{v_[1] / 64:.1f} h)")
    print(f"  A5-1b total {st2_tot[0]:.0f}–{st2_tot[1]:.0f} core-h")
    print(f"  A5-1b actual (15 Nov): " + ", ".join(f"{k_} {v_:.0f}" for k_, v_ in A5_1B_ACTUAL.items())
          + f" = {SPENT['A5_1b']:.0f} core-h; spent {spent_upto('A5_1b'):.0f}")
    print(f"  A3-1 slab relaxation actual {SPENT['A3_relax']:.1f} core-h; spent {spent_upto('A3_relax'):.0f}")
    print(f"  A5-2 job 1 actual (15 Nov; estimate 230–390): " + ", ".join(f"{k_} {v_:.0f}" for k_, v_ in A5_2_JOB1_ACTUAL.items())
          + f" = {SPENT['A5_2_job1']:.0f} core-h; A3-2 PBE twin + vacuum check {SPENT['A3_2_pbe']:.1f}; spent "
          f"{spent_upto('A3_2_pbe'):.0f}")
    out["after_A5_relax"]["A5_2_job1_actual"] = A5_2_JOB1_ACTUAL
    out["after_A5_relax"]["A3_2_pbe_actual"] = A3_2_PBE_ACTUAL
    print(f"  A5-2 job 2 + slab PBE Wannier90 actual {SPENT['A5_2_job2']:.0f} core-h (estimate 70–225); A3-2 HSE06 "
          f"{SPENT['A3_2_hse']:.0f} (A3 estimate 145–270 for A3-2); A5-3 r110_8 {SPENT['A5_3_r110_8']:.0f} (registered "
          f"450–960, + 10–20 % on two nodes); spent {spent_upto('A5_3_r110_8'):.0f}")
    print(f"  A3-3 slab HSE Wannier90 actual {SPENT['A3_3_hw90']:.1f} core-h; spent {spent_upto('A3_3_hw90'):.0f}")
    print(f"  A5-3 r100_6 actual {SPENT['A5_3_r100_6']:.0f} core-h (registered 480–1,030); spent "
          f"{spent_upto('A5_3_r100_6'):.0f}")
    out["after_A5_relax"].update(A5_2_job2_actual=A5_2_JOB2_ACTUAL, A3_2_hse_actual=A3_2_HSE_ACTUAL,
                                 A5_3_actual=A5_3_ACTUAL, A3_3_hw90_actual=A3_3_HW90_ACTUAL,
                                 spent=spent_upto("A5_3_r100_6"))
    (HERE / "cost_a4.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
