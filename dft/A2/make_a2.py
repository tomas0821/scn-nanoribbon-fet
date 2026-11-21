#!/usr/bin/env python3
"""Generate the A2 inputs: the M-ScN sheet (P4/nmm two-plane rocksalt-(001) bilayer, 2 f.u.), NC set at 80 Ry,
HSE06 + Wannier90 and the PBE twin (plan §4 A2). Settings follow A1, including its lessons.

Geometry: the A0 PBE (PAW 90) 2Dxy vc-relax (a = 3.1179 Å, Sc/N planes 2.174/2.177 Å apart), scaled by the bulk
ratio a_exp / a_PBE = 4.501 / 4.5134, the same convention as bulk A1 (experimental a). Vacuum: c = 22.5 Å (A0).
A1 lessons applied:
  * 16×16 Wannier grid (bulk 8³ failed the interpolation criterion, 16³ fixed it); route O (IBZ + open_grid.x);
  * HSE q 8×8×1 (bulk q 4 → 8 moved the gaps 12–15 meV); a q 4 vs q 8 check at k 8×8 here;
  * pw2wannier90 on ≤ 16 ranks; srun < /dev/null in the tag loop.

Stage A2-1 (run_a2.sh, one 64-rank job, restartable tags):
  pbe_scf        PBE SCF 16×16×1, m-v 0.005 Ry, nbnd 56 (density, E_F)
  pbe_pot        pp.x: bare + Hartree potential (plot_num 11) as a cube → vacuum level (frozen window must sit below it)
  pbe_bands      PBE bands Γ–X̄–M̄–Γ (reference for the PBE-twin interpolation criterion)
  pbe_nscf_k16   PBE NSCF on the explicit full 16×16×1 list, nosym/noinv (pw2wannier90, PBE twin)
  pbe_c30        PBE SCF at c = 30 Å, 12×12×1, nbnd 24 (vacuum convergence of the gaps)
  hse_k8_q4      HSE06 SCF, 8×8×1 IBZ, q 4×4×1, nbnd 24 (q check)
  hse_k8_q8      HSE06 SCF, 8×8×1 IBZ, q 8×8×1, nbnd 24 (q check)
  hse_ibz_k16    HSE06 SCF, 16×16×1 IBZ, q 8×8×1, nbnd 56 (production)
  hse_ibz_k16_open  open_grid.x → the full 16×16×1 grid for pw2wannier90
Stage A2-2 (Wannier90, after reading the A2-1 eigenvalues and the vacuum level): `win()` below.
"""
from pathlib import Path

import numpy as np

BOHR = 0.529177210903
PSEUDO_DIR = "/path/to/scn_fet/pseudo"
A_PBE, C = 3.117899, 22.5                     # A0 sheet_paw90_vcrelax (Å)
Z_PBE = {"Sc1": 0.4516957094, "N1": 0.4516288717, "N2": 0.5483711283, "Sc2": 0.5483042906}
SCALE = 4.501 / 4.5134                        # bulk a_exp / a_PBE(PAW 90, A0)
A = A_PBE * SCALE
ECUT, DUAL, FOCK = 80, 4, 160
NBND, NBND_SMALL = 56, 24
KW = 16
HERE = Path(__file__).resolve().parent
OUT = HERE / "inputs"
PATH = [("G", (0.0, 0.0, 0.0)), ("X", (0.5, 0.0, 0.0)), ("M", (0.5, 0.5, 0.0)), ("G", (0.0, 0.0, 0.0))]


def zfrac(c=C):
    """Fractional z of the four atoms in a cell of height c: PBE offsets from the mid-plane, scaled by SCALE."""
    return {k: 0.5 + (v - 0.5) * C * SCALE / c for k, v in Z_PBE.items()}


def species(c=C):
    z = zfrac(c)
    return ("ATOMIC_SPECIES\n  Sc 44.955908 Sc_ONCV_PBE-1.0.oncvpsp.upf\n  N  14.007 N.oncvpsp.upf\n"
            "ATOMIC_POSITIONS crystal\n"
            f"  Sc 0.0 0.0 {z['Sc1']:.10f}\n  N  0.5 0.5 {z['N1']:.10f}\n"
            f"  N  0.0 0.0 {z['N2']:.10f}\n  Sc 0.5 0.5 {z['Sc2']:.10f}\n")


def full_grid(n=KW):
    return [(i / n, j / n, 0.0) for i in range(n) for j in range(n)]


def kblock_list(n=KW):
    ks = full_grid(n)
    return f"K_POINTS crystal\n{len(ks)}\n" + "".join(f"  {a:.8f} {b:.8f} {c:.8f} {1 / len(ks):.10f}\n" for a, b, c in ks)


def control(calc, prefix, extra=""):
    return (f"&control\n  calculation='{calc}', prefix='{prefix}', outdir='./tmp_{prefix}', pseudo_dir='{PSEUDO_DIR}',\n"
            f"  verbosity='high', max_seconds=80000{extra}\n/\n")


def system(occ, nbnd, extra="", c=C):
    return (f"&system\n  ibrav=6, celldm(1)={A / BOHR:.10f}, celldm(3)={c / A:.10f}, nat=4, ntyp=2,\n"
            f"  ecutwfc={ECUT}, ecutrho={ECUT * DUAL},\n  {occ}, nbnd={nbnd}{extra}\n/\n")


SMEAR = "occupations='smearing', smearing='mv', degauss=0.005"
FIXED = "occupations='fixed'"


def hse(nq):
    return f",\n  input_dft='hse', nqx1={nq}, nqx2={nq}, nqx3=1, ecutfock={FOCK}"


def electrons(thr="1.0d-10", extra=""):
    return f"&electrons\n  conv_thr={thr}, mixing_beta=0.4{extra}\n/\n"


def stage1():
    OUT.mkdir(exist_ok=True)
    tags = []

    def w(tag, text):
        (OUT / f"{tag}.in").write_text(text)
        tags.append(tag)

    w("pbe_scf", control("scf", "pbe_scf", ", tstress=.true., tprnfor=.true.") + system(SMEAR, NBND) + electrons()
      + species() + f"K_POINTS automatic\n  {KW} {KW} 1 0 0 0\n")
    (OUT / "pbe_pot.in").write_text(
        "&inputpp\n  prefix='pbe_scf', outdir='./tmp_pbe_scf', filplot='pbe_scf.pot', plot_num=11\n/\n"
        "&plot\n  nfile=1, filepp(1)='pbe_scf.pot', weight(1)=1.0, iflag=3, output_format=6, fileout='pbe_scf_pot.cube'\n/\n")
    tags.append("pbe_pot")
    seg = "".join(f"  {x:.4f} {y:.4f} {z:.4f} {60 if i < len(PATH) - 1 else 1}\n" for i, (_, (x, y, z)) in enumerate(PATH))
    w("pbe_bands", control("bands", "pbe_scf") + system(SMEAR, NBND) + electrons(extra=", diago_full_acc=.true.")
      + species() + f"K_POINTS crystal_b\n  {len(PATH)}\n" + seg)
    w("pbe_nscf_k16", control("nscf", "pbe_scf") + system(SMEAR, NBND, ",\n  nosym=.true., noinv=.true.")
      + electrons(extra=", diago_full_acc=.true.") + species() + kblock_list())
    w("pbe_c30", control("scf", "pbe_c30") + system(SMEAR, NBND_SMALL, c=30.0) + electrons()
      + species(30.0) + "K_POINTS automatic\n  12 12 1 0 0 0\n")
    for nq in (4, 8):
        w(f"hse_k8_q{nq}", control("scf", f"hse_k8_q{nq}") + system(FIXED, NBND_SMALL, hse(nq)) + electrons("1.0d-9")
          + species() + "K_POINTS automatic\n  8 8 1 0 0 0\n")
    w("hse_ibz_k16", control("scf", "hse_ibz_k16") + system(FIXED, NBND, hse(8))
      + electrons("1.0d-9", ", diago_full_acc=.true.") + species() + f"K_POINTS automatic\n  {KW} {KW} 1 0 0 0\n")
    (OUT / "hse_ibz_k16_open.in").write_text("&inputpp\n  outdir='./tmp_hse_ibz_k16', prefix='hse_ibz_k16'\n/\n")
    tags.append("hse_ibz_k16_open")
    (HERE / "tags.txt").write_text("\n".join(tags) + "\n")
    z = zfrac()
    print(f"wrote {len(tags)} tags to {OUT} | a = {A:.5f} Å, Sc–Sc plane gap {(z['Sc2'] - z['Sc1']) * C:.4f} Å, "
          f"N–N {(z['N2'] - z['N1']) * C:.4f} Å, in-plane Sc–N {A / np.sqrt(2):.4f} Å")


def win(seedname, kpts, dis_froz_max, dis_win_max, nbnd=NBND):
    """A2-2 .win: 16 projection-only WFs (2 × Sc d + 2 × N p); bands 1–10 excluded (2 × (Sc 3s + 3p), 2 × N 2s)."""
    z = zfrac()
    lines = ["num_wann  = 16", f"num_bands = {nbnd - 10}", "exclude_bands = 1-10",
             f"dis_froz_max = {dis_froz_max:.4f}", f"dis_win_max = {dis_win_max:.4f}",
             "dis_num_iter = 4000", "dis_conv_tol = 1.0e-10", "num_iter = 0     ! projection-only, as A1",
             "use_ws_distance = .false.", "write_hr = .true.", "write_xyz = .true.", "write_u_matrices = .true.",
             "bands_plot = .true.", "bands_num_points = 120", "", "begin projections", "Sc:d", "N:p",
             "end projections", "", "begin unit_cell_cart", "ang", f"  {A:.8f} 0 0", f"  0 {A:.8f} 0", f"  0 0 {C:.8f}",
             "end unit_cell_cart", "", "begin atoms_frac",
             f"Sc 0.0 0.0 {z['Sc1']:.10f}", f"N  0.5 0.5 {z['N1']:.10f}", f"N  0.0 0.0 {z['N2']:.10f}",
             f"Sc 0.5 0.5 {z['Sc2']:.10f}", "end atoms_frac", "", "begin kpoint_path"]
    for (l1, c1), (l2, c2) in zip(PATH[:-1], PATH[1:]):
        lines.append(f"{l1} {c1[0]:.6f} {c1[1]:.6f} {c1[2]:.6f}  {l2} {c2[0]:.6f} {c2[1]:.6f} {c2[2]:.6f}")
    lines += ["end kpoint_path", "", f"mp_grid = {KW} {KW} 1", "", "begin kpoints"]
    lines += [f"  {x:.8f} {y:.8f} {z_:.8f}" for x, y, z_ in kpts] + ["end kpoints", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    stage1()
