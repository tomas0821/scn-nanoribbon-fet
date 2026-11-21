#!/usr/bin/env python3
"""Generate the A1 inputs: bulk rocksalt ScN, NC set at 80 Ry, HSE06 + Wannier90 and the PBE twin (plan §4 A1).

Cutoffs from A0 (2026-11-11): NC 80/320 Ry (Sc SG15 + N Dojo, 11 + 5 valence e⁻). Lattice: experimental
a = 4.501 Å (decision 2026-11-11; the literature HSE06 references use 4.48–4.50 Å).

Stage A1-1 (run_a1.sh, one 64-rank job, restartable tags):
  pbe_scf        PBE SCF, 12³, m-v 0.01 Ry, nbnd 30 (reference E_F and charge density)
  pbe_scf_bands  PBE bands along L–Γ–X–W–K–Γ (Wannier-quality reference for the PBE twin)
  pbe_nscf       PBE NSCF on the explicit full 8×8×8 list, nosym/noinv (pw2wannier90 input for the PBE twin)
  hse_ibz        HSE06 SCF, 8×8×8 IBZ, q 4×4×4, ecutfock 160, fixed occupations, nbnd 30
  hse_ibz_open   open_grid.x: unfold hse_ibz to the full grid (route O, to be tested against route N)
  hse_nosym      HSE06 SCF on the explicit full 8×8×8 list, nosym/noinv, q 4×4×4 (route N)
  hse_q8         HSE06 SCF, IBZ, q 8×8×8, nbnd 12 (q convergence of the gaps)
  hse_fock320    HSE06 SCF, IBZ, q 4×4×4, ecutfock 320, nbnd 12 (ecutfock convergence)
Stage A1-2 (Wannier90; windows set from the A1-1 eigenvalues): `win()` below writes the .win/.pw2wan files.
"""
from pathlib import Path

BOHR = 0.529177210903
PSEUDO_DIR = "/path/to/scn_fet/pseudo"
A_EXP = 4.501                        # Å, experiment (B_targets table A)
ECUT, DUAL, FOCK = 80, 4, 160
NBND = 30                            # 5 semicore/N 2s + 3 N p + 5 Sc d + Sc s/p and margin, for disentanglement
KW = 8                               # Wannier grid
HERE = Path(__file__).resolve().parent
OUT = HERE / "inputs"

SPECIES = ("ATOMIC_SPECIES\n  Sc 44.955908 Sc_ONCV_PBE-1.0.oncvpsp.upf\n  N  14.007 N.oncvpsp.upf\n"
           "ATOMIC_POSITIONS crystal\n  Sc 0.0 0.0 0.0\n  N  0.5 0.5 0.5\n")
PATH = [("L", (0.5, 0.5, 0.5)), ("G", (0.0, 0.0, 0.0)), ("X", (0.0, 0.0, 1.0)), ("W", (0.5, 0.0, 1.0)),
        ("K", (0.75, 0.0, 0.75)), ("G", (0.0, 0.0, 0.0))]      # Cartesian, 2π/a (QE tpiba_b for ibrav = 2)


def full_grid(n=KW):
    """kmesh.pl order: crystal coordinates of the reciprocal basis, x slowest."""
    return [(i / n, j / n, k / n) for i in range(n) for j in range(n) for k in range(n)]


def kblock_list(n=KW):
    ks = full_grid(n)
    return f"K_POINTS crystal\n{len(ks)}\n" + "".join(f"  {a:.8f} {b:.8f} {c:.8f} {1 / len(ks):.10f}\n" for a, b, c in ks)


def control(calc, prefix, extra=""):
    return (f"&control\n  calculation='{calc}', prefix='{prefix}', outdir='./tmp_{prefix}', pseudo_dir='{PSEUDO_DIR}',\n"
            f"  verbosity='high', max_seconds=80000{extra}\n/\n")


def system(occ, nbnd, extra=""):
    return (f"&system\n  ibrav=2, celldm(1)={A_EXP / BOHR:.10f}, nat=2, ntyp=2,\n  ecutwfc={ECUT}, ecutrho={ECUT * DUAL},\n"
            f"  {occ}, nbnd={nbnd}{extra}\n/\n")


SMEAR = "occupations='smearing', smearing='mv', degauss=0.01"
FIXED = "occupations='fixed'"


def hse(nq, fock=FOCK, nosym=False):
    s = f",\n  input_dft='hse', nqx1={nq}, nqx2={nq}, nqx3={nq}, ecutfock={fock}"
    return s + (",\n  nosym=.true., noinv=.true." if nosym else "")


def electrons(thr="1.0d-10", extra=""):
    return f"&electrons\n  conv_thr={thr}, mixing_beta=0.5{extra}\n/\n"


def stage1():
    OUT.mkdir(exist_ok=True)
    tags = []

    def w(tag, text):
        (OUT / f"{tag}.in").write_text(text)
        tags.append(tag)

    w("pbe_scf", control("scf", "pbe_scf", ", tstress=.true.") + system(SMEAR, NBND) + electrons()
      + SPECIES + "K_POINTS automatic\n  12 12 12 0 0 0\n")
    path = "".join(f"  {x:.4f} {y:.4f} {z:.4f} {40 if i < len(PATH) - 1 else 1}\n" for i, (_, (x, y, z)) in enumerate(PATH))
    w("pbe_scf_bands", control("bands", "pbe_scf") + system(SMEAR, NBND) + electrons(extra=", diago_full_acc=.true.")
      + SPECIES + f"K_POINTS tpiba_b\n  {len(PATH)}\n" + path)
    w("pbe_nscf", control("nscf", "pbe_scf") + system(SMEAR, NBND, ",\n  nosym=.true., noinv=.true.")
      + electrons(extra=", diago_full_acc=.true.") + SPECIES + kblock_list())
    w("hse_ibz", control("scf", "hse_ibz") + system(FIXED, NBND, hse(4)) + electrons("1.0d-9", ", diago_full_acc=.true.")
      + SPECIES + f"K_POINTS automatic\n  {KW} {KW} {KW} 0 0 0\n")
    (OUT / "hse_ibz_open.in").write_text("&inputpp\n  outdir='./tmp_hse_ibz', prefix='hse_ibz'\n/\n")
    tags.append("hse_ibz_open")
    w("hse_nosym", control("scf", "hse_nosym") + system(FIXED, NBND, hse(4, nosym=True))
      + electrons("1.0d-9", ", diago_full_acc=.true.") + SPECIES + kblock_list())
    w("hse_q8", control("scf", "hse_q8") + system(FIXED, 12, hse(8)) + electrons("1.0d-9")
      + SPECIES + f"K_POINTS automatic\n  {KW} {KW} {KW} 0 0 0\n")
    w("hse_fock320", control("scf", "hse_fock320") + system(FIXED, 12, hse(4, fock=320)) + electrons("1.0d-9")
      + SPECIES + f"K_POINTS automatic\n  {KW} {KW} {KW} 0 0 0\n")
    (HERE / "tags.txt").write_text("\n".join(tags) + "\n")
    print(f"wrote {len(tags)} inputs to {OUT}")


def win(seedname, kpts, dis_froz_max, dis_win_max, dis_win_min=None, nbnd=NBND):
    """Stage A1-2 .win: 8 projection-only WFs (Sc d + N p), semicore and N 2s excluded (bands 1–5).
    kpts: the k list in the order of the run that pw2wannier90 reads (full_grid() for N, open_grid's list for O)."""
    a = A_EXP
    cell = [(-a / 2, 0, a / 2), (0, a / 2, a / 2), (-a / 2, a / 2, 0)]          # QE ibrav = 2
    lines = [f"num_wann  = 8", f"num_bands = {nbnd - 5}", "exclude_bands = 1-5",
             f"dis_froz_max = {dis_froz_max:.4f}", f"dis_win_max = {dis_win_max:.4f}"]
    if dis_win_min is not None:
        lines.append(f"dis_win_min = {dis_win_min:.4f}")
    lines += ["dis_num_iter = 4000", "dis_conv_tol = 1.0e-10", "num_iter = 0     ! projection-only: SK-compatible character",
              "use_ws_distance = .false.   ! WS shifts are irrelevant here (A1 diagnostic: pbe_ws = pbe); scnfet.w90.apply_wsvec reads them if switched on",
              "write_hr = .true.", "write_xyz = .true.", "write_u_matrices = .true.", "bands_plot = .true.",
              "bands_num_points = 120", "", "begin projections", "Sc:d", "N:p", "end projections", "",
              "begin unit_cell_cart", "ang"] + [f"  {x:.8f} {y:.8f} {z:.8f}" for x, y, z in cell] + [
              "end unit_cell_cart", "", "begin atoms_frac", "Sc 0.0 0.0 0.0", "N  0.5 0.5 0.5", "end atoms_frac", "",
              "begin kpoint_path"]
    # W90 kpoint_path in fractional coordinates of the reciprocal basis
    import numpy as np
    A = np.array(cell)
    for (l1, c1), (l2, c2) in zip(PATH[:-1], PATH[1:]):
        f1 = A @ np.array(c1) / a
        f2 = A @ np.array(c2) / a
        lines.append(f"{l1} {f1[0]:.6f} {f1[1]:.6f} {f1[2]:.6f}  {l2} {f2[0]:.6f} {f2[1]:.6f} {f2[2]:.6f}")
    lines += ["end kpoint_path", "", f"mp_grid = {KW} {KW} {KW}", "", "begin kpoints"]
    lines += [f"  {x:.8f} {y:.8f} {z:.8f}" for x, y, z in kpts] + ["end kpoints", ""]
    return "\n".join(lines)


def open_grid_list(path):
    """The k list open_grid.x prints for the .win (crystal coordinates), in the order of its wavefunctions."""
    lines = Path(path).read_text().splitlines()
    i = next(j for j, l in enumerate(lines) if "List to be put in the .win file of wannier90" in l)
    ks = []
    for l in lines[i + 1:]:
        f = l.split()
        if len(f) != 4:
            break
        ks.append(tuple(float(v) for v in f[:3]))
    return ks


def stage2():
    """A1-2: HSE q 8³ IBZ (production, criterion 2 failed at q 4³) + open_grid, and Wannier90 for four cases:
    pbe (PBE twin, from pbe_nscf), hse_n (q4 route N), hse_o (q4 route O), hse_q8 (q8 route O, production).
    Windows from the A1-1 survey (a1_stage1.json): frozen ≤ VBM + 6.5 eV (HSE) / + 5.4 eV (PBE), just below the
    first non-model band (band 14 minimum: +7.02 / +5.88 eV); outer ≤ VBM + 13.5 / + 11.8 eV, below band 15."""
    import json
    st = json.loads((HERE / "a1_stage1.json").read_text())
    v_hse, v_pbe = st["survey"]["hse_nosym"]["vbm"], st["survey"]["pbe_nscf"]["vbm"]
    tags = []

    def w(tag, text):
        (OUT / f"{tag}.in").write_text(text)
        tags.append(tag)

    w("hse_ibz_q8", control("scf", "hse_ibz_q8") + system(FIXED, NBND, hse(8)) + electrons("1.0d-9", ", diago_full_acc=.true.")
      + SPECIES + f"K_POINTS automatic\n  {KW} {KW} {KW} 0 0 0\n")
    (OUT / "hse_ibz_q8_open.in").write_text("&inputpp\n  outdir='./tmp_hse_ibz_q8', prefix='hse_ibz_q8'\n/\n")
    tags.append("hse_ibz_q8_open")
    cases = {  # seed: (outdir, prefix, frozen max, outer max, k list or None = from open_grid in the job)
        "pbe": ("tmp_pbe_scf", "pbe_scf", v_pbe + 5.4, v_pbe + 11.8, full_grid()),
        "hse_n": ("tmp_hse_nosym", "hse_nosym", v_hse + 6.5, v_hse + 13.5, full_grid()),
        "hse_o": ("tmp_hse_ibz", "hse_ibz_open", v_hse + 6.5, v_hse + 13.5, open_grid_list(HERE / "outputs" / "hse_ibz_open.out")),
        "hse_q8": ("tmp_hse_ibz_q8", "hse_ibz_q8_open", v_hse + 6.5, v_hse + 13.5, None),
    }
    for seed, (od, pre, fz, wm, ks) in cases.items():
        d = HERE / "w90" / seed
        d.mkdir(parents=True, exist_ok=True)
        text = win(seed, ks if ks is not None else [], fz, wm)
        if ks is None:      # the job appends the open_grid list between these markers
            text = text.replace("begin kpoints\nend kpoints", "begin kpoints\n@KLIST@\nend kpoints")
        (d / f"{seed}.win").write_text(text)
        (d / f"{seed}.pw2wan").write_text(f"&inputpp\n  outdir='../../{od}', prefix='{pre}', seedname='{seed}',\n"
                                         f"  write_amn=.true., write_mmn=.true., write_unk=.false.\n/\n")
        tags.append(f"w90_{seed}")
    (HERE / "tags_w90.txt").write_text("\n".join(tags) + "\n")
    print("stage 2:", tags, f"| windows HSE frozen {v_hse + 6.5:.3f} outer {v_hse + 13.5:.3f}; PBE frozen {v_pbe + 5.4:.3f} outer {v_pbe + 11.8:.3f}")


if __name__ == "__main__":
    import sys
    stage2() if "stage2" in sys.argv else stage1()
