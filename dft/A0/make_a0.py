#!/usr/bin/env python3
"""Generate the A0 cutoff / pseudopotential test inputs (plan §4, A0; the project DFT settings).

What A0 answers:
  * the PAW set (SSSP-precision 1.3.0: Sc.pbe-spn-kjpaw_psl.0.2.3 + N.oncvpsp Dojo), at ecutwfc 60/75/90/110 Ry and
    ecutrho = 8x. This is the PBE relaxation set;
  * the all-NC set (Sc_ONCV_PBE-1.0 SG15 + N.oncvpsp Dojo), at 80/90 Ry plus a 110 Ry reference and
    ecutrho = 4x. This is the HSE06 set; the plan fixes HSE at 80 Ry with ecutfock 160;
  * rocksalt ScN: lattice constant from an E(V) fit, and the sign and size of the PBE Gamma-X gap
    (target about -0.04 eV);
  * the M-ScN sheet (P4/nmm, 2 f.u.): in-plane a from a 2Dxy vc-relax (PAW only, since only PAW relaxes), and
    the PBE gap at a fixed geometry for both sets (target 0.388 eV, Alexandria agm2000022730). The same geometry
    gives a PAW-vs-NC consistency check.

Run it here; it writes inputs/*.in plus tags.txt. Copy dft/A0/ to the cluster and sbatch run_a0.sh.
"""
from pathlib import Path

BOHR = 0.529177210903
PSEUDO_DIR = "/path/to/scn_fet/pseudo"
MASS = {"Sc": 44.955908, "N": 14.007}
SETS = {  # name: (Sc file, N file, ecutwfc, dual)
    **{f"paw{e}": ("Sc.pbe-spn-kjpaw_psl.0.2.3.UPF", "N.oncvpsp.upf", e, 8) for e in (60, 75, 90, 110)},
    **{f"nc{e}": ("Sc_ONCV_PBE-1.0.oncvpsp.upf", "N.oncvpsp.upf", e, 4) for e in (80, 90, 110)},
}
A_BULK = [4.42, 4.45, 4.48, 4.51, 4.54, 4.57, 4.60]   # Å, centred on the PBE literature value (≈ 4.51-4.52)
A_FIX_BULK = 4.50                                    # Å, fixed-a gap comparisons (the plan's a ≈ 4.50 Å)
SHEET_A, SHEET_D, SHEET_C = 3.11, 2.18, 22.5         # Å: Liang SI S1 (4.40/sqrt2; interlayer 2.18); vacuum ≈ 20.3 Å
K_BULK, K_SHEET = 12, 12                             # even Γ-centred meshes contain Γ, X (bulk) and Γ, X̄, M̄ (sheet)

out = Path(__file__).resolve().parent / "inputs"
out.mkdir(exist_ok=True)
tags = []


def species(s):
    sc, n, *_ = SETS[s]
    return f"ATOMIC_SPECIES\n  Sc {MASS['Sc']} {sc}\n  N  {MASS['N']} {n}\n"


def control(calc, tag, extra=""):
    return (f"&control\n  calculation='{calc}', prefix='{tag}', outdir='./tmp_{tag}', pseudo_dir='{PSEUDO_DIR}',\n"
            f"  tstress=.true., tprnfor=.true., verbosity='high', max_seconds=80000{extra}\n/\n")


def system(s, ibrav, celldm, nat, nbnd, degauss):
    _, _, e, dual = SETS[s]
    cd = ", ".join(f"celldm({i})={v:.10f}" for i, v in celldm.items())
    return (f"&system\n  ibrav={ibrav}, {cd}, nat={nat}, ntyp=2,\n  ecutwfc={e}, ecutrho={e * dual},\n"
            f"  occupations='smearing', smearing='mv', degauss={degauss}, nbnd={nbnd}\n/\n")


def write(tag, text):
    (out / f"{tag}.in").write_text(text)
    tags.append(tag)


def bulk(s, a, k=K_BULK, tag=None):
    tag = tag or f"bulk_{s}_a{a:.2f}"
    write(tag, control("scf", tag) + system(s, 2, {1: a / BOHR}, 2, 16, 0.01)
          + "&electrons\n  conv_thr=1.0d-10, mixing_beta=0.5\n/\n" + species(s)
          + "ATOMIC_POSITIONS crystal\n  Sc 0.0 0.0 0.0\n  N  0.5 0.5 0.5\n"
          + f"K_POINTS automatic\n  {k} {k} {k} 0 0 0\n")


def sheet_cell(a):
    z1, z2 = 0.5 - SHEET_D / (2 * SHEET_C), 0.5 + SHEET_D / (2 * SHEET_C)
    # rocksalt-(001) bilayer: Sc over N between the planes; P4/nmm, centrosymmetric, so no dipole correction is needed
    pos = (f"ATOMIC_POSITIONS crystal\n  Sc 0.0 0.0 {z1:.8f}\n  N  0.5 0.5 {z1:.8f}\n"
           f"  N  0.0 0.0 {z2:.8f}\n  Sc 0.5 0.5 {z2:.8f}\n")
    return {1: a / BOHR, 3: SHEET_C / a}, pos


def sheet_scf(s):
    tag = f"sheet_{s}_fixed"
    cd, pos = sheet_cell(SHEET_A)
    write(tag, control("scf", tag) + system(s, 6, cd, 4, 24, 0.005)
          + "&electrons\n  conv_thr=1.0d-10, mixing_beta=0.4\n/\n" + species(s) + pos
          + f"K_POINTS automatic\n  {K_SHEET} {K_SHEET} 1 0 0 0\n")
    return tag


def sheet_bands(s, scf_tag):
    # Γ–X̄–M̄–Γ, 30 points per segment; reuses the scf save (same prefix/outdir)
    tag = f"{scf_tag}_bands"
    cd, pos = sheet_cell(SHEET_A)
    ctrl = (f"&control\n  calculation='bands', prefix='{scf_tag}', outdir='./tmp_{scf_tag}', "
            f"pseudo_dir='{PSEUDO_DIR}', verbosity='high'\n/\n")
    write(tag, ctrl + system(s, 6, cd, 4, 24, 0.005) + "&electrons\n  diago_full_acc=.true.\n/\n"
          + species(s) + pos
          + "K_POINTS tpiba_b\n  4\n  0.0 0.0 0.0 30\n  0.5 0.0 0.0 30\n  0.5 0.5 0.0 30\n  0.0 0.0 0.0 1\n")


def sheet_relax(s):
    tag = f"sheet_{s}_vcrelax"
    cd, pos = sheet_cell(SHEET_A)
    write(tag, control("vc-relax", tag, ",\n  forc_conv_thr=1.0d-4, etot_conv_thr=1.0d-6, nstep=150")
          + system(s, 6, cd, 4, 24, 0.005)
          + "&electrons\n  conv_thr=1.0d-10, mixing_beta=0.4\n/\n&ions\n  ion_dynamics='bfgs'\n/\n"
          + "&cell\n  cell_dynamics='bfgs', cell_dofree='2Dxy', press_conv_thr=0.2\n/\n" + species(s) + pos
          + f"K_POINTS automatic\n  {K_SHEET} {K_SHEET} 1 0 0 0\n")


# 1. bulk E(V) for every set (cheap: 2 atoms)
for s in SETS:
    for a in A_BULK:
        bulk(s, a)
# 2. bulk at the fixed a = 4.50 Å, with a k-convergence ladder for the production cutoffs
for s in SETS:
    bulk(s, A_FIX_BULK)
for s in ("paw90", "nc80"):
    for k in (8, 16, 20):
        bulk(s, A_FIX_BULK, k=k, tag=f"bulk_{s}_a{A_FIX_BULK:.2f}_k{k}")
# 3. sheet at a fixed geometry, every set; bands for the two production sets
for s in SETS:
    t = sheet_scf(s)
    if s in ("paw90", "nc80"):
        sheet_bands(s, t)
# 4. sheet 2Dxy vc-relax, PAW only (relaxations are PAW by design)
for s in ("paw60", "paw75", "paw90", "paw110"):
    sheet_relax(s)

(out.parent / "tags.txt").write_text("\n".join(tags) + "\n")
print(f"{len(tags)} inputs -> {out}")
