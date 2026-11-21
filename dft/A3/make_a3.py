#!/usr/bin/env python3
"""A3 inputs: the 3-plane slab hold-out (README.md, written before any A3 input existed).

Stage A3-1 (`make_a3.py`): one PAW PBE relaxation, internal coordinates at fixed a_p and c, as A4a / A5 otherwise
(make_a4a.py: SSSP precision PAW, 90/720 Ry, gaussian 0.005 Ry, BFGS, forc 1e-3, etot 1e-5), k 12 × 12 × 1 (A0's
sheet grid), nspin 1.

Stage A3-2 (`make_a3.py stage2`, after relax_slab3 is copied back): NC PBE + HSE06 as A2, at the relaxed geometry
(NC set and cards from make_a4b; A2's smearing m-v 0.005 Ry, conv_thr 1e-10 / 1e-9, mixing 0.4; nbnd 84 = A2's 28 per
f.u.). Three jobs, all in scn_fet/A3:
  PBE twin (tags_slab3.txt; A4b-1's run_a4b.sh unchanged, A4B_DIR=…/A3, R = slab3; the tag names are A4b-1's):
    scf_slab3 (16 × 16 × 1), pot_slab3 (pp.x vacuum level), bands_slab3 (Γ–X̄–M̄–Γ, 60 per segment, A2's path),
    projwfc_slab3 (on the path: the W2 character test), nscf_slab3 (the full 16 × 16 × 1 list, nosym), w90pp_slab3,
    p2w_slab3_B (Sc d + N p), p2w_slab3_S (+ Sc s on every Sc; amn only)
  HSE06 (tags_hse_slab3.txt; run_a3_hse.sh = run_a5_hse.sh's three steps + open_grid and Wannier90 as run_a4b2.sh):
    hse_slab3 (IBZ 16 × 16 × 1, q 8 × 8 × 1, fixed occupations, ecutfock 160), hpot_slab3, hproj_slab3, hopen_slab3,
    hw90pp_slab3 (open_grid k list into the H / HS .win), hp2w_slab3 (H: amn mmn eig; HS: amn)
  vacuum check (tags_vac.txt; run_a3.sh): vac_slab3, the relaxed slab centred in 27.5 Å vacuum, PBE 12 × 12 × 1
The Wannier90 variants (A3-3) follow from the windows.

    ~/.venvs/scnfet/bin/python -P dft/A3/make_a3.py [stage2]          (from the project root)
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
_spec = importlib.util.spec_from_file_location("make_a4a", ROOT / "dft" / "A4a" / "make_a4a.py")
MK = importlib.util.module_from_spec(_spec)
_argv, sys.argv = sys.argv, ["make_a4a.py"]
_spec.loader.exec_module(MK)
sys.argv = _argv
from scnfet import structures as S          # noqa: E402

NPLANES, VAC, K_RELAX = 3, 20.0, 12
OUT = HERE / "inputs"


def structure():
    return S.slab(NPLANES, a_p=MK.A_P, d=MK.D_PLANE, vacuum=VAC)


def write(tag, text):
    OUT.mkdir(exist_ok=True)
    (OUT / f"{tag}.in").write_text(text)
    return tag


def relax_input():
    st = structure()
    body, order = MK.cards(st, list(st.species), f"K_POINTS automatic\n  {K_RELAX} {K_RELAX} 1 0 0 0\n")
    return (MK.control("relax", "relax_slab3", ",\n  forc_conv_thr=1.0d-3, etot_conv_thr=1.0d-5, nstep=300")
            + MK.system(len(st.species), len(order)) + MK.ELECTRONS + "&ions\n  ion_dynamics='bfgs'\n/\n" + body)


def stage1():
    st = structure()
    tag = write("relax_slab3", relax_input())
    (HERE / "tags_relax.txt").write_text(tag + "\n")
    info = dict(nplanes=NPLANES, vacuum=VAC, nat=len(st.species), species=list(st.species),
                cell=[round(float(x), 4) for x in np.diag(st.cell)], k_relax=K_RELAX)
    (HERE / "a3_inputs.json").write_text(json.dumps(info, indent=1))
    print(f"{tag}: {len(st.species)} atoms {st.species}, cell {np.diag(st.cell).round(4)}, z {st.pos[:, 2].round(4)}")


# ------------------------------------------------------------------------------------------------- A3-2
import re                                    # noqa: E402

KEY, NBND, KW, NQ, VAC_CHECK, NBND_VAC = "slab3", 84, 16, 8, 27.5, 36
PATH = [("G", (0.0, 0.0, 0.0)), ("X", (0.5, 0.0, 0.0)), ("M", (0.5, 0.5, 0.0)), ("G", (0.0, 0.0, 0.0))]   # A2
SMEAR_A2, FIXED = "occupations='smearing', smearing='mv', degauss=0.005", "occupations='fixed'"
W90 = HERE / "w90"


def _mb():
    spec = importlib.util.spec_from_file_location("make_a4b", ROOT / "dft" / "A4b" / "make_a4b.py")
    mb = importlib.util.module_from_spec(spec)
    argv, sys.argv[:] = list(sys.argv), ["make_a4b.py"]
    spec.loader.exec_module(mb)
    sys.argv[:] = argv
    return mb


def relaxed():
    """Final positions from outputs/relax_slab3.out (as make_a4a.relaxed)."""
    txt = (HERE / "outputs" / "relax_slab3.out").read_text()
    assert "bfgs converged" in txt and "End final coordinates" in txt, "relax_slab3 not converged"
    blk = txt[txt.rindex("Begin final coordinates"):txt.rindex("End final coordinates")]
    rows = re.findall(r"^\s*([A-Za-z]+\d*)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", blk.split("ATOMIC_POSITIONS")[1], re.M)
    st = structure()
    pos = np.array([[float(v) for v in r[1:]] for r in rows])
    assert len(pos) == len(st.species) and [r[0] for r in rows] == list(st.species)
    return S.Structure(st.cell, list(st.species), pos, st.pbc, st.frame)


def control(calc, prefix, extra=""):
    return (f"&control\n  calculation='{calc}', prefix='{prefix}', outdir='./tmp_{prefix}', pseudo_dir='{MK.PSEUDO_DIR}',\n"
            f"  verbosity='high', max_seconds=80000{extra}\n/\n")


def system(st, occ, nbnd, extra=""):
    return (f"&system\n  ibrav=0, nat={len(st.species)}, ntyp={len(set(st.species))},\n  ecutwfc=80, ecutrho=320,\n"
            f"  {occ}, nbnd={nbnd}{extra}\n/\n")


def electrons(thr="1.0d-10", full_acc=False):
    return (f"&electrons\n  conv_thr={thr}, mixing_beta=0.4" + (", diago_full_acc=.true." if full_acc else "") + "\n/\n")


def full_grid(n=KW):
    return [(i / n, j / n, 0.0) for i in range(n) for j in range(n)]


def klist(n=KW):
    ks = full_grid(n)
    return f"K_POINTS crystal\n{len(ks)}\n" + "".join(f"  {a:.8f} {b:.8f} {c:.8f} {1 / len(ks):.10f}\n" for a, b, c in ks)


def kpath():
    seg = "".join(f"  {x:.4f} {y:.4f} {z:.4f} {60 if i < len(PATH) - 1 else 1}\n" for i, (_, (x, y, z)) in enumerate(PATH))
    return f"K_POINTS crystal_b\n  {len(PATH)}\n" + seg


def nexcl(st):
    return 4 * st.species.count("Sc") + st.species.count("N")


def win(st, seed, klines):
    """The slab .win, as make_a4b.win in 2D: projection-only, Sc 3s 3p + N 2s excluded; placeholder windows (the -pp
    run ignores them; A3-3 sets the real ones). seed B: Sc d + N p; S: + Sc s on every Sc."""
    proj = ["Sc:d", "N:p"] + (["Sc:s"] if seed in ("S", "HS") else [])
    nw = 5 * st.species.count("Sc") + 3 * st.species.count("N") + (st.species.count("Sc") if seed in ("S", "HS") else 0)
    ne = nexcl(st)
    L = [f"num_wann  = {nw}", f"num_bands = {NBND - ne}", f"exclude_bands = 1-{ne}", "dis_froz_max = 5.0000",
         "dis_win_max = 30.0000", "dis_num_iter = 4000", "dis_conv_tol = 1.0e-10", "num_iter = 0     ! projection-only, as A1/A2",
         "use_ws_distance = .false.", "write_hr = .true.", "write_xyz = .true.", "", "begin projections", *proj,
         "end projections", "", "begin unit_cell_cart", "ang"]
    L += ["  " + " ".join(f"{v:.10f}" for v in row) for row in st.cell]
    L += ["end unit_cell_cart", "", "begin atoms_cart", "ang"]
    L += [f"{x:3s} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}" for x, p in zip(st.species, st.pos)]
    L += ["end atoms_cart", "", f"mp_grid = {KW} {KW} 1", "", "begin kpoints", *klines, "end kpoints", ""]
    return "\n".join(L), nw


def stage2():
    mb = _mb()
    st = relaxed()
    cards = mb.cards
    pbe = []

    def w(tag, text, lst):
        write(tag, text)
        lst.append(tag)

    k16 = f"K_POINTS automatic\n  {KW} {KW} 1 0 0 0\n"
    w(f"scf_{KEY}", control("scf", f"scf_{KEY}", ", tstress=.true., tprnfor=.true.") + system(st, SMEAR_A2, NBND)
      + electrons() + cards(st, k16), pbe)
    w(f"pot_{KEY}", f"&inputpp\n  prefix='scf_{KEY}', outdir='./tmp_scf_{KEY}', filplot='pot_{KEY}.dat', plot_num=11\n/\n"
                    f"&plot\n  nfile=1, filepp(1)='pot_{KEY}.dat', weight(1)=1.0, iflag=3, output_format=6, "
                    f"fileout='pot_{KEY}.cube'\n/\n", pbe)
    w(f"bands_{KEY}", control("bands", f"scf_{KEY}") + system(st, SMEAR_A2, NBND) + electrons(full_acc=True)
      + cards(st, kpath()), pbe)
    w(f"projwfc_{KEY}", f"&projwfc\n  prefix='scf_{KEY}', outdir='./tmp_scf_{KEY}',\n"
                        f"  filpdos='pdos_{KEY}', degauss=0.005, DeltaE=0.01, kresolveddos=.true.\n/\n", pbe)
    w(f"nscf_{KEY}", control("nscf", f"scf_{KEY}") + system(st, SMEAR_A2, NBND, ",\n  nosym=.true., noinv=.true.")
      + electrons(full_acc=True) + cards(st, klist()), pbe)
    d = W90 / KEY
    d.mkdir(parents=True, exist_ok=True)
    kl = [f"  {a:.8f} {b:.8f} {c:.8f}" for a, b, c in full_grid()]
    nw = {}
    for seed in ("B", "S"):
        text, nw[seed] = win(st, seed, kl)
        (d / f"{KEY}_{seed}.win").write_text(text)
        (d / f"{KEY}_{seed}.pw2wan").write_text(
            f"&inputpp\n  outdir='../../tmp_scf_{KEY}', prefix='scf_{KEY}', seedname='{KEY}_{seed}',\n"
            f"  write_amn=.true., write_mmn={'.true.' if seed == 'B' else '.false.'}, write_unk=.false.\n/\n")
    pbe += [f"w90pp_{KEY}", f"p2w_{KEY}_B", f"p2w_{KEY}_S"]
    (HERE / f"tags_{KEY}.txt").write_text("\n".join(pbe) + "\n")
    # HSE06
    hse = []
    w(f"hse_{KEY}", control("scf", f"hse_{KEY}") + system(st, FIXED, NBND, f",\n  input_dft='hse', nqx1={NQ}, nqx2={NQ}, "
                                                                        f"nqx3=1, ecutfock=160")
      + electrons("1.0d-9", full_acc=True) + cards(st, k16), hse)
    w(f"hpot_{KEY}", f"&inputpp\n  prefix='hse_{KEY}', outdir='./tmp_hse_{KEY}', filplot='hpot_{KEY}.dat', plot_num=11\n/\n"
                     f"&plot\n  nfile=1, filepp(1)='hpot_{KEY}.dat', weight(1)=1.0, iflag=3, output_format=6, "
                     f"fileout='hpot_{KEY}.cube'\n/\n", hse)
    w(f"hproj_{KEY}", f"&projwfc\n  prefix='hse_{KEY}', outdir='./tmp_hse_{KEY}',\n"
                      f"  filpdos='hpdos_{KEY}', degauss=0.005, DeltaE=0.01, kresolveddos=.true.\n/\n", hse)
    w(f"hopen_{KEY}", f"&inputpp\n  outdir='./tmp_hse_{KEY}', prefix='hse_{KEY}'\n/\n", hse)
    dh = d / "H"
    dh.mkdir(exist_ok=True)
    for seed in ("H", "HS"):
        text, _ = win(st, seed, ["@KLIST@"])
        (dh / f"{KEY}_{seed}.win").write_text(text)
        (dh / f"{KEY}_{seed}.pw2wan").write_text(
            f"&inputpp\n  outdir='../../../tmp_hse_{KEY}', prefix='hse_{KEY}_open', seedname='{KEY}_{seed}',\n"
            f"  write_amn=.true., write_mmn={'.true.' if seed == 'H' else '.false.'}, write_unk=.false.\n/\n")
    hse += [f"hw90pp_{KEY}", f"hp2w_{KEY}"]
    (HERE / f"tags_hse_{KEY}.txt").write_text("\n".join(hse) + "\n")
    # vacuum check: the relaxed slab centred in a cell with 27.5 Å vacuum
    big = S.slab(NPLANES, a_p=MK.A_P, d=MK.D_PLANE, vacuum=VAC_CHECK)
    shift = np.array([0.0, 0.0, (big.cell[2, 2] - st.cell[2, 2]) / 2])
    stv = S.Structure(big.cell, list(st.species), st.pos + shift, st.pbc, st.frame)
    write(f"vac_{KEY}", control("scf", f"vac_{KEY}") + system(stv, SMEAR_A2, NBND_VAC) + electrons()
          + cards(stv, f"K_POINTS automatic\n  12 12 1 0 0 0\n"))
    (HERE / "tags_vac.txt").write_text(f"vac_{KEY}\n")
    z = st.pos[:, 2]
    info = json.loads((HERE / "a3_inputs.json").read_text())
    info["relaxed"] = dict(z=[round(float(v), 4) for v in z], species=list(st.species),
                           spacing_outer_middle=[round(float(z[2] - z[0]), 4), round(float(z[3] - z[1]), 4)],
                           rumpling_outer_N_minus_Sc=round(float(z[0] - z[1]), 4))
    info["stage2"] = dict(nbnd=NBND, nexcl=nexcl(st), num_wann_B=nw["B"], num_wann_S=nw["S"], k=KW, nq=NQ,
                          vac_check=VAC_CHECK, c_vac=round(float(big.cell[2, 2]), 4), nbnd_vac=NBND_VAC,
                          tags_pbe=pbe, tags_hse=hse)
    (HERE / "a3_inputs.json").write_text(json.dumps(info, indent=1))
    print(f"PBE twin: {', '.join(pbe)}\nHSE: {', '.join(hse)}\nvacuum: vac_{KEY} (c {big.cell[2, 2]:.3f} Å)")
    print(f"nbnd {NBND}, excluded {nexcl(st)}, num_wann B {nw['B']} / S {nw['S']}; relaxed z {z.round(4)}")


if __name__ == "__main__":
    if "stage2" in sys.argv:
        stage2()
    else:
        stage1()
