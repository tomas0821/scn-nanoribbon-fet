#!/usr/bin/env python3
"""A4b-1 inputs: the PBE twin and the edge Wannier recipe test on the four calibration ribbons (README.md).

Geometry: the A4a PAW relaxations, centred in the A4b vacuum (make_a4a.recentred): [100] n = 3 bare 15 Å; [100] n = 3 H,
[110] n = 4 bare and H 20 Å (decisions of 12–13 Nov). The 20 Å cells are exactly those of the A4a 20 Å checks.
Set: the NC set of A1/A2 (SG15 Sc, Dojo N, SG15 H), 80/320 Ry, gaussian 0.005 Ry (as A4a; the ribbons are gapped).

Per ribbon, one 64-rank job (run_a4b.sh <ribbon>), tags in this order:
  scf_<r>       PBE SCF, k 15 along x
  pot_<r>       pp.x bare + Hartree potential (plot_num 11) as a cube → vacuum level (window rule P1)
  bands_<r>     bands on the 31-k Γ–X̄ path (reference for P2–P5), up to ≈ VBM + 4 eV
  projwfc_<r>   Löwdin projections on the path → first non-model state (variant W2)
  nscf_<r>      NSCF on the full W90 grid (12 k [100], 16 k [110]), nosym/noinv, nbnd ≥ the outer window
  w90pp_<r>     wannier90.x -pp for both projection sets
  p2w_<r>_B     pw2wannier90, base set B (amn + mmn + eig)
  p2w_<r>_S     pw2wannier90, set S = B + Sc s on every Sc (diagnostic W3; amn only, mmn/eig linked from B)
The Wannier90 runs themselves (W1, W2 on B; W3 on S) follow in a second job once the windows are set from these outputs
(parse_a4b.py windows), as for A2.

    ~/.venvs/scnfet/bin/python -P dft/A4b/make_a4b.py          (from the project root)
    ~/.venvs/scnfet/bin/python -P dft/A4b/make_a4b.py hse      A4b-2: HSE06 inputs (tags_hse_<r>.txt, w90/<r>/H/)
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

PSEUDO_DIR = MK.PSEUDO_DIR
PP = {"Sc": ("44.955908", "Sc_ONCV_PBE-1.0.oncvpsp.upf"), "N": ("14.007", "N.oncvpsp.upf"),
      "H": ("1.008", "H_ONCV_PBE-1.0.oncvpsp.upf")}
ECUT, DUAL = 80, 4
SMEAR = MK.SMEAR
ZV = {"Sc": 11, "N": 5, "H": 1}
RIBBONS = ("r100_3", "r100_3h", "r110_4", "r110_4h")
VAC = {"r100_3": 15.0, "r100_3h": 20.0, "r110_4": 20.0, "r110_4h": 20.0}
NK_W90 = {"100": 12, "110": 16}                  # the sheet's density (16 per 3.11 Å), as cost_a4.kgrid
# bands run: up to ≈ VBM + 4 eV (A4a PAW path: 129 / 125 / 89 / 93 bands reach it at 15 / 15 / 15 / 20 Å), + room for
# the 20 Å vacuum states. NSCF: the outer window VBM + 10.23 eV (PBE twin, cost test T2) needs N(L) bands; measured
# 265 ([100] bare, 15 Å) and 201 ([110] H, 20 Å); estimated ≈ 265 ([100] H) and ≈ 200 ([110] bare) at 20 Å.
NB_BANDS = {"r100_3": 140, "r100_3h": 150, "r110_4": 105, "r110_4h": 105}
NB_NSCF = {"r100_3": 290, "r100_3h": 300, "r110_4": 225, "r110_4h": 225}
EDGE_S = {"r110_4"}                              # edge Sc get an s projection (Gate M: the Sc-s CB edge pair)
OUT = HERE / "inputs"
W90 = HERE / "w90"


def direction(key):
    return MK.ALL[key][0]


def geometry(key):
    """The A4a relaxed ribbon in the A4b cell; for r100_3 and r110_4h the relaxed cell already is it (shift 0)."""
    return MK.recentred(key, VAC[key])


def nexcl(st):
    """Sc 3s 3p (4 per Sc) + N 2s (1 per N): the lowest bands, below the model window (A2: 10 for 2 f.u.)."""
    return 4 * st.species.count("Sc") + st.species.count("N")


def nocc(st):
    return sum(ZV[s] for s in st.species) // 2


def control(calc, prefix, extra=""):
    return (f"&control\n  calculation='{calc}', prefix='{prefix}', outdir='./tmp_{prefix}', pseudo_dir='{PSEUDO_DIR}',\n"
            f"  verbosity='high', max_seconds=80000{extra}\n/\n")


def system(st, nbnd=None, extra=""):
    s = (f"&system\n  ibrav=0, nat={len(st.species)}, ntyp={len(set(st.species))},\n"
         f"  ecutwfc={ECUT}, ecutrho={ECUT * DUAL},\n  {SMEAR}")
    if nbnd:
        s += f", nbnd={nbnd}"
    return s + extra + "\n/\n"


def electrons(full_acc=False):
    return ("&electrons\n  conv_thr=1.0d-9, mixing_beta=0.3, mixing_mode='local-TF', electron_maxstep=300"
            + (", diago_full_acc=.true." if full_acc else "") + "\n/\n")


def cards(st, kpts):
    order = list(dict.fromkeys(st.species))
    s = "ATOMIC_SPECIES\n" + "".join(f"  {x:3s} {PP[x][0]} {PP[x][1]}\n" for x in order)
    s += "CELL_PARAMETERS angstrom\n" + "".join("  " + " ".join(f"{v:.10f}" for v in row) + "\n" for row in st.cell)
    s += "ATOMIC_POSITIONS angstrom\n" + "".join(f"  {x:3s} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}\n"
                                                 for x, p in zip(st.species, st.pos))
    return s + kpts


def grid(nk):
    return [(i / nk, 0.0, 0.0) for i in range(nk)]


def klist(nk):
    return f"K_POINTS crystal\n{nk}\n" + "".join(f"  {a:.10f} {b:.10f} {c:.10f} {1 / nk:.12f}\n" for a, b, c in grid(nk))


def projections(key, st, seed):
    """Projection block lines and the WF map [(atom index, species, orbital)] in W90 order (lines in order; within a
    species line, atoms in atoms-block order; d: z², xz, yz, x²−y², xy; p: z, x, y)."""
    sp = list(st.species)
    left, right = edge_rows(key, st)
    lines, wf = ["ang", "Sc:d", "N:p"], []
    wf += [(i, "Sc", o) for i, s in enumerate(sp) if s == "Sc" for o in ("dz2", "dxz", "dyz", "dx2-y2", "dxy")]
    wf += [(i, "N", o) for i, s in enumerate(sp) if s == "N" for o in ("pz", "px", "py")]
    if "H" in sp:
        lines.append("H:s")
        wf += [(i, "H", "s") for i, s in enumerate(sp) if s == "H"]
    if seed == "S":
        lines.append("Sc:s")
        wf += [(i, "Sc", "s") for i, s in enumerate(sp) if s == "Sc"]
    elif key in EDGE_S:
        for i, s in enumerate(sp):
            if s == "Sc" and (left[i] or right[i]):
                x, y, z = st.pos[i]
                lines.append(f"c={x:.8f},{y:.8f},{z:.8f}:s")
                wf.append((i, "Sc", "s"))
    return lines, wf


def edge_rows(key, st):
    """Outermost host rows (|y − y_edge| < 0.1 Å in the ideal cut, both planes), as parse_a4a.edge_masks; atom order
    is unchanged by the relaxation, so the masks apply to the relaxed structure."""
    left, right, _ = _edge_masks(key)
    assert len(left) == len(st.species)
    return left, right


def _edge_masks(key):
    st0 = MK.structure(key)
    sp = np.array(st0.species)
    host = sp != "H"
    y = st0.pos[:, 1]
    ylo, yhi = y[host].min(), y[host].max()
    return host & (np.abs(y - ylo) < 0.1), host & (np.abs(y - yhi) < 0.1), sp == "H"


def win(key, st, seed, froz=5.0, outer=30.0, nbnd=None, klist_placeholder=False):
    """The .win for one ribbon and projection set. Windows are placeholders for `-pp` (it ignores them); the variant
    runs get the real ones from parse_a4b.py windows. nbnd: total bands of the run (default the A4b-1 NSCF);
    klist_placeholder: @KLIST@ instead of the k list, filled on the cluster from open_grid.x (A4b-2, as A2)."""
    nk = NK_W90[direction(key)]
    lines_p, wf = projections(key, st, seed)
    ne = nexcl(st)
    L = [f"num_wann  = {len(wf)}", f"num_bands = {(nbnd or NB_NSCF[key]) - ne}", f"exclude_bands = 1-{ne}",
         f"dis_froz_max = {froz:.4f}", f"dis_win_max = {outer:.4f}", "dis_num_iter = 4000", "dis_conv_tol = 1.0e-10",
         "num_iter = 0     ! projection-only, as A1/A2", "use_ws_distance = .false.", "write_hr = .true.",
         "write_xyz = .true.", "", "begin projections", *lines_p, "end projections", "", "begin unit_cell_cart", "ang"]
    L += ["  " + " ".join(f"{v:.10f}" for v in row) for row in st.cell]
    L += ["end unit_cell_cart", "", "begin atoms_cart", "ang"]
    L += [f"{s:3s} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}" for s, p in zip(st.species, st.pos)]
    L += ["end atoms_cart", "", f"mp_grid = {nk} 1 1", "", "begin kpoints"]
    L += (["@KLIST@"] if klist_placeholder else [f"  {a:.10f} {b:.10f} {c:.10f}" for a, b, c in grid(nk)])
    L += ["end kpoints", ""]
    return "\n".join(L), wf


def write(tag, text):
    OUT.mkdir(exist_ok=True)
    (OUT / f"{tag}.in").write_text(text)
    return tag


def ribbon_inputs(key):
    st = geometry(key)
    nk = NK_W90[direction(key)]
    tags = [write(f"scf_{key}", control("scf", f"scf_{key}") + system(st) + electrons() + cards(st, MK.kauto(MK.K15)))]
    tags.append(write(f"pot_{key}",
                      f"&inputpp\n  prefix='scf_{key}', outdir='./tmp_scf_{key}', filplot='pot_{key}.dat', plot_num=11\n/\n"
                      f"&plot\n  nfile=1, filepp(1)='pot_{key}.dat', weight(1)=1.0, iflag=3, output_format=6, "
                      f"fileout='pot_{key}.cube'\n/\n"))
    tags.append(write(f"bands_{key}", control("bands", f"scf_{key}") + system(st, NB_BANDS[key]) + electrons(True)
                      + cards(st, MK.kpath())))
    tags.append(write(f"projwfc_{key}", f"&projwfc\n  prefix='scf_{key}', outdir='./tmp_scf_{key}',\n"
                                        f"  filpdos='pdos_{key}', degauss=0.005, DeltaE=0.01, kresolveddos=.true.\n/\n"))
    tags.append(write(f"nscf_{key}", control("nscf", f"scf_{key}")
                      + system(st, NB_NSCF[key], ",\n  nosym=.true., noinv=.true.") + electrons(True)
                      + cards(st, klist(nk))))
    d = W90 / key
    d.mkdir(parents=True, exist_ok=True)
    maps = {}
    for seed in ("B", "S"):
        text, wf = win(key, st, seed)
        (d / f"{key}_{seed}.win").write_text(text)
        (d / f"{key}_{seed}.pw2wan").write_text(
            f"&inputpp\n  outdir='../../tmp_scf_{key}', prefix='scf_{key}', seedname='{key}_{seed}',\n"
            f"  write_amn=.true., write_mmn={'.true.' if seed == 'B' else '.false.'}, write_unk=.false.\n/\n")
        maps[seed] = [dict(atom=a, species=s, orbital=o) for a, s, o in wf]
    (d / "wf_map.json").write_text(json.dumps(maps, indent=0))
    tags += [f"w90pp_{key}", f"p2w_{key}_B", f"p2w_{key}_S"]
    (HERE / f"tags_{key}.txt").write_text("\n".join(tags) + "\n")
    return st, maps


# ------------------------------------------------------------------------------------------------- A4b-2 (HSE06)
# README "HSE stage" and "A4b-2 staging". Recipe W2 for all four ribbons (A4b-1 decision, 13 Nov), 20 Å kept
# (decision, 13 Nov). HSE06 as A2: fixed occupations, IBZ SCF + open_grid.x, ecutfock 160, default divergence treatment.
NQ = {"100": 6, "110": 8}                        # README recipe item 1 (q along x; 1 × 1 across)
NK_IBZ = {"100": 7, "110": 9}                    # 12 / 16 k with time reversal: k = 0 … ½
HSE_OUTER = 11.05                                # cost test T2: HSE outer window, eV above the VBM
OUTER_MARGIN = 0.5
# nbnd: smallest N whose minimum over the A4b-1 PBE NSCF grid is ≥ VBM + 10.23 + 0.5 eV (parse_a4b.py hsurvey). The
# 0.5 eV margin covers vacuum-like states that HSE lifts less than the +0.82 eV the T2 rule assumes (A2 sheet, band
# minima: +0.69 to +1.40 eV); without it a few such states could cap the HSE outer window (H1 → rerun).
NB_HSE = {"r100_3": 267, "r100_3h": 273, "r110_4": 214, "r110_4h": 205}


def hse_inputs(key):
    st = geometry(key)
    d_ = direction(key)
    nk, nq = NK_W90[d_], NQ[d_]
    sysb = (f"&system\n  ibrav=0, nat={len(st.species)}, ntyp={len(set(st.species))},\n"
            f"  ecutwfc={ECUT}, ecutrho={ECUT * DUAL},\n  occupations='fixed', nbnd={NB_HSE[key]},\n"
            f"  input_dft='hse', nqx1={nq}, nqx2=1, nqx3=1, ecutfock={2 * ECUT}\n/\n")
    tags = [write(f"hse_{key}", control("scf", f"hse_{key}") + sysb + electrons(True)
                  + cards(st, f"K_POINTS automatic\n  {nk} 1 1 0 0 0\n"))]
    tags.append(write(f"hpot_{key}",
                      f"&inputpp\n  prefix='hse_{key}', outdir='./tmp_hse_{key}', filplot='hpot_{key}.dat', plot_num=11\n/\n"
                      f"&plot\n  nfile=1, filepp(1)='hpot_{key}.dat', weight(1)=1.0, iflag=3, output_format=6, "
                      f"fileout='hpot_{key}.cube'\n/\n"))
    tags.append(write(f"hproj_{key}", f"&projwfc\n  prefix='hse_{key}', outdir='./tmp_hse_{key}',\n"
                                      f"  filpdos='hpdos_{key}', degauss=0.005, DeltaE=0.01, kresolveddos=.true.\n/\n"))
    tags.append(write(f"hopen_{key}", f"&inputpp\n  outdir='./tmp_hse_{key}', prefix='hse_{key}'\n/\n"))
    d = W90 / key / "H"
    d.mkdir(parents=True, exist_ok=True)
    text, wf = win(key, st, "B", nbnd=NB_HSE[key], klist_placeholder=True)
    (d / f"{key}_H.win").write_text(text)
    (d / f"{key}_H.pw2wan").write_text(
        f"&inputpp\n  outdir='../../../tmp_hse_{key}', prefix='hse_{key}_open', seedname='{key}_H',\n"
        f"  write_amn=.true., write_mmn=.true., write_unk=.false.\n/\n")
    tags += [f"hw90pp_{key}", f"hp2w_{key}"]
    (HERE / f"tags_hse_{key}.txt").write_text("\n".join(tags) + "\n")
    return st, wf


def hse_set_s(key):
    """Amendment 2026-11-14: set S (B + Sc s on every Sc) on the HSE run, for a ribbon whose H3 fired. The k list is
    taken from the cluster-filled w90/<r>/H/<r>_H.win (open_grid order); .mmn/.eig are linked from <r>_H."""
    st = geometry(key)
    d = W90 / key / "H"
    filled = (d / f"{key}_H.win").read_text()
    assert "@KLIST@" not in filled, "copy back the filled <r>_H.win first"
    kl = filled.split("begin kpoints\n")[1].split("end kpoints")[0].rstrip("\n")
    text, wf = win(key, st, "S", nbnd=NB_HSE[key], klist_placeholder=True)
    (d / f"{key}_HS.win").write_text(text.replace("@KLIST@", kl))
    (d / f"{key}_HS.pw2wan").write_text(
        f"&inputpp\n  outdir='../../../tmp_hse_{key}', prefix='hse_{key}_open', seedname='{key}_HS',\n"
        f"  write_amn=.true., write_mmn=.false., write_unk=.false.\n/\n")
    return len(wf)


def main_hse():
    info = json.loads((HERE / "a4b_inputs.json").read_text())
    for key in RIBBONS:
        st, wf = hse_inputs(key)
        d_ = direction(key)
        assert len(wf) == info[key]["num_wann_B"]
        info[key].update(nbnd_hse=NB_HSE[key], nq_hse=NQ[d_], nk_ibz_hse=NK_IBZ[d_])
        print(f"{key}: HSE nbnd {NB_HSE[key]}, k {NK_W90[d_]} (IBZ {NK_IBZ[d_]}), q {NQ[d_]}, num_wann {len(wf)}, "
              f"num_bands {NB_HSE[key] - nexcl(st)}")
    (HERE / "a4b_inputs.json").write_text(json.dumps(info, indent=1))


def main():
    info = {}
    for key in RIBBONS:
        st, maps = ribbon_inputs(key)
        ne, no = nexcl(st), nocc(st)
        info[key] = dict(vacuum=VAC[key], nat=len(st.species), nocc=no, nexcl=ne, nk_w90=NK_W90[direction(key)],
                         nbnd_bands=NB_BANDS[key], nbnd_nscf=NB_NSCF[key], num_wann_B=len(maps["B"]),
                         num_wann_S=len(maps["S"]), cell=np.diag(st.cell).round(4).tolist())
        print(f"{key}: {len(st.species)} at, cell {np.diag(st.cell).round(3)}, nocc {no}, exclude 1-{ne}, "
              f"W90 k {NK_W90[direction(key)]}, num_wann B {len(maps['B'])} / S {len(maps['S'])}, "
              f"num_bands {NB_NSCF[key] - ne}, occupied in window {no - ne}")
    (HERE / "a4b_inputs.json").write_text(json.dumps(info, indent=1))


if __name__ == "__main__":
    if sys.argv[1:2] == ["hse_s"]:
        for k in sys.argv[2:]:
            print(f"{k}: set S on the HSE run, num_wann {hse_set_s(k)}")
    else:
        main_hse() if sys.argv[1:] == ["hse"] else main()
