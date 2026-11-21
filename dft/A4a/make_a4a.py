#!/usr/bin/env python3
"""A4a inputs: the calibration ribbons and Gate M (plan §4 A4a, §3 D4; criteria in README.md).

Ribbons (scnfet.structures.ribbon; W edge to edge between nuclei):
  r100_3   [100] n = 3 bare, 24 atoms, W 1.10 nm   (non-polar edges, mixed Sc/N rows)
  r110_4   [110] n = 4 bare, 16 atoms, W 1.09 nm   (polar edges, opposite polarity in the two planes)
  r100_3h  [100] n = 3 + 8 H (one H per edge atom along the edge normal), W 1.47 nm: the mandatory H-terminated
           ribbon. [100] because each edge atom has exactly one cut bond; the [110] H termination stays open.

Geometry convention = A2: the A0 PBE sheet scaled by the bulk a_exp / a_PBE (a_p 3.1093 Å, planes 2.1685 Å), so the
[110] period is a_p and the [100] period a_p √2. A4a relaxes the internal coordinates at that fixed period; the relaxed
structures are then used unchanged for A4b (HSE-W90). 15 Å vacuum in y and z (checked by `vac20`).

Settings: PBE, SSSP precision as A0 (Sc PAW psl 0.2.3 spn, N Dojo, H SG15), 90/720 Ry (A0 rule); gaussian smearing
0.005 Ry for every run, so NM and magnetic energies compare; k 15 along x (Γ-centred, 8 irreducible); k-test at 30.

Stage A4a-1 (this script, `python3 make_a4a.py`): three relaxations, one 64-rank job.
  relax_r100_3, relax_r110_4   relax from a ferromagnetic start (nspin 2): the moment survives or collapses
  relax_r100_3h                relax nspin 1
Stage A4a-2 (`python3 make_a4a.py stage2`, after A4a-1 is copied back to outputs/): single points at the relaxed
geometries: NM, FM, AFM across the edges, AFM between the planes (bare); NM + FM check (H); k 30 for NM, FM and AFM
edges (bare); bands Γ–X̄ (31 k, 28 bands per f.u.: also the band survey for A4b's nbnd) and projwfc for NM; vacuum check
(r110_4 NM at 20 Å).
"""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from scnfet import structures as S          # noqa: E402

PSEUDO_DIR = "/path/to/scn_fet/pseudo"
PP = {"Sc": ("44.955908", "Sc.pbe-spn-kjpaw_psl.0.2.3.UPF"), "N": ("14.007", "N.oncvpsp.upf"),
      "H": ("1.008", "H_ONCV_PBE-1.0.oncvpsp.upf")}
A_P = 3.117899 * 4.501 / 4.5134
D_PLANE = 2.1685
VAC = 15.0
ECUT, DUAL = 90, 8
K15, K30 = 15, 30
SMEAR = "occupations='smearing', smearing='gaussian', degauss=0.005"
MAG0 = {"Sc": 0.10, "N": 0.20, "H": 0.0}       # starting_magnetization: ≈ 1.1 μB per Sc, 1.0 per N
RIBBONS = {"r100_3": ("100", 3, False), "r110_4": ("110", 4, False), "r100_3h": ("100", 3, True)}
# Added after Gate M (decisions 2026-11-12): the [110] H ribbon, one H per edge atom on the edge normal, at the
# [110] vacuum (20 Å, criterion 5 failed for [110]). Kept out of RIBBONS so the stage-1/2 loops stay as run.
EXTRA = {"r110_4h": ("110", 4, True)}
VAC_OF = {"r110_4h": 20.0}
ALL = {**RIBBONS, **EXTRA}
OUT = HERE / "inputs"


def structure(key, vacuum=None):
    d, n, h = ALL[key]
    return S.ribbon(d, n, a_p=A_P, d=D_PLANE, vacuum=vacuum or VAC_OF.get(key, VAC), passivate=h)


def control(calc, tag, extra=""):
    return (f"&control\n  calculation='{calc}', prefix='{tag}', outdir='./tmp_{tag}', pseudo_dir='{PSEUDO_DIR}',\n"
            f"  verbosity='high',{' tprnfor=.true.,' if calc in ('relax', 'scf') else ''} max_seconds=82000{extra}\n/\n")


def system(nat, ntyp, spin=None, nbnd=None):
    """spin: None (nspin 1) or {species_label: starting_magnetization}."""
    s = (f"&system\n  ibrav=0, nat={nat}, ntyp={ntyp},\n  ecutwfc={ECUT}, ecutrho={ECUT * DUAL},\n  {SMEAR}")
    if nbnd:
        s += f", nbnd={nbnd}"
    if spin:
        s += ",\n  nspin=2, " + ", ".join(f"starting_magnetization({i + 1})={m:+.2f}" for i, m in enumerate(spin.values()))
    return s + "\n/\n"


ELECTRONS = "&electrons\n  conv_thr=1.0d-8, mixing_beta=0.3, mixing_mode='local-TF', electron_maxstep=300\n/\n"


def cards(st, labels, kpts):
    """labels: per-atom species label (may be split, e.g. 'Sc1'); ATOMIC_SPECIES lists the labels in first-seen order."""
    order = list(dict.fromkeys(labels))
    s = "ATOMIC_SPECIES\n" + "".join(f"  {l:4s} {PP[re.sub(r'[0-9]', '', l)][0]} {PP[re.sub(r'[0-9]', '', l)][1]}\n"
                                      for l in order)
    s += "CELL_PARAMETERS angstrom\n" + "".join("  " + " ".join(f"{v:.10f}" for v in row) + "\n" for row in st.cell)
    s += "ATOMIC_POSITIONS angstrom\n" + "".join(f"  {l:4s} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}\n"
                                                 for l, p in zip(labels, st.pos))
    s += kpts
    return s, order


def kauto(n):
    return f"K_POINTS automatic\n  {n} 1 1 0 0 0\n"


def kpath(n=31):
    return "K_POINTS crystal_b\n2\n  0.0 0.0 0.0 30\n  0.5 0.0 0.0 1\n"


def write(tag, text):
    OUT.mkdir(exist_ok=True)
    (OUT / f"{tag}.in").write_text(text)
    return tag


def spin_labels(st, mode):
    """Species labels and starting magnetisations for a spin configuration.
    fm: all up. afm_edges: atoms left of the ribbon centre (y) up, right down. afm_planes: lower plane up, upper down."""
    sp, pos = list(st.species), st.pos
    yc, zc = pos[:, 1].mean(), pos[:, 2].mean()
    labels, mag = [], {}
    for s, p in zip(sp, pos):
        if mode == "fm" or s == "H":
            lab, sign = s, 1
        elif mode == "afm_edges":
            lab, sign = (s + "1", 1) if p[1] < yc else (s + "2", -1)
        elif mode == "afm_planes":
            lab, sign = (s + "1", 1) if p[2] < zc else (s + "2", -1)
        else:
            raise ValueError(mode)
        labels.append(lab)
        mag[lab] = sign * MAG0[s]
    order = list(dict.fromkeys(labels))
    return labels, {l: mag[l] for l in order}


def relax_input(key):
    st = structure(key)
    bare = not ALL[key][2]
    if bare:
        labels, mag = spin_labels(st, "fm")
    else:
        labels, mag = list(st.species), None
    body, order = cards(st, labels, kauto(K15))
    return (control("relax", f"relax_{key}", ",\n  forc_conv_thr=1.0d-3, etot_conv_thr=1.0d-5, nstep=200")
            + system(len(st.species), len(order), mag) + ELECTRONS + "&ions\n  ion_dynamics='bfgs'\n/\n" + body)


def stage1():
    tags = [write(f"relax_{k}", relax_input(k)) for k in RIBBONS]
    (HERE / "tags_stage1.txt").write_text("\n".join(tags) + "\n")
    for k in RIBBONS:
        st = structure(k)
        print(f"{k}: {len(st.species)} atoms, W {S.width(st) / 10:.2f} nm, cell {np.diag(st.cell).round(3)}")
    print("wrote", ", ".join(tags))


def relaxed(key):
    """Final positions from outputs/relax_<key>.out (Begin/End final coordinates, angstrom), species as written."""
    txt = (HERE / "outputs" / f"relax_{key}.out").read_text()
    blk = txt[txt.rindex("Begin final coordinates"):txt.rindex("End final coordinates")]
    rows = re.findall(r"^\s*([A-Za-z]+\d*)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", blk.split("ATOMIC_POSITIONS")[1], re.M)
    st = structure(key)
    pos = np.array([[float(v) for v in r[1:]] for r in rows])
    assert len(pos) == len(st.species), "relaxed atom count differs"
    assert [re.sub(r"\d", "", r[0]) for r in rows] == list(st.species), "species order differs"
    return S.Structure(st.cell, list(st.species), pos, st.pbc, st.frame)


def single_point(tag, st, spin_mode, k, nbnd=None, calc="scf"):
    if spin_mode == "nm":
        labels, mag = list(st.species), None
    else:
        labels, mag = spin_labels(st, spin_mode)
    body, order = cards(st, labels, k)
    return write(tag, control(calc, tag) + system(len(st.species), len(order), mag, nbnd) + ELECTRONS + body)


def stage2_ribbon(key):
    """Gate M single points at a relaxed geometry (stage 2; also the [110] H ribbon after its relaxation)."""
    h = ALL[key][2]
    tags = []
    st = relaxed(key)
    tags.append(single_point(f"nm_{key}", st, "nm", kauto(K15)))
    tags.append(single_point(f"fm_{key}", st, "fm", kauto(K15)))
    if not h:
        for mode in ("afm_edges", "afm_planes"):
            tags.append(single_point(f"{mode}_{key}", st, mode, kauto(K15)))
        for mode in ("nm", "fm", "afm_edges"):
            tags.append(single_point(f"{mode}_k30_{key}", st, mode, kauto(K30)))
    # bands: reuse the NM SCF charge density (same prefix/outdir), so write with the NM prefix. 28 bands per f.u.
    # (+ 1 per H) also survey where the W90 outer window (VBM + 11.05 eV HSE, cost test T2) falls: A4b's nbnd.
    nbnd = 28 * st.species.count("Sc") + st.species.count("H")
    # more vacuum than 15 Å (r110_4h at 20 Å): scale the empty bands by the cell volume, since the extra bands below the
    # survey level are vacuum states (added 2026-11-13; leaves the 15 Å ribbons' nbnd unchanged)
    vol = lambda c: c[1, 1] * c[2, 2]                                     # noqa: E731
    f_vac = vol(st.cell) / vol(structure(key, vacuum=VAC).cell)
    if f_vac > 1 + 1e-6:
        nocc = sum(dict(Sc=11, N=5, H=1)[x] for x in st.species) // 2
        nbnd = nocc + int(round((nbnd - nocc) * f_vac))
    txt = single_point(f"bands_{key}", st, "nm", kpath(), nbnd=nbnd, calc="bands")
    p = OUT / f"bands_{key}.in"
    p.write_text(p.read_text().replace(f"prefix='bands_{key}', outdir='./tmp_bands_{key}'",
                                       f"prefix='nm_{key}', outdir='./tmp_nm_{key}'"))
    tags.append(txt)
    tags.append(write(f"projwfc_{key}", f"&projwfc\n  prefix='nm_{key}', outdir='./tmp_nm_{key}',\n"
                                        f"  filpdos='pdos_{key}', degauss=0.005, DeltaE=0.01, kresolveddos=.true.\n/\n"))
    return tags


def recentred(key, vac):
    """The relaxed ribbon centred in the cell S.ribbon builds for vacuum `vac` (shift in y and z only)."""
    st = relaxed(key)
    d, n, h = ALL[key]
    big = S.ribbon(d, n, a_p=A_P, d=D_PLANE, vacuum=vac, passivate=h)
    shift = (np.diag(big.cell) - np.diag(st.cell)) / 2 * np.array([0, 1, 1])
    return S.Structure(big.cell, list(st.species), st.pos + shift, st.pbc, st.frame)


def vac20(key, vac=20.0):
    """NM SCF at k 15 with the relaxed ribbon centred in a cell with `vac` Å vacuum (Gate M criterion 5; 25 Å added
    2026-11-13)."""
    return single_point(f"vac{vac:.0f}_{key}", recentred(key, vac), "nm", kauto(K15))


def stage2():
    tags = []
    for key in RIBBONS:
        tags += stage2_ribbon(key)
    st = relaxed("r110_4")
    big = S.ribbon("110", 4, a_p=A_P, d=D_PLANE, vacuum=20.0)
    shift = (np.diag(big.cell) - np.diag(st.cell)) / 2 * np.array([0, 1, 1])
    st20 = S.Structure(big.cell, list(st.species), st.pos + shift, st.pbc, st.frame)
    tags.append(single_point("vac20_r110_4", st20, "nm", kauto(K15)))
    (HERE / "tags_stage2.txt").write_text("\n".join(tags) + "\n")
    # one job per ribbon (bands/projwfc reuse that ribbon's tmp_nm_ dir, so a ribbon's tags stay in one job, in order)
    for key in RIBBONS:
        mine = [t for t in tags if t.endswith(f"_{key}")]
        (HERE / f"tags_stage2_{key}.txt").write_text("\n".join(mine) + "\n")
    print("wrote", len(tags), "stage-2 inputs;", ", ".join(f"tags_stage2_{k}.txt" for k in RIBBONS))


def stage_vac20():
    """Decision 2026-11-12 after criterion 5 failed: the same 20 Å check on the [100] ribbons (bare and H); a
    ribbon whose NM gap moves > 5 meV also goes to 20 Å in A4b/A5. [110] ribbons use 20 Å."""
    tags = [vac20(k) for k in ("r100_3", "r100_3h")]
    (HERE / "tags_vac20.txt").write_text("\n".join(tags) + "\n")
    print("wrote", ", ".join(tags))


def stage_vac25():
    """The 25 Å check (2026-11-13, with A4b-1): do the two ribbons that failed 15 → 20 Å pass 20 → 25 Å (same 5 meV)?"""
    tags = [vac20(k, 25.0) for k in ("r100_3h", "r110_4")]
    (HERE / "tags_vac25.txt").write_text("\n".join(tags) + "\n")
    print("wrote", ", ".join(tags))


def stage_h110(step):
    """The [110] n = 4 H ribbon (LR3, decision 2026-11-12). step 1: relax (nspin 1, as r100_3h);
    step 2: NM + FM check, bands + projwfc (stage2_ribbon)."""
    key = "r110_4h"
    if step == 1:
        tags = [write(f"relax_{key}", relax_input(key))]
        st = structure(key)
        print(f"{key}: {len(st.species)} atoms, W {S.width(st) / 10:.2f} nm, cell {np.diag(st.cell).round(3)}")
    else:
        tags = stage2_ribbon(key)
    (HERE / f"tags_h110_{step}.txt").write_text("\n".join(tags) + "\n")
    print("wrote", ", ".join(tags))


if __name__ == "__main__":
    if "stage2" in sys.argv:
        stage2()
    elif "vac20" in sys.argv:
        stage_vac20()
    elif "vac25" in sys.argv:
        stage_vac25()
    elif "h110" in sys.argv:
        stage_h110(2 if "2" in sys.argv else 1)
    else:
        stage1()
