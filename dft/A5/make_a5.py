#!/usr/bin/env python3
"""A5 inputs: the two pre-registered hold-out ribbons (README.md, written before any A5 input existed).

Stage A5-1 (`make_a5.py`): one PAW PBE relaxation per hold-out, nspin 1, from the ideal cut, as A4a otherwise
(make_a4a.py: SSSP precision PAW, 90/720 Ry, gaussian 0.005 Ry, k 15, BFGS, forc 1e-3, etot 1e-5).
  r110_8   [110] n = 8 bare, 20 Å vacuum
  r100_6   [100] n = 6 bare, 15 Å vacuum
Stage A5-1b (`make_a5.py stage2`, after the relaxations are copied back to outputs/): at the relaxed geometry, as
A4a stage 2 (make_a4a.stage2_ribbon) but only what README.md registers: NM, FM, AFM-across-the-edges single points at
k 15 (Gate M criterion 1); bands on the 31-k Γ–X̄ path with A4a's band survey count (28 per f.u., empty bands scaled
by the cell's yz area beyond 15 Å: also the band survey for A5-2's NSCF nbnd, as A4a's was for A4b-1); projwfc on it.
One tag list per ribbon (bands/projwfc reuse that ribbon's tmp_nm_ dir).
Stage A5-2 (`make_a5.py twin`, after A5-1b is copied back): the NC PBE twin as A4b-1 (README.md "A5-2"), with
make_a4b.ribbon_inputs unchanged, pointed at the A5 relaxed geometry (already in its A5 cell) and written to twin/
(inputs/, w90/<r>/, tags_<r>.txt), so no A5-1 file is touched; on the cluster it runs in scn_fet/A5/twin with A4b-1's
run_a4b.sh (A4B_DIR). Set B gets Sc s on the edge Sc of bare [110] (EDGE_S), as r110_4. Band counts from the A5-1b
PAW survey (parse_a5.py → a5_stage2.json), with A4b-1's rules: bands = N(VBM + 4 eV) + 8 %; NSCF = N(L + 1 eV) + 3 %,
both rounded up to a multiple of 5.
Stage A5-3 (`make_a5.py hse`): the HSE06 gap single points as A4b-2 (README.md "A5-3"): make_a4b.hse_inputs
unchanged, nbnd = nocc + 32, only its first three steps kept (hse_<r> SCF on the IBZ of the 16 / 12-k grid with q 8 / 6,
hpot_<r> vacuum level, hproj_<r> projections); no open_grid / Wannier90 (pruned after generation). Written to hse/;
run with hse/run_a5_hse.sh (run_a4b2.sh's three steps; 63 ranks per node, one pool per IBZ k).

    ~/.venvs/scnfet/bin/python -P dft/A5/make_a5.py [stage2 | twin | hse]          (from the project root)
"""
import importlib.util
import json
import re
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

HOLDOUTS = {"r110_8": ("110", 8, False), "r100_6": ("100", 6, False)}
VAC = {"r110_8": 20.0, "r100_6": 15.0}
OUT = HERE / "inputs"


def structure(key):
    d, n, h = HOLDOUTS[key]
    return S.ribbon(d, n, a_p=MK.A_P, d=MK.D_PLANE, vacuum=VAC[key], passivate=h)


def write(tag, text):
    OUT.mkdir(exist_ok=True)
    (OUT / f"{tag}.in").write_text(text)
    return tag


def relax_input(key):
    st = structure(key)
    body, order = MK.cards(st, list(st.species), MK.kauto(MK.K15))
    return (MK.control("relax", f"relax_{key}", ",\n  forc_conv_thr=1.0d-3, etot_conv_thr=1.0d-5, nstep=300")
            + MK.system(len(st.species), len(order)) + MK.ELECTRONS + "&ions\n  ion_dynamics='bfgs'\n/\n" + body)


def stage1():
    info = {}
    for key in HOLDOUTS:
        tag = write(f"relax_{key}", relax_input(key))
        st = structure(key)
        info[key] = dict(vacuum=VAC[key], nat=len(st.species), nSc=st.species.count("Sc"),
                         W_nm=round(S.width(st) / 10, 3), cell=[round(float(x), 4) for x in st.cell.diagonal()])
        print(f"{tag}: {len(st.species)} atoms, W {S.width(st) / 10:.3f} nm, cell {st.cell.diagonal().round(3)}")
    (HERE / "a5_inputs.json").write_text(json.dumps(info, indent=1))


def relaxed(key):
    """Final positions from outputs/relax_<key>.out (as make_a4a.relaxed, with the A5 structures)."""
    txt = (HERE / "outputs" / f"relax_{key}.out").read_text()
    assert "End of BFGS Geometry Optimization" in txt and "bfgs converged" in txt, f"relax_{key} not converged"
    blk = txt[txt.rindex("Begin final coordinates"):txt.rindex("End final coordinates")]
    rows = re.findall(r"^\s*([A-Za-z]+\d*)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", blk.split("ATOMIC_POSITIONS")[1], re.M)
    st = structure(key)
    pos = np.array([[float(v) for v in r[1:]] for r in rows])
    assert len(pos) == len(st.species), "relaxed atom count differs"
    assert [re.sub(r"\d", "", r[0]) for r in rows] == list(st.species), "species order differs"
    return S.Structure(st.cell, list(st.species), pos, st.pbc, st.frame)


def geometry(key):
    """Descriptive: displacement from the ideal cut (edge row vs interior; edge row = |y - y_edge| < 0.1 Å, as Gate M),
    the relaxed width, and the energy gained."""
    st0, st = structure(key), relaxed(key)
    d = st.pos - st0.pos
    d[:, 0] -= st.cell[0, 0] * np.round(d[:, 0] / st.cell[0, 0])
    y0 = st0.pos[:, 1]
    edge = (np.abs(y0 - y0.min()) < 0.1) | (np.abs(y0 - y0.max()) < 0.1)
    txt = (HERE / "outputs" / f"relax_{key}.out").read_text()
    e = [float(x) for x in re.findall(r"^!\s+total energy\s+=\s+(-?[\d.]+) Ry", txt, re.M)]
    nrm = np.linalg.norm(d, axis=1)
    return dict(max_disp_edge=round(float(nrm[edge].max()), 4), max_disp_interior=round(float(nrm[~edge].max()), 4),
                mean_disp_interior=round(float(nrm[~edge].mean()), 4), n_edge_atoms=int(edge.sum()),
                W_ideal_nm=round(S.width(st0) / 10, 4), W_relaxed_nm=round(S.width(st) / 10, 4),
                dE_relax_eV=round((e[-1] - e[0]) * 13.605693123, 4), n_scf=len(e),
                max_dz_edge=round(float(np.abs(d[edge, 2]).max()), 4), max_dy_edge=round(float(np.abs(d[edge, 1]).max()), 4))


def stage2_ribbon(key):
    MK.OUT = OUT                                   # make_a4a.write writes to its module OUT: point it here
    st = relaxed(key)
    k15 = MK.kauto(MK.K15)
    tags = [MK.single_point(f"{m}_{key}", st, m, k15) for m in ("nm", "fm", "afm_edges")]
    # bands: the A4a survey count (make_a4a.stage2_ribbon): 28 per f.u., empty bands scaled by the yz area beyond 15 Å
    nbnd = 28 * st.species.count("Sc")
    f_vac = (st.cell[1, 1] * st.cell[2, 2]) / np.prod(np.diag(S.ribbon(*HOLDOUTS[key][:2], a_p=MK.A_P, d=MK.D_PLANE,
                                                                       vacuum=15.0).cell)[1:])
    nocc = sum(dict(Sc=11, N=5, H=1)[x] for x in st.species) // 2
    if f_vac > 1 + 1e-6:
        nbnd = nocc + int(round((nbnd - nocc) * f_vac))
    MK.single_point(f"bands_{key}", st, "nm", MK.kpath(), nbnd=nbnd, calc="bands")
    p = OUT / f"bands_{key}.in"
    p.write_text(p.read_text().replace(f"prefix='bands_{key}', outdir='./tmp_bands_{key}'",
                                       f"prefix='nm_{key}', outdir='./tmp_nm_{key}'"))
    tags.append(f"bands_{key}")
    tags.append(MK.write(f"projwfc_{key}", f"&projwfc\n  prefix='nm_{key}', outdir='./tmp_nm_{key}',\n"
                                           f"  filpdos='pdos_{key}', degauss=0.005, DeltaE=0.01, kresolveddos=.true.\n/\n"))
    (HERE / f"tags_stage2_{key}.txt").write_text("\n".join(tags) + "\n")
    return tags, dict(nbnd_bands=nbnd, nocc=nocc, f_vac=round(float(f_vac), 4), **geometry(key))


def stage2():
    info = json.loads((HERE / "a5_inputs.json").read_text())
    for key in HOLDOUTS:
        tags, g = stage2_ribbon(key)
        info[key]["stage2"] = g
        print(f"{key}: {', '.join(tags)}")
        print("   ", g)
    (HERE / "a5_inputs.json").write_text(json.dumps(info, indent=1))


TWIN = HERE / "twin"


def _ceil5(x):
    return int(-(-x // 5) * 5)


def twin_counts():
    """A5-2 band counts from the A5-1b PAW band survey (a5_stage2.json and outputs/bands_<r>.out)."""
    import importlib.util as iu
    spec = iu.spec_from_file_location("parse_a5", HERE / "parse_a5.py")
    p5 = iu.module_from_spec(spec)
    spec.loader.exec_module(p5)
    res = json.loads((HERE / "a5_stage2.json").read_text())
    nb = {}
    for key in HOLDOUTS:
        nocc = res[key]["nocc"]
        kp, ep, _ = p5.PA.P.bands(p5.PA.out(f"bands_{key}"))
        emin = ep.min(axis=0) - ep[:, nocc - 1].max()
        n4 = int(np.where(emin >= 4.0)[0][0]) + 1
        nl1 = res[key]["survey"]["levels"]["L+1.0"]["nbnd"]
        nb[key] = dict(N_vbm4=n4, N_L1=nl1, bands=_ceil5(n4 * 1.08), nscf=_ceil5(nl1 * 1.03))
    return nb


def _make_a4b(where):
    """make_a4b, unchanged, pointed at the A5 hold-outs: relaxed geometry (already in the A5 cells, 20 / 15 Å), the A5
    ideal cuts for the edge rows, and output under `where` (inputs/, w90/, tags_*.txt)."""
    import types
    spec = importlib.util.spec_from_file_location("make_a4b", ROOT / "dft" / "A4b" / "make_a4b.py")
    MB = importlib.util.module_from_spec(spec)
    argv, sys.argv[:] = list(sys.argv), ["make_a4b.py"]
    spec.loader.exec_module(MB)
    sys.argv[:] = argv
    MB.MK = types.SimpleNamespace(structure=structure, ALL=HOLDOUTS, kauto=MK.kauto, K15=MK.K15, kpath=MK.kpath)
    MB.geometry = relaxed
    MB.HERE, MB.OUT, MB.W90 = where, where / "inputs", where / "w90"
    return MB


def twin():
    nb = twin_counts()
    MB = _make_a4b(TWIN)
    MB.EDGE_S = set(MB.EDGE_S) | {"r110_8"}          # bare [110]: Sc s on the edge Sc (set B), as r110_4
    for key in HOLDOUTS:
        MB.NB_BANDS[key], MB.NB_NSCF[key] = nb[key]["bands"], nb[key]["nscf"]
    TWIN.mkdir(exist_ok=True)
    info = json.loads((HERE / "a5_inputs.json").read_text())
    for key in HOLDOUTS:
        st, maps = MB.ribbon_inputs(key)
        ne, no = MB.nexcl(st), MB.nocc(st)
        info[key]["twin"] = dict(nb[key], nexcl=ne, nocc_nc=no, nk_w90=MB.NK_W90[HOLDOUTS[key][0]],
                                 num_wann_B=len(maps["B"]), num_wann_S=len(maps["S"]), num_bands=nb[key]["nscf"] - ne,
                                 edge_s=sum(1 for m in maps["B"] if m["orbital"] == "s"))
        print(f"{key}: bands {nb[key]['bands']} (N(VBM+4) {nb[key]['N_vbm4']}), NSCF {nb[key]['nscf']} (N(L+1) "
              f"{nb[key]['N_L1']}), k {MB.NK_W90[HOLDOUTS[key][0]]}, num_wann B {len(maps['B'])} / S {len(maps['S'])}, "
              f"num_bands {nb[key]['nscf'] - ne}, excluded {ne}, occupied {no}")
    (HERE / "a5_inputs.json").write_text(json.dumps(info, indent=1))


HSE = HERE / "hse"


def hse():
    import shutil
    MB = _make_a4b(HSE)
    HSE.mkdir(exist_ok=True)
    info = json.loads((HERE / "a5_inputs.json").read_text())
    for key in HOLDOUTS:
        st = relaxed(key)
        MB.NB_HSE[key] = MB.nocc(st) + 32
        MB.hse_inputs(key)
        tags = [f"hse_{key}", f"hpot_{key}", f"hproj_{key}"]
        (HSE / f"tags_hse_{key}.txt").write_text("\n".join(tags) + "\n")
        (HSE / "inputs" / f"hopen_{key}.in").unlink()
        d_ = HOLDOUTS[key][0]
        info[key]["hse"] = dict(nbnd=MB.NB_HSE[key], nocc=MB.nocc(st), nk=MB.NK_W90[d_], nk_ibz=MB.NK_IBZ[d_], nq=MB.NQ[d_])
        print(f"{key}: HSE nbnd {MB.NB_HSE[key]} (nocc {MB.nocc(st)} + 32), k {MB.NK_W90[d_]} (IBZ {MB.NK_IBZ[d_]}), "
              f"q {MB.NQ[d_]}; tags {', '.join(tags)}")
    shutil.rmtree(HSE / "w90")                       # hse_inputs also writes the A4b-2 Wannier inputs: not used here
    (HERE / "a5_inputs.json").write_text(json.dumps(info, indent=1))


if __name__ == "__main__":
    if "stage2" in sys.argv:
        stage2()
    elif "twin" in sys.argv:
        twin()
    elif "hse" in sys.argv:
        hse()
    else:
        stage1()
