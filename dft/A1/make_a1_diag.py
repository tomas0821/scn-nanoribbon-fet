#!/usr/bin/env python3
"""A1 Wannier diagnostic (criterion 4 failed): which change fixes the PBE-twin interpolation?

Cases, all from the existing PBE NSCF (tmp_pbe_scf) and the same windows as `pbe`:
  pbe_gi   8 WFs, use_ws_distance = .false.  (the original setup; consistency check against parse_a1_w90.py)
  pbe_ws   8 WFs, use_ws_distance = .true.   (cause i)
  pbe9     9 WFs, + Sc:s, use_ws_distance = .false.  (cause ii)
  pbe9_ws  9 WFs, + Sc:s, use_ws_distance = .true.   (both)
Each run uses Wannier90's geninterp to evaluate its bands at exactly the 201 k-points of pbe_scf_bands, so the
comparison with Quantum ESPRESSO needs no wsvec reader on our side.
"""
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
A_RED = np.array([[-0.5, 0, 0.5], [0, 0.5, 0.5], [-0.5, 0.5, 0]])        # QE ibrav = 2 rows a_i / a


def qe_path_k():
    txt = (HERE / "outputs" / "pbe_scf_bands.out").read_text()
    i = txt.rfind("End of band structure calculation")
    ks = re.findall(r"k =\s*(-?[\d.]+)\s*(-?[\d.]+)\s*(-?[\d.]+) \(\s*\d+ PWs\)", txt[i:])
    return np.array(ks, float)                                            # Cartesian, 2π/a


def main():
    base = (HERE / "w90" / "pbe" / "pbe.win").read_text()
    kc = qe_path_k()
    frac = kc @ A_RED.T                                                   # f_i = a_i · k (2π/a units)
    gk = "path k-points of pbe_scf_bands (crystal)\ncrystal\n" + f"{len(frac)}\n" + "".join(
        f"{i + 1} {a:.10f} {b:.10f} {c:.10f}\n" for i, (a, b, c) in enumerate(frac))
    cases = {"pbe_gi": (8, False), "pbe_ws": (8, True), "pbe9": (9, False), "pbe9_ws": (9, True)}
    for seed, (nw, ws) in cases.items():
        w = base
        w = re.sub(r"num_wann\s*=\s*8", f"num_wann  = {nw}", w)
        if nw == 9:
            w = w.replace("Sc:d", "Sc:s;d")
            # band 14 (Sc 4s) reaches VBM + 13.11 eV: the outer window must hold ≥ 9 bands at every k
            import json
            vbm = json.loads((HERE / "a1_stage1.json").read_text())["survey"]["pbe_nscf"]["vbm"]
            w = re.sub(r"dis_win_max = [\d.]+", f"dis_win_max = {vbm + 14.0:.4f}", w)
        w = re.sub(r"use_ws_distance = \.false\.[^\n]*", f"use_ws_distance = {'.true.' if ws else '.false.'}", w)
        w = w.replace("bands_plot = .true.", "bands_plot = .true.\ngeninterp = .true.")
        d = HERE / "w90" / seed
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{seed}.win").write_text(w)
        (d / f"{seed}_geninterp.kpt").write_text(gk)
        if seed == "pbe9":
            (d / "pbe9.pw2wan").write_text("&inputpp\n  outdir='../../tmp_pbe_scf', prefix='pbe_scf', seedname='pbe9',\n"
                                           "  write_amn=.true., write_mmn=.true., write_unk=.false.\n/\n")
    print("wrote", list(cases), "| geninterp points:", len(frac))


if __name__ == "__main__":
    main()
