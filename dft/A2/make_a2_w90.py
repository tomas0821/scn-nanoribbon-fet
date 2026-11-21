#!/usr/bin/env python3
"""A2-2: Wannier90 for the sheet (PBE twin `pbe_k16`, HSE06 `hse_k16`), 16 projection-only WFs on 16×16×1.

Windows by the pre-registered rule (README criterion 3), from a2_stage1.json:
  frozen top = min(E_vac − 0.5 eV, first non-model band minimum − 0.1 eV)   [eV, absolute]
  outer top  = VBM + 13.5 (HSE) / + 11.8 (PBE), as bulk A1, capped below the lowest top-band energy (nbnd 56)
E_vac for PBE is the planar-averaged bare + Hartree potential in the vacuum from `pbe_scf` (pp.x). For HSE it is
estimated from the same PBE potential (same density to first order) and **checked** in this job by `hse_pot` (pp.x
on the HSE density); parse_a2.py reports the margin, and the HSE Wannier run is redone if it is under 0.5 eV.
"""
import json
from pathlib import Path

import importlib.util

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("make_a2", HERE / "make_a2.py")
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)
OUTER = {"pbe_k16": 11.8, "hse_k16": 13.5}


def main():
    s = json.loads((HERE / "a2_stage1.json").read_text())
    e_vac = s["vacuum"]["E_vac_eV"]
    rows = {}
    for seed, run in (("pbe_k16", "pbe_nscf_k16"), ("hse_k16", "hse_ibz_k16")):
        vbm = s[run]["vbm"]
        sv = s[f"survey_{run}"]
        froz = min(e_vac - 0.5, vbm + sv["first_nonmodel_min"] - 0.1)
        outer = min(vbm + OUTER[seed], vbm + sv["top_band_min"] - 0.05)
        rows[seed] = dict(vbm=vbm, froz=froz, froz_rel=froz - vbm, outer=outer, outer_rel=outer - vbm,
                          vac_margin=e_vac - froz)
        d = HERE / "w90" / seed
        d.mkdir(parents=True, exist_ok=True)
        if seed == "pbe_k16":
            text = M.win(seed, M.full_grid(), froz, outer)
            od, pre = "tmp_pbe_scf", "pbe_scf"
        else:
            text = M.win(seed, [], froz, outer).replace("begin kpoints\nend kpoints", "begin kpoints\n@KLIST@\nend kpoints")
            od, pre = "tmp_hse_ibz_k16", "hse_ibz_k16_open"
        (d / f"{seed}.win").write_text(text)
        (d / f"{seed}.pw2wan").write_text(f"&inputpp\n  outdir='../../{od}', prefix='{pre}', seedname='{seed}',\n"
                                         "  write_amn=.true., write_mmn=.true., write_unk=.false.\n/\n")
    (M.OUT / "hse_pot.in").write_text(
        "&inputpp\n  prefix='hse_ibz_k16', outdir='./tmp_hse_ibz_k16', filplot='hse_ibz_k16.pot', plot_num=11\n/\n"
        "&plot\n  nfile=1, filepp(1)='hse_ibz_k16.pot', weight(1)=1.0, iflag=3, output_format=6, fileout='hse_pot.cube'\n/\n")
    (HERE / "a2_windows.json").write_text(json.dumps(rows, indent=1))
    print(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
