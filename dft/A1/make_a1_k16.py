#!/usr/bin/env python3
"""A1 production HSE06 Wannier on a 16³ grid (decision 2026-11-11, after the round-2 diagnostic showed the 8³
grid limits the interpolation). Route O: HSE06 SCF on the 16³ IBZ (q 8³ divides 16), open_grid.x, Wannier90 with
the same windows and projection-only WFs as hse_q8. The PBE twin on 16³ is w90/pbe_k16 (already done)."""
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("make_a1", HERE / "make_a1.py")
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)

v = json.loads((HERE / "a1_stage1.json").read_text())["survey"]["hse_nosym"]["vbm"]
(M.OUT / "hse_ibz_k16.in").write_text(
    M.control("scf", "hse_ibz_k16") + M.system(M.FIXED, M.NBND, M.hse(8)) + M.electrons("1.0d-9", ", diago_full_acc=.true.")
    + M.SPECIES + "K_POINTS automatic\n  16 16 16 0 0 0\n")
(M.OUT / "hse_ibz_k16_open.in").write_text("&inputpp\n  outdir='./tmp_hse_ibz_k16', prefix='hse_ibz_k16'\n/\n")
d = HERE / "w90" / "hse_k16"
d.mkdir(parents=True, exist_ok=True)
w = M.win("hse_k16", [], v + 6.5, v + 13.5).replace("mp_grid = 8 8 8", "mp_grid = 16 16 16")
w = w.replace("begin kpoints\nend kpoints", "begin kpoints\n@KLIST@\nend kpoints")
(d / "hse_k16.win").write_text(w)
(d / "hse_k16.pw2wan").write_text("&inputpp\n  outdir='../../tmp_hse_ibz_k16', prefix='hse_ibz_k16_open', seedname='hse_k16',\n"
                                  "  write_amn=.true., write_mmn=.true., write_unk=.false.\n/\n")
print(f"windows: frozen {v + 6.5:.4f}, outer {v + 13.5:.4f}")
