#!/usr/bin/env python3
"""A1 Wannier diagnostic, round 2 (criterion 4 still failed after round 1). PBE only, same windows as `pbe`:
  pbe_mlwf  8³ grid, the same 8 WFs maximally localised (num_iter = 200): reuses pbe.amn/.mmn/.eig
  pbe_k12   12³ grid (new NSCF, 1728 k), projection-only
  pbe_k16   16³ grid (new NSCF, 4096 k), projection-only
Each NSCF runs in its own copy of tmp_pbe_scf (prefix pbe_scf), so the 8³ wavefunctions stay untouched.
Scored by parse_a1_diag.py (criterion 4 on the exact QE path).
"""
import re
import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("make_a1", HERE / "make_a1.py")
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)


def main():
    base_nscf = (HERE / "inputs" / "pbe_nscf.in").read_text()
    head = base_nscf[:base_nscf.index("K_POINTS")]
    pbe_win = (HERE / "w90" / "pbe" / "pbe.win").read_text()
    froz = float(re.search(r"dis_froz_max = ([\d.]+)", pbe_win).group(1))
    wmax = float(re.search(r"dis_win_max = ([\d.]+)", pbe_win).group(1))
    gk = (HERE / "w90" / "pbe_gi" / "pbe_gi_geninterp.kpt").read_text()       # the QE path, for scoring
    for n in (12, 16):
        txt = head.replace("outdir='./tmp_pbe_scf'", f"outdir='./tmp_pbe_k{n}'") + M.kblock_list(n)
        (HERE / "inputs" / f"pbe_nscf_k{n}.in").write_text(txt)
        seed = f"pbe_k{n}"
        d = HERE / "w90" / seed
        d.mkdir(parents=True, exist_ok=True)
        w = M.win(seed, M.full_grid(n), froz, wmax).replace("mp_grid = 8 8 8", f"mp_grid = {n} {n} {n}")
        (d / f"{seed}.win").write_text(w)
        (d / f"{seed}.pw2wan").write_text(f"&inputpp\n  outdir='../../tmp_pbe_k{n}', prefix='pbe_scf', seedname='{seed}',\n"
                                          "  write_amn=.true., write_mmn=.true., write_unk=.false.\n/\n")
        (d / f"{seed}_geninterp.kpt").write_text(gk)
    d = HERE / "w90" / "pbe_mlwf"
    d.mkdir(parents=True, exist_ok=True)
    w = pbe_win.replace("num_iter = 0     ! projection-only: SK-compatible character",
                        "num_iter = 200   ! maximally localised (round-2 diagnostic)\nconv_tol = 1.0e-10\nconv_window = 5")
    assert "num_iter = 200" in w
    (d / "pbe_mlwf.win").write_text(w)
    (d / "pbe_mlwf_geninterp.kpt").write_text(gk)
    print("wrote pbe_nscf_k12/k16 inputs and w90/pbe_k12, pbe_k16, pbe_mlwf | windows", froz, wmax)


if __name__ == "__main__":
    main()
