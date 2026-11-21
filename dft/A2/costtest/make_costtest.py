#!/usr/bin/env python3
"""HSE cost tests on the A2 sheet: T1 (smearing vs fixed EXX cost) and T2 (lower Wannier outer window).
Criteria and decisions: README.md (written before any output). Writes inputs/, w90/<seed>_nb<N>/ and windows.json."""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
A2 = HERE.parent
WIN = json.loads((A2 / "a2_windows.json").read_text())
NEXCL = 10
NWANN = 16
OUTER_REL = {"hse_k16": 13.5, "pbe_k16": 11.8}       # the A2 rule (bulk A1 values)
NBND = (34, 36, 40, 44)
T1 = {"ct_fix24": ("fixed", 24), "ct_fix40": ("fixed", 40), "ct_sm24": ("smearing", 24), "ct_sm40": ("smearing", 40)}
SMEAR = "occupations='smearing', smearing='gaussian', degauss=0.005"


def t1():
    base = (A2 / "inputs" / "hse_k8_q4.in").read_text()
    (HERE / "inputs").mkdir(exist_ok=True)
    for tag, (occ, nb) in T1.items():
        s = base.replace("prefix='hse_k8_q4', outdir='./tmp_hse_k8_q4'", f"prefix='{tag}', outdir='./tmp_{tag}'")
        s = s.replace("occupations='fixed', nbnd=24", ("occupations='fixed'" if occ == "fixed" else SMEAR) + f", nbnd={nb}")
        assert f"prefix='{tag}'" in s and f"nbnd={nb}" in s, tag
        (HERE / "inputs" / f"{tag}.in").write_text(s)
    (HERE / "tags_t1.txt").write_text("\n".join(T1) + "\n")


def eig(seed):
    d = np.loadtxt(HERE / "ref" / f"{seed}.eig")        # copied from the cluster's A2/w90/<seed>/
    nb, nk = int(d[:, 0].max()), int(d[:, 1].max())
    return d[:, 2].reshape(nk, nb)                      # W90 band b is total band b + NEXCL


def t2():
    rows = {}
    for seed in ("hse_k16", "pbe_k16"):
        E = eig(seed)
        vbm, froz, base = WIN[seed]["vbm"], WIN[seed]["froz"], WIN[seed]["outer"]
        win = (A2 / "w90" / seed / f"{seed}.win").read_text()
        assert f"dis_win_max = {base:.4f}" in win
        for N in NBND:
            outer = min(vbm + OUTER_REL[seed], E[:, N - NEXCL - 1].min() - 0.05)
            if abs(outer - base) < 1e-6:
                continue                                 # identical to the A2 baseline
            n_in = (E <= outer).sum(1)
            assert n_in.min() >= NWANN and outer > froz + 1.0, (seed, N, n_in.min(), outer)
            v = f"{seed}_nb{N}"
            d = HERE / "w90" / v
            d.mkdir(parents=True, exist_ok=True)
            (d / f"{seed}.win").write_text(win.replace(f"dis_win_max = {base:.4f}", f"dis_win_max = {outer:.4f}"))
            rows[v] = dict(seed=seed, nbnd=N, per_fu=N / 2, outer=outer, outer_rel_vbm=outer - vbm,
                           bands_in_outer_min=int(n_in.min()), bands_in_outer_max=int(n_in.max()))
    (HERE / "windows.json").write_text(json.dumps(rows, indent=1))
    for v, r in rows.items():
        print(f"{v:14s} outer {r['outer']:7.3f} (VBM + {r['outer_rel_vbm']:5.2f})  bands/k {r['bands_in_outer_min']}-"
              f"{r['bands_in_outer_max']}")


if __name__ == "__main__":
    t1()
    t2()
