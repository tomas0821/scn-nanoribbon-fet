"""LR1 test 1 (decision 2026-11-11): add the Sc 4s orbital to the bulk SK model (plan §10 LR1 mitigation).

Model: the plan's 13 parameters + Es, sp_s (NN Sc s–N p σ), ss1_s and sd1_s (2NN), as an "Ss" pseudo-site at the Sc
position (sk.with_s, sk.s_shells). Targets: as scripts/a1_sk_fit.py, but the DFT grid now covers all 9 bands
(QE 6–14) up to the outer window (13.5 / 11.8 eV), below the first band outside the manifold. Variants, HSE06:
  S   balanced weighting (far σ 0.1 eV), as H
  Sn  near-gap only (far σ 0.5 eV, masses ×3), as Hn
  PS  the PBE twin, balanced
Seeds: H / Hn / P fits with Es = 8, sp_s = 1.5, ss1_s = −0.3, sd1_s = −0.3, plus multi-start.
Scored on gate (a) exactly as the 8-band fits.

Run from the project root: PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/a1_sk_fit_s.py
Writes fits/A1_bulk/sk_fit_s.json.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from scnfet import sk

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("a1_sk_fit", ROOT / "scripts" / "a1_sk_fit.py")
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)
A, OUT = F.A, F.OUT
S0 = dict(Es=8.0, sp_s=1.5, ss1_s=-0.3, sd1_s=-0.3)


def main():
    only = sys.argv[1:] or None
    st = sk.with_s(sk.rocksalt(A))
    model = sk.build(st, shells=sk.rocksalt_shells(A) + sk.s_shells(A), names=sk.PARAM_NAMES + sk.S_PARAM_NAMES)
    base = json.loads((OUT / "sk_fits.json").read_text())["results"]
    ext = json.loads((OUT / "sk_fit_extensions.json").read_text())
    res = {}
    for name, lab, near_only, seeds in (("S", "HSE06", False, ("H", "Hn")), ("Sn", "HSE06", True, ("Hn", "H")),
                                        ("PS", "PBE", False, ("P",))):
        if only and name not in only:
            continue
        case = F.Case(lab, model)
        T = case.targets()
        if near_only:
            for mt in T.masses:
                mt.rel_sigma /= 3
            T.sigma = np.where(~T.floor & (T.sigma > F.SIG_NEAR), 0.5, T.sigma)
        allp = {**base, **ext}
        x0s = {s: model.x_from({**allp[s]["params"], **S0}) for s in seeds}
        rec, _ = F.run(model, case, T, ("Ep",), x0s, name)
        rec["n_targets"] = dict(n=int(len(T.E)), n_k=int(len(T.k)), floors=int(T.floor.sum()))
        res[name] = rec
    f = OUT / "sk_fit_s.json"
    old = json.loads(f.read_text()) if f.exists() else {}
    f.write_text(json.dumps({**old, **res}, indent=1, default=float))
    print("wrote", f)


if __name__ == "__main__":
    main()
