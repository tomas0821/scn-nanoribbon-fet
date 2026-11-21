"""LR1 test 2 (decision 2026-11-11): a three-centre Sc–N–Sc term on the 2NN Sc–Sc d–d block (sk.ThreeCentre).

Why: the A1 Wannier 2NN Sc–Sc block is not two-centre (its x²−y² element is −0.23 eV where two-centre theory
predicts +0.39). With the three-centre term the projection reproduces that block exactly (max |Δ| 0.0 meV) and the
in-model residual falls from 105 to 7 meV (HSE06). Not in the plan: adopting it would be a model change (a
separate decision).

Model: the plan's 13 parameters + t3_s, t3_ppp, t3_psp. t3_pss is fixed at 0: σ–σ through N needs ĉ_ik·ĉ_kj ≠ 0, and
Sc–N–Sc is a 90° path in rocksalt, so the term vanishes identically in bulk (nearly so in the sheet).
Targets, weights and the gate (a) score: exactly as scripts/a1_sk_fit.py (8 bands). Variants:
  T   HSE06, balanced (as H)        Tn  HSE06, near-gap only (as Hn)        PT  PBE twin, balanced (as P)
  Tm  as Tn with mass weights ×10 instead of ×3: can the masses reach ±10 % while the ±1 eV RMS stays ≤ 30 meV?
  TSm5 / TSm10  as TSn (three-centre + Sc 4s, near-gap) with mass weights ×5 / ×10, seeded from TSn
  [TS / TSn: three-centre + Sc 4s, run only if asked on the command line]
Seeds: the matching 13-parameter fit + projected t3 values, and the projection of the extended model.

Run from the project root: PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/a1_sk_fit_t3.py [T Tn PT TS TSn]
Writes fits/A1_bulk/sk_fit_t3.json.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from scnfet import sk, skfit

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("a1_sk_fit", ROOT / "scripts" / "a1_sk_fit.py")
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)
A, OUT = F.A, F.OUT
FIXED = ("Ep", "t3_pss")


def make_model(with_s=False):
    tc = sk.three_centre_rocksalt(A)
    st = sk.with_s(sk.rocksalt(A)) if with_s else sk.rocksalt(A)
    shells = sk.rocksalt_shells(A) + (sk.s_shells(A) if with_s else ())
    names = sk.PARAM_NAMES + (sk.S_PARAM_NAMES if with_s else ()) + sk.T3_PARAM_NAMES
    return sk.build(st, shells=shells, names=names, three_centre=tc)


def projected(model, case):
    if model.norb != 8:
        return None
    x, _ = skfit.project_hr(model, case.hm)
    p = dict(zip(model.names, x))
    ep = p["Ep"]
    for k in ("Ep", "Et2g", "Eeg"):
        p[k] -= ep
    p["t3_pss"] = 0.0
    return p


def main():
    only = sys.argv[1:] or ["T", "Tn", "PT"]
    prev = {**json.loads((OUT / "sk_fits.json").read_text())["results"],
            **json.loads((OUT / "sk_fit_extensions.json").read_text())}
    s_fits = json.loads((OUT / "sk_fit_s.json").read_text()) if (OUT / "sk_fit_s.json").exists() else {}
    res = {}
    for name, lab, near_only, with_s, base in (("T", "HSE06", False, False, "H"), ("Tn", "HSE06", True, False, "Hn"),
                                               ("PT", "PBE", False, False, "P"), ("TS", "HSE06", False, True, "S"),
                                               ("TSn", "HSE06", True, True, "Sn"), ("Tm", "HSE06", True, False, "Hn"),
                                               ("TSm5", "HSE06", True, True, "Sn"), ("TSm10", "HSE06", True, True, "Sn")):
        if name not in only:
            continue
        model = make_model(with_s)
        case = F.Case(lab, model)
        T = case.targets()
        if near_only:
            for mt in T.masses:
                mt.rel_sigma /= {"Tm": 10, "TSm5": 5, "TSm10": 10}.get(name, 3)
            T.sigma = np.where(~T.floor & (T.sigma > F.SIG_NEAR), 0.5, T.sigma)
        pj = projected(model, F.Case(lab, make_model(False))) if not with_s else None
        t3 = {k: v for k, v in (pj or {}).items() if k.startswith("t3")} or dict(t3_s=0.45, t3_ppp=0.1, t3_psp=-0.15)
        bp = (s_fits if with_s else prev)[base]["params"]
        if name in ("Tm", "TSm5", "TSm10"):
            bp = json.loads((OUT / "sk_fit_t3.json").read_text())["Tn" if name == "Tm" else "TSn"]["params"]
        x0s = ({f"{base}+t3": model.x_from({**bp, **t3, "t3_pss": 0.0})} if name not in ("Tm", "TSm5", "TSm10")
               else {"Tn" if name == "Tm" else "TSn": model.x_from(bp)})
        if pj is not None:
            x0s["proj"] = model.x_from(pj)
        rec, _ = F.run(model, case, T, FIXED, x0s, name)
        rec["n_targets"] = dict(n=int(len(T.E)), n_k=int(len(T.k)), floors=int(T.floor.sum()))
        res[name] = rec
    f = OUT / "sk_fit_t3.json"
    old = json.loads(f.read_text()) if f.exists() else {}
    f.write_text(json.dumps({**old, **res}, indent=1, default=float))
    print("wrote", f)


if __name__ == "__main__":
    main()
