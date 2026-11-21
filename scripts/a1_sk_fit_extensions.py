"""Diagnostic for the SK basis decision: which extra shells would let the
bulk HSE06 fit pass gate (a)? Evidence only. The model specification (plan D3) was not changed by this diagnostic.

Base = the 13-parameter model of a1_sk_fit.py (variant H), same targets and weights. Extensions, each a two-centre
shell with free signs:
  pp2  N–N at a (6 bonds; the largest dropped Wannier hops, up to 124 meV)        +2
  pd2  Sc–N at a√3/2 = 3.90 Å (8 bonds; dropped up to 53 meV)                       +2
  dd3  Sc–Sc at a√(3/2) = 5.51 Å (24 bonds; shifts the d_xy Γ–X spacing by ≈ 0.5 eV) +3
  all  the three together                                                          +7
  Hm   base model, mass weights ×3 (rel σ 0.01 / 0.017): can the base reach the masses at an energy cost?
  Hn   base model, near-gap only: far targets σ 0.5 eV (from 0.1), mass weights ×3: gate (a) scores only ±1 eV
       around the edges and the masses, so this is the most gate-friendly weighting of the plan's model
  Hn_all  as Hn with all three extra shells
Seeds: the H fit (extra parameters 0) and the Wannier projection onto the extended model.

Run from the project root: PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/a1_sk_fit_extensions.py
Writes fits/A1_bulk/sk_fit_extensions.json.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np

from scnfet import sk, skfit

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("a1_sk_fit", ROOT / "scripts" / "a1_sk_fit.py")
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)
A = F.A
OUT = F.OUT


def shells_ext(which):
    base = sk.rocksalt_shells(A)
    nn, n2, n3, n4, n5 = A / 2, A / np.sqrt(2), A * np.sqrt(3) / 2, A, A * np.sqrt(1.5)
    extra = {"pd2": sk.Shell("pd2", ("N", "Sc"), n3, (nn + n3) / 2, (n3 + A * np.sqrt(5) / 2) / 2, 3.5, ("pd2_s", "pd2_p")),
             "pp2": sk.Shell("pp2", ("N", "N"), n4, (n2 + n4) / 2, (n4 + n5) / 2, 2.0, ("pp2_s", "pp2_p")),
             "dd3": sk.Shell("dd3", ("Sc", "Sc"), n5, (n4 + n5) / 2, (n5 + A * np.sqrt(2)) / 2, 5.0,
                             ("dd3_s", "dd3_p", "dd3_d"))}
    sel = [extra[w] for w in which]
    return base + tuple(sel), sk.PARAM_NAMES + tuple(p for s in sel for p in s.params)


def main():
    import sys
    only = sys.argv[1:] or None
    prev = json.loads((OUT / "sk_fits.json").read_text())["results"]["H"]["params"]
    res = {}
    for name, which, mass_scale in (("pp2", ["pp2"], 1), ("pd2", ["pd2"], 1), ("dd3", ["dd3"], 1),
                                    ("all", ["pp2", "pd2", "dd3"], 1), ("Hm", [], 3), ("Hn", [], 3),
                                    ("Hn_all", ["pp2", "pd2", "dd3"], 3)):
        if only and name not in only:
            continue
        shells, names = shells_ext(which)
        model = sk.build(sk.rocksalt(A), shells=shells, names=names)
        case = F.Case("HSE06", model)
        T = case.targets()
        for mt in T.masses:
            mt.rel_sigma /= mass_scale
        if name.startswith("Hn"):
            T.sigma = np.where(~T.floor & (T.sigma > F.SIG_NEAR), 0.5, T.sigma)
        xp, _ = skfit.project_hr(model, case.hm)
        p = dict(zip(model.names, xp))
        ep = p["Ep"]
        for k in ("Ep", "Et2g", "Eeg"):
            p[k] -= ep
        x0s = {"H-fit": model.x_from(prev), "proj": model.x_from(p)}
        rec, _ = F.run(model, case, T, ("Ep",), x0s, name)
        rec["n_params_free"] = len(model.names) - 1
        res[name] = rec
    f = OUT / "sk_fit_extensions.json"
    old = json.loads(f.read_text()) if f.exists() else {}
    f.write_text(json.dumps({**old, **res}, indent=1, default=float))
    print("wrote", OUT / "sk_fit_extensions.json")


if __name__ == "__main__":
    main()
