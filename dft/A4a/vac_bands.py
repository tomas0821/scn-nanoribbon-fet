"""Gate M criterion 5 follow-up: per-band eigenvalue shift at Γ between nm_<key> (15 Å) and vac20_<key> (20 Å).
Run from the project root:  PYTHONPATH=. ~/.venvs/scnfet/bin/python -P dft/A4a/vac_bands.py r100_3h r100_3 r110_4
"""
import sys, numpy as np, importlib.util
from pathlib import Path
spec = importlib.util.spec_from_file_location("pa", "dft/A4a/parse_a4a.py"); pa = importlib.util.module_from_spec(spec); spec.loader.exec_module(pa)
P = pa.P
for key in sys.argv[1:]:
    n = pa.nocc_of(f"nm_{key}")
    E = {}
    for tag in (f"nm_{key}", f"vac20_{key}"):
        k, e, _ = P.bands(pa.out(tag))
        g = int(np.argmin(np.linalg.norm(k, axis=1))); assert np.linalg.norm(k[g]) < 1e-6
        E[tag] = e[g]
    a, b = E[f"nm_{key}"], E[f"vac20_{key}"]
    print(f"== {key}  nocc {n}  nbnd {len(a)}")
    d = b - a
    refs = {"VBM": d[n-1], "lowest band": d[0], "mean of lowest 12": d[:12].mean(), "median all occ": np.median(d[:n])}
    for r, v in refs.items(): print(f"   ref {r:18s} shift {1e3*v:+8.1f} meV")
    print(f"   spread of occ-band shifts (bands 1..n-3): std {1e3*d[:n-2].std():.1f} meV, min {1e3*d[:n-2].min():+.1f} max {1e3*d[:n-2].max():+.1f}")
    ref = np.median(d[:n])
    for i in list(range(n-6, n+3)):
        print(f"   band {i+1:4d} E15−VBM {a[i]-a[n-1]:+7.3f}  Δ(median-aligned) {1e3*(d[i]-ref):+7.1f} meV   Δ(VBM-aligned) {1e3*(d[i]-d[n-1]):+7.1f}")
