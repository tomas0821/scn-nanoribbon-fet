"""Band-by-band character near the gap at Γ and X̄ from a Gate M projwfc run: E − VBM, edge-row and H weight, |ψ|²,
and the orbital shares (Sc s/p/d, N s/p, H s) of the printed Löwdin weight. Descriptive (Gate M criterion 4 follow-up).
Run from the project root:  PYTHONPATH=. ~/.venvs/scnfet/bin/python -P dft/A4a/edge_bands.py r110_4h [r110_4 ...]
"""
import importlib.util
import re
import sys

import numpy as np

spec = importlib.util.spec_from_file_location("pa", "dft/A4a/parse_a4a.py")
pa = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pa)
LNAME = "spdf"


def states(tag):
    txt = pa.out(tag).read_text()
    return {int(s): (sp.strip(), LNAME[int(l)]) for s, a, sp, l in
            re.findall(r"state #\s*(\d+): atom\s+(\d+) \((\w+)\s*\), wfc\s+\d+ \(l=(\d)", txt)}


def comps(tag):
    """per k: list of (E, {(species, l): weight}, {atom: weight}, psi2)."""
    txt = pa.out(tag).read_text()
    stmap = states(tag)
    atom = {int(s): int(a) - 1 for s, a in re.findall(r"state #\s*(\d+): atom\s+(\d+)", txt)}
    res = []
    for kb in re.split(r"\n k =", txt.split("Lowdin Charges")[0])[1:]:
        st = []
        for m in re.finditer(r"==== e\(\s*\d+\) =\s*(-?[\d.]+) eV ====(.*?)\|psi\|\^2 = ([\d.]+)", kb, re.S):
            wo, wa = {}, {}
            for c, s in re.findall(r"([\d.]+)\*\[#\s*(\d+)\]", m.group(2)):
                key = stmap[int(s)]
                wo[key] = wo.get(key, 0.0) + float(c)
                wa[atom[int(s)]] = wa.get(atom[int(s)], 0.0) + float(c)
            st.append((float(m.group(1)), wo, wa, float(m.group(3))))
        res.append(st)
    return res


for key in sys.argv[1:]:
    nocc = pa.nocc_of(f"nm_{key}")
    left, right, hm = pa.edge_masks(key)
    kp, ep, _ = pa.P.bands(pa.out(f"bands_{key}"))
    vbm = ep[:, nocc - 1].max()
    cw = comps(f"projwfc_{key}")
    print(f"=== {key}  nocc {nocc}  VBM {vbm:.3f} eV  (edge = outermost host rows; shares of the printed weight)")
    for ik, lab in ((0, "Γ"), (len(kp) - 1, "X̄")):
        print(f"  -- {lab}")
        print("   band  E−VBM   edge    H   edge+H  |ψ|²   Sc s  Sc p  Sc d   N s   N p   H s")
        for b in range(nocc - 8, nocc + 6):
            E, wo, wa, p2 = cw[ik][b]
            assert abs(E - ep[ik, b]) < 2e-3
            ew = pa.edge_weight(wa, left, right, hm)
            tot = sum(wo.values())
            sh = lambda s, l: wo.get((s, l), 0.0) / tot          # noqa: E731
            mark = " ←VBM" if b == nocc - 1 else (" ←CB1" if b == nocc else "")
            print(f"   {b + 1:4d} {E - vbm:+7.3f}  {ew['edge']:.2f}  {ew['H']:.2f}  {ew['edge_plus_H']:.2f}   {p2:.2f}   "
                  f"{sh('Sc', 's'):.2f}  {sh('Sc', 'p'):.2f}  {sh('Sc', 'd'):.2f}  {sh('N', 's'):.2f}  {sh('N', 'p'):.2f}  "
                  f"{sh('H', 's'):.2f}{mark}")
