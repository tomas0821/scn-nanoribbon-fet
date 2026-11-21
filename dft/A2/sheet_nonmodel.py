"""A4b amendment 2026-11-14, item 6 (diagnostic only): does the M-ScN sheet have a vacuum-like (non-model) empty state
near its CBM? projwfc.x on the A2 HSE06 SCF (IBZ 16×16) and PBE NSCF (16×16), run_a2_proj.sh.

Same measure as A4b: s_out = |ψ|² × (Löwdin share outside Sc d + N p) + (1 − |ψ|²); non-model if ≥ 0.25 (0.20 /
0.30 reported). Reported, never used to change an A2 or gate (a) result.

    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P dft/A2/sheet_nonmodel.py      → dft/A2/a2_sheet_nonmodel.json
"""
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location("parse_a4b", ROOT / "dft" / "A4b" / "parse_a4b.py")
PA = importlib.util.module_from_spec(spec)
argv, sys.argv = sys.argv, ["parse_a4b.py"]
spec.loader.exec_module(PA)
sys.argv = argv

NOCC = 16                                    # 2 f.u. × (Sc 11 + N 5) / 2


def analyse(tag):
    S = PA.projwfc_states(tag, HERE / "outputs" / f"{tag}.out")
    vbm = max(st[NOCC - 1][0] for st in S)
    cbm = min(st[NOCC][0] for st in S)
    rows = []
    for ik, states in enumerate(S):
        for b in range(NOCC, len(states)):
            E, w, p2 = states[b]
            tot = sum(w.values())
            out = sum(v for (a, sp, l), v in w.items() if (sp, l) not in (("Sc", "d"), ("N", "p"))) / tot if tot else 1.0
            byl = {}
            for (a, sp, l), v in w.items():
                byl[sp + " " + l] = byl.get(sp + " " + l, 0.0) + v / (tot or 1.0)
            rows.append(dict(E_rel_cbm=E - cbm, k_index=ik, band=b + 1, psi2=p2, s_out=p2 * out + (1 - p2),
                             top=dict(sorted(((k, round(v, 3)) for k, v in byl.items()), key=lambda x: -x[1])[:4])))
    rows.sort(key=lambda r: r["E_rel_cbm"])
    first = {f"{t:.2f}": next((r for r in rows if r["s_out"] >= t), None) for t in (0.20, 0.25, 0.30)}
    return dict(nk=len(S), vbm=vbm, cbm=cbm, gap=cbm - vbm, first_nonmodel=first,
                lowest_empty=[r for r in rows if r["E_rel_cbm"] <= 1.0][:12])


def main():
    res = {}
    for tag, lab in (("hproj_sheet", "HSE06"), ("pproj_sheet", "PBE")):
        p = HERE / "outputs" / f"{tag}.out"
        if not p.exists():
            print(f"{lab}: {p.name} missing")
            continue
        r = res[lab] = analyse(tag)
        print(f"=== sheet {lab}: {r['nk']} k, gap {r['gap']:.3f} eV")
        for t, fn in r["first_nonmodel"].items():
            print(f"   first non-model (s_out ≥ {t}): " + ("none" if fn is None else
                  f"CBM + {fn['E_rel_cbm']:.3f} eV, k {fn['k_index']}, band {fn['band']}, |ψ|² {fn['psi2']:.2f}, {fn['top']}"))
    (HERE / "a2_sheet_nonmodel.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
