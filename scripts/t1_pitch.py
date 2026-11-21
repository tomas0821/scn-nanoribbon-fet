"""Pitch-normalised I_ON (fits/PREREG_t1b.md §1). Post-processing of fits/t1/corr/ only; no device solve.

    python3 scripts/t1_pitch.py
Writes fits/t1/pitch/pitch.json and fits/t1/pitch/summary.md.

I_fp = I_ON(W) × W / (W + s), s edge to edge. s_min = 2.0 nm (the DFT-cell gap); s_GAA = 2 t_HfO2 at EOT 0.6
(oxide shells touching). GAA is reported at s_GAA only; parallel plate and the single gate at both spacings.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORR = ROOT / "fits" / "t1" / "corr"
OUT = ROOT / "fits" / "t1" / "pitch"
PREREG = ROOT / "fits" / "PREREG_t1b.md"
IRDS = 1.979
EOT = 0.6
T_HFO2 = EOT * 25.0 / 3.9                       # nm, the scnfet/cg.py stack
SPACINGS = {"s2": 2.0, "sGAA": 2 * T_HFO2}
GEOMS = {"s2": ("parallel_plate", "single"), "sGAA": ("parallel_plate", "single", "gaa")}
PER_W_GEOMS = ("parallel_plate", "single", "gaa")
FAMS = ("100", "110", "110h")
LAB = {"100": "[100] bare", "110": "[110] bare", "110h": "[110] H (in-sample only)",
       "100h": "[100] H (DFT point)"}
KINDS = ("ballistic", "B2d")


def load(fam):
    rs = [json.loads(p.read_text()) for p in CORR.glob(f"{fam}_n*.json")]
    return sorted(rs, key=lambda r: r["n"])


def point(W, carriers):
    """Per-W and per-footprint I_ON for one ribbon (all carriers, geometries, kinds)."""
    out = {"per_W": {}, "fp": {}}
    for c, x in carriers.items():
        out["per_W"][c] = {g: {k: x["I_on"][f"{g} EOT 0.6"][k] for k in KINDS} for g in PER_W_GEOMS}
        out["fp"][c] = {}
        for sk, s in SPACINGS.items():
            ff = W / (W + s)
            out["fp"][c][sk] = {"fill": ff, **{g: {k: x["I_on"][f"{g} EOT 0.6"][k] * ff for k in KINDS}
                                               for g in GEOMS[sk]}}
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    md5 = hashlib.md5(PREREG.read_bytes()).hexdigest()
    res = {"prereg_md5": md5, "EOT": EOT, "t_HfO2_nm": T_HFO2, "spacings_nm": SPACINGS, "families": {}, "dft": {}}
    for f in FAMS:
        res["families"][f] = [{"n": r["n"], "W_nm": r["W_nm"], **point(r["W_nm"], r["carriers"])} for r in load(f)]
    dft = json.loads((CORR / "dft.json").read_text())["points"]
    for p, d in dft.items():
        res["dft"][p] = {"key": d["key"], "W_nm": d["W_nm"], **point(d["W_nm"], d["carriers"])}

    # pass counts: per W (registered) and per footprint
    counts = {}
    for f, rs in res["families"].items():
        counts[f] = {}
        for c in ("n", "p"):
            row = {}
            for g in PER_W_GEOMS:
                for k in KINDS:
                    row[f"per_W {g} {k}"] = sum(r["per_W"][c][g][k] >= IRDS for r in rs)
            for sk in SPACINGS:
                for g in GEOMS[sk]:
                    for k in KINDS:
                        row[f"{sk} {g} {k}"] = sum(r["fp"][c][sk][g][k] >= IRDS for r in rs)
            row["total"] = len(rs)
            counts[f][c] = row
    res["irds_pass_counts"] = counts
    (OUT / "pitch.json").write_text(json.dumps(res, indent=1))

    # summary.md
    sg = SPACINGS["sGAA"]
    L = ["# T1b: pitch-normalised I_ON (`fits/PREREG_t1b.md` §1)", "",
         f"I_fp = I_ON × W / (W + s), mA per μm of footprint, EOT 0.6. s = 2.0 nm (DFT-cell gap) and s = {sg:.2f} nm "
         f"(2 t_HfO2; oxide shells touching; GAA only here). × B = 2D Fröhlich in vacuum. IRDS ≥ {IRDS}.", "",
         "## IRDS pass counts (ribbons passing / total), × B", "",
         "| family | carrier | per W: pp / single / GAA | s = 2 nm: pp / single | "
         f"s = {sg:.1f} nm: pp / single / GAA | of |", "|---|---|---|---|---|---|"]
    for f in FAMS:
        for c in ("n", "p"):
            k = counts[f][c]
            L.append(f"| {LAB[f]} | {c} | {k['per_W parallel_plate B2d']} / {k['per_W single B2d']} / "
                     f"{k['per_W gaa B2d']} | {k['s2 parallel_plate B2d']} / {k['s2 single B2d']} | "
                     f"{k['sGAA parallel_plate B2d']} / {k['sGAA single B2d']} / {k['sGAA gaa B2d']} | {k['total']} |")
    L += ["", "## IRDS pass counts, ballistic (before B)", "",
          "| family | carrier | per W: pp / single / GAA | s = 2 nm: pp / single | "
          f"s = {sg:.1f} nm: pp / single / GAA | of |", "|---|---|---|---|---|---|"]
    for f in FAMS:
        for c in ("n", "p"):
            k = counts[f][c]
            L.append(f"| {LAB[f]} | {c} | {k['per_W parallel_plate ballistic']} / {k['per_W single ballistic']} / "
                     f"{k['per_W gaa ballistic']} | {k['s2 parallel_plate ballistic']} / {k['s2 single ballistic']} | "
                     f"{k['sGAA parallel_plate ballistic']} / {k['sGAA single ballistic']} / "
                     f"{k['sGAA gaa ballistic']} | {k['total']} |")
    L.append("")
    for f in FAMS:
        L += [f"## {LAB[f]}: × B, mA/μm", "",
              "| n | W (nm) | carrier | per W pp / single / GAA | fill s 2 / s 7.7 | s = 2: pp / single | "
              "s = 7.7: pp / single / GAA |", "|---|---|---|---|---|---|---|"]
        for r in res["families"][f]:
            for c in ("n", "p"):
                w, a, b = r["per_W"][c], r["fp"][c]["s2"], r["fp"][c]["sGAA"]
                L.append(f"| {r['n']} | {r['W_nm']:.2f} | {c} | {w['parallel_plate']['B2d']:.2f} / "
                         f"{w['single']['B2d']:.2f} / {w['gaa']['B2d']:.2f} | {a['fill']:.2f} / {b['fill']:.2f} | "
                         f"{a['parallel_plate']['B2d']:.2f} / {a['single']['B2d']:.2f} | "
                         f"{b['parallel_plate']['B2d']:.2f} / {b['single']['B2d']:.2f} / {b['gaa']['B2d']:.2f} |")
        L.append("")
    L += ["## Calibration DFT-Wannier ribbons, × B", "",
          "| ribbon | W (nm) | carrier | per W pp / single | s = 2: pp / single | s = 7.7: pp / single / GAA |",
          "|---|---|---|---|---|---|"]
    for p, d in res["dft"].items():
        for c in ("n", "p"):
            w, a, b = d["per_W"][c], d["fp"][c]["s2"], d["fp"][c]["sGAA"]
            L.append(f"| {d['key']} | {d['W_nm']:.2f} | {c} | {w['parallel_plate']['B2d']:.2f} / {w['single']['B2d']:.2f} | "
                     f"{a['parallel_plate']['B2d']:.2f} / {a['single']['B2d']:.2f} | {b['parallel_plate']['B2d']:.2f} / "
                     f"{b['single']['B2d']:.2f} / {b['gaa']['B2d']:.2f} |")
    (OUT / "summary.md").write_text("\n".join(L) + "\n")
    print(f"wrote {OUT / 'pitch.json'} and summary.md (prereg md5 {md5})")


if __name__ == "__main__":
    main()
