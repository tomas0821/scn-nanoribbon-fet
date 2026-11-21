"""Array-C_G results × Fröhlich B, per W and per footprint (fits/PREREG_array.md). Post-processing only.

    python3 scripts/t1_array_post.py
Reads fits/t1/cg_array/{s2,sGAA}/ and fits/t1/corr/; writes fits/t1/cg_array/array.json and summary.md.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARR = ROOT / "fits" / "t1" / "cg_array"
CORR = ROOT / "fits" / "t1" / "corr"
IRDS = 1.979
FAMS = ("100", "110", "110h")
LAB = {"100": "[100] bare", "110": "[110] bare", "110h": "[110] H (in-sample only)"}


def load(d, fam):
    return {json.loads(p.read_text())["n"]: json.loads(p.read_text()) for p in d.glob(f"{fam}_n*.json")}


def main():
    md5 = hashlib.md5((ROOT / "fits" / "PREREG_array.md").read_bytes()).hexdigest()
    tab = json.loads((ARR / "cg_array_table.json").read_text())
    S = tab["spacings_nm"]
    res = dict(prereg_md5=md5, spacings_nm=S, families={}, counts={})
    for f in FAMS:
        corr = load(CORR, f)
        arrs = {k: load(ARR / k, f) for k in S}
        rows = []
        for n in sorted(corr):
            c = corr[n]
            W = c["W_nm"]
            row = dict(n=n, W_nm=W)
            for k, s in S.items():
                a = arrs[k][n]
                ff = W / (W + s)
                row[k] = {"fill": ff, "C_ratio_pp": a["tob"]["n EOT 0.6"]["C_ratio_pp"],
                          "C_ratio_iso": a["tob"]["n EOT 0.6"]["C_ratio_iso"]}
                for car in ("n", "p"):
                    I = a["tob"][f"{car} EOT 0.6"]["I_on"]
                    B2, B3 = c["carriers"][car]["frohlich"]["2d"]["B"], c["carriers"][car]["frohlich"]["3d"]["B"]
                    row[k][car] = dict(per_W=dict(ballistic=I, B2d=I * B2, B3d=I * B3),
                                       fp=dict(ballistic=I * ff, B2d=I * B2 * ff, B3d=I * B3 * ff))
            rows.append(row)
        res["families"][f] = rows
        res["counts"][f] = {f"{k} {car} {norm} {kind}": sum(r[k][car][norm][kind] >= IRDS for r in rows)
                            for k in S for car in ("n", "p") for norm in ("per_W", "fp")
                            for kind in ("ballistic", "B2d", "B3d")}
        res["counts"][f]["total"] = len(rows)
    (ARR / "array.json").write_text(json.dumps(res, indent=1))

    sg = S["sGAA"]
    L = ["# Array C_G (single gate, periodic): I_ON × B per W and per footprint (`fits/PREREG_array.md`)", "",
         f"EOT 0.6; s = 2.0 nm and s = {sg:.2f} nm (edge to edge). × B = 2D Fröhlich in vacuum (3D in the JSON). "
         f"IRDS ≥ {IRDS}. GAA at s = {sg:.2f} nm equals the isolated GAA (`fits/t1/pitch/`).", "",
         "## IRDS pass counts (ribbons passing / total)", "",
         "| family | carrier | s = 2: per W ball. / × B | s = 2: per footprint ball. / × B | "
         f"s = {sg:.1f}: per W ball. / × B | s = {sg:.1f}: per footprint ball. / × B | of |",
         "|---|---|---|---|---|---|---|"]
    for f in FAMS:
        k = res["counts"][f]
        for car in ("n", "p"):
            L.append(f"| {LAB[f]} | {car} | " + " | ".join(
                f"{k[f'{s} {car} {norm} ballistic']} / {k[f'{s} {car} {norm} B2d']}"
                for s in S for norm in ("per_W", "fp")) + f" | {k['total']} |")
    L.append("")
    for f in FAMS:
        L += [f"## {LAB[f]}", "",
              "| n | W (nm) | carrier | s = 2: C/C_pp, C/C_iso | per W × B | per footprint × B | "
              f"s = {sg:.1f}: C/C_pp, C/C_iso | per W × B | per footprint × B |", "|---|---|---|---|---|---|---|---|---|"]
        for r in res["families"][f]:
            for car in ("n", "p"):
                a, b = r["s2"], r["sGAA"]
                L.append(f"| {r['n']} | {r['W_nm']:.2f} | {car} | {a['C_ratio_pp']:.2f}, {a['C_ratio_iso']:.2f} | "
                         f"{a[car]['per_W']['B2d']:.2f} | {a[car]['fp']['B2d']:.2f} | "
                         f"{b['C_ratio_pp']:.2f}, {b['C_ratio_iso']:.2f} | {b[car]['per_W']['B2d']:.2f} | "
                         f"{b[car]['fp']['B2d']:.2f} |")
        L.append("")
    (ARR / "summary.md").write_text("\n".join(L) + "\n")
    print("wrote fits/t1/cg_array/array.json and summary.md")


if __name__ == "__main__":
    main()
