"""Parse the A0 cutoff/pseudopotential tests and apply the criteria pre-registered in README.md (2026-11-10).

Run from the project root (laptop; the outputs were copied from the cluster):
    ~/.venvs/scnfet/bin/python -P dft/A0/parse_a0.py
Writes dft/A0/a0_summary.json and figures/A0_cutoff_convergence.png, and prints the verdict tables.

Each criterion compares a cutoff with the highest cutoff of the same set (paw110, nc110):
  bulk BM3 a0 ±0.002 Å; E(V) shape ±1 meV/atom (energies relative to the a = 4.51 Å point); bulk Γ–X gap at
  a = 4.50 Å ±10 meV; sheet gap at fixed geometry ±10 meV; PAW sheet 2Dxy vc-relax a ±0.003 Å; pressure at
  a = 4.50 Å ±1 kbar. Production = the lowest passing cutoff.
"""
import json
import re
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "outputs"
RY = 13.605693122994
BOHR = 0.529177210903
GPA = 160.21766208            # eV/Å³ → GPa
A_LIST = [4.42, 4.45, 4.48, 4.51, 4.54, 4.57, 4.60]
SETS = {"paw": [60, 75, 90, 110], "nc": [80, 90, 110]}
TOL = dict(a0=0.002, shape=1.0, gap_bulk=10.0, gap_sheet=10.0, a_vc=0.003, P=1.0)   # Å, meV/atom, meV, meV, Å, kbar


def parse(tag):
    txt = (OUT / f"{tag}.out").read_text()
    if "JOB DONE" not in txt:
        raise RuntimeError(f"{tag}: no JOB DONE")
    if "convergence NOT achieved" in txt:
        raise RuntimeError(f"{tag}: SCF not converged")
    r = dict(tag=tag)
    alat = float(re.search(r"lattice parameter \(alat\)\s*=\s*([\d.]+)", txt).group(1))
    r["alat_bohr"] = alat
    r["nelec"] = float(re.search(r"number of electrons\s*=\s*([\d.]+)", txt).group(1))
    E = re.findall(r"^!\s+total energy\s+=\s+(-?[\d.]+) Ry", txt, re.M)
    r["E_Ry"] = float(E[-1]) if E else None
    P = re.findall(r"P=\s*(-?[\d.]+)", txt)
    r["P_kbar"] = float(P[-1]) if P else None
    # eigenvalues of the last SCF / band block
    i = max(txt.rfind("End of self-consistent calculation"), txt.rfind("End of band structure calculation"))
    block = txt[i:]
    ks, eig = [], []
    for m in re.finditer(r"k =\s*(-?[\d.]+)\s*(-?[\d.]+)\s*(-?[\d.]+) \(\s*\d+ PWs\)\s+bands \(ev\):\s+(.*?)(?:\n\s*\n\s*(?:occupation|k =|the Fermi|highest|Writing|\s*$))",
                         block, re.S):
        ks.append([float(m.group(j)) for j in (1, 2, 3)])
        eig.append([float(v) for v in re.findall(r"-?\d+\.\d+", m.group(4))])
    r["k"] = np.array(ks)
    r["eig"] = np.array(eig)
    # final cell of a vc-relax (alat units)
    if "Begin final coordinates" in txt:
        fin = txt[txt.rfind("Begin final coordinates"):]
        m = re.search(r"CELL_PARAMETERS \(alat=\s*([\d.]+)\)\s*\n(.*)\n(.*)\n(.*)\n", fin)
        a_alat = float(m.group(1))
        cell = np.array([[float(v) for v in m.group(j).split()] for j in (2, 3, 4)]) * a_alat * BOHR
        r["cell_A"] = cell
    return r


def bm3(V, E0, V0, B0, Bp):
    x = (V0 / V) ** (2 / 3)
    return E0 + 9 * V0 * B0 / 16 * ((x - 1) ** 3 * Bp + (x - 1) ** 2 * (6 - 4 * x))


def fit_bm3(a, E_eV):
    V = np.asarray(a) ** 3 / 4                       # fcc primitive cell
    c = np.polyfit(V, E_eV, 2)
    V0 = -c[1] / (2 * c[0])
    p0 = [np.polyval(c, V0), V0, 2 * c[0] * V0, 4.0]
    p, _ = curve_fit(bm3, V, E_eV, p0=p0, maxfev=20000)
    return dict(a0=(4 * p[1]) ** (1 / 3), B0_GPa=p[2] * GPA, Bp=p[3], rms_meV=1e3 * float(np.sqrt(np.mean((bm3(V, *p) - E_eV) ** 2))))


def bulk_gap(r):
    """Γ15 → X3 for 16 electrons: band 8 at Γ (index 7) to band 9 at X (index 8). k in 2π/alat Cartesian."""
    nocc = int(round(r["nelec"] / 2))
    k, e = r["k"], r["eig"]
    iG = np.argmin(np.linalg.norm(k, axis=1))
    isX = np.isclose(np.sort(np.abs(k), axis=1), [0, 0, 1], atol=1e-4).all(axis=1)
    if not isX.any():
        raise RuntimeError(f"{r['tag']}: X not in the k list")
    return dict(gap_GX=e[isX, nocc].min() - e[iG, nocc - 1], vbm_at_G=bool(np.isclose(e[:, nocc - 1].max(), e[iG, nocc - 1])),
                cbm_at_X=bool(np.isclose(e[:, nocc].min(), e[isX, nocc].min())))


def sheet_gap(r):
    nocc = int(round(r["nelec"] / 2))
    e, k = r["eig"], r["k"]
    vbm, cbm = e[:, nocc - 1].max(), e[:, nocc].min()
    iG = np.argmin(np.linalg.norm(k, axis=1))
    return dict(gap=cbm - vbm, gap_G=e[iG, nocc] - e[iG, nocc - 1], vbm_k=k[np.argmax(e[:, nocc - 1])].tolist(),
                cbm_k=k[np.argmin(e[:, nocc])].tolist())


def main():
    res = {}
    for fam, cuts in SETS.items():
        for c in cuts:
            s = f"{fam}{c}"
            ev = [parse(f"bulk_{s}_a{a:.2f}") for a in A_LIST]
            E = np.array([q["E_Ry"] for q in ev]) * RY
            a_out = [q["alat_bohr"] * BOHR for q in ev]
            assert np.allclose(a_out, A_LIST, atol=1e-4), (s, a_out)
            fit = fit_bm3(A_LIST, E)
            b = parse(f"bulk_{s}_a4.50")
            sh = parse(f"sheet_{s}_fixed")
            rec = dict(set=s, family=fam, ecut=c, **fit, E_rel_meV_per_atom=(1e3 * (E - E[3]) / 2).tolist(),
                       P_kbar=b["P_kbar"], **bulk_gap(b), sheet=sheet_gap(sh), sheet_P_kbar=sh["P_kbar"])
            if fam == "paw":
                vc = parse(f"sheet_{s}_vcrelax")
                rec["a_vc"] = float(np.linalg.norm(vc["cell_A"][0]))
            res[s] = rec
    # k ladder at a = 4.50
    ladder = {}
    for s in ("paw90", "nc80"):
        ladder[s] = {}
        for kk, tag in ((8, f"bulk_{s}_a4.50_k8"), (12, f"bulk_{s}_a4.50"), (16, f"bulk_{s}_a4.50_k16"), (20, f"bulk_{s}_a4.50_k20")):
            q = parse(tag)
            ladder[s][kk] = dict(E_Ry=q["E_Ry"], P_kbar=q["P_kbar"], **bulk_gap(q))
    # sheet band paths
    paths = {}
    for s in ("paw90", "nc80"):
        q = parse(f"sheet_{s}_fixed_bands")
        nocc = int(round(q["nelec"] / 2))
        paths[s] = dict(k=q["k"].tolist(), vb=q["eig"][:, nocc - 1].tolist(), cb=q["eig"][:, nocc].tolist(),
                        gap=float(q["eig"][:, nocc].min() - q["eig"][:, nocc - 1].max()))

    # criteria vs the highest cutoff of each family
    verdict = {}
    for fam, cuts in SETS.items():
        ref = res[f"{fam}{cuts[-1]}"]
        for c in cuts:
            r = res[f"{fam}{c}"]
            d = dict(a0=abs(r["a0"] - ref["a0"]),
                     shape=float(np.max(np.abs(np.array(r["E_rel_meV_per_atom"]) - np.array(ref["E_rel_meV_per_atom"])))),
                     gap_bulk=1e3 * abs(r["gap_GX"] - ref["gap_GX"]), gap_sheet=1e3 * abs(r["sheet"]["gap"] - ref["sheet"]["gap"]),
                     P=abs(r["P_kbar"] - ref["P_kbar"]))
            if fam == "paw":
                d["a_vc"] = abs(r["a_vc"] - ref["a_vc"])
            ok = {k: v <= TOL[k] + 1e-12 for k, v in d.items()}
            verdict[r["set"]] = dict(dev=d, ok=ok, passes=all(ok.values()), reference=(c == cuts[-1]))
    choice = {}
    for fam, cuts in SETS.items():
        passing = [c for c in cuts[:-1] if verdict[f"{fam}{c}"]["passes"]]
        choice[fam] = min(passing) if passing else cuts[-1]
    physics = dict(
        bulk_a0_PAW_ref=res["paw110"]["a0"], bulk_a0_NC_ref=res["nc110"]["a0"],
        gap_GX_PAW_ref=res["paw110"]["gap_GX"], gap_GX_NC_ref=res["nc110"]["gap_GX"],
        sheet_gap_PAW_ref=res["paw110"]["sheet"]["gap"], sheet_gap_NC_ref=res["nc110"]["sheet"]["gap"],
        paw_vs_nc_sheet_gap_meV=1e3 * abs(res["paw110"]["sheet"]["gap"] - res["nc110"]["sheet"]["gap"]),
        k_ladder_gap_k12_vs_k20_meV={s: 1e3 * abs(ladder[s][12]["gap_GX"] - ladder[s][20]["gap_GX"]) for s in ladder},
        k_ladder_E_k12_vs_k20_meV_per_atom={s: 1e3 * abs(ladder[s][12]["E_Ry"] - ladder[s][20]["E_Ry"]) * RY / 2 for s in ladder})

    summary = dict(results=res, verdict=verdict, choice=choice, physics=physics, ladder=ladder, sheet_paths=paths, tolerances=TOL)
    (HERE / "a0_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    report(res, verdict, choice, physics, ladder)
    plot(res, verdict)


def report(res, verdict, choice, physics, ladder):
    print("set     a0(Å)    B0(GPa) fitRMS  ΓX(eV)  sheetEg  P(kbar)  a_vc(Å) | Δa0     shape   ΔgapB  ΔgapS  ΔP    Δa_vc  pass")
    for s, r in res.items():
        v = verdict[s]
        d = v["dev"]
        print(f"{s:7s} {r['a0']:.4f}  {r['B0_GPa']:6.1f}  {r['rms_meV']:5.2f}  {r['gap_GX']:+.3f}  {r['sheet']['gap']:.3f}  "
              f"{r['P_kbar']:7.2f}  {r.get('a_vc', float('nan')):.4f} | {d['a0']:.4f}  {d['shape']:5.2f}  {d['gap_bulk']:5.1f}  "
              f"{d['gap_sheet']:5.1f}  {d['P']:5.2f}  {d.get('a_vc', float('nan')):.4f}  "
              f"{'REF' if v['reference'] else ('yes' if v['passes'] else 'no: ' + ','.join(k for k, o in v['ok'].items() if not o))}")
    print("choice:", choice)
    print("physics:", json.dumps(physics, indent=1, default=float))
    print("k ladder Γ–X (eV):", {s: {k: round(q["gap_GX"], 4) for k, q in L.items()} for s, L in ladder.items()})
    print("VBM at Γ / CBM at X:", {s: (r["vbm_at_G"], r["cbm_at_X"]) for s, r in res.items()})
    print("sheet VBM/CBM k:", {s: (r["sheet"]["vbm_k"], r["sheet"]["cbm_k"]) for s, r in res.items() if s in ("paw110", "nc110")})


def plot(res, verdict):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, ink2, grid, band = "#17181a", "#4d4f53", "#e6e6e4", "#ecebe7"
    col = {"paw": "#2a78d6", "nc": "#eb6834"}
    lab = {"paw": "PAW (Sc psl 0.2.3), dual 8", "nc": "NC (Sc SG15), dual 4"}
    plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": ink2, "xtick.color": ink2, "ytick.color": ink2,
                         "axes.labelcolor": ink, "axes.linewidth": 0.8})
    panels = [("a0", "bulk a₀ (BM3), Δ vs reference (mÅ)", 1e3, TOL["a0"] * 1e3),
              ("shape", "E(V) shape, max |Δ| (meV/atom)", 1, TOL["shape"]),
              ("gap_bulk", "bulk Γ–X gap at 4.50 Å, |Δ| (meV)", 1, TOL["gap_bulk"]),
              ("gap_sheet", "sheet gap, fixed geometry, |Δ| (meV)", 1, TOL["gap_sheet"]),
              ("P", "pressure at 4.50 Å, |Δ| (kbar)", 1, TOL["P"]),
              ("a_vc", "sheet a, 2Dxy vc-relax (PAW), |Δ| (mÅ)", 1e3, TOL["a_vc"] * 1e3)]
    fig, axes = plt.subplots(2, 3, figsize=(10, 5.6), constrained_layout=True)
    for ax, (key, title, scale, tol) in zip(axes.ravel(), panels):
        ymax = tol * 1.6
        for fam, cuts in SETS.items():
            if key == "a_vc" and fam == "nc":
                continue
            y = [verdict[f"{fam}{c}"]["dev"][key] * scale for c in cuts]
            ymax = max(ymax, max(y) * 1.15)
            ax.plot(cuts, y, "-o", color=col[fam], lw=2, ms=6, mec="white", mew=1.2, label=lab[fam], zorder=3)
        ax.axhspan(0, tol, color=band, zorder=0)
        ax.axhline(tol, color=ink2, lw=0.8, ls=(0, (3, 3)), zorder=1)
        ax.text(0.98, tol, f" tolerance {tol:g}", transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=7.5, color=ink2)
        ax.set_ylim(0, ymax)
        ax.set_xticks([60, 75, 80, 90, 110])
        ax.set_xlabel("ecutwfc (Ry)")
        ax.set_title(title, loc="left", fontsize=9, color=ink)
        ax.grid(axis="y", color=grid, lw=0.6)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0, 0].legend(frameon=False, fontsize=7.5, loc="upper right")
    fig.suptitle("A0: cutoff convergence against the highest cutoff of each set (110 Ry); shaded = pre-registered tolerance",
                 x=0.01, ha="left", fontsize=9.5, color=ink)
    out = ROOT / "figures" / "A0_cutoff_convergence.png"
    fig.savefig(out, dpi=200, facecolor="white")
    print("wrote", out)


if __name__ == "__main__":
    main()
