"""A2-1 analysis: the pre-registered criteria 1–2 (dft/A2/README.md), the vacuum level, the band survey that sets the
Wannier windows (criterion 3), and the physics flags (gaps, M̄ valley offset).

    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P dft/A2/parse_a2.py        (after copying outputs/ and pbe_scf_pot.cube)
Writes dft/A2/a2_stage1.json.
"""
import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RY = 13.605693122994
BOHR = 0.529177210903
NV = 16          # occupied bands: 2 × (Sc 3s + 3p) + 2 × N 2s + 6 N p
NEXCL = 10       # excluded from W90 (semicore + N 2s)
NMODEL = 16      # 2 × Sc d + 2 × N p


def alat_of(txt):
    return float(re.search(r"lattice parameter \(alat\)\s*=\s*([\d.]+)", txt).group(1)) * BOHR


def bands(path):
    txt = Path(path).read_text()
    a = alat_of(txt)
    i = max(txt.rfind("End of self-consistent calculation"), txt.rfind("End of band structure calculation"))
    ks, es = [], []
    for m in re.finditer(r"k =\s*(-?[\d.]+)\s*(-?[\d.]+)\s*(-?[\d.]+) \(\s*\d+ PWs\)\s+bands \(ev\):\s+(.*?)(?=\n\s*\n\s*(?:occupation|k =|the Fermi|highest|Writing|$))", txt[i:], re.S):
        ks.append([float(m.group(j)) for j in (1, 2, 3)])
        es.append([float(v) for v in re.findall(r"-?\d+\.\d+", m.group(4))])
    return np.array(ks) * 2 * np.pi / a, np.array(es), a


def levels(path):
    """Valence and conduction edges at Γ, X̄ and M̄ (whichever are in the k list) and the global gap."""
    k, e, a = bands(path)
    g = np.pi / a
    vbm, cbm = e[:, NV - 1].max(), e[:, NV].min()
    out = dict(vbm=float(vbm), gap=float(cbm - vbm), k_vbm=(k[np.argmax(e[:, NV - 1])] / g).round(3).tolist(),
               k_cbm=(k[np.argmin(e[:, NV])] / g).round(3).tolist())
    for lab, kk in (("G", (0, 0, 0)), ("X", (1, 0, 0)), ("M", (1, 1, 0))):
        d = np.linalg.norm(np.abs(k[:, :2] / g) - np.array(kk[:2]), axis=1)
        if d.min() < 1e-3:
            j = int(np.argmin(d))
            out[f"CB_{lab}"] = float(e[j, NV] - vbm)
            out[f"VB_{lab}"] = float(e[j, NV - 1] - vbm)
    if "CB_M" in out and "CB_G" in out:
        out["Mbar_offset"] = out["CB_M"] - out["CB_G"]
    return out


def vacuum_level(cube):
    lines = Path(cube).read_text().splitlines()
    nat = abs(int(lines[2].split()[0]))
    n = [int(lines[3 + i].split()[0]) for i in range(3)]
    vals = np.array(" ".join(lines[6 + nat:]).split(), float).reshape(n)        # x slowest, z fastest
    vz = vals.mean(axis=(0, 1)) * RY
    nz = n[2]
    edge = np.r_[vz[: nz // 20], vz[-nz // 20:]]
    return float(edge.mean()), float(edge.max() - edge.min()), vz


def survey(path, vbm):
    """Window survey (eV from the VBM): the model bands are NEXCL+1 … NEXCL+NMODEL; report the first non-model band's
    minimum and how far the band list reaches."""
    _, e, _ = bands(path)
    first_out = e[:, NEXCL + NMODEL] - vbm
    return dict(model_band_max=float(e[:, NEXCL + NMODEL - 1].max() - vbm), first_nonmodel_min=float(first_out.min()),
                top_band_min=float(e[:, -1].min() - vbm), nbnd=int(e.shape[1]), nk=int(e.shape[0]))


def main():
    o = HERE / "outputs"
    out = {}
    for tag in ("pbe_scf", "pbe_bands", "pbe_nscf_k16", "pbe_c30", "hse_k8_q4", "hse_k8_q8", "hse_ibz_k16"):
        f = o / f"{tag}.out"
        if f.exists() and "JOB DONE" in f.read_text():
            out[tag] = levels(f)
    # criterion 1: q 4 → 8 at k 8×8 (gap at Γ, M̄ offset, X̄ CB)
    q4, q8 = out.get("hse_k8_q4"), out.get("hse_k8_q8")
    if q4 and q8:
        d = {k: 1e3 * (q8[k] - q4[k]) for k in ("gap", "CB_G", "CB_X", "CB_M", "Mbar_offset") if k in q4 and k in q8}
        out["criterion1_q"] = dict(delta_meV=d, ok=bool(all(abs(v) <= 10 for v in d.values())))
    # criterion 2: vacuum 22.5 → 30 Å (PBE gap and M̄ offset; both on 12×12 vs 16×16 grids, Γ and M̄ on both)
    a, b = out.get("pbe_scf"), out.get("pbe_c30")
    if a and b:
        d = {k: 1e3 * (b[k] - a[k]) for k in ("gap", "CB_G", "CB_M", "Mbar_offset") if k in a and k in b}
        out["criterion2_vacuum"] = dict(delta_meV=d, ok=bool(all(abs(v) <= 5 for v in d.values())))
    # vacuum level and window survey (criterion 3 input)
    cube = HERE / "pbe_scf_pot.cube"
    if cube.exists() and "pbe_scf" in out:
        ev, spread, _ = vacuum_level(cube)
        out["vacuum"] = dict(E_vac_eV=ev, plateau_spread_meV=1e3 * spread, E_vac_minus_pbe_vbm=ev - out["pbe_scf"]["vbm"])
    for tag in ("pbe_nscf_k16", "hse_ibz_k16"):
        if tag in out:
            out[f"survey_{tag}"] = survey(o / f"{tag}.out", out[tag]["vbm"])
    (HERE / "a2_stage1.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
