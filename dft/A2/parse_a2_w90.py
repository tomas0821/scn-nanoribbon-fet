"""A2-2 analysis: pre-registered criteria 3–5 (dft/A2/README.md) for the sheet Wannier models.

  3. frozen top ≥ 0.5 eV below the vacuum level: PBE from pbe_scf_pot.cube (set by rule), HSE checked here against
     hse_pot.cube (pp.x on the HSE density);
  4. PBE twin `pbe_k16` against `pbe_bands` (Γ–X̄–M̄–Γ) inside the frozen window: max ≤ 10 meV, RMS ≤ 3 meV;
     WF centres within 0.05 Å of their atoms;
  5. HSE `hse_k16` against the HSE 16×16 IBZ eigenvalues inside the frozen window: ≤ 1 meV.

    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P dft/A2/parse_a2_w90.py
Writes dft/A2/a2_w90.json.
"""
import importlib.util
import json
import re
from pathlib import Path

import numpy as np

from scnfet import w90

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("parse_a2", HERE / "parse_a2.py")
P = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P)
spec = importlib.util.spec_from_file_location("make_a2", HERE / "make_a2.py")
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)
CELL = np.diag([M.A, M.A, M.C])


def centres(seed):
    txt = (HERE / "w90" / seed / f"{seed}.wout").read_text()
    blk = txt[txt.rfind("Final State"):]
    c = np.array([[float(v) for v in m] for m in re.findall(
        r"WF centre and spread\s+\d+\s+\(\s*(-?[\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+)\s*\)\s+([\d.]+)", blk)])
    om = float(re.search(r"Omega Total\s*=\s*([\d.]+)", blk).group(1))
    z = M.zfrac()
    atoms = np.array([(0, 0, z["Sc1"]), (0.5, 0.5, z["N1"]), (0, 0, z["N2"]), (0.5, 0.5, z["Sc2"])]) @ CELL
    inv = np.linalg.inv(CELL)
    off, dz = [], []
    for r in c[:, :3]:
        d = []
        for a in atoms:
            f = (r - a) @ inv
            f -= np.round(f)                     # wrap all three: W90 reports z modulo c (e.g. -10.14 for 12.36)
            d.append(f @ CELL)
        j = int(np.argmin([np.linalg.norm(v) for v in d]))
        off.append(float(np.linalg.norm(d[j])))
        dz.append(float(d[j][2]))
    return dict(max_offset_A=float(max(off)), offsets_A=[round(v, 4) for v in off], dz_A=[round(v, 4) for v in dz],
                spreads_A2=c[:, 3].round(4).tolist(), omega_total=om)


def interp(seed, kcart):
    hr = w90.read_hr(HERE / "w90" / seed / f"{seed}_hr.dat")
    return np.linalg.eigvalsh(w90.hk_frac(hr, kcart @ CELL.T / (2 * np.pi)))


def main():
    s1 = json.loads((HERE / "a2_stage1.json").read_text())
    win = json.loads((HERE / "a2_windows.json").read_text())
    out = {}
    # criterion 3, HSE: vacuum level from the HSE density
    cube = HERE / "hse_pot.cube"
    if cube.exists():
        ev, spread, _ = P.vacuum_level(cube)
        out["criterion3_hse_vacuum"] = dict(E_vac_eV=ev, E_vac_minus_hse_vbm=ev - s1["hse_ibz_k16"]["vbm"],
                                            margin_eV=ev - win["hse_k16"]["froz"],
                                            ok=bool(ev - win["hse_k16"]["froz"] >= 0.5))
    out["criterion3_pbe_vacuum"] = dict(margin_eV=win["pbe_k16"]["vac_margin"], ok=bool(win["pbe_k16"]["vac_margin"] >= 0.5 - 1e-9))
    # criterion 4: PBE twin on the PBE band path
    k, e, _ = P.bands(HERE / "outputs" / "pbe_bands.out")
    ref = e[:, P.NEXCL:P.NEXCL + P.NMODEL]
    E = interp("pbe_k16", k)
    inwin = ref <= win["pbe_k16"]["froz"]
    err = np.abs(E - ref)[inwin]
    vbm = s1["pbe_bands"]["vbm"]
    worst = np.unravel_index(np.argmax(np.where(inwin, np.abs(E - ref), 0)), ref.shape)
    cen = {sd: centres(sd) for sd in ("pbe_k16", "hse_k16")}
    out["criterion4"] = dict(max_meV=1e3 * float(err.max()), rms_meV=1e3 * float(np.sqrt(np.mean(err ** 2))), n=int(err.size),
                             worst=dict(k=int(worst[0]), band=int(worst[1]) + P.NEXCL + 1,
                                        E_rel_vbm=float(ref[worst] - vbm)),
                             centres_ok=bool(all(c["max_offset_A"] <= 0.05 for c in cen.values())),
                             ok=bool(err.max() <= 0.010 and np.sqrt(np.mean(err ** 2)) <= 0.003))
    for top, lab in ((0.01, "vb"), (1.4, "cb_le_cbm_plus_1"), (2.6, "le_2p6")):
        m = inwin & (ref - vbm <= top)
        out["criterion4"][f"max_{lab}_meV"] = 1e3 * float(np.abs(E - ref)[m].max())
    # criterion 5: HSE Wannier on the HSE grid
    k, e, _ = P.bands(HERE / "outputs" / "hse_ibz_k16.out")
    ref = e[:, P.NEXCL:P.NEXCL + P.NMODEL]
    E = interp("hse_k16", k)
    inwin = ref <= win["hse_k16"]["froz"]
    err = np.abs(E - ref)[inwin]
    out["criterion5"] = dict(max_meV=1e3 * float(err.max()), n=int(err.size), ok=bool(err.max() <= 0.001))
    out["centres"] = cen
    (HERE / "a2_w90.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "centres"}, indent=1))
    print({k: (v["max_offset_A"], v["omega_total"]) for k, v in cen.items()})


if __name__ == "__main__":
    main()
