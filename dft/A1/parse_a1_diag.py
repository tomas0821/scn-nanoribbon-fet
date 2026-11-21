"""Score the A1 Wannier diagnostic: criterion 4 (PBE-twin interpolation inside the frozen window, max ≤ 10 meV,
RMS ≤ 3 meV) for pbe_gi, pbe_ws, pbe9 and pbe9_ws, from Wannier90's geninterp output on the QE band path.

    ~/.venvs/scnfet/bin/python -P dft/A1/parse_a1_diag.py      (after copying w90/<seed>/<seed>_geninterp.dat back)
"""
import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FROZ = 5.4


def qe_bands():
    txt = (HERE / "outputs" / "pbe_scf_bands.out").read_text()
    i = txt.rfind("End of band structure calculation")
    es = [[float(v) for v in re.findall(r"-?\d+\.\d+", m)] for m in
          re.findall(r"bands \(ev\):\s+(.*?)(?=\n\s*\n\s*(?:k =|Writing|$))", txt[i:], re.S)]
    return np.array(es)


def geninterp(seed):
    """Bands of `seed` at the QE path k-points, from _hr.dat (+ _wsvec.dat when use_ws_distance is on).
    (geninterp is a postw90.x module, so wannier90.x ignored it; scnfet.w90 reproduces Wannier90's own band output
    to 0.04 meV, the _hr.dat print precision.)"""
    import sys
    sys.path.insert(0, str(HERE.parents[1]))
    from scnfet import w90
    d = HERE / "w90" / seed
    hr = w90.read_hr(d / f"{seed}_hr.dat")
    if "use_ws_distance = .true." in (d / f"{seed}.win").read_text():
        hr = w90.apply_wsvec(hr, w90.read_wsvec(d / f"{seed}_wsvec.dat"))
    kf = np.array([[float(v) for v in l.split()[1:4]] for l in (d / f"{seed}_geninterp.kpt").read_text().splitlines()[3:]])
    return np.linalg.eigvalsh(w90.hk_frac(hr, kf))


def main():
    eq = qe_bands()
    vbm = json.loads((HERE / "a1_stage1.json").read_text())["survey"]["pbe_nscf"]["vbm"]
    ref = eq[:, 5:13]                                   # QE bands 6–13 = the 8 Sc d + N p bands
    inwin = ref <= vbm + FROZ
    out = {}
    seeds = [x for x in ("pbe_gi", "pbe_ws", "pbe9", "pbe9_ws", "pbe_mlwf", "pbe_k12", "pbe_k16")
             if (HERE / "w90" / x / f"{x}_hr.dat").exists()]
    for seed in seeds:
        try:
            E = geninterp(seed)
        except FileNotFoundError:
            out[seed] = "missing"
            continue
        err = np.abs(E[:, :8] - ref)[inwin]
        worst = np.unravel_index(np.argmax(np.where(inwin, np.abs(E[:, :8] - ref), 0)), ref.shape)
        out[seed] = dict(max_meV=round(1e3 * float(err.max()), 1), rms_meV=round(1e3 * float(np.sqrt(np.mean(err ** 2))), 2),
                         n_over_20meV=int((err > 0.02).sum()), worst=dict(k=int(worst[0]), band=int(worst[1]) + 6),
                         passes=bool(err.max() <= 0.010 and np.sqrt(np.mean(err ** 2)) <= 0.003))
        for top, lab in ((0.01, "vb"), (2.5, "le_2p5"), (4.0, "le_4p0")):
            m = inwin & (ref <= vbm + top) & ((ref > vbm + 0.01) if lab == "le_2p5" else True)
            e = np.abs(E[:, :8] - ref)[m]
            out[seed][f"max_{lab}_meV"] = round(1e3 * float(e.max()), 1)
        if seed == "pbe9":
            s_band = E[:, 8] - vbm
            out[seed]["ninth_band_rel_vbm_min_max"] = [round(float(s_band.min()), 2), round(float(s_band.max()), 2)]
    (HERE / "a1_diag.json").write_text(json.dumps(out, indent=1))
    for k, v in out.items():
        print(k, v)


if __name__ == "__main__":
    main()
