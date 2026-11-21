"""Post hoc, descriptive, not pre-registered (2026-11-20): the on-current share of the nearly-free-electron (NFE) band
of the bare [100] n = 3 calibration ribbon.

In HSE06 a vacuum-reaching NFE state sits 0.35 eV above the conduction-band minimum of [100] n = 3, inside the device
window (dft/A4b/README.md, the 14 Nov amendment). Its DFT device Hamiltonian is HW3 (Sc s on every Sc), which contains
the band; HW2 stops the frozen window 0.1 eV below it and leaves the band out. The n-FET on-current of the two
Hamiltonians, with the same ToB settings and gate capacitances as the registered DFT points (parallel plate, and the
Laplace single gate and GAA of fits/t1/cg/cg_table.json), measures the band's share. Read-only: it imports the frozen
scripts/t1_bb_widths.py without changing it and writes only fits/posthoc_2026-11-20/nfe_share.json.

    PYTHONPATH=. python scripts/posthoc_nfe_share.py
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scnfet import tob, w90 as W  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, [str(path)]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


def main():
    T = _load("t1", ROOT / "scripts" / "t1_bb_widths.py")
    kc = "r100_3"
    st = T.MB.geometry(kc)
    nv, Lx, nat = T.nval(st), st.cell[0, 0], T.n_host(st)
    dft = json.loads((ROOT / "fits/t1/dft.json").read_text())["calibration"]["100"]
    assert dft["key"] == kc and dft["variant"] == "HW3"
    W_A = dft["W_A"]
    cgt = json.loads((ROOT / "fits/t1/cg/cg_table.json").read_text())["widths"]["dft_100"]
    out = dict(key=kc, W_A=W_A, note=__doc__.split("\n\n")[0], variants={})
    for var in ("HW3", "HW2"):
        hr = W.read_hr(T.G.hr_path(T.MB, kc, var))
        E = T.G.w_eigs(hr, T.KG)
        Ef = T.G.w_eigs(hr, T.KF)
        cbm = Ef[:, nv].min()
        # conduction bands at Gamma within 0.6 eV of the CBM (the NFE state is the one HW2 lacks)
        g = int(np.argmin(np.abs(T.KF[:, 0])))
        cb_gamma = sorted(float(e - cbm) for e in Ef[g, nv:] if e - cbm < 0.6)
        ch = tob.Channel(E=np.sort(E, axis=1), a=Lx * 1e-10, n_val=nv)
        res = dict(gap=float(cbm - Ef[:, nv - 1].max()), n_bands=int(E.shape[1]), cb_gamma_rel_cbm_eV=cb_gamma)
        for geom in ("parallel_plate", "single", "gaa"):
            for eot in (1.0, 0.6):
                Cg = cgt[f"EOT {eot}"][geom]["C_F_per_m"]
                f = tob.ToBFET(ch, W=W_A * 1e-10, EOT=eot * 1e-9, carrier="n", n_atoms_per_cell=nat, C_G=Cg)
                m = f.metrics(span=T.G.TOB_SPAN, n_sweep=T.G.TOB_NSWEEP)
                res[f"n {geom} EOT {eot}"] = dict(I_on=float(m["I_on"]), v_inj=float(m["v_inj"]))
        out["variants"][var] = res
        print(f"{var}: {E.shape[1]} bands, gap {res['gap']:.4f} eV, CB at Gamma (rel. CBM, < 0.6 eV): "
              + ", ".join(f"{x:+.3f}" for x in cb_gamma), flush=True)
    chk = out["variants"]["HW3"]["n parallel_plate EOT 0.6"]["I_on"]
    print(f"check: HW3 parallel-plate EOT 0.6 I_ON {chk:.4f} vs stored dft.json {dft['HSE']['tob']['n EOT 0.6']['I_on']:.4f}")
    out["share"] = {}
    for k in out["variants"]["HW3"]:
        if not k.startswith("n "):
            continue
        a, b = out["variants"]["HW3"][k]["I_on"], out["variants"]["HW2"][k]["I_on"]
        out["share"][k] = dict(HW3=a, HW2=b, share=(a - b) / a)
        print(f"  {k:28s} HW3 {a:.3f}  HW2 {b:.3f} mA/um  ->  NFE share {100 * (a - b) / a:+.1f} %")
    d = ROOT / "fits/posthoc_2026-11-20"
    d.mkdir(parents=True, exist_ok=True)
    (d / "nfe_share.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
