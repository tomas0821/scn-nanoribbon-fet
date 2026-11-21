"""A4b-1 analysis (README.md). Two steps:

  windows   after job 1 (outputs/, pot_<r>.cube.gz copied back): vacuum level, first non-model state, the frozen and
            outer windows of W1 / W2 / W3, the band-count checks, and the variant .win files (w90/<r>/<v>/).
  criteria  after job 2 (the variant .wout and _hr.dat copied back): P1–P6 and the decision rule → a4b_criteria.json.
  hsurvey   A4b-2 staging: HSE nbnd from the A4b-1 NSCF grid (VBM + 10.23 + 0.5 eV).
  hwindows  after A4b-2 job 1 (outputs/hse_*, hproj_*, hpot_<r>.cube.gz, w90/<r>/H/<r>_H.win copied back): H1, the H3
            stop flag, the W1/W2 .win files (w90/<r>/H/<v>/) → a4b2_windows.json.
  hcriteria after A4b-2 job 2: H2 and the descriptive HSE results → a4b2_criteria.json.

    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P dft/A4b/parse_a4b.py windows [ribbon ...]
Writes dft/A4b/a4b_windows.json.

Rules (README, written before any A4b output):
  W1  frozen top = min(E_vac − 0.5 eV, minimum of the first band beyond the model count − 0.1 eV)   (A2's rule as run)
  W2  frozen top = min(W1, first non-model state − 0.1 eV); non-model: s_out = |ψ|² × (Löwdin share outside the
      projection set) + (1 − |ψ|²) ≥ 0.25, empty states on the 31-k path; also reported at 0.20 and 0.30
  W3  as W1, projection set S (B + Sc s on every Sc)
  outer = min(VBM + 10.23 eV, lowest energy of the top NSCF band − 0.05 eV); a cap below VBM + 10.23 is flagged
  E_vac: the bare + Hartree potential averaged along x, then over the (y, z) points whose distance to the nearest atom
  is ≥ 90 % of the largest such distance in the cell; the spread is reported.
"""
import gzip
import importlib.util
import json
import os
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, [str(path)]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


MB = _load("make_a4b", HERE / "make_a4b.py")
P = _load("parse_a2", ROOT / "dft" / "A2" / "parse_a2.py")
RY, BOHR = 13.605693122994, 0.529177210903
OUTER_PBE = 10.23            # cost test T2: PBE-twin outer window, eV above the VBM
VAC_MARGIN, NONMODEL_MARGIN = 0.5, 0.1
S_OUT = 0.25
CHECK_THR = (0.20, 0.30)
LNAME = "spdf"


def out(tag):
    return HERE / "outputs" / f"{tag}.out"


def vacuum_level(cube_gz):
    with gzip.open(cube_gz, "rt") as f:
        lines = f.read().splitlines()
    nat = abs(int(lines[2].split()[0]))
    n, vec = [], []
    for i in range(3):
        w = lines[3 + i].split()
        n.append(int(w[0]))
        vec.append([float(x) for x in w[1:4]])
    vec = np.array(vec) * BOHR                                  # step vectors, Å (orthorhombic cells here)
    atoms = np.array([[float(x) for x in lines[6 + j].split()[2:5]] for j in range(nat)]) * BOHR
    vals = np.array(" ".join(lines[6 + nat:]).split(), float).reshape(n) * RY      # x slowest, z fastest
    vyz = vals.mean(axis=0)                                     # average along the ribbon axis
    Ly, Lz = n[1] * vec[1, 1], n[2] * vec[2, 2]
    y = np.arange(n[1]) * vec[1, 1]
    z = np.arange(n[2]) * vec[2, 2]
    Y, Z = np.meshgrid(y, z, indexing="ij")
    dmin = np.full(Y.shape, np.inf)
    for a in atoms:
        dy = Y - a[1]
        dy -= Ly * np.round(dy / Ly)
        dz = Z - a[2]
        dz -= Lz * np.round(dz / Lz)
        dmin = np.minimum(dmin, np.hypot(dy, dz))
    far = dmin >= 0.9 * dmin.max()
    v = vyz[far]
    return dict(E_vac_eV=float(v.mean()), spread_meV=1e3 * float(v.max() - v.min()), n_points=int(far.sum()),
                max_dist_A=float(dmin.max()))


def projwfc_states(tag, path=None):
    """per path k: [(E, {(atom, species, l): weight}, |ψ|²)]"""
    txt = Path(path or out(tag)).read_text()
    stm = {int(s): (int(a) - 1, sp.strip(), LNAME[int(l)]) for s, a, sp, l in
           re.findall(r"state #\s*(\d+): atom\s+(\d+) \((\w+)\s*\), wfc\s+\d+ \(l=(\d)", txt)}
    res = []
    for kb in re.split(r"\n k =", txt.split("Lowdin Charges")[0])[1:]:
        st = []
        for m in re.finditer(r"==== e\(\s*\d+\) =\s*(-?[\d.]+) eV ====(.*?)\|psi\|\^2 = ([\d.]+)", kb, re.S):
            w = {}
            for c, s in re.findall(r"([\d.]+)\*\[#\s*(\d+)\]", m.group(2)):
                w[stm[int(s)]] = w.get(stm[int(s)], 0.0) + float(c)
            st.append((float(m.group(1)), w, float(m.group(3))))
        res.append(st)
    return res


def in_model_fn(key, st, seed):
    left, right = MB.edge_rows(key, st)
    edge = left | right

    def f(atom, sp, l):
        if (sp, l) in (("Sc", "d"), ("N", "p"), ("H", "s")):
            return True
        if sp == "Sc" and l == "s":
            return seed == "S" or (key in MB.EDGE_S and bool(edge[atom]))
        return False
    return f


def first_nonmodel(key, st, nocc, vbm, seed="B", path=None):
    inm = in_model_fn(key, st, seed)
    rows = []
    for ik, states in enumerate(projwfc_states(f"projwfc_{key}", path)):
        for b in range(nocc, len(states)):
            E, w, p2 = states[b]
            tot = sum(w.values())
            outside = sum(v for (a, sp, l), v in w.items() if not inm(a, sp, l)) / tot if tot > 0 else 1.0
            rows.append((E, p2 * outside + (1 - p2), ik, b + 1, p2))
    rows.sort()
    res = {}
    for thr in (S_OUT, *CHECK_THR):
        hit = next((r for r in rows if r[1] >= thr), None)
        res[f"{thr:.2f}"] = None if hit is None else dict(E_eV=hit[0], E_rel_vbm=hit[0] - vbm, s_out=round(hit[1], 3),
                                                          k_index=hit[2], band=hit[3], psi2=hit[4])
    return res


def ribbon_windows(key):
    info = json.loads((HERE / "a4b_inputs.json").read_text())[key]
    st = MB.geometry(key)
    nocc, ne = info["nocc"], info["nexcl"]
    kp, ep, _ = P.bands(out(f"bands_{key}"))
    kg, eg, _ = P.bands(out(f"nscf_{key}"))
    assert eg.shape == (info["nk_w90"], info["nbnd_nscf"]), f"{key}: NSCF grid {eg.shape}"
    vbm = float(max(ep[:, nocc - 1].max(), eg[:, nocc - 1].max()))
    cbm = float(min(ep[:, nocc].min(), eg[:, nocc].min()))
    vac = vacuum_level(HERE / f"pot_{key}.cube.gz")
    nm = {s: first_nonmodel(key, st, nocc, vbm, s) for s in ("B", "S")}
    r = dict(vbm=vbm, cbm=cbm, gap_path=float(ep[:, nocc].min() - ep[:, nocc - 1].max()), vacuum=vac,
             E_vac_rel_vbm=vac["E_vac_eV"] - vbm, path_top_band_min_rel_vbm=float(ep[:, -1].min() - vbm),
             nonmodel=nm)
    top_min = float(eg[:, -1].min())
    outer = min(vbm + OUTER_PBE, top_min - 0.05)
    r.update(outer=outer, outer_rel_vbm=outer - vbm, outer_capped=bool(outer < vbm + OUTER_PBE - 1e-9),
             nscf_top_band_min_rel_vbm=top_min - vbm)
    # band survey at this vacuum (for A4b-2's HSE nbnd): smallest N with grid band minimum ≥ VBM + 10.28 / + 11.28
    emin = eg.min(axis=0) - vbm
    r["survey_grid"] = {lab: (int(np.where(emin >= L)[0][0]) + 1 if (emin >= L).any() else None)
                        for lab, L in (("L", OUTER_PBE + 0.05), ("L+1.0", OUTER_PBE + 1.05))}
    var = {}
    for v, seed in (("W1", "B"), ("W2", "B"), ("W3", "S")):
        nw = info[f"num_wann_{seed}"]
        first_idx = float(eg[:, ne + nw].min()) if ne + nw < eg.shape[1] else np.inf    # first band beyond the model count
        froz = min(vac["E_vac_eV"] - VAC_MARGIN, first_idx - NONMODEL_MARGIN)
        binding = "vacuum" if froz == vac["E_vac_eV"] - VAC_MARGIN else "band count"
        if v == "W2":
            fn = nm["B"][f"{S_OUT:.2f}"]
            if fn is not None and fn["E_eV"] - NONMODEL_MARGIN < froz:
                froz, binding = fn["E_eV"] - NONMODEL_MARGIN, "non-model state"
        win = eg[:, ne:]
        n_froz = (win <= froz).sum(axis=1)
        n_outer = (win <= outer).sum(axis=1)
        var[v] = dict(seed=seed, num_wann=nw, froz=froz, froz_rel_vbm=froz - vbm, froz_rel_cbm=froz - cbm, binding=binding,
                      max_frozen_bands=int(n_froz.max()), min_outer_bands=int(n_outer.min()),
                      ok_counts=bool(n_froz.max() <= nw <= n_outer.min()),
                      vac_margin=vac["E_vac_eV"] - froz)
        d = MB.W90 / key / v
        d.mkdir(parents=True, exist_ok=True)
        base = (MB.W90 / key / f"{key}_{seed}.win").read_text()
        text = re.sub(r"dis_froz_max = .*", f"dis_froz_max = {froz:.4f}", base)
        text = re.sub(r"dis_win_max = .*", f"dis_win_max = {outer:.4f}", text)
        (d / f"{key}_{v}.win").write_text(text)
    r["variants"] = var
    # the path reference must hold every band below the highest frozen top (P3), with 0.1 eV to spare
    r["path_covers_frozen"] = bool(ep[:, -1].min() >= max(w["froz"] for w in var.values()) + 0.1)
    return r


def windows(keys):
    res = {}
    for key in keys:
        if not (out(f"nscf_{key}").exists() and (HERE / f"pot_{key}.cube.gz").exists()):
            print(f"{key}: outputs missing, skipped")
            continue
        r = res[key] = ribbon_windows(key)
        if not r["path_covers_frozen"]:
            print(f"   WARNING {key}: the band path stops below the highest frozen top + 0.1 eV (raise NB_BANDS)")
        print(f"=== {key}: VBM {r['vbm']:.3f}, gap {r['gap_path']:.3f} eV; E_vac − VBM {r['E_vac_rel_vbm']:.3f} eV "
              f"(spread {r['vacuum']['spread_meV']:.2f} meV, {r['vacuum']['n_points']} pts); outer VBM + "
              f"{r['outer_rel_vbm']:.2f}{' CAPPED' if r['outer_capped'] else ''}; survey N(L) {r['survey_grid']}")
        for thr, fn in r["nonmodel"]["B"].items():
            print(f"   first non-model (B, s_out ≥ {thr}): " + ("none" if fn is None else
                  f"VBM + {fn['E_rel_vbm']:.3f} (CBM + {fn['E_eV'] - r['cbm']:.3f}), band {fn['band']}, |ψ|² {fn['psi2']:.2f}"))
        for v, w in r["variants"].items():
            print(f"   {v} ({w['seed']}, {w['num_wann']} WFs): frozen VBM + {w['froz_rel_vbm']:.3f} = CBM + "
                  f"{w['froz_rel_cbm']:.3f} eV [{w['binding']}], vacuum margin {w['vac_margin']:.2f}; bands frozen ≤ "
                  f"{w['max_frozen_bands']}, in outer ≥ {w['min_outer_bands']} → {'ok' if w['ok_counts'] else 'FAIL'}")
    p = HERE / "a4b_windows.json"
    old = json.loads(p.read_text()) if p.exists() else {}
    old.update(res)
    p.write_text(json.dumps(old, indent=1, default=float))


# ---------------------------------------------------------------------------------------------------------- criteria
LIM = dict(device_max_meV=5.0, frozen_rms_meV=3.0, gap_meV=3.0, mass_rel=0.03, centre_ref_A=0.10)
DEVICE = (0.5, 0.5)                      # device window: [VBM − 0.5, CBM + 0.5] eV
HB2M = 3.80998212                        # ħ²/2m_e, eV Å²
MASS_WIN = 0.030                         # eV: path points within 30 meV of the extremum, at least 3


def wout_centres(path):
    txt = Path(path).read_text()
    blk = txt[txt.rfind("Final State"):]
    c = np.array([[float(v) for v in m] for m in re.findall(
        r"WF centre and spread\s+\d+\s+\(\s*(-?[\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+)\s*\)\s+([\d.]+)", blk)])
    om = float(re.search(r"Omega Total\s*=\s*([\d.]+)", blk).group(1))
    return c[:, :3], c[:, 3], om


def edge_mass(kx, E, band, kind):
    """Parabolic fit E = E0 + c (k − k0)² at the band extremum on the path (k0 = the extremal path point), over the points
    within 30 meV of it, at least the 3 nearest. Returns (m*/m_e signed by curvature, k0 index, n points)."""
    e = E[:, band]
    i0 = int(np.argmax(e) if kind == "v" else np.argmin(e))
    d = np.abs(e - e[i0])
    idx = np.where(d <= MASS_WIN)[0]
    if len(idx) < 3:
        idx = np.argsort(np.abs(kx - kx[i0]))[:3]
    dk2 = (kx[idx] - kx[i0]) ** 2
    A = np.vstack([np.ones_like(dk2), dk2]).T
    c = np.linalg.lstsq(A, e[idx], rcond=None)[0][1]
    return HB2M / c, i0, len(idx), idx


def variant_criteria(key, v, info, win, w90dir):
    from scnfet import w90 as W
    st = MB.geometry(key)
    nocc, ne = info["nocc"], info["nexcl"]
    var = win["variants"][v]
    nw = var["num_wann"]
    seed = f"{key}_{v}"
    kp, ep, _ = P.bands(out(f"bands_{key}"))
    Lx = st.cell[0, 0]
    kx = kp[:, 0]                                            # Å⁻¹ (P.bands: 2π/alat units × 2π/alat)
    hr = W.read_hr(w90dir / f"{seed}_hr.dat")
    assert hr.num_wann == nw
    Ew = np.linalg.eigvalsh(W.hk_1d(hr, kx * Lx / (2 * np.pi), axis=0))          # (nk, nw), sorted
    ref = ep[:, ne:ne + nw]                                  # DFT bands matched by index above the excluded ones
    nb = min(ref.shape[1], nw)
    ref, Ew = ref[:, :nb], Ew[:, :nb]
    vbm, cbm = win["vbm"], win["cbm"]
    err = np.abs(Ew - ref)
    dev = (ref >= vbm - DEVICE[0]) & (ref <= cbm + DEVICE[1])
    froz = ref <= var["froz"]
    iv, ic = nocc - 1 - ne, nocc - ne                         # VBM / CBM band within the window
    gap_dft = ref[:, ic].min() - ref[:, iv].max()
    gap_w = Ew[:, ic].min() - Ew[:, iv].max()
    masses = {}
    for lab, b, kind in (("VB", iv, "v"), ("CB", ic, "c")):
        m_d, i0, n, idx = edge_mass(kx, ref, b, kind)
        c_w = np.linalg.lstsq(np.vstack([np.ones(len(idx)), (kx[idx] - kx[i0]) ** 2]).T, Ew[idx, b], rcond=None)[0][1]
        m_w = HB2M / c_w
        masses[lab] = dict(m_dft=float(m_d), m_w90=float(m_w), rel=float(m_w / m_d - 1), k_index=i0, n_points=n)
    c, spr, om = wout_centres(w90dir / f"{seed}.wout")
    wfmap = json.loads((MB.W90 / key / "wf_map.json").read_text())[var["seed"]]
    assert len(wfmap) == len(c) == nw
    cell = np.diag(st.cell)
    offs = []
    for (r_, s_), m in zip(zip(c, spr), wfmap):
        d = r_ - st.pos[m["atom"]]
        d -= cell * np.round(d / cell)                       # W90 reports centres modulo the cell
        offs.append(dict(atom=m["atom"], species=m["species"], orbital=m["orbital"], offset_A=float(np.linalg.norm(d)),
                         d_A=d.round(4).tolist(), spread_A2=float(s_)))
    host = [o for o in offs if not (o["species"] == "Sc" and o["orbital"] in ("dz2", "s")) and o["species"] != "H"]
    res = dict(
        P1=dict(vac_margin_eV=var["vac_margin"], counts_ok=var["ok_counts"], outer_capped=win["outer_capped"],
                ok=bool(var["vac_margin"] >= VAC_MARGIN - 1e-9 and var["ok_counts"] and not win["outer_capped"])),
        P2=dict(max_meV=1e3 * float(err[dev].max()), n=int(dev.sum()), ok=bool(1e3 * err[dev].max() <= LIM["device_max_meV"])),
        P3=dict(rms_meV=1e3 * float(np.sqrt(np.mean(err[froz] ** 2))), max_meV=1e3 * float(err[froz].max()),
                n=int(froz.sum()), ok=bool(1e3 * np.sqrt(np.mean(err[froz] ** 2)) <= LIM["frozen_rms_meV"])),
        P4=dict(gap_dft=float(gap_dft), gap_w90=float(gap_w), d_meV=1e3 * float(gap_w - gap_dft),
                ok=bool(1e3 * abs(gap_w - gap_dft) <= LIM["gap_meV"])),
        P5=dict(masses=masses, ok=bool(all(abs(m["rel"]) <= LIM["mass_rel"] for m in masses.values()))),
        P6=dict(omega_total=om, host_nonz2_max_offset_A=max(o["offset_A"] for o in host),
                host_nonz2_over_ref=sum(o["offset_A"] > LIM["centre_ref_A"] for o in host),
                z2_offsets_A=[round(o["offset_A"], 3) for o in offs if o["orbital"] == "dz2"],
                z2_spreads_A2=[round(o["spread_A2"], 2) for o in offs if o["orbital"] == "dz2"],
                extra=[dict(atom=o["atom"], species=o["species"], orbital=o["orbital"], offset_A=round(o["offset_A"], 3),
                            spread_A2=round(o["spread_A2"], 2)) for o in offs if o["orbital"] == "s"]),
        wf=offs)
    res["pass_P1_P5"] = bool(all(res[k]["ok"] for k in ("P1", "P2", "P3", "P4", "P5")))
    return res


def criteria(keys):
    info_all = json.loads((HERE / "a4b_inputs.json").read_text())
    wins = json.loads((HERE / "a4b_windows.json").read_text())
    res = {}
    for key in keys:
        if key not in wins:
            print(f"{key}: no windows, skipped")
            continue
        r = res[key] = {}
        for v in ("W1", "W2", "W3"):
            d = MB.W90 / key / v
            if not ((d / f"{key}_{v}_hr.dat").exists() and (d / f"{key}_{v}.wout").exists()):
                print(f"{key} {v}: Wannier output missing")
                continue
            r[v] = variant_criteria(key, v, info_all[key], wins[key], d)
        # decision rule (README): W2 if it passes P1–P5, else W1, else back for review (LR6); W3 is reported only
        if r.get("W2", {}).get("pass_P1_P5"):
            r["decision"] = "W2"
        elif r.get("W1", {}).get("pass_P1_P5"):
            r["decision"] = "W1"
        elif "W1" in r and "W2" in r:
            r["decision"] = "neither: LR6, back for review" + (" (W3 alone passes: needs a decision)" if r.get("W3", {}).get("pass_P1_P5") else "")
        else:
            r["decision"] = "incomplete"
        print(f"=== {key}: decision {r['decision']}")
        for v in ("W1", "W2", "W3"):
            if v not in r:
                continue
            c = r[v]
            print(f"   {v}: P1 {c['P1']['ok']} | P2 device max {c['P2']['max_meV']:.2f} meV {c['P2']['ok']} | P3 frozen RMS "
                  f"{c['P3']['rms_meV']:.2f} (max {c['P3']['max_meV']:.1f}) {c['P3']['ok']} | P4 gap {c['P4']['d_meV']:+.2f} meV "
                  f"{c['P4']['ok']} | P5 m* VB {c['P5']['masses']['VB']['rel']:+.1%} CB {c['P5']['masses']['CB']['rel']:+.1%} "
                  f"{c['P5']['ok']} | P6 host max {c['P6']['host_nonz2_max_offset_A']:.3f} Å, z² {c['P6']['z2_offsets_A'][:4]}…, "
                  f"Ω {c['P6']['omega_total']:.1f}")
    p = HERE / "a4b_criteria.json"
    old = json.loads(p.read_text()) if p.exists() else {}
    old.update(res)
    p.write_text(json.dumps(old, indent=1, default=float))


# ------------------------------------------------------------------------------------------- A4b-2 (HSE06): H1–H3
# README "HSE stage (A4b-2): criteria registered now" (frozen 13 Nov):
#   H1  windows by the same rule on the HSE eigenvalues: E_vac from the HSE potential, non-model states from projwfc.x
#       on the HSE SCF (fallback, if that fails: band index from the PBE twin), outer = VBM + 11.05 eV, not capped
#   H2  Wannier90 reproduces the HSE grid eigenvalues inside the frozen window to ≤ 1 meV (A2 criterion 5)
#   H3  stop flag: any HSE non-model state inside the device window [VBM − 0.5, CBM + 0.5] → A4b stops for a decision
# The recipe is W2 (A4b-1 decision); W1 is run alongside and reported. Gap, edge masses and centres are descriptive
# (no HSE band path exists to test them against; the PBE twin did that).
# A4B_LOCAL (tests only): read outputs/, hpot cubes and w90/<r>/H/ from another directory.
BASE = Path(os.environ.get("A4B_LOCAL", HERE))
H2_MEV = 1.0
H_SEED = "B"
HSE_OUTER = MB.HSE_OUTER


def hout(tag):
    return BASE / "outputs" / f"{tag}.out"


def hdir(key, v=None):
    d = BASE / "w90" / key / "H"
    return d / v if v else d


def ribbon_hwindows(key):
    info = json.loads((HERE / "a4b_inputs.json").read_text())[key]
    st = MB.geometry(key)
    nocc, ne, nw = info["nocc"], info["nexcl"], info["num_wann_B"]
    txt = hout(f"hse_{key}").read_text()
    # the converged HSE run ends with QE's "!!" total energy; a max_seconds stop also prints JOB DONE, so check
    if "!!   total energy" not in txt or re.search(r"Maximum CPU time exceeded", txt):
        raise SystemExit(f"{key}: the HSE SCF did not converge (no final '!!' energy, or a max_seconds stop); not used")
    dexx = [float(x) for x in re.findall(r"est\. exchange err \(dexx\)\s*=\s*([-\d.E+]+)", txt)]
    kh, eh, _ = P.bands(hout(f"hse_{key}"))                   # last eigenvalue block: the final HSE one (IBZ)
    vbm, cbm = float(eh[:, nocc - 1].max()), float(eh[:, nocc].min())
    vac = vacuum_level(BASE / f"hpot_{key}.cube.gz")
    try:
        nm = first_nonmodel(key, st, nocc, vbm, H_SEED, hout(f"hproj_{key}"))
        nm_source = "projwfc on the HSE SCF"
    except Exception as e:                                    # H1 fallback, registered: band index from the PBE twin
        raise SystemExit(f"{key}: projwfc on the HSE SCF did not parse ({e}); apply the registered fallback by hand")
    r = dict(nbnd=int(eh.shape[1]), nk_ibz=int(eh.shape[0]), exx_outer_steps=len(dexx), dexx_final_Ry=dexx[-1],
             vbm=vbm, cbm=cbm, gap_grid=cbm - vbm, vacuum=vac,
             E_vac_rel_vbm=vac["E_vac_eV"] - vbm, nonmodel=nm, nonmodel_source=nm_source)
    top_min = float(eh[:, -1].min())
    outer = min(vbm + HSE_OUTER, top_min - 0.05)
    r.update(outer=outer, outer_rel_vbm=outer - vbm, outer_capped=bool(outer < vbm + HSE_OUTER - 1e-9),
             top_band_min_rel_vbm=top_min - vbm)
    fn = nm[f"{S_OUT:.2f}"]
    dev_top = cbm + DEVICE[1]
    r["H3"] = dict(device_top=dev_top, first_nonmodel_E=None if fn is None else fn["E_eV"],
                   first_nonmodel_rel_cbm=None if fn is None else fn["E_eV"] - cbm,
                   stop=bool(fn is not None and fn["E_eV"] <= dev_top),
                   stop_at_0p20=bool(nm["0.20"] is not None and nm["0.20"]["E_eV"] <= dev_top),
                   stop_at_0p30=bool(nm["0.30"] is not None and nm["0.30"]["E_eV"] <= dev_top))
    bases = {"B": (hdir(key) / f"{key}_H.win").read_text()}
    if "@KLIST@" in bases["B"]:
        raise SystemExit(f"{key}: w90/{key}/H/{key}_H.win still has @KLIST@; copy back the cluster copy (hw90pp fills it)")
    variants = [("W1", "B"), ("W2", "B")]
    if (hdir(key) / f"{key}_HS.win").exists():                 # amendment 2026-11-14: W3 (set S) when H3 fired
        bases["S"] = (hdir(key) / f"{key}_HS.win").read_text()
        variants.append(("W3", "S"))
    var = {}
    for v, seed in variants:
        nw = info[f"num_wann_{seed}"]
        base = bases[seed]
        first_idx = float(eh[:, ne + nw].min()) if ne + nw < eh.shape[1] else np.inf
        froz = min(vac["E_vac_eV"] - VAC_MARGIN, first_idx - NONMODEL_MARGIN)
        binding = "vacuum" if froz == vac["E_vac_eV"] - VAC_MARGIN else "band count"
        if v == "W2" and fn is not None and fn["E_eV"] - NONMODEL_MARGIN < froz:
            froz, binding = fn["E_eV"] - NONMODEL_MARGIN, "non-model state"
        win = eh[:, ne:]
        n_froz, n_outer = (win <= froz).sum(axis=1), (win <= outer).sum(axis=1)
        var[v] = dict(seed=seed, num_wann=nw, froz=froz, froz_rel_vbm=froz - vbm, froz_rel_cbm=froz - cbm, binding=binding,
                      device_inside_frozen=bool(froz >= cbm + DEVICE[1]),
                      max_frozen_bands=int(n_froz.max()), min_outer_bands=int(n_outer.min()),
                      ok_counts=bool(n_froz.max() <= nw <= n_outer.min()), vac_margin=vac["E_vac_eV"] - froz)
        var[v]["H1_ok"] = bool(var[v]["vac_margin"] >= VAC_MARGIN - 1e-9 and var[v]["ok_counts"] and not r["outer_capped"])
        d = hdir(key, v)
        d.mkdir(parents=True, exist_ok=True)
        text = re.sub(r"dis_froz_max = .*", f"dis_froz_max = {froz:.4f}", base)
        text = re.sub(r"dis_win_max = .*", f"dis_win_max = {outer:.4f}", text)
        (d / f"{key}_H{v}.win").write_text(text)
    r["variants"] = var
    return r


def hwindows(keys):
    pbe = json.loads((HERE / "a4b_windows.json").read_text())
    res = {}
    for key in keys:
        if not (hout(f"hse_{key}").exists() and (BASE / f"hpot_{key}.cube.gz").exists()):
            print(f"{key}: HSE outputs missing, skipped")
            continue
        r = res[key] = ribbon_hwindows(key)
        q = pbe.get(key, {})
        print(f"=== {key}: HSE gap (IBZ grid) {r['gap_grid']:.3f} eV (PBE twin {q.get('gap_path', float('nan')):.3f}); "
              f"E_vac − VBM {r['E_vac_rel_vbm']:.3f} (PBE {q.get('E_vac_rel_vbm', float('nan')):.3f}); outer VBM + "
              f"{r['outer_rel_vbm']:.2f}{' CAPPED' if r['outer_capped'] else ''} (top band min VBM + {r['top_band_min_rel_vbm']:.2f})")
        for thr, fn in r["nonmodel"].items():
            print(f"   first non-model (s_out ≥ {thr}): " + ("none" if fn is None else
                  f"CBM + {fn['E_eV'] - r['cbm']:.3f} eV, band {fn['band']}, |ψ|² {fn['psi2']:.2f}"))
        h3 = r["H3"]
        print(f"   H3 {'**STOP: a non-model state lies inside the device window**' if h3['stop'] else 'clear'} "
              f"(device top CBM + 0.5; at 0.20 / 0.30: {'stop' if h3['stop_at_0p20'] else 'clear'} / "
              f"{'stop' if h3['stop_at_0p30'] else 'clear'})")
        for v, w in r["variants"].items():
            print(f"   {v}: frozen CBM + {w['froz_rel_cbm']:.3f} eV [{w['binding']}], vacuum margin {w['vac_margin']:.2f}; "
                  f"bands frozen ≤ {w['max_frozen_bands']}, in outer ≥ {w['min_outer_bands']} (set {w['seed']}, num_wann "
                  f"{w['num_wann']}; device window inside frozen: {w['device_inside_frozen']}) → H1 {'ok' if w['H1_ok'] else 'FAIL'}")
    p = BASE / "a4b2_windows.json"
    old = json.loads(p.read_text()) if p.exists() else {}
    old.update(res)
    p.write_text(json.dumps(old, indent=1, default=float))


def hvariant(key, v, info, win):
    from scnfet import w90 as W
    st = MB.geometry(key)
    pset = win["variants"][v].get("seed", H_SEED)
    nocc, ne, nw = info["nocc"], info["nexcl"], info[f"num_wann_{pset}"]
    seed = f"{key}_H{v}"
    d = hdir(key, v)
    kh, eh, _ = P.bands(hout(f"hse_{key}"))
    Lx = st.cell[0, 0]
    hr = W.read_hr(d / f"{seed}_hr.dat")
    assert hr.num_wann == nw
    ref = eh[:, ne:ne + nw]
    Ew = np.linalg.eigvalsh(W.hk_1d(hr, kh[:, 0] * Lx / (2 * np.pi), axis=0))
    froz = ref <= win["variants"][v]["froz"]
    err = np.abs(Ew - ref)[froz]
    # descriptive: the Wannier band structure on a fine Γ–X̄ path (gap, edges, masses)
    kf = np.linspace(0, np.pi / Lx, 301)
    Ef = np.linalg.eigvalsh(W.hk_1d(hr, kf * Lx / (2 * np.pi), axis=0))
    iv, ic = nocc - 1 - ne, nocc - ne
    vb, cb = Ef[:, iv].max(), Ef[:, ic].min()
    masses = {}
    for lab, b, kind in (("VB", iv, "v"), ("CB", ic, "c")):
        m, i0, n, _ = edge_mass(kf, Ef, b, kind)
        masses[lab] = dict(m=float(m), k_frac=float(kf[i0] / kf[-1]), n_points=n)
    c, spr, om = wout_centres(d / f"{seed}.wout")
    wfmap = json.loads((MB.W90 / key / "wf_map.json").read_text())[pset]
    cell = np.diag(st.cell)
    offs = []
    for (r_, s_), m in zip(zip(c, spr), wfmap):
        dd = r_ - st.pos[m["atom"]]
        dd -= cell * np.round(dd / cell)
        offs.append(dict(atom=m["atom"], species=m["species"], orbital=m["orbital"], offset_A=float(np.linalg.norm(dd)),
                         spread_A2=float(s_)))
    z2 = [o for o in offs if o["orbital"] == "dz2"]
    return dict(H2=dict(max_meV=1e3 * float(err.max()), n=int(err.size), ok=bool(1e3 * err.max() <= H2_MEV)),
                gap_w90_path=float(cb - vb), vbm_k_frac=float(kf[np.argmax(Ef[:, iv])] / kf[-1]),
                cbm_k_frac=float(kf[np.argmin(Ef[:, ic])] / kf[-1]), masses=masses, omega_total=om,
                z2_offset_A=[min(o["offset_A"] for o in z2), max(o["offset_A"] for o in z2)],
                z2_spread_A2=[min(o["spread_A2"] for o in z2), max(o["spread_A2"] for o in z2)],
                extra=[dict(atom=o["atom"], species=o["species"], offset_A=round(o["offset_A"], 3),
                            spread_A2=round(o["spread_A2"], 2)) for o in offs if o["orbital"] == "s"])


H4_MEV = 5.0
# Dated decisions that replace the rule's pick (README "Amendment 2026-11-14 (2)").
DECISION_OVERRIDE = {
    "r110_4": ("W1 device Hamiltonian and SK blocks (decision 2026-11-14)",
               "W3 passes H1/H2 on the grid but differs from W1 by up to 222 meV between grid points (interior Sc s "
               "WFs delocalised, 168–195 Å²; disentanglement at the cap); W2's frozen top lies below the edge CB."),
}


def h4_pairs(key, info, win, variants, nk=601):
    """H4 (amendment 2): every pair of variants agrees within 5 meV over the part of the device window both freeze,
    on a 601-point Γ–X̄ path (bands index-matched above the excluded ones)."""
    from scnfet import w90 as W
    st = MB.geometry(key)
    t = np.linspace(0, 0.5, nk)
    E = {}
    for v in variants:
        hr = W.read_hr(hdir(key, v) / f"{key}_H{v}_hr.dat")
        E[v] = np.linalg.eigvalsh(W.hk_1d(hr, t, axis=0))
    vbm, cbm = win["vbm"], win["cbm"]
    res = {}
    for i, a in enumerate(variants):
        for b in variants[i + 1:]:
            n = min(E[a].shape[1], E[b].shape[1])
            top = min(cbm + DEVICE[1], win["variants"][a]["froz"], win["variants"][b]["froz"])
            m = (E[a][:, :n] >= vbm - DEVICE[0]) & (E[a][:, :n] <= top)
            d = 1e3 * float(np.abs(E[a][:, :n] - E[b][:, :n])[m].max()) if m.any() else float("nan")
            res[f"{a}-{b}"] = dict(max_meV=d, top_rel_cbm=top - cbm, ok=bool(d <= H4_MEV))
    return res


def hcriteria(keys):
    info_all = json.loads((HERE / "a4b_inputs.json").read_text())
    wins = json.loads((BASE / "a4b2_windows.json").read_text())
    res = {}
    for key in keys:
        if key not in wins:
            print(f"{key}: no HSE windows, skipped")
            continue
        r = res[key] = {}
        w, h3 = wins[key], wins[key]["H3"]["stop"]
        for v in w["variants"]:
            d = hdir(key, v)
            if (d / f"{key}_H{v}_hr.dat").exists() and (d / f"{key}_H{v}.wout").exists():
                r[v] = hvariant(key, v, info_all[key], wins[key])
            else:
                print(f"{key} {v}: Wannier output missing")
        ok2 = "W2" in r and w["variants"]["W2"]["H1_ok"] and r["W2"]["H2"]["ok"]
        if h3:      # amendment 2026-11-14: W3 for the device Hamiltonian if it passes H1, H2 and covers the device window
            w3 = w["variants"].get("W3")
            if w3 is None or "W3" not in r:
                r["decision"] = "STOP (H3): W3 pending"
            elif w3["H1_ok"] and r["W3"]["H2"]["ok"] and w3["device_inside_frozen"]:
                r["decision"] = "W3 device Hamiltonian (amendment 2026-11-14); W2 for the SK blocks" + ("" if ok2 else " (W2 fails H1/H2!)")
            else:
                r["decision"] = "W3 fails: W2 + the band reported as a limitation (amendment fallback)"
        else:
            r["decision"] = ("W2" if ok2 else "incomplete" if "W2" not in r else "W2 fails H1/H2: back for review")
        r["H4"] = h4_pairs(key, info_all[key], w, [v for v in ("W1", "W2", "W3") if v in r])
        if key in DECISION_OVERRIDE:
            r["decision_by_rule"] = r["decision"]
            r["decision"], r["decision_note"] = DECISION_OVERRIDE[key]
        print(f"=== {key}: {r['decision']}" + (f"  [by rule: {r['decision_by_rule']}]" if "decision_by_rule" in r else ""))
        for pair, c in r["H4"].items():
            print(f"   H4 {pair}: max |ΔE| {c['max_meV']:.1f} meV up to CBM + {c['top_rel_cbm']:.2f} eV → {'ok' if c['ok'] else 'FAIL'}")
        for v in ("W1", "W2", "W3"):
            if v in r:
                c = r[v]
                print(f"   {v}: H1 {w['variants'][v]['H1_ok']} | H2 grid max {c['H2']['max_meV']:.3f} meV ({c['H2']['n']} values) "
                      f"{c['H2']['ok']} | gap {c['gap_w90_path']:.3f} eV (VBM k {c['vbm_k_frac']:.2f}, CBM k {c['cbm_k_frac']:.2f} "
                      f"of Γ–X̄) | m* VB {c['masses']['VB']['m']:.3f} CB {c['masses']['CB']['m']:.3f} | z² "
                      f"{c['z2_offset_A'][0]:.2f}–{c['z2_offset_A'][1]:.2f} Å | Ω {c['omega_total']:.1f}")
    p = BASE / "a4b2_criteria.json"
    old = json.loads(p.read_text()) if p.exists() else {}
    old.update(res)
    p.write_text(json.dumps(old, indent=1, default=float))


def hsurvey(keys, margin=0.5):
    """A4b-2 nbnd from the A4b-1 PBE NSCF grid: smallest N whose grid minimum is ≥ VBM + 10.23 + margin eV."""
    wins = json.loads((HERE / "a4b_windows.json").read_text())
    for key in keys:
        kg, eg, _ = P.bands(out(f"nscf_{key}"))
        emin = eg.min(axis=0) - wins[key]["vbm"]
        idx = np.where(emin >= OUTER_PBE + margin)[0]
        print(f"{key}: N = {idx[0] + 1 if len(idx) else 'beyond ' + str(len(emin))} (make_a4b.NB_HSE {MB.NB_HSE[key]})")


if __name__ == "__main__":
    modes = dict(windows=windows, criteria=criteria, hwindows=hwindows, hcriteria=hcriteria, hsurvey=hsurvey)
    if len(sys.argv) > 1 and sys.argv[1] in modes:
        modes[sys.argv[1]](sys.argv[2:] or list(MB.RIBBONS))
    else:
        print(__doc__)
