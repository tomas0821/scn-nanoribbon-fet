"""The pre-registered joint bulk + sheet SK fit and gate (a) (fits/PREREG_sk_bulk_sheet.md, written 2026-11-11 before
any A2 output). Model: the plan's 13 parameters + the three-centre Sc–N–Sc term (t3_pss fixed at 0) + sheet on-site
shifts by C4v irreducible representations (sh_N_xy fixed at 0: with the sheet referenced to its own VBM, a uniform
shift of all sheet sites is a null direction). Weights: the "Tn" recipe.

Data:
  * bulk: as scripts/a1_sk_fit.py (dft/A1, hse_k16 / pbe_k16, 16³ grid);
  * sheet: dft/A2, W90 `hse_k16` / `pbe_k16` (16 WFs, 16×16 grid) and the DFT grid outputs.
The sheet Wannier functions are mapped to the model's sites by their centres (`.wout` final state) and to orbitals
by the W90 real-harmonic order (w90.W90_ORDER); `w90.regauge` then moves H(R) into the model's gauge.

Run from the project root: PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/a2_sk_fit_bulk_sheet.py [HSE06|PBE]
Writes fits/A2_joint/sk_joint_<label>.json and figures/A2_SK_joint_fit_<label>.png.
"""
import importlib.util
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

from scnfet import sk, skfit, w90
from scnfet.structures import R45

ROOT = Path(__file__).resolve().parents[1]
A2 = ROOT / "dft" / "A2"
OUT = ROOT / "fits" / "A2_joint"
A_BULK = 4.501
E_SPLIT, SIG_NEAR, SIG_FAR, MASS_SCALE = 4.0, 0.02, 0.5, 3        # the "Tn" recipe (pre-registered)
N_STARTS = 16
NV_SHEET = 6                                                    # N p valence bands in the 2 f.u. sheet
SHEET = {"HSE06": dict(seed="hse_k16", grid="hse_ibz_k16.out"), "PBE": dict(seed="pbe_k16", grid="pbe_nscf_k16.out")}
FIXED = ("Ep", "t3_pss", "sh_N_xy")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


F = _load("a1_sk_fit", ROOT / "scripts" / "a1_sk_fit.py")
T3 = _load("a1_sk_fit_t3", ROOT / "scripts" / "a1_sk_fit_t3.py")
M2 = _load("make_a2", A2 / "make_a2.py")

# C4v on-site patterns, model orbital order (sk.ORBITALS): Sc dxy, dyz, dzx, dx2-y2, dz2 | N px, py, pz
SHIFTS = {"sh_Sc_xy": ("Sc", np.diag([1.0, 0, 0, 0, 0])), "sh_Sc_xzyz": ("Sc", np.diag([0, 1.0, 1, 0, 0])),
          "sh_Sc_x2y2": ("Sc", np.diag([0, 0, 0, 1.0, 0])), "sh_Sc_z2": ("Sc", np.diag([0, 0, 0, 0, 1.0])),
          "sh_N_z": ("N", np.diag([0, 0, 1.0])), "sh_N_xy": ("N", np.diag([1.0, 1, 0]))}


def sheet_structure():
    """The A2 sheet exactly as computed (make_a2.py): atoms in the .win order Sc1, N1, N2, Sc2."""
    z = M2.zfrac()
    cell = np.diag([M2.A, M2.A, M2.C])
    frac = [(0, 0, z["Sc1"]), (0.5, 0.5, z["N1"]), (0, 0, z["N2"]), (0.5, 0.5, z["Sc2"])]
    # Sc–N bonds along (1, ±1, 0)/√2: the on-site t2g/e_g split is taken on the cubic axes (structures.R45)
    return sk.Structure(cell, ["Sc", "N", "N", "Sc"], np.array(frac) @ cell, pbc=(True, True, False), frame=R45)


def sheet_model():
    st = sheet_structure()
    terms = {n: {i: M for i, s in enumerate(st.species) if s == sp} for n, (sp, M) in SHIFTS.items()}
    return sk.build(st, shells=sk.rocksalt_shells(A_BULK), names=sk.PARAM_NAMES + sk.T3_PARAM_NAMES,
                    onsite_terms=terms, three_centre=sk.three_centre_rocksalt(A_BULK))


def wout_centres(path):
    txt = Path(path).read_text()
    blk = txt[txt.rfind("Final State"):]
    return np.array([[float(v) for v in m] for m in re.findall(
        r"WF centre and spread\s+\d+\s+\(\s*(-?[\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+)\s*\)", blk)])


def w90_sites(model, centres, cell):
    """Nearest model site of each WF centre, wrapping all three directions (W90 reports z modulo c, e.g. −10.14 Å for
    a WF at 12.36 Å). The assignment must be unambiguous (nearest < half the second nearest); the offsets are the
    WF-quality numbers of A2 criterion 4 and are only reported here. (A2: the Sc z² WFs sit 0.8–0.94 Å out into the
    vacuum, unambiguously on their Sc.)"""
    st = model.structure
    inv = np.linalg.inv(cell)
    site_of, off = [], []
    for c in centres:
        d = []
        for p in st.pos:
            f = (c - p) @ inv
            f -= np.round(f)
            d.append(np.linalg.norm(f @ cell))
        d = np.array(d)
        j = int(np.argmin(d))
        if d[j] >= 0.5 * np.sort(d)[1]:
            raise ValueError(f"WF centre {c} has no unambiguous site (distances {np.round(d, 3)})")
        site_of.append(j)
        off.append(float(d[j]))
    return site_of, off


def w90_perm(model, site_of):
    """Model orbital i ← W90 WF perm[i]; the WFs on one site keep W90's real-harmonic order for that l."""
    perm = []
    for i, lab in enumerate(model.orb_label):
        site, orb = int(lab.split(":")[0][2:] if lab.startswith("Sc") else lab.split(":")[0][1:]), lab.split(":")[1]
        wfs = [j for j, s in enumerate(site_of) if s == site]
        order = w90.W90_ORDER["d" if orb.startswith("d") else "p"]
        perm.append(wfs[order.index(orb.replace("dzx", "dxz"))])
    assert sorted(perm) == list(range(len(perm)))
    return perm


def qe_bands(path, alat):
    txt = Path(path).read_text()
    i = max(txt.rfind("End of self-consistent calculation"), txt.rfind("End of band structure calculation"))
    ks, es = [], []
    for m in re.finditer(r"k =\s*(-?[\d.]+)\s*(-?[\d.]+)\s*(-?[\d.]+) \(\s*\d+ PWs\)\s+bands \(ev\):\s+(.*?)(?=\n\s*\n\s*(?:occupation|k =|the Fermi|highest|Writing|$))", txt[i:], re.S):
        ks.append([float(m.group(j)) for j in (1, 2, 3)])
        es.append([float(v) for v in re.findall(r"-?\d+\.\d+", m.group(4))])
    return np.array(ks) * 2 * np.pi / alat, np.array(es)


class SheetCase:
    def __init__(self, label, model):
        c = SHEET[label]
        self.label, self.model = label, model
        d = A2 / "w90" / c["seed"]
        hr = w90.read_hr(d / f"{c['seed']}_hr.dat")
        cell = np.diag([M2.A, M2.A, M2.C])
        cen_w = wout_centres(d / f"{c['seed']}.wout")
        site_of, self.wf_offsets = w90_sites(model, cen_w, cell)
        perm = w90_perm(model, site_of)
        st = model.structure
        cen_m = st.pos[model.orb_site]
        # projection-only WFs live in W90's home cell at their projection site; the hr.dat R convention refers to
        # that site, not to the (Berry-phase) centre, so the gauge positions are the atoms
        self.hm = w90.regauge(hr, cell, st.pos[site_of], st.cell, cen_m, perm=perm)
        self.st = st
        self.froz = float(re.search(r"dis_froz_max = ([\d.]+)", (d / f"{c['seed']}.win").read_text()).group(1))
        kg, eg = qe_bands(A2 / "outputs" / c["grid"], M2.A)
        self.kg = kg
        eg = eg[:, 10:26]
        self.vbm_dft = eg[:, NV_SHEET - 1].max()
        self.froz -= self.vbm_dft
        self.eg = eg - self.vbm_dft
        kp, _, _ = self.path(200)
        ep = self.E_w90(kp)
        self.vbm = ep[:, NV_SHEET - 1].max()
        self.k_vbm = kp[np.argmax(ep[:, NV_SHEET - 1])]
        self.k_cbm = kp[np.argmin(ep[:, NV_SHEET])]
        self.masses = []
        for kind, k0, edge in (("e", self.k_cbm, NV_SHEET), ("h", self.k_vbm, NV_SHEET - 1)):
            e0 = self.E_w90(k0[None])[0]
            bands = [b for b in range(16) if abs(e0[b] - e0[edge]) < 1e-3]
            for u, ul in (([1, 0, 0], "100"), ([1, 1, 0], "110")):
                for b in bands:
                    self.masses.append((f"{kind}{b}_{ul}", k0, u, b, 0.03 / MASS_SCALE if kind == "e" else 0.05 / MASS_SCALE))
        self.mass_ref = {m[0]: self.mass_w90(*m[1:4]) for m in self.masses}

    def E_w90(self, k):
        return np.linalg.eigvalsh(w90.hk_frac(self.hm, self.st.kfrac(k)))

    def mass_w90(self, k0, u, band, h=0.004):
        u = np.asarray(u, float) / np.linalg.norm(u)
        E = self.E_w90(np.array([k0 - h * u, k0, k0 + h * u]))[:, band]
        return sk.HB2M / ((E[0] + E[2] - 2 * E[1]) / h ** 2)

    def path(self, n):
        g = np.pi / M2.A
        P = {"G": np.zeros(3), "X": np.array([g, 0, 0]), "M": np.array([g, g, 0])}
        return sk.path(P, ["G", "X", "M", "G"], n=n)

    def targets(self):
        kp, _, _ = self.path(30)
        ep = self.E_w90(kp) - self.vbm
        cbm = min(ep[:, NV_SHEET].min(), self.eg[:, NV_SHEET].min())
        ks, ik, band, E, sig, grp, flo = [], [], [], [], [], [], []
        for kk, ee, src in ((kp, ep, "path"), (self.kg, self.eg, "grid")):
            for i in range(len(kk)):
                n0 = len(ks)
                ks.append(kk[i])
                for b in range(16):
                    e = ee[i, b]
                    if src == "path" and e > E_SPLIT:
                        continue
                    near = (b < NV_SHEET and e >= -1.5) or (b >= NV_SHEET and e <= cbm + 2.5)
                    ik.append(n0)
                    band.append(b)
                    if src == "grid" and e > self.froz:
                        E.append(self.froz), sig.append(0.05), grp.append("floor"), flo.append(True)
                    else:
                        E.append(e), sig.append(SIG_NEAR if near else SIG_FAR)
                        grp.append(("VB" if b < NV_SHEET else "CB") + ("_near" if near else "_far")), flo.append(False)
        masses = [skfit.MassTarget(n, np.asarray(k0), np.asarray(u, float), b, self.mass_ref[n], rs)
                  for n, k0, u, b, rs in self.masses]
        return skfit.Targets(f"{self.label} sheet", np.array(ks), np.array(ik), np.array(band), np.array(E),
                             np.array(sig), masses=masses, ref=(self.k_vbm, NV_SHEET - 1), group=np.array(grp),
                             floor=np.array(flo))

    def score(self, model, x):
        kp, _, _ = self.path(80)
        out = {}
        e_vbm = model.eigvals(x, self.k_vbm[None])[0, NV_SHEET - 1]
        for nm, kk, ref in (("path", kp, self.E_w90(kp) - self.vbm), ("grid", self.kg, self.eg)):
            e = model.eigvals(x, kk) - e_vbm
            cbm = ref[:, NV_SHEET].min()
            m = np.zeros_like(ref, bool)
            m[:, :NV_SHEET] = ref[:, :NV_SHEET] >= -1.0
            m[:, NV_SHEET:] = ref[:, NV_SHEET:] <= cbm + 1.0
            err = (e - ref)[m]
            out[nm] = dict(rms_pm1eV_meV=1e3 * float(np.sqrt(np.mean(err ** 2))),
                           max_pm1eV_meV=1e3 * float(np.max(np.abs(err))))
        out["masses"] = {}
        for n, k0, u, b, _ in self.masses:
            mm = float(model.mass(x, k0, u, b))
            out["masses"][n] = dict(sk=mm, ref=float(self.mass_ref[n]), rel=float(mm / self.mass_ref[n] - 1))
        g = np.pi / M2.A
        eps = {}
        for lab, k in (("G", np.zeros(3)), ("X", np.array([g, 0, 0])), ("M", np.array([g, g, 0]))):
            eps[lab] = (model.eigvals(x, k[None])[0] - e_vbm, self.E_w90(k[None])[0] - self.vbm)
        out["levels"] = {f"CB_{lab}": dict(sk=float(s[NV_SHEET]), ref=float(r[NV_SHEET])) for lab, (s, r) in eps.items()}
        out["levels"]["Mbar_offset"] = dict(sk=float(eps["M"][0][NV_SHEET] - eps["G"][0][NV_SHEET]),
                                            ref=float(eps["M"][1][NV_SHEET] - eps["G"][1][NV_SHEET]))
        out["gate_a"] = dict(rms_ok=bool(max(out["path"]["rms_pm1eV_meV"], out["grid"]["rms_pm1eV_meV"]) <= 30),
                             masses_ok=bool(all(abs(v["rel"]) <= 0.10 for v in out["masses"].values())))
        return out


def seed_shifts(case_sheet, bulk_onsite):
    """Sheet on-site shifts that put the model's on-site diagonal (bulk crystal field on the cubic axes, frame R45,
    plus the shifts) on the Wannier on-site energies; both relative to N1 px, whose shift sh_N_xy is the fixed gauge.
    Frame-agnostic: the bulk part comes from the model itself. (Shift names are in the sheet-cell frame: there dxy
    points at the in-plane N, i.e. it is the cubic-frame e_g x²−y², and dx2-y2 is the cubic t2g xy.)"""
    ms, hm = case_sheet.model, case_sheet.hm
    lab = ms.orb_label
    H0w = hm.H[np.where((hm.R == 0).all(1))[0][0]].real
    H0m = ms.HR_blocks(ms.x_from({k: bulk_onsite[k] for k in ("Ep", "Et2g", "Eeg")}))[np.where((ms.R == 0).all(1))[0][0]]
    j = lab.index("N1:px")
    d = {l: (H0w[i, i] - H0w[j, j]) - (H0m[i, i] - H0m[j, j]) for i, l in enumerate(lab)}
    return {"sh_Sc_xy": d["Sc0:dxy"], "sh_Sc_xzyz": d["Sc0:dyz"], "sh_Sc_x2y2": d["Sc0:dx2-y2"],
            "sh_Sc_z2": d["Sc0:dz2"], "sh_N_z": d["N1:pz"], "sh_N_xy": 0.0}


def plot(rec, cb, cs, mb, ms):
    """Main-text Fig. 4 at 0.95 \\linewidth of sn-jnl (scripts/figstyle.py; re-rendered 2026-11-19 without the in-image
    title): (a) bulk, (b) sheet; grey HSE06 (or PBE) Wannier bands, orange SK, shading the ±1 eV windows scored."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import figstyle as FS

    FS.use(7.5)
    grid, ref, col = "#dcdad3", "#9b9a96", "#eb6834"
    fig, axes = plt.subplots(1, 2, figsize=(FS.width(0.95), 3.0), sharey=True,
                             gridspec_kw=dict(width_ratios=[1.35, 1]))
    names = ms.names
    x = np.array([rec["params"][n] for n in names])
    xb = x[[names.index(n) for n in mb.names]]
    lab = rec["label"]
    for i, (ax, (case, model, xx, nb, nv, title)) in enumerate(zip(axes, (
            (cb, mb, xb, 8, 3, "bulk ScN"), (cs, ms, x, 16, NV_SHEET, "M-ScN sheet")))):
        kp, sv, ticks = case.path(120)
        ew = case.E_w90(kp) - case.vbm
        if case is cs:
            e = model.eigvals(xx, kp) - model.eigvals(xx, cs.k_vbm[None])[0, NV_SHEET - 1]
        else:
            e = model.eigvals(xx, kp) - model.eigvals(xx, np.zeros((1, 3)))[0, 2]
        cbm = ew[:, nv].min()
        ax.axhspan(-1.0, 0.0, color=grid, alpha=0.6, lw=0, zorder=0)
        ax.axhspan(cbm, cbm + 1.0, color=grid, alpha=0.6, lw=0, zorder=0)
        for b in range(nb):
            ax.plot(sv, ew[:, b], color=ref, lw=2.6, alpha=0.55, zorder=2, solid_capstyle="round")
        for b in range(nb):
            ax.plot(sv, e[:, b], color=col, lw=0.9, zorder=3)
        for t, l in ticks:
            ax.axvline(t, color=grid, lw=0.6, zorder=0)
        ax.set_xticks([t for t, _ in ticks], ["Γ" if l == "G" else (l + "̄" if case is cs else l) for _, l in ticks])
        ax.set_xlim(ticks[0][0], ticks[-1][0])
        ax.set_ylim(-6.5, 8.0)
        FS.letter(ax, "ab"[i], title)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("$E - E_\\mathrm{VBM}$ (eV)")
    h = [Line2D([], [], color=ref, lw=2.6, alpha=0.55), Line2D([], [], color=col, lw=0.9),
         Patch(color=grid, alpha=0.6, lw=0)]
    fig.legend(h, [f"{lab} Wannier bands", "SK joint fit", "±1 eV fit-test windows"],
               loc="outside lower center", ncol=3)
    fig_path = ROOT / "figures" / f"A2_SK_joint_fit_{lab}.png"
    fig_path.parent.mkdir(exist_ok=True)
    FS.save(fig, fig_path)


def main():
    label = sys.argv[1] if len(sys.argv) > 1 else "HSE06"
    mb = T3.make_model(False)
    ms = sheet_model()
    names = ms.names
    assert set(mb.names) <= set(names)
    cb, cs = F.Case(label, mb), SheetCase(label, ms)
    if "--plot-only" in sys.argv:
        plot(json.loads((OUT / f"sk_joint_{label}.json").read_text()), cb, cs, mb, ms)
        return
    Tb, Ts = cb.targets(), cs.targets()
    for T in (Tb,):
        for mt in T.masses:
            mt.rel_sigma /= MASS_SCALE
        T.sigma = np.where(~T.floor & (T.sigma > SIG_NEAR), SIG_FAR, T.sigma)
    t3 = json.loads((ROOT / "fits" / "A1_bulk" / "sk_fit_t3.json").read_text())
    base = t3["Tn" if label == "HSE06" else "PT"]["params"]
    proj = json.loads((ROOT / "fits" / "A1_bulk" / "sk_projection.json").read_text())[label]["projected_params"]
    shifts = seed_shifts(cs, proj)
    free = [n for n in names if n not in FIXED]
    lo, hi = skfit.sign_bounds(names, free)
    x0s = {"bulk-fit + 0": np.array([{**base}.get(n, 0.0) for n in names]),
           "bulk-fit + W90 shifts": np.array([{**base, **shifts}.get(n, 0.0) for n in names])}
    for n in FIXED:
        for x in x0s.values():
            x[names.index(n)] = 0.0
    t = time.time()
    best = None
    for sname, x0 in x0s.items():
        r = skfit.fit_joint([(mb, Tb), (ms, Ts)], names, x0, fixed=FIXED, bounds=(lo, hi), n_starts=N_STARTS, seed=1)
        print(f"seed {sname}: cost {r.cost:.1f}", flush=True)
        if best is None or r.cost < best.cost:
            best, bseed = r, sname
    x = best.x
    xb = x[[names.index(n) for n in mb.names]]
    sb, ss = cb.score(mb, xb), cs.score(ms, x)
    gate = dict(bulk=sb["gate_a"], sheet=ss["gate_a"],
                passes=bool(all(sb["gate_a"].values()) and all(ss["gate_a"].values())))
    rec = dict(label=label, prereg="fits/PREREG_sk_bulk_sheet.md", seconds=round(time.time() - t, 1), cost=float(best.cost),
               best_seed=bseed, params={k: float(v) for k, v in best.params().items()},
               sigma={k: (None if not np.isfinite(v) else float(v)) for k, v in best.sigma().items()},
               at_bound=best.at_bound, bulk=sb, sheet=ss, gate_a=gate,
               sheet_mass_refs={k: float(v) for k, v in cs.mass_ref.items()}, sheet_froz_rel_vbm=cs.froz,
               seed_shifts={k: float(v) for k, v in shifts.items()}, sheet_wf_offsets_A=cs.wf_offsets,
               notes=["sheet on-site t2g/e_g split on the cubic axes (structures.R45); shift names in the sheet-cell "
                      "frame: sh_Sc_xy = cubic e_g x2-y2 (points at in-plane N), sh_Sc_x2y2 = cubic t2g xy",
                      "Sc z2 WFs 0.8 (HSE) / 0.94 (PBE) A off-site into the vacuum (A2 criterion 4); site mapping "
                      "by nearest atom, gauge positions = atoms"])
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"sk_joint_{label}.json").write_text(json.dumps(rec, indent=1, default=float))
    print(json.dumps(dict(gate=gate, bulk_rms=[sb["path"]["rms_pm1eV_meV"], sb["grid"]["rms_pm1eV_meV"]],
                          sheet_rms=[ss["path"]["rms_pm1eV_meV"], ss["grid"]["rms_pm1eV_meV"]],
                          bulk_masses={k: round(100 * v["rel"], 1) for k, v in sb["masses"].items()},
                          sheet_masses={k: round(100 * v["rel"], 1) for k, v in ss["masses"].items()},
                          sheet_levels=ss["levels"]), indent=1))
    plot(rec, cb, cs, mb, ms)


if __name__ == "__main__":
    main()
