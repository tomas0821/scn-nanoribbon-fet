"""Fitting the linear SK model: band and mass targets, multi-start least squares, and the Wannier projection.

Plan §3 D3: initial guess from a linear least-squares projection of the Wannier H(R) onto the SK forms
(`project_hr`), then multi-start refinement against bands with mass weights (`fit`). D4: error bands from the fit
covariance and from cross-validation; `FitResult.cov` and the near-optimal spread across starts feed those.

Residuals are dimensionless: (E_model − E_target)/σ_E for energies, and (m_target · c/(ħ²/m_e) − 1)/σ_rel for
masses, c being the band curvature. The inverse-mass form stays finite when a trial band goes flat.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from scnfet.sk import HB2M, CANONICAL_SIGN


@dataclass
class MassTarget:
    label: str
    k0: np.ndarray        # Cartesian 1/Å
    u: np.ndarray         # direction
    band: int             # sorted band index
    m: float              # target mass, m_e (negative for a hole band maximum)
    rel_sigma: float = 0.05


@dataclass
class CharTarget:
    """Band `band` at k must carry at least `wmin` of its weight on orbitals whose label ends in one of `orbitals`
    (e.g. (":dxy",)). Energies alone cannot tell a d_xy band from a d_z² band, and a fit without character
    information can land in a minimum with the characters swapped (it did: CBM at X became d_z²). The default σ is
    strict: at σ = 0.05 the swapped minimum still won on the Deng targets."""
    label: str
    k: np.ndarray
    band: int
    orbitals: tuple
    wmin: float = 0.8
    sigma: float = 0.01


@dataclass
class Targets:
    label: str
    k: np.ndarray                 # (nk, 3) unique Cartesian k-points
    ik: np.ndarray                # (n,) index into k
    band: np.ndarray              # (n,) sorted band index
    E: np.ndarray                 # (n,) target energy, eV
    sigma: np.ndarray             # (n,) eV
    masses: list = field(default_factory=list)
    ref: tuple = None             # (k (3,), band): energies are relative to this model eigenvalue (e.g. the VBM)
    group: np.ndarray = None      # (n,) str labels for reporting (e.g. "VB", "CB")
    floor: np.ndarray = None      # (n,) bool: one-sided target, model E ≥ E (penalised only when below)
    chars: list = field(default_factory=list)    # CharTarget list

    def __post_init__(self):
        if self.group is None:
            self.group = np.array(["all"] * len(self.E))
        if self.floor is None:
            self.floor = np.zeros(len(self.E), bool)


def model_values(model, x, T):
    """Model energies at the target (k, band) pairs (referenced like the targets) and model masses."""
    k = T.k if T.ref is None else np.vstack([T.k, T.ref[0][None]])
    E = model.eigvals(x, k)
    vals = E[T.ik, T.band]
    if T.ref is not None:
        vals = vals - E[-1, T.ref[1]]
    curv = np.array([model.curvature(x, mt.k0, mt.u, mt.band) for mt in T.masses])
    return vals, curv


def char_weights(model, x, T):
    """Weight of each CharTarget's band on its orbital set."""
    out = []
    for ct in T.chars:
        _, v = np.linalg.eigh(model.hk(x, np.asarray(ct.k)[None])[0])
        sel = [i for i, l in enumerate(model.orb_label) if l.endswith(tuple(ct.orbitals))]
        out.append(float(np.sum(np.abs(v[sel, ct.band]) ** 2)))
    return np.array(out)


def residuals(model, x, T):
    vals, curv = model_values(model, x, T)
    d = vals - T.E
    r_e = np.where(T.floor, np.minimum(d, 0.0), d) / T.sigma
    r_m = np.array([(mt.m * c / HB2M - 1.0) / mt.rel_sigma for mt, c in zip(T.masses, curv)])
    w = char_weights(model, x, T)
    r_c = np.array([max(ct.wmin - wi, 0.0) / ct.sigma for ct, wi in zip(T.chars, w)])
    return np.concatenate([r_e, r_m, r_c])


@dataclass
class FitResult:
    names: tuple
    x: np.ndarray
    free: np.ndarray              # bool mask of fitted parameters
    cost: float                   # ½ Σ r²
    n_res: int
    cov: np.ndarray               # (n_free, n_free) covariance of the free parameters, scaled by the reduced χ²
    starts: list                  # [(cost, x)] for every start, sorted by cost
    at_bound: list                # names of free parameters sitting on a bound
    near_best_spread: dict        # name → (min, max) over starts within 2 % of the best cost

    def params(self):
        return dict(zip(self.names, self.x))

    def sigma(self):
        s = np.full(len(self.names), np.nan)
        s[self.free] = np.sqrt(np.clip(np.diag(self.cov), 0, None))
        return dict(zip(self.names, s))


def sign_bounds(names, free_names, big=50.0):
    lo, hi = [], []
    for n in free_names:
        s = CANONICAL_SIGN.get(n, 0)
        lo.append(0.0 if s > 0 else -big)
        hi.append(0.0 if s < 0 else big)
    return np.array(lo), np.array(hi)


def fit(model, T, x0, fixed=(), bounds="canonical", n_starts=24, spread=0.4, seed=0, max_nfev=4000):
    """Multi-start least squares. `fixed` names keep their x0 values (e.g. "Ep" to fix the energy gauge when the
    targets are VBM-referenced). bounds: "canonical" (Harrison signs), None, or (lo, hi) over the free names."""
    names = model.names
    x0 = np.asarray(x0, float)
    free = np.array([n not in fixed for n in names])
    fnames = [n for n, f in zip(names, free) if f]
    if bounds == "canonical":
        lo, hi = sign_bounds(names, fnames)
    elif bounds is None:
        lo, hi = np.full(free.sum(), -np.inf), np.full(free.sum(), np.inf)
    else:
        lo, hi = map(np.asarray, bounds)
    eps = 1e-9
    def clip(v):
        return np.minimum(np.maximum(v, lo + eps * (np.isfinite(lo))), hi - eps * (np.isfinite(hi)))

    def full(xf):
        x = x0.copy()
        x[free] = xf
        return x

    fun = lambda xf: residuals(model, full(xf), T)
    rng = np.random.default_rng(seed)
    xf0 = clip(x0[free])
    starts = [xf0] + [clip(xf0 * (1 + spread * rng.standard_normal(xf0.size)) + 0.1 * rng.standard_normal(xf0.size))
                      for _ in range(n_starts - 1)]
    runs = []
    for s in starts:
        r = least_squares(fun, s, bounds=(lo, hi), x_scale="jac", max_nfev=max_nfev)
        runs.append((r.cost, r.x, r.jac))
    runs.sort(key=lambda t: t[0])
    cost, xb, J = runs[0]
    n_res = J.shape[0]
    dof = max(n_res - xb.size, 1)
    cov = (2 * cost / dof) * np.linalg.pinv(J.T @ J)
    tol = 1e-6 * (np.abs(lo) + np.abs(hi) + 1)
    at_bound = [n for n, v, a, b in zip(fnames, xb, lo, hi) if (np.isfinite(a) and v - a < 1e-6)
                or (np.isfinite(b) and b - v < 1e-6)]
    near = np.array([x for c, x, _ in runs if c <= 1.02 * cost + 1e-12])
    spread_d = {n: (near[:, i].min(), near[:, i].max()) for i, n in enumerate(fnames)}
    return FitResult(names, full(xb), free, cost, n_res, cov, [(c, full(x)) for c, x, _ in runs], at_bound,
                     spread_d)


def report(model, x, T):
    """RMS energy error per group (meV) and each mass against its target."""
    vals, curv = model_values(model, x, T)
    out = {"rms_meV": {}, "max_meV": {}, "masses": {},
           "floor_violations": int(np.sum(T.floor & (vals < T.E - 1e-6))),
           "floor_worst_meV": 1e3 * float(np.max(np.r_[0.0, (T.E - vals)[T.floor]]))}
    for g in np.unique(T.group[~T.floor]):
        m = (T.group == g) & ~T.floor
        err = vals[m] - T.E[m]
        out["rms_meV"][g] = 1e3 * float(np.sqrt(np.mean(err ** 2)))
        out["max_meV"][g] = 1e3 * float(np.max(np.abs(err)))
    for mt, c in zip(T.masses, curv):
        out["masses"][mt.label] = (float(HB2M / c), mt.m)
    out["chars"] = {ct.label: (w, ct.wmin) for ct, w in zip(T.chars, char_weights(model, x, T))}
    return out


# ---------------------------------------------------------------- Wannier → SK linear projection
def project_hr(model, hr, weights=None, rcond=None):
    """Least-squares SK parameters from a Wannier H(R) in the model's orbital order and lattice gauge.

    Only the (R, m, n) elements the model can represent enter the fit; the norm of the remaining elements
    (longer range, or symmetry-forbidden in the SK form) is returned as `out_of_model`, the projection's measure of
    what the SK form cannot carry. hr.H is divided by hr.deg, as in H(k).
    Returns (x, info) with info = dict(rms_in_model, out_of_model, n_elements).
    """
    rmap = {tuple(r): i for i, r in enumerate(hr.R)}
    Hw = hr.H.real / hr.deg[:, None, None]
    A_rows, b, w = [], [], []
    used = np.zeros(Hw.shape, bool)
    for r_i, R in enumerate(model.R):
        j = rmap.get(tuple(R))
        if j is None:
            continue
        mask = (np.abs(model.G[:, r_i]).sum(0) > 0) | (np.abs(model.C[r_i]) > 0)
        A_rows.append(model.G[:, r_i][:, mask].T)
        b.append(Hw[j][mask] - model.C[r_i][mask])
        w.append(np.ones(mask.sum()) if weights is None else weights[r_i][mask])
        used[j][mask] = True
    A, b, w = np.vstack(A_rows), np.concatenate(b), np.concatenate(w)
    x, *_ = np.linalg.lstsq(A * w[:, None], b * w, rcond=rcond)
    res = A @ x - b
    return x, dict(rms_in_model=float(np.sqrt(np.mean(res ** 2))),
                   out_of_model=float(np.sqrt(np.sum(Hw[~used] ** 2))), n_elements=int(b.size))


# ---------------------------------------------------------------- joint fits over several structures
def fit_joint(parts, names, x0, fixed=(), bounds=None, n_starts=24, spread=0.4, seed=0, max_nfev=4000):
    """Multi-start least squares over several (model, Targets) parts that share parameters by name (e.g. bulk and
    sheet; the sheet model carries extra on-site shift parameters). `names` is the global parameter list; each
    part's model.names must be a subset. bounds: (lo, hi) over the free names, or None. Returns a FitResult over
    `names` (cov from the stacked Jacobian)."""
    names = tuple(names)
    x0 = np.asarray(x0, float)
    free = np.array([n not in fixed for n in names])
    idx = [np.array([names.index(n) for n in m.names]) for m, _ in parts]
    if bounds is None:
        lo, hi = np.full(free.sum(), -np.inf), np.full(free.sum(), np.inf)
    else:
        lo, hi = map(np.asarray, bounds)
    eps = 1e-9

    def clip(v):
        return np.minimum(np.maximum(v, lo + eps * np.isfinite(lo)), hi - eps * np.isfinite(hi))

    def full(xf):
        x = x0.copy()
        x[free] = xf
        return x

    def fun(xf):
        x = full(xf)
        return np.concatenate([residuals(m, x[i], T) for (m, T), i in zip(parts, idx)])

    rng = np.random.default_rng(seed)
    xf0 = clip(x0[free])
    starts = [xf0] + [clip(xf0 * (1 + spread * rng.standard_normal(xf0.size)) + 0.1 * rng.standard_normal(xf0.size))
                      for _ in range(n_starts - 1)]
    runs = []
    for s in starts:
        r = least_squares(fun, s, bounds=(lo, hi), x_scale="jac", max_nfev=max_nfev)
        runs.append((r.cost, r.x, r.jac))
    runs.sort(key=lambda t: t[0])
    cost, xb, J = runs[0]
    fnames = [n for n, f in zip(names, free) if f]
    dof = max(J.shape[0] - xb.size, 1)
    cov = (2 * cost / dof) * np.linalg.pinv(J.T @ J)
    at_bound = [n for n, v, a, b in zip(fnames, xb, lo, hi) if (np.isfinite(a) and v - a < 1e-6)
                or (np.isfinite(b) and b - v < 1e-6)]
    near = np.array([x for c, x, _ in runs if c <= 1.02 * cost + 1e-12])
    spread_d = {n: (near[:, i].min(), near[:, i].max()) for i, n in enumerate(fnames)}
    return FitResult(names, full(xb), free, cost, J.shape[0], cov, [(c, full(x)) for c, x, _ in runs], at_bound,
                     spread_d)
