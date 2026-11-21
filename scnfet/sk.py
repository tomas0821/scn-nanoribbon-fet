"""Slater–Koster two-centre tight binding for Sc d + N p models of ScN: bulk rocksalt, sheets and ribbons.

Model (plan §3 D3):
  * orthogonal basis, Sc d5 (xy, yz, zx, x²−y², 3z²−r²) + N p3 (x, y, z);
  * shells: NN Sc–N pd; 2NN Sc–Sc dd and N–N pp; 3NN Sc–Sc dd at d = a ("dd2"). The 3NN shell is what unlocks
    m_l: with 2NN dd only, the d_xy band along Γ–X_z is an exact cosine and E(Γ25′) − E(X3) = 8(ddπ + ddδ) fixes
    m_l = 8ħ²/(m a² [E(Γ25′) − E(X3)]). The Sc–N shell at a√3/2 is not in the model;
  * Harrison scaling V(d) = V(d_ref) (d_ref/d)^η per shell, with d_ref the bulk shell distance.

For a fixed geometry and fixed exponents, H(R) = C(R) + Σ_p x_p G_p(R) is linear in the parameters. `LinearModel`
stores the G_p, so a parameter update costs one tensor contraction, and projecting a Wannier H(R) onto the SK forms
(the fit's initial guess, plan D3) is an ordinary linear least-squares problem (`scnfet.skfit.project_hr`).

Units: Å and eV. k vectors are Cartesian in 1/Å unless a name says `frac`.
"""
from dataclasses import dataclass, field
import itertools

import numpy as np

from scnfet import w90

S3 = np.sqrt(3.0)
HB2M = 7.619964  # ħ²/m_e in eV Å²

ORBITALS = {"Sc": ("dxy", "dyz", "dzx", "dx2-y2", "dz2"), "N": ("px", "py", "pz"), "Ss": ("s",)}
ANGULAR = {"Sc": 2, "N": 1, "Ss": 0}
# "Ss" is the Sc 4s orbital as a pseudo-site at the Sc position (LR1 test, 2026-11-11; `with_s`). The builder handles
# one angular momentum per site, so a d + s atom is two co-located sites; their on-site s–e_g coupling is zero in
# cubic symmetry and same-position pairs are skipped.

PARAM_NAMES = ("Ep", "Et2g", "Eeg",
               "pd_s", "pd_p",
               "dd1_s", "dd1_p", "dd1_d",
               "pp_s", "pp_p",
               "dd2_s", "dd2_p", "dd2_d")

# Canonical signs (Harrison) for the bond integrals, used as optional fit bounds: +1 ≥ 0, −1 ≤ 0, 0 = free.
# Only pd and pp keep Harrison signs. The Sc–Sc dd integrals are effective and left free (T0 debug fits on Deng,
# 2026-11-11):
#   * unlocking m_l requires dd2_d > 0, since m_l = ħ²/(a²[(dd1_p + dd1_d) − 2 dd2_d]) and the 2NN-only value
#     (1.06–1.08) is lighter than HSE (1.43);
#   * at Γ and X the p–d coupling vanishes by symmetry, so the CB levels there are pure fcc dd tight binding.
#     Keeping X3 = d_xy as the CBM while the other d levels at X stay above ~5 eV needs dd1_d > 0 (fits: +0.27–0.29).
S_PARAM_NAMES = ("Es", "sp_s", "ss1_s", "sd1_s")    # Sc 4s: on-site, NN Sc s–N p σ, 2NN Sc s–Sc s σ, 2NN Sc s–Sc d σ

CANONICAL_SIGN = {"pd_s": -1, "pd_p": +1, "dd1_s": 0, "dd1_p": 0, "dd1_d": 0, "pp_s": +1, "pp_p": -1,
                  "dd2_s": 0, "dd2_p": 0, "dd2_d": 0}


# ---------------------------------------------------------------- Slater–Koster tables (Slater & Koster 1954)
def sk_pp(l, m, n, s, p):
    c = np.array([l, m, n])
    return np.outer(c, c) * (s - p) + np.eye(3) * p


def sk_pd(l, m, n, s, p):
    """p at the origin, d at R = d (l, m, n). Rows px, py, pz; columns xy, yz, zx, x²−y², 3z²−r²."""
    M = np.empty((3, 5))
    M[0, 0] = S3 * l * l * m * s + m * (1 - 2 * l * l) * p
    M[0, 1] = S3 * l * m * n * s - 2 * l * m * n * p
    M[0, 2] = S3 * l * l * n * s + n * (1 - 2 * l * l) * p
    M[0, 3] = S3 / 2 * l * (l * l - m * m) * s + l * (1 - l * l + m * m) * p
    M[0, 4] = l * (n * n - (l * l + m * m) / 2) * s - S3 * l * n * n * p
    M[1, 0] = S3 * m * m * l * s + l * (1 - 2 * m * m) * p
    M[1, 1] = S3 * m * m * n * s + n * (1 - 2 * m * m) * p
    M[1, 2] = S3 * l * m * n * s - 2 * l * m * n * p
    M[1, 3] = S3 / 2 * m * (l * l - m * m) * s - m * (1 + l * l - m * m) * p
    M[1, 4] = m * (n * n - (l * l + m * m) / 2) * s - S3 * m * n * n * p
    M[2, 0] = S3 * l * m * n * s - 2 * l * m * n * p
    M[2, 1] = S3 * n * n * m * s + m * (1 - 2 * n * n) * p
    M[2, 2] = S3 * n * n * l * s + l * (1 - 2 * n * n) * p
    M[2, 3] = S3 / 2 * n * (l * l - m * m) * s - n * (l * l - m * m) * p
    M[2, 4] = n * (n * n - (l * l + m * m) / 2) * s + S3 * n * (l * l + m * m) * p
    return M


def sk_dd(l, m, n, s, p, d):
    l2, m2, n2 = l * l, m * m, n * n
    M = np.empty((5, 5))
    M[0, 0] = 3 * l2 * m2 * s + (l2 + m2 - 4 * l2 * m2) * p + (n2 + l2 * m2) * d
    M[1, 1] = 3 * m2 * n2 * s + (m2 + n2 - 4 * m2 * n2) * p + (l2 + m2 * n2) * d
    M[2, 2] = 3 * n2 * l2 * s + (n2 + l2 - 4 * n2 * l2) * p + (m2 + n2 * l2) * d
    M[0, 1] = 3 * l * m2 * n * s + l * n * (1 - 4 * m2) * p + l * n * (m2 - 1) * d
    M[1, 2] = 3 * m * n2 * l * s + m * l * (1 - 4 * n2) * p + m * l * (n2 - 1) * d
    M[0, 2] = 3 * l2 * m * n * s + m * n * (1 - 4 * l2) * p + m * n * (l2 - 1) * d
    M[0, 3] = 1.5 * l * m * (l2 - m2) * s + 2 * l * m * (m2 - l2) * p + 0.5 * l * m * (l2 - m2) * d
    M[1, 3] = 1.5 * m * n * (l2 - m2) * s - m * n * (1 + 2 * (l2 - m2)) * p + m * n * (1 + (l2 - m2) / 2) * d
    M[2, 3] = 1.5 * n * l * (l2 - m2) * s + n * l * (1 - 2 * (l2 - m2)) * p - n * l * (1 - (l2 - m2) / 2) * d
    M[0, 4] = S3 * (l * m * (n2 - (l2 + m2) / 2) * s - 2 * l * m * n2 * p + 0.5 * l * m * (1 + n2) * d)
    M[1, 4] = S3 * (m * n * (n2 - (l2 + m2) / 2) * s + m * n * (l2 + m2 - n2) * p - 0.5 * m * n * (l2 + m2) * d)
    M[2, 4] = S3 * (l * n * (n2 - (l2 + m2) / 2) * s + l * n * (l2 + m2 - n2) * p - 0.5 * l * n * (l2 + m2) * d)
    M[3, 3] = 0.75 * (l2 - m2) ** 2 * s + (l2 + m2 - (l2 - m2) ** 2) * p + (n2 + (l2 - m2) ** 2 / 4) * d
    M[3, 4] = S3 * ((l2 - m2) * (n2 - (l2 + m2) / 2) / 2 * s + n2 * (m2 - l2) * p + (1 + n2) * (l2 - m2) / 4 * d)
    M[4, 4] = (n2 - (l2 + m2) / 2) ** 2 * s + 3 * n2 * (l2 + m2) * p + 0.75 * (l2 + m2) ** 2 * d
    for i in range(5):
        for j in range(i):
            M[i, j] = M[j, i]
    return M


def sk_sp(l, m, n, s):
    """s at the origin, p at R = d (l, m, n): one row, columns px, py, pz."""
    return np.array([[l, m, n]], float) * s


def sk_sd(l, m, n, s):
    """s at the origin, d at R: one row, columns xy, yz, zx, x²−y², 3z²−r²."""
    return np.array([[S3 * l * m, S3 * m * n, S3 * n * l, S3 / 2 * (l * l - m * m), n * n - (l * l + m * m) / 2]]) * s


def sk_block(li, lj, c, V):
    """Hopping block from an orbital of angular momentum li at the origin to lj at R, with R/|R| = c.

    V holds the bond integrals: (σ,) for ss/sp/sd, (σ, π) for pp/pd, (σ, π, δ) for dd. The reversed blocks use
    E_ji(c) = E_ij(−c)ᵀ.
    """
    l, m, n = c
    if (li, lj) == (0, 0):
        return np.array([[V[0]]], float)
    if (li, lj) == (0, 1):
        return sk_sp(l, m, n, *V)
    if (li, lj) == (1, 0):
        return sk_sp(-l, -m, -n, *V).T
    if (li, lj) == (0, 2):
        return sk_sd(l, m, n, *V)
    if (li, lj) == (2, 0):
        return sk_sd(-l, -m, -n, *V).T
    if (li, lj) == (1, 1):
        return sk_pp(l, m, n, *V)
    if (li, lj) == (1, 2):
        return sk_pd(l, m, n, *V)
    if (li, lj) == (2, 1):
        return sk_pd(-l, -m, -n, *V).T
    if (li, lj) == (2, 2):
        return sk_dd(l, m, n, *V)
    raise ValueError(f"no SK table for l = {li}, {lj}")


# ---------------------------------------------------------------- structures and shells
@dataclass
class Structure:
    cell: np.ndarray           # (3, 3) rows are lattice vectors, Å
    species: list              # e.g. ["Sc", "N", ...]
    pos: np.ndarray            # (nat, 3) Cartesian, Å
    pbc: tuple = (True, True, True)
    frame: np.ndarray = None   # (3, 3) columns: the rocksalt cubic axes in this structure's Cartesian frame (None = same)

    def __post_init__(self):
        self.cell = np.asarray(self.cell, float)
        self.pos = np.asarray(self.pos, float)
        if self.frame is not None:
            self.frame = np.asarray(self.frame, float)
            assert np.allclose(self.frame.T @ self.frame, np.eye(3), atol=1e-10), "frame must be orthonormal"

    def kfrac(self, kcart):
        """Cartesian k (1/Å) → fractional coordinates of the reciprocal basis, so that k·R = 2π k_frac·n."""
        return np.atleast_2d(kcart) @ self.cell.T / (2 * np.pi)


# real d orbitals (ORBITALS["Sc"] order) as quadratic forms r^T M r, orthonormal under 2·tr(A B)
D_FORMS = [np.array(m, float) for m in (
    [[0, .5, 0], [.5, 0, 0], [0, 0, 0]], [[0, 0, 0], [0, 0, .5], [0, .5, 0]], [[0, 0, .5], [0, 0, 0], [.5, 0, 0]],
    [[.5, 0, 0], [0, -.5, 0], [0, 0, 0]])] + [np.diag([-1, -1, 2]) / (2 * np.sqrt(3))]


def d_rep(Q):
    """D[γ, α] with d_α(Q s) = Σ_γ D[γ, α] d_γ(s) (real d orbitals, ORBITALS["Sc"] order)."""
    return np.array([[2 * np.trace(Mg @ Q.T @ Ma @ Q) for Ma in D_FORMS] for Mg in D_FORMS])


def crystal_field_projectors(frame=None):
    """(P_t2g, P_eg) in the structure frame. The octahedral t2g/e_g split is defined on the rocksalt cubic axes; for a
    structure whose Cartesian frame is rotated from them (e.g. the M-ScN primitive cell, Sc–N bonds along [110]),
    the cubic-frame orbitals are d_α(Qᵀ s) with Q = frame, so P = D Π Dᵀ with D = d_rep(Qᵀ)."""
    if frame is None:
        return np.diag([1.0, 1, 1, 0, 0]), np.diag([0, 0, 0, 1.0, 1])
    D = d_rep(np.asarray(frame, float).T)
    return D @ np.diag([1.0, 1, 1, 0, 0]) @ D.T, D @ np.diag([0, 0, 0, 1.0, 1]) @ D.T


def rocksalt(a=4.50):
    """Primitive fcc cell of rocksalt ScN with Sc at the origin."""
    cell = 0.5 * a * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], float)
    return Structure(cell, ["Sc", "N"], np.array([[0, 0, 0], [a / 2, 0, 0]]))


@dataclass(frozen=True)
class Shell:
    name: str
    pair: tuple       # sorted species pair
    d_ref: float      # reference (bulk) distance, Å
    d_min: float      # distances in [d_min, d_max) belong to this shell
    d_max: float
    eta: float        # Harrison exponent
    params: tuple     # bond-integral parameter names (σ, π[, δ])


DEFAULT_ETA = {"pd": 3.5, "dd1": 5.0, "dd2": 5.0, "pp": 2.0}   # Harrison 1980 values; to be fitted (plan D3)


def rocksalt_shells(a0=4.50, eta=None):
    """Shell windows from the bulk distances: NN a/2, 2NN a/√2, 3NN Sc–Sc a. Each window ends halfway to the next
    shell of that species pair, so strained sheets and ribbons (Sc–N 2.18–2.20 Å, Sc–Sc 3.10–3.11 Å) classify."""
    e = {**DEFAULT_ETA, **(eta or {})}
    nn, n2, n3, n4 = a0 / 2, a0 / np.sqrt(2), a0 * np.sqrt(3) / 2, a0
    return (Shell("pd", ("N", "Sc"), nn, 0.75 * nn, (nn + n3) / 2, e["pd"], ("pd_s", "pd_p")),
            Shell("dd1", ("Sc", "Sc"), n2, 0.75 * n2, (n2 + n4) / 2, e["dd1"], ("dd1_s", "dd1_p", "dd1_d")),
            Shell("dd2", ("Sc", "Sc"), n4, (n2 + n4) / 2, (n4 + a0 * np.sqrt(1.5)) / 2, e["dd2"],
                  ("dd2_s", "dd2_p", "dd2_d")),
            Shell("pp", ("N", "N"), n2, 0.75 * n2, (n2 + n4) / 2, e["pp"], ("pp_s", "pp_p")))


S_ETA = {"sp": 2.0, "ss1": 2.0, "sd1": 3.5}      # Harrison 1980 (s–p, s–s ∝ d⁻², s–d ∝ d⁻³·⁵)


def s_shells(a0=4.50, eta=None):
    """Shells for the Sc 4s pseudo-site "Ss": NN Ss–N (sp σ), 2NN Ss–Ss (ss σ) and 2NN Sc–Ss (sd σ); same windows as
    the pd and dd1 shells."""
    e = {**S_ETA, **(eta or {})}
    nn, n2, n3, n4 = a0 / 2, a0 / np.sqrt(2), a0 * np.sqrt(3) / 2, a0
    return (Shell("sp", ("N", "Ss"), nn, 0.75 * nn, (nn + n3) / 2, e["sp"], ("sp_s",)),
            Shell("ss1", ("Ss", "Ss"), n2, 0.75 * n2, (n2 + n4) / 2, e["ss1"], ("ss1_s",)),
            Shell("sd1", ("Sc", "Ss"), n2, 0.75 * n2, (n2 + n4) / 2, e["sd1"], ("sd1_s",)))


def with_s(st):
    """The structure with an "Ss" (Sc 4s) pseudo-site added at every Sc position (appended after the real sites)."""
    sc = [i for i, s in enumerate(st.species) if s == "Sc"]
    return Structure(st.cell.copy(), list(st.species) + ["Ss"] * len(sc), np.vstack([st.pos, st.pos[sc]]), st.pbc,
                     st.frame)


# ---------------------------------------------------------------- the linear model
@dataclass
class LinearModel:
    structure: Structure
    names: tuple                 # parameter names, the order of x
    R: np.ndarray                # (nR, 3) integer lattice vectors
    G: np.ndarray                # (n_par, nR, norb, norb) basis matrices
    C: np.ndarray                # (nR, norb, norb) constant part (fixed on-site blocks)
    orb_site: np.ndarray         # (norb,) site index of each orbital
    orb_label: list              # e.g. "Sc0:dxy"
    bonds: dict = field(default_factory=dict, repr=False)   # shell name → number of directed bonds (diagnostic)

    @property
    def norb(self):
        return self.C.shape[1]

    def x_from(self, params):
        """Parameter vector from a dict (missing names → 0)."""
        return np.array([params.get(n, 0.0) for n in self.names], float)

    def HR_blocks(self, x):
        return self.C + np.tensordot(np.asarray(x, float), self.G, axes=1)

    def hr(self, x):
        """The model as a Wannier90-style H(R) (lattice gauge, deg = 1)."""
        return w90.HR(R=self.R.copy(), deg=np.ones(len(self.R)), H=self.HR_blocks(x).astype(complex))

    def hk(self, x, kcart):
        kf = self.structure.kfrac(kcart)
        phase = np.exp(2j * np.pi * kf @ self.R.T)                         # (nk, nR)
        return np.einsum("kr,rmn->kmn", phase, self.HR_blocks(x))

    def eigvals(self, x, kcart):
        return np.linalg.eigvalsh(self.hk(x, kcart))

    def curvature(self, x, k0, u, band, h=0.004):
        """d²E/dk² (eV Å²) of `band` (sorted index) at k0 along the unit vector u, by central differences."""
        u = np.asarray(u, float) / np.linalg.norm(u)
        E = self.eigvals(x, np.array([k0 - h * u, k0, k0 + h * u]))[:, band]
        return (E[0] + E[2] - 2 * E[1]) / h ** 2

    def mass(self, x, k0, u, band, h=0.004):
        """Effective mass in m_e (signed: negative for a band maximum)."""
        return HB2M / self.curvature(x, k0, u, band, h)


def _translations(st, rmax):
    inv = np.linalg.inv(st.cell)                 # columns: reciprocal vectors / 2π
    ranges = []
    for i in range(3):
        if st.pbc[i]:
            spacing = 1.0 / np.linalg.norm(inv[:, i])
            ranges.append(range(-int(np.ceil(rmax / spacing)) - 1, int(np.ceil(rmax / spacing)) + 2))
        else:
            ranges.append([0])
    return np.array(list(itertools.product(*ranges)), int)


@dataclass(frozen=True)
class ThreeCentre:
    """Three-centre correction to the Sc–Sc 2NN d–d block through each N bonded to both Sc (LR1 test 2, 2026-11-11).

    The A1 Wannier 2NN Sc–Sc block breaks the two-centre form (x²−y² element −0.23 eV where two-centre theory gives
    +0.39). Model: for every Sc pair (i, j) in `pair_window` and every N site k within `via_window` of both,
        H_ij += (s_ik s_kj)^η-scaled Σ_p x_p B_p(i→k→j),
    with B built from SK blocks of a fictitious channel at k, so every term has the full point symmetry:
        t3_s    s channel:  sd(ĉ_ik) · sd(ĉ_kj)                      (d σ-lobes towards k, both sides)
        t3_pss  p channel, σ on both legs:   pd_σ(ĉ_ik) · pd_σ(ĉ_kj)
        t3_ppp  p channel, π on both legs
        t3_psp  p channel, σ on one leg and π on the other (symmetrised)
    Pairs at an edge share fewer N atoms, so the term is environment-dependent by construction. η scales each leg as
    (d_ref/d)^η, as the pd shell does."""
    pair: tuple = ("Sc", "Sc")
    via: str = "N"
    pair_window: tuple = (2.4, 3.84)
    via_window: tuple = (1.6, 3.07)
    d_ref: float = 2.25
    eta: float = 3.5
    params: tuple = ("t3_s", "t3_pss", "t3_ppp", "t3_psp")


T3_PARAM_NAMES = ThreeCentre().params


def three_centre_rocksalt(a0=4.50, eta=3.5):
    nn, n2, n3, n4 = a0 / 2, a0 / np.sqrt(2), a0 * np.sqrt(3) / 2, a0
    return ThreeCentre(pair_window=(0.75 * n2, (n2 + n4) / 2), via_window=(0.75 * nn, (nn + n3) / 2), d_ref=nn, eta=eta)


def _t3_blocks(c1, c2):
    """The four three-centre patterns for d at the origin → via site (direction c1) → d (direction c2 from via)."""
    s1, s2 = sk_block(2, 0, c1, (1.0,)), sk_block(0, 2, c2, (1.0,))
    a_s, a_p = sk_block(2, 1, c1, (1.0, 0.0)), sk_block(2, 1, c1, (0.0, 1.0))
    b_s, b_p = sk_block(1, 2, c2, (1.0, 0.0)), sk_block(1, 2, c2, (0.0, 1.0))
    return (s1 @ s2, a_s @ b_s, a_p @ b_p, a_s @ b_p + a_p @ b_s)


def build(st, shells=None, names=PARAM_NAMES, onsite_terms=None, onsite_const=None, three_centre=None):
    """Assemble the linear SK model for structure `st`.

    shells: tuple of `Shell` (default `rocksalt_shells(4.50)`).
    three_centre: an optional `ThreeCentre` (its params must be in `names`).
    onsite_terms: {param_name: {site: (n, n) matrix}}, extra *fitted* on-site parameters (surface or edge shifts,
        crystal fields), each a matrix pattern scaled by its parameter.
    onsite_const: {site: (n, n) matrix}, fixed on-site additions (e.g. edge blocks read from Wannier H(R)).
    """
    shells = shells or rocksalt_shells()
    onsite_terms = onsite_terms or {}
    names = tuple(names) + tuple(n for n in onsite_terms if n not in names)
    idx = {n: i for i, n in enumerate(names)}
    nat = len(st.species)
    nb = [len(ORBITALS[s]) for s in st.species]
    off = np.concatenate([[0], np.cumsum(nb)])
    norb = int(off[-1])
    rmax = max(sh.d_max for sh in shells)
    T = _translations(st, rmax)
    Tc = T @ st.cell

    blocks = {}          # R tuple -> index
    def ridx(n):
        key = tuple(int(v) for v in n)
        if key not in blocks:
            blocks[key] = len(blocks)
        return blocks[key]
    ridx((0, 0, 0))
    entries = []         # (param index, R index, i0, j0, matrix)
    nbond = {sh.name: 0 for sh in shells}
    by_pair = {}
    for sh in shells:
        by_pair.setdefault(sh.pair, []).append(sh)

    for i in range(nat):
        for j in range(nat):
            pair = tuple(sorted((st.species[i], st.species[j])))
            cand = by_pair.get(pair)
            if not cand:
                continue
            dvec = st.pos[j] + Tc - st.pos[i]
            dist = np.linalg.norm(dvec, axis=1)
            for t in np.where((dist > 1e-6) & (dist < rmax))[0]:
                d = dist[t]
                sh = next((s for s in cand if s.d_min <= d < s.d_max), None)
                if sh is None:
                    if d < min(s.d_min for s in cand):
                        raise ValueError(f"atoms {i},{j} only {d:.3f} Å apart: below every shell window")
                    continue
                c = dvec[t] / d
                scale = (sh.d_ref / d) ** sh.eta
                r = ridx(T[t])
                li, lj = ANGULAR[st.species[i]], ANGULAR[st.species[j]]
                for k, pname in enumerate(sh.params):
                    V = np.zeros(len(sh.params))
                    V[k] = scale
                    entries.append((idx[pname], r, off[i], off[j], sk_block(li, lj, c, V)))
                nbond[sh.name] += 1

    if three_centre is not None:
        tc = three_centre
        rmax3 = max(rmax, tc.pair_window[1], tc.via_window[1])
        T3 = _translations(st, 2 * rmax3)
        T3c = T3 @ st.cell
        nbond["t3"] = 0
        vias = [k for k, s in enumerate(st.species) if s == tc.via]
        for i in range(nat):
            for j in range(nat):
                if tuple(sorted((st.species[i], st.species[j]))) != tuple(sorted(tc.pair)):
                    continue
                dij = st.pos[j] + T3c - st.pos[i]
                lij = np.linalg.norm(dij, axis=1)
                for t in np.where((lij >= tc.pair_window[0]) & (lij < tc.pair_window[1]))[0]:
                    rj = st.pos[j] + T3c[t]
                    r = None
                    for k in vias:
                        rk = st.pos[k] + T3c - st.pos[i]                 # via sites relative to atom i
                        d1 = np.linalg.norm(rk, axis=1)
                        d2 = np.linalg.norm(st.pos[i] + rk - rj, axis=1)
                        for u in np.where((d1 >= tc.via_window[0]) & (d1 < tc.via_window[1])
                                          & (d2 >= tc.via_window[0]) & (d2 < tc.via_window[1]))[0]:
                            c1 = rk[u] / d1[u]
                            c2 = (rj - (st.pos[i] + rk[u])) / d2[u]
                            scale = (tc.d_ref / d1[u]) ** tc.eta * (tc.d_ref / d2[u]) ** tc.eta
                            r = ridx(T3[t]) if r is None else r
                            for pname, M in zip(tc.params, _t3_blocks(c1, c2)):
                                entries.append((idx[pname], r, off[i], off[j], scale * M))
                            nbond["t3"] += 1

    nR = len(blocks)
    G = np.zeros((len(names), nR, norb, norb))
    C = np.zeros((nR, norb, norb))
    for p, r, i0, j0, M in entries:
        G[p, r, i0:i0 + M.shape[0], j0:j0 + M.shape[1]] += M
    P_t2g, P_eg = crystal_field_projectors(st.frame)
    for i, s in enumerate(st.species):
        o = off[i]
        if s == "Sc":
            G[idx["Et2g"], 0, o:o + 5, o:o + 5] = P_t2g
            G[idx["Eeg"], 0, o:o + 5, o:o + 5] = P_eg
        elif s == "N":
            G[idx["Ep"], 0, o + np.arange(3), o + np.arange(3)] = 1.0
        elif s == "Ss":
            G[idx["Es"], 0, o, o] = 1.0
        else:
            raise ValueError(f"species {s!r} not in the Sc d + N p model")
    for pname, sites in onsite_terms.items():
        for site, M in sites.items():
            o = off[site]
            G[idx[pname], 0, o:o + nb[site], o:o + nb[site]] += np.asarray(M, float)
    for site, M in (onsite_const or {}).items():
        o = off[site]
        C[0, o:o + nb[site], o:o + nb[site]] += np.asarray(M, float)

    R = np.array(sorted(blocks, key=blocks.get), int)
    labels = [f"{s}{i}:{o}" for i, s in enumerate(st.species) for o in ORBITALS[s]]
    orb_site = np.repeat(np.arange(nat), nb)
    return LinearModel(st, names, R, G, C, orb_site, labels, nbond)


# ---------------------------------------------------------------- fcc k-points (Cartesian, 1/Å)
def fcc_points(a):
    """High-symmetry points of the fcc BZ in Cartesian 1/Å. X is X_z, so the Δ line Γ–X runs along z."""
    u = 2 * np.pi / a
    return {"G": np.zeros(3), "X": u * np.array([0, 0, 1.0]), "L": u * np.array([0.5, 0.5, 0.5]),
            "U": u * np.array([0.25, 0.25, 1.0]), "K": u * np.array([0.75, 0.75, 0.0]),
            "W": u * np.array([0.5, 0.0, 1.0])}


def path(points, labels, n=60):
    """Piecewise-linear path through labelled points. A label "U|K" style jump is given as two consecutive
    points with `None` between them. Returns (kcart (nk, 3), s (nk,) path length, ticks [(s, label)])."""
    ks, ss, ticks = [], [], []
    s0 = 0.0
    seg_start, start_in = None, False
    for lab in labels:
        if lab is None:
            seg_start = None
            continue
        p = points[lab]
        if seg_start is None:
            ticks.append((s0, lab))
            seg_start, start_in = p, False
            continue
        L = np.linalg.norm(p - seg_start)
        t = np.linspace(0, 1, n, endpoint=False)
        if start_in:
            t = t[1:]                      # the start point is the previous segment's end point
        start_in = True
        ks.append(seg_start + np.outer(t, p - seg_start))
        ss.append(s0 + t * L)
        ks.append(p[None])
        ss.append([s0 + L])
        s0 += L
        ticks.append((s0, lab))
        seg_start = p
    return np.vstack(ks), np.concatenate(ss), ticks
