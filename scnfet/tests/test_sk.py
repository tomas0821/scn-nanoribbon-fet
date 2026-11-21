"""Exact checks of scnfet.sk and scnfet.skfit (no external data).

Run: ~/.venvs/scnfet/bin/python -m pytest scnfet/tests -q   (from the project root)
"""
import itertools

import numpy as np
import pytest

from scnfet import sk, skfit, w90

TOY = dict(Ep=-2.02, Et2g=2.71, Eeg=4.93, pd_s=-2.40, pd_p=1.38, dd1_s=-0.12, dd1_p=0.58, dd1_d=-0.22,
           pp_s=0.80, pp_p=-0.18, dd2_s=-0.03, dd2_p=0.05, dd2_d=-0.02)   # A_skmodel.md toy + small 3NN


def rand_rotation(rng):
    q, r = np.linalg.qr(rng.standard_normal((3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    return q if np.linalg.det(q) > 0 else -q


d_rep = sk.d_rep


@pytest.mark.parametrize("seed", range(5))
def test_sk_tables_rotation_covariant(seed):
    """h(QR) = Dᵀ h(R) D for every block; p transforms with D = Qᵀ. Catches any typo in the SK tables."""
    rng = np.random.default_rng(seed)
    Q = rand_rotation(rng)
    c = rng.standard_normal(3)
    c /= np.linalg.norm(c)
    Dp, Dd = Q.T, d_rep(Q)
    assert np.allclose(Dd.T @ Dd, np.eye(5), atol=1e-12)
    Vpp, Vpd, Vdd = rng.standard_normal(2), rng.standard_normal(2), rng.standard_normal(3)
    Vs, Ds = rng.standard_normal(1), np.eye(1)
    for (li, lj, V, Di, Dj) in [(1, 1, Vpp, Dp, Dp), (1, 2, Vpd, Dp, Dd), (2, 1, Vpd, Dd, Dp), (2, 2, Vdd, Dd, Dd),
                                (0, 0, Vs, Ds, Ds), (0, 1, Vs, Ds, Dp), (1, 0, Vs, Dp, Ds), (0, 2, Vs, Ds, Dd),
                                (2, 0, Vs, Dd, Ds)]:
        lhs = sk.sk_block(li, lj, Q @ c, V)
        rhs = Di.T @ sk.sk_block(li, lj, c, V) @ Dj
        assert np.allclose(lhs, rhs, atol=1e-12), (li, lj)


def test_parity_and_transpose():
    c = np.array([0.3, -0.5, 0.81])
    c /= np.linalg.norm(c)
    V = (-1.3, 0.7)
    assert np.allclose(sk.sk_block(1, 2, -c, V), -sk.sk_block(1, 2, c, V))           # p odd, d even
    assert np.allclose(sk.sk_block(2, 1, c, V), sk.sk_block(1, 2, -c, V).T)
    assert np.allclose(sk.sk_block(2, 2, c, (1, 2, 3)), sk.sk_block(2, 2, c, (1, 2, 3)).T)


@pytest.fixture(scope="module")
def bulk():
    m = sk.build(sk.rocksalt(4.50))
    return m, m.x_from(TOY)


def test_bulk_shells_and_hermiticity(bulk):
    m, x = bulk
    # directed bonds per primitive cell: Sc–N 6+6, Sc–Sc 12 (2NN) and 6 (3NN), N–N 12
    assert m.bonds == {"pd": 12, "dd1": 12, "dd2": 6, "pp": 12}
    k = np.random.default_rng(1).standard_normal((7, 3))
    Hk = m.hk(x, k)
    assert np.allclose(Hk, np.conj(np.transpose(Hk, (0, 2, 1))), atol=1e-12)


def test_bulk_degeneracies_and_cubic_symmetry(bulk):
    m, x = bulk
    P = sk.fcc_points(4.50)
    eG = m.eigvals(x, P["G"][None])[0]
    assert np.allclose(eG[:3], eG[0]) and np.allclose(eG[3:6], eG[3]) and np.allclose(eG[6:8], eG[6])
    assert eG[2] < eG[3]                                   # 3-fold p (Γ15) below 3-fold t2g (Γ25′)
    k = np.array([0.31, -0.17, 0.52])
    e0 = m.eigvals(x, k[None])[0]
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product([1, -1], repeat=3):
            e = m.eigvals(x, (np.array(signs) * k[list(perm)])[None])[0]
            assert np.allclose(e, e0, atol=1e-10)


def test_dxy_band_closed_form_and_ml_lock(bulk):
    """Along Γ–X_z, d_xy mixes with nothing, so its band is a sum of exact cosines.

    In-plane neighbours give constants: 4(¾σ1 + ¼δ1) from the 2NN and 4π2 from the 3NN. Out-of-plane 2NN give
    4(π1 + δ1) cos(q a/2); 3NN along z give 2δ2 cos(q a). Hence
    E(Γ25′) − E(X3) = 8(π1 + δ1) and m_l = ħ²/(a²[(π1 + δ1) − 2δ2]).
    """
    m, x = bulk
    a = 4.50
    q = np.linspace(0, 2 * np.pi / a, 41)
    k = np.column_stack([0 * q, 0 * q, q])
    Hk = m.hk(x, k)
    dxy = Hk[:, 0, 0].real                    # Sc0:dxy diagonal element
    assert np.allclose(np.abs(Hk[:, 0, 1:]), 0, atol=1e-12)    # decoupled from every other orbital
    p1, d1, s1, d2, p2 = TOY["dd1_p"], TOY["dd1_d"], TOY["dd1_s"], TOY["dd2_d"], TOY["dd2_p"]
    const = TOY["Et2g"] + 4 * (0.75 * s1 + 0.25 * d1) + 4 * p2
    ref = const + 4 * (p1 + d1) * np.cos(q * a / 2) + 2 * d2 * np.cos(q * a)
    assert np.allclose(dxy, ref, atol=1e-12)
    assert dxy[0] - dxy[-1] == pytest.approx(8 * (p1 + d1), abs=1e-12)
    X = sk.fcc_points(a)["X"]
    # the d_xy band is the CBM at X in this toy: band index 3
    ml = m.mass(x, X, [0, 0, 1], 3, h=1e-3)
    assert ml == pytest.approx(sk.HB2M / (a * a * ((p1 + d1) - 2 * d2)), rel=1e-5)


def test_hr_export_feeds_tob_1d():
    """A structure periodic along x only: LinearModel.hr() must give the same bands through w90.bands_1d."""
    a = 4.40                     # one Sc–N row of the two-plane sheet, along the conventional [100] axis
    st = sk.Structure(cell=np.diag([a, 30.0, 30.0]), species=["Sc", "N", "N", "Sc"],
                      pos=[[0, 0, 0], [a / 2, 0, 0], [0, 0, 2.18], [a / 2, 0, 2.18]], pbc=(True, False, False))
    m = sk.build(st)
    x = m.x_from(TOY)
    k, E = w90.bands_1d(m.hr(x), nk=64)
    E2 = m.eigvals(x, np.column_stack([2 * np.pi * k / a, 0 * k, 0 * k]))
    assert np.allclose(E, np.sort(E2, axis=1), atol=1e-10)


def test_projection_recovers_parameters_from_own_hr():
    """Round trip on a strained two-plane sheet: project_hr(model.hr(x)) must return x exactly."""
    a = 3.11
    st = sk.Structure(cell=np.array([[a, 0, 0], [0, a, 0], [0, 0, 25.0]]), species=["Sc", "N", "N", "Sc"],
                      pos=[[0, 0, 10], [a / 2, a / 2, 10], [0, 0, 12.18], [a / 2, a / 2, 12.18]],
                      pbc=(True, True, False))
    terms = {"dSc_surf": {0: np.eye(5), 3: np.eye(5)}, "dN_surf": {1: np.eye(3), 2: np.eye(3)}}
    m = sk.build(st, onsite_terms=terms)
    x = m.x_from({**TOY, "dSc_surf": -0.6, "dN_surf": 0.3})
    xp, info = skfit.project_hr(m, m.hr(x))
    # the Sc and N surface shifts are degenerate with Et2g/Eeg and Ep in a sheet where every atom is surface:
    # compare the Hamiltonians, and the bond integrals directly
    assert np.allclose(m.HR_blocks(xp), m.HR_blocks(x), atol=1e-10)
    bonds = [i for i, n in enumerate(m.names) if n in sk.CANONICAL_SIGN]     # the 10 bond integrals (any sign)
    assert np.allclose(xp[bonds], x[bonds], atol=1e-10)
    assert info["out_of_model"] < 1e-12


def test_fit_recovers_synthetic_targets(bulk):
    """Targets generated by the model itself, fitted from a perturbed start: cost → 0, parameters recovered."""
    m, x = bulk
    P = sk.fcc_points(4.50)
    kp, s, _ = sk.path(P, ["L", "G", "X", "U", None, "K", "G"], n=12)
    E = m.eigvals(x, kp)
    vbm = E[np.argmin(np.linalg.norm(kp, axis=1)), 2]
    ik, band = np.meshgrid(np.arange(len(kp)), np.arange(8), indexing="ij")
    T = skfit.Targets("synthetic", kp, ik.ravel(), band.ravel(), (E - vbm).ravel(), np.full(ik.size, 0.01),
                      masses=[skfit.MassTarget("ml", P["X"], np.array([0, 0, 1.0]), 3, m.mass(x, P["X"], [0, 0, 1], 3)),
                              skfit.MassTarget("mt", P["X"], np.array([1.0, 0, 0]), 3, m.mass(x, P["X"], [1, 0, 0], 3))],
                      ref=(P["G"], 2))
    rng = np.random.default_rng(3)
    x0 = x * (1 + 0.15 * rng.standard_normal(x.size))
    x0[m.names.index("Ep")] = x[m.names.index("Ep")]
    res = skfit.fit(m, T, x0, fixed=("Ep",), n_starts=4, seed=1)
    assert res.cost < 1e-8
    assert np.allclose(res.x, x, atol=1e-5)


def test_character_and_floor_targets(bulk):
    """CharTarget reads the orbital weight of a band; floor targets penalise only from below."""
    m, x = bulk
    P = sk.fcc_points(4.50)
    T = skfit.Targets("t", P["X"][None], np.array([0, 0]), np.array([3, 4]), np.array([0.0, 99.0]), np.ones(2),
                      floor=np.array([True, True]),
                      chars=[skfit.CharTarget("X3 = d_xy", P["X"], 3, (":dxy",), 0.9),
                             skfit.CharTarget("X3 = d_z2 (wrong)", P["X"], 3, (":dz2",), 0.9)])
    w = skfit.char_weights(m, x, T)
    assert w[0] > 0.999 and w[1] < 1e-6                 # the toy has the physical X3 = d_xy
    r = skfit.residuals(m, x, T)
    E = m.eigvals(x, P["X"][None])[0]
    assert r[0] == 0.0                                   # band 3 at X lies above the 0 eV floor
    assert r[1] == pytest.approx(E[4] - 99.0)            # band 4 lies below the 99 eV floor: negative residual
    assert r[2] == 0.0 and r[3] == pytest.approx(0.9 / 0.01)


def test_s_pseudosite_bulk():
    """Rocksalt with Sc 4s: 9 bands, Hermitian, cubic (Γ1 s singlet at Es + 12 ssσ when sd vanishes at Γ by symmetry),
    and with every s coupling at 0 the 8 Sc d + N p bands are unchanged."""
    a = 4.50
    st = sk.with_s(sk.rocksalt(a))
    assert st.species == ["Sc", "N", "Ss"]
    m = sk.build(st, shells=sk.rocksalt_shells(a) + sk.s_shells(a), names=sk.PARAM_NAMES + sk.S_PARAM_NAMES)
    assert m.norb == 9 and m.bonds["sp"] == 12 and m.bonds["ss1"] == 12 and m.bonds["sd1"] == 24
    p = dict(TOY, Es=7.0, sp_s=1.4, ss1_s=-0.3, sd1_s=-0.4)
    x = m.x_from(p)
    rng = np.random.default_rng(3)
    k = rng.standard_normal((6, 3))
    Hk = m.hk(x, k)
    assert np.allclose(Hk, np.conj(np.transpose(Hk, (0, 2, 1))))
    eG = m.eigvals(x, np.zeros((1, 3)))[0]
    assert np.isclose(eG, 7.0 + 12 * -0.3).sum() == 1        # Γ1: s couples to neither p (odd) nor d (Σ sd = 0) at Γ
    Q = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], float)     # cubic C3
    assert np.allclose(m.eigvals(x, k), m.eigvals(x, k @ Q.T))
    m8 = sk.build(sk.rocksalt(a))
    x0 = m.x_from(dict(TOY, Es=50.0))
    assert np.allclose(m.eigvals(x0, k)[:, :8], m8.eigvals(m8.x_from(TOY), k))


def test_three_centre_term_bulk():
    """Three-centre Sc–N–Sc term on rocksalt: each 2NN Sc pair shares exactly 2 N (12 pairs × 2 = 24 paths per Sc), the
    model stays Hermitian and cubic, and with the t3 parameters at 0 nothing changes."""
    a = 4.50
    tc = sk.three_centre_rocksalt(a)
    m = sk.build(sk.rocksalt(a), names=sk.PARAM_NAMES + sk.T3_PARAM_NAMES, three_centre=tc)
    assert m.bonds["t3"] == 24
    x = m.x_from(dict(TOY, t3_s=0.3, t3_pss=-0.2, t3_ppp=0.15, t3_psp=0.1))
    rng = np.random.default_rng(5)
    k = rng.standard_normal((6, 3))
    Hk = m.hk(x, k)
    assert np.allclose(Hk, np.conj(np.transpose(Hk, (0, 2, 1))))
    for Q in (np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], float), np.diag([-1.0, 1, 1]), np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1.0]])):
        assert np.allclose(m.eigvals(x, k), m.eigvals(x, k @ Q.T))
    m0 = sk.build(sk.rocksalt(a))
    assert np.allclose(m.eigvals(m.x_from(TOY), k), m0.eigvals(m0.x_from(TOY), k))


def test_fit_joint_recovers_shared_and_sheet_parameters():
    """Bulk + 2-plane sheet with a sheet-only on-site shift: synthetic targets from known parameters are recovered by
    the joint fit (shared parameters from both parts, the shift from the sheet alone)."""
    from scnfet import structures
    a = 4.50
    mb = sk.build(sk.rocksalt(a), shells=sk.rocksalt_shells(a))
    st = structures.slab(2, a_p=a / np.sqrt(2), d=a / 2)
    sc = [i for i, s in enumerate(st.species) if s == "Sc"]
    ms = sk.build(st, shells=sk.rocksalt_shells(a),
                  onsite_terms={"sh_z2": {i: np.diag([0, 0, 0, 0, 1.0]) for i in sc}})
    names = ms.names
    truth = dict(TOY, sh_z2=0.35)
    rng = np.random.default_rng(2)

    def targets(m, k):
        x = m.x_from(truth)
        E = m.eigvals(x, k)
        nb = E.shape[1]
        ik, band = np.repeat(np.arange(len(k)), nb), np.tile(np.arange(nb), len(k))
        return skfit.Targets("syn", k, ik, band, E[ik, band], np.full(ik.size, 0.01))

    kb = rng.standard_normal((12, 3))
    ks = np.c_[rng.standard_normal((12, 2)), np.zeros(12)]
    parts = [(mb, targets(mb, kb)), (ms, targets(ms, ks))]
    x0 = np.array([truth[n] for n in names]) * (1 + 0.05 * rng.standard_normal(len(names)))
    r = skfit.fit_joint(parts, names, x0, n_starts=2)
    got = r.params()
    assert abs(got["sh_z2"] - 0.35) < 1e-4 and all(abs(got[n] - TOY[n]) < 1e-4 for n in TOY)


@pytest.mark.parametrize("seed", range(3))
def test_crystal_field_frame_covariant(seed):
    """Bulk built in a rotated Cartesian frame (frame = Q) has the bands of the cubic-frame bulk at the rotated k.
    Covers the on-site t2g/e_g projectors together with the two- and three-centre terms."""
    Q = rand_rotation(np.random.default_rng(10 + seed))
    st0 = sk.rocksalt(4.50)
    st1 = sk.Structure(st0.cell @ Q.T, st0.species, st0.pos @ Q.T, frame=Q)
    t3 = sk.three_centre_rocksalt(4.50)
    names = sk.PARAM_NAMES + sk.T3_PARAM_NAMES
    m0 = sk.build(st0, shells=sk.rocksalt_shells(4.50), names=names, three_centre=t3)
    m1 = sk.build(st1, shells=sk.rocksalt_shells(4.50), names=names, three_centre=t3)
    p = dict(TOY, t3_s=0.33, t3_ppp=0.16, t3_psp=-0.48)
    k = np.random.default_rng(seed).standard_normal((6, 3))
    assert np.allclose(m0.eigvals(m0.x_from(p), k), m1.eigvals(m1.x_from(p), k @ Q.T), atol=1e-10)
    # a rotated frame without `frame` is a different (wrong) crystal field
    m2 = sk.build(sk.Structure(st1.cell, st1.species, st1.pos), shells=sk.rocksalt_shells(4.50), names=names,
                  three_centre=t3)
    assert not np.allclose(m0.eigvals(m0.x_from(p), k), m2.eigvals(m2.x_from(p), k @ Q.T), atol=1e-3)


def test_sheet_primitive_vs_conventional_cell():
    """The sheet in its primitive cell (bonds along [110], frame R45) folds onto the conventional cell (bonds along
    x, y): conventional Γ = primitive Γ ∪ primitive M̄."""
    from scnfet import structures as S
    a = 4.50
    kw = dict(shells=sk.rocksalt_shells(a), names=sk.PARAM_NAMES + sk.T3_PARAM_NAMES,
              three_centre=sk.three_centre_rocksalt(a))
    p = dict(TOY, t3_s=0.33, t3_ppp=0.16, t3_psp=-0.48)
    mp = sk.build(S.slab(2, a_p=a / np.sqrt(2), d=a / 2, cell="primitive"), **kw)
    mc = sk.build(S.slab(2, a_p=a / np.sqrt(2), d=a / 2, cell="conventional"), **kw)
    g = np.pi / (a / np.sqrt(2))
    ep = np.sort(np.r_[mp.eigvals(mp.x_from(p), np.zeros((1, 3)))[0],
                       mp.eigvals(mp.x_from(p), np.array([[g, g, 0]]))[0]])
    ec = mc.eigvals(mc.x_from(p), np.zeros((1, 3)))[0]
    assert np.allclose(ep, ec, atol=1e-10)
    P_t2g, _ = sk.crystal_field_projectors(S.R45)
    assert np.allclose(np.diag(P_t2g), [0, 1, 1, 1, 0], atol=1e-12)     # dx2-y2 (between the bonds) is t2g here
