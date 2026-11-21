import numpy as np

from scnfet import frohlich as F
from scnfet import guard as G


def test_3d_rate_reproduces_mu():
    """Mu et al. PRB 104, 075118 (2021): ab initio τ 104–114 fs below ħω_LO with their LDA inputs."""
    f = F.Frohlich(0.39, eps_inf=9.8, eps_0=9.8 * (75.8 / 33.4) ** 2, hw=0.0758)
    assert 104 <= 1e15 / f.rate_3d(0.005) <= 118


def test_2d_rate_against_brute_force():
    """Eq. 14 by the closed angular form against a Cartesian q-grid with a Gaussian-broadened delta."""
    f = F.Frohlich(0.2)
    E = 0.03
    k = f.k(E)
    q = np.arange(-0.35, 0.35, 0.0004)
    qx, qy = np.meshgrid(q, q, indexing="ij")
    Ef = F.H2M / f.m * ((k + qx) ** 2 + qy ** 2)
    g2 = f.K_A / (f.eps12 + f.r * np.hypot(qx, qy)) ** 2
    sig = 0.002
    tot = 0.0
    for D, N in ((f.hw, f.N0), (-f.hw, f.N0 + 1)):
        if E + D > 0:
            d = np.exp(-0.5 * ((Ef - E - D) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
            tot += N * (g2 * d).sum() * 0.0004 ** 2
    brute = 2 * np.pi / F.HBAR * tot / (2 * np.pi) ** 2
    assert abs(f.rate_2d(E) / brute - 1) < 0.03


def test_2d_large_q_limit_is_bulk():
    """Sohier: for r_eff q >> eps12 the 2D coupling equals the bulk one (V = A t)."""
    f = F.Frohlich(0.3)
    bulk = 4 * np.pi * F.E2 * f.hw / 2 / f.t * (1 / f.eps_inf - 1 / f.eps_0)   # |g_3D|² A q²
    assert abs(f.K_A / f.r ** 2 / bulk - 1) < 1e-12


def test_ballisticity_bounds():
    assert F.ballisticity(15.0) == 0.5 and 0 < F.ballisticity(1.0) < F.ballisticity(100.0) < 1


def test_square_barrier_transmission():
    m, V0, a, E = 0.3, 0.3, 30.0, np.array([0.05, 0.1, 0.2])
    U = np.full(300, V0)
    T = G.t_qm(E, U, a / 300, m, 0.0)
    kap = np.sqrt((V0 - E) * m / G.H2M)
    exact = 1 / (1 + V0 ** 2 * np.sinh(kap * a) ** 2 / (4 * E * (V0 - E)))
    assert np.allclose(T, exact, rtol=1e-6)
    assert np.allclose(G.t_qm(E, np.zeros(50), 1.0, m, 0.0), 1.0)


def test_sdt_limits():
    heavy = G.sdt(5.0, 2.0, 0.30, 0.05)[0]
    light = G.sdt(0.1, 2.0, 0.30, 0.05)[0]
    wider = G.sdt(0.1, 3.0, 0.30, 0.05)[0]
    assert abs(heavy - 1) < 0.05 and light > 1.05 and wider > light


def test_gidl_closed_window():
    I, w = G.gidl_bound(0.2, 1.5, 2.0, 0.30, 0.05)
    assert w < 0 and I < 1e-30
