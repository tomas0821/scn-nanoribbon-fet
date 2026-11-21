"""scnfet.w90: Wigner–Seitz shift folding (use_ws_distance) and 3D H(k).

The real-data check (apply_wsvec reproduces Wannier90's own _band.dat to 0.04 meV for the A1 bulk runs) lives in
dft/A1; these are self-contained identities.
"""
import numpy as np

from scnfet import w90


def two_orbital_hr():
    R = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
    H = np.zeros((3, 2, 2), complex)
    H[0] = [[0.3, 0.1], [0.1, -0.2]]
    H[1] = [[-1.0, 0.4], [0.2, -0.5]]
    H[2] = np.conj(H[1].T)
    return w90.HR(R=R, deg=np.ones(3), H=H)


def test_trivial_wsvec_changes_nothing():
    hr = two_orbital_hr()
    ws = {(tuple(R), m, n): np.zeros((1, 3), int) for R in hr.R for m in range(2) for n in range(2)}
    k = np.random.default_rng(0).random((9, 3))
    assert np.allclose(w90.hk_frac(w90.apply_wsvec(hr, ws), k), w90.hk_frac(hr, k))


def test_wsvec_shift_averages_the_phase():
    """One element split over two equivalent lattice vectors R and R + T: H(k) = h e^{2πik·R}(1 + e^{2πik·T})/2."""
    R = np.array([[0, 0, 0], [2, 0, 0], [-2, 0, 0]])
    H = np.zeros((3, 1, 1), complex)
    H[1, 0, 0] = H[2, 0, 0] = 0.6
    hr = w90.HR(R=R, deg=np.ones(3), H=H)
    ws = {((2, 0, 0), 0, 0): np.array([[0, 0, 0], [-4, 0, 0]]), ((-2, 0, 0), 0, 0): np.array([[0, 0, 0], [4, 0, 0]])}
    k = np.array([[0.13, 0, 0], [0.31, 0, 0]])
    got = w90.hk_frac(w90.apply_wsvec(hr, ws), k)[:, 0, 0].real
    x = 2 * np.pi * k[:, 0]
    want = 0.6 * (np.cos(2 * x) + np.cos(-2 * x))      # both halves: (e^{i2x} + e^{-i2x})/2 + (e^{-i2x} + e^{i2x})/2
    assert np.allclose(got, want)


def test_hk_frac_matches_hk_1d_along_axis():
    hr = two_orbital_hr()
    k = np.linspace(-0.5, 0.5, 11)
    assert np.allclose(w90.hk_frac(hr, np.c_[k, 0 * k, 0 * k]), w90.hk_1d(hr, k, axis=0))


def test_regauge_preserves_bands_and_bond_vectors():
    """A 2-site fcc model written in one cell/home-cell/order is moved to another; the bands at the same Cartesian
    k must agree, and every element must keep its Cartesian bond vector."""
    a = 4.5
    cA = 0.5 * a * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], float)
    cB = 0.5 * a * np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]], float)
    posA = np.array([[0, 0, 0], [a / 2, 0, 0]])                 # orbitals: site 0, site 1
    posB = np.array([[-a / 2, a / 2, a / 2], [0, 0, 0]])        # same sites, other home cell, swapped order
    rng = np.random.default_rng(1)
    R = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, -1], [0, -1, 1]])
    H = np.zeros((5, 2, 2), complex)
    H[0] = [[0.2, 0.7], [0.7, -0.4]]
    for i, j in ((1, 2), (3, 4)):
        H[i] = rng.normal(size=(2, 2))
        H[j] = H[i].T
    hrA = w90.HR(R=R, deg=np.ones(5), H=H)
    hrB = w90.regauge(hrA, cA, posA, cB, posB, perm=[1, 0])
    k = rng.normal(size=(6, 3))
    eA = np.linalg.eigvalsh(w90.hk_frac(hrA, k @ cA.T / (2 * np.pi)))
    eB = np.linalg.eigvalsh(w90.hk_frac(hrB, k @ cB.T / (2 * np.pi)))
    assert np.allclose(eA, eB)
    back = w90.regauge(hrB, cB, posB, cA, posA, perm=[1, 0])
    order = np.lexsort(back.R.T[::-1])
    nz = np.abs(back.H[order]).sum(axis=(1, 2)) > 0
    assert np.array_equal(back.R[order][nz], R[np.lexsort(R.T[::-1])]) and np.allclose(
        back.H[order][nz], H[np.lexsort(R.T[::-1])])
