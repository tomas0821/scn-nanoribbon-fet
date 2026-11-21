import numpy as np

from scnfet import bb


def rz(deg):
    t = np.radians(deg)
    return np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1.0]])


def test_orb_rot_orthogonal_and_composes():
    for names in (list(bb.ORB["Sc"]), list(bb.ORB["N"])):
        A, B = bb.orb_rot(rz(30), names), bb.orb_rot(rz(15), names)
        assert np.allclose(A.T @ A, np.eye(len(names)))
        assert np.allclose(bb.orb_rot(rz(30) @ rz(15), names), A @ B)


def test_orb_rot_90_degrees():
    """φ_i(Qᵀ r) for Q = 90° about z: (x, y) → (y, −x), so xz → yz, yz → −xz, xy → −xy, x² − y² → −(x² − y²)."""
    d = list(bb.ORB["Sc"])
    C = bb.orb_rot(rz(90), d)
    e = {o: i for i, o in enumerate(d)}
    assert np.isclose(C[e["dyz"], e["dxz"]], 1) and np.isclose(C[e["dxz"], e["dyz"]], -1)
    assert np.isclose(C[e["dxy"], e["dxy"]], -1) and np.isclose(C[e["dx2-y2"], e["dx2-y2"]], -1)
    assert np.isclose(C[e["dz2"], e["dz2"]], 1)
    p = list(bb.ORB["N"])
    Cp = bb.orb_rot(rz(90), p)
    q = {o: i for i, o in enumerate(p)}
    assert np.isclose(Cp[q["py"], q["px"]], 1) and np.isclose(Cp[q["px"], q["py"]], -1)


def test_xkey_wraps():
    assert bb._xkey(4.39999, 4.4) == 0.0 and bb._xkey(-1e-9, 4.4) == 0.0 and bb._xkey(2.2, 4.4) == 2.2
