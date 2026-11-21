"""Gate capacitance of a periodic ribbon array under one planar gate (fits/PREREG_array.md).

An infinite array of identical ribbons (width W, thickness d, edge-to-edge gap s, pitch P = W + s), all at the same
potential, at height t above a grounded planar gate, encapsulated in one dielectric (eps_below / eps_above as in
`cg.cg_single`). By symmetry the planes x = ±P/2 midway between ribbons carry no normal field, so one period with
zero-flux side boundaries is exact. Returns C_G per ribbon per unit length (F/m), the quantity the ToB code takes.

Limits used as tests: s = 0 is a continuous sheet (C = ε W / t exactly, the parallel-plate value); s → ∞ is the
isolated ribbon of `cg.cg_single`.

GAA needs no array solve: for s ≥ 2t each ribbon's gate box encloses it, so the array value equals `cg.cg_gaa`.
Uses the same finite-volume solver and grid rules as `cg.py`. Lengths in nm.
"""
import numpy as np

from scnfet.cg import EPS0, _axis, _solve


def cg_single_array(W, s, t, d=0.217, eps_below=25.0, eps_above=25.0, h=0.02, extent=None):
    """Single planar gate, periodic array of pitch W + s. F/m per ribbon."""
    P = W + s
    L = extent or 40.0 * max(P, t)
    x = _axis([-P / 2, -W / 2, 0.0, W / 2, P / 2], h, -P / 2, P / 2)
    y = _axis([0.0, t, t + d], h, 0.0, t + d + L)
    yc = (y[:-1] + y[1:]) / 2
    eps = np.where(yc[None, :] < t, eps_below, eps_above) * np.ones((len(x) - 1, 1))
    X, Y = np.meshgrid(x, y, indexing="ij")
    tol = 1e-9
    cond = (np.abs(X) <= W / 2 + tol) & (Y >= t - tol) & (Y <= t + d + tol)
    gate = np.abs(Y) < tol
    return EPS0 * _solve(x, y, eps, cond, gate)
