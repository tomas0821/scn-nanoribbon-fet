"""Geometric gate capacitance C_G(W) of a nanoribbon by the 2D Laplace equation (plan D6, referee item K3).

Cross-section only (the ribbon is long along x), so C_G is per unit length (F/m). The ribbon is a perfect conductor
of width W (edge to edge, as everywhere in the project) and thickness d; the gate is grounded. Two geometries:

  * "single": a planar gate at y = 0, the ribbon's bottom face at y = t; dielectric eps_below for y < t and
    eps_above for y ≥ t (eps_above = eps_below: the ribbon encapsulated). Outer boundaries are zero-flux (Neumann),
    so every field line from the ribbon ends on the gate and C is the two-conductor capacitance.
  * "gaa": a rectangular gate at distance t from every face of the ribbon, dielectric eps inside.

Finite volumes on a rectilinear grid with step h at every corner and interface, growing geometrically (×1.15)
towards segment middles and outside. C = ε0 · Σ_edges g_e (V_a − V_b)², the discrete field energy at unit voltage, which equals the
charge for the finite-volume solution. Lengths in nm.

    parallel_plate(W, EOT) = ε0 · 3.9 · W / EOT        (the signed-off convention; the EOT reference is SiO2)
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

EPS0 = 8.8541878128e-12
EPS_SIO2 = 3.9


def parallel_plate(W_nm, eot_nm):
    """F/m (W and EOT in nm; the ratio makes the units cancel)."""
    return EPS0 * EPS_SIO2 * W_nm / eot_nm


def _segment(a, b, h, growth):
    """Points in (a, b]: step h at both ends, growing geometrically towards the middle."""
    L = b - a
    left, s, x = [], h, 0.0
    while x + s < L / 2:
        x += s
        left.append(x)
        s *= growth
    inner = np.array(left)
    pts = np.concatenate([inner, [L / 2] if L > 2 * h else [], L - inner[::-1], [L]])
    pts = np.unique(np.round(pts, 12))
    return a + pts[pts > 1e-12]


def _axis(fixed, h, lo, hi, growth=1.15):
    """Nodes: every fixed point; graded steps between consecutive fixed points (h at each fixed point, growing
    towards segment middles); then geometric growth from the outermost fixed points out to lo and hi."""
    f = np.unique(np.round(np.asarray(fixed, float), 10))
    pts = [f[0]]
    for a, b in zip(f[:-1], f[1:]):
        pts += list(_segment(a, b, h, growth))
    out = []
    step, x = h, f[-1]
    while x < hi - 1e-12:
        step *= growth
        x = min(x + step, hi)
        out.append(x)
    inn = []
    step, x = h, f[0]
    while x > lo + 1e-12:
        step *= growth
        x = max(x - step, lo)
        inn.append(x)
    return np.array(sorted(inn) + pts + out)


def _solve(x, y, eps_cell, conductor, gate):
    """eps_cell: (nx−1, ny−1) permittivity of each cell; conductor / gate: (nx, ny) boolean node masks (V = 1 / 0).
    Neumann elsewhere on the boundary. Returns C/ε0 (dimensionless, per unit length)."""
    nx, ny = len(x), len(y)
    dx, dy = np.diff(x), np.diff(y)
    # edge conductances: x-edges between (i, j) and (i+1, j); the dual face spans the half cells above and below
    ex = np.zeros((nx - 1, ny))
    ex[:, :-1] += eps_cell * dy[None, :] / 2
    ex[:, 1:] += eps_cell * dy[None, :] / 2
    gx = ex / dx[:, None]
    ey = np.zeros((nx, ny - 1))
    ey[:-1, :] += eps_cell * dx[:, None] / 2
    ey[1:, :] += eps_cell * dx[:, None] / 2
    gy = ey / dy[None, :]
    fixed = conductor | gate
    V = np.where(conductor, 1.0, 0.0)
    idx = -np.ones((nx, ny), int)
    free = ~fixed
    idx[free] = np.arange(free.sum())
    rows, cols, vals = [], [], []
    b = np.zeros(free.sum())
    diag = np.zeros(free.sum())

    def couple(ia, ja, ib, jb, g):
        a, c = idx[ia, ja], idx[ib, jb]
        ma, mc = a >= 0, c >= 0
        # a free, c free: off-diagonals both ways; a free, c fixed: rhs
        np.add.at(diag, a[ma], g[ma])
        np.add.at(diag, c[mc], g[mc])
        both = ma & mc
        rows.extend([a[both], c[both]])
        cols.extend([c[both], a[both]])
        vals.extend([-g[both], -g[both]])
        af = ma & ~mc
        np.add.at(b, a[af], g[af] * V[ib, jb][af])
        cf = mc & ~ma
        np.add.at(b, c[cf], g[cf] * V[ia, ja][cf])

    I, J = np.meshgrid(np.arange(nx - 1), np.arange(ny), indexing="ij")
    couple(I.ravel(), J.ravel(), I.ravel() + 1, J.ravel(), gx.ravel())
    I, J = np.meshgrid(np.arange(nx), np.arange(ny - 1), indexing="ij")
    couple(I.ravel(), J.ravel(), I.ravel(), J.ravel() + 1, gy.ravel())
    n = free.sum()
    A = sp.csr_matrix((np.concatenate(vals + [diag]), (np.concatenate(rows + [np.arange(n)]),
                                                     np.concatenate(cols + [np.arange(n)]))), shape=(n, n))
    V[free] = spl.spsolve(A.tocsc(), b)
    return float((gx * np.diff(V, axis=0) ** 2).sum() + (gy * np.diff(V, axis=1) ** 2).sum())


def cg_single(W, t, d=0.217, eps_below=25.0, eps_above=25.0, h=0.02, extent=None):
    """Single planar gate. Returns C_G per unit length (F/m). W, t, d in nm."""
    L = extent or 40.0 * max(W, t)
    x = _axis([-W / 2, 0.0, W / 2], h, -L, L)
    y = _axis([0.0, t, t + d], h, 0.0, t + d + L)
    xc, yc = (x[:-1] + x[1:]) / 2, (y[:-1] + y[1:]) / 2
    eps = np.where(yc[None, :] < t, eps_below, eps_above) * np.ones((len(xc), 1))
    X, Y = np.meshgrid(x, y, indexing="ij")
    tol = 1e-9
    cond = (np.abs(X) <= W / 2 + tol) & (Y >= t - tol) & (Y <= t + d + tol)
    gate = np.abs(Y) < tol
    return EPS0 * _solve(x, y, eps, cond, gate)


def cg_gaa(W, t, d=0.217, eps=25.0, h=0.02):
    """Rectangular gate all around, at distance t from every face. Returns F/m."""
    x = _axis([-W / 2 - t, -W / 2, 0.0, W / 2, W / 2 + t], h, -W / 2 - t, W / 2 + t)
    y = _axis([-t, 0.0, d, d + t], h, -t, d + t)
    X, Y = np.meshgrid(x, y, indexing="ij")
    tol = 1e-9
    cond = (np.abs(X) <= W / 2 + tol) & (Y >= -tol) & (Y <= d + tol)
    gate = (np.abs(np.abs(X) - (W / 2 + t)) < tol) | (np.abs(Y + t) < tol) | (np.abs(Y - d - t) < tol)
    epsc = np.full((len(x) - 1, len(y) - 1), float(eps))
    return EPS0 * _solve(x, y, epsc, cond, gate)


def microstrip_hammerstad(W, h_sub, eps_r):
    """Reference for the tests: zero-thickness strip at height h_sub over a ground plane in a homogeneous medium
    (Hammerstad–Jensen / Wheeler closed form, accurate to ≈ 1 %). Returns F/m."""
    u = W / h_sub
    eta0 = 376.730313668
    if u <= 1:
        Z0 = 60.0 / np.sqrt(eps_r) * np.log(8 / u + u / 4)
    else:
        Z0 = eta0 / np.sqrt(eps_r) / (u + 1.393 + 0.667 * np.log(u + 1.444))
    return np.sqrt(eps_r) / (299792458.0 * Z0)
