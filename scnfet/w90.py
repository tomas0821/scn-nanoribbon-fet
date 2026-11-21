"""Wannier90 `_hr.dat` reader and 1D band structures along one lattice direction.

H(k) = sum_R exp(2*pi*i k R_t) H(R) / deg(R), with k in fractional units of the transport reciprocal vector. Any
R components in the confined directions get phase 0 (k_perp = 0). That is right for an isolated ribbon cell,
whose Wannier k-mesh is 1 along the confined directions, so those R components are zero anyway.
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class HR:
    R: np.ndarray      # (nR, 3) int lattice vectors
    deg: np.ndarray    # (nR,) Wigner–Seitz degeneracies
    H: np.ndarray      # (nR, nw, nw) complex, eV

    @property
    def num_wann(self):
        return self.H.shape[1]


def read_hr(path):
    with open(path) as f:
        f.readline()                                  # header / date
        nw = int(f.readline())
        nR = int(f.readline())
        deg = []
        while len(deg) < nR:
            deg += [int(x) for x in f.readline().split()]
        data = np.loadtxt(f)
    if data.shape[0] != nR * nw * nw:
        raise ValueError(f"{path}: expected {nR * nw * nw} matrix lines, found {data.shape[0]}")
    data = data.reshape(nR, nw * nw, 7)
    R = data[:, 0, :3].astype(int)
    m = data[:, :, 3].astype(int) - 1
    n = data[:, :, 4].astype(int) - 1
    H = np.zeros((nR, nw, nw), complex)
    for i in range(nR):
        H[i, m[i], n[i]] = data[i, :, 5] + 1j * data[i, :, 6]
    return HR(R=R, deg=np.asarray(deg, float), H=H)


def write_hr(path, hr, header="written by scnfet.w90"):
    """Inverse of read_hr (used by the tests to build synthetic models)."""
    nR, nw = len(hr.R), hr.num_wann
    with open(path, "w") as f:
        f.write(header + "\n" + f"{nw}\n{nR}\n")
        d = [str(int(x)) for x in hr.deg]
        for i in range(0, nR, 15):
            f.write(" ".join(d[i:i + 15]) + "\n")
        for i in range(nR):
            for n in range(nw):                       # W90 order: m fastest
                for m in range(nw):
                    z = hr.H[i, m, n]
                    f.write(f"{hr.R[i, 0]:5d}{hr.R[i, 1]:5d}{hr.R[i, 2]:5d}{m + 1:5d}{n + 1:5d}"
                            f"{z.real:16.10f}{z.imag:16.10f}\n")


def read_wsvec(path):
    """Wannier90 `_wsvec.dat` (use_ws_distance): {(R tuple, m, n): (nT, 3) int shifts}, with m, n 0-based."""
    out = {}
    with open(path) as f:
        f.readline()                                  # header
        while True:
            line = f.readline()
            if not line.strip():
                break
            r1, r2, r3, m, n = (int(v) for v in line.split())
            nT = int(f.readline())
            T = np.array([[int(v) for v in f.readline().split()] for _ in range(nT)])
            out[((r1, r2, r3), m - 1, n - 1)] = T
    return out


def apply_wsvec(hr, ws):
    """Fold the Wigner–Seitz shifts into an equivalent H(R) with deg = 1:
    H'(R + T)_mn += H(R)_mn / (deg(R) · N_T) for every shift T of (R, m, n). The result interpolates exactly as
    Wannier90 does with use_ws_distance = .true., so hk_1d / bands_1d / the ToB code need no special case."""
    nw = hr.num_wann
    acc = {}
    for i, R in enumerate(hr.R):
        Rt = tuple(int(v) for v in R)
        for m in range(nw):
            for n in range(nw):
                h = hr.H[i, m, n] / hr.deg[i]
                if h == 0:
                    continue
                T = ws.get((Rt, m, n))
                if T is None:
                    T = np.zeros((1, 3), int)
                for t in T:
                    key = (Rt[0] + int(t[0]), Rt[1] + int(t[1]), Rt[2] + int(t[2]))
                    acc.setdefault(key, np.zeros((nw, nw), complex))[m, n] += h / len(T)
    keys = sorted(acc)
    return HR(R=np.array(keys, int), deg=np.ones(len(keys)), H=np.array([acc[k] for k in keys]))


def hk_frac(hr, kfrac):
    """H(k) for 3D fractional k (nk, 3): sum_R exp(2πi k·R) H(R) / deg(R), Hermitised."""
    kfrac = np.atleast_2d(kfrac)
    phase = np.exp(2j * np.pi * kfrac @ hr.R.T) / hr.deg
    Hk = np.einsum("kr,rmn->kmn", phase, hr.H)
    return 0.5 * (Hk + np.conj(np.transpose(Hk, (0, 2, 1))))


def hk_1d(hr, k, axis=0):
    """H(k) for fractional k (array) along lattice `axis`. Returns (nk, nw, nw)."""
    k = np.atleast_1d(k)
    phase = np.exp(2j * np.pi * np.outer(k, hr.R[:, axis])) / hr.deg     # (nk, nR)
    Hk = np.einsum("kr,rmn->kmn", phase, hr.H)
    return 0.5 * (Hk + np.conj(np.transpose(Hk, (0, 2, 1))))


def bands_1d(hr, nk=4000, axis=0, chunk=500):
    """Sorted eigenvalues on a uniform periodic grid k in [-1/2, 1/2). Returns (k, E) with E (nk, nw)."""
    k = (np.arange(nk) / nk) - 0.5
    E = np.empty((nk, hr.num_wann))
    for i in range(0, nk, chunk):
        E[i:i + chunk] = np.linalg.eigvalsh(hk_1d(hr, k[i:i + chunk], axis))
    return k, E


# W90 real-harmonic order for `X:d` / `X:p` projections (user guide, Tables 3.1-3.2)
W90_ORDER = {"d": ("dz2", "dxz", "dyz", "dx2-y2", "dxy"), "p": ("pz", "px", "py")}


def regauge(hr, cell_from, pos_from, cell_to, pos_to, perm=None, signs=None):
    """Re-express H(R) in another lattice basis, home-cell choice and orbital order.

    Element (R, perm[i], perm[j]) of `hr`, a hop from orbital perm[i] at pos_from[perm[i]] to orbital perm[j] at
    pos_from[perm[j]] + R·cell_from, becomes element (R', i, j) with the same Cartesian bond vector:
    pos_to[j] + R'·cell_to − pos_to[i] = pos_from[perm[j]] + R·cell_from − pos_from[perm[i]].
    pos_* are per-orbital centres (Å, rows), cell_* rows are lattice vectors. `signs` (±1 per new orbital)
    flips orbital phases. H is divided by deg (the result has deg = 1); the map is one-to-one for fixed (i, j),
    so no two elements collide. Raises if R' is not an integer vector (the two cells don't describe one lattice).
    """
    nw = hr.num_wann
    perm = np.arange(nw) if perm is None else np.asarray(perm)
    s = np.ones(nw) if signs is None else np.asarray(signs, float)
    cell_from, cell_to = np.asarray(cell_from, float), np.asarray(cell_to, float)
    pf, pt = np.asarray(pos_from, float)[perm], np.asarray(pos_to, float)
    inv = np.linalg.inv(cell_to)
    Rc = (hr.R @ cell_from) @ inv                                                     # (nR, 3)
    shift = ((pf[None, :, :] - pf[:, None, :]) - (pt[None, :, :] - pt[:, None, :])) @ inv   # (nw, nw, 3)
    Rn = Rc[:, None, None, :] + shift[None]                                          # (nR, nw, nw, 3)
    Ri = np.rint(Rn).astype(int)
    if np.max(np.abs(Rn - Ri)) > 1e-4:
        raise ValueError(f"non-integer R' (max {np.max(np.abs(Rn - Ri)):.3g}): cells or centres inconsistent")
    Hs = hr.H[:, perm][:, :, perm] / hr.deg[:, None, None] * np.outer(s, s)[None]
    keys, inv_idx = np.unique(Ri.reshape(-1, 3), axis=0, return_inverse=True)
    H = np.zeros((len(keys), nw, nw), complex)
    ii, jj = np.meshgrid(np.arange(nw), np.arange(nw), indexing="ij")
    rr = inv_idx.reshape(len(hr.R), nw, nw)
    for r in range(len(hr.R)):
        H[rr[r], ii, jj] = Hs[r]
    return HR(R=keys.astype(int), deg=np.ones(len(keys)), H=H)
