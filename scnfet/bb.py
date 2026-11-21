"""Building-block Wannier ribbons (after Lee, Nardelli & Marzari, PRL 95, 076804 (2005)).

A wide ribbon's Wannier Hamiltonian is assembled from two smaller Wannier calculations, with no fitting:
  * the **edge source**, a narrow calibration ribbon with the same direction and edge chemistry: every matrix element
    whose two orbitals both sit within the calibration ribbon's half-width of the same edge is copied from the
    matched pair of the calibration ribbon (matched by row counted from that edge, plane and position along x);
  * the **interior source**, the periodic sheet: every other matrix element is the sheet's element for the same
    species, planes and bond vector, with the orbitals rotated from the sheet's axes into the ribbon's.
The two calculations have different energy zeros, so each source's on-site energies are referred to its own vacuum
level (planar-averaged bare + Hartree potential in the vacuum, as the A2/A4b window rules used).

Everything is in the Wannier90 real-harmonic order (d: z², xz, yz, x²−y², xy; p: z, x, y), so the assembled
Hamiltonian is directly comparable with a DFT-Wannier ribbon. Units: Å, eV.
"""
from dataclasses import dataclass

import numpy as np

from scnfet import w90 as W

ORB = {"Sc": ("dz2", "dxz", "dyz", "dx2-y2", "dxy"), "N": ("pz", "px", "py"), "H": ("s",)}
_DF = {"dz2": np.diag([-1.0, -1, 2]) / (2 * np.sqrt(3)),
       "dxz": np.array([[0, 0, .5], [0, 0, 0], [.5, 0, 0]]),
       "dyz": np.array([[0, 0, 0], [0, 0, .5], [0, .5, 0]]),
       "dx2-y2": np.diag([.5, -.5, 0]),
       "dxy": np.array([[0, .5, 0], [.5, 0, 0], [0, 0, 0]])}
_PI = {"px": 0, "py": 1, "pz": 2}


def orb_rot(Q, names):
    """C with φ_i(Qᵀ r) = Σ_j C[j, i] φ_j(r) for the real orbitals `names` (one l). Orbitals defined on axes rotated
    by Q (a vector with components v there has components Q v here) are φ'_i = Σ_j C[j, i] φ_j, so H' = Cᵀ H C."""
    Q = np.asarray(Q, float)
    if names[0].startswith("d"):
        return np.array([[2 * np.trace(_DF[j] @ Q @ _DF[i] @ Q.T) for i in names] for j in names])
    if names[0].startswith("p"):
        return np.array([[Q[_PI[j], _PI[i]] for i in names] for j in names])
    return np.eye(len(names))


@dataclass
class Source:
    """One Wannier calculation, flattened to deg = 1, with each WF labelled by its atom and orbital.

    pos: the atom positions the WF centres sit on (the .win positions); cell rows are lattice vectors; E_vac: the
    vacuum level in the calculation's own energy zero."""
    hr: W.HR
    cell: np.ndarray
    pos: np.ndarray
    species: list
    wf: list                 # [(atom, orbital name)]
    E_vac: float

    def __post_init__(self):
        if not np.allclose(self.hr.deg, 1):
            self.hr = W.HR(R=self.hr.R, deg=np.ones(len(self.hr.R)), H=self.hr.H / self.hr.deg[:, None, None])
        self.idx = {}
        for j, (a, o) in enumerate(self.wf):
            self.idx.setdefault(a, {})[o] = j
        self.iR = {tuple(int(v) for v in R): i for i, R in enumerate(self.hr.R)}

    def block(self, a, b, R, oa, ob):
        """⟨a, oa; 0 | H | b, ob; R⟩ for orbital-name lists oa, ob (zero rows/columns for orbitals a source lacks;
        zero if R is outside the source's range). On-site diagonals are referred to the vacuum level."""
        i = self.iR.get(tuple(R))
        out = np.zeros((len(oa), len(ob)))
        if i is None:
            return out
        ra = [(p, self.idx[a][o]) for p, o in enumerate(oa) if o in self.idx.get(a, {})]
        rb = [(q, self.idx[b][o]) for q, o in enumerate(ob) if o in self.idx.get(b, {})]
        if not ra or not rb:
            return out
        Hs = self.hr.H[i]
        blk = Hs[np.ix_([j for _, j in ra], [j for _, j in rb])]
        # real gauge expected (time reversal); r110_4h HW1 carries 0.4 meV of numerical noise, the others < 0.1 meV
        assert np.abs(blk.imag).max() < 1e-3, "complex hopping in a real (time-reversal) Wannier gauge"
        out[np.ix_([p for p, _ in ra], [q for q, _ in rb])] = blk.real
        if a == b and not any(R):
            for p, o in enumerate(oa):
                if o in self.idx.get(a, {}) and o in ob:
                    out[p, ob.index(o)] -= self.E_vac
        return out


def wf_list(species, s_atoms=()):
    """The WF list in make_a4b order: Sc d (atom order), N p, H s, then Sc s on `s_atoms`."""
    wf = [(i, o) for i, s in enumerate(species) if s == "Sc" for o in ORB["Sc"]]
    wf += [(i, o) for i, s in enumerate(species) if s == "N" for o in ORB["N"]]
    wf += [(i, "s") for i, s in enumerate(species) if s == "H"]
    wf += [(i, "s") for i in s_atoms]
    return wf


@dataclass
class Sheet:
    """The periodic sheet as the interior source: lookup by (species, plane) of both atoms and the in-plane bond in
    units of a_p/2 (every M-ScN site is on that half lattice)."""
    src: Source
    a_p: float
    plane_z: float           # z separating the two planes

    def __post_init__(self):
        s = self.src
        self.cls = [(sp, int(z > self.plane_z)) for sp, z in zip(s.species, s.pos[:, 2])]
        assert len(set(self.cls)) == len(self.cls), "the sheet needs one atom per (species, plane)"
        self.atom_of = {c: i for i, c in enumerate(self.cls)}
        self.table = {}
        for R in s.hr.R:
            for m in range(len(s.species)):
                for n in range(len(s.species)):
                    d = s.pos[n] + R @ s.cell - s.pos[m]
                    u = 2 * d[:2] / self.a_p
                    ui = np.rint(u).astype(int)
                    assert np.abs(u - ui).max() < 1e-3
                    key = (self.cls[m], self.cls[n], int(ui[0]), int(ui[1]))
                    assert key not in self.table
                    self.table[key] = (m, n, tuple(int(v) for v in R), float(np.linalg.norm(d[:2])))

    def block(self, ca, cb, d_sheet, oa, ob):
        """Sheet element for atoms of class ca → cb at in-plane bond d_sheet (sheet axes); None if out of range."""
        u = 2 * np.asarray(d_sheet[:2]) / self.a_p
        ui = np.rint(u).astype(int)
        if np.abs(u - ui).max() > 0.05:
            raise ValueError(f"bond {d_sheet} is not on the sheet's half lattice")
        hit = self.table.get((ca, cb, int(ui[0]), int(ui[1])))
        if hit is None:
            return None
        m, n, R, _ = hit
        return self.src.block(m, n, R, oa, ob)


def rows_planes(ideal_pos, species):
    """Row index by ideal y (H included: an H row is its own outermost row) and plane index by ideal z."""
    y, z = np.round(ideal_pos[:, 1], 3), np.round(ideal_pos[:, 2], 3)
    host = np.array([s != "H" for s in species])
    uy = np.unique(y)
    uz = np.unique(z[host])
    assert len(uz) == 2, f"{len(uz)} host planes"
    plane = np.where(host, np.searchsorted(uz, z), np.argmin(np.abs(z[:, None] - uz[None]), axis=1))
    return np.searchsorted(uy, y), len(uy), plane


@dataclass
class Ribbon:
    """Ideal-cut description of a ribbon: transport period Lx (along x), positions, species, rows and planes."""
    pos: np.ndarray
    species: list
    Lx: float

    def __post_init__(self):
        self.row, self.nrow, self.plane = rows_planes(self.pos, self.species)


def _xkey(x, L):
    v = float(x) % L
    return 0.0 if L - v < 1e-3 else round(v, 2)


def assemble(target, cal, cal_src, sheet, Q_sheet, rcut, s_atoms=(), h=None, sheet_shift=0.0):
    """Building-block H(R) for `target` (Ribbon), from the calibration ribbon `cal` (Ribbon, its ideal cut) whose
    Wannier Hamiltonian is `cal_src` (Source, same atom order), and the `sheet` (Sheet). Q_sheet maps target-frame
    vectors to sheet-frame vectors (d_sheet = Q_sheet d). Elements with |d| > rcut are zero.

    Matching: target atom a at row r maps to the calibration atom with the same plane, the same x (mod Lx) and row r
    (left edge) or r − (N_T − N_C) (right edge). An element is copied from the calibration ribbon when both atoms
    lie in the same edge region (row < N_C / 2 from that edge); every other element comes from the sheet.
    s_atoms: target atoms that carry an Sc s WF (the calibration ribbon's edge s, matched the same way).
    h: edge-region depth in rows (default N_C / 2, the whole calibration half; smaller values are in-sample tests).
    sheet_shift: a constant added to the sheet's on-site energies (eV), aligning its zero to the calibration
    ribbon's (see `centre_shift`); 0 keeps the two vacuum levels aligned.
    Returns (HR with deg 1, WF list, source map {(a, b, Rx): "L" | "R" | "S"})."""
    T, C = target, cal
    assert abs(T.Lx - C.Lx) < 1e-3, "different transport periods"
    h = C.nrow // 2 if h is None else h
    assert 1 <= h <= C.nrow // 2
    dN = T.nrow - C.nrow
    assert dN >= 0 and dN % 2 == 0
    key_c = {(int(C.row[i]), int(C.plane[i]), _xkey(C.pos[i, 0], C.Lx), C.species[i]): i
             for i in range(len(C.species))}

    def cal_atom(a, side):
        r = int(T.row[a]) if side == "L" else int(T.row[a]) - dN
        return key_c[(r, int(T.plane[a]), _xkey(T.pos[a, 0], T.Lx), T.species[a])]

    side = np.where(T.row < h, "L", np.where(T.row >= T.nrow - h, "R", "I"))
    wf = wf_list(T.species, s_atoms)
    orbs = {}
    for j, (a, o) in enumerate(wf):
        orbs.setdefault(a, []).append((j, o))
    cls_t = [(sp, int(p)) for sp, p in zip(T.species, T.plane)]
    rot = {}
    nmax = int(np.ceil(rcut / T.Lx)) + 1
    Rs = list(range(-nmax, nmax + 1))
    nw = len(wf)
    H = np.zeros((len(Rs), nw, nw))
    used = {}
    for ir, Rx in enumerate(Rs):
        for a in range(len(T.species)):
            ia, oa = zip(*orbs[a])
            for b in range(len(T.species)):
                d = T.pos[b] + np.array([Rx * T.Lx, 0, 0]) - T.pos[a]
                if np.linalg.norm(d) > rcut:
                    continue
                ib, ob = zip(*orbs[b])
                sa, sb = side[a], side[b]
                if sa == sb and sa in ("L", "R"):
                    a_c, b_c = cal_atom(a, sa), cal_atom(b, sb)
                    dc = cal_src.pos[b_c] - cal_src.pos[a_c]
                    Rc = int(np.rint((d[0] - dc[0]) / C.Lx))
                    # sanity only (atoms are matched by ideal row): relaxed H moves up to 0.7 Å in y (A4a)
                    assert abs(d[0] - dc[0] - Rc * C.Lx) < C.Lx / 4 and abs(d[1] - dc[1]) < 2.0, "atom match off"
                    blk = cal_src.block(a_c, b_c, (Rc, 0, 0), list(oa), list(ob))
                    used[(a, b, Rx)] = sa
                else:
                    if T.species[a] == "H" or T.species[b] == "H":
                        continue              # H only ever pairs inside an edge region (no H in the sheet)
                    ds = Q_sheet @ d
                    oa_r, ob_r = [o for o in oa if o != "s"], [o for o in ob if o != "s"]
                    blk = np.zeros((len(oa), len(ob)))
                    if oa_r and ob_r:
                        # sheet orbitals are on sheet axes: rotate into the target's (H' = Cᵀ H C)
                        for o_ in (tuple(oa_r), tuple(ob_r)):
                            if o_ not in rot:
                                rot[o_] = orb_rot(Q_sheet, list(o_))
                        Ca, Cb = rot[tuple(oa_r)], rot[tuple(ob_r)]
                        bs = sheet.block(cls_t[a], cls_t[b], ds, oa_r, ob_r)
                        if bs is not None:
                            pa = [oa.index(o) for o in oa_r]
                            pb = [ob.index(o) for o in ob_r]
                            blk[np.ix_(pa, pb)] = Ca.T @ bs @ Cb
                            if a == b and Rx == 0:
                                blk[pa, pa] += sheet_shift
                    used[(a, b, Rx)] = "S"
                H[ir][np.ix_(ia, ib)] = blk
    Rv = np.array([[r, 0, 0] for r in Rs])
    hr = W.HR(R=Rv, deg=np.ones(len(Rs)), H=H.astype(complex))
    herm = max(np.abs(H[i] - H[len(Rs) - 1 - i].T).max() for i in range(len(Rs)))
    assert herm < 1e-8, f"assembled H(R) not Hermitian ({herm:.2e})"
    return hr, wf, used


def centre_shift(cal, cal_src, sheet):
    """Sheet → calibration alignment (Lee/Nardelli/Marzari: match the on-site energies where both describe the same
    environment): the mean on-site diagonal of the calibration ribbon's two centre rows (all their atoms and
    orbitals) minus the sheet's for the same atoms. Both sources are vacuum-referred, so this is the residual offset."""
    C = cal
    cen = [a for a in range(len(C.species)) if C.row[a] in (C.nrow // 2 - 1, C.nrow // 2) and C.species[a] != "H"]
    mine, ref = [], []
    for a in cen:
        o = list(ORB[C.species[a]])
        cls = (C.species[a], int(C.plane[a]))
        mine += list(np.diag(cal_src.block(a, a, (0, 0, 0), o, o)))
        ref += list(np.diag(sheet.block(cls, cls, np.zeros(3), o, o)))
    return float(np.mean(mine) - np.mean(ref))


def edge_map(target, cal, h=None):
    """For each target atom: (side, calibration atom) if it lies in an edge region ("L" / "R", as `assemble` matches
    it), else ("I", None). Added for T1 (geometry and bookkeeping); `assemble` is unchanged."""
    T, C = target, cal
    h = C.nrow // 2 if h is None else h
    dN = T.nrow - C.nrow
    key_c = {(int(C.row[i]), int(C.plane[i]), _xkey(C.pos[i, 0], C.Lx), C.species[i]): i
             for i in range(len(C.species))}
    out = []
    for a in range(len(T.species)):
        r = int(T.row[a])
        if r < h:
            out.append(("L", key_c[(r, int(T.plane[a]), _xkey(T.pos[a, 0], T.Lx), T.species[a])]))
        elif r >= T.nrow - h:
            out.append(("R", key_c[(r - dN, int(T.plane[a]), _xkey(T.pos[a, 0], T.Lx), T.species[a])]))
        else:
            out.append(("I", None))
    return out
