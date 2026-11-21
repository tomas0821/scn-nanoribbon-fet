"""Structure generators: the M-ScN sheet and n-plane rocksalt-(001) slabs, and [100] / [110] nanoribbons.

Geometry (all Å). The sheet is the rocksalt-(001) bilayer, P4/nmm: in each plane Sc and N alternate on a square
lattice of spacing a_p/√2 (the in-plane Sc–N bond), and each plane swaps species relative to the one below (rocksalt
stacking, planes `d` apart). a_p is the primitive (2-f.u.) lattice constant, a_c = a_p √2 the 4-f.u. cell.

Ribbon conventions (plan §3 D4; signed-off device conventions):
  * [100]: transport along the Sc–N bond (cubic [100]), period a_c. n counts 4-f.u. cells across, so a ribbon has
    2n rows of spacing a_c/2, mixed Sc/N rows (non-polar edges) and 8n atoms per cell (two planes).
    n = 3 → W = 1.1 nm, n = 6 → 2.4 nm, n = 14 → 5.9 nm.
  * [110]: transport along the primitive vector (cubic [110]), period a_p. n counts primitive cells across: 2n rows
    of spacing a_p/2, each row a single species within a plane (polar edges, opposite polarity in the two planes),
    and 4n atoms per cell. n = 4 → 1.1 nm, n = 8 → 2.3 nm, n = 20 → 6.1 nm.
  * W is edge to edge between the outermost nuclei, edge H included (signed off 2026-11-11). `width()` measures it
    from the structure, so it stays right for relaxed or H-terminated edges.
  * Transport is along x, the width along y, the planes along z. Ribbons are periodic along x only.

Lattice constants: Liang's PBE-D3 values (a_p 3.11, d 2.18) are the defaults, as in A0. Our PBE relaxation (A0,
PAW 90 Ry) gives a_p 3.118, d 2.177 (`PBE_A0`).
"""
import numpy as np

from scnfet.sk import Structure

LIANG = dict(a_p=3.11, d=2.18)
PBE_A0 = dict(a_p=3.118, d=2.177)
BOND_H = {"Sc": 1.83, "N": 1.02}        # Å, starting Sc–H / N–H bond lengths for H-terminated edges (to be relaxed)
# The rocksalt cubic axes (columns) in a frame whose x, y run along the primitive vectors, i.e. with the Sc–N bonds
# along (1, ±1, 0)/√2: the primitive sheet/slab cell and the [110] ribbon. The SK on-site t2g/e_g split is defined on
# these axes (sk.crystal_field_projectors); the conventional cell and the [100] ribbon are already on them.
R45 = np.array([[1, -1, 0], [1, 1, 0], [0, 0, np.sqrt(2)]]) / np.sqrt(2)


def _species(parity):
    return "Sc" if parity % 2 == 0 else "N"


def slab(nplanes=2, a_p=LIANG["a_p"], d=LIANG["d"], vacuum=20.0, cell="primitive"):
    """n-plane rocksalt-(001) slab, periodic in plane. cell: "primitive" (2 atoms per plane, a_p × a_p) or
    "conventional" (4 per plane, a_c × a_c, axes along the Sc–N bonds)."""
    c = (nplanes - 1) * d + vacuum
    z0 = vacuum / 2
    pos, sp = [], []
    if cell == "primitive":
        L = np.array([[a_p, 0, 0], [0, a_p, 0], [0, 0, c]])
        sites = [(0.0, 0.0, 0), (a_p / 2, a_p / 2, 1)]                     # (x, y, sublattice parity)
    elif cell == "conventional":
        a_c = a_p * np.sqrt(2)
        L = np.array([[a_c, 0, 0], [0, a_c, 0], [0, 0, c]])
        sites = [(i * a_c / 2, j * a_c / 2, i + j) for i in range(2) for j in range(2)]
    else:
        raise ValueError(cell)
    for p in range(nplanes):
        for x, y, par in sites:
            pos.append([x, y, z0 + p * d])
            sp.append(_species(par + p))
    return Structure(L, sp, np.array(pos), pbc=(True, True, False), frame=R45 if cell == "primitive" else None)


def ribbon(direction, n, nplanes=2, a_p=LIANG["a_p"], d=LIANG["d"], vacuum=15.0, passivate=False,
           periodic_width=False):
    """Bare (or H-terminated) nanoribbon along x. vacuum: added both across the width (y) and out of plane (z).
    periodic_width=True makes y periodic with period equal to the 2n rows' extent: a sheet supercell, for tests."""
    a_c = a_p * np.sqrt(2)
    if direction == "100":
        period, dy = a_c, a_c / 2
        row = lambda j: [(0.0, j), (a_c / 2, j + 1)]                       # (x, parity) of the 2 atoms in row j
    elif direction == "110":
        period, dy = a_p, a_p / 2
        row = lambda j: [((j % 2) * a_p / 2, j % 2)]                       # one atom per row; species alternates by row
    else:
        raise ValueError("direction must be '100' or '110'")
    nrow = 2 * n
    span_y = nrow * dy if periodic_width else (nrow - 1) * dy + vacuum
    span_z = (nplanes - 1) * d + vacuum
    y0, z0 = (0.0 if periodic_width else vacuum / 2), vacuum / 2
    pos, sp = [], []
    for p in range(nplanes):
        for j in range(nrow):
            for x, par in row(j):
                pos.append([x, y0 + j * dy, z0 + p * d])
                sp.append(_species(par + p))
    pos = np.array(pos)
    if passivate:
        if periodic_width:
            raise ValueError("a periodic-width ribbon has no edges to passivate")
        pos, sp = _passivate(pos, sp, y0, y0 + (nrow - 1) * dy)
    L = np.array([[period, 0, 0], [0, span_y, 0], [0, 0, span_z]])
    return Structure(L, list(sp), pos, pbc=(True, bool(periodic_width), False),
                     frame=R45 if direction == "110" else None)


def _passivate(pos, sp, ylo, yhi):
    """One H per edge atom, along the edge normal (±y), at BOND_H. A starting geometry for the DFT relaxation (A4a),
    not a relaxed one; the SK model has no H basis yet.
    [100] edges: each edge atom has one cut in-plane bond, along ±y, so this saturates it.
    [110] edges: each edge atom has two cut bonds, diagonal to the edge. One H per cut bond would put neighbouring
    H atoms 0.5 Å apart, so the single H on the normal is a placeholder; the [110] termination is a chemistry
    decision for A4a."""
    pos, sp = list(map(np.array, pos)), list(sp)
    add_p, add_s = [], []
    for r, s in zip(pos, sp):
        for ye, sgn in ((ylo, -1), (yhi, +1)):
            if abs(r[1] - ye) < 1e-6:
                add_p.append(r + np.array([0, sgn * BOND_H[s], 0]))
                add_s.append("H")
    return np.array(pos + add_p), sp + add_s


def width(st, axis=1):
    """Edge-to-edge width (Å): the spread of the nuclei across the ribbon, H included."""
    return float(np.ptp(st.pos[:, axis]))


def counts(st):
    return {s: st.species.count(s) for s in sorted(set(st.species))}


def width_table(directions=("100", "110"), ns=None, **kw):
    """W (nm) for a list of n per direction (bare ribbons)."""
    ns = ns or {"100": [3, 4, 6, 8, 10, 14], "110": [4, 6, 8, 10, 14, 20]}
    return {dr: {n: width(ribbon(dr, n, **kw)) / 10 for n in ns[dr]} for dr in directions}
