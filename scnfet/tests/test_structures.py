"""Checks of scnfet.structures: counts, widths, bonding, and band folding through the SK builder."""
import numpy as np
import pytest

from scnfet import sk, structures as S

TOY = dict(Ep=-2.02, Et2g=2.71, Eeg=4.93, pd_s=-2.40, pd_p=1.38, dd1_s=-0.12, dd1_p=0.58, dd1_d=-0.22,
           pp_s=0.80, pp_p=-0.18, dd2_s=-0.03, dd2_p=0.05, dd2_d=-0.02)


def neighbours(st, i, rmax):
    """(distance, species) of all atoms within rmax of atom i, periodic images along the periodic axes."""
    out = []
    rng = [range(-1, 2) if p else [0] for p in st.pbc]
    for a in rng[0]:
        for b in rng[1]:
            for c in rng[2]:
                T = np.array([a, b, c]) @ st.cell
                d = np.linalg.norm(st.pos + T - st.pos[i], axis=1)
                out += [(float(x), st.species[j]) for j, x in enumerate(d) if 1e-6 < x < rmax]
    return out


@pytest.mark.parametrize("direction, n, natoms, W_nm", [("100", 3, 24, 1.1), ("100", 6, 48, 2.4), ("100", 14, 112, 5.9),
                                                        ("110", 4, 16, 1.1), ("110", 8, 32, 2.3), ("110", 20, 80, 6.1)])
def test_counts_and_widths_match_the_plan(direction, n, natoms, W_nm):
    r = S.ribbon(direction, n)
    assert len(r.species) == natoms
    assert S.counts(r) == {"N": natoms // 2, "Sc": natoms // 2}            # stoichiometric
    assert S.width(r) / 10 == pytest.approx(W_nm, abs=0.05)
    a_c = S.LIANG["a_p"] * np.sqrt(2)
    step = a_c / 2 if direction == "100" else S.LIANG["a_p"] / 2
    assert S.width(r) == pytest.approx((2 * n - 1) * step)


@pytest.mark.parametrize("direction", ["100", "110"])
def test_rocksalt_bonding(direction):
    """Every nearest neighbour is the other species: 2.20 Å in plane, 2.18 Å between planes. Interior atoms are
    5-fold. Bare edge atoms are 4-fold on [100] (one in-plane bond cut) but 3-fold on [110], where the in-plane bonds
    run diagonally to the edge and two are cut. [110] edge rows are single-species within a plane."""
    r = S.ribbon(direction, 4)
    a_in = S.LIANG["a_p"] / np.sqrt(2)
    y = r.pos[:, 1]
    for i, s in enumerate(r.species):
        nb = neighbours(r, i, 2.5)
        assert all(sp != s for _, sp in nb)
        assert all(min(abs(dd - a_in), abs(dd - S.LIANG["d"])) < 1e-9 for dd, _ in nb)
        edge = abs(y[i] - y.min()) < 1e-6 or abs(y[i] - y.max()) < 1e-6
        assert len(nb) == ((4 if direction == "100" else 3) if edge else 5)
    if direction == "110":
        for p in range(2):
            z = sorted(set(np.round(r.pos[:, 2], 6)))[p]
            row = [s for s, q in zip(r.species, r.pos) if abs(q[1] - y.min()) < 1e-6 and abs(q[2] - z) < 1e-6]
            assert len(set(row)) == 1
        lower = [s for s, q in zip(r.species, r.pos) if abs(q[1] - y.min()) < 1e-6]
        assert set(lower) == {"Sc", "N"}                                     # opposite polarity in the two planes


def test_passivation_adds_one_H_per_edge_atom():
    r = S.ribbon("100", 3, passivate=True)
    assert S.counts(r)["H"] == 2 * 2 * 2                                     # 2 edges × 2 atoms per row × 2 planes
    bare = S.ribbon("100", 3)
    # each [100] edge row holds both species, so the outermost nucleus on each side is the Sc–H hydrogen
    assert S.width(r) == pytest.approx(S.width(bare) + 2 * S.BOND_H["Sc"])


@pytest.mark.parametrize("direction, n", [("110", 3), ("100", 2)])
def test_periodic_ribbon_is_a_folded_sheet(direction, n):
    """A width-periodic ribbon is a sheet supercell: at (k_x, 0) its SK eigenvalues equal the sheet's at the n
    folded k_y = 2π m / (n · period_y). Checks geometry, species pattern and the builder together."""
    rib = S.ribbon(direction, n, periodic_width=True)
    sheet = S.slab(2, cell="primitive" if direction == "110" else "conventional")
    mr, ms = sk.build(rib), sk.build(sheet)
    xr, xs = mr.x_from(TOY), ms.x_from(TOY)
    Py = rib.cell[1, 1]
    for kx in (0.0, 0.37, 0.81):
        er = np.sort(mr.eigvals(xr, np.array([[kx, 0, 0]]))[0])
        ks = np.array([[kx, 2 * np.pi * m / Py, 0] for m in range(n)])
        es = np.sort(ms.eigvals(xs, ks).ravel())
        assert np.allclose(er, es, atol=1e-9)


def test_sheet_slabs():
    for npl in (2, 3):
        st = S.slab(npl)
        assert len(st.species) == 2 * npl
        z = np.unique(np.round(st.pos[:, 2], 6))
        assert np.allclose(np.diff(z), S.LIANG["d"])
    conv = S.slab(2, cell="conventional")
    assert np.isclose(conv.cell[0, 0], S.LIANG["a_p"] * np.sqrt(2))


def test_width_table():
    t = S.width_table()
    assert t["100"][6] == pytest.approx(2.42, abs=0.01) and t["110"][8] == pytest.approx(2.33, abs=0.01)
