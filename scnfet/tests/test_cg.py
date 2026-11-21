import numpy as np

from scnfet import cg


def test_microstrip_matches_hammerstad():
    """Zero-thickness strip over a ground plane, homogeneous medium: within 2 % of the closed form (≈ 1 % accurate)."""
    for W in (0.5, 1.0, 3.0, 6.0):
        c = cg.cg_single(W, 1.0, d=0.0, eps_below=1.0, eps_above=1.0, h=0.02)
        assert abs(c / cg.microstrip_hammerstad(W, 1.0, 1.0) - 1) < 0.02


def test_wide_limits():
    pp = cg.parallel_plate(200.0, 1.0)
    single = cg.cg_single(200.0, 1.0, d=0.0, eps_below=3.9, eps_above=3.9, h=0.05) / pp
    gaa = cg.cg_gaa(200.0, 1.0, d=0.0, eps=3.9, h=0.05) / (2 * pp)
    assert 1.0 < single < 1.05 and 1.0 < gaa < 1.02


def test_eot_scaling():
    """HfO2 of thickness EOT·25/3.9 tends to the parallel-plate value for wide ribbons."""
    t = 1.0 * 25 / 3.9
    r = cg.cg_single(400.0, t, eps_below=25.0, eps_above=25.0, h=0.1) / cg.parallel_plate(400.0, 1.0)
    assert 1.0 < r < 1.1
