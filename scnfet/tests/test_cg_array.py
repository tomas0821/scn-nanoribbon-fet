"""Periodic-array C_G (scnfet/cg_array.py): exact limits and monotonicity."""
import numpy as np

from scnfet import cg
from scnfet.cg_array import cg_single_array

T = 0.6 * 25 / 3.9


def test_zero_gap_is_parallel_plate():
    # a continuous sheet: C = eps W / t exactly, which is the parallel-plate convention at this EOT
    for W in (1.09, 5.93):
        assert abs(cg_single_array(W, 0.0, T) / cg.parallel_plate(W, 0.6) - 1) < 2e-3


def test_wide_gap_tends_to_isolated():
    W = 2.41
    assert abs(cg_single_array(W, 100.0, T) / cg.cg_single(W, T) - 1) < 1e-2


def test_monotonic_between_limits():
    W = 1.09
    c = np.array([cg_single_array(W, s, T) for s in (0.0, 1.0, 2.0, 4.0, 7.69, 16.0)])
    assert np.all(np.diff(c) > 0)
    assert c[-1] < cg.cg_single(W, T)
