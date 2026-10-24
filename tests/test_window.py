"""T33: a reading's force and acceleration are averages over one window.

The theorem's example: F(t) = 5(1 + sin 2t) in a fixed direction, m = 2. Over the window
[0, 1] the pair gives exactly 2; with the acceleration window shifted by half a window,
1.935; with it doubled, 2.417. Noise-free, so the exact branch reads the pair.
"""
import numpy as np
import pytest

from met.full import exact_mass

U = np.array([1.0, 2.0, 2.0]) / 3.0          # any fixed direction


def mean_push(t0, t1):
    """Window average of 5(1 + sin 2t): 5 + 5 (cos 2t0 - cos 2t1) / (2 (t1 - t0))."""
    return 5.0 + 5.0 * (np.cos(2 * t0) - np.cos(2 * t1)) / (2.0 * (t1 - t0))


@pytest.mark.parametrize("acc_window, expected", [
    ((0.0, 1.0), 2.000),     # same window: Newton II holds exactly
    ((0.5, 1.5), 1.935),     # shifted by half
    ((0.0, 2.0), 2.417),     # doubled
])
def test_t33_window_mismatch(acc_window, expected):
    m = 2.0
    force = mean_push(0.0, 1.0) * U
    acceleration = mean_push(*acc_window) / m * U
    kind, value = exact_mass(force, acceleration)
    assert kind == "one"
    assert value == pytest.approx(expected, abs=5e-4)

