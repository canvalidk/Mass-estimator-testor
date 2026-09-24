"""Single-reading readouts: exact null law, known values, symmetries, convergence."""

import math

import numpy as np
import pytest

from met.analyst import estimate

ALL = ["ratio_of_means", "median", "geometric", "reciprocal_root", "log_sd", "interval_95"]


@pytest.mark.parametrize("d", [1, 2, 3])
@pytest.mark.parametrize("sf,sa", [(1.0, 1.0), (2.0, 0.5), (0.3, 3.0)])
def test_null_law_is_half_cauchy_at_the_instrument_ratio(d, sf, sa):
    s = sf / sa
    r = estimate(np.zeros(d), np.zeros(d), sf, sa, ALL)
    assert r["ratio_of_means"] == pytest.approx(s, rel=1e-12)
    assert r["median"] == pytest.approx(s, rel=1e-10)
    assert r["geometric"] == pytest.approx(s, rel=1e-10)
    assert r["reciprocal_root"] == pytest.approx(s, rel=1e-10)
    assert r["log_sd"] == pytest.approx(math.pi / 2, rel=1e-12)
    expected = [s * math.tan(0.025 * math.pi / 2), s * math.tan(0.975 * math.pi / 2)]
    assert r["interval_95"] == pytest.approx(expected, rel=1e-7)


def test_equation_readme_example():
    """Equation README (21 Sept): 4.993754 kg, 95% [4.344146, 5.740363], log SD 0.071035.

    The README's interval and log SD carry that implementation's stated rtol 2e-4.
    Our values converge to [4.3442404, 5.7402376] and 0.0710234 (see convergence test).
    """
    r = estimate([10, 0, 0], [2, 0.1, 0], 0.5, 0.1, ["ratio_of_means", "interval_95", "log_sd"])
    assert r["ratio_of_means"] == pytest.approx(4.993754, abs=5e-7)
    assert r["interval_95"] == pytest.approx([4.344146, 5.740363], rel=2e-4)
    assert r["log_sd"] == pytest.approx(0.071035, rel=2e-4)
    assert r["interval_95"] == pytest.approx([4.3442404, 5.7402376], rel=2e-7)


def test_contribution_11_single_trials():
    t1 = estimate([2, 0.2, 0], [1, 0, 0], 0.5, 0.4)
    t2 = estimate([5, 0, 0.1], [2.1, 0.2, 0], 1, 0.6)
    assert t1["ratio_of_means"] == pytest.approx(1.993026469631, rel=1e-11)
    assert t2["ratio_of_means"] == pytest.approx(2.373589604907, rel=1e-11)


@pytest.mark.parametrize("d", [1, 2, 3])
def test_reciprocal_symmetry(d):
    rng = np.random.default_rng(10 + d)
    for _ in range(4):
        f, a = rng.normal(size=d) * 3, rng.normal(size=d)
        sf, sa = rng.uniform(0.3, 2, size=2)
        r = estimate(f, a, sf, sa, ALL)
        w = estimate(a, f, sa, sf, ALL)
        assert r["ratio_of_means"] * w["ratio_of_means"] == pytest.approx(1, rel=1e-9)
        assert r["median"] * w["median"] == pytest.approx(1, rel=1e-9)
        assert r["geometric"] * w["geometric"] == pytest.approx(1, rel=1e-9)
        assert r["reciprocal_root"] * w["reciprocal_root"] == pytest.approx(1, rel=1e-9)
        assert r["log_sd"] == pytest.approx(w["log_sd"], rel=1e-9)
        assert r["interval_95"] * w["interval_95"][::-1] == pytest.approx([1, 1], rel=1e-7)


@pytest.mark.parametrize("d", [1, 2, 3])
def test_force_unit_change(d):
    rng = np.random.default_rng(20 + d)
    f, a = rng.normal(size=d) * 3, rng.normal(size=d)
    r = estimate(f, a, 0.7, 1.3, ALL)
    c = 1000.0
    w = estimate(f * c, a, 0.7 * c, 1.3, ALL)
    for key in ("ratio_of_means", "median", "geometric", "reciprocal_root"):
        assert w[key] == pytest.approx(c * r[key], rel=1e-9)
    assert w["interval_95"] == pytest.approx(c * r["interval_95"], rel=1e-7)


@pytest.mark.parametrize("d", [2, 3])
def test_rotation_invariance(d):
    rng = np.random.default_rng(30 + d)
    q, _ = np.linalg.qr(rng.normal(size=(d, d)))
    f, a = rng.normal(size=d) * 2, rng.normal(size=d)
    r = estimate(f, a, 0.8, 0.6, ALL)
    w = estimate(q @ f, q @ a, 0.8, 0.6, ALL)
    for key in ALL:
        assert np.allclose(w[key], r[key], rtol=1e-9)


def test_direction_sign_in_one_dimension():
    """In 1D both directions are allowed: flipping both signs changes nothing."""
    r = estimate([3.0], [-1.0], 1, 1, ALL)
    w = estimate([-3.0], [1.0], 1, 1, ALL)
    for key in ALL:
        assert np.allclose(w[key], r[key], rtol=1e-12)


def test_grid_convergence():
    for f, a in (([10, 0, 0], [2, 0.1, 0]), ([0.3, -0.2, 0], [0.1, 0.4, 0]), ([40, 1, 0], [8, 0, 0.2])):
        coarse = estimate(f, a, 1.0, 0.5, ALL)
        fine = estimate(f, a, 1.0, 0.5, ALL, numerics={"grid_points": 16001, "cutoff": 50.0})
        for key in ALL:
            assert np.allclose(coarse[key], fine[key], rtol=2e-7), key
