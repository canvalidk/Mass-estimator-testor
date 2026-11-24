"""What combining results in any order does and does not force on what is carried (T83; L25,
R-2026-11-24-1). Harvested 24 Nov 2026 from the exploration branch `explore/e17-l25-and-g25`
(its E17-L1 and E17-L2) and checked by Claude on main's own code.

Three things are measured. A normalised law with one fixed starting law also combines in any
order, so the strong clause of L25 does not by itself exclude it. A count of readings cannot be
carried continuously through a zero reading: main's scale-free test on scalar readings is left
alone by an exact zero and changes its verdict beside one reading a thousandth of a noise unit
in size. On vector readings the same test is steady beside tiny readings."""
import numpy as np
import pytest

from met.conditional_lr import conditional_test
from met.evidence import Curve
from met.full import isotropic

GRID = np.geomspace(0.02, 400.0, 20001)
START = 1.0 / (1.0 + GRID ** 2)                      # a half-Cauchy starting law, centre 1


def _area(values):
    return float(np.sum(0.5 * (values[1:] + values[:-1]) * np.diff(GRID)))


def _law(readings):
    raw = START * np.exp(-0.5 * Curve(readings, scale=1.0)(GRID))
    return raw / _area(raw)


def _join(first, second):
    raw = first * second / START                     # the starting law is in both; divide it out once
    return raw / _area(raw)


def _apart(first, second):
    return 0.5 * _area(np.abs(first - second))


def test_a_normalised_law_with_one_fixed_starting_law_also_combines_in_any_order():
    a = [isotropic([6.0], [2.0], 1.0, 1.0)]
    b = [isotropic([9.0], [3.0], 1.0, 1.0)]
    c = [isotropic([2.0], [1.5], 1.0, 1.0)]
    nothing = START / _area(START)                   # the all-zero experiment's law
    assert _apart(_join(_law(a), _law(b)), _law(a + b)) < 1e-12              # joined = pooled
    assert _apart(_join(_law(a), _law(b)), _join(_law(b), _law(a))) < 1e-12  # any order
    assert _apart(_join(_join(_law(a), _law(b)), _law(c)), _join(_law(a), _join(_law(b), _law(c)))) < 1e-12
    assert _apart(_join(_law(a), nothing), _law(a)) < 1e-12                  # the all-zero experiment changes nothing
    # what does fail is multiplying two normalised laws as they stand (the starting law counted
    # twice), which is the check on the exploration branch
    product = _law(a) * _law(b)
    assert _apart(product / _area(product), _law(a + b)) == pytest.approx(0.198, abs=0.005)
    # and unnormalised curves add, with the all-zero experiment's flat curve as identity (T73)
    masses = np.array([0.5, 1.5, 3.0, 8.0])
    assert Curve(a + b, scale=1.0)(masses) == pytest.approx(Curve(a, scale=1.0)(masses) + Curve(b, scale=1.0)(masses), rel=1e-12)


@pytest.mark.parametrize("which", ["five equal readings", "five different readings"])
def test_a_zero_reading_changes_nothing_and_a_tiny_one_changes_the_verdict_on_scalar_readings(which):
    if which == "five equal readings":
        series, mass = [isotropic([5.0], [1.0], 1.0, 1.0)] * 5, 1.0
    else:
        rng = np.random.default_rng(5)
        series = [isotropic([5.0 + rng.normal()], [1.0 + 0.3 * rng.normal()], 1.0, 1.0) for _ in range(5)]
        mass = 2.5

    def p(extra):
        return conditional_test(series + extra, mass, null="spherical", scale=1.0)["p_value"]

    assert p([]) == pytest.approx(2 / 32)                                    # kept at 95%
    assert p([isotropic([0.0], [0.0], 1.0, 1.0)]) == pytest.approx(2 / 32)   # an exact zero: nothing
    assert p([isotropic([1e-3], [0.0], 1.0, 1.0)]) == pytest.approx(2 / 64)  # a thousandth of a noise unit: rejected
    assert p([isotropic([1e-1], [0.0], 1.0, 1.0)]) == pytest.approx(2 / 64)
    # In exact arithmetic the halving happens for every size above zero: a sign test on n readings
    # cannot go below 2 / 2^n, and the tiny reading is a sixth reading. met.conditional_lr counts
    # differences below 1e-9 of the statistic as ties, so at 1e-6 it still returns 2/32.
    assert p([isotropic([1e-6], [0.0], 1.0, 1.0)]) == pytest.approx(2 / 32)


def test_on_vector_readings_the_scale_free_test_is_steady_beside_tiny_readings():
    series = [isotropic([5.0, 0.2], [1.0, 0.1], 1.0, 1.0), isotropic([4.8, -0.3], [1.0, -0.05], 1.0, 1.0),
              isotropic([5.2, 0.1], [1.0, 0.02], 1.0, 1.0)]
    tiny = [isotropic([1e-3, -1e-3], [1e-3, 1e-3], 1.0, 1.0)] * 20
    seen = []
    for mass in (4.6, 4.8, 5.0):
        alone = conditional_test(series, mass, null="spherical", scale=1.0, draws=1999, seed=1)["p_value"]
        beside = conditional_test(series + tiny, mass, null="spherical", scale=1.0, draws=1999, seed=1)["p_value"]
        assert abs(alone - beside) < 0.04
        seen.append(alone)
    assert min(seen) < 0.05 < max(seen)                  # the masses straddle the verdict
