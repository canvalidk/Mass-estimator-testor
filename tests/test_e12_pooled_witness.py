"""Finite identities and boundary cases for the E12 scalar pooled proof."""
import math

import numpy as np
import pytest
from numpy.polynomial.hermite import hermgauss

from explore.e12_scalar_witness import effective_count, pooled_witness
from met.gaussian_information import gaussian_angle_quantile


@pytest.mark.parametrize("rate", [1.0, 2.0])
@pytest.mark.parametrize("radius", [0.0, 0.75, 3.0])
def test_gaussian_count_moment_by_hermite_quadrature(rate, radius):
    nodes, weights = hermgauss(80)
    x, y = np.meshgrid(math.sqrt(2) * nodes, math.sqrt(2) * nodes)
    t = (x + radius) ** 2 + y ** 2
    h = np.array([effective_count([value], rate=rate) for value in t.ravel()]).reshape(t.shape)
    actual = np.sum(h * np.outer(weights, weights)) / math.pi
    expected = 1 + (math.exp(-rate * radius ** 2 / (2 * (1 + rate)))
                    - math.exp(-rate * radius ** 2 / (1 + 2 * rate))) / rate
    assert actual == pytest.approx(expected, abs=4e-12)
    assert expected >= 1.0


def test_exact_zero_padding_and_arbitrarily_small_rows():
    rows = np.array([[3., 3.], [2.5, -2.5]]) / math.sqrt(2)
    reference = pooled_witness(rows, alpha=0.125)
    assert pooled_witness(np.vstack([[0., 0.], rows, [0., 0.]]), alpha=0.125) == reference
    near = pooled_witness(np.vstack([rows, [1e-7, 2e-7]]), alpha=0.125)
    assert near.precision == pytest.approx(reference.precision, abs=3e-14)
    assert effective_count([0.]) == 0.
    assert effective_count([1e-34]) == pytest.approx(2e-34, rel=1e-14, abs=0.)


def test_swap_and_independent_row_signs_preserve_concentration():
    rows = np.array([[4., 1.], [1., 2.], [0., -1.]])
    base = pooled_witness(rows, alpha=0.125)
    swapped = pooled_witness(rows[::-1, ::-1] * [[-1], [1], [-1]], alpha=0.125)
    assert swapped.angle == pytest.approx(math.pi/2 - base.angle)
    assert swapped.positive_residual == pytest.approx(base.positive_residual)
    assert swapped.precision == pytest.approx(base.precision)


def test_positive_arc_moves_the_best_negative_line_to_an_axis():
    rows = np.array([[1., -3.], [2., -6.]])
    result = pooled_witness(rows, alpha=0.125)
    assert result.positive_residual == pytest.approx(5.)
    assert result.largest == pytest.approx(50.)
    assert result.angle < 0.
    assert result.precision > 0.


@pytest.mark.parametrize("alpha", [0.0625, 0.125, 1.])
def test_incompatible_shrinking_ray_converges_to_stated_nonzero_tilt(alpha):
    rows = np.column_stack(([.2, 1., 1.4], np.ones(3)))
    # The shape has lambda_1=5.6, eigengap 5.2 and least positive residual 0.4.
    limit = 3 / .4 * (5.6 / .4) ** alpha * 5.2 ** 2 / 5.6
    for scale in [50., 100.]:
        fit = pooled_witness(scale * rows, alpha=alpha)
        assert fit.precision == pytest.approx(limit, rel=2e-13)
        assert 0. < fit.discount < 1.


def test_isotropic_and_axis_exceptions_and_compatible_interior_line():
    assert pooled_witness([], alpha=1.).precision == 0.
    assert pooled_witness(np.zeros((3, 2)), alpha=1.).count == 0.
    isotropic = pooled_witness(100 * np.eye(2), alpha=1.)
    assert isotropic.gap == isotropic.precision == 0.
    assert isotropic.positive_residual > 0.
    axis = pooled_witness([[0., 100.]], alpha=1.)
    assert axis.positive_residual == 0. and axis.discount == 1.
    assert axis.precision == pytest.approx(10000.)
    aligned = pooled_witness([[20., 20.], [10., 10.]], alpha=1.)
    assert aligned.positive_residual == pytest.approx(0., abs=1e-22)
    assert aligned.precision == pytest.approx(1000.)


def test_count_convention_changes_the_finite_law():
    rows = np.array([[3., 3.], [2.5, -2.5]]) / math.sqrt(2)
    first = pooled_witness(rows, alpha=0.125, count_rate=1.)
    second = pooled_witness(rows, alpha=0.125, count_rate=2.)
    assert first.precision == pytest.approx(0.2960530326)
    assert second.precision == pytest.approx(0.2818607553)
    assert first.angle == pytest.approx(second.angle)


@pytest.mark.parametrize("alpha,far_discount,far_half_width", [
    (1., 1., .13944), (.125, .12840736, .38914), (.0625, .07210028, .51934),
])
def test_two_ruled_scenes_at_the_declared_instrument(alpha, far_discount, far_half_width):
    def report(forces):
        rows = np.column_stack((np.array(forces) / .1, [500., 500.]))
        fit = pooled_witness(rows, alpha=alpha)
        masses = [50 * math.tan(gaussian_angle_quantile(p, fitted_angle=fit.angle,
                                                        precision=fit.precision))
                  for p in [.025, .975]]
        return fit.discount, (masses[1] - masses[0]) / 2

    assert report([5.5, 5.6]) == pytest.approx((1., .13944), abs=1e-5)
    assert report([5., 6.]) == pytest.approx((far_discount, far_half_width), abs=1e-5)


def test_far_discount_switches_off_in_the_known_acceleration_limit():
    # Unbounded positive alpha makes the ratio lambda_max/Q arbitrarily large.
    rows = np.column_stack(([50., 60.], [1e15, 1e15]))
    fit = pooled_witness(rows, alpha=.125)
    assert fit.positive_residual == pytest.approx(50.)
    assert fit.discount == 1.


def test_two_row_gap_inflation_coverage_in_original_independent_coordinates():
    # Quadrature conditional on W1=x integrates over W2 in the original
    # coordinates. The formula computes the same event after rotating to Z,G.
    from scipy.integrate import quad
    from scipy.special import ndtr, ndtri

    p, cutoff, multiplier = .95, 2., 1.5
    z = ndtri((1 + p) / 2)
    root_two = math.sqrt(2)

    def joint_quiet(x, width):
        lo = max(x - root_two * cutoff, -x - root_two * width)
        hi = min(x + root_two * cutoff, -x + root_two * width)
        return max(0., ndtr(hi) - ndtr(lo))

    def integrand(x):
        expanded = ndtr(-x + root_two * multiplier * z) - ndtr(-x - root_two * multiplier * z)
        conditional = expanded + joint_quiet(x, z) - joint_quiet(x, multiplier * z)
        return conditional * math.exp(-x*x/2) / math.sqrt(2*math.pi)

    corners = [sign * (cutoff + relative * w) / root_two
               for sign in [-1, 1] for relative in [-1, 1] for w in [z, multiplier*z]]
    direct = quad(integrand, -12., 12., points=sorted(corners), epsabs=1e-12)[0]
    rotated = p + 2 * ndtr(-cutoff) * (2 * ndtr(multiplier * z) - 1 - p)
    assert direct == pytest.approx(rotated, abs=1e-12)
    assert direct > p + .002


def test_fixed_n_amplification_switches_off_the_extra_discount():
    # A fixed realized, large residual in the transverse direction.
    u, v = np.array([1., 1.]) / math.sqrt(2), np.array([1., -1.]) / math.sqrt(2)
    rows = 1e8 * u + np.array([3., -3.])[:, None] * v
    fit = pooled_witness(rows, alpha=.125)
    assert fit.positive_residual == pytest.approx(18., abs=1e-6)
    assert fit.discount == 1. and fit.precision == pytest.approx(2e16)


@pytest.mark.parametrize("rows,alpha,rate", [
    ([[1., 2.]], 0., 1.), ([[1., 2.]], 1., float("nan")),
    ([[float("inf"), 2.]], 1., 1.), ([[1., 2., 3.]], 1., 1.),
])
def test_invalid_inputs(rows, alpha, rate):
    with pytest.raises(ValueError):
        pooled_witness(rows, alpha=alpha, count_rate=rate)
