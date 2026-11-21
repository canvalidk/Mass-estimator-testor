"""Check the positive-mass common fit of bounded per-reading residuals."""
import math

import numpy as np
import pytest
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize_scalar

from explore.e12_residual_join import residual_log_evidence, residual_minimum
from met.full import Reading, _h_terms


ANGLES = np.linspace(0., math.pi / 2, 101)


def row_grams(rows):
    rows = np.asarray(rows, dtype=float)
    return rows[:, :, None] * rows[:, None, :]


def response(cost, dimension=1):
    scale = dimension / (1 - 2**(-dimension / 2))
    return -scale * np.expm1(-.5 * np.asarray(cost))


def test_distinct_finite_laws_for_two_positive_incompatible_readings():
    grams = row_grams([(2, 1), (1, 2)])
    # The midline minimizes the bounded cost here as well as the raw cost.
    expected = 2 * response(.5)
    minimum = residual_minimum(grams)
    assert minimum == pytest.approx(expected, abs=1e-8)
    excess = minimum - response(1.)
    assert excess > .1
    penalty = 36/11**3 * excess/(1+minimum)
    angles = np.arctan([1., .5])
    contrasts = [np.subtract(*residual_log_evidence(grams, angles, tau=t)) for t in (0., 64.)]
    difference = (64/9)/20 * (1 - 1/(1 + 64*penalty))
    assert contrasts[0] - contrasts[1] == pytest.approx(difference, abs=1e-9)


@pytest.mark.parametrize("gram,dimension,least", [
    ([[4., 2.], [2., 1.]], 1, 0.),
    ([[4., -2.], [-2., 1.]], 1, 1.),
    ([[7., 2.], [2., 4.]], 2, 3.),
    ([[5., -2.], [-2., 4.]], 2, 4.),
    ([[2., 0.], [0., 2.]], 2, 2.),
])
def test_singleton_uses_exact_constrained_residual_and_has_no_tau_dependence(gram, dimension, least):
    grams = np.array([gram])
    assert residual_minimum(grams, dimension=dimension) == pytest.approx(response(least, dimension), abs=1e-8)
    base = residual_log_evidence(grams, ANGLES, dimension=dimension)
    perturbed = residual_log_evidence(grams, ANGLES, dimension=dimension, tau=1e4)
    np.testing.assert_allclose(base, perturbed, atol=1e-8, rtol=0.)


def test_block_concatenation_preserves_readings_and_zero_and_channel_symmetries():
    first, second, third = [row_grams(rows) for rows in (
        [(3, 1), (-1, 4)], [(2, 1)], [(1, 3), (1, 2)],
    )]
    zeros = np.zeros((2, 2, 2))
    left = np.concatenate((np.concatenate((first, second)), third))
    right = np.concatenate((third, np.concatenate((zeros, second, first))))
    assert residual_minimum(left) == pytest.approx(residual_minimum(right), abs=1e-8)
    np.testing.assert_allclose(residual_log_evidence(left, ANGLES, tau=64),
                               residual_log_evidence(right, ANGLES, tau=64), atol=1e-8)
    swapped = left[:, ::-1, ::-1]
    np.testing.assert_allclose(residual_log_evidence(left, ANGLES, tau=64),
                               residual_log_evidence(swapped, math.pi/2 - ANGLES, tau=64), atol=1e-8)
    assert residual_minimum(zeros) == 0.
    np.testing.assert_allclose(residual_log_evidence(zeros, ANGLES), 0., atol=0.)


def test_narrow_competing_mass_fits_have_a_global_minimum_between_individual_lines():
    # Four close but distinct lines and three copies of a distant competing line.
    # The objective is almost flat between the two clusters.
    magnitude, centre, other = 6000., .74312371, 1.235
    offsets = np.array([-.48, -.26, .44, .59])
    axes = np.r_[centre + offsets/magnitude, [other]*3]
    grams = row_grams(magnitude * np.column_stack((np.sin(axes), np.cos(axes))))

    def direct_cost(offset):
        theta = centre + offset/magnitude
        return response((magnitude * np.sin(axes-theta))**2).sum()

    best = minimize_scalar(direct_cost, bounds=(-.5, .6), method="bounded", options={"xatol": 1e-13})
    assert best.success
    assert min(direct_cost(x) for x in offsets) - best.fun > .1
    assert direct_cost((other-centre)*magnitude) - best.fun > 1.
    assert residual_minimum(grams) == pytest.approx(best.fun, rel=0., abs=1e-7)


def test_fixed_positive_conflict_counts_readings_outside_the_best_supported_exact_mass():
    # Three readings agree on 2, two on 1/3, one on 4. A measured zero fits each.
    rows = [(2, 1), (4, 2), (6, 3), (1, 3), (2, 6), (4, 1), (0, 0)]
    grams = row_grams(rows)
    limit = 3 / (1 - math.sqrt(.5))
    assert 0 < residual_minimum(grams) < limit - 1.
    assert residual_minimum(10000 * grams) == pytest.approx(limit, abs=1e-7)
    assert residual_minimum(10000 * row_grams([(2, 1), (-4, -2), (0, 0)])) == pytest.approx(0., abs=1e-8)


def test_negative_slopes_have_residual_at_every_allowed_mass_including_both_axes():
    # Each absolute residual exceeds cos(theta) + sin(theta) away from pi/2.
    grams = row_grams([(-4, 1), (2, -1)])
    assert residual_minimum(grams) == pytest.approx(2 * response(1.), abs=1e-8)
    assert residual_minimum(10000 * grams) == pytest.approx(2 / (1 - math.sqrt(.5)), abs=1e-8)


def test_scalar_outer_product_keeps_its_exact_fit_at_large_signal():
    # Rounding a scalar Gram at this scale can leave an order-one positive
    # smaller eigenvalue even though the declared reading lies on one line.
    grams = 1e16 * row_grams([(1.54, 1.79)])
    assert residual_minimum(grams) == pytest.approx(0., abs=1e-8)


@pytest.mark.parametrize("tau,expected", [
    (0., [.307038, 1.528611, 15.258568]),
    (64., [.270048, 1.510203, 16.247667]),
])
def test_joint_mass_median_and_central_interval_under_uniform_angular_prior(tau, expected):
    grams = row_grams([(3, 1), (1, 2)])
    angles = np.linspace(0, math.pi/2, 100001)
    log_weights = residual_log_evidence(grams, angles, tau=tau)
    cdf = cumulative_trapezoid(np.exp(log_weights - log_weights.max()), angles, initial=0.)
    mass_quantiles = np.tan(np.interp([.025, .5, .975], cdf / cdf[-1], angles))
    np.testing.assert_allclose(mass_quantiles, expected, atol=1e-6, rtol=0.)


@pytest.mark.parametrize("k", [1, 2])
def test_full_t36_residual_threshold_is_a_degree_2k_polynomial_inequality(k):
    # Work in a proper k-plane with correlated force/acceleration errors.
    rng = np.random.default_rng(1100 + k)
    d = k + 1
    plane, _ = np.linalg.qr(rng.normal(size=(d, k)))
    matrix = rng.normal(size=(2*d, 2*d))
    covariance = matrix @ matrix.T + np.eye(2*d)
    observation = rng.normal(size=2*d)
    reading = Reading(observation[:d], observation[d:], covariance, plane)
    pre = reading.pre()

    # The T36 residual excludes the fixed component outside the declared span.
    root = np.linalg.cholesky(covariance)
    zero = np.zeros_like(plane)
    span = np.linalg.solve(root, np.block([[plane, zero], [zero, plane]]))
    basis, _ = np.linalg.qr(span)
    projected = basis @ (basis.T @ np.linalg.solve(root, observation))
    assert projected @ projected == pytest.approx(pre["energy"], abs=1e-12)

    masses = np.array([0., .04, .22, .61, .99, 1.52, 2.36, 2.88, 3.])
    costs = pre["energy"] - _h_terms(pre, 1., np.arctan(masses))
    threshold = (max(costs) + min(costs)) / 2

    def determinants(mass):
        design = mass*span[:, :k] + span[:, k:]
        augmented = np.column_stack((design, projected))
        bordered = augmented.T @ augmented
        bordered[-1, -1] -= threshold
        return np.linalg.det(design.T @ design), np.linalg.det(bordered)

    denominator, numerator = np.array([determinants(m) for m in masses]).T
    assert np.all(denominator > 0)
    np.testing.assert_allclose(numerator, denominator * (costs-threshold), atol=1e-12)
    np.testing.assert_array_equal(numerator < 0, costs < threshold)

    # Reconstruct with 2k+1 other points and check the degree bound off-grid.
    nodes = 1.5 * (1 + np.polynomial.chebyshev.chebpts1(2*k + 1))
    values = [determinants(m)[1] for m in nodes]
    polynomial = np.polynomial.Chebyshev.fit(nodes, values, 2*k)
    np.testing.assert_allclose(polynomial(masses), numerator, atol=1e-12)
