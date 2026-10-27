"""O17: narrow scalar laws can coexist with a nonzero tensor/scalar mismatch.

These are deterministic permitted-input families, not a claim that Gaussian noise
is exactly zero in a typical generated experiment. They test the implementation
against derived limits; they do not invent an adequacy requirement or a no-mass gate.
"""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import chi2, norm

from met.analyst import (bartlett_weight, estimate, joint_log_density,
                         misfit_weights, prior_log_density)
from met.law import ReadingSet


GOLDEN = (1 + math.sqrt(5)) / 2
THETA = math.atan(GOLDEN)
Q0 = (7 - 3 * math.sqrt(5)) / 4
CART = {"reference": "cartesian2", "nuisance": "radius"}
OURS29 = {"rule": "likelihood_product",
          "prior": [{"uniform_angle": {"center": 1.0}}],
          "misfit": "kink", "calibrate": "bartlett"}
READOUTS = ["ratio_of_means_rss", "median", "interval_95"]


def _family(name, n, sd):
    if name == "oblique":
        a = np.tile(np.array([1., 1., 0.]) / math.sqrt(2), (n, 1))
    else:
        assert name == "axes" and n % 2 == 0
        a = np.tile(np.eye(3)[:2], (n // 2, 1))
    f = a @ np.diag([1., 2., 4.])
    readings = ReadingSet(f[None], a[None], np.full((1, n), sd), np.full((1, n), sd))
    return f, a, readings


def _curvature_per_reading(name):
    return 45 + 21 * math.sqrt(5) if name == "oblique" else 225 / 4


@pytest.mark.parametrize("name", ["oblique", "axes"])
def test_tensor_family_has_no_common_exact_scalar_despite_positive_pairings(name):
    f, a, _ = _family(name, 8, 0.03)
    m_least_squares = np.sum(f * a) / np.sum(a * a)
    assert np.all(np.sum(f * a, axis=1) > 0)
    assert np.linalg.norm(f - m_least_squares * a) > 0.1
    # Independent stacked SVD profiles the unrestricted common scalar line.
    singular = np.linalg.svd(np.column_stack((f.ravel(), a.ravel())), compute_uv=False)
    assert singular[-1] ** 2 == pytest.approx(8 * Q0, rel=1e-12)


@pytest.mark.parametrize("name", ["oblique", "axes"])
def test_individual_misfits_and_bartlett_power_distinguish_two_obstructions(name):
    n, sd = 8, 0.03
    f, a, readings = _family(name, n, sd)
    weights = misfit_weights(readings, OURS29)
    power = bartlett_weight(readings, OURS29)[0]
    if name == "oblique":
        np.testing.assert_allclose(weights, 2 * sd ** 2 / Q0, rtol=2e-13)
        assert power == pytest.approx(1.0)
    else:
        np.testing.assert_array_equal(weights, np.ones((1, n)))
        assert power == pytest.approx((15 * math.sqrt(5) / 2) * sd ** 2, rel=1e-12)
    total_cost = np.linalg.svd(np.column_stack((f.ravel(), a.ravel())), compute_uv=False)[-1] ** 2 / sd ** 2
    assert total_cost == pytest.approx(n * Q0 / sd ** 2, rel=1e-12)
    assert total_cost > 100


@pytest.mark.parametrize("name", ["oblique", "axes"])
def test_quoted_angle_kernel_freezes_while_absolute_cost_grows(name):
    u = np.linspace(-2, 3, 401)[None]
    theta = np.arctan(np.exp(u[0]))
    n, curves, costs = 8, [], []
    for sd in (0.03, 0.01):
        f, a, readings = _family(name, n, sd)
        lp = joint_log_density(readings, u, CART, OURS29, alpha=False)[0][0]
        evidence = lp - prior_log_density(OURS29["prior"], u, readings)[0]
        evidence -= evidence[0]
        expected = n * _curvature_per_reading(name) / 4 * (
            np.cos(2 * (theta - THETA)) - np.cos(2 * (theta[0] - THETA)))
        np.testing.assert_allclose(evidence, expected, rtol=2e-11, atol=2e-8)
        curves.append(lp - lp[0])
        costs.append(np.linalg.svd(np.column_stack((f.ravel(), a.ravel())), compute_uv=False)[-1] ** 2 / sd ** 2)
    np.testing.assert_allclose(curves[0], curves[1], rtol=2e-11, atol=2e-8)
    assert costs[1] / costs[0] == pytest.approx(9, rel=1e-12)


@pytest.fixture(scope="module", params=["oblique", "axes"])
def family_summaries(request):
    name = request.param
    results = []
    for n in (8, 32, 128):
        f, a, _ = _family(name, n, 0.03)
        out = estimate(f, a, 0.03, 0.03, READOUTS, law=CART, combine=OURS29)
        results.append((n, out))
    return name, results


def test_scalar_law_concentrates_at_a_pseudomass_despite_the_exact_obstruction(family_summaries):
    name, results = family_summaries
    scaled_widths = []
    for n, out in results:
        assert not out["unresolved"]
        assert out["median"] == pytest.approx(GOLDEN, rel=2e-4)
        assert out["ratio_of_means_rss"] == pytest.approx(GOLDEN, rel=2e-3)
        lo, hi = out["interval_95"]
        assert lo < GOLDEN < hi
        scaled_widths.append(math.sqrt(n) * math.log(hi / lo))
    expected = 2 * norm.ppf(0.975) * math.sqrt(5 / _curvature_per_reading(name))
    assert scaled_widths[-1] == pytest.approx(expected, rel=0.02)
    assert scaled_widths[-1] / scaled_widths[-2] == pytest.approx(1, rel=0.02)


def test_independent_compact_angle_integral_checks_the_reported_interval(family_summaries):
    name, results = family_summaries
    n, out = results[1]
    coefficient = n * _curvature_per_reading(name) / 4

    def kernel(theta):
        return math.exp(coefficient * (math.cos(2 * (theta - THETA)) - 1))

    area = quad(kernel, 0, math.pi / 2, points=[THETA], epsabs=1e-12, epsrel=1e-11)[0]
    for mass, expected in zip(out["interval_95"], (0.025, 0.975)):
        bound = math.atan(mass)
        points = [THETA] if THETA < bound else []
        cdf = quad(kernel, 0, bound, points=points, epsabs=1e-12, epsrel=1e-11)[0] / area
        assert cdf == pytest.approx(expected, abs=3e-5)


def test_common_profile_cost_is_bounded_by_the_cost_at_any_true_scalar():
    """The conservative chi-square diagnostic uses 3N, with no fitted df claim.

    Under a true scalar m0, K(m0) is chi-square(3N). Since Cmin<=K(m0),
    sf(Cmin)>=sf(K(m0)), giving a superuniform diagnostic under that model.
    This fixed-data check verifies the ordering, not a fitted-data p-value law.
    """
    m0, sd = 1.7, 0.05
    latent = np.vstack((np.eye(3), [1., -1., 1.]))
    a = latent + np.array([[.02, -.03, .01], [.01, .04, -.02],
                           [-.04, .03, .01], [.02, .01, -.03]])
    f = m0 * latent + np.array([[.1, -.04, .02], [.03, -.1, .06],
                                [-.02, .1, -.03], [.05, -.03, .08]])
    assert np.sum(f * a) > 0  # The projective minimum lies at positive slope here.
    cmin = np.linalg.svd(np.column_stack((f.ravel(), a.ravel())), compute_uv=False)[-1] ** 2 / sd ** 2
    fixed_cost = np.sum((f - m0 * a) ** 2) / (sd ** 2 * (1 + m0 ** 2))
    assert cmin <= fixed_cost
    assert chi2.sf(cmin, 3 * len(a)) >= chi2.sf(fixed_cost, 3 * len(a))
