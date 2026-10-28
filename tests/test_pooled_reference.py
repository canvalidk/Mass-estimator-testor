"""O14 moment identities and the separate point movement caused by scale learning."""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import roots_hermitenorm
from scipy.stats import gamma

import met.full as full
from met.pooled_reference import _integrated_moments, moments


@pytest.mark.parametrize("k", [2, 3])
def test_null_first_and_second_moments_match_independent_integrals(k):
    result = moments(k, 0)
    mean, second = _integrated_moments(k, 0.0, 96)
    assert mean == pytest.approx(result["mean"], abs=7e-4)
    assert second - mean ** 2 == pytest.approx(result["variance"], abs=3e-3)


@pytest.mark.parametrize("k", [2, 3])
def test_general_moments_increase_toward_strong_signal_chi_square(k):
    results = [moments(k, strength) for strength in (0, 1, 8, 1000)]
    means = [r["mean"] for r in results]
    assert all(a < b < k - 1 for a, b in zip(means, means[1:]))
    assert means[-1] == pytest.approx(k - 1, rel=0.005)
    assert results[-1]["variance"] == pytest.approx(2 * (k - 1), rel=0.025)
    assert all(r["mean_error_estimate"] >= 0 for r in results)


@pytest.mark.parametrize("k,lam,nodes", [(1, 0, 32), (True, 0, 32), (2, -1, 32),
                                                (2, True, 32), (2, float("inf"), 32), (2, 1, 2)])
def test_reference_rejects_inputs_outside_its_scope(k, lam, nodes):
    with pytest.raises(ValueError):
        moments(k, lam, nodes)


def _null_padding_case():
    base = [full.isotropic([2, 0], [0, 0], 1, 1),
            full.isotropic([2, 0], [0, 2], 1, 1)]
    zero = full.isotropic([0, 0], [0, 0], 1, 1)
    return base, zero


def _pooled_precision(readings):
    states = [full.reading_state(r, 1.0) for r in readings]
    return sum(s["nu"] for s in states) / sum(s["Q"] for s in states)


def test_zero_padding_changes_pooled_precision_even_with_a_cap():
    base, zero = _null_padding_case()
    before, after = _pooled_precision(base), _pooled_precision(base + [zero])
    assert (before, after) == pytest.approx((0.5, 0.75))
    assert min(1.0, before) != min(1.0, after)
    # Existing Bartlett itself ignores the zero; this failure is the new pool.
    weights = [1.0, 0.25]
    assert full.bartlett_eta(base, 1.0, weights)[0] == pytest.approx(1.0)
    assert full.bartlett_eta(base + [zero], 1.0, weights + [1.0])[0] == pytest.approx(1.0)


def test_scale_power_moves_the_point_with_the_original_cloud_frozen():
    """Hold A fixed to exclude O2's extra-cloud mechanism entirely.

    L(m)=2m²/(1+m²) is strictly increasing, so the derivative of this point
    is Cov_(A p_t)(m,L)>0. The same law's median moves too.
    """
    base, zero = _null_padding_case()
    lm = np.linspace(-18, 18, 2401)
    mass = np.exp(lm)
    size = full._stacked_mean_length(mass, base)
    evidence = 2 * mass ** 2 / (1 + mass ** 2)
    points, medians, covariances = [], [], []
    for readings in (base, base + [zero]):
        precision = _pooled_precision(readings)
        # Both original readings retain their old g; the zero contributes nothing.
        weights = [1.0, 0.25] + [1.0] * (len(readings) - len(base))
        lp = full._log_law(readings, 1.0, weights, precision, lm)
        (median,), probability = full._quantiles(lm, lp, [0.5])
        w = size * probability
        w /= w.sum()
        point = np.sum(mass * w)
        covariance = np.sum(w * (mass - point) * (evidence - np.sum(w * evidence)))
        points.append(point)
        medians.append(median)
        covariances.append(covariance)
    assert points[1] > points[0] + 0.05
    assert medians[1] > medians[0] + 0.05
    assert min(covariances) > 0


def test_independent_compact_integrals_confirm_scale_point_movement():
    """An explicit fixed A=1/sqrt(1+m²), checked by two variable changes."""
    points = []
    for t in (0.5, 0.75):
        compact = math.exp(2 * t) * quad(lambda u: math.exp(-2 * t * u * u), 0, 1)[0]
        compact /= quad(lambda u: math.exp(2 * t * u * u), 0, 1)[0]

        def angular(phi, numerator):
            m = math.tan(phi)
            return (math.sin(phi) if numerator else math.cos(phi)) * math.exp(2 * t * m * m / (1 + m * m))

        direct = quad(lambda phi: angular(phi, True), 0, math.pi / 2)[0]
        direct /= quad(lambda phi: angular(phi, False), 0, math.pi / 2)[0]
        assert compact == pytest.approx(direct, abs=1e-12)
        points.append(compact)
    assert points[1] > points[0]


@pytest.mark.parametrize("n", [3, 5, 20])
def test_three_dimensional_null_pooled_precision_has_inverse_gamma_moments(n):
    """At k=3 null Q is exactly Exp(1), hence sum Q is Gamma(N,1).

    The reciprocal is not unbiased for the true precision 1: even at N->infinity
    its mean tends to 2. Moments fail to exist at N=1 (mean) and N<=2 (variance).
    """
    first = quad(lambda z: 2 * n / z * gamma.pdf(z, n), 0, np.inf, epsabs=1e-9)[0]
    second = quad(lambda z: (2 * n / z) ** 2 * gamma.pdf(z, n), 0, np.inf, epsabs=1e-8)[0]
    assert first == pytest.approx(2 * n / (n - 1), rel=1e-9)
    assert second - first * first == pytest.approx(4 * n * n / ((n - 1) ** 2 * (n - 2)), rel=1e-9)


@pytest.mark.parametrize("k", [2, 3])
@pytest.mark.parametrize("lam", [0.0, 2.0, 8.0])
def test_unweighted_eta_reference_uses_signal_and_score_noise_moments(k, lam):
    """At the true line, curvature has mean lambda, score² has mean lambda+k.

    Hermite integration is exact for these Gaussian polynomial moments. Cross-
    coordinate score products have zero mean, so their second moments add.
    """
    nodes, weights = roots_hermitenorm(5)
    weights /= math.sqrt(2 * math.pi)
    curvature, score_second = 0.0, 0.0
    for axis in range(k):
        parallel = nodes[:, None] + (math.sqrt(lam) if axis == 0 else 0.0)
        perpendicular = nodes[None, :]
        joint = weights[:, None] * weights[None, :]
        curvature += np.sum(joint * (parallel ** 2 - perpendicular ** 2))
        score_second += np.sum(joint * parallel ** 2 * perpendicular ** 2)
    assert curvature == pytest.approx(lam, abs=1e-12)
    assert score_second == pytest.approx(lam + k, abs=1e-12)
