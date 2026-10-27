"""O13's 27 Oct candidate, kept local: production ours29 is unchanged.

Exact Gaussian moments and a finite lower-bound counterexample judge the proposed
numerator. Other checks separate that failure from algebraic properties of the
positive-cost weight. No random study or fitted true errors are used.
"""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, norm

import met.full as full
from met.misalignment import isotropic_minima


NULL_MOMENTS = {
    1: (0.0, 0.5 - 1 / math.pi, 0.5 - 1 / math.pi),
    2: (2 - math.pi / 2, math.pi / 4 - 0.5, 1.5 - math.pi / 4),
    3: (1.0, 1 - 2 / math.pi, 2 - 2 / math.pi),
}


def _kappa(k):
    return 0.5 - 1 / math.pi if k == 1 else k - 1


def _weight(reading):
    qplus = isotropic_minima(reading)["positive_min_cost"]
    return 1.0 if qplus == 0 else min(1.0, _kappa(reading.k) / qplus)


def _candidate_law(readings, c, lm):
    weights = [_weight(r) for r in readings]
    eta = full.bartlett_eta(readings, c, weights)[0] if readings else 1.0
    return full._log_law(readings, c, weights, eta, lm), eta


@pytest.mark.parametrize("k", [1, 2, 3])
def test_null_expectations_from_independent_radial_integrals(k):
    """Central eigen-shape integration and E min(A,B)=integral P(A>t)^2 dt.

    T=trace has mean 2k. For k>1, v=gap/T=sin(t) has density
    (k-1) sin(t) cos(t)^(k-2) dt; k=1 is rank one, Q=0.
    Sign symmetry makes E Qplus half E Q plus half E min(A,B).
    """
    q = 0.0 if k == 1 else quad(
        lambda t: k * (1 - math.sin(t)) * (k - 1)
        * math.sin(t) * math.cos(t) ** (k - 2), 0, math.pi / 2,
        epsabs=1e-11, epsrel=1e-11)[0]
    minimum = quad(lambda t: chi2.sf(t, k) ** 2, 0, np.inf,
                   epsabs=1e-10, epsrel=1e-10)[0]
    qplus = (q + minimum) / 2
    np.testing.assert_allclose((q, qplus - q, qplus), NULL_MOMENTS[k],
                               rtol=2e-10, atol=2e-11)


@pytest.mark.parametrize("k", [1, 2, 3])
def test_finite_unequal_channel_counterexample_exceeds_proposed_numerator(k):
    """O13 proves this lower bound for true means (32 e1, .5 e1), mass 64.

    This is a bound at finite positive channel strengths, not just a boundary
    experiment or a Monte Carlo discrepancy. The Gaussian negative-part moment
    is also checked by direct integration independently of its closed form.
    """
    r, b = 32.0, 0.5
    h = (1 + b * b) * norm.cdf(-b) - b * norm.pdf(b)
    direct = quad(lambda z: z * z * norm.pdf(z, loc=b), -np.inf, 0,
                  epsabs=1e-12)[0]
    assert h == pytest.approx(direct, abs=1e-12)
    eb, eb2 = k + b * b, (k + b * b) ** 2 + 2 * (k + 2 * b * b)
    lower = k - 1 + h - 4 * eb2 / r ** 2 - eb * norm.cdf(-r / 2)
    assert lower > _kappa(k) + 0.01


@pytest.mark.parametrize("k", [1, 2, 3])
def test_pointwise_bound_against_independent_latent_fits(k):
    """Validate 0 <= M-Qplus <= B^2/A, the finite-bound proof's geometry.

    LS fits [sin(phi)I; cos(phi)I], without using the Gram eigenvalue or
    isotropic_minima to obtain the oracle cost. Endpoints are infima on m>0.
    """
    for size in (0.2, 2.0, 32.0):
        for cosine in (-1.0, -0.4, 0.0, 0.6, 1.0):
            if k == 1 and abs(cosine) != 1:
                continue
            x, y = np.zeros(k), np.zeros(k)
            x[0], y[0] = size, cosine
            if k > 1:
                y[1] = math.sqrt(1 - cosine ** 2)
            observed = np.r_[x, y]

            def cost(phi):
                design = np.vstack((math.sin(phi) * np.eye(k),
                                    math.cos(phi) * np.eye(k)))
                fit = np.linalg.lstsq(design, observed, rcond=None)[0]
                error = observed - design @ fit
                return float(error @ error)

            interior = minimize_scalar(cost, bounds=(0, math.pi / 2), method="bounded",
                                       options={"xatol": 1e-13})
            qplus = min(cost(0), cost(math.pi / 2), interior.fun)
            a, b = x @ x, y @ y
            m = b - max(float(x @ y) / math.sqrt(a), 0) ** 2
            assert -1e-10 <= m - qplus <= b * b / a + 1e-10
            actual = isotropic_minima(full.isotropic(x, y, 1, 1))["positive_min_cost"]
            assert actual == pytest.approx(qplus, abs=1e-8)


@pytest.mark.parametrize("k", [2, 3])
def test_supremum_numerator_cannot_preserve_all_acute_weights(k):
    """O13 proves k-.5 <= S_k <= k for the supremum over all positive means.

    An acute pair with Q=2k keeps both caps inactive. Every numerator allowed
    by those proven bounds gives a strictly larger weight than (k-1)/Q.
    """
    x, y = np.zeros(k), np.zeros(k)
    x[0] = math.sqrt(4 * k)
    y[:2] = math.sqrt(4 * k) * np.array([0.5, math.sqrt(3) / 2])
    q = isotropic_minima(full.isotropic(x, y, 1, 1))["positive_min_cost"]
    assert q == pytest.approx(2 * k)
    assert (k - 0.5) / q > (k - 1) / q
    assert k / q < 1


@pytest.mark.parametrize("k", [1, 2, 3])
def test_candidate_bounds_antiparallel_evidence(k):
    """With Qplus>0 fixed-pair costs scale as SD^-2 and the candidate cancels it."""
    f, a = np.zeros(k), np.zeros(k)
    f[0], a[0] = 2.0, -1.0
    lm = np.linspace(-14, 14, 1201)
    spans = []
    for sd in (0.1, 0.03, 0.01):
        law, eta = _candidate_law([full.isotropic(f, a, sd, sd)], 1.0, lm)
        spans.append(np.ptp(law - full._log_prior(lm, 1.0)))
        assert eta == pytest.approx(1.0)
    np.testing.assert_allclose(spans, np.full(3, spans[0]), rtol=2e-10, atol=1e-10)
    assert spans[0] > 0


def test_candidate_does_not_bound_a_zero_cost_endpoint_case():
    """A force-only pair has Qplus=0: no finite positive fit, but zero infimum.

    The candidate's bounded-tail argument requires Qplus>0. It must not be
    extended silently to every pair lacking a finite positive exact solution.
    """
    lm = np.linspace(-14, 14, 1201)
    spans = []
    for sd in (0.1, 0.01):
        reading = full.isotropic([1.0], [0.0], sd, sd)
        assert _weight(reading) == 1.0
        law, _ = _candidate_law([reading], 1.0, lm)
        spans.append(np.ptp(law - full._log_prior(lm, 1.0)))
    assert spans[1] / spans[0] == pytest.approx(100.0, rel=1e-10)


@pytest.mark.parametrize("k", [1, 2, 3])
def test_candidate_has_no_gate_and_leaves_acute_weights_alone(k):
    for magnitude in (0.0, 0.1, 1.0, 100.0):
        x, y = np.zeros(k), np.zeros(k)
        x[0], y[0] = 2 * magnitude, -magnitude
        assert 0 < _weight(full.isotropic(x, y, 1, 1)) <= 1
    x, y = np.zeros(k), np.zeros(k)
    x[0], y[0] = 4, 2
    if k > 1:
        y[1] = 3
    reading = full.isotropic(x, y, 1, 1)
    state = full.reading_state(reading, 1.0)
    assert _weight(reading) == pytest.approx(full.misfit_weight(state["Q"], k - 1))


@pytest.mark.parametrize("k", [1, 2, 3])
def test_candidate_commutes_with_swap_and_units(k):
    x, y = np.zeros(k), np.zeros(k)
    x[0], y[0] = 2, -1
    if k > 1:
        y[1] = 0.4
    original = full.isotropic(x, y, 0.7, 1.3)
    swapped = full.isotropic(y, x, 1.3, 0.7)
    units = full.isotropic(3 * x, 5 * y, 2.1, 6.5)
    assert _weight(swapped) == pytest.approx(_weight(original), rel=1e-13)
    assert _weight(units) == pytest.approx(_weight(original), rel=1e-13)
    lm = np.linspace(-8, 8, 401)
    law = _candidate_law([original], 2, lm)[0]
    swap_law = _candidate_law([swapped], 0.5, -lm)[0]
    unit_law = _candidate_law([units], 1.2, lm + math.log(3 / 5))[0]
    np.testing.assert_allclose(swap_law, law, atol=2e-11, rtol=1e-12)
    np.testing.assert_allclose(unit_law, law, atol=2e-11, rtol=1e-12)


@pytest.mark.parametrize("k", [1, 2, 3])
def test_exact_zeros_leave_candidate_law_and_eta_alone(k):
    f, a = np.zeros(k), np.zeros(k)
    f[0], a[0] = 3, -1
    readings = [full.isotropic(f, a, 0.4, 0.7), full.isotropic(-f, a, 0.4, 0.7)]
    zero = full.isotropic(np.zeros(k), np.zeros(k), 0.4, 0.7)
    assert _weight(zero) == 1.0
    lm = np.linspace(-8, 8, 401)
    base, eta = _candidate_law(readings, 2.0, lm)
    padded, padded_eta = _candidate_law(readings + [zero] * 4, 2.0, lm)
    assert padded_eta == pytest.approx(eta, abs=1e-12)
    np.testing.assert_allclose(padded, base, atol=1e-12)


def test_candidate_keeps_the_known_all_zero_cloud_point_counterexample():
    """g(0)=1 preserves the law, but cannot repair O2's cloud point on null input."""
    c = 10.0
    lm = np.linspace(-18, 18, 1801)
    zero = full.isotropic([0.0, 0.0], [0.0, 0.0], 1, 1)
    base, _ = _candidate_law([], c, lm)
    padded, _ = _candidate_law([zero], c, lm)
    np.testing.assert_array_equal(base, padded)
    _, probability = full._quantiles(lm, padded, [0.5])
    size = full._stacked_mean_length(np.exp(lm), [zero])
    point = np.sum(np.exp(lm) * size * probability) / np.sum(size * probability)
    assert point == pytest.approx(4.91319551, rel=3e-6)
    assert abs(math.log(point / c)) > 0.5
