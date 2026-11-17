"""E2: finite checks for score combination and its null-block limit.

These deterministic identities and fixtures check ingredients of the argument;
they do not simulate coverage or prove the asymptotic quantifiers.
"""

import math

import numpy as np
import pytest
from scipy import integrate, optimize, stats

from met.conditional_score import gaussian_score_accepts, gaussian_score_test
from met.full import Reading


def _scalar_readings(rows):
    return [Reading([row[0]], [row[1]], np.eye(2)) for row in rows]


def _full_blocks():
    factor = np.array([[1.4, 0., 0., 0.], [.2, 1.1, 0., 0.],
                       [-.3, .1, 1.3, 0.], [.4, -.2, .3, .9]])
    restricted = Reading([2., -.7], [.5, 1.1], factor @ factor.T,
                         np.array([[.6], [.8]]))
    unrestricted = Reading([.4, 1.7], [-.6, .9], np.diag([.8, 2., 1.3, .7]))
    correlated = Reading([3.], [.2], np.array([[2., -.4], [-.4, 1.]]))
    return [[restricted], [unrestricted, correlated]]


@pytest.mark.parametrize("mass", [.13, 1., 4.7])
def test_retained_score_and_variance_reproduce_full_type_pooling(mass):
    left, right = _full_blocks()
    a = gaussian_score_test(left, mass)
    b = gaussian_score_test(right, mass)
    pooled = gaussian_score_test(left + right, mass)
    signed = a["score"] + b["score"]
    variance = a["variance"] + b["variance"]
    assert pooled["score"] == pytest.approx(signed, abs=1e-13)
    assert pooled["variance"] == pytest.approx(variance, abs=1e-13)
    assert pooled["statistic"] == pytest.approx(signed**2 / variance, abs=1e-13)
    # Conditional on each block's own projections, this is weighted normal
    # combination. The weights must be retained with the signed pivots.
    combined = (math.sqrt(a["variance"]) * a["normal_score"]
                + math.sqrt(b["variance"]) * b["normal_score"]) / math.sqrt(variance)
    assert pooled["normal_score"] == pytest.approx(combined, abs=1e-13)


@pytest.mark.parametrize("mass", [.1, 1., 8.])
def test_score_state_combination_ignores_zero_readings_and_order(mass):
    left, right = _full_blocks()
    zero = Reading([0.], [0.], np.array([[2., .3], [.3, 1.]]))
    readings = left + right
    original = gaussian_score_test(readings, mass)
    assert gaussian_score_test([zero] + readings[::-1] + [zero], mass) == original


def test_identical_entire_pvalue_reports_do_not_determine_the_combined_set():
    # Both positive scatters t I have zero score at every mass, strictly
    # positive variance, and p-value one. This uses no zero-denominator rule.
    blocks = [_scalar_readings(math.sqrt(t) * np.eye(2)) for t in (1., 3.)]
    masses = np.geomspace(.01, 100., 21)
    for block in blocks:
        for mass in masses:
            report = gaussian_score_test(block, mass)
            assert report["statistic"] < 1e-27
            assert report["p_value"] == pytest.approx(1., abs=2e-14)
            assert report["variance"] > 0
        assert gaussian_score_accepts(block, masses, probability=.95).all()

    extra = _scalar_readings([[3., -1.]])
    # At m=1 the extra block has S=4 and C=2. The first two reports
    # contribute S=0 but C=1 and C=3 respectively.
    first = gaussian_score_test(blocks[0] + extra, 1.)
    second = gaussian_score_test(blocks[1] + extra, 1.)
    assert first["statistic"] == pytest.approx(16/3)
    assert second["statistic"] == pytest.approx(16/5)
    assert not gaussian_score_accepts(blocks[0] + extra, 1., probability=.95)
    assert gaussian_score_accepts(blocks[1] + extra, 1., probability=.95)


def test_null_scatter_harmonics_have_the_stated_covariance_by_quadrature():
    nodes, weights = np.polynomial.hermite.hermgauss(3)
    xx, yy = np.meshgrid(math.sqrt(2)*nodes, math.sqrt(2)*nodes, indexing="ij")
    ww = (np.outer(weights, weights)/math.pi).ravel()
    x, y = xx.ravel(), yy.ravel()
    # Individual summands are not Gaussian or independent. The CLT uses
    # these means/covariances, giving independent Gaussian limits.
    terms = np.column_stack(((x*x-y*y)/2, x*y, (x*x+y*y-2)/2))
    np.testing.assert_allclose(ww @ terms, 0., atol=1e-14)
    np.testing.assert_allclose(terms.T @ (ww[:, None]*terms), np.eye(3), atol=1e-14)
    for theta, phi in ((.13, .8), (.4, 1.2), (.7, .7)):
        direct = (np.cos(theta)*x-np.sin(theta)*y)*(np.sin(theta)*x+np.cos(theta)*y)
        harmonic = terms[:, 0]*math.sin(2*theta)+terms[:, 1]*math.cos(2*theta)
        other = terms[:, 0]*math.sin(2*phi)+terms[:, 1]*math.cos(2*phi)
        np.testing.assert_allclose(direct, harmonic, atol=1e-14)
        assert ww @ (direct*other) == pytest.approx(math.cos(2*(theta-phi)), abs=1e-14)


def test_deterministic_scatter_sequence_loses_a_fixed_external_block():
    # A finite ingredient check: realize N I + sqrt(N) A by two rows.
    # This is not a claim that this deterministic path is a null sample.
    perturbation = np.array([[2., 1.2], [1.2, -.4]])
    external_rows = np.array([[6., 3.]])
    external = _scalar_readings(external_rows)
    external_gram = external_rows.T @ external_rows
    angles = np.array([.1, .3, .7, 1.1, 1.4])
    limit = 1.2*np.sin(2*angles)+1.2*np.cos(2*angles)
    errors, differences = [], []
    for n in (100., 10000., 1000000.):
        gram = n*np.eye(2)+math.sqrt(n)*perturbation
        rows = np.linalg.cholesky(gram).T
        np.testing.assert_allclose(rows.T @ rows, gram, rtol=1e-14, atol=1e-12)
        null_block = _scalar_readings(rows)
        pooled_values, null_values = [], []
        for theta in angles:
            u = np.array([math.sin(theta), math.cos(theta)])
            v = np.array([math.cos(theta), -math.sin(theta)])
            leading = v @ perturbation @ u
            expected = ((leading + (v @ external_gram @ u)/math.sqrt(n))
                        / math.sqrt(1+(u @ perturbation @ u)/math.sqrt(n)
                                    +(u @ external_gram @ u)/n))
            pooled = gaussian_score_test(external + null_block, math.tan(theta))["normal_score"]
            pure = gaussian_score_test(null_block, math.tan(theta))["normal_score"]
            assert pooled == pytest.approx(expected, abs=5e-10)
            pooled_values.append(pooled)
            null_values.append(pure)
        errors.append(np.max(np.abs(np.asarray(pooled_values)-limit)))
        differences.append(np.max(np.abs(np.asarray(pooled_values)-null_values)))
    assert errors[2] < errors[1] < errors[0]
    assert differences[2] < differences[1] < differences[0]
    assert errors[2] < .03
    assert differences[2] < differences[0]/50
    baseline = [gaussian_score_test(external, math.tan(t))["normal_score"] for t in angles]
    assert np.max(np.abs(np.asarray(pooled_values)-baseline)) > 1


@pytest.mark.parametrize("z1,z2", [(1., 0.), (0., 1.), (2., -.4), (-3., .2)])
def test_limit_curve_supremum_is_the_two_harmonic_radius(z1, z2):
    # The closed arc has the same supremum as the open positive-mass arc.
    maximum_angle = (math.atan2(z1, z2) % math.pi)/2
    value = z1*math.sin(2*maximum_angle)+z2*math.cos(2*maximum_angle)
    assert value*value == pytest.approx(z1*z1+z2*z2, abs=1e-14)
    grid = np.linspace(0., math.pi/2, 1001)
    sampled = z1*np.sin(2*grid)+z2*np.cos(2*grid)
    assert np.max(sampled*sampled) <= z1*z1+z2*z2+1e-13


@pytest.mark.parametrize("probability", [.5, .95, .99])
def test_limit_whole_domain_probability_by_independent_normal_integral(probability):
    q = stats.chi2.ppf(probability, 1)
    radius = math.sqrt(q)

    def integrand(angle):
        x = radius*math.sin(angle)
        height = radius*math.cos(angle)
        return stats.norm.pdf(x)*(stats.norm.cdf(height)-stats.norm.cdf(-height))*height

    measured, error = integrate.quad(integrand, -math.pi/2, math.pi/2,
                                     epsabs=1e-12, epsrel=1e-12)
    whole = 1-math.exp(-q/2)
    assert measured == pytest.approx(whole, abs=1e-12)
    assert error < 1e-10
    assert whole < probability
    if probability == .95:
        assert whole == pytest.approx(.8534999355139157, abs=1e-13)


def _own_projection_weighted_score(readings, mass, *, damped):
    # On this scalar, unit-noise slice with angle derivative, each C_i is
    # precisely ||U_i||^2. On the full type a coordinate-dependent C_i must
    # not be silently substituted for that invariant projected norm.
    reports = [gaussian_score_test([reading], mass) for reading in readings]
    weights = [1/(1+report["variance"]) if damped else 1. for report in reports]
    score = math.fsum(weight*report["score"] for weight, report in zip(weights, reports))
    variance = math.fsum(weight*weight*report["variance"] for weight, report in zip(weights, reports))
    return score, variance, score/math.sqrt(variance) if variance else 0.


def test_reading_own_weights_have_the_exact_conditional_variance_by_quadrature():
    parallel, mass = np.array([.5, 2.]), 1.7
    angle = math.atan(mass)
    u = np.array([math.sin(angle), math.cos(angle)])
    v = np.array([math.cos(angle), -math.sin(angle)])
    nodes, weights = np.polynomial.hermite.hermgauss(3)
    for damped in (False, True):
        h = 1/(1+parallel**2) if damped else np.ones_like(parallel)
        expected_variance = np.sum(h*h*parallel**2)
        mean, second_moment = 0., 0.
        for i, first in enumerate(math.sqrt(2)*nodes):
            for j, second in enumerate(math.sqrt(2)*nodes):
                data = parallel[:, None]*u + np.array([first, second])[:, None]*v
                score, variance, _ = _own_projection_weighted_score(_scalar_readings(data), mass, damped=damped)
                probability_weight = weights[i]*weights[j]/math.pi
                assert variance == pytest.approx(expected_variance, abs=1e-13)
                mean += probability_weight*score
                second_moment += probability_weight*score*score
        assert mean == pytest.approx(0., abs=1e-13)
        assert second_moment == pytest.approx(expected_variance, abs=1e-13)


def test_reading_own_weights_keep_conditional_coverage_by_independent_integration():
    parallel, mass, probability = np.array([.5, 2.]), 1.7, .95
    angle = math.atan(mass)
    u = np.array([math.sin(angle), math.cos(angle)])
    v = np.array([math.cos(angle), -math.sin(angle)])
    threshold = stats.norm.ppf((1+probability)/2)
    for damped in (False, True):
        def statistic(first, second):
            data = parallel[:, None]*u + np.array([first, second])[:, None]*v
            return _own_projection_weighted_score(_scalar_readings(data), mass, damped=damped)[2]

        def conditional(first):
            roots = [optimize.brentq(lambda second: statistic(first, second)-bound, -100., 100.)
                     for bound in (-threshold, threshold)]
            return stats.norm.pdf(first)*(stats.norm.cdf(roots[1])-stats.norm.cdf(roots[0]))

        measured, error = integrate.quad(conditional, -10., 10., epsabs=1e-11, epsrel=1e-11)
        assert measured == pytest.approx(probability, abs=2e-11)
        assert error < 1e-9


def test_distinct_own_projection_weights_give_different_reachable_verdicts():
    # At theta=pi/4 these rows have parallel coordinates (3,1) and
    # residual coordinates (1,-4). Both exact conditional procedures exist,
    # but their realized 95% verdicts differ.
    data = np.array([[4., 2.], [-3., 5.]])/math.sqrt(2)
    readings = _scalar_readings(data)
    flat = _own_projection_weighted_score(readings, 1., damped=False)
    damped = _own_projection_weighted_score(readings, 1., damped=True)
    assert flat[0] == pytest.approx(-1.)
    assert flat[1] == pytest.approx(10.)
    assert flat[2]**2 == pytest.approx(.1)
    assert damped[0] == pytest.approx(-1.7)
    assert damped[1] == pytest.approx(.34)
    assert damped[2]**2 == pytest.approx(8.5)
    critical = stats.chi2.ppf(.95, 1)
    assert flat[2]**2 < critical < damped[2]**2
    zero = _scalar_readings([[0., 0.]])
    for weighted in (False, True):
        original = _own_projection_weighted_score(readings, 1., damped=weighted)
        reordered = _own_projection_weighted_score(zero + readings[::-1], 1., damped=weighted)
        assert reordered == original
