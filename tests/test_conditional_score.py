"""Independent geometry and distribution checks for T45's research tests."""
import math

import numpy as np
import pytest
from scipy import integrate, optimize, stats

from met.conditional_score import gaussian_score_accepts, gaussian_score_test
from met.full import Reading


def _scalar_readings(data):
    return [Reading([row[0]], [row[1]], np.eye(2)) for row in data]


def _example():
    matrix = np.array([[1.4, .2, -.1, .3], [.1, 1.1, .4, 0.],
                       [.3, -.2, 1.7, .2], [0., .1, -.1, 1.2]])
    return Reading([2., -.7], [.5, 1.1], matrix@matrix.T,
                   np.array([[.6], [.8]]))


@pytest.mark.parametrize("mass", [.13, 1., 4.7])
def test_full_geometry_score_matches_independent_profile_finite_difference(mass):
    reading = _example()
    data = np.r_[reading.force, reading.acceleration]
    inverse = np.linalg.inv(reading.cov)
    angle = math.atan(mass)

    def half_profile(phi):
        basis = np.vstack([math.sin(phi)*reading.subspace,
                           math.cos(phi)*reading.subspace])
        b = basis.T @ inverse @ data
        return .5*float(b @ np.linalg.solve(basis.T@inverse@basis, b))

    step = 1e-4
    derivative = (half_profile(angle-2*step)-8*half_profile(angle-step)
                  +8*half_profile(angle+step)-half_profile(angle+2*step))/(12*step)
    actual = gaussian_score_test([reading], mass)
    assert actual["score"] == pytest.approx(derivative, abs=1e-11)


@pytest.mark.parametrize("parallel", [[.3, 1.4], [-3., .2], [2., -1.3]])
@pytest.mark.parametrize("probability", [.5, .95])
def test_exact_conditional_coverage_by_independent_residual_quadrature(parallel, probability):
    # Fix the realized parallel observations, rather than declaring a law
    # for signal amplitudes. Integrate independent residual Gaussians.
    parallel = np.asarray(parallel)
    mass, angle = 2., math.atan(2.)
    u = np.array([math.sin(angle), math.cos(angle)])
    v = np.array([math.cos(angle), -math.sin(angle)])
    z = stats.norm.ppf((1+probability)/2)

    def score(w0, w1):
        observations = parallel[:, None]*u + np.array([w0, w1])[:, None]*v
        return gaussian_score_test(_scalar_readings(observations), mass)["normal_score"]

    def conditional(w0):
        bracket = (abs(parallel[0]*w0)+z*np.linalg.norm(parallel)+1)/abs(parallel[1])
        roots = sorted(optimize.brentq(lambda w1: score(w0, w1)-bound, -bracket, bracket)
                       for bound in (-z, z))
        return stats.norm.pdf(w0)*(stats.norm.cdf(roots[1])-stats.norm.cdf(roots[0]))

    value, error = integrate.quad(conditional, -12., 12., epsabs=1e-11)
    assert value == pytest.approx(probability, abs=2e-11)
    assert error < 1e-9


def test_full_covariance_conditional_cf_from_gaussian_residual_integration():
    # Independently construct whitened U,V using an SVD, not the helper's
    # QR. The remaining three Gaussian coordinates are integrated by GH.
    base, mass = _example(), 1.7
    chol = np.linalg.cholesky(base.cov)
    n = base.subspace
    basis = np.linalg.solve(chol, np.vstack([mass*n, n]))
    rotation, _, _ = np.linalg.svd(basis, full_matrices=True)
    observed_u = .83*rotation[:, 0]
    nodes, weights = np.polynomial.hermite.hermgauss(12)
    characteristic = 0j
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            for k in range(len(nodes)):
                residual = math.sqrt(2)*rotation[:, 1:]@nodes[[i, j, k]]
                data = chol@(observed_u+residual)
                reading = Reading(data[:2], data[2:], base.cov, n)
                normal = gaussian_score_test([reading], mass)["normal_score"]
                characteristic += weights[i]*weights[j]*weights[k]*np.exp(.7j*normal)/math.pi**1.5
    assert characteristic.real == pytest.approx(math.exp(-.7**2/2), abs=1e-12)
    assert characteristic.imag == pytest.approx(0., abs=1e-13)


def test_full_type_zero_padding_frame_basis_units_and_swap_invariance():
    reading, mass = _example(), 2.3
    base = gaussian_score_test([reading], mass)
    angle = .7
    rotation = np.array([[math.cos(angle), -math.sin(angle)],
                         [math.sin(angle), math.cos(angle)]])
    transform = np.zeros((4, 4))
    transform[:2, :2] = 7*rotation
    transform[2:, 2:] = .2*rotation
    data = transform@np.r_[reading.force, reading.acceleration]
    changed = Reading(data[:2], data[2:], transform@reading.cov@transform.T,
                      -rotation@reading.subspace)
    assert gaussian_score_test([changed], mass*7/.2)["statistic"] == pytest.approx(base["statistic"], rel=2e-13)
    swapped = Reading(reading.acceleration, reading.force,
                      reading.cov[np.ix_([2, 3, 0, 1], [2, 3, 0, 1])], reading.subspace)
    assert gaussian_score_test([swapped], 1/mass)["statistic"] == pytest.approx(base["statistic"], rel=2e-13)
    reversed_sign = Reading(-reading.force, -reading.acceleration, reading.cov, reading.subspace)
    assert gaussian_score_test([reversed_sign], mass)["statistic"] == pytest.approx(base["statistic"])
    zero = Reading([0.], [0.], np.array([[2., .3], [.3, 1.]]))
    assert gaussian_score_test([reading, zero], mass) == base


def test_multiple_subspace_dimensions_use_full_conditional_score_variance():
    data = np.array([[2., -.3], [.2, 1.1]])
    vector = Reading(data[0], data[1], np.eye(4))
    components = _scalar_readings(data.T)
    for mass in (.3, 2., 8.):
        full = gaussian_score_test([vector], mass)
        split = gaussian_score_test(components, mass)
        assert full["statistic"] == pytest.approx(split["statistic"], abs=2e-13)
        assert full["variance"] == pytest.approx(split["variance"], abs=2e-13)


def test_isotropic_score_curve_and_maximum_against_scatter_eigenvalues():
    data = np.array([[4., .2], [-.3, 1.4], [1., -.8]])
    scatter = data.T@data
    eigenvalues = np.linalg.eigvalsh(scatter)
    results = []
    for theta in np.linspace(-math.pi/2+.001, math.pi/2-.001, 401):
        # Full projective formula checked independently. The public helper
        # accepts positive masses, so verify it on that part of the curve.
        u = np.array([math.sin(theta), math.cos(theta)])
        v = np.array([math.cos(theta), -math.sin(theta)])
        along, score = u@scatter@u, u@scatter@v
        value = score**2/along
        results.append(value)
        if theta > 0:
            actual = gaussian_score_test(_scalar_readings(data), math.tan(theta))
            assert actual["statistic"] == pytest.approx(value, abs=4e-14)
    maximum = (math.sqrt(eigenvalues[1])-math.sqrt(eigenvalues[0]))**2
    assert max(results) <= maximum*(1+1e-13)
    assert max(results) == pytest.approx(maximum, rel=2e-4)


def test_confidence_set_can_be_disconnected_and_include_bad_stationary_line():
    # Principal line at -0.1 radians, minor line at pi/2-0.1. The positive
    # domain contains an accepted lower tail AND the wrong minor-line branch.
    phi = -.1
    u = np.array([math.sin(phi), math.cos(phi)])
    v = np.array([math.cos(phi), -math.sin(phi)])
    data = np.array([10*u, v])
    angles = np.array([.01, .6, math.pi/2-.1, math.pi/2-.01])
    accepted = gaussian_score_accepts(_scalar_readings(data), np.tan(angles), probability=.95)
    assert accepted.tolist() == [True, False, True, False]
    assert gaussian_score_test(_scalar_readings(data), math.tan(angles[2]))["statistic"] < 1e-25


def test_full_domain_acceptance_is_not_a_proper_law_central_interval():
    readings = _scalar_readings(np.eye(2))
    assert gaussian_score_accepts(readings, [1e-8, 1., 1e8], probability=.1).all()
    assert gaussian_score_accepts(readings, 2., probability=.95)
    # The half-Cauchy starting law has finite central quantiles, while the
    # exact score set is the whole positive domain for this scatter.
    low, high = stats.halfcauchy.ppf([.025, .975])
    assert 0 < low < high < math.inf
    assert gaussian_score_accepts(readings, high*100, probability=.95)


def test_empty_and_all_zero_series_have_an_explicit_degenerate_convention():
    for readings in ([], _scalar_readings(np.zeros((2, 2)))):
        value = gaussian_score_test(readings, 2.)
        assert value["degenerate"] and value["statistic"] == 0 and value["p_value"] == 1


@pytest.mark.parametrize("scale", [1., 3., 1e4])
def test_nonzero_observations_with_zero_candidate_projection_follow_declared_convention(scale):
    reading = Reading([scale], [-2*scale], np.eye(2))
    actual = gaussian_score_test([reading], 2.)
    assert actual["degenerate"] and actual["score"] == actual["variance"] == 0.
    assert actual["statistic"] == 0. and actual["p_value"] == 1.
    assert gaussian_score_accepts([reading], 2., probability=.95)


def test_optional_zero_guard_does_not_require_a_representable_precision_vector():
    # Cholesky whitening is finite even though Sigma^-1 y overflows. The
    # optional guard must leave the usable QR statistic available.
    reading = Reading([1e10], [1.], np.diag([1e-310, 1.]))
    actual = gaussian_score_test([reading], 1e100)
    assert math.isfinite(actual["statistic"]) and not actual["degenerate"]
    assert actual["statistic"] == pytest.approx(1.)


def test_extreme_statistic_saturates_to_rejection_without_a_power_overflow():
    acceleration = np.nextafter(np.nextafter(-2e160, 0.), 0.)
    reading = Reading([1e160], [acceleration], np.eye(2))
    actual = gaussian_score_test([reading], 2.)
    assert not actual["degenerate"] and actual["statistic"] == math.inf
    assert actual["p_value"] == 0.


@pytest.mark.parametrize("mass", [0., -1., math.inf, math.nan])
def test_invalid_candidate_mass_is_rejected(mass):
    with pytest.raises(ValueError):
        gaussian_score_test([_example()], mass)


def test_membership_validation_and_explicit_probability():
    with pytest.raises(ValueError):
        gaussian_score_accepts([], [0.], probability=.95)
    with pytest.raises(ValueError):
        gaussian_score_accepts([], [1.], probability=1.)
    with pytest.raises(ValueError):
        gaussian_score_test([object()], 1.)
