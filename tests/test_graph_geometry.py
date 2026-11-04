"""T49: independent projector, Jacobian, curve-volume and restriction checks."""
import math

import numpy as np
import pytest
from scipy import integrate

from met.graph_geometry import (graph_chart_log_volume, graph_metric,
                                 graph_projector, graph_projector_derivative,
                                 scalar_chart_substitution_density, scalar_curve_volume)


def _covariance(dimension):
    base = np.eye(dimension)+np.arange(dimension**2).reshape(dimension, dimension)/(5*dimension**2)
    return base@base.T


@pytest.mark.parametrize("operator", [np.array([[1.7]]), np.array([[1.4, .2], [-.3, .8]]),
                                      np.array([[1., .3], [-.4, .8], [.2, -.1]])])
def test_graph_projector_is_the_declared_mean_subspace(operator):
    p, q = operator.shape
    covariance = _covariance(p+q)
    chol = np.linalg.cholesky(covariance)
    basis = np.linalg.solve(chol, np.vstack([operator, np.eye(q)]))
    projector = graph_projector(operator, covariance=covariance)
    np.testing.assert_allclose(projector@basis, basis, atol=2e-14)
    np.testing.assert_allclose(projector@projector, projector, atol=2e-14)
    assert np.trace(projector) == pytest.approx(q, abs=2e-14)


def test_full_covariance_projection_cost_matches_independent_physical_residual():
    operator = np.array([[1.8, .3], [-.2, 1.2]])
    covariance = _covariance(4)
    data = np.array([2., -.4, 1.1, .3])
    z = np.linalg.solve(np.linalg.cholesky(covariance), data)
    residual = z-graph_projector(operator, covariance=covariance)@z
    constraint = np.c_[np.eye(2), -operator]
    physical = constraint@data
    expected = physical@np.linalg.solve(constraint@covariance@constraint.T, physical)
    assert residual@residual == pytest.approx(expected, rel=2e-14)


@pytest.mark.parametrize("use_covariance", [False, True])
def test_projector_derivative_matches_five_point_numerical_geometry(use_covariance):
    operator = np.array([[1.3, .2], [-.4, .7]])
    direction = np.array([[.2, -.3], [.7, .1]])
    covariance = _covariance(4) if use_covariance else None
    step = 1e-4
    project = lambda t: graph_projector(operator+t*direction, covariance=covariance)
    numeric = (project(-2*step)-8*project(-step)+8*project(step)-project(2*step))/(12*step)
    actual = graph_projector_derivative(operator, direction, covariance=covariance)
    np.testing.assert_allclose(actual, numeric, atol=3e-12)
    assert graph_metric(operator, direction, covariance=covariance) == pytest.approx(np.sum(numeric**2)/2, abs=2e-12)


def test_identity_metric_matches_the_independent_graph_chart_formula():
    matrix = np.array([[1., .3], [-.4, .8], [.2, -.1]])
    first = np.array([[.1, .2], [-.1, .3], [.4, .2]])
    second = np.array([[.3, -.2], [.2, .1], [-.1, .4]])
    left = np.linalg.inv(np.eye(3)+matrix@matrix.T)
    right = np.linalg.inv(np.eye(2)+matrix.T@matrix)
    expected = np.trace(left@first@right@second.T)
    assert graph_metric(matrix, first, second) == pytest.approx(expected, abs=2e-14)


@pytest.mark.parametrize("matrix", [np.array([[1.2], [.3]]), np.array([[1.1, .2], [-.3, .8]]),
                                    np.array([[1., .3], [-.4, .8], [.2, -.1]])])
@pytest.mark.parametrize("use_covariance", [False, True])
def test_chart_volume_matches_independently_differentiated_projector_gram_determinant(matrix, use_covariance):
    covariance = _covariance(sum(matrix.shape)) if use_covariance else None
    step, derivatives = 1e-4, []
    for coordinate in range(matrix.size):
        basis = np.zeros_like(matrix)
        basis.flat[coordinate] = 1.
        project = lambda t: graph_projector(matrix+t*basis, covariance=covariance)
        derivatives.append((project(-2*step)-8*project(-step)+8*project(step)-project(2*step))/(12*step))
    gram = np.array([[np.sum(first*second)/2 for second in derivatives] for first in derivatives])
    sign, log_determinant = np.linalg.slogdet(gram)
    assert sign == 1
    assert graph_chart_log_volume(matrix, covariance=covariance) == pytest.approx(log_determinant/2, abs=2e-10)


@pytest.mark.parametrize("dimension", [1, 2, 3])
@pytest.mark.parametrize("mass", [.2, 1., 4.])
def test_scalar_restriction_has_intrinsic_curve_volume_independent_of_dimension_after_normalizing(dimension, mass):
    ratio = 2.3
    covariance = np.diag([ratio**2]*dimension+[1.]*dimension)
    metric = graph_metric(mass*np.eye(dimension), np.eye(dimension), covariance=covariance)
    speed = scalar_curve_volume(mass, dimension=dimension, noise_ratio=ratio)
    assert math.sqrt(metric) == pytest.approx(speed, rel=3e-14)
    normalized = speed/(math.sqrt(dimension)*math.pi/2)
    assert normalized == pytest.approx(2*ratio/(math.pi*(ratio**2+mass**2)))


def test_positive_curve_volume_is_finite_and_its_median_is_the_declared_noise_ratio():
    for dimension, ratio in ((1, .7), (3, 2.)):
        speed = lambda m: scalar_curve_volume(m, dimension=dimension, noise_ratio=ratio)
        total = integrate.quad(speed, 0, np.inf, epsabs=2e-12)[0]
        lower = integrate.quad(speed, 0, ratio, epsabs=2e-12)[0]
        assert total == pytest.approx(math.sqrt(dimension)*math.pi/2, abs=2e-12)
        assert lower/total == pytest.approx(.5, abs=2e-13)


@pytest.mark.parametrize("rho", [-.6, .6])
def test_correlated_noise_curve_volume_is_shifted_truncated_cauchy_with_same_median_ratio(rho):
    ratio = 2.
    covariance = np.array([[ratio**2, rho*ratio], [rho*ratio, 1.]])
    def speed(mass):
        return math.sqrt(graph_metric(np.array([[mass]]), np.ones((1, 1)), covariance=covariance))
    for mass in (.3, 1., 4.):
        explicit = ratio*math.sqrt(1-rho*rho)/(mass*mass-2*rho*ratio*mass+ratio*ratio)
        assert speed(mass) == pytest.approx(explicit, rel=3e-14)
        reciprocal = ratio*ratio/mass
        assert speed(mass) == pytest.approx(speed(reciprocal)*ratio*ratio/mass**2, rel=4e-14)
    total = integrate.quad(speed, 0, np.inf, epsabs=2e-11)[0]
    lower = integrate.quad(speed, 0, ratio, epsabs=2e-11)[0]
    assert total == pytest.approx(math.pi/2+math.asin(rho), abs=2e-11)
    assert lower/total == pytest.approx(.5, abs=2e-12)


def test_ambient_density_substitution_is_a_different_normalized_law_and_breaks_reciprocity():
    mass = .5
    for dimension in (1, 2):
        naive = lambda m: scalar_chart_substitution_density(m, dimension=dimension, noise_ratio=1.)
        assert integrate.quad(naive, 0, np.inf)[0] == pytest.approx(1., abs=2e-12)
        transformed = naive(1/mass)/mass**2
        if dimension == 1:
            assert naive(mass) == pytest.approx(transformed)
        else:
            assert abs(naive(mass)-transformed) > .3
            assert integrate.quad(naive, 0, 1.)[0] > .9
    # Evaluating the ambient matrix density on mI gives the exponent d^2.
    d = 2
    assert graph_chart_log_volume(mass*np.eye(d)) == pytest.approx(-d*d*math.log1p(mass*mass))


def test_two_transverse_tube_conventions_produce_distinct_scalar_limits():
    # Independently integrate ambient density over 3D traceless transverse
    # balls. Their finite-radius weights must approach different limits;
    # this does not just evaluate the proposed limiting density on the curve.
    radial, rw = np.polynomial.legendre.leggauss(6)
    cosine, cw = np.polynomial.legendre.leggauss(8)
    radial, rw = (radial+1)/2, rw/2
    directions = [np.diag([1., -1.])/math.sqrt(2),
                  np.array([[0., 1.], [0., 0.]]), np.array([[0., 0.], [1., 0.]])]

    def tube_weight(mass, epsilon, canonical_radius):
        radius = epsilon*(1+mass*mass) if canonical_radius else epsilon
        total = 0.
        for r, wr in zip(radial, rw):
            for c, wc in zip(cosine, cw):
                sine = math.sqrt(1-c*c)
                for azimuth in np.arange(16)*2*math.pi/16:
                    coordinates = [sine*math.cos(azimuth), sine*math.sin(azimuth), c]
                    offset = sum(x*basis for x, basis in zip(coordinates, directions))
                    matrix = mass*np.eye(2)+radius*r*offset
                    density = np.linalg.det(np.eye(2)+matrix.T@matrix)**-2
                    total += wr*wc*r*r*density*(2*math.pi/16)
        return radius**3*total

    scale = (1+1.5**2)/(1+.5**2)
    for canonical in (False, True):
        expected = scale if canonical else scale**4
        errors = []
        for epsilon in (.02, .01):
            actual = tube_weight(.5, epsilon, canonical)/tube_weight(1.5, epsilon, canonical)
            errors.append(abs(actual-expected))
        assert errors[1] < errors[0]
        assert errors[1]/expected < 1e-3


def test_covariance_scale_leaves_hypothesis_geometry_unchanged_but_does_not_select_fiber_volume():
    matrix = np.array([[1.4, .2], [-.3, .9]])
    covariance = _covariance(4)
    direction = np.eye(2)
    for scale in (.2, 7.):
        np.testing.assert_allclose(graph_projector(matrix, covariance=scale*covariance),
                                   graph_projector(matrix, covariance=covariance), atol=2e-14)
        assert graph_metric(matrix, direction, covariance=scale*covariance) == pytest.approx(graph_metric(matrix, direction, covariance=covariance))
        assert graph_chart_log_volume(matrix, covariance=scale*covariance) == pytest.approx(graph_chart_log_volume(matrix, covariance=covariance), abs=2e-14)


def test_automatic_product_geometry_prior_changes_on_zero_padding_with_new_instrument_geometry():
    # This hypothetical automatic prior rule uses declarations, not data.
    # Appending a zero reading adds no data factor but adds another g_i.
    g = lambda m, ratio: ratio**2/(m*m+ratio**2)**2
    original = lambda m: math.sqrt(g(m, 1.))
    appended = lambda m: math.sqrt(g(m, 1.)+g(m, 2.))
    old_z = integrate.quad(original, 0, np.inf)[0]
    new_z = integrate.quad(appended, 0, np.inf)[0]
    old_cdf = integrate.quad(original, 0, 1.)[0]/old_z
    new_cdf = integrate.quad(appended, 0, 1.)[0]/new_z
    assert old_cdf == pytest.approx(.5)
    assert abs(new_cdf-old_cdf) > .03


def test_geometry_input_validation():
    with pytest.raises(ValueError):
        graph_projector([1., 2.])
    with pytest.raises(ValueError):
        graph_projector(np.eye(2), covariance=np.diag([1., 1., 1., 0.]))
    with pytest.raises(ValueError):
        graph_metric(np.eye(2), np.ones((1, 2)))
    with pytest.raises(ValueError):
        scalar_curve_volume(1., dimension=0, noise_ratio=1.)
