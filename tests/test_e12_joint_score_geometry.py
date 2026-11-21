"""Finite checks of the shared-accelerometer score tensor and integrability."""
import numpy as np
import pytest
from numpy.polynomial.hermite import hermgauss
from scipy.integrate import quad


def graph(masses, whitening):
    basis = whitening[:, :2]
    raw = whitening @ np.r_[masses, 1.]
    size = np.linalg.norm(raw)
    u = raw/size
    du = (basis-np.outer(u, u @ basis))/size
    return u, du


def normal_quadrature(dim, order=4):
    nodes, weights = hermgauss(order)
    axes = np.meshgrid(*([np.sqrt(2)*nodes]*dim), indexing="ij")
    node_weights = np.meshgrid(*([weights/np.sqrt(np.pi)]*dim), indexing="ij")
    return np.array(axes).reshape(dim, -1).T, np.prod(node_weights, axis=0).ravel()


def test_exact_score_covariance_expected_hessian_and_tensor_rotation():
    m = np.array([.8, 1.5])
    whitening = np.array([[1., .2, -.1], [.1, .7, .3], [0., -.1, 1.3]])
    u, du = graph(m, whitening)
    metric = du.T @ du
    a = -.87
    exact = a*(whitening @ np.r_[m, 1.])
    r = u @ exact
    noise, weights = normal_quadrature(3)
    ys = exact + noise
    scores = ((ys @ u)[:, None])*(ys @ du)
    np.testing.assert_allclose(weights @ scores, 0, atol=1e-14)
    np.testing.assert_allclose(scores.T @ (weights[:, None]*scores), (1+r**2)*metric)
    score_step = 1e-5
    hessian = np.empty((2, 2))
    def mean_score(candidate):
        unit, deriv = graph(candidate, whitening)
        return weights @ ((ys @ unit)[:, None]*(ys @ deriv))
    for j in range(2):
        half_step = np.eye(2)[j]*score_step
        hessian[:, j] = (mean_score(m+half_step)-mean_score(m-half_step))/(2*score_step)
    np.testing.assert_allclose(-hessian, r**2*metric, atol=1e-9)
    rotation = np.array([[0., -1., 0.], [0., 0., 1.], [1., 0., 0.]])
    unit_rotated, derivative_rotated = graph(m, rotation @ whitening)
    np.testing.assert_allclose(derivative_rotated.T @ derivative_rotated, metric)
    np.testing.assert_allclose((ys @ rotation.T @ unit_rotated)[:, None] *
                               (ys @ rotation.T @ derivative_rotated), scores)


def test_correlated_physical_sensors_require_distinct_marginal_discounts():
    m = np.ones(2)
    subtract_acceleration = np.array([[1., 0., -1.], [0., 1., -1.], [0., 0., 1.]])
    metrics, rank_weights = [], []
    for sd, acceleration in [([.5, 1., 1.], 1.), ([1., .5, 1.], 2.)]:
        whitening = np.diag(1/np.array(sd)) @ subtract_acceleration
        inv = np.linalg.inv(whitening)
        physical_covariance = inv @ inv.T
        assert physical_covariance[0, 2] != 0.
        u, du = graph(m, whitening)
        rank_weights.append((acceleration*(u @ whitening @ np.r_[m, 1.]))**2)
        metrics.append(du.T @ du)
    np.testing.assert_allclose(metrics, [np.diag([4., 1.]), np.diag([1., 4.])])
    a = sum(w*metric for w, metric in zip(rank_weights, metrics))
    b = a + np.sum(metrics, axis=0)
    np.testing.assert_allclose(a, np.diag([8., 17.]))
    np.testing.assert_allclose(b, np.diag([13., 22.]))
    sandwich = np.linalg.inv(a) @ b @ np.linalg.inv(a)
    eta = np.diag(np.linalg.inv(a))/np.diag(sandwich)
    np.testing.assert_allclose(eta, [8/13, 17/22])
    assert np.sqrt(eta[1]/eta[0]) == pytest.approx(1.12057, abs=.001)


def test_conditional_two_mass_score_has_chi_square_first_two_moments():
    # Fixed observed longitudinal components, including a negative push and a zero.
    m = [.6, 1.3]
    whitening = [np.eye(3), np.diag([.6, 1.1, 1.5]), np.diag([1., 1.4, 2.])]
    longitudinal = [1.32, -.48, 0.]
    noise, weights = normal_quadrature(4)
    scores = np.zeros((len(weights), 2))
    conditional_covariance = np.zeros((2, 2))
    for i, (measurement, observed) in enumerate(zip(whitening, longitudinal)):
        u, deriv = graph(m, measurement)
        conditional_covariance += observed**2 * (deriv.T @ deriv)
        if observed:
            _, _, tangent = np.linalg.svd(u[None, :], full_matrices=True)
            normal = noise[:, 2*i:2*i+2] @ tangent[1:]
            scores += observed * (normal @ deriv)
    pivot = np.einsum("ij,ij->i", scores, np.linalg.solve(conditional_covariance, scores.T).T)
    np.testing.assert_allclose(scores.T @ (weights[:, None]*scores), conditional_covariance)
    assert weights @ pivot == pytest.approx(2.)
    assert weights @ (pivot-2.)**2 == pytest.approx(4.)


def cost_curve(mass, observations):
    u, du = graph(mass, np.eye(3))
    along = observations @ u
    costs = np.sum((observations-along[:, None]*u)**2, axis=1)
    gradients = -2*along[:, None]*(observations @ du)
    return costs, gradients


def integrate_rectangle(center, size, observations, row_weights):
    vertices = center + size*np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
    circulation = 0.
    for start, end in zip(vertices, np.roll(vertices, -1, axis=0)):
        step = end-start
        def rate(t):
            costs, gradients = cost_curve(start+t*step, observations)
            return row_weights(costs) @ (gradients @ step)
        circulation += quad(rate, 0., 1., epsabs=1e-13)[0]
    return circulation


def test_nonintegrable_symmetric_discount_has_circulation_on_actual_mass_rectangle():
    center = np.array([1., 1.])
    u, _ = graph(center, np.eye(3))
    e = np.array([1., -1., 0.])/np.sqrt(2)
    f = np.cross(u, e)
    observations = np.array([u+e, u+np.sqrt(2)*f])
    k, grads = cost_curve(center, observations)
    np.testing.assert_allclose(k, [1., 2.])
    expected_curl = (k[0]-k[1])*np.linalg.det(grads)/(1+np.prod(k))**2
    assert expected_curl < -.1
    def invalid(k):
        return np.ones(2)/(1+np.prod(k))
    for size, tolerance in [(.01, .04), (.001, .004)]:
        result = integrate_rectangle(center, size, observations, invalid)
        assert result/size**2 == pytest.approx(expected_curl, rel=tolerance)
        for valid in [lambda z: 1/(1+z), lambda z: np.ones(2)/(1+sum(z))]:
            assert abs(integrate_rectangle(center, size, observations, valid)) < 1e-12


def test_bounded_conservative_weights_obey_parallel_line_bound():
    k = np.array([.13, 12., .6])
    l = np.array([.7, 1., 3.])
    for step in [.05, 10., 1000.]:
        weight = lambda s: sum(1/(1+k+s)-1/(1+l+s))
        integral = quad(weight, 0., step)[0]
        potential = lambda x: np.log1p(x).sum()
        assert integral == pytest.approx(potential(k+step)-potential(k)
                                         -potential(l+step)+potential(l))
        assert abs(integral) <= np.abs(k-l).sum()
