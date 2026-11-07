"""Independent profile, geometry and conditional distribution checks for T50."""
import math

import numpy as np
import pytest
from scipy import integrate, optimize, stats

from met.conditional_score import gaussian_score_test
from met.full import Reading
from met.operator_score import (gaussian_operator_profile,
                                gaussian_operator_score_test,
                                operator_tangent_basis)


def _example():
    factor = np.array([[1.4, .2, -.1, .3], [.1, 1.1, .4, 0.],
                       [.3, -.2, 1.7, .2], [0., .1, -.1, 1.2]])
    return [Reading([2., -.7], [.5, 1.1], factor@factor.T),
            Reading([-.3, 1.7], [1.2, -.2], np.diag([.8, 1.3, 1.5, .7]))]


def _test(readings, matrix, tangents=None):
    if tangents is None:
        tangents = operator_tangent_basis(len(matrix), symmetric=False)
    return gaussian_operator_score_test(readings, matrix, tangents=tangents,
                                       rank_rtol=1e-11, rank_atol=0.)


def _profile_normal_equations(readings, matrix):
    # Independent raw-coordinate precision calculation, without whitened QR.
    basis = np.vstack([matrix, np.eye(len(matrix))])
    result = 0.
    for reading in readings:
        data = np.r_[reading.force, reading.acceleration]
        precision = np.linalg.inv(reading.cov)
        amplitude = np.linalg.solve(basis.T@precision@basis, basis.T@precision@data)
        residual = data-basis@amplitude
        result -= .5*residual@precision@residual
    return float(result)


def test_full_covariance_profile_and_vector_score_against_independent_differences():
    readings = _example()
    matrix = np.array([[1.7, -.2], [.4, .9]])
    tangents = operator_tangent_basis(2, symmetric=False)
    actual = _test(readings, matrix, tangents)
    assert gaussian_operator_profile(readings, matrix) == pytest.approx(
        _profile_normal_equations(readings, matrix), abs=2e-14)
    step = 1e-4
    finite_difference = []
    for direction in tangents:
        values = [_profile_normal_equations(readings, matrix+t*step*direction)
                  for t in (-2, -1, 1, 2)]
        finite_difference.append((values[0]-8*values[1]+8*values[2]-values[3])/(12*step))
    assert actual["score"] == pytest.approx(finite_difference, abs=1e-11)

    # Differentiate the projector itself, independently of QR tangents.
    expected = np.zeros((4, 4))
    for reading in readings:
        chol = np.linalg.cholesky(reading.cov)
        y = np.linalg.solve(chol, np.r_[reading.force, reading.acceleration])

        def projector(candidate):
            basis = np.linalg.solve(chol, np.vstack([candidate, np.eye(2)]))
            return basis@np.linalg.solve(basis.T@basis, basis.T)

        p = projector(matrix)
        columns = []
        for direction in tangents:
            derivative = (projector(matrix+step*direction)-projector(matrix-step*direction))/(2*step)
            columns.append((np.eye(4)-p)@derivative@p@y)
        g = np.column_stack(columns)
        expected += g.T@g
    assert actual["covariance"] == pytest.approx(expected, rel=2e-8, abs=2e-9)


@pytest.mark.parametrize("count,rank", [(1, 2), (2, 4), (3, 4)])
def test_conditional_gaussian_moments_by_independent_hermite_quadrature(count, rank):
    matrix = np.array([[1.3, -.2], [.4, .8]])
    declarations = (_example()+[Reading([.8, -1.1], [.4, .6], np.diag([1.3, .9, .7, 1.1]))])[:count]
    geometries = []
    for i, reading in enumerate(declarations):
        chol = np.linalg.cholesky(reading.cov)
        basis = np.linalg.solve(chol, np.vstack([matrix, np.eye(2)]))
        rotation, _, _ = np.linalg.svd(basis, full_matrices=True)
        u = rotation[:, :2]@np.array([[.7, -1.2], [1.1, .4], [-.3, .9]][i])
        geometries.append((chol, rotation[:, 2:], u, reading.cov))
    nodes, weights = np.polynomial.hermite.hermgauss(3)
    means = np.zeros(4)
    second = np.zeros((4, 4))
    first_statistic, second_statistic = 0., 0.
    covariance = None
    for indices in np.ndindex(*(3,)*(2*count)):
        w = math.prod(weights[list(indices)])/math.pi**count
        readings = []
        for i, (chol, perpendicular, u, cov) in enumerate(geometries):
            data = chol@(u+math.sqrt(2)*perpendicular@nodes[list(indices[2*i:2*i+2])])
            readings.append(Reading(data[:2], data[2:], cov))
        result = _test(readings, matrix)
        assert result["rank"] == rank
        means += w*result["score"]
        second += w*np.outer(result["score"], result["score"])
        first_statistic += w*result["statistic"]
        second_statistic += w*result["statistic"]**2
        if covariance is None:
            covariance = result["covariance"]
        else:
            assert result["covariance"] == pytest.approx(covariance, abs=2e-14)
    assert means == pytest.approx(np.zeros(4), abs=2e-15)
    assert second == pytest.approx(covariance, abs=2e-14)
    assert first_statistic == pytest.approx(rank, abs=2e-13)
    assert second_statistic == pytest.approx(rank*(rank+2), abs=2e-12)


def test_rank_two_conditional_acceptance_probability_by_gaussian_quadrature():
    matrix = np.array([[1.3, -.2], [.4, .8]])
    reading = _example()[0]
    chol = np.linalg.cholesky(reading.cov)
    basis = np.linalg.solve(chol, np.vstack([matrix, np.eye(2)]))
    rotation, _, _ = np.linalg.svd(basis, full_matrices=True)
    u = rotation[:, :2]@np.array([.7, -1.2])
    threshold = stats.chi2.ppf(.9, 2)

    def statistic(first, second):
        data = chol@(u+rotation[:, 2:]@np.array([first, second]))
        return _test([Reading(data[:2], data[2:], reading.cov)], matrix)["statistic"]

    def conditional(phi):
        # This change of integration variable removes the square-root
        # boundary singularity without prescribing the score's CDF.
        first = math.sqrt(threshold)*math.sin(phi)
        positive = optimize.brentq(lambda second: statistic(first, second)-threshold,
                                   0., 2*math.sqrt(threshold), xtol=1e-13)
        negative = optimize.brentq(lambda second: statistic(first, second)-threshold,
                                   -2*math.sqrt(threshold), 0., xtol=1e-13)
        return (math.sqrt(threshold)*math.cos(phi)*stats.norm.pdf(first)
                *(stats.norm.cdf(positive)-stats.norm.cdf(negative)))

    value, error = integrate.quad(conditional, -math.pi/2, math.pi/2,
                                  epsabs=1e-10, epsrel=1e-10)
    assert value == pytest.approx(.9, abs=2e-10)
    assert error < 1e-8


def test_units_frames_swap_and_invertible_tangent_reparameterization():
    matrix = np.array([[1.7, -.2], [.4, .9]])
    basis = operator_tangent_basis(2, symmetric=False)
    readings = _example()
    base = _test(readings, matrix, basis)
    # Independent invertible force/acceleration coordinate maps encompass
    # ordinary unit changes and simultaneous orthogonal frame changes.
    force_map = np.array([[3., .4], [-.2, 2.]])
    acceleration_map = np.array([[.3, -.1], [.05, .5]])
    transform = np.zeros((4, 4))
    transform[:2, :2], transform[2:, 2:] = force_map, acceleration_map
    inverse_acceleration = np.linalg.inv(acceleration_map)
    changed_matrix = force_map@matrix@inverse_acceleration
    changed_tangents = np.array([force_map@e@inverse_acceleration for e in basis])
    changed = []
    swapped = []
    permutation = [2, 3, 0, 1]
    for reading in readings:
        data = transform@np.r_[reading.force, reading.acceleration]
        changed.append(Reading(data[:2], data[2:], transform@reading.cov@transform.T))
        swapped.append(Reading(reading.acceleration, reading.force,
                               reading.cov[np.ix_(permutation, permutation)]))
    assert _test(changed, changed_matrix, changed_tangents)["statistic"] == pytest.approx(base["statistic"], abs=2e-13)
    inverse = np.linalg.inv(matrix)
    swap_tangents = np.array([-inverse@e@inverse for e in basis])
    assert _test(swapped, inverse, swap_tangents)["statistic"] == pytest.approx(base["statistic"], abs=2e-13)
    jacobian = np.array([[1., .2, 0., .1], [0., 1.2, -.2, 0.],
                         [.1, 0., .8, .3], [0., .1, 0., 1.1]])
    transformed_basis = np.einsum("ba,bij->aij", jacobian, basis)
    reparameterized = _test(readings, matrix, transformed_basis)
    assert reparameterized["score"] == pytest.approx(jacobian.T@base["score"], abs=2e-14)
    assert reparameterized["covariance"] == pytest.approx(jacobian.T@base["covariance"]@jacobian, abs=2e-14)
    assert reparameterized["statistic"] == pytest.approx(base["statistic"], abs=2e-13)


@pytest.mark.parametrize("mass", [.2, 2., 7.])
def test_scalar_matrix_restriction_is_t45_and_not_full_matrix_test(mass):
    readings = _example()
    scalar = gaussian_score_test(readings, mass)
    restricted = _test(readings, mass*np.eye(2), np.eye(2)[None])
    assert restricted["rank"] == 1
    assert restricted["statistic"] == pytest.approx(scalar["statistic"], abs=5e-14)
    assert restricted["score"][0]*(1+mass**2) == pytest.approx(scalar["score"], abs=2e-14)
    assert restricted["covariance"][0, 0]*(1+mass**2)**2 == pytest.approx(scalar["variance"], abs=2e-13)
    full = _test(readings, mass*np.eye(2))
    assert full["rank"] == 4
    assert abs(full["statistic"]-restricted["statistic"]) > .1


def test_isotropic_profile_is_bingham_tls_and_scalar_trace_scatter():
    data = np.array([[3., .2, 1., -.4], [.7, -2., .1, 1.2],
                     [1.1, .3, -.5, .7], [-.2, .6, 1., .4]])
    readings = [Reading(y[:2], y[2:], np.eye(4)) for y in data]
    scatter = data.T@data
    matrix = np.array([[1.1, -.2], [.3, .7]])
    graph = np.vstack([matrix, np.eye(2)])
    projector = graph@np.linalg.solve(graph.T@graph, graph.T)
    expected = -.5*np.trace(scatter)+.5*np.trace(scatter@projector)
    assert gaussian_operator_profile(readings, matrix) == pytest.approx(expected, abs=2e-14)
    eigenvalues, eigenvectors = np.linalg.eigh(scatter)
    plane = eigenvectors[:, -2:]
    tls = np.linalg.solve(plane[2:].T, plane[:2].T).T
    assert gaussian_operator_profile(readings, tls) == pytest.approx(-.5*sum(eigenvalues[:2]), abs=2e-14)
    for mass in (.2, 1., 3.):
        reduced = np.array([[np.trace(scatter[:2, :2]), np.trace(scatter[:2, 2:])],
                            [np.trace(scatter[2:, :2]), np.trace(scatter[2:, 2:])]])
        u = np.array([mass, 1.])/math.hypot(mass, 1.)
        scalar_expected = -.5*np.trace(scatter)+.5*u@reduced@u
        assert gaussian_operator_profile(readings, mass*np.eye(2)) == pytest.approx(scalar_expected, abs=2e-14)


def test_zero_padding_rank_zero_and_recognized_nonzero_zero_projection():
    matrix = np.array([[1.7, -.2], [.4, .9]])
    readings = _example()
    base = _test(readings, matrix)
    padded = _test(readings+[Reading([0., 0.], [0., 0.], np.eye(4))], matrix)
    assert padded["rank"] == base["rank"]
    assert padded["score"] == pytest.approx(base["score"], abs=0.)
    assert padded["covariance"] == pytest.approx(base["covariance"], abs=0.)
    assert padded["statistic"] == pytest.approx(base["statistic"], abs=2e-14)
    for data in ([], [Reading([0., 0.], [0., 0.], np.eye(4))],
                 [Reading([1., 2.], [-2., -4.], np.eye(4))]):
        result = _test(data, 2*np.eye(2))
        assert result["degenerate"] and result["rank"] == 0
        assert result["statistic"] == 0. and result["p_value"] == 1.


def test_symmetric_basis_and_explicit_numerical_rank_tolerance():
    symmetric = operator_tangent_basis(2, symmetric=True)
    assert symmetric.shape == (3, 2, 2)
    assert np.einsum("aij,bij->ab", symmetric, symmetric) == pytest.approx(np.eye(3), abs=1e-15)
    assert _test(_example(), np.array([[1., .2], [.2, .8]]), symmetric)["rank"] == 3
    basis = operator_tangent_basis(2, symmetric=False)
    basis[-1] *= 1e-14
    result = _test(_example(), np.eye(2), basis)
    assert result["rank"] == 3
    assert result["rank_cutoff"] > 0.


@pytest.mark.parametrize("count,full_rank,symmetric_rank", [(1, 3, 3), (2, 6, 5), (3, 9, 6)])
def test_three_dimensional_rank_counts_match_the_unexamined_operator_directions(count, full_rank, symmetric_rank):
    matrix = np.diag([1., 1.5, 2.])
    latent = np.array([[1., .2, -.3], [-.1, .7, .4], [.2, -.4, 1.1]])
    readings = [Reading(matrix@vector, vector, np.eye(6)) for vector in latent[:count]]
    assert _test(readings, matrix)["rank"] == full_rank
    symmetric = operator_tangent_basis(3, symmetric=True)
    assert _test(readings, matrix, symmetric)["rank"] == symmetric_rank


def test_validation_rejects_restricted_subspaces_and_invalid_declarations():
    restricted = Reading([1., 2.], [.3, .4], np.eye(4), np.array([[1.], [0.]]))
    with pytest.raises(ValueError, match="full subspace"):
        _test([restricted], np.eye(2))
    with pytest.raises(ValueError):
        _test(_example(), np.ones((2, 3)))
    with pytest.raises(ValueError):
        _test(_example(), np.eye(2), np.zeros((0, 2, 2)))
    with pytest.raises(ValueError):
        gaussian_operator_score_test([], np.eye(2), tangents=np.eye(2)[None],
                                     rank_rtol=-1., rank_atol=0.)
    with pytest.raises(ValueError):
        operator_tangent_basis(0, symmetric=False)


def test_svd_statistic_overflow_obeys_the_explicit_range_rejection():
    reading = Reading([1e200], [1e-200], np.eye(2))
    with pytest.raises(ValueError, match="statistic exceeds the floating-point range"):
        _test([reading], np.zeros((1, 1)), np.ones((1, 1, 1)))
