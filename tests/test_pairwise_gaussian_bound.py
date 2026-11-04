"""T68: independent Gaussian cubature of pairwise profile improvements.

Every graph projector is constructed from the physical mean relation and
declared covariance here. No pairwise-anchor runtime helper is imported.
Polynomial moments use exact-degree Gauss--Hermite rules. Positive-part
and Laplace cubatures are numerical checks with explicit refinement;
they are not proofs of a uniform bound or limiting interval coverage.
"""
from functools import lru_cache
import math

import numpy as np
import pytest


IDENTITY = np.eye(2)
CORRELATED = np.array([[1., .9], [.9, 1.]])
DIAGONAL = np.diag([.1, 2.])

# One or two physical acquisitions, each with its own free fixed signal.
# Radii are signed whitened signal lengths, used only to specify worlds.
CASES = [
    (.6, 3., (0.,), (IDENTITY,)),
    (1.2, .3, (.25,), (IDENTITY,)),
    (1.2, 4/15, (4.,), (CORRELATED,)),
    (.5, 5., (15.,), (DIAGONAL,)),
    (1.2, 1.2, (10.,), (CORRELATED,)),
    (.8, 2., (0., 0.), (IDENTITY, CORRELATED)),
    (.6, 4., (0., 3.), (IDENTITY, DIAGONAL)),
    (1.2, 4/15, (2., -4.), (CORRELATED, DIAGONAL)),
    (2.3, .7, (.1, .5),
     (np.array([[1., -.4], [-.4, 1.]]), np.diag([4., .25]))),
]


def _physical_experiment(true_mass, candidate_mass, radii, covariances):
    dimension = 2*len(radii)
    difference = np.zeros((dimension, dimension))
    true_projector = np.zeros_like(difference)
    mean = np.zeros(dimension)
    for index, (radius, covariance) in enumerate(zip(radii, covariances)):
        whitening = np.linalg.solve(np.linalg.cholesky(covariance), IDENTITY)
        true_column = whitening@np.array([true_mass, 1.])
        candidate_column = whitening@np.array([candidate_mass, 1.])
        p0 = np.outer(true_column, true_column)/(true_column@true_column)
        pt = np.outer(candidate_column, candidate_column)/(candidate_column@candidate_column)
        block = slice(2*index, 2*index+2)
        true_projector[block, block] = p0
        difference[block, block] = pt-p0
        mean[block] = radius*true_column/np.linalg.norm(true_column)
    return difference, mean, true_projector


@lru_cache(maxsize=None)
def _normal_rule(order, dimension):
    nodes, weights = np.polynomial.hermite.hermgauss(order)
    grid = np.indices((order,)*dimension).reshape(dimension, -1).T
    points = math.sqrt(2)*nodes[grid]
    probabilities = np.prod(weights[grid]/math.sqrt(math.pi), axis=1)
    return points, probabilities


def _fields(difference, mean, order):
    nodes, weights = _normal_rule(order, len(mean))
    observations = nodes+mean
    projected = observations@difference.T
    improvement = np.einsum("ni,ni->n", observations, projected)/2
    energy = np.einsum("ni,ni->n", projected, projected)
    return improvement, energy, weights


def _moments(difference, mean):
    square = difference@difference
    a = float(mean@square@mean)
    trace = float(np.trace(square))
    return a, trace, a+trace


@pytest.mark.parametrize("case", CASES)
def test_gaussian_profile_difference_has_physical_mean_variance_and_energy(case):
    difference, mean, true_projector = _physical_experiment(*case)
    np.testing.assert_allclose(true_projector@mean, mean, rtol=3e-14, atol=3e-14)
    assert abs(np.trace(difference)) < 2e-14
    assert np.linalg.norm(difference, 2) <= 1+3e-14
    a, trace, b = _moments(difference, mean)
    # The negative mean identity needs the true mean in the true graph;
    # it would be false for an arbitrary mean and symmetric difference.
    assert mean@difference@mean == pytest.approx(-a, rel=5e-13, abs=2e-13)
    improvement, energy, weights = _fields(difference, mean, 3)
    expected = -a/2
    assert weights@improvement == pytest.approx(expected, rel=5e-13, abs=2e-13)
    assert weights@((improvement-expected)**2) == pytest.approx(
        a+trace/2, rel=8e-13, abs=2e-13)
    assert weights@energy == pytest.approx(b, rel=5e-13, abs=2e-13)


@pytest.mark.parametrize("case", CASES)
def test_centered_fourth_moment_matches_exact_gaussian_formula_and_bound(case):
    difference, mean, _ = _physical_experiment(*case)
    a, trace, b = _moments(difference, mean)
    fourth_operator = np.linalg.matrix_power(difference, 4)
    exact_fourth = (3*(a+trace/2)**2+12*mean@fourth_operator@mean
                    +3*np.trace(fourth_operator))
    # Degree eight in Gaussian coordinates: five Hermite nodes per
    # coordinate integrate this polynomial exactly before roundoff.
    improvement, _, weights = _fields(difference, mean, 5)
    observed = float(weights@((improvement+a/2)**4))
    assert observed == pytest.approx(exact_fourth, rel=2e-12, abs=2e-12)
    assert exact_fourth <= 3*b*b+12*b+2e-10
    assert exact_fourth <= 24*(b*b+b)+2e-10


def test_fourth_moment_bound_is_centered_not_a_raw_profile_moment():
    difference, mean, _ = _physical_experiment(1.2, 4/15, (100.,), (CORRELATED,))
    _, _, b = _moments(difference, mean)
    improvement, _, weights = _fields(difference, mean, 5)
    assert weights@(improvement**4) > 1000*24*(b*b+b)


def _exact_log_laplace(difference, mean, parameter):
    square = difference@difference
    denominator = np.eye(len(mean))+2*parameter*square
    sign, log_determinant = np.linalg.slogdet(denominator)
    assert sign == 1
    return (-log_determinant/2
            -parameter*mean@square@np.linalg.solve(denominator, mean))


@pytest.mark.parametrize("case", CASES)
def test_exact_laplace_transform_has_the_uniform_exponential_bound(case):
    difference, mean, _ = _physical_experiment(*case)
    _, _, b = _moments(difference, mean)
    for parameter in (1e-4, .05, .25, 1., 10., 1000.):
        logarithm = _exact_log_laplace(difference, mean, parameter)
        bound = -parameter*b/(1+2*parameter)
        assert logarithm <= bound+3e-11
    # Exponential Markov at s=1/4: P(V <= b/2) <= exp(-b/24).
    tilted_log_bound = b/8-b/6
    assert tilted_log_bound == pytest.approx(-b/24, rel=2e-15, abs=1e-14)


@pytest.mark.parametrize("case", [CASES[i] for i in (0, 1, 2, 5, 6, 7, 8)])
def test_laplace_identity_agrees_with_refined_physical_gaussian_cubature(case):
    difference, mean, _ = _physical_experiment(*case)
    for parameter in (.05, .25):
        target = math.exp(_exact_log_laplace(difference, mean, parameter))
        estimates = []
        for order in (15, 23):
            _, energy, weights = _fields(difference, mean, order)
            estimates.append(float(weights@np.exp(-parameter*energy)))
        assert estimates[1] == pytest.approx(target, rel=2e-9, abs=2e-12)
        assert estimates[0] == pytest.approx(estimates[1], rel=2e-6, abs=2e-9)


@pytest.mark.parametrize("case", CASES)
def test_positive_self_normalized_improvement_checks_uniform_moment_constants(case):
    difference, mean, _ = _physical_experiment(*case)
    estimates = []
    for order in ((43, 73) if len(mean) == 2 else (15, 23)):
        improvement, energy, weights = _fields(difference, mean, order)
        positive = np.maximum(improvement, 0)/np.sqrt(1+energy)
        first, second = float(weights@positive), float(weights@(positive**2))
        assert first <= math.sqrt(18)
        assert second <= 100
        estimates.append((first, second))
    # The positive part is nonsmooth. This explicit numerical refinement
    # is descriptive checking, not polynomial exactness or a proof over
    # all instruments and signals. Constants above have ample margin.
    np.testing.assert_allclose(estimates[0], estimates[1], rtol=.03, atol=.004)
