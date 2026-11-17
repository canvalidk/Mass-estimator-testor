"""B27 finite Gaussian identities. No sampled coverage or rate simulations."""

from itertools import product

import numpy as np
import pytest
from numpy.polynomial.hermite import hermgauss


def gaussian_nodes(dimension):
    """Integrate degree-four polynomials against independent unit normals."""
    nodes, weights = hermgauss(3)
    indices = np.array(list(product(range(len(nodes)), repeat=dimension)))
    return np.sqrt(2) * nodes[indices], np.prod(weights[indices], axis=1) / np.pi**(dimension / 2)


def pair_direction(mass, ratio):
    return np.array([mass, ratio]) / np.hypot(mass, ratio)


def projector(basis):
    frame, _ = np.linalg.qr(basis, mode="reduced")
    return frame @ frame.T


def test_vector_pooling_and_common_gaussian_scatter_moments():
    rng = np.random.default_rng(271)
    groups = [rng.normal(size=(size, 2)) for size in (1, 4, 2)]
    flat = np.concatenate(groups)
    np.testing.assert_allclose(sum(group.T @ group for group in groups), flat.T @ flat)

    direction = pair_direction(1.7, .8)
    nodes, weights = gaussian_nodes(2)
    for radius in (0., -.5, 2.):
        gaussian_pairs = radius * direction + nodes
        scatter = gaussian_pairs[:, :, None] * gaussian_pairs[:, None, :]
        mean = np.eye(2) + radius**2 * np.outer(direction, direction)
        np.testing.assert_allclose(np.einsum("i,ijk->jk", weights, scatter), mean, atol=1e-14)
        squared_error = np.sum((scatter - mean)**2, axis=(1, 2))
        assert weights @ squared_error == pytest.approx(6 * (1 + radius**2))


@pytest.mark.parametrize("ratio", [.01, 1., 7., 200.])
def test_unequal_ratio_chord_and_log_mass_score(ratio):
    m, truth, radius = 3., .5, 1.2
    u0 = pair_direction(truth, ratio)
    u = pair_direction(m, ratio)
    perpendicular = np.array([u[1], -u[0]])
    speed = truth * ratio / (ratio**2 + truth**2)
    chord = abs(perpendicular @ u0)
    assert chord == pytest.approx(abs(m-truth)*ratio / (np.hypot(m, ratio)*np.hypot(truth, ratio)))
    reference_speed = ratio / (1 + ratio**2)
    assert .5 * reference_speed <= speed <= 2 * reference_speed

    # At truth, the profile quadratic score is the product of independent
    # parallel N(radius, 1) and perpendicular N(0, 1) coordinates.
    nodes, weights = gaussian_nodes(2)
    score = -2 * speed * (radius + nodes[:, 0]) * nodes[:, 1]
    assert weights @ score == pytest.approx(0., abs=1e-14)
    assert weights @ score**2 == pytest.approx(4 * speed**2 * (1 + radius**2))


def test_unbounded_ratios_separate_unweighted_and_weighted_benchmarks():
    indices = np.arange(1., 1001.)
    ratios = indices**2
    reference_speed_squared = (ratios / (1 + ratios**2))**2
    # L is bounded by sum i^-4, which is below the integral bound 4/3.
    assert np.sum(reference_speed_squared) < 4/3
    # With r_i^2=1, total whitened energy grows like D. Pairwise KL of the
    # fixed-signal Gaussian experiments is bounded by half the sum i^-4.
    u1 = np.array([np.ones_like(ratios), ratios]) / np.hypot(1., ratios)
    u2 = np.array([2 * np.ones_like(ratios), ratios]) / np.hypot(2., ratios)
    assert np.sum((u1-u2)**2) / 2 < 2/3
    # With r_i^2=i^3, weighted signal energy grows like the harmonic sum.
    weighted_signal = reference_speed_squared * indices**3
    np.testing.assert_array_less(1 / (5 * indices), weighted_signal)
    np.testing.assert_array_less(weighted_signal, 1 / indices)


def test_proper_subspace_compression_and_general_profile_moments():
    # B selects a shared one-dimensional physical subspace from two channels.
    basis = np.eye(4)[:, [0, 2]]
    outside = np.eye(4)[:, [1, 3]]
    factor = np.array([[2., 0., 0., 0.], [.3, 1.1, 0., 0.],
                       [-.5, .4, 1.3, 0.], [.6, -.2, .2, .8]])
    covariance = factor @ factor.T
    precision = np.linalg.inv(covariance)
    compressed = np.linalg.inv(basis.T @ precision @ basis)
    cross = basis.T @ covariance @ outside
    conditional_map = basis.T - cross @ np.linalg.solve(outside.T @ covariance @ outside, outside.T)
    np.testing.assert_allclose(compressed @ basis.T @ precision, conditional_map, atol=1e-14)
    np.testing.assert_allclose(conditional_map @ covariance @ conditional_map.T, compressed, atol=1e-14)
    np.testing.assert_allclose(conditional_map @ covariance @ outside, 0., atol=1e-14)

    truth, candidate, acceleration = 1.3, .7, .9
    whitened_basis = np.linalg.solve(factor, basis)
    true_graph = whitened_basis @ np.array([[truth], [1.]])
    candidate_graph = whitened_basis @ np.array([[candidate], [1.]])
    mean = acceleration * true_graph[:, 0]
    contrast = projector(true_graph) - projector(candidate_graph)
    row = np.array([1., -candidate])
    gap = ((candidate-truth) * acceleration)**2 / (row @ compressed @ row)
    assert mean @ contrast @ mean == pytest.approx(gap)
    nodes, weights = gaussian_nodes(4)
    readings = mean + nodes
    values = np.einsum("ni,ij,nj->n", readings, contrast, readings)
    assert weights @ values == pytest.approx(gap)
    assert weights @ (values-gap)**2 == pytest.approx(4*gap + 2*np.sum(contrast**2))


def test_correlated_histories_equal_information_and_energy_but_bounded_kl():
    masses = (1., 2.)
    energy, information = np.zeros(2), np.zeros(2)
    kl = 0.
    added_information = 0.
    for pair in range(1, 13):
        epsilon = 2.**(-pair)
        for label, mass in enumerate(masses):
            whitener = np.array([[1., -mass], [0., mass*epsilon]])
            graph = whitener @ np.array([[mass], [1.]])
            # Each truth puts its only push in the matching reading of the pair.
            means = np.zeros((2, 2))
            means[label] = graph[:, 0] / mass
            tangent = whitener @ np.array([1., 0.])
            log_mass_score = (np.eye(2) - projector(graph)) @ tangent
            information[label] += log_mass_score @ log_mass_score
            energy += np.sum(means**2, axis=1)
            difference = means[0] - means[1]
            kl += (difference @ difference) / 2

        # Two uncorrelated readings add unit energy in each world. Reciprocal
        # ratios swap the angular speeds at the two masses, preserving equality
        # of the exact information oracles while their directions become close.
        speeds = []
        for ratio in (1/epsilon, 2*epsilon):
            directions = np.array([pair_direction(mass, ratio) for mass in masses])
            energy += np.sum(directions**2, axis=1)
            speed = np.array([mass*ratio / (mass**2 + ratio**2) for mass in masses])
            speeds.append(speed)
            information += speed**2
            difference = directions[0] - directions[1]
            assert np.linalg.norm(difference) <= epsilon
            kl += (difference @ difference) / 2
        np.testing.assert_allclose(speeds[0], speeds[1][::-1])
        added_information += sum(speed[0]**2 for speed in speeds)
        np.testing.assert_allclose(information, pair + added_information)
        np.testing.assert_allclose(energy, 2*pair + (1.-4.**(-pair)) / 3)
        assert kl <= 2*(1.-4.**(-pair)) / 3
        assert kl < 2/3
