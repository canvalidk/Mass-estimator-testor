"""B14: local modified-profile geometry and obstructions to a global mass law.

All candidate factors here are test-local mathematical constructions. The pole
example continues an unrestricted whole-line local formula into positive mass;
it does not call that continuation a regular likelihood or change production.
"""
import itertools
import math

import numpy as np
import pytest
from scipy import integrate, optimize


def _u(theta, dimension=1):
    return np.vstack([math.sin(theta) * np.eye(dimension),
                      math.cos(theta) * np.eye(dimension)])


def _line_mle(data):
    """One k=1 unrestricted whole-line fit, in its positive-amplitude chart."""
    if np.linalg.norm(data) == 0:
        raise ValueError("the direction is unidentified at measured zero")
    return math.atan2(data[0], data[1]), np.linalg.norm(data)


def _formal_factor(theta, data):
    """Formal continuation: Lp / |cos(theta-theta_hat)|, omitting data constants."""
    theta_hat, _ = _line_mle(data)
    along = _u(theta)[:, 0] @ data
    residual_energy = data @ data - along * along
    return math.exp(-residual_energy / 2) / abs(math.cos(theta - theta_hat))


def test_exact_saturated_chart_has_inverse_sample_space_jacobian():
    """Vary amplitude at fixed direction MLE; numerically invert the constrained fit."""
    theta_hat, amplitude_hat = 0.35, 1.7
    direction_hat = _u(theta_hat)[:, 0]
    for theta in [0.2, 0.7, 1.2]:
        direction = _u(theta)[:, 0]

        def constrained_fit(amplitude):
            return direction @ (amplitude * direction_hat)

        fitted = constrained_fit(amplitude_hat)
        h = 1e-5
        inverse = [optimize.brentq(lambda a: constrained_fit(a) - target, 0.1, 4.0)
                   for target in [fitted - h, fitted + h]]
        inverse_derivative = (inverse[1] - inverse[0]) / (2 * h)
        expected = 1 / math.cos(theta - theta_hat)
        assert inverse_derivative == pytest.approx(expected, rel=1e-9)
        if theta == 1.2:
            # Reversing the Jacobian would give cos, a different modification.
            assert abs(inverse_derivative - math.cos(theta - theta_hat)) > 0.5


def _gaussian_score_cross_covariance(theta, theta_hat, dimension):
    """Independent Gaussian quadrature of centered finite-difference nuisance scores."""
    nodes, weights = np.polynomial.hermite.hermgauss(3)
    indices = np.array(list(itertools.product(range(3), repeat=2 * dimension)))
    quadrature_weights = np.prod(weights[indices], axis=1) / math.pi**dimension
    noise = math.sqrt(2) * nodes[indices]
    basis, fitted_basis = _u(theta, dimension), _u(theta_hat, dimension)
    fitted_amplitude = 0.4 * np.arange(1, dimension + 1)
    mean = fitted_basis @ fitted_amplitude
    data = mean + noise
    constrained_amplitude = basis.T @ mean

    def scores(matrix, amplitude):
        out = []
        h = 1e-5
        for j in range(dimension):
            delta = np.zeros(dimension)
            delta[j] = h
            plus = -0.5 * np.sum((data - matrix @ (amplitude + delta))**2, axis=1)
            minus = -0.5 * np.sum((data - matrix @ (amplitude - delta))**2, axis=1)
            out.append((plus - minus) / (2 * h))
        score = np.array(out).T
        return score - quadrature_weights @ score

    first = scores(basis, constrained_amplitude)
    second = scores(fitted_basis, fitted_amplitude)
    return first.T @ (quadrature_weights[:, None] * second)


@pytest.mark.parametrize("dimension", [1, 3])
def test_expected_nuisance_score_covariance_has_cosine_block(dimension):
    theta, theta_hat = 0.27, 0.81
    covariance = _gaussian_score_cross_covariance(theta, theta_hat, dimension)
    np.testing.assert_allclose(covariance, math.cos(theta - theta_hat) * np.eye(dimension),
                               rtol=1e-9, atol=1e-9)
    assert abs(np.linalg.det(covariance)) == pytest.approx(
        abs(math.cos(theta - theta_hat))**dimension, rel=1e-9)


def test_formal_positive_mass_continuation_has_nonintegrable_pole():
    """Valid data have a negative unrestricted fit and a positive-mass continuation pole.

    The positive-mass constrained optimum is a boundary case, so the global
    continuation is not claimed to satisfy regular BN construction assumptions.
    """
    data = np.array([-1.0, 2.0])
    theta_hat, _ = _line_mle(data)
    pole = math.atan(2.0)
    assert theta_hat == pytest.approx(-math.atan(0.5))
    assert pole - theta_hat == pytest.approx(math.pi / 2)
    # K at the pole is |data|^2=5, and |cos| vanishes linearly.
    for h in [1e-3, 1e-4, 1e-5]:
        for sign in [-1, 1]:
            assert h * _formal_factor(pole + sign * h, data) == pytest.approx(
                math.exp(-2.5), rel=3e-6)

    # The normalized half-Cauchy prior is uniform in positive angle.
    density = lambda theta: (2 / math.pi) * _formal_factor(theta, data)
    cutoffs = [1e-2, 1e-3, 1e-4]
    integrals = [integrate.quad(density, 0, pole - epsilon, epsabs=1e-10)[0] +
                 integrate.quad(density, pole + epsilon, math.pi / 2, epsabs=1e-10)[0]
                 for epsilon in cutoffs]
    assert np.all(np.diff(integrals) > 0)
    # Each decade adds a nonzero constant, rather than converging to a normalizer.
    expected_growth = (4 / math.pi) * math.exp(-2.5) * math.log(10)
    assert integrals[-1] - integrals[-2] == pytest.approx(expected_growth, rel=1e-5)


def test_measured_zero_has_no_direction_independent_formal_extension():
    with pytest.raises(ValueError, match="unidentified"):
        _formal_factor(0.6, np.zeros(2))
    candidates = [0.3, 0.8]
    limiting_ratios = []
    for direction in [0.2, 0.9]:
        tiny_data = 1e-8 * _u(direction)[:, 0]
        ratio = _formal_factor(candidates[0], tiny_data) / _formal_factor(candidates[1], tiny_data)
        limit = abs(math.cos(candidates[1] - direction) / math.cos(candidates[0] - direction))
        assert ratio == pytest.approx(limit, rel=1e-13)
        assert abs(ratio - 1) > 0.1
        limiting_ratios.append(ratio)
    assert abs(limiting_ratios[0] - limiting_ratios[1]) > 0.1


def test_physical_nuisance_reparameterization_changes_only_data_constant():
    """j_alpha,alpha and the sample Jacobian cancel candidate-dependent scaling."""
    fitted_mass, fitted_alpha = 0.5, 0.8
    fitted_candidate = np.array([fitted_mass, 1.0])
    fitted_gram = fitted_candidate @ fitted_candidate
    theta_hat = math.atan(fitted_mass)
    ratios = []
    for mass in [0.2, 0.7, 2.0, 8.0]:
        candidate = np.array([mass, 1.0])
        gram = candidate @ candidate

        def constrained_alpha(alpha):
            return candidate @ (alpha * fitted_candidate) / gram

        h = 1e-5
        forward = (constrained_alpha(fitted_alpha + h) -
                   constrained_alpha(fitted_alpha - h)) / (2 * h)
        physical_adjustment = 1 / (abs(forward) * math.sqrt(gram))
        arc_adjustment = 1 / abs(math.cos(math.atan(mass) - theta_hat))
        ratios.append(physical_adjustment / arc_adjustment)
    np.testing.assert_allclose(ratios, 1 / math.sqrt(fitted_gram), rtol=1e-9)


def _formal_physical_mass_density(force, acceleration, sf, sa, mass):
    centre = sf / sa
    data = np.array([force / sf, acceleration / sa])
    theta = math.atan(mass / centre)
    prior = 2 * centre / (math.pi * (centre**2 + mass**2))
    return prior * _formal_factor(theta, data)


def test_formal_factor_respects_units_and_swap_without_becoming_proper():
    force, acceleration, sf, sa = 1.7, 0.8, 2.0, 0.7
    force_units, acceleration_units = 3.0, 0.2
    scale = force_units / acceleration_units
    for mass in [0.3, 0.8, 2.0, 8.0]:
        original = _formal_physical_mass_density(force, acceleration, sf, sa, mass)
        scaled = _formal_physical_mass_density(force_units * force, acceleration_units * acceleration,
                                               force_units * sf, acceleration_units * sa, scale * mass)
        assert scaled == pytest.approx(original / scale, rel=1e-13)
        swapped = _formal_physical_mass_density(acceleration, force, sa, sf, 1 / mass)
        assert swapped / mass**2 == pytest.approx(original, rel=1e-13)


def test_private_profile_corrections_are_not_the_common_mass_correction():
    """Expected cross-score determinants use the joint unrestricted mass estimate."""
    observations = [np.array([1.0, 1.0]), np.array([0.0, 2.0])]
    scatter = sum(np.outer(z, z) for z in observations)
    _, eigenvectors = np.linalg.eigh(scatter)
    fitted_direction = eigenvectors[:, -1]
    if fitted_direction[1] < 0:
        fitted_direction *= -1
    shared_hat = math.atan2(*fitted_direction)
    assert shared_hat == pytest.approx(0.5 * math.atan(0.5))
    private_hats = [_line_mle(z)[0] for z in observations]
    correction_ratios = []
    for theta in [0.2, 0.4, 0.9]:
        direction = _u(theta)[:, 0]
        shared_cross_determinant = (direction @ fitted_direction)**2
        private_cross_determinant = math.prod(direction @ _u(hat)[:, 0] for hat in private_hats)
        shared_correction = 1 / abs(shared_cross_determinant)
        private_product = 1 / abs(private_cross_determinant)
        correction_ratios.append(shared_correction / private_product)
    # The difference varies with candidate mass, so a data-only factor cannot remove it.
    assert np.ptp(correction_ratios) > 0.1
