"""B17's stopping counterexample: a declared metric does not fix scale weights.

These are checks of T48's explicitly conditional integration, not another
estimator. The brief's proposed forcing premise fails; its consequence and
forward-model studies are therefore not pursued.
"""
import math

import numpy as np
import pytest
from scipy import integrate, special, stats


def _density(cost, k, nu, b):
    """T48 section 2, including its absolute (not mass-law) constant."""
    alpha = nu / 2
    exponent = alpha + k / 2 - b
    return math.exp(b * math.log(alpha) + special.gammaln(exponent)
                    - special.gammaln(alpha) - k / 2 * math.log(math.pi * nu)
                    - exponent * math.log1p(cost / nu))


def _radial_integral(density, k):
    area = 2 * math.pi ** (k / 2) / special.gamma(k / 2)
    value, error = integrate.quad(
        lambda radius: area * radius ** (k - 1) * density(radius * radius),
        0, np.inf, epsabs=2e-10, epsrel=2e-10, limit=300)
    assert error < 2e-8
    return value


@pytest.mark.parametrize("b", [-1.0, -0.5, 0.0, 0.5, 1.0])
@pytest.mark.parametrize("k", [1, 3])
def test_residual_total_is_noise_size_moment_not_only_b_zero(k, b):
    """Independent residual-space quadrature catches the brief's false claim."""
    nu = 10.0
    actual = _radial_integral(lambda cost: _density(cost, k, nu, b), k)
    moment = math.exp(b * math.log(5) + special.gammaln(5 - b) - special.gammaln(5))
    assert actual == pytest.approx(moment, rel=2e-9)
    if b in (-1.0, 0.0):
        assert actual == pytest.approx(1.0, abs=2e-9)
    else:
        assert abs(actual - 1) > 0.01


@pytest.mark.parametrize("nu", [0.5, 10.0])
@pytest.mark.parametrize("k", [1, 3])
def test_b_minus_one_uses_same_explicit_gaussian_and_ig_but_another_normalized_density(k, nu):
    residual = np.arange(1, k + 1, dtype=float) / (k + 1)
    cost = float(residual @ residual)

    def integrand(v):
        log_gaussian = -k / 2 * math.log(2 * math.pi * v) - cost / (2 * v)
        return math.exp(log_gaussian - math.log(v)
                        + stats.invgamma.logpdf(v, nu / 2, scale=nu / 2))

    mixed, error = integrate.quad(integrand, 0, np.inf,
                                  epsabs=2e-11, epsrel=2e-10, limit=300)
    actual = _density(cost, k, nu, -1.0)
    assert error < 1e-8
    assert actual == pytest.approx(mixed, rel=2e-8)
    equivalent = stats.multivariate_t.pdf(
        residual, shape=np.eye(k) * nu / (nu + 2), df=nu + 2)
    assert actual == pytest.approx(equivalent, rel=2e-14)
    assert _radial_integral(lambda c: _density(c, k, nu, -1.0), k) == pytest.approx(1, abs=2e-9)
    # Different relative evidence, so this is not a removable mass-free constant.
    ratio = _density(9.0, k, nu, -1) / _density(0.0, k, nu, -1)
    b_zero_ratio = _density(9.0, k, nu, 0) / _density(0.0, k, nu, 0)
    assert ratio == pytest.approx(b_zero_ratio / (1 + 9 / nu), rel=2e-14)


@pytest.mark.parametrize("b", [-1.0, -0.5, 0.0, 0.5, 1.0])
def test_normalization_at_each_fixed_scale_is_the_additional_condition(b):
    """An already mixed density can be normalized without this condition."""
    v, k = 4.0, 3
    total = _radial_integral(
        lambda cost: v ** b * math.exp(-k / 2 * math.log(2 * math.pi * v) - cost / (2 * v)), k)
    assert total == pytest.approx(v ** b, rel=2e-12)
    assert (abs(total - 1) < 1e-10) == (b == 0)


def _cost(force, acceleration, covariance, mass):
    """Physical least squares, with no mass-law or exponent implementation."""
    k = len(force)
    chol = np.linalg.cholesky(covariance)
    observed = np.linalg.solve(chol, np.r_[force, acceleration])
    design = np.linalg.solve(chol, np.vstack([mass * np.eye(k), np.eye(k)]))
    residual = observed - design @ np.linalg.lstsq(design, observed, rcond=None)[0]
    return float(residual @ residual)


def test_b_minus_one_respects_declared_geometry_units_swap_frames_and_zero():
    rng = np.random.default_rng(1717)
    matrix = rng.normal(size=(6, 6))
    covariance = matrix @ matrix.T + np.eye(6)
    force, acceleration = np.array([2., -.3, .7]), np.array([1., .4, -.2])
    masses, nu, k = [0.2, 1.0, 4.0], 10.0, 3
    lam, mu = 3., 0.7
    units = np.diag([lam] * k + [mu] * k)
    swap = np.r_[np.arange(k, 2 * k), np.arange(k)]
    rotation, _ = np.linalg.qr(rng.normal(size=(k, k)))
    frame = np.kron(np.eye(2), rotation)
    for mass in masses:
        cost = _cost(force, acceleration, covariance, mass)
        transformed = [
            _cost(lam * force, mu * acceleration, units @ covariance @ units, lam / mu * mass),
            _cost(acceleration, force, covariance[np.ix_(swap, swap)], 1 / mass),
            _cost(rotation @ force, rotation @ acceleration, frame @ covariance @ frame.T, mass),
        ]
        assert transformed == pytest.approx([cost] * 3, rel=2e-13)
        for changed in transformed:
            assert _density(changed, k, nu, -1) == pytest.approx(_density(cost, k, nu, -1), rel=2e-13)
        assert _density(_cost(np.zeros(k), np.zeros(k), covariance, mass), k, nu, -1) == _density(0, k, nu, -1)
