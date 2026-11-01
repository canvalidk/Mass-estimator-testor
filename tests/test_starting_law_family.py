"""B08: algebra of admissible starting densities, not a general-prior met API.

met.full.ours30 accepts only a centre c; _log_prior hard-codes the half-Cauchy.
These tests therefore check the three requested members as probability densities
and the all-zero Gaussian update algebra independently, without monkeypatching or
claiming that production can accept a log-normal or log-logistic starting law.
"""
import math

import numpy as np
import pytest
from scipy import integrate, stats


FAMILIES = ("half_cauchy", "log_normal", "log_logistic")
CENTRES = (0.3, 1.0, 7.0)


def _log_shape(family, x):
    """Log of a nonnegative, normalized, even density g of x=log(m/c)."""
    if family == "half_cauchy":
        # 1 / (pi cosh x), evaluated without overflowing cosh.
        ax = np.abs(x)
        return math.log(2 / math.pi) - ax - np.logaddexp(0.0, -2 * ax)
    if family == "log_normal":
        return stats.norm.logpdf(x, scale=0.8)
    if family == "log_logistic":
        return stats.logistic.logpdf(x, scale=1.3)
    raise ValueError(family)


def _density(family, m, centre):
    x = np.log(m) - math.log(centre)
    return np.exp(_log_shape(family, x) - np.log(m))


def _reference(family, centre):
    if family == "half_cauchy":
        return stats.halfcauchy(scale=centre)
    if family == "log_normal":
        return stats.lognorm(s=0.8, scale=centre)
    return stats.fisk(c=1 / 1.3, scale=centre)


@pytest.mark.parametrize("family", FAMILIES)
@pytest.mark.parametrize("centre", CENTRES)
def test_starting_density_is_normalized_log_pushforward_with_centre_median(family, centre):
    """The 1/m Jacobian gives the standard three physical-mass distributions."""
    x = np.linspace(-8, 8, 101)
    m = centre * np.exp(x)
    density = _density(family, m, centre)
    np.testing.assert_allclose(density, _reference(family, centre).pdf(m), rtol=1e-13)
    np.testing.assert_allclose(m * density, np.exp(_log_shape(family, x)), rtol=1e-13)
    np.testing.assert_allclose(_log_shape(family, x), _log_shape(family, -x), atol=2e-14)

    total = integrate.quad(lambda mass: _density(family, mass, centre), 0, np.inf,
                           epsabs=1e-10, epsrel=1e-10)[0]
    below_c = integrate.quad(lambda mass: _density(family, mass, centre), 0, centre,
                             epsabs=1e-11, epsrel=1e-11)[0]
    assert total == pytest.approx(1.0, abs=2e-10)
    assert below_c == pytest.approx(0.5, abs=2e-11)


@pytest.mark.parametrize("family", FAMILIES)
@pytest.mark.parametrize("centre", CENTRES)
def test_starting_density_obeys_units_and_reciprocal_jacobians(family, centre):
    m = centre * np.exp(np.linspace(-5, 5, 37))
    density = _density(family, m, centre)
    # Image under m -> c^2/m includes the reciprocal map's density Jacobian.
    reciprocal_image = _density(family, centre**2 / m, centre) * centre**2 / m**2
    np.testing.assert_allclose(reciprocal_image, density, rtol=2e-13)
    for force_units, acceleration_units in [(0.2, 3.0), (4.0, 0.3)]:
        scale = force_units / acceleration_units
        np.testing.assert_allclose(
            _density(family, scale * m, scale * centre), density / scale, rtol=2e-13)


def _null_gaussian_factor(m):
    """Integrate one zero reading in noise volume, independently of met.full.

    In physical latent acceleration alpha, the zero-data Gaussian kernel is
    exp(-G(m) alpha^2/2). Noise volume dv=sqrt(G(m)) d alpha cancels its scale.
    The omitted Gaussian likelihood normalizer is constant in candidate mass.
    """
    covariance = np.array([[1.6, 0.4], [0.4, 0.9]])
    candidate = np.array([m, 1.0])
    gram = candidate @ np.linalg.solve(covariance, candidate)
    physical_integral = integrate.quad(
        lambda alpha: math.exp(-0.5 * gram * alpha**2), -np.inf, np.inf,
        epsabs=1e-11, epsrel=1e-11)[0]
    return math.sqrt(gram) * physical_integral


@pytest.mark.parametrize("family", FAMILIES)
@pytest.mark.parametrize("centre", CENTRES)
def test_all_zero_starting_family_algebra_returns_pi_and_median_c(family, centre):
    """Constant null evidence leaves each pi unchanged; this is not an API test."""
    masses = centre * np.array([0.125, 0.5, 1.0, 2.0, 8.0])
    factors = np.array([_null_gaussian_factor(m) for m in masses])
    np.testing.assert_allclose(factors, math.sqrt(2 * math.pi), rtol=1e-12)
    # Any number of independent null readings multiplies only a common constant.
    for readings in [0, 1, 4]:
        null_evidence = factors**readings
        normalizer = math.sqrt(2 * math.pi)**readings
        updated = _density(family, masses, centre) * null_evidence / normalizer
        np.testing.assert_allclose(updated, _density(family, masses, centre), rtol=4e-12)
    # Normalize the physical density, then verify its actual integrated median.
    total = integrate.quad(lambda m: _density(family, m, centre), 0, np.inf,
                           epsabs=1e-10, epsrel=1e-10)[0]
    null_cdf_at_c = integrate.quad(lambda m: _density(family, m, centre) / total,
                                   0, centre, epsabs=1e-11, epsrel=1e-11)[0]
    assert null_cdf_at_c == pytest.approx(0.5, abs=2e-10)


@pytest.mark.parametrize("centre", CENTRES)
def test_half_cauchy_is_uniform_in_angle(centre):
    theta = np.linspace(0.01, math.pi / 2 - 0.01, 73)
    mass = centre * np.tan(theta)
    jacobian = centre / np.cos(theta)**2
    np.testing.assert_allclose(
        _density("half_cauchy", mass, centre) * jacobian, 2 / math.pi, rtol=1e-13)


def test_even_proper_log_density_can_have_nonunique_mass_median():
    """L3/L4/properness give c as a median, but do not force uniqueness."""
    # Equal uniform pieces on [-2,-1] and [1,2] leave a central gap.
    def log_cdf(x):
        return 0.5 * np.clip(x + 2, 0, 1) + 0.5 * np.clip(x - 1, 0, 1)

    assert log_cdf(-np.inf) == 0
    assert log_cdf(np.inf) == 1
    # In mass coordinates every point in [c/e, c*e] has probability 1/2 below it.
    centre = 3.0
    medians = centre * np.exp(np.array([-1.0, -0.4, 0.0, 0.6, 1.0]))
    np.testing.assert_allclose(log_cdf(np.log(medians / centre)), 0.5, atol=1e-15)
    assert log_cdf(-1 - 1e-6) < 0.5
    # The lower CDF-quantile convention chooses c/e, not c, for this law.
    assert centre * math.exp(-1) != centre


@pytest.mark.parametrize("log_atoms,weights", [([0.0], [1.0]),
                                              ([-math.log(2), math.log(2)], [0.5, 0.5])],
                         ids=["point_mass", "reciprocal_pair"])
def test_reciprocal_atomic_laws_are_admissible_without_density(log_atoms, weights):
    """The declared law type does not require absolute continuity in dm.

    These symmetric log-laws obey properness, units covariance and the reciprocal
    swap, but positive mass on singletons rules out an ordinary density g.
    """
    centre = 1.7
    atoms = centre * np.exp(log_atoms)
    weights = np.asarray(weights)
    assert np.all(atoms > 0)
    assert np.all(weights > 0)
    assert weights.sum() == 1
    image_order = np.argsort(centre**2 / atoms)
    np.testing.assert_allclose((centre**2 / atoms)[image_order], atoms, rtol=1e-14)
    np.testing.assert_array_equal(weights[image_order], weights)
    assert weights[atoms <= centre].sum() >= 0.5
    assert weights[atoms >= centre].sum() >= 0.5
    # Positive rescaling acts on the whole law and its centre together.
    scale = 0.3 / 4.0
    np.testing.assert_allclose(np.log((scale * atoms) / (scale * centre)), log_atoms)
    assert weights[atoms == atoms[0]].sum() > 0  # a singleton has positive probability


@pytest.mark.parametrize("omega", [0.2, 1.0, 4.0])
def test_ours28_gaussian_size_marginal_fisher_is_theta_omega_orthogonal(omega):
    """Check the marginal Gaussian size model, not the fixed-signal (theta,r) model.

    In whitened coordinates, alpha~N(0,omega^2) gives covariance
    V=I+omega^2 u u^T. The standard Gaussian trace formula gives the Fisher
    matrix after integrating alpha out. This is not a reference-prior assertion.
    """
    def covariance(theta, size):
        direction = np.array([math.sin(theta), math.cos(theta)])
        return np.eye(2) + size**2 * np.outer(direction, direction)

    expected = np.diag([omega**4 / (1 + omega**2),
                        2 * omega**2 / (1 + omega**2)**2])
    for theta in [0.17, 0.81, 1.3]:
        htheta = 1e-5
        homega = 1e-5 * max(1.0, omega)
        derivatives = [
            (covariance(theta + htheta, omega) - covariance(theta - htheta, omega)) /
            (2 * htheta),
            (covariance(theta, omega + homega) - covariance(theta, omega - homega)) /
            (2 * homega),
        ]
        inverse_products = [np.linalg.solve(covariance(theta, omega), derivative)
                            for derivative in derivatives]
        fisher = np.array([[0.5 * np.trace(a @ b) for b in inverse_products]
                           for a in inverse_products])
        np.testing.assert_allclose(fisher, expected, rtol=2e-8, atol=2e-10)
