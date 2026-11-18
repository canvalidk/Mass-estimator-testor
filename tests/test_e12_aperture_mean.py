"""Analytic and deterministic checks for the aperture-mean comparison."""

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import norm


def arc_cost(gram, beta):
    low, high = np.linalg.eigvalsh(gram)
    principal = np.linalg.eigh(gram)[1][:, -1]
    angle = math.atan2(principal[1], principal[0]) % math.pi
    distance = min(angle - beta, math.pi - angle) if angle > beta else 0.
    return low + (high - low) * math.sin(distance) ** 2


@pytest.mark.parametrize("beta", [math.pi / 3, .4 * math.pi, math.pi / 2, .8 * math.pi])
def test_reflection_bounds_mean_pointwise_for_every_pair_of_eigenvalues(beta):
    reflection = np.diag([1., -1.])
    for angle in np.linspace(0., math.pi, 169):
        u = np.array([math.cos(angle), math.sin(angle)])
        for low, high in [(0., 0.), (0., 9.), (.7, 2.4), (5., 5.)]:
            gram = low * np.eye(2) + (high - low) * np.outer(u, u)
            reflected = reflection @ gram @ reflection
            left = arc_cost(gram, beta) + arc_cost(reflected, beta)
            assert left <= low + gram[1, 1] + 1e-13


def test_reflection_inequality_has_a_sharp_geometric_threshold():
    beta = math.pi / 3 - .025
    angle = (math.pi - beta) / 2
    u = np.array([math.cos(angle), math.sin(angle)])
    gram = np.outer(u, u)
    reflection = np.diag([1., -1.])
    total = arc_cost(gram, beta) + arc_cost(reflection @ gram @ reflection, beta)
    assert total > gram[1, 1] + 1e-5


def gain_lower_bound(h, strength):
    r = strength
    return (norm.pdf(r + 2 / r) / r**3
            * norm.cdf((h * (r + 1 / r) - 1) / math.sqrt(1 - h**2)))


def loss_upper_bound(h, strength):
    v = strength / math.sqrt(2 * (1 - h))
    return norm.sf(v) + (1 - h) * v * norm.pdf(v) / 2


@pytest.mark.parametrize("h,strength", [(.9, 2.1), (.75, 7.), (.65, 10.), (.55, 18.)])
def test_gain_exceeds_even_the_total_loss_bound(h, strength):
    gain = gain_lower_bound(h, strength)
    loss = loss_upper_bound(h, strength)
    assert gain > 100 * loss > 0
    if h == .75:
        assert gain / loss > 1e7
    # Independently integrate the conditional second moment on the losing half-space.
    tau = math.sqrt(2 * (1 - h))
    v = strength / tau
    integral = quad(lambda z: norm.pdf(z) * ((tau * z) ** 2 / 4 + (1 + h) / 2),
                    v, v + 10, epsabs=loss * 1e-10, epsrel=1e-10)[0]
    assert integral == pytest.approx(loss, rel=1e-9, abs=0.)


def scalar_mean_by_conditioning(rho, force_mean, accel_mean):
    """Integrate min(X²,Y²) on opposite signs in physical channel coordinates."""
    sd = math.sqrt(1 - rho**2)

    def integrand(x):
        a = abs(x)
        # Flip both channels on the negative half so eligible Y is always negative.
        mu = (accel_mean + rho * (x - force_mean)) * np.sign(x)
        bottom, top = (-a - mu) / sd, -mu / sd
        p = norm.cdf(top) - norm.cdf(bottom)
        lower_density, upper_density = norm.pdf(bottom), norm.pdf(top)
        square = ((mu**2 + sd**2) * p + 2 * mu * sd * (lower_density - upper_density)
                  + sd**2 * (bottom * lower_density - top * upper_density))
        return norm.pdf(x - force_mean) * (square + a**2 * norm.cdf(bottom))

    return sum(quad(integrand, low, high, epsabs=2e-11, epsrel=2e-11)[0]
               for low, high in [(-12., 0.), (0., force_mean + 12.)])


def test_strict_excess_survives_an_interior_mass_and_has_an_observable_margin():
    # The boundary proof has an independent quantitative lower bound.
    margin = gain_lower_bound(.9, 2.1) - loss_upper_bound(.9, 2.1)
    assert scalar_mean_by_conditioning(-.9, 2.1, 0.) > .5 + margin
    assert scalar_mean_by_conditioning(-.9, 2.1, .001) > .5


@pytest.mark.parametrize("rho", [-.5, -.25, 0., .5])
def test_boundary_scalar_means_stay_below_half_on_the_wider_apertures(rho):
    for strength in [0., 1., 2.5]:
        assert scalar_mean_by_conditioning(rho, strength, 0.) < .5
