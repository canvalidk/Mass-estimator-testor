"""The hierarchical excitation rule: a shared Gaussian prior on every reading's excitation."""

import math

import numpy as np
import pytest
from scipy.integrate import quad

from met.analyst import (Estimator, _log_lower_gamma_integral, joint_log_density, posterior_summaries)
from met.law import ReadingSet, cartesian1_mean_radius
from met.world import generate, series_rngs

CART = {"reference": "cartesian1", "nuisance": "radius"}
PRIOR = [{"uniform_angle": {"center": "first_reading"}}]


def _hier(b=1.0):
    return {"rule": "hierarchical", "prior": PRIOR, "excitation": {"beta_b": b}}


@pytest.mark.parametrize("c", [0.8, 3.5, 76.0])
@pytest.mark.parametrize("z", [0.0, 0.3, 2.0, 40.0, 75.0, 77.0, 900.0])
def test_beta_integral_closed_form(c, z):
    ref = quad(lambda g: g ** (c - 1) * math.exp(-g * z), 0, 1, epsabs=0, epsrel=1e-13, limit=200)[0]
    assert float(_log_lower_gamma_integral(c, np.array([z]))[0]) == pytest.approx(math.log(ref), abs=1e-9)


def test_law_against_brute_force():
    """1D, two readings: integrate v_i under N(0, omega^2) and beta ~ Beta(1, b) numerically."""
    rng = np.random.default_rng(2)
    x, y = rng.normal(size=2) * 1.5 + 1, rng.normal(size=2) + 0.5
    b = 1.7
    readings = ReadingSet(x[None, :, None], y[None, :, None], np.ones((1, 2)), np.ones((1, 2)))
    thetas = np.linspace(0.15, 1.35, 6)
    u = np.log(np.tan(thetas))[None]
    lp, _ = joint_log_density(readings, u, CART, _hier(b))

    def marginal(t):
        s, c = math.sin(t), math.cos(t)

        def given_beta(beta):
            omega2 = beta / (1 - beta)
            total = 1.0
            for xi, yi in zip(x, y):
                f = lambda v: (math.exp(-v * v / (2 * omega2) - (xi - v * s) ** 2 / 2 - (yi - v * c) ** 2 / 2)
                               / math.sqrt(2 * math.pi * omega2))
                total *= quad(f, -60, 60, epsabs=0, epsrel=1e-11)[0]
            return total * b * (1 - beta) ** (b - 1)
        return math.log(quad(given_beta, 1e-9, 1 - 1e-9, epsabs=0, epsrel=1e-9, limit=200)[0])

    brute = np.array([marginal(t) for t in thetas])
    prior = -np.log(np.cosh(u[0]))   # uniform-angle factor at s = 1, up to a constant
    diff = brute + prior - lp[0]
    assert np.ptp(diff) < 1e-6


def test_conditional_acceleration_quadrature():
    """E[alpha | theta] from 24 quantile nodes against direct integration over beta."""
    rng = np.random.default_rng(3)
    f, a = rng.normal(size=(1, 4, 3)) * 2, rng.normal(size=(1, 4, 3))
    readings = ReadingSet(f, a, np.full((1, 4), 1.0), np.full((1, 4), 1.0))
    theta = 0.7
    u = np.array([[math.log(math.tan(theta))]])
    _, acc = joint_log_density(readings, u, CART, _hier(1.0))
    x, y = readings.force[0], readings.acceleration[0]
    h = np.linalg.norm(x * math.sin(theta) + y * math.cos(theta), axis=1)
    n, c = 12, 12 / 2 + 1.0
    z = 0.5 * np.sum(h * h)
    w = lambda beta: (1 - beta) ** (c - 1) * math.exp(beta * z - z)
    norm = quad(w, 0, 1, epsabs=0, epsrel=1e-12)[0]
    val = quad(lambda beta: w(beta) * sum(math.sqrt(beta) * float(cartesian1_mean_radius(math.sqrt(beta) * hi, 3))
                                          for hi in h), 0, 1, epsabs=0, epsrel=1e-12)[0] / norm
    assert acc[0, 0] == pytest.approx(math.cos(theta) * val, rel=2e-3)


def test_strong_excitation_approaches_the_flat_excitation_prior():
    """With excitation far above the noise, beta -> 1 and the law approaches the cartesian product."""
    rng = np.random.default_rng(4)
    v = rng.normal(size=(1, 20, 3)) * 30
    readings = ReadingSet(2 * v + rng.normal(size=v.shape), v + rng.normal(size=v.shape),
                          np.ones((1, 20)), np.ones((1, 20)))
    flat = posterior_summaries(readings, CART, {"rule": "likelihood_product", "prior": PRIOR},
                               ["median", "log_sd"])
    hier = posterior_summaries(readings, CART, _hier(1.0), ["median", "log_sd"])
    assert hier["median"][0] == pytest.approx(flat["median"][0], rel=1e-3)
    assert hier["log_sd"][0] == pytest.approx(flat["log_sd"][0], rel=2e-2)


def test_weak_excitation_widens_the_law():
    rng = np.random.default_rng(5)
    v = rng.normal(size=(1, 30, 3)) * 0.8
    readings = ReadingSet(v + rng.normal(size=v.shape), v + rng.normal(size=v.shape),
                          np.ones((1, 30)), np.ones((1, 30)))
    flat = posterior_summaries(readings, CART, {"rule": "likelihood_product", "prior": PRIOR}, ["log_sd"])
    hier = posterior_summaries(readings, CART, _hier(1.0), ["log_sd"])
    assert hier["log_sd"][0] > 1.2 * flat["log_sd"][0]


def test_validation():
    base = {"name": "h", "reduce": "none", "readouts": ["median"]}
    with pytest.raises(ValueError, match="cartesian1"):
        Estimator(dict(base, law={"reference": "flat", "nuisance": "radius"}, combine=_hier()))
    with pytest.raises(ValueError, match="beta_b"):
        Estimator(dict(base, law=CART, combine={"rule": "hierarchical", "prior": PRIOR, "excitation": {"beta_b": 0}}))
    est = Estimator(dict(base, law=CART, combine=_hier()))
    est.check_world({"readings": 50, "design": "new_excitation", "dimension": 3})


def test_gaussian_excitation_world():
    world = {"dimension": 3, "design": "new_excitation", "mass": 2.0, "acceleration_snr": {"normal_rms": 3.0},
             "direction": "random", "readings": 50, "noise": {"force_sd": 1.0, "acceleration_sd": 0.5},
             "supplied_noise": "exact"}
    _, truth = generate(world, series_rngs(1, world, 200))
    rms = math.sqrt(np.mean(np.sum(truth["true_acceleration"] ** 2, axis=2)))
    assert rms == pytest.approx(3.0 * 0.5, rel=0.03)
    from met.world import validate_cell
    with pytest.raises(ValueError, match="direction: random"):
        validate_cell(dict(world, direction="fixed"))
