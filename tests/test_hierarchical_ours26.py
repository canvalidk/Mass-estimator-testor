"""The hierarchical rule under ours26 (27 September): a learned push scale in noise units.

Every standardised push, each in its own noise pair's units, gets one shared prior
N(0, omega^2 I), beta = omega^2 / (1 + omega^2). Integrating the pushes out raises
ours26's likelihood to the power beta (the tempering identity), so the peak stays where
ours26 has it and the law widens by 1/beta. The readout is the stacked one: given theta
and beta the stacked push is N(beta W, beta I).

Checks, each against a computation done a different way: the tempering identity; the
shrunk stacked weight against Monte Carlo, and the two implementations of it against each
other; the beta average against direct quadrature over beta; beta -> 1 gives ours26. Then
what should hold: a vector series and its axes agree, the reciprocal symmetry, units,
and the refusals.
"""

import math

import numpy as np
import pytest
from scipy.integrate import quad

from met.analyst import (Estimator, _hierarchical_rss_alpha, joint_log_density, posterior_summaries,
                         prior_log_density)
from met.law import ReadingSet

CART1 = {"reference": "cartesian1", "nuisance": "radius"}
PRIOR = [{"uniform_angle": {"center": "first_reading_geometric"}}]
READOUTS = ["ratio_of_means_rss", "median", "interval_95"]
SF, SA = np.array([0.9, 0.4, 1.5]), np.array([0.5, 0.5, 0.2])      # three different ratios


def _hier(b=1.0, fixed=None):
    out = {"rule": "hierarchical", "prior": PRIOR, "excitation": {"beta_b": b}}
    if fixed is not None:
        out["_fixed_beta"] = fixed
    return out


FLAT = {"rule": "likelihood_product", "prior": PRIOR}


def _per_axis(seed, n, d=3, scale=1.0):
    rng = np.random.default_rng(seed)
    f = rng.normal(1.0, 1.0, (1, n, d)) * scale
    a = rng.normal(0.4, 0.6, (1, n, d)) * scale
    return ReadingSet(f, a, np.broadcast_to(SF[:d], (1, n, d)), np.broadcast_to(SA[:d], (1, n, d)))


def _isotropic(seed, n, d=3, sf=0.7, sa=0.3):
    rng = np.random.default_rng(seed)
    f, a = rng.normal(1.0, 1.0, (1, n, d)), rng.normal(0.4, 0.6, (1, n, d))
    return ReadingSet(f, a, np.full((1, n), sf), np.full((1, n), sa))


def _angles(readings, u):
    """Per-component theta_g at log masses u (S,), with the components' standardised data."""
    sf, sa = readings.force_sd_axes[0].ravel(), readings.acceleration_sd_axes[0].ravel()
    x = (readings.force[0] / readings.force_sd_axes[0]).ravel()
    y = (readings.acceleration[0] / readings.acceleration_sd_axes[0]).ravel()
    th = np.arctan(np.exp(u)[:, None] / (sf / sa)[None, :])          # (S, K)
    return th, x, y, sa


# ---------------------------------------------------------------- the law

@pytest.mark.parametrize("beta", [0.2, 0.7, 0.97])
def test_tempering_identity(beta):
    """With beta fixed, the law's log-likelihood is beta times ours26's, up to a constant."""
    readings = _per_axis(1, 4)
    u = np.linspace(-2.5, 2.5, 41)[None]
    prior = prior_log_density(PRIOR, u, readings)
    lp_h, _ = joint_log_density(readings, u, CART1, _hier(fixed=beta))
    lp_f, _ = joint_log_density(readings, u, CART1, FLAT)
    diff = (lp_h - prior) - beta * (lp_f - prior)
    assert np.ptp(diff) < 1e-10


def test_peak_does_not_move():
    readings = _per_axis(2, 6)
    u = np.linspace(-1.5, 2.5, 4001)[None]
    prior = prior_log_density(PRIOR, u, readings)
    flat = np.argmax(joint_log_density(readings, u, CART1, FLAT)[0] - prior)
    for beta in (0.05, 0.5, 0.99):
        lp, _ = joint_log_density(readings, u, CART1, _hier(fixed=beta))
        assert np.argmax(lp - prior) == flat


def test_law_matches_brute_force_per_axis():
    """Two 2D readings with per-axis noise: integrate every standardised push under
    N(0, omega^2) and beta ~ Beta(1, b) numerically, component by component."""
    readings = _per_axis(3, 2, d=2)
    b = 1.4
    us = np.log(np.array([0.3, 0.9, 2.0, 5.0]))
    lp, _ = joint_log_density(readings, us[None], CART1, _hier(b))
    th, x, y, _ = _angles(readings, us)

    def marginal(row):
        def given_beta(beta):
            omega2 = beta / (1 - beta)
            total = 1.0
            for xi, yi, t in zip(x, y, th[row]):
                s, c = math.sin(t), math.cos(t)
                f = lambda v: (math.exp(-v * v / (2 * omega2) - (xi - v * s) ** 2 / 2 - (yi - v * c) ** 2 / 2)
                               / math.sqrt(2 * math.pi * omega2))
                total *= quad(f, -60, 60, epsabs=0, epsrel=1e-11)[0]
            return total * b * (1 - beta) ** (b - 1)
        return math.log(quad(given_beta, 1e-9, 1 - 1e-9, epsabs=0, epsrel=1e-10, limit=200)[0])

    brute = np.array([marginal(i) for i in range(len(us))])
    prior = prior_log_density(PRIOR, us[None], readings)[0]
    assert np.ptp(brute + prior - lp[0]) < 1e-6


# ---------------------------------------------------------------- the stacked weight

@pytest.mark.parametrize("beta", [0.15, 0.6])
def test_shrunk_weight_against_monte_carlo(beta):
    readings = _per_axis(4, 3)
    u = np.array([[math.log(1.3)]])
    got = readings.stacked_cond_alpha_at(u, np.array([[beta]]))[0, 0]
    th, x, y, sa = _angles(readings, u[0])
    w = x * np.sin(th[0]) + y * np.cos(th[0])
    rng = np.random.default_rng(9)
    v = beta * w + math.sqrt(beta) * rng.normal(size=(1_000_000, w.size))
    mc = np.linalg.norm(v * (sa * np.cos(th[0]))[None, :], axis=1)
    assert got == pytest.approx(mc.mean(), abs=4 * mc.std() / math.sqrt(mc.size))


def test_isotropic_shortcut_equals_weighted_path():
    readings = _isotropic(5, 4)
    u = np.linspace(-1.0, 2.0, 7)[None]
    beta = np.linspace(0.1, 0.95, 7)[None]
    fast = readings.stacked_cond_alpha_at(u, beta)
    slow = readings._weighted_cond_alpha_at(u, beta)
    assert np.allclose(fast, slow, rtol=1e-9, atol=0)


@pytest.mark.parametrize("make", [_per_axis, _isotropic])
@pytest.mark.parametrize("b", [1.0, 3.0])
def test_beta_average_against_quadrature(make, b):
    """E|a*_stack| given theta, averaged over beta | theta by tanh-sinh nodes in its quantile,
    against direct quadrature over beta."""
    readings = make(6, 3)
    us = np.log(np.array([0.4, 1.1, 3.0]))
    got = _hierarchical_rss_alpha(readings, us[None], _hier(b))[0]
    like, _, _ = readings.reading_terms(us[None], "radius", reference="cartesian1")
    z = np.sum(like, axis=1)[0]
    c = 0.5 * readings.readings * readings.dimension + b
    for i, ui in enumerate(us):
        w = lambda beta: (1 - beta) ** (c - 1) * math.exp((beta - 1) * z[i])
        norm = quad(w, 0, 1, epsabs=0, epsrel=1e-12)[0]
        val = quad(lambda beta: w(beta) * readings.stacked_cond_alpha_at(np.array([[ui]]), np.array([[beta]]))[0, 0],
                   0, 1, epsabs=0, epsrel=1e-10)[0] / norm
        assert got[i] == pytest.approx(val, rel=1e-9)


# ---------------------------------------------------------------- whole estimator

def test_beta_to_one_is_ours26():
    readings = _per_axis(7, 5)
    flat = posterior_summaries(readings, CART1, FLAT, READOUTS)
    near = posterior_summaries(readings, CART1, _hier(fixed=1 - 1e-12), READOUTS)
    for key in READOUTS:
        assert np.allclose(near[key], flat[key], rtol=1e-8, atol=0), key


@pytest.mark.parametrize("n", [1, 4])
def test_axis_splitting(n):
    readings = _per_axis(8, n)
    whole = posterior_summaries(readings, CART1, _hier(), READOUTS)
    axes = posterior_summaries(readings.split_axes(), CART1, _hier(), READOUTS)
    for key in READOUTS:
        assert np.allclose(whole[key], axes[key], rtol=1e-8, atol=0), key


def test_learned_scale_widens_but_keeps_the_centre_close():
    readings = _per_axis(10, 30, scale=0.6)
    flat = posterior_summaries(readings, CART1, FLAT, ["median", "log_sd"])
    hier = posterior_summaries(readings, CART1, _hier(), ["median", "log_sd"])
    assert hier["log_sd"][0] > 1.2 * flat["log_sd"][0]


def test_reciprocal_symmetry():
    readings = _per_axis(11, 4)
    swapped = ReadingSet(readings.acceleration, readings.force,
                         readings.acceleration_sd_axes, readings.force_sd_axes)
    one = posterior_summaries(readings, CART1, _hier(), ["ratio_of_means_rss", "median"])
    two = posterior_summaries(swapped, CART1, _hier(), ["ratio_of_means_rss", "median"])
    for key in ("ratio_of_means_rss", "median"):
        assert one[key][0] * two[key][0] == pytest.approx(1.0, rel=2e-6), key


def test_units():
    readings = _per_axis(12, 4)
    lam = 3.7
    scaled = ReadingSet(lam * readings.force, readings.acceleration,
                        lam * readings.force_sd_axes, readings.acceleration_sd_axes)
    one = posterior_summaries(readings, CART1, _hier(), ["ratio_of_means_rss", "median"])
    two = posterior_summaries(scaled, CART1, _hier(), ["ratio_of_means_rss", "median"])
    for key in ("ratio_of_means_rss", "median"):
        assert two[key][0] == pytest.approx(lam * one[key][0], rel=2e-6), key


# ---------------------------------------------------------------- settings

def test_validation():
    world = {"readings": 5, "design": "new_excitation", "dimension": 3,
             "noise": {"force_sd": {"per_axis": [1, 1, 1]}, "acceleration_sd": {"per_axis": [0.2, 0.5, 1.25]}}}
    base = {"name": "h", "reduce": "none", "law": CART1, "combine": _hier()}
    Estimator(dict(base, readouts=READOUTS)).check_world(world)
    Estimator(dict(base, reduce="split_axes", readouts=READOUTS)).check_world(world)
    with pytest.raises(ValueError, match="plain sum"):
        Estimator(dict(base, readouts=["ratio_of_means"])).check_world(world)
    with pytest.raises(ValueError, match="first_reading"):
        combine = {"rule": "hierarchical", "prior": [{"uniform_angle": {"center": "first_reading"}}],
                   "excitation": {"beta_b": 1}}
        Estimator(dict(base, combine=combine, readouts=READOUTS)).check_world(world)
