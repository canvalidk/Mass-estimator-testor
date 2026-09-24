"""The cartesian reference: flat in the components of the latent vector."""

import math

import numpy as np
import pytest
from scipy.integrate import quad, dblquad

from met.analyst import estimate, posterior_summaries
from met.law import ReadingSet, cartesian_mean_radius

READOUTS = ["ratio_of_means", "median", "log_sd", "interval_95"]
CART = {"reference": "cartesian", "nuisance": "radius"}
FLAT = {"reference": "flat", "nuisance": "radius"}


@pytest.mark.parametrize("d", [1, 2, 3])
@pytest.mark.parametrize("h", [0.0, 0.4, 1.5, 4.0])
def test_mean_radius_is_the_noncentral_chi_mean(d, h):
    """E|v| for v ~ N(w, I_d), |w| = h, by direct integration in polar form about w's axis."""
    if d == 1:
        ref = quad(lambda v: abs(v) * math.exp(-(v - h) ** 2 / 2), -np.inf, np.inf)[0] / math.sqrt(2 * math.pi)
    elif d == 2:
        ref = dblquad(lambda phi, r: r * r * math.exp(-(r * r - 2 * r * h * math.cos(phi) + h * h) / 2),
                      0, 12 + h, 0, 2 * math.pi)[0] / (2 * math.pi)
    else:
        ref = dblquad(lambda c, r: r ** 3 * math.exp(-(r * r - 2 * r * h * c + h * h) / 2),
                      0, 12 + h, -1, 1)[0] * 2 * math.pi / (2 * math.pi) ** 1.5
    assert float(cartesian_mean_radius(h, d)) == pytest.approx(ref, rel=1e-7)


def test_cartesian_equals_flat_in_two_dimensions():
    for seed in range(3):
        rng = np.random.default_rng(seed)
        f, a = rng.normal(size=2) * 3, rng.normal(size=2)
        law_c = estimate(f, a, 0.8, 1.1, READOUTS, law=CART)
        law_f = estimate(f, a, 0.8, 1.1, READOUTS, law=FLAT)
        for key in READOUTS:
            assert np.allclose(law_c[key], law_f[key], rtol=1e-10), key


@pytest.mark.parametrize("d", [1, 3])
def test_cartesian_differs_from_flat_outside_two_dimensions(d):
    f, a = np.full(d, 1.0), np.full(d, 0.6)
    assert estimate(f, a, 1, 1, law=CART)["ratio_of_means"] != pytest.approx(
        estimate(f, a, 1, 1, law=FLAT)["ratio_of_means"], rel=1e-6)


@pytest.mark.parametrize("d", [1, 2, 3])
def test_null_law_is_still_the_half_cauchy(d):
    r = estimate(np.zeros(d), np.zeros(d), 2.0, 0.5, READOUTS, law=CART)
    assert r["ratio_of_means"] == pytest.approx(4.0, rel=1e-12)
    assert r["log_sd"] == pytest.approx(math.pi / 2, rel=1e-12)


@pytest.mark.parametrize("d", [1, 3])
def test_single_reading_law_by_brute_force(d):
    """Integrate exp(-|x - v sin t|^2/2 - |y - v cos t|^2/2) over v in R^d directly, on a theta grid."""
    rng = np.random.default_rng(40 + d)
    x, y = rng.normal(size=d) * 2 + 1, rng.normal(size=d) + 0.5
    theta = np.linspace(0.05, math.pi / 2 - 0.05, 7)
    brute = []
    for t in theta:
        w = x * math.sin(t) + y * math.cos(t)
        # the v-integral factorises over components; do each 1D integral numerically
        val = 1.0
        for k in range(d):
            val *= quad(lambda v: math.exp(-(x[k] - v * math.sin(t)) ** 2 / 2 - (y[k] - v * math.cos(t)) ** 2 / 2),
                        -np.inf, np.inf)[0]
        brute.append(math.log(val))
    readings = ReadingSet(x[None, None], y[None, None], [[1.0]], [[1.0]])
    u = np.log(np.tan(theta))[None]
    like, _, _ = readings.reading_terms(u, "radius", "cartesian")
    diff = np.array(brute) - like[0, 0]
    assert np.ptp(diff) < 1e-8


@pytest.mark.parametrize("d", [1, 3])
def test_vonmises_state_is_exact_under_cartesian_in_every_dimension(d):
    rng = np.random.default_rng(50 + d)
    n, b = 7, 3
    v = rng.normal(size=(b, n, d))
    v /= np.linalg.norm(v, axis=2, keepdims=True)
    r = ReadingSet(2.0 * v + rng.normal(size=(b, n, d)), v + 0.6 * rng.normal(size=(b, n, d)),
                   np.full((b, n), 1.0), np.full((b, n), 0.6))
    old = {"rule": "likelihood_product", "prior": [{"uniform_angle": {"center": "first_reading"}}]}
    a = posterior_summaries(r, CART, {"rule": "sequential", "old": old, "carry": "curve_only"}, READOUTS)
    b_ = posterior_summaries(r, CART, {"rule": "sequential", "old": old, "carry": "vonmises_state"}, READOUTS)
    for key in READOUTS:
        assert np.allclose(a[key], b_[key], rtol=1e-9), key
