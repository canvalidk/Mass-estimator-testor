"""our25: the sandwich-calibrated law (likelihood_product with calibrate: sandwich)."""

import numpy as np
import pytest

from met.analyst import Estimator, joint_log_density, sandwich_weight
from met.law import ReadingSet

CART = {"reference": "cartesian", "nuisance": "radius"}
PRIOR = [{"uniform_angle": {"center": "first_reading"}}]


def _world(rng, n, d, m, snr):
    u = rng.normal(size=(n, d))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    a_star = snr * u
    return ReadingSet((m * a_star + rng.normal(size=(n, d)))[None], (a_star + rng.normal(size=(n, d)))[None],
                      np.ones((1, n)), np.ones((1, n)))


@pytest.mark.parametrize("d,m,snr", [(1, 1.0, 1.0), (1, 4.0, 1.0), (3, 1.0, 1.0), (3, 0.25, 1.0), (3, 4.0, 1.0)])
def test_weight_matches_the_exact_large_n_factor(d, m, snr):
    """For the cartesian law H/J -> r*^2 / (d + r*^2) per reading. One estimate from
    3000 readings has an SD of about 5 percent, so eight series are averaged."""
    rng = np.random.default_rng(100 + d)
    combine = {"rule": "likelihood_product", "prior": PRIOR, "calibrate": "sandwich"}
    w = [sandwich_weight(_world(rng, 3000, d, m, snr), CART, combine)[0] for _ in range(8)]
    r2 = snr**2 * (1 + m * m)
    assert np.mean(w) == pytest.approx(r2 / (d + r2), rel=0.05)


def test_calibrated_density_is_the_tempered_likelihood_times_the_prior():
    readings = _world(np.random.default_rng(7), 30, 3, 1.0, 1.0)
    plain = {"rule": "likelihood_product", "prior": PRIOR}
    calibrated = dict(plain, calibrate="sandwich")
    u = np.linspace(-3, 3, 121)[None]
    lp0, acc0 = joint_log_density(readings, u, CART, plain)
    lp1, acc1 = joint_log_density(readings, u, CART, calibrated)
    w = sandwich_weight(readings, CART, calibrated)[0]
    prior = joint_log_density(readings, u, CART, plain)[0] - np.sum(
        readings.reading_terms(u, "radius", reference="cartesian")[0], axis=1)
    assert 0 < w <= 1
    assert np.allclose(lp1, w * (lp0 - prior) + prior, atol=1e-9)
    assert np.array_equal(acc0, acc1)


def test_calibrate_is_refused_where_it_cannot_apply():
    base = {"name": "x", "reduce": "none", "readouts": ["median"]}
    est = Estimator(dict(base, law=CART, combine={"rule": "likelihood_product", "prior": PRIOR, "calibrate": "sandwich"}))
    with pytest.raises(ValueError, match="at least two"):
        est.check_world({"readings": 1, "design": "new_excitation", "dimension": 3})
    est.check_world({"readings": 5, "design": "new_excitation", "dimension": 3})
    acc = Estimator(dict(base, law={"reference": "cartesian", "nuisance": "acceleration"},
                         combine={"rule": "likelihood_product", "prior": ["flat_mass"], "calibrate": "sandwich"}))
    with pytest.raises(ValueError, match="nuisance radius"):
        acc.check_world({"readings": 5, "design": "new_excitation", "dimension": 3})
    with pytest.raises(ValueError, match="calibrate must be one of"):
        Estimator(dict(base, law=CART, combine={"rule": "likelihood_product", "prior": PRIOR, "calibrate": "magic"}))
    with pytest.raises(ValueError):
        Estimator(dict(base, law=CART, combine={"rule": "single", "prior": [], "calibrate": "sandwich"}))
    with pytest.raises(ValueError, match="sequential"):
        Estimator(dict(base, law=CART, combine={"rule": "sequential", "carry": "exact",
                                                "old": {"rule": "likelihood_product", "prior": PRIOR,
                                                        "calibrate": "sandwich"}}))
