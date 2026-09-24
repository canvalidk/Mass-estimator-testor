"""The sequential rule: old information carried into the prior slot."""

import math

import numpy as np
import pytest

from met.analyst import Estimator, posterior_summaries, tilt_lambda_for_sd, tilt_sd
from met.law import ReadingSet

FLAT = {"reference": "flat", "nuisance": "radius"}
OLD = {"rule": "likelihood_product", "prior": [{"uniform_angle": {"center": "first_reading"}}]}
READOUTS = ["ratio_of_means", "median", "log_sd", "interval_95"]


def _readings(seed, n, b=3, scale=1.0):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(b, n, 3))
    v /= np.linalg.norm(v, axis=2, keepdims=True)
    return ReadingSet(2.0 * v * scale + rng.normal(size=(b, n, 3)), v * scale + 0.6 * rng.normal(size=(b, n, 3)),
                      np.full((b, n), 1.0), np.full((b, n), 0.6))


def _seq(carry):
    return {"rule": "sequential", "old": OLD, "carry": carry}


def test_exact_carry_is_the_rule_applied_to_all_readings():
    r = _readings(1, 6)
    a = posterior_summaries(r, FLAT, OLD, READOUTS)
    b = posterior_summaries(r, FLAT, _seq("exact"), READOUTS)
    for key in READOUTS:
        assert np.allclose(a[key], b[key], rtol=1e-9), key


def test_curve_only_changes_only_the_ratio_of_means():
    r = _readings(2, 6)
    a = posterior_summaries(r, FLAT, _seq("exact"), READOUTS)
    b = posterior_summaries(r, FLAT, _seq("curve_only"), READOUTS)
    for key in ("median", "log_sd", "interval_95"):
        assert np.allclose(a[key], b[key], rtol=1e-9), key
    assert not np.allclose(a["ratio_of_means"], b["ratio_of_means"], rtol=1e-6)


def test_tilt_inversion():
    for value in (-30.0, -4.0, 0.0, 0.5, 4.0, 64.0, 1e3, 1e5, 1e6):
        sd = tilt_sd(value, points=200001)
        back = float(tilt_lambda_for_sd(sd))
        assert back == pytest.approx(value, rel=2e-3, abs=2e-3)
    assert float(tilt_lambda_for_sd(math.pi / 2)) == pytest.approx(0.0, abs=1e-6)


def test_tilt_fit_is_exact_when_the_old_law_is_the_null_law():
    """Zero old readings leave the half-Cauchy at s, which is the tilt family at lambda = 0."""
    rng = np.random.default_rng(3)
    f = np.concatenate([np.zeros((2, 4, 3)), rng.normal(size=(2, 1, 3)) + 2], axis=1)
    a = np.concatenate([np.zeros((2, 4, 3)), rng.normal(size=(2, 1, 3)) + 1], axis=1)
    r = ReadingSet(f, a, np.full((2, 5), 1.3), np.full((2, 5), 0.8))
    exact = posterior_summaries(r, FLAT, _seq("curve_only"), READOUTS)
    tilt = posterior_summaries(r, FLAT, _seq("tilt_fit"), READOUTS)
    for key in READOUTS:
        assert np.allclose(exact[key], tilt[key], rtol=1e-6), key


def test_lognormal_fit_is_exact_for_a_nearly_normal_old_law():
    """With very precise old readings the old law is close to normal in log m."""
    r = _readings(4, 5, scale=200.0)
    exact = posterior_summaries(r, FLAT, _seq("curve_only"), READOUTS)
    fit = posterior_summaries(r, FLAT, _seq("lognormal_fit"), READOUTS)
    assert np.allclose(exact["median"], fit["median"], rtol=1e-4)
    assert np.allclose(exact["log_sd"], fit["log_sd"], rtol=1e-2)


def test_validation():
    base = {"name": "s", "reduce": "none", "law": FLAT, "readouts": ["median"]}
    with pytest.raises(ValueError, match="carry"):
        Estimator(dict(base, combine={"rule": "sequential", "old": OLD, "carry": "points"}))
    with pytest.raises(ValueError, match="prior belongs"):
        Estimator(dict(base, combine={"rule": "sequential", "old": OLD, "carry": "exact", "prior": []}))
    est = Estimator(dict(base, combine=_seq("lognormal_fit")))
    with pytest.raises(ValueError, match="at least two"):
        est.check_world({"readings": 1, "design": "new_excitation"})
    bad = Estimator(dict(base, combine={"rule": "sequential", "carry": "exact",
                                        "old": {"rule": "likelihood_product", "prior": ["flat_mass"]}}))
    with pytest.raises(ValueError, match="old readings"):
        bad.check_world({"readings": 3, "design": "new_excitation"})
