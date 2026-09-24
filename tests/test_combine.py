"""Combination rules and priors: identities, contribution 11/12 values, improper laws."""

import math

import numpy as np
import pytest

from met.analyst import joint_log_density, posterior_summaries
from met.law import ReadingSet

FLAT = {"reference": "flat", "nuisance": "radius"}
SYMMETRIC = {"rule": "likelihood_product", "prior": [{"uniform_angle": {"center": "first_reading"}}]}
EQ18 = {"rule": "posterior_product", "coordinate": "mass", "prior": []}
EQ20 = {"rule": "posterior_product", "coordinate": "log_mass", "prior": []}
READOUTS = ["ratio_of_means", "median", "log_sd", "interval_95"]


def _series(force, acceleration, sf, sa):
    force, acceleration = np.asarray(force, float), np.asarray(acceleration, float)
    n = force.shape[0]
    return ReadingSet(force[None], acceleration[None], np.broadcast_to(sf, (1, n)), np.broadcast_to(sa, (1, n)))


def _same(a, b, rtol=1e-9):
    for key in READOUTS:
        assert np.allclose(a[key], b[key], rtol=rtol), key


def test_all_rules_agree_for_one_reading():
    readings = _series([[3.0, 0.5, 0]], [[1.0, -0.2, 0.1]], 0.8, 0.5)
    base = posterior_summaries(readings, FLAT, {"rule": "single", "prior": []}, READOUTS)
    variants = [
        (FLAT, SYMMETRIC), (FLAT, EQ18), (FLAT, EQ20),
        ({"reference": "flat", "nuisance": "acceleration"}, {"rule": "likelihood_product", "prior": ["flat_mass"]}),
        ({"reference": "flat", "nuisance": "force"}, {"rule": "likelihood_product", "prior": ["flat_inverse_mass"]}),
    ]
    for law, combine in variants:
        _same(base, posterior_summaries(readings, law, combine, READOUTS))


def test_eq18_is_acceleration_nuisance_with_flat_mass_prior():
    rng = np.random.default_rng(3)
    readings = ReadingSet(rng.normal(size=(3, 4, 3)) * 2 + 3, rng.normal(size=(3, 4, 3)) + 1,
                          rng.uniform(0.5, 2, size=(3, 4)), rng.uniform(0.5, 2, size=(3, 4)))
    a = posterior_summaries(readings, FLAT, EQ18, READOUTS)
    b = posterior_summaries(readings, {"reference": "flat", "nuisance": "acceleration"},
                            {"rule": "likelihood_product", "prior": ["flat_mass"]}, READOUTS)
    _same(a, b)


def test_angle_identities_with_common_instrument_ratio():
    """Contribution 12 §7: eq18 = sym * cos^{2(N-1)}, eq20 = sym * (sin 2θ / 2)^{N-1}."""
    rng = np.random.default_rng(4)
    n, s = 5, 1.7
    readings = ReadingSet(rng.normal(size=(1, n, 3)) + 2, rng.normal(size=(1, n, 3)) + 1,
                          np.full((1, n), s), np.full((1, n), 1.0))
    u = np.linspace(-5, 5, 201)[None]
    theta = np.arctan(np.exp(u) / s)
    sym, _ = joint_log_density(readings, u, FLAT, SYMMETRIC)
    e18, _ = joint_log_density(readings, u, FLAT, EQ18)
    e20, _ = joint_log_density(readings, u, FLAT, EQ20)
    d18 = e18 - sym - 2 * (n - 1) * np.log(np.cos(theta))
    d20 = e20 - sym - (n - 1) * np.log(np.sin(2 * theta) / 2)
    assert np.ptp(d18) < 1e-9 and np.ptp(d20) < 1e-9


def test_contribution_11_combined_readout():
    """Contribution 11 §7: eq. (28) under eq. (18) = 2.008902161134; trial 1 scaled x10 -> 1.997267503687."""
    f = [[2, 0.2, 0], [5, 0, 0.1]]
    a = [[1, 0, 0], [2.1, 0.2, 0]]
    readings = _series(f, a, np.array([[0.5, 1.0]]), np.array([[0.4, 0.6]]))
    out = posterior_summaries(readings, FLAT, EQ18, ["ratio_of_means"])
    assert out["ratio_of_means"][0] == pytest.approx(2.008902161134, rel=1e-10)
    scaled = _series([[20, 2, 0], [5, 0, 0.1]], [[10, 0, 0], [2.1, 0.2, 0]],
                     np.array([[5.0, 1.0]]), np.array([[4.0, 0.6]]))
    out = posterior_summaries(scaled, FLAT, EQ18, ["ratio_of_means"])
    assert out["ratio_of_means"][0] == pytest.approx(1.997267503687, rel=1e-10)


def test_zero_readings():
    """Contribution 11 §4 and 12 §7: exact zero observations."""
    s = 2.0
    for n in (1, 2, 5):
        zeros = _series(np.zeros((n, 3)), np.zeros((n, 3)), s, 1.0)
        sym = posterior_summaries(zeros, FLAT, SYMMETRIC, READOUTS)
        assert sym["log_sd"][0] == pytest.approx(math.pi / 2, rel=1e-12)
        assert sym["ratio_of_means"][0] == pytest.approx(s, rel=1e-12)
    two = _series(np.zeros((2, 3)), np.zeros((2, 3)), s, 1.0)
    assert posterior_summaries(two, FLAT, EQ18, ["ratio_of_means"])["ratio_of_means"][0] == pytest.approx(s / 2, rel=1e-10)
    assert posterior_summaries(two, FLAT, EQ20, ["ratio_of_means"])["ratio_of_means"][0] == pytest.approx(s, rel=1e-10)


def test_sech_tilt_null_log_sds():
    """Contribution 12 §4: the tilt reproduces the closed-form null law's log SDs."""
    zero = _series(np.zeros((1, 3)), np.zeros((1, 3)), 1.0, 1.0)
    expected = {-4: 2.852883, 0: 1.570796, 4: 0.666264, 16: 0.263564, 64: 0.126513}
    for lam, sd in expected.items():
        combine = {"rule": "single", "prior": [{"sech_tilt": {"lambda": lam, "center": "first_reading"}}]}
        out = posterior_summaries(zero, FLAT, combine, ["log_sd", "ratio_of_means"])
        assert out["log_sd"][0] == pytest.approx(sd, abs=1e-6)
        assert out["ratio_of_means"][0] == pytest.approx(1.0, rel=1e-10)


def test_improper_law_is_flagged_not_reported():
    readings = _series([[1.0, 0, 0]], [[1.0, 0, 0]], 1.0, 1.0)
    out = posterior_summaries(readings, FLAT, {"rule": "likelihood_product", "prior": ["flat_mass"]}, READOUTS)
    assert out["unresolved"][0]
    assert np.isnan(out["ratio_of_means"][0]) and np.all(np.isnan(out["interval_95"][0]))


def test_pooling_a_repeated_pair_equals_one_precise_reading():
    rng = np.random.default_rng(5)
    f = rng.normal(size=(1, 20, 3)) + [3, 0, 0]
    a = rng.normal(size=(1, 20, 3)) + [1, 0, 0]
    readings = ReadingSet(f, a, np.full((1, 20), 1.0), np.full((1, 20), 1.0))
    pooled = readings.pooled()
    direct = _series(f.mean(axis=1), a.mean(axis=1), 1 / math.sqrt(20), 1 / math.sqrt(20))
    _same(posterior_summaries(pooled, FLAT, {"rule": "single", "prior": []}, READOUTS),
          posterior_summaries(direct, FLAT, {"rule": "single", "prior": []}, READOUTS), rtol=1e-12)


# ---------------------------------------------------------------- added after the independent review

import zlib

from met.analyst import tail_slopes


@pytest.mark.parametrize("nuisance,rule,prior,n", [
    (nu, rule, prior, n)
    for nu in ("radius", "acceleration", "force")
    for rule in ("single", "likelihood_product", "posterior_mass", "posterior_log")
    for prior in ([], ["flat_mass"], ["flat_log_mass"], ["flat_inverse_mass"],
                  [{"uniform_angle": {"center": "first_reading"}}],
                  [{"sech_tilt": {"lambda": 2.0, "center": 1.0}}])
    for n in (1, 3)
    if not (rule == "single" and n != 1)
])
def test_static_tail_slopes_match_the_computed_density(nuisance, rule, prior, n):
    """The properness check reads tail slopes from a table; compare with the real density."""
    combine = {"rule": {"posterior_mass": "posterior_product", "posterior_log": "posterior_product"}.get(rule, rule),
               "prior": prior}
    if rule == "posterior_mass":
        combine["coordinate"] = "mass"
    if rule == "posterior_log":
        combine["coordinate"] = "log_mass"
    rng = np.random.default_rng(zlib.crc32(f"{nuisance}{rule}{n}{prior}".encode()))
    readings = ReadingSet(rng.normal(size=(1, n, 3)) * 2 + 1, rng.normal(size=(1, n, 3)) + 0.5,
                          rng.uniform(0.5, 2, size=(1, n)), rng.uniform(0.5, 2, size=(1, n)))
    law = {"reference": "flat", "nuisance": nuisance}
    u = np.array([[-80.0, -70.0, 70.0, 80.0]])
    lp, _ = joint_log_density(readings, u, law, combine)
    minus = (lp[0, 1] - lp[0, 0]) / 10
    plus = (lp[0, 3] - lp[0, 2]) / 10
    expected_plus, expected_minus, _ = tail_slopes(law, combine, n)
    assert plus == pytest.approx(expected_plus, abs=1e-6)
    assert minus == pytest.approx(expected_minus, abs=1e-6)


@pytest.mark.parametrize("snr", [3e4, 1e5, 1e6])
def test_very_precise_readings_are_resolved(snr):
    """At high SNR the law of log m is close to normal with sd sqrt(1/snr_F^2 + 1/snr_a^2)."""
    readings = _series([[3.0 * snr, 0, 0]], [[snr, 0, 0]], 1.0, 1.0)
    out = posterior_summaries(readings, FLAT, {"rule": "single", "prior": []},
                              ["ratio_of_means", "median", "log_sd", "interval_95"])
    sd = math.sqrt(1 / (3 * snr) ** 2 + 1 / snr**2)
    assert not out["unresolved"][0]
    assert out["log_sd"][0] == pytest.approx(sd, rel=1e-3)
    assert out["ratio_of_means"][0] == pytest.approx(3.0, abs=0.05 * sd * 3)
    lo, hi = np.log(out["interval_95"][0] / out["median"][0])
    assert -lo == pytest.approx(1.959964 * sd, rel=2e-3)
    assert hi == pytest.approx(1.959964 * sd, rel=2e-3)


def test_pooling_equals_the_same_pair_likelihood_by_brute_force():
    """1D: the N-reading same-pair law, integrated directly over (f, alpha, direction),
    against the estimator applied to the pooled pair."""
    rng = np.random.default_rng(6)
    n, sf, sa = 4, 1.0, 0.7
    f_obs = 2.4 + sf * rng.standard_normal(n)
    a_obs = 1.1 + sa * rng.standard_normal(n)
    f = np.linspace(1e-6, 8.0, 1601)[:, None]
    a = np.linspace(1e-6, 4.0, 1601)[None, :]
    logs = []
    for sign in (1.0, -1.0):
        q = sum((fo - sign * f) ** 2 / sf**2 + (ao - sign * a) ** 2 / sa**2 for fo, ao in zip(f_obs, a_obs))
        logs.append(-q / 2)
    top = max(np.max(x) for x in logs)
    w = sum(np.exp(x - top) for x in logs)
    brute = float(np.sum(w * f) / np.sum(w * a))
    readings = ReadingSet(f_obs[None, :, None], a_obs[None, :, None], np.full((1, n), sf), np.full((1, n), sa))
    pooled = posterior_summaries(readings.pooled(), FLAT, {"rule": "single", "prior": []}, ["ratio_of_means"])
    assert pooled["ratio_of_means"][0] == pytest.approx(brute, rel=2e-5)
