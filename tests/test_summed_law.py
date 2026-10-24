"""The likelihood product from three sums (ReadingSet.summed_log_like) and the lazy E[alpha | m].

Two speed changes of 24 Oct, checked against the per-reading computation they replace:

- reading_terms(alpha=False) skips E[alpha | m] (erf, Bessel); log_like and log_ref must be
  the same bit for bit, under every reference and nuisance.
- Under either cartesian reference a reading's radius-nuisance log-likelihood is h_i^2 / 2,
  linear in (P_i, Q_i, D_i) at a common angle, so with one s per series the (weighted) sum
  over readings is one curve from three sums. It must equal the per-reading sum to rounding,
  and every readout built on it must agree with the per-reading route to rounding.
- The sandwich (ours25) keeps the per-reading sum: its weight on a series whose best line is
  at negative mass is set by rounding noise (sandwich_weight's docstring), so a different
  summation order would move its recorded numbers.
"""

import numpy as np
import pytest

from met.analyst import joint_log_density, posterior_summaries, sandwich_weight
from met.law import ReadingSet

PRIOR = [{"uniform_angle": {"center": "first_reading"}}]


def _readings(rng, b, n, d, signal, same_s=True, per_axis=False):
    m = np.exp(rng.normal(size=(b, 1, 1)))
    alpha = signal * rng.normal(size=(b, n, d))
    # SDs differ between readings; with same_s their ratio is bitwise common (powers of 2
    # scale exactly), which is what summed_log_like requires
    scale = 2.0 ** rng.integers(-2, 3, size=(b, n))
    sa = 0.7 * scale
    sf = 0.91 * scale if same_s else rng.uniform(0.2, 2.0, size=(b, n))
    f = m * alpha + sf[..., None] * rng.normal(size=(b, n, d))
    a = alpha + sa[..., None] * rng.normal(size=(b, n, d))
    if per_axis:
        sf = sf[..., None] * rng.uniform(0.5, 2.0, size=(b, n, d))
        sa = np.broadcast_to(sa[..., None], (b, n, d)).copy()
    return ReadingSet(f, a, sf, sa)


@pytest.mark.parametrize("reference", ["flat", "cartesian1", "cartesian2"])
@pytest.mark.parametrize("nuisance", ["radius", "acceleration", "force"])
@pytest.mark.parametrize("d", [1, 2, 3])
def test_skipping_alpha_leaves_the_law_bit_for_bit(reference, nuisance, d):
    r = _readings(np.random.default_rng(d), 3, 5, d, 1.0, same_s=False)
    u = np.linspace(-4, 4, 301)[None].repeat(3, axis=0)
    full = r.reading_terms(u, nuisance, reference=reference)
    lazy = r.reading_terms(u, nuisance, reference=reference, alpha=False)
    assert np.array_equal(full[0], lazy[0]) and np.array_equal(full[1], lazy[1])
    assert full[2] is not None and lazy[2] is None


@pytest.mark.parametrize("reference", ["cartesian1", "cartesian2"])
@pytest.mark.parametrize("d", [1, 2, 3])
@pytest.mark.parametrize("signal", [0.0, 0.5, 4.0])
def test_three_sums_equal_the_per_reading_sum(reference, d, signal):
    rng = np.random.default_rng(int(10 * signal) + d)
    r = _readings(rng, 4, 40, d, signal)
    u = np.log(r.s[:, :1]) + np.linspace(-30, 30, 1201)[None]
    per = r.reading_terms(u, "radius", reference=reference)[0]
    g = rng.uniform(0.1, 1.0, size=r.P.shape)
    for w in (None, g):
        want = np.sum(per if w is None else w[:, :, None] * per, axis=1)
        got = r.summed_log_like(u, reference, w)
        scale = np.sum(r.P + r.Q, axis=1)[:, None]
        assert np.all(np.abs(got - want) <= 1e-13 * scale + 1e-300)


def test_three_sums_refuse_where_the_sum_does_not_reduce():
    rng = np.random.default_rng(3)
    u = np.linspace(-3, 3, 11)[None].repeat(2, axis=0)
    same = _readings(rng, 2, 4, 3, 1.0)
    assert same.summed_log_like(u, "flat") is None
    assert _readings(rng, 2, 4, 3, 1.0, same_s=False).summed_log_like(u, "cartesian2") is None
    assert _readings(rng, 2, 4, 3, 1.0, per_axis=True).summed_log_like(u, "cartesian2") is None
    assert same.summed_log_like(u, "cartesian2") is not None


COMBINES = [
    {"rule": "likelihood_product", "prior": PRIOR},
    {"rule": "likelihood_product", "prior": PRIOR, "misfit": "kink"},
    {"rule": "likelihood_product", "prior": PRIOR, "misfit": "kink", "calibrate": "bartlett"},
    {"rule": "likelihood_product", "prior": PRIOR, "calibrate": "sandwich"},
]


@pytest.mark.parametrize("combine", COMBINES, ids=["ours27", "ours29_law", "ours29", "ours25"])
@pytest.mark.parametrize("d", [1, 3])
def test_readouts_agree_with_the_per_reading_route(combine, d, monkeypatch):
    law = {"reference": "cartesian2", "nuisance": "radius"}
    readouts = ["ratio_of_means_rss", "median", "interval_95", "geometric"]
    r = _readings(np.random.default_rng(20 + d), 6, 30, d, 0.8)
    fast = posterior_summaries(r, law, combine, readouts)
    r2 = _readings(np.random.default_rng(20 + d), 6, 30, d, 0.8)
    monkeypatch.setattr(ReadingSet, "summed_log_like", lambda self, *a, **k: None)
    slow = posterior_summaries(r2, law, combine, readouts)
    for name in readouts:
        np.testing.assert_allclose(fast[name], slow[name], rtol=1e-11, atol=0)


def test_alpha_flag_changes_only_whether_alpha_comes_back():
    law = {"reference": "cartesian2", "nuisance": "radius"}
    r = _readings(np.random.default_rng(9), 3, 12, 3, 1.0)
    u = np.log(r.s[:, :1]) + np.linspace(-5, 5, 401)[None]
    for combine in COMBINES:
        lp1, acc = joint_log_density(r, u, law, combine)
        lp0, none = joint_log_density(r, u, law, combine, alpha=False)
        assert acc is not None and none is None
        np.testing.assert_allclose(lp0, lp1, rtol=1e-12, atol=1e-12 * np.max(np.abs(lp1)))


def test_the_sandwich_search_stays_reading_by_reading(monkeypatch):
    """ours25's weight must not go through the three sums (see the module docstring)."""
    law = {"reference": "cartesian2", "nuisance": "radius"}
    combine = {"rule": "likelihood_product", "prior": PRIOR, "calibrate": "sandwich"}
    def refuse(self, *a, **k):
        raise AssertionError("sandwich_weight used summed_log_like")
    monkeypatch.setattr(ReadingSet, "summed_log_like", refuse)
    w = sandwich_weight(_readings(np.random.default_rng(4), 3, 20, 3, 1.0), law, combine)
    assert w.shape == (3,)
