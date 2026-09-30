"""ratio_of_means_rss (ours25b): the many-reading readout pooled as the stacked vector's length.

Each check is against a computation done a different way: the noncentral chi
mean by quadrature over the component along w and the chi-square of the rest;
the whole readout by an independent quadrature in the angle theta (not the
log-mass grid). Then the properties the readout was built for: a 3D reading and
its three axes agree, the recording frame does not matter, and the series enters
only through three sums.
"""

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import roots_legendre
from scipy.stats import chi2
from scipy.spatial.transform import Rotation

from met.analyst import Estimator, estimate
from met.law import ReadingSet, cartesian2_mean_radius, stacked_mean_radius, _noncentral_chi_mean_series

CART2 = {"reference": "cartesian2", "nuisance": "radius"}
CART1 = {"reference": "cartesian1", "nuisance": "radius"}
PRODUCT = {"rule": "likelihood_product", "prior": [{"uniform_angle": {"center": "first_reading"}}]}
BOTH = ["ratio_of_means", "ratio_of_means_rss", "median"]


def _rss(force, acceleration, sf=1.0, sa=1.0, law=CART2, readouts=BOTH):
    return estimate(force, acceleration, sf, sa, readouts, law=law, combine=PRODUCT)


def _split(force, acceleration):
    return np.asarray(force, float).reshape(-1, 1), np.asarray(acceleration, float).reshape(-1, 1)


# ---------------------------------------------------------------- the noncentral chi mean

@pytest.mark.parametrize("k", [1, 2, 3])
def test_series_matches_closed_forms(k):
    h = np.array([0.0, 1e-3, 0.4, 1.7, 5.0, 23.0, 140.0])
    assert _noncentral_chi_mean_series(k, h) == pytest.approx(cartesian2_mean_radius(h, k), rel=1e-10)


@pytest.mark.parametrize("k", [4, 7, 30, 150])
@pytest.mark.parametrize("h", [0.0, 0.8, 4.0, 25.0])
def test_stacked_mean_radius_by_quadrature(k, h):
    """E sqrt((h + z)^2 + R^2), z ~ N(0, 1) along w, R^2 ~ chi^2_{k-1} across it."""
    def inner(z):
        return quad(lambda q: math.sqrt((h + z) ** 2 + q) * chi2.pdf(q, k - 1), 0, np.inf, limit=200)[0]
    value = quad(lambda z: inner(z) * math.exp(-0.5 * z * z) / math.sqrt(2 * math.pi), -12, 12, limit=200)[0]
    assert stacked_mean_radius(np.array([h]), k)[0] == pytest.approx(value, rel=1e-8)


def test_stacked_table_extends_and_stays_exact():
    h = np.array([3.0, 900.0, 2500.0])
    assert stacked_mean_radius(h, 12) == pytest.approx(_noncentral_chi_mean_series(12, h), rel=1e-8)


# ---------------------------------------------------------------- the readout

def _theta_quadrature(force, acceleration, sf, sa, points=4000):
    """The ours24 law and both poolings, computed in theta with its own code
    (Gauss-Legendre on (0, pi/2))."""
    x, y = np.asarray(force, float) / sf, np.asarray(acceleration, float) / sa
    n, d = x.shape
    nodes, weights = roots_legendre(points)
    t = (nodes + 1) * math.pi / 4
    sn, cs = np.sin(t), np.cos(t)
    lp = -0.5 * np.sum((x[:, :, None] * cs - y[:, :, None] * sn) ** 2, axis=(0, 1))
    w = weights * np.exp(lp - lp.max())
    h2 = np.sum((x[:, :, None] * sn + y[:, :, None] * cs) ** 2, axis=1)          # (N, T)
    g_sum = np.sum(cartesian2_mean_radius(np.sqrt(h2), d), axis=0)
    g_rss = stacked_mean_radius(np.sqrt(np.sum(h2, axis=0)), n * d)
    s = sf / sa
    return {name: s * np.sum(w * sn * g) / np.sum(w * cs * g)
            for name, g in (("ratio_of_means", g_sum), ("ratio_of_means_rss", g_rss))}


@pytest.mark.parametrize("n,d", [(1, 3), (3, 1), (4, 2), (5, 3)])
def test_readout_against_theta_quadrature(n, d):
    rng = np.random.default_rng(100 * n + d)
    f, a = rng.normal(1.0, 1.0, (n, d)), rng.normal(0.4, 0.5, (n, d))
    got = _rss(f, a, 0.7, 0.4)
    ref = _theta_quadrature(f, a, 0.7, 0.4)
    assert got["ratio_of_means"] == pytest.approx(ref["ratio_of_means"], rel=1e-9)
    assert got["ratio_of_means_rss"] == pytest.approx(ref["ratio_of_means_rss"], rel=1e-8)


@pytest.mark.parametrize("d", [1, 2, 3])
def test_one_reading_is_the_ratio_of_means(d):
    rng = np.random.default_rng(d)
    for _ in range(3):
        r = _rss(rng.normal(size=d) * 2, rng.normal(size=d), 1.3, 0.6)
        assert r["ratio_of_means_rss"] == pytest.approx(r["ratio_of_means"], rel=1e-12)


@pytest.mark.parametrize("n", [1, 3, 20])
def test_zero_readings_give_the_instrument_ratio(n):
    r = _rss(np.zeros((n, 3)), np.zeros((n, 3)), 2.0, 0.5)
    assert r["ratio_of_means_rss"] == pytest.approx(4.0, rel=1e-10)


def test_reciprocal_symmetry():
    rng = np.random.default_rng(3)
    f, a = rng.normal(size=(4, 3)) * 2, rng.normal(size=(4, 3))
    r = _rss(f, a, 1.5, 0.5)
    w = _rss(a, f, 0.5, 1.5)
    assert r["ratio_of_means_rss"] * w["ratio_of_means_rss"] == pytest.approx(1, rel=1e-9)


def test_cartesian1_and_cartesian2_agree():
    rng = np.random.default_rng(4)
    f, a = rng.normal(size=(3, 3)), rng.normal(size=(3, 3))
    assert _rss(f, a, law=CART1)["ratio_of_means_rss"] == pytest.approx(
        _rss(f, a, law=CART2)["ratio_of_means_rss"], rel=1e-10)


# ---------------------------------------------------------------- what it was built for

@pytest.mark.parametrize("n", [1, 4])
def test_axis_splitting_holds_for_the_readout(n):
    """A 3D series and the same series as 3N one-dimensional readings: ours25b agrees
    exactly; ours24 (plain sum of sizes) does not."""
    rng = np.random.default_rng(20 + n)
    f, a = rng.normal(0.8, 1.0, (n, 3)), rng.normal(0.3, 1.0, (n, 3))
    vector, axes = _rss(f, a, 0.9, 0.6), _rss(*_split(f, a), 0.9, 0.6)
    assert axes["ratio_of_means_rss"] == pytest.approx(vector["ratio_of_means_rss"], rel=1e-9)
    assert axes["median"] == pytest.approx(vector["median"], rel=1e-9)
    assert abs(math.log(axes["ratio_of_means"] / vector["ratio_of_means"])) > 1e-4


def test_split_axes_reduce_matches_splitting_by_hand():
    rng = np.random.default_rng(5)
    f, a = rng.normal(size=(2, 4, 3)), rng.normal(size=(2, 4, 3))
    readings = ReadingSet(f, a, np.full((2, 4), 0.8), np.full((2, 4), 0.5))
    spec = {"name": "e", "reduce": "split_axes", "law": CART2, "combine": PRODUCT, "readouts": BOTH}
    out = Estimator(spec).evaluate(readings)
    for b in range(2):
        by_hand = _rss(*_split(f[b], a[b]), 0.8, 0.5)
        for name in BOTH:
            assert out[name][b] == pytest.approx(by_hand[name], rel=1e-12)


def test_recording_frame_does_not_matter():
    """Recorded as three 1D readings, ours25b does not depend on the orientation of
    the lab axes; ours24 does."""
    rng = np.random.default_rng(6)
    f, a = rng.normal(1.0, 1.0, (1, 3)), rng.normal(0.5, 1.0, (1, 3))
    rss, plain = [], []
    for rot in Rotation.random(6, random_state=7).as_matrix():
        r = _rss(*_split(f @ rot.T, a @ rot.T))
        rss.append(r["ratio_of_means_rss"])
        plain.append(r["ratio_of_means"])
    assert np.ptp(np.log(rss)) < 1e-9
    assert np.ptp(np.log(plain)) > 1e-3


def test_series_enters_only_through_three_sums():
    """Two different 2 x 2D series with the same |X|^2, |Y|^2 and X.Y: one ours25b answer."""
    x1 = np.array([[2.0, 0.0], [0.0, 0.5]])
    y1 = np.array([[1.5, 0.3], [0.2, 0.4]])
    stacked = np.stack([x1.ravel(), y1.ravel()], axis=1)
    q, _ = np.linalg.qr(np.random.default_rng(8).normal(size=(4, 4)))
    x2, y2 = (q @ stacked)[:, 0].reshape(2, 2), (q @ stacked)[:, 1].reshape(2, 2)
    r1, r2 = _rss(x1, y1), _rss(x2, y2)
    assert r2["ratio_of_means_rss"] == pytest.approx(r1["ratio_of_means_rss"], rel=1e-9)
    assert r2["median"] == pytest.approx(r1["median"], rel=1e-9)
    assert abs(math.log(r2["ratio_of_means"] / r1["ratio_of_means"])) > 1e-4


# ---------------------------------------------------------------- refusals

def _spec(**changes):
    spec = {"name": "e", "reduce": "none", "law": CART2, "combine": PRODUCT, "readouts": ["ratio_of_means_rss"]}
    spec.update(changes)
    return spec


def test_refuses_the_flat_reference():
    with pytest.raises(ValueError, match="cartesian"):
        Estimator(_spec(law={"reference": "flat", "nuisance": "radius"}))


def test_refuses_sequential_rule_accepts_hierarchical():
    """The hierarchical rule's stacked readout was added on 27 September
    (tests/test_hierarchical_ours26.py); the sequential rule without stacked_state is refused."""
    Estimator(_spec(law=CART1, combine={"rule": "hierarchical", "prior": [{"uniform_angle": {
        "center": "first_reading"}}], "excitation": {"beta_b": 1}}))
    with pytest.raises(ValueError, match="ratio_of_means_rss needs rule"):
        Estimator(_spec(combine={"rule": "sequential", "carry": "exact", "old": PRODUCT}))


def test_unequal_sds_within_a_series_are_ours26():
    """ours25b refused SDs that differ between readings of a series; ours26 computes them
    (the weighted stacked mean, tested in test_ours26.py). One reading's SDs agree, so the
    closed form still holds for it alone."""
    readings = ReadingSet(np.ones((1, 2, 3)), np.ones((1, 2, 3)), np.array([[1.0, 2.0]]), np.array([[1.0, 1.0]]))
    out = Estimator(_spec(combine={"rule": "likelihood_product", "prior": [{"uniform_angle": {
        "center": 1.5}}]})).evaluate(readings)
    assert np.isfinite(out["ratio_of_means_rss"][0])
