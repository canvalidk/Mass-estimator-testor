"""Point-level composition on data, and the parametrisation-free point (O2, 28 Oct, Claude).

1. Lumping on data: bodies 1 and 2 share each acceleration reading, and the lump's force
   reading is F1 + F2 (its noise adds). L5 asks R(lump) = R(1) + R(2). The cloud readout (T4)
   and the law's median fail it by about the same amount at finite noise; both pass as the
   signal grows (L11).
2. T35: a point that commutes with every increasing reparametrisation of mass is a quantile,
   and the swap makes it the median.
"""
import math
import os

import numpy as np
import pytest

from met.full import isotropic, ours29

SF, SA = 0.2, 0.1
M1, M2 = 1.0, 2.0


def _lumping_defects(size, n_readings, n_series, seed, priors=(1.0, 2.0, 3.0)):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_series):
        a = rng.standard_normal((n_readings, 3))
        a *= size / np.linalg.norm(a, axis=1, keepdims=True)
        at = a + SA * rng.standard_normal((n_readings, 3))
        f1 = M1 * a + SF * rng.standard_normal((n_readings, 3))
        f2 = M2 * a + SF * rng.standard_normal((n_readings, 3))
        one = ours29([isotropic(f, x, SF, SA) for f, x in zip(f1, at)], priors[0])
        two = ours29([isotropic(f, x, SF, SA) for f, x in zip(f2, at)], priors[1])
        lump = ours29([isotropic(f, x, SF * math.sqrt(2), SA) for f, x in zip(f1 + f2, at)], priors[2])
        out.append([abs(lump[k] - one[k] - two[k]) / (one[k] + two[k]) for k in ("readout", "median")])
    return np.array(out)


@pytest.mark.parametrize("size, low, high", [(2.0, 0.0, 0.01), (0.1, 0.05, 1.0)])
def test_cloud_readout_and_median_fail_lumping_on_data_alike(size, low, high):
    """Median over series of |R(lump) - R(1) - R(2)| / (R(1) + R(2)); priors compose (1 + 2 = 3).
    Recorded 28 Oct, 30 series, seed 7: strong 0.273% (cloud) vs 0.273% (median); weak 20.96%
    vs 23.25%. Default run: 8 series; MET_SLOW: 30."""
    n = 30 if os.environ.get("MET_SLOW") else 8
    d = np.median(_lumping_defects(size, 4, n, seed=7), axis=0)
    assert low <= d[0] <= high and low <= d[1] <= high
    assert 0.6 < d[0] / d[1] < 1.4


def test_lumping_defect_shrinks_with_signal():
    """L5 holds in the limit for both (consistency, L11), not at finite noise."""
    weak, mid, strong = (np.median(_lumping_defects(s, 4, 6, seed=11), axis=0) for s in (0.1, 0.4, 2.0))
    assert (weak > mid).all() and (mid > strong).all()


# --- T35: the parametrisation-free point --------------------------------------------------

def _quantile(lm, law, q):
    w = law / law.sum()
    cdf = np.cumsum(w) - 0.5 * w            # midpoint rule: symmetric when the order is reversed
    return float(np.interp(q, cdf, lm))


@pytest.mark.parametrize("transform", [lambda x: -x, lambda x: 2 * x + 0.3, lambda x: x ** 3, np.sinh])
def test_t35_quantiles_commute_with_increasing_maps_of_log_mass(transform):
    """Any increasing g: the q-quantile of g(ln m) is g(q-quantile of ln m); a decreasing g
    (the swap, -x) sends q to 1 - q. Checked on an ours29 law."""
    r = ours29([isotropic([0.4, 0.1, 0.0], [0.2, 0.05, 0.0], SF, SA)], 2.0, readout=False)
    lm, law = r["lm"], r["law"]
    for q in (0.25, 0.5, 0.8):
        y = transform(lm)
        order = np.argsort(y)
        direct = _quantile(y[order], law[order], q)
        increasing = transform(1.0) > transform(0.0)
        through = transform(_quantile(lm, law, q if increasing else 1 - q))
        assert direct == pytest.approx(through, abs=2e-3 * (1 + abs(through)))


def test_t35_geometric_mean_is_not_parametrisation_free():
    """exp E[ln m] commutes with scaling and the swap but not with m -> m + const in ln m's
    cube: it fails the general increasing map, so the swap and units alone do not force the median."""
    r = ours29([isotropic([0.4, 0.1, 0.0], [0.2, 0.05, 0.0], SF, SA)], 2.0, readout=False)
    lm, law = r["lm"], r["law"] / r["law"].sum()
    mean_log = float((lm * law).sum())
    cubed = float(((lm ** 3) * law).sum())
    assert abs(cubed - mean_log ** 3) > 1e-3


# --- the point in met/ (28 Oct) --------------------------------------------------------------

def _null(n):
    return [isotropic([0.0] * 3, [0.0] * 3, SF, SA) for _ in range(n)]


def test_point_is_the_median_and_meets_z_exactly():
    sig = [isotropic([0.4 * v for v in f], [0.4 * v for v in x], SF, SA)
           for f, x in zip([[4.0, 0.3, 0.0], [6.1, -0.2, 0.1]], [[2.0, 0.0, 0.1], [2.0, 0.1, 0.0]])]
    base = ours29(sig, 2.0, readout=False)
    assert base["point"] == base["median"]
    padded = ours29(sig + _null(50), 2.0, readout=False)
    assert padded["point"] == pytest.approx(base["point"], rel=1e-9)


@pytest.mark.parametrize("c", [0.5, 2.0, 7.0])
def test_point_returns_the_centre_on_the_all_zero_experiment(c):
    """For every c, not only c = s (the cloud readout gives R(c), O2)."""
    r = ours29(_null(3), c, readout=False)
    assert r["point"] == pytest.approx(c, rel=1e-6)


def test_point_swap():
    fwd = ours29([isotropic([0.5, 0.2, 0.0], [0.2, 0.0, 0.1], SF, SA)], 1.3, readout=False)
    bwd = ours29([isotropic([0.2, 0.0, 0.1], [0.5, 0.2, 0.0], SA, SF)], 1 / 1.3, readout=False)
    assert fwd["point"] * bwd["point"] == pytest.approx(1, rel=1e-6)
