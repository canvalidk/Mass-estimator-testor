"""ours32, the probe: ours31's reading level with eta put back (met/ours32.py; problems/O18.md, 1 Nov).

Construction checks and the premise rows ours32 shares with ours31. Coverage is tools/ours32_coverage.py.
"""
import math

import numpy as np
import pytest

from met.full import isotropic
from met.ours31 import _per_reading, ours31
from met.ours32 import _l_terms, ours32


def _iso(f, a, sd=1.0):
    return isotropic(np.atleast_1d(f), np.atleast_1d(a), sd, sd)


def _series(seed=1, n=6, k=3, signal=1.0):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        x = signal * np.eye(k)[0]
        out.append(isotropic(x + rng.normal(size=k), x + rng.normal(size=k), 1.0, 1.0))
    return out


def test_without_eta_is_ours31():
    rd = _series()
    a, b = ours31(rd, 1.0, 10.0), ours32(rd, 1.0, 10.0, calibrate=False)
    assert b["median"] == pytest.approx(a["median"], rel=1e-12)
    assert b["interval"] == pytest.approx(a["interval"], rel=1e-12)


def test_l_derivatives_match_finite_differences():
    rd = _series(seed=2)
    st = _per_reading(rd, 1.0, [10.0] * len(rd))
    phi, h = np.array([0.3, -1.0, 1.2]), 1e-5
    for i, p in enumerate(_l_terms(st, 1.0, phi, 2)):
        up, dn = _l_terms(st, 1.0, phi + h, 0)[i][0], _l_terms(st, 1.0, phi - h, 0)[i][0]
        assert np.allclose((up - dn) / (2 * h), p[1], atol=1e-8)
        assert np.allclose((up - 2 * p[0] + dn) / h ** 2, p[2], atol=1e-3)


def test_l_is_finite_on_the_whole_line():
    """H <= E gives H+ - H >= -Q+, so the log's argument is at least nu/(nu + Q+) > 0."""
    rd = [_iso([2.0, 0, 0], [-1.0, 0, 0]), _iso(-2.0, 1.0)] + _series(seed=3)
    st = _per_reading(rd, 1.0, [0.5] * len(rd))
    phi = np.linspace(-math.pi / 2 + 1e-9, math.pi / 2, 2001)
    for p in _l_terms(st, 1.0, phi, 2):
        assert all(np.all(np.isfinite(t)) for t in p)


def test_z_and_all_zero():
    rd, z = _series(seed=4), _iso([0.0, 0, 0], [0.0, 0, 0])
    a, b = ours32(rd, 1.0, 10.0), ours32(rd + [z, z], 1.0, 10.0)
    assert b["eta"] == pytest.approx(a["eta"], rel=1e-9)
    assert b["median"] == pytest.approx(a["median"], rel=1e-9)
    for c in (0.3, 1.0, 7.0):
        assert ours32([z] * 3, c, 10.0)["median"] == pytest.approx(c, rel=1e-6)


@pytest.mark.parametrize("nu_eff", [0.5, 4.0, 50.0])
def test_l2_totality_on_b11s_inputs(nu_eff):
    cases = [[_iso(2.0, 1.0)], [_iso(0.0, 1.0)], [_iso(1.0, 0.0)], [_iso(-2.0, 1.0)],
             [_iso([2.0, 0, 0], [-1.0, 0, 0])], [_iso(0.0, 0.0)] * 3, [_iso(2.0, 1.0)] * 5]
    for rd in cases:
        r = ours32(rd, 1.0, nu_eff)
        lo, hi = r["interval"]
        assert 0 < lo < r["median"] < hi < math.inf and 0 < r["eta"] <= 1


def test_units_and_swap():
    rd = _series(seed=5)
    base = ours32(rd, 1.0, 10.0)["median"]
    lam = 3.0
    scaled = [isotropic(r.force * lam, r.acceleration, lam, 1.0) for r in rd]
    assert ours32(scaled, lam, 10.0)["median"] == pytest.approx(lam * base, rel=1e-4)
    swapped = [isotropic(r.acceleration, r.force, 1.0, 1.0) for r in rd]
    assert ours32(swapped, 1.0, 10.0)["median"] == pytest.approx(1 / base, rel=1e-4)
