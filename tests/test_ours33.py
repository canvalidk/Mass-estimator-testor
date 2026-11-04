"""ours33, the probe (met/ours33.py; problems/O18.md, 4 Nov): construction checks and its premise report.

Each test is one row of the report: passes / fails / unknown per premise. ours33 is T48's exponent on
ours32's residual with T47's normal-score warp in place of the power eta. Coverage is
tools/ours33_coverage.py.
"""
import math

import numpy as np
import pytest

from met.full import isotropic
from met.ours31 import _per_reading
from met.ours32 import ours32
from met.ours32 import _l_terms as _l_terms_32
from met.ours33 import _l_terms, _tail_quantiles, ours33, warp_probabilities


def _iso(f, a, sd=1.0):
    return isotropic(np.atleast_1d(f), np.atleast_1d(a), sd, sd)


def _series(seed=1, n=6, k=3, signal=1.0):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        x = signal * np.eye(k)[0]
        out.append(isotropic(x + rng.normal(size=k), x + rng.normal(size=k), 1.0, 1.0))
    return out


def _t39_pair(radius, nu=10.0, d=3.0):
    """T39 part 5: two 1D readings at true mass 2 with opposite perpendicular errors of size d."""
    th0 = math.atan(2.0)
    u = np.array([math.sin(th0), math.cos(th0)])
    v = np.array([math.cos(th0), -math.sin(th0)])
    plus, minus = radius * u + d * v, radius * u - d * v
    return [_iso(plus[0], plus[1]), _iso(minus[0], minus[1])]


# --- construction -------------------------------------------------------------------------------

def test_exponent_is_t48s_one_more_residual_dimension():
    """ours33's log factor is ours32's scaled by (nu + k)/(nu + k - 1): the same residual, exponent (nu+k)/2."""
    rd = _series(seed=2)
    st = _per_reading(rd, 1.0, [10.0] * len(rd))
    phi = np.array([0.3, -1.0, 1.2])
    for p33, p32, (r, ne, _, _) in zip(_l_terms(st, 1.0, phi, 2), _l_terms_32(st, 1.0, phi, 2), st):
        ratio = (ne + r.k) / (ne + r.k - 1)
        for a, b in zip(p33, p32):
            assert np.allclose(a, ratio * b, rtol=1e-12, atol=1e-12)


def test_derivatives_match_finite_differences():
    rd = _series(seed=3)
    st = _per_reading(rd, 1.0, [10.0] * len(rd))
    phi, h = np.array([0.2, -0.9, 1.1]), 1e-5
    for i, p in enumerate(_l_terms(st, 1.0, phi, 2)):
        up, dn = _l_terms(st, 1.0, phi + h, 0)[i][0], _l_terms(st, 1.0, phi - h, 0)[i][0]
        assert np.allclose(p[1], (up - dn) / (2 * h), rtol=1e-5, atol=1e-6)
        assert np.allclose(p[2], (up - 2 * p[0] + dn) / h ** 2, rtol=1e-3, atol=1e-4)


def test_warp_fixes_the_median_and_is_the_identity_at_eta_one():
    assert warp_probabilities(1.0, [0.025, 0.5, 0.975]) == pytest.approx([0.025, 0.5, 0.975])
    for eta in (0.7, 0.2, 0.01):
        lo, med, hi = warp_probabilities(eta, [0.025, 0.5, 0.975])
        assert med == 0.5 and lo < 0.025 and hi > 0.975 and lo == pytest.approx(1 - hi)


def test_tail_quantiles_resolve_ranks_below_floating_point():
    lm = np.linspace(-14.0, 14.0, 28001)
    lp = -0.5 * lm ** 2                                   # a standard normal law in ln m
    lo, hi = _tail_quantiles(lm, lp, 0.025)
    assert math.log(lo) == pytest.approx(-1.95996, abs=2e-3) and math.log(hi) == pytest.approx(1.95996, abs=2e-3)
    lo, hi = _tail_quantiles(lm, lp, 1e-30)              # 1 - 1e-30 == 1.0 in float64: the upper tail is
    assert math.log(hi) == pytest.approx(11.4638, abs=2e-3)   # accumulated from its own end
    assert math.log(lo) == pytest.approx(-11.4638, abs=2e-3)


def test_point_is_the_raw_median_whatever_eta():
    rd = _series(seed=4)
    a, b = ours33(rd, 1.0, 10.0), ours33(rd, 1.0, 10.0, calibrate=False)
    assert a["eta"] < 1.0 and a["median"] == b["median"] and a["interval"] != a["raw_interval"]
    assert b["interval"] == b["raw_interval"]
    assert a["interval"][0] < a["raw_interval"][0] and a["interval"][1] > a["raw_interval"][1]


# --- the two T39 breaks of ours32, by construction -----------------------------------------------

@pytest.mark.parametrize("radius", [1000.0, 10000.0])
def test_t39_nu_below_one_single_reading_concentrates_under_t48(radius):
    """T39 part 4: one 1D reading at nu = 1/2, true mass 2, amplitude growing: ours32's median tends to
    1.82; T48's exponent (nu + k)/2 = 3/4 gives an integrable tail and the raw law concentrates."""
    th0 = math.atan(2.0)
    y = radius * np.array([math.sin(th0), math.cos(th0)]) + np.array([0.3, -0.2])
    rd = [_iso(y[0], y[1])]
    assert ours32(rd, 1.0, 0.5)["median"] < 1.85
    assert ours33(rd, 1.0, 0.5)["median"] == pytest.approx(2.0, rel=0.02)


def test_t39_opposed_pair_at_nu_ten_concentrates_under_the_warp():
    """T39 part 5: eta -> 1/90 makes ours32's powered law non-concentrating (median -> 1.31). The warp
    leaves the raw law's median and concentration alone; eta only widens the quoted interval."""
    rd = _t39_pair(1000.0)
    a, b = ours32(rd, 1.0, 10.0), ours33(rd, 1.0, 10.0)
    assert a["median"] == pytest.approx(1.3022, abs=0.01)
    assert b["median"] == pytest.approx(2.0, abs=1e-3)
    assert b["raw_interval"][0] > 1.95 and b["raw_interval"][1] < 2.05
    assert b["eta"] < 0.02 and b["interval"][0] <= 2.0 <= b["interval"][1]


def test_the_warped_width_at_small_eta_is_the_starting_laws_tail():
    """Hole 1 of O18 (4 Nov), measured: on T39's pair the warp asks for the raw law's 1e-80 ranks, which
    lie in pi's Cauchy tails; the quoted interval spans tens of decades while both readings sit at 2."""
    b = ours33(_t39_pair(1000.0), 1.0, 10.0)
    lo, hi = b["interval"]
    assert math.log10(hi / lo) > 60


# --- premise rows shared with ours31/ours32 ---------------------------------------------------------

@pytest.mark.parametrize("nu_eff", [0.5, 4.0, 50.0])
def test_l2_totality_on_b11s_inputs(nu_eff):
    cases = [[_iso(2.0, 1.0)], [_iso(0.0, 1.0)], [_iso(1.0, 0.0)], [_iso(-2.0, 1.0)],
             [_iso([2.0, 0, 0], [-1.0, 0, 0])], [_iso(0.0, 0.0)] * 3, [_iso(2.0, 1.0)] * 5]
    for rd in cases:
        r = ours33(rd, 1.0, nu_eff)
        lo, hi = r["interval"]
        assert 0 < lo < r["median"] < hi < math.inf and 0 < r["eta"] <= 1


def test_z_and_all_zero():
    rd = [_iso([2.0, 0.1, 0], [1.0, 0, 0.1]), _iso([1.0, 0.3, 0], [0.6, -0.2, 0])]
    z = _iso([0.0, 0, 0], [0.0, 0, 0])
    base, padded = ours33(rd, 1.0, 10.0), ours33(rd + [z] * 4, 1.0, 10.0)
    assert padded["eta"] == base["eta"] and padded["median"] == base["median"]
    assert np.allclose(base["law"], padded["law"], rtol=0, atol=1e-15)
    for c in (0.3, 1.0, 7.0):
        r = ours33([isotropic([0.0, 0, 0], [0.0, 0, 0], 1.0, 1.0 / c)] * 3, c, 10.0)
        assert r["eta"] == 1.0 and r["median"] == pytest.approx(c, rel=1e-6)


def test_units_and_swap():
    rd = _series(seed=5)
    base = ours33(rd, 1.0, 10.0)
    lam = 3.0
    scaled = ours33([isotropic(r.force * lam, r.acceleration, lam, 1.0) for r in rd], lam, 10.0)
    assert scaled["median"] == pytest.approx(lam * base["median"], rel=1e-4)
    assert scaled["eta"] == pytest.approx(base["eta"], rel=1e-9)
    swapped = ours33([isotropic(r.acceleration, r.force, 1.0, 1.0) for r in rd], 1.0, 10.0)
    assert swapped["median"] == pytest.approx(1 / base["median"], rel=1e-4)
    assert swapped["interval"][0] == pytest.approx(1 / base["interval"][1], rel=1e-3)


@pytest.mark.parametrize("pair", ["misaligned 30 deg, 3D", "anti-parallel, 1D", "anti-parallel, 3D"])
def test_t30_bounded_as_the_declared_noise_shrinks(pair):
    """T30 for the raw law (T48 section 3): a fixed pair no positive mass fits keeps bounded,
    non-vanishing relative evidence as sigma -> 0; the raw interval converges to a positive width."""
    th = math.radians(30)
    f, a = {"misaligned 30 deg, 3D": ([2 * math.cos(th), 2 * math.sin(th), 0], [1.0, 0, 0]),
            "anti-parallel, 1D": ([-2.0], [1.0]),
            "anti-parallel, 3D": ([-2.0, 0, 0], [1.0, 0, 0])}[pair]
    widths = []
    for sd in (1e-1, 1e-2, 1e-3):
        lo, hi = ours33([_iso(f, a, sd)], 1.0, 10.0, calibrate=False)["raw_interval"]
        widths.append(math.log(hi / lo))
    assert widths[-1] > 0.05 and abs(widths[-1] - widths[-2]) < 0.05 * widths[-1]
