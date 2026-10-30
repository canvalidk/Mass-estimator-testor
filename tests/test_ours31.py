"""ours31, the probe built on T37 (met/ours31.py; problems/O18.md): its premise report.

Each test is one row of the probe's report: passes / fails / unknown per premise. A failing
premise is recorded as a strict xfail naming what fails, not hidden.
"""
import math

import numpy as np
import pytest

from met import admissible
from met.full import Reading, isotropic, ours30
from met.ours31 import _log_factors, _per_reading, ours31


def _iso(f, a, sd=1.0):
    return isotropic(np.atleast_1d(f), np.atleast_1d(a), sd, sd)


# --- L2: never undefined on measured inputs, including B11's counterexamples ----------------

@pytest.mark.parametrize("nu_eff", [0.5, 4.0, 50.0])
def test_l2_totality_on_b11s_inputs(nu_eff):
    """The exact 1D fit (2, 1), the boundary fits (0, 1) and (1, 0), an anti-parallel pair and the
    all-zero series: every factor is in (0, 1], so the law is at most pi and proper."""
    cases = [[_iso(2.0, 1.0)], [_iso(0.0, 1.0)], [_iso(1.0, 0.0)], [_iso(-2.0, 1.0)],
             [_iso([2.0, 0, 0], [-1.0, 0, 0])], [_iso(0.0, 0.0)] * 3, [_iso(2.0, 1.0)] * 5]
    lm = np.linspace(math.log(1e-8), math.log(1e8), 4001)
    for rd in cases:
        lf = _log_factors(_per_reading(rd, 1.0, [nu_eff] * len(rd)), 1.0, lm)
        assert np.all(np.isfinite(lf)) and lf.max() <= 1e-12
        r = ours31(rd, 1.0, nu_eff)
        lo, hi = r["interval"]
        assert 0 < lo < r["median"] < hi < math.inf


# --- Z (for the law), and the all-zero experiment --------------------------------------------

def test_z_zero_readings_leave_the_law_alone():
    rd = [_iso([2.0, 0.1, 0], [1.0, 0, 0.1]), _iso([1.0, 0.3, 0], [0.6, -0.2, 0])]
    base = ours31(rd, 1.0, 10.0)
    padded = ours31(rd + [_iso([0.0, 0, 0], [0.0, 0, 0])] * 5, 1.0, 10.0)
    assert np.array_equal(base["lm"], padded["lm"]) and np.allclose(base["law"], padded["law"], rtol=0, atol=1e-15)
    assert padded["median"] == base["median"]


@pytest.mark.parametrize("c", [0.3, 1.0, 7.0])
def test_all_zero_returns_the_centre(c):
    r = ours31([isotropic([0.0, 0, 0], [0.0, 0, 0], 1.0, 1.0 / c)] * 4, c, 10.0)
    assert r["median"] == pytest.approx(c, rel=1e-6)


# --- T30: bounded, non-vanishing evidence on a pair no positive mass fits as sigma -> 0 ------

@pytest.mark.parametrize("pair", ["misaligned 30 deg, 3D", "anti-parallel, 1D", "anti-parallel, 3D"])
def test_t30_bounded_as_the_declared_noise_shrinks(pair):
    th = math.radians(30)
    f, a = {"misaligned 30 deg, 3D": ([2 * math.cos(th), 2 * math.sin(th), 0], [1.0, 0, 0]),
            "anti-parallel, 1D": ([-2.0], [1.0]),
            "anti-parallel, 3D": ([-2.0, 0, 0], [1.0, 0, 0])}[pair]
    widths = []
    for sd in (1e-3, 1e-5, 1e-7):
        lo, hi = ours31([_iso(f, a, sd)], 1.0, 10.0)["interval"]
        widths.append(math.log(hi / lo))
    assert widths[2] > 0.05                                  # not collapsing
    assert widths[2] == pytest.approx(widths[1], rel=1e-3)   # and converged: bounded


def test_t30_fails_for_the_fixed_size_gaussian():
    """nu_eff = infinity is the fixed-size Gaussian: confidence on a misaligned pair is unbounded."""
    th = math.radians(30)
    f, a = [2 * math.cos(th), 2 * math.sin(th), 0], [1.0, 0, 0]
    w = [math.log(np.divide(*ours31([_iso(f, a, sd)], 1.0, math.inf)["interval"][::-1])) for sd in (1e-2, 1e-4)]
    assert w[1] < 0.02 * w[0]


# --- L3 (units), the swap (L4 / T5) ----------------------------------------------------------

def test_units_and_swap():
    rng = np.random.default_rng(31)
    rd = [Reading(rng.normal(size=3) + [2, 0, 0], rng.normal(size=3) + [1, 0, 0], np.diag([0.5, 0.5, 0.5, 0.2, 0.2, 0.2]))
          for _ in range(6)]
    base = ours31(rd, 1.6, 8.0)["median"]
    lam, mu = 3.0, 0.5                                       # force x3, acceleration x0.5: m x6
    sc = np.diag([lam] * 3 + [mu] * 3)
    scaled = [Reading(r.force * lam, r.acceleration * mu, sc @ r.cov @ sc) for r in rd]
    assert ours31(scaled, 1.6 * lam / mu, 8.0)["median"] == pytest.approx(base * lam / mu, rel=1e-6)
    p = np.block([[np.zeros((3, 3)), np.eye(3)], [np.eye(3), np.zeros((3, 3))]])
    swapped = [Reading(r.acceleration, r.force, p @ r.cov @ p) for r in rd]
    assert ours31(swapped, 1 / 1.6, 8.0)["median"] == pytest.approx(1 / base, rel=1e-5)


# --- L24 (candidate): the reading's precision is a cut posterior of a proper noise-size law ----

@pytest.mark.parametrize("k", [1, 2, 3])
def test_l24_the_peak_precision_is_b12s_inverse_gamma_control(k):
    """At its own best line the factor's curvature is (nu_eff + nu)/(nu_eff + Q+): B12's inverse-gamma
    control in the precision convention, which B12 classifies as admissible."""
    nu_eff, nu = 10.0, k - 1
    verdict = admissible.check_form("invgamma", nu, alpha=nu_eff / 2, normalization="precision")
    assert verdict["cm"] == "pass" and verdict["proper"] == "pass" and verdict["precision"] == "pass"
    for q in (0.0, 1.0, 7.5, 40.0):
        assert admissible.weight("invgamma", q, nu, alpha=nu_eff / 2, normalization="precision") == \
            pytest.approx((nu_eff + nu) / (nu_eff + q), rel=1e-12)


# --- consistency with the old law ------------------------------------------------------------

def test_infinite_nu_eff_is_ours27s_raw_law():
    rng = np.random.default_rng(27)
    rd = [_iso(rng.normal(size=3) + [2, 0, 0], rng.normal(size=3) + [1, 0, 0]) for _ in range(8)]
    a = ours31(rd, 1.0, math.inf)["median"]
    b = ours30(rd, 1.0, misfit=False, calibrate=False, readout=False)["median"]
    assert a == pytest.approx(b, rel=1e-9)
