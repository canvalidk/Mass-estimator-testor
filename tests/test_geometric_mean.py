"""The comparison readout exp E[ln m] beside the median, on ours30, ours31, ours32 and ours33 (O2, 6 Nov; ours33 7 Nov).

Reported, not the point: G6 is a sort between the reflection-equivariant location functionals of
ln m, and these are the L1 (median) and L2 (geometric mean) members. The tests check that the
geometric mean is a law functional with the same equivariances as the median, that the two agree
where the law concentrates, and record how far apart they are where it does not.
"""
import math

import numpy as np
import pytest

from met.full import Reading, geometric_mean, isotropic, ours30
from met.ours31 import ours31
from met.ours32 import ours32
from met.ours33 import ours33

C = 2.0
NU = 10.0


def _series(rng, m, n, signal, sf=0.2, sa=0.1, d=3):
    out = []
    for _ in range(n):
        a = rng.normal(size=d)
        a *= signal * sa / np.linalg.norm(a)
        f = m * a
        out.append(isotropic(f + rng.normal(scale=sf, size=d), a + rng.normal(scale=sa, size=d), sf, sa))
    return out


def _run(name, readings, **kw):
    if name == "ours30":
        return ours30(readings, C, readout=False, **kw)
    if name == "ours31":
        return ours31(readings, C, NU, **kw)
    if name == "ours32":
        return ours32(readings, C, NU, **kw)
    return ours33(readings, C, NU, **kw)


EQUATIONS = ["ours30", "ours31", "ours32", "ours33"]


@pytest.mark.parametrize("name", EQUATIONS)
def test_geometric_mean_is_the_law_functional_exp_mean_log_m(name):
    rng = np.random.default_rng(6)
    out = _run(name, _series(rng, 3.0, 4, 1.0))
    direct = math.exp(float(np.dot(out["law"], out["lm"])))
    assert out["geometric_mean"] == pytest.approx(direct, rel=1e-12)
    assert geometric_mean(out["lm"], out["law"]) == out["geometric_mean"]
    assert math.isfinite(out["geometric_mean"]) and out["geometric_mean"] > 0


@pytest.mark.parametrize("name", EQUATIONS)
def test_all_zero_experiment_gives_c_for_both_readouts(name):
    """T24: the all-zero law is pi, symmetric in ln m about ln c, so both centres are c."""
    readings = [isotropic([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.2, 0.1) for _ in range(3)]
    out = _run(name, readings)
    assert out["median"] == pytest.approx(C, rel=1e-6)
    assert out["geometric_mean"] == pytest.approx(C, rel=1e-6)


@pytest.mark.parametrize("name", EQUATIONS)
def test_a_zero_reading_leaves_the_geometric_mean_alone(name):
    """Z at the point: every law functional inherits it from the law (T24)."""
    rng = np.random.default_rng(7)
    readings = _series(rng, 3.0, 4, 1.0)
    padded = readings + [isotropic([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.2, 0.1)]
    a, b = _run(name, readings, grid=(8.0, 2001)), _run(name, padded, grid=(8.0, 2001))
    assert b["geometric_mean"] == pytest.approx(a["geometric_mean"], rel=1e-9)
    assert b["median"] == pytest.approx(a["median"], rel=1e-9)


@pytest.mark.parametrize("name", EQUATIONS)
def test_units_and_the_swap_act_on_the_geometric_mean_as_on_the_median(name):
    """L3: lambda/mu scales both; the same-centre swap (T5) sends both to c^2 / value."""
    rng = np.random.default_rng(8)
    readings = _series(rng, 3.0, 4, 1.0)
    base = _run(name, readings, grid=(8.0, 2001))

    lam, mu = 3.0, 0.5   # force in new units, acceleration in new units
    scaled = [Reading(lam * r.force, mu * r.acceleration,
                      np.diag(np.r_[np.full(3, lam**2), np.full(3, mu**2)]) * r.cov)
              for r in readings]
    c2 = C * lam / mu
    if name == "ours30":
        out = ours30(scaled, c2, readout=False, grid=(8.0, 2001))
    elif name == "ours31":
        out = ours31(scaled, c2, NU, grid=(8.0, 2001))
    elif name == "ours32":
        out = ours32(scaled, c2, NU, grid=(8.0, 2001))
    else:
        out = ours33(scaled, c2, NU, grid=(8.0, 2001))
    assert out["geometric_mean"] == pytest.approx(base["geometric_mean"] * lam / mu, rel=1e-9)
    assert out["median"] == pytest.approx(base["median"] * lam / mu, rel=1e-9)

    transform = np.block([[np.zeros((3, 3)), C * np.eye(3)], [np.eye(3) / C, np.zeros((3, 3))]])
    swapped = []
    for r in readings:
        y = transform @ np.r_[r.force, r.acceleration]
        swapped.append(Reading(y[:3], y[3:], transform @ r.cov @ transform.T))
    out = _run(name, swapped, grid=(8.0, 2001))
    assert out["geometric_mean"] == pytest.approx(C**2 / base["geometric_mean"], rel=1e-7)
    assert out["median"] == pytest.approx(C**2 / base["median"], rel=1e-7)


@pytest.mark.parametrize("name", EQUATIONS)
def test_the_two_readouts_agree_where_the_law_concentrates(name):
    rng = np.random.default_rng(9)
    out = _run(name, _series(rng, 3.0, 8, 20.0))
    assert out["geometric_mean"] == pytest.approx(out["median"], rel=2e-3)


@pytest.mark.parametrize("name", EQUATIONS)
def test_the_gap_at_weak_signal_is_measured_not_assumed(name):
    """The fork's size: where the law is lopsided in ln m the two centres part. Recorded as a
    bound, so a change that closes or widens the gap is seen; the numbers are in O2 (6 Nov)."""
    rng = np.random.default_rng(10)
    gaps = []
    for _ in range(20):
        out = _run(name, _series(rng, 3.0, 5, 0.7))
        gaps.append(abs(math.log(out["geometric_mean"] / out["median"])))
    gaps = np.array(gaps)
    assert gaps.max() > 0.01          # they do differ at weak signal
    assert np.median(gaps) < 0.5      # and by less than a factor e^0.5 on typical series
