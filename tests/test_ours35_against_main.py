"""ours35 (met/ours35.py, taken across from explore/e16-ours35-door-a on 24 Nov 2026, as it is;
R-2026-11-24-1) against main's own code on the full reading type: its curve against met.evidence,
its score against met.conditional_score, its point against the curve's least, its kept set against
main's score test mass by mass. Then three scenes read against the rows. Claude's check, 24 Nov;
its two-body sets (lump, common_force) are not checked here."""
import math

import numpy as np
import pytest

from met.conditional_score import gaussian_score_accepts, gaussian_score_test
from met.evidence import Curve, leftover
from met.full import Reading
from met.ours35 import State, ours35, scalar


def _series(rng):
    d = int(rng.integers(1, 4))
    out = []
    for _ in range(int(rng.integers(2, 6))):
        k = int(rng.integers(1, d + 1))
        q, _ = np.linalg.qr(rng.normal(size=(d, d)))
        a = q[:, :k] @ rng.normal(size=k) * 4.0
        mix = rng.normal(size=(2 * d, 2 * d))
        cov = 0.2 * (mix @ mix.T) / (2 * d) + 0.05 * np.eye(2 * d)
        y = np.concatenate([2.0 * a, a]) + np.linalg.cholesky(cov) @ rng.normal(size=2 * d)
        out.append(Reading(y[:d], y[d:], 0.5 * (cov + cov.T), subspace=q[:, :k]))
    return out


def _inside(intervals, m):
    return any(lo <= m <= hi for lo, hi in intervals)


def _near_an_end(intervals, m, rel=1e-3):
    ends = [e for pair in intervals for e in pair if 0.0 < e < math.inf]
    return any(abs(m - e) <= rel * e for e in ends)


def test_ours35_curve_score_point_and_set_agree_with_main_on_the_full_type():
    rng = np.random.default_rng(35)
    masses = np.geomspace(0.05, 80.0, 61)
    for _ in range(12):
        readings = _series(rng)
        state, curve = State(readings), Curve(readings)
        for m in (0.3, 1.0, 2.0, 6.0):
            k, s, c = state.ksc_m(m)
            assert k == pytest.approx(curve(m), rel=1e-8, abs=1e-8)
            ref = gaussian_score_test(readings, m)
            assert s / math.sqrt(c) == pytest.approx(ref["normal_score"], rel=1e-6, abs=1e-7)
            assert state.p_value(m) == pytest.approx(ref["p_value"], rel=1e-5, abs=1e-12)
        out = ours35(readings)
        assert out["point"] == pytest.approx(leftover(readings)["mass"], rel=1e-5)
        assert out["leftover"] == pytest.approx(leftover(readings)["least"], rel=1e-7, abs=1e-9)
        kept = gaussian_score_accepts(readings, masses, probability=0.95)
        for m, keep in zip(masses, kept):
            if not _near_an_end(out["set"], m):
                assert _inside(out["set"], m) == bool(keep)


def test_ours35_keeps_the_far_pair_as_narrowly_as_the_near_pair():
    """Can's two scenes of 17 Nov on two direct scalar readings, each +- 0.1. Both pairs are kept
    with the same width (0.2772); what differs is the leftover beside the set, 0.5 against 50.
    R-2026-11-17-3 asks for the far pair to be believed less in the result; this probe is built
    not to (it is the experiment without the weighting), so this row is a known fail."""
    near = ours35([scalar(5.5, 1.0, 0.1, 1e-6), scalar(5.6, 1.0, 0.1, 1e-6)])
    far = ours35([scalar(5.0, 1.0, 0.1, 1e-6), scalar(6.0, 1.0, 0.1, 1e-6)])
    (nlo, nhi), (flo, fhi) = near["set"][0], far["set"][0]
    assert (nlo, nhi) == pytest.approx((5.4114, 5.6886), abs=1e-4)
    assert (flo, fhi) == pytest.approx((5.3614, 5.6386), abs=1e-4)
    assert nhi - nlo == pytest.approx(fhi - flo, abs=1e-8)
    assert near["leftover"] == pytest.approx(0.5, abs=1e-6)
    assert far["leftover"] == pytest.approx(50.0, abs=1e-5)


def test_on_readings_against_themselves_the_set_leaves_out_the_point_and_keeps_the_worst_mass():
    """Two scalar readings with force against acceleration, unit noise: (10, -3) and (10, -3.1).
    The curve is least at infinite mass (18.61), so that is the point (T81). The score test
    rejects it (p = 1.6e-5) and keeps only a sliver round the mass the curve is HIGHEST at
    (0.305, where the curve is 218.6). Recomputed on main's curve and score test: the same.
    This answers T81's "not shown" (whether the end lies inside the kept set when the least is
    at an end): on this input it does not. From exploration E15 (its E15-3), checked here."""
    readings = [scalar(10.0, -3.0, 1.0, 1.0), scalar(10.0, -3.1, 1.0, 1.0)]
    out = ours35(readings)
    assert out["point"] == math.inf and out["leftover"] == pytest.approx(18.61, abs=1e-9)
    assert len(out["set"]) == 1
    lo, hi = out["set"][0]
    assert 0.304 < lo < 0.305 < hi < 0.306
    assert out["state"].p_value(math.inf) == pytest.approx(1.608e-5, rel=1e-3)

    main = [Reading([10.0], [-3.0], np.eye(2)), Reading([10.0], [-3.1], np.eye(2))]
    curve = Curve(main)
    top = max(float(curve(m)) for m in np.geomspace(1e-3, 1e3, 4001))
    assert curve(0.305) == pytest.approx(218.605, abs=1e-2) and curve(0.305) >= top - 1e-3
    assert curve(1e8) == pytest.approx(18.61, abs=1e-4)
    assert leftover(main)["at_end"] and leftover(main)["least"] == pytest.approx(18.61, abs=1e-6)
    assert gaussian_score_test(main, 0.305)["p_value"] > 0.95
    assert gaussian_score_test(main, 1e8)["p_value"] == pytest.approx(1.608e-5, rel=1e-3)
    assert gaussian_score_test(main, 1.0)["p_value"] < 1e-30
