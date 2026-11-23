"""ours34 (met/ours34.py, taken across from explore/e14-ours34-walls on 23 Nov 2026, as it is)
against main's own code on the full reading type: its curve against met.evidence, its Gaussian
score against met.conditional_score, its point against the curve's least. Claude's check; the
probe itself is not otherwise checked beyond its own tests (tests/test_ours34.py)."""
import math

import numpy as np
import pytest

from met.conditional_score import gaussian_score_test
from met.evidence import Curve, leftover
from met.full import Reading
from met.ours34 import State, ours34


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


def test_ours34_curve_score_and_point_agree_with_main_on_the_full_type():
    rng = np.random.default_rng(34)
    for _ in range(12):
        readings = _series(rng)
        state, curve = State(readings), Curve(readings)
        for m in (0.3, 1.0, 2.0, 6.0):
            k, s, c = state.ksc_m(m)
            assert k == pytest.approx(curve(m), rel=1e-8, abs=1e-8)
            ref = gaussian_score_test(readings, m)
            assert s / math.sqrt(c) == pytest.approx(ref["normal_score"], rel=1e-6, abs=1e-7)
        assert ours34(readings)["point"] == pytest.approx(leftover(readings)["mass"], rel=1e-5)


def test_on_two_scalar_readings_ours34_widens_the_agreeing_pair_too():
    """Can's two scenes of 17 Nov as he gave them: two scalar readings of one mass, each +- 0.1.
    The near pair (5.5, 5.6), which he called convincing, is kept as 4.91 to 6.19, where the
    declared Gaussian gives 5.41 to 5.69. The far pair (5, 6) is kept from 0 to 11.85, with a
    second piece at very large masses. A finding against R-2026-11-17-3's near scene."""
    from met.ours34 import scalar
    near = ours34([scalar(5.5, 1.0, 0.1, 1e-6), scalar(5.6, 1.0, 0.1, 1e-6)])
    assert near["point"] == pytest.approx(5.55, abs=1e-6)
    assert near["set"][0] == pytest.approx((4.9147, 6.1853), abs=1e-3)
    assert near["raw_set"][0] == pytest.approx((5.4114, 5.6886), abs=1e-3)
    far = ours34([scalar(5.0, 1.0, 0.1, 1e-6), scalar(6.0, 1.0, 0.1, 1e-6)])
    assert far["set"][0] == pytest.approx((0.0, 11.853), abs=1e-2) and len(far["set"]) == 2
