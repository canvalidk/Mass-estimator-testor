"""What has to be carried from a block of readings so that blocks can be combined: past one
instrument the evidence curve alone is not enough (T77 (6); problems/O20.md).

Harvested 23 Nov 2026 from the exploration branch `explore/e14-ours34-walls` (Antigravity with Can,
22 Nov; its Theorem E14-3) and checked here with met.conditional_score and met.evidence, which
share no code with that branch.
"""
import math

import numpy as np
import pytest
from scipy import stats

from met.conditional_score import gaussian_score_test
from met.evidence import Curve
from met.full import Reading, isotropic


def test_two_blocks_with_one_curve_give_opposite_verdicts_beside_the_same_reading():
    """Forty readings on each of two instruments (noise ratio 1 and 2), each block a force-only and an
    acceleration-only reading of equal size in noise units, twenty times. Both curves are the
    constant 20 at every mass. Beside one further reading (F = 5, a = 1, unit noise), the exact
    Gaussian score test of mass 1 keeps it with the first block and rejects it with the second."""
    first = [isotropic([1.0], [0.0], 1.0, 1.0), isotropic([0.0], [1.0], 1.0, 1.0)] * 20
    second = [isotropic([2.0], [0.0], 2.0, 1.0), isotropic([0.0], [1.0], 2.0, 1.0)] * 20
    new = [isotropic([5.0], [1.0], 1.0, 1.0)]
    masses = np.array([0.0, 0.1, 0.5, 1.0, 2.0, 7.0, 50.0, math.inf])
    assert Curve(first, scale=1.0)(masses) == pytest.approx(20.0, abs=1e-12)
    assert Curve(second, scale=1.0)(masses) == pytest.approx(20.0, abs=1e-12)

    cut = stats.chi2.ppf(0.95, 1)
    with_first = gaussian_score_test(first + new, 1.0)["normal_score"] ** 2
    with_second = gaussian_score_test(second + new, 1.0)["normal_score"] ** 2
    assert with_first == pytest.approx(3.7895, abs=1e-4) and with_first < cut
    assert with_second == pytest.approx(4.6753, abs=1e-4) and with_second > cut
    # the reading on its own rejects mass 1 outright: both blocks, whose curves are flat, dilute it
    assert gaussian_score_test(new, 1.0)["normal_score"] ** 2 == pytest.approx(8.0)
    # what differs is the variance of the score, which the curve does not show
    assert gaussian_score_test(first, 1.0)["variance"] == pytest.approx(20.0)
    assert gaussian_score_test(second, 1.0)["variance"] == pytest.approx(12.8)


def test_the_score_is_minus_half_the_slope_of_the_curve_on_the_full_type():
    # any covariance, declared subspaces of every dimension, 1D to 3D: so the curve carries the score
    rng = np.random.default_rng(3)
    for _ in range(25):
        d = int(rng.integers(1, 4))
        readings = []
        for _ in range(int(rng.integers(1, 5))):
            k = int(rng.integers(1, d + 1))
            q, _ = np.linalg.qr(rng.normal(size=(d, d)))
            mix = rng.normal(size=(2 * d, 2 * d))
            cov = mix @ mix.T + 0.3 * np.eye(2 * d)
            readings.append(Reading(2 * rng.normal(size=d), rng.normal(size=d), 0.5 * (cov + cov.T), subspace=q[:, :k]))
        m = float(np.exp(rng.normal()))
        u, h = math.atan(m), 1e-5
        curve = Curve(readings, scale=1.0)
        slope = (curve.at_angle(u + h) - curve.at_angle(u - h)) / (2 * h)
        score = gaussian_score_test(readings, m)["score"]
        assert score == pytest.approx(-0.5 * slope, rel=1e-6, abs=1e-6)


def test_on_one_instrument_the_scores_variance_is_the_energy_less_the_curve():
    """In the instrument's own angle the conditional variance of the score is E - K, the part of the
    readings' energy the candidate does explain, and E is twice the mean level of the curve. So one
    instrument's curve, its level included, carries the variance; the noise ratio enters only
    through the change of coordinate, which is why two instruments need two curves."""
    rng = np.random.default_rng(4)
    sf, sa = 0.3, 0.7
    readings = [isotropic(rng.normal(size=3), rng.normal(size=3), sf, sa) for _ in range(5)]
    energy = sum(float(r.force @ r.force) / sf ** 2 + float(r.acceleration @ r.acceleration) / sa ** 2 for r in readings)
    s = sf / sa
    curve = Curve(readings, scale=s)
    assert 2 * float(np.mean(curve.at_angle(np.linspace(0.0, math.pi, 2000, endpoint=False)))) == pytest.approx(energy, rel=1e-9)
    for m in (0.2, 1.0, 3.0):
        jacobian = s / (s * s + m * m) * (1 + m * m)      # d(instrument angle) / d(atan m)
        assert gaussian_score_test(readings, m)["variance"] == pytest.approx((energy - curve(m)) * jacobian ** 2, rel=1e-9)
