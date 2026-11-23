"""The Gaussian score set widened by the pooled weight from pairs: scenes (problems/O14.md, 23 Nov;
tools/pooled_weight_on_the_score_set.py). Measurements, not a rule."""
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from contradiction_in_pairs import pair_weight          # noqa: E402

from met.conditional_score import gaussian_score_test   # noqa: E402
from met.full import isotropic                          # noqa: E402

CUT = stats.chi2.ppf(0.95, 1)


def _weight(pairs, sf, sa):
    return pair_weight([(np.atleast_1d(np.asarray(f, float)) / sf, np.atleast_1d(np.asarray(a, float)) / sa) for f, a in pairs])


def _kept(pairs, sf, sa, w, grid):
    rs = [isotropic(np.atleast_1d(f), np.atleast_1d(a), sf, sa) for f, a in pairs]
    return grid[np.array([gaussian_score_test(rs, m)["statistic"] for m in grid]) <= CUT / w]


def test_two_direct_readings_the_near_pair_untouched_the_far_pair_about_half_the_gap_per_sigma():
    grid = np.linspace(4.0, 7.0, 3001)
    near, far = [(5.5, 1.0), (5.6, 1.0)], [(5.0, 1.0), (6.0, 1.0)]
    assert _weight(near, 0.1, 1e-4) == 1.0
    w = _weight(far, 0.1, 1e-4)
    assert w == pytest.approx(0.0200, abs=3e-4)
    k = _kept(far, 0.1, 1e-4, w, grid)
    assert (k.min(), k.max()) == pytest.approx((4.521, 6.479), abs=3e-3)
    plain = _kept(far, 0.1, 1e-4, 1.0, grid)
    assert (plain.min(), plain.max()) == pytest.approx((5.362, 5.638), abs=3e-3)


def test_the_weight_is_unmoved_by_exact_zeros_and_by_tiny_readings_and_stops_following_sigma():
    far = [(5.0, 1.0), (6.0, 1.0)]
    w = _weight(far, 0.1, 0.1)
    assert _weight(far + [(0.0, 0.0)] * 40, 0.1, 0.1) == w
    assert _weight(far + [(1e-9, 1e-9)] * 40, 0.1, 0.1) == pytest.approx(w, rel=1e-9)
    assert _weight(far, 0.01, 0.01) == pytest.approx(w / 100, rel=0.02)      # so the widened set is the same set


def test_a_stray_reading_drags_the_plain_set_off_the_others_and_the_widened_set_still_holds_them():
    rng = np.random.default_rng(99)
    a = 1.0 + 0.02 * rng.normal(size=99)
    many = [(5.0 + 0.1 * rng.normal(), x) for x in a] + [(10.0, 1.0)]
    grid = np.linspace(4.5, 5.8, 1301)
    w = _weight(many, 0.1, 0.02)
    assert w < 0.05
    assert _kept(many, 0.1, 0.02, 1.0, grid).min() > 5.0 and _kept(many, 0.1, 0.02, w, grid).min() < 5.0


def test_a_reading_with_force_against_acceleration_keeps_only_the_mass_it_fits_worst():
    """The sign condition (18 Nov): no pair is contradicted, so the weight is 1; and the exact Gaussian score
    test rejects every positive mass except the one where the curve is highest, where the score is zero
    (E14's hilltop). A statement that keeps only the worst mass is not an answer."""
    one = [(10.0, -10.0)]
    assert _weight(one, 1.0, 1.0) == 1.0
    grid = np.exp(np.linspace(math.log(1e-3), math.log(1e3), 601))
    kept = _kept(one, 1.0, 1.0, 1.0, grid)
    assert len(kept) <= 1 and all(abs(m - 1.0) < 1e-9 for m in kept)
    from met.evidence import Curve
    curve = Curve([isotropic([10.0], [-10.0], 1.0, 1.0)], scale=1.0)
    assert curve(1.0) == pytest.approx(200.0) and curve(0.0) == pytest.approx(100.0) and curve(math.inf) == pytest.approx(100.0)
