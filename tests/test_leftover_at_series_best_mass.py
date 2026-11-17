"""Where a reading's leftover is read (Claude, 17 Nov 2026; problems/O14.md). One-number readings."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import leftover_at_series_best_mass as lb  # noqa: E402


def _half(result):
    return 0.5 * (result[2] - result[1])


def test_the_sign_problem_comes_from_reading_the_leftover_at_the_readings_own_best_positive_mass():
    # truth 2, signal 1 per reading: ours30's recorded 1.73 (G22) is reproduced by its way of reading the
    # leftover, and reading it at the series' best mass keeps the estimate at 2
    plain, at_series, own = lb.sign_check(1.0, n=600_000)
    assert plain == pytest.approx(2.0, abs=0.02) and at_series == pytest.approx(2.0, abs=0.02)
    assert own == pytest.approx(1.73, abs=0.03)


def test_the_rule_at_every_candidate_mass_does_not_leave_agreeing_readings_at_face_value():
    near = lb.cases()["5.5, 5.6"]
    assert _half(lb.law(near, kind="gauss")) == pytest.approx(0.139, abs=0.003)
    assert _half(lb.law(near, kind="every candidate")) > 0.5
    for kind in ("pooled", "per reading"):
        assert _half(lb.law(near, kind=kind)) == pytest.approx(0.139, abs=0.003)


def test_pooled_and_per_reading_agree_on_two_readings_and_part_ways_on_a_stray_among_many():
    c = lb.cases()
    for kind in ("pooled", "per reading"):
        assert _half(lb.law(c["5, 6"], kind=kind)) == pytest.approx(0.69, abs=0.02)
    pooled, per = lb.law(c["the same 99 and one at 9"], kind="pooled"), lb.law(c["the same 99 and one at 9"], kind="per reading")
    assert _half(pooled) > 3 * _half(per)
    assert _half(per) < 1.2 * _half(lb.law(c["99 ordinary near 5"], kind="gauss"))
