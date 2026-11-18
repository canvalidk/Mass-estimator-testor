"""The weak-signal share and the fit; the score law's false peak (Claude, 18 Nov 2026; problems/O14.md, 18 Nov,
afternoon; tools/share_and_the_fit.py). Measurements, no rule adopted."""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import share_and_the_fit as sf  # noqa: E402


def test_at_the_true_line_the_share_has_no_count_and_a_zero_reading_changes_nothing():
    # E[A] = Delta and E[T] = Delta + M at the true line, whatever the pushes
    rng = np.random.default_rng(3)
    k, n, a_sum, t_sum, d_sum = 3, 12, 0.0, 0.0, 0.0
    for _ in range(4000):
        F, A, delta = sf.series(rng, k, n, 1.0)
        t, r = (F + A) / math.sqrt(2), (F - A) / math.sqrt(2)
        a_sum += (t * t).sum() - (r * r).sum(); t_sum += (t * t).sum(); d_sum += delta
    assert a_sum / d_sum == pytest.approx(1.0, abs=0.03)
    assert t_sum / (d_sum + 4000 * n * k) == pytest.approx(1.0, abs=0.01)
    F, A, _ = sf.series(rng, k, n, 1.0)
    padded = sf.phasor(np.vstack([F, np.zeros((5, k))]), np.vstack([A, np.zeros((5, k))]))
    assert padded == pytest.approx(sf.phasor(F, A))


def test_the_gram_share_is_never_below_the_true_line_share_on_any_series():
    # the fit puts A at its largest, |Z|, and x / (E/2 + x/2) increases with x: the overshoot is one-sided
    rng = np.random.default_rng(8)
    for k, n, snr in ((1, 5, 0.5), (3, 20, 1.0), (1, 100, 2.0)):
        for _ in range(300):
            F, A, _ = sf.series(rng, k, n, snr)
            absZ, _, energy = sf.phasor(F, A)
            assert sf.gram_share(absZ, energy) >= sf.true_line_share(F, A) - 1e-12


def test_the_recorded_gram_means_are_reproduced_and_the_first_order_correction_removes_most_of_the_overshoot():
    weak = sf.share_cell(0, 1, 5, 0.5)                    # the 7 Nov table's first row: Gram 0.658, target 0.189
    assert weak["target"] == pytest.approx(0.189, abs=0.002) and weak["gram"] == pytest.approx(0.658, abs=0.002)
    assert weak["first_order"] == pytest.approx(0.32, abs=0.02) and weak["under_the_law"] == pytest.approx(0.28, abs=0.02)
    assert weak["rms_first_order"] < 0.7 * weak["rms_gram"]
    mid = sf.share_cell(5, 1, 20, 1.0)                    # at moderate signal the correction costs scatter
    assert mid["gram"] == pytest.approx(0.542, abs=0.002) and mid["rms_first_order"] > mid["rms_gram"]
    assert weak["several_solutions"] >= 1                 # the self-consistent form is not always single-valued


def test_the_score_laws_false_peak_goes_when_signal_energy_cannot_be_negative():
    rows = sf.false_peak(2000, n_series=12, grid=100001)
    false_side = {name: np.mean([r["false_side"] for r in v]) for name, v in rows.items()}
    assert false_side["score law"] == pytest.approx(0.75, abs=0.03)          # T56's 0.7536
    assert false_side["signal energy >= 0"] < 1e-6 and false_side["gaussian"] < 1e-6
    med = np.median([r["median"] for r in rows["signal energy >= 0"]])
    assert med == pytest.approx(sf.M0, abs=0.02)
    # the local width is the Godambe one: wider than the Gaussian's by about sqrt((Delta + M)/Delta), r0^2 = 28/19
    width = lambda name: np.median([math.log(r["hi"] / r["lo"]) for r in rows[name]])
    assert width("signal energy >= 0") / width("gaussian") == pytest.approx(math.sqrt((28 / 19 + 1) / (28 / 19)), abs=0.08)
