"""Measurements on ways to read the leftover (Claude, 18 Nov 2026; problems/O14.md, 18 Nov; tools/discount_lab.py,
tools/discount_scalar_slice.py). They record what each body does on the scenes; no rule is adopted."""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import discount_lab as lab  # noqa: E402
import discount_scalar_slice as ss  # noqa: E402

GRID = 40001
SCENES = lab.scenes()


def _law(name, sigma, body):
    readings, _, c, ratio = SCENES[name]
    return lab.law(readings, sigma, sigma * ratio, c, body, grid=GRID)


def _bounded(name, body):
    a, b = _law(name, 0.05, body), _law(name, 0.01, body)
    return b["hi"] / b["lo"] == pytest.approx(a["hi"] / a["lo"], rel=0.05) and b["point"] == pytest.approx(a["point"], rel=0.05)


C, C1, C2 = (n for n in SCENES if n.startswith("C"))


def test_with_the_fit_corrected_count_the_far_pair_is_the_scatter_of_the_two_readings():
    # expected leftover at the fitted mass is N - 1, not N: weight (s1^2 + s2^2)/gap^2 = 0.02, half-width 0.98, where
    # the 17 Nov table (count N) has 0.04 and 0.69
    r = _law("B: 5, 6", 0.1, "one weight, by count")
    assert r["w"][0] == pytest.approx(0.0202, abs=0.0003)
    assert 0.5 * (r["hi"] - r["lo"]) == pytest.approx(0.98, abs=0.01)
    x, s = np.array([5.0, 6.0]), np.array([0.05, 0.1])
    assert ss.one_weight(x, s)[2] == pytest.approx(0.0125) and ss.per_reading(x, s)[2] == pytest.approx([0.0125, 0.0125])
    assert ss.one_weight(x, s, count="N")[2] == pytest.approx(0.0250)


def test_the_expected_leftover_at_the_fitted_mass_is_one_minus_the_leverage():
    sim, formula = ss.expected_leftover_at_the_fit(reps=200000)
    assert sim == pytest.approx(formula, abs=0.01)


@pytest.mark.parametrize("body", ("gauss", "eta as on the pages"))
def test_the_gaussian_and_eta_are_not_bounded_on_any_right_angle_scene(body):
    assert _law(C, 0.01, body)["point"] > 200 and _law(C1, 0.01, body)["point"] > 300
    r = _law(C2, 0.01, body)
    assert r["hi"] / r["lo"] < 1.05


@pytest.mark.parametrize("body", ("one weight, by count", "one weight, by energy", "one weight, by signal energy"))
def test_these_one_weight_bodies_hold_every_right_angle_scene(body):
    for name in (C, C1, C2, "one 2D reading F = (2, 0), a = (0, 1)", "anti-parallel scalar reading F = 1, a = -1"):
        assert _bounded(name, body)


def test_by_leverage_is_blind_to_a_reading_at_right_angles_and_undefined_on_the_pair():
    assert _law(C1, 0.01, "one weight, by leverage")["w"][0] == 1.0
    assert _law(C2, 0.01, "one weight, by leverage")["w"][0] == 1.0
    assert _law(C, 0.2, "one weight, by leverage")["undefined"]


def test_per_reading_leaves_the_agreeing_readings_at_face_value_on_the_right_angle_scenes():
    for name in (C1, C2):
        assert _law(name, 0.01, "one weight per reading")["w"][:2] == pytest.approx([1.0, 1.0])
    r = _law(C2, 0.01, "one weight per reading")
    assert r["hi"] / r["lo"] < 1.05


def test_noise_only_readings_dilute_the_count_and_exact_zeros_dilute_it_more():
    base = _law("B: 5, 6", 0.1, "one weight, by count")["w"][0]
    assert _law("B and 40 noise-only readings", 0.1, "one weight, by count")["w"][0] > 15 * base
    assert _law("B and 40 exactly zero readings", 0.1, "one weight, by count")["w"][0] > 35 * base


@pytest.mark.parametrize("body", ("one weight, by leverage", "one weight, by energy", "one weight, by signal energy"))
def test_the_weighted_bodies_are_not_diluted_by_noise_only_or_zero_readings(body):
    base = _law("B: 5, 6", 0.1, body)["w"][0]
    assert _law("B and 40 noise-only readings", 0.1, body)["w"][0] == pytest.approx(base, rel=0.02)
    assert _law("B and 40 exactly zero readings", 0.1, body)["w"][0] == pytest.approx(base, rel=1e-9)


def test_by_signal_energy_keeps_ordinary_scenes_at_face_value_and_an_aligned_pair_concentrating():
    for name in ("A: 5.5, 5.6 (each +- 0.1)", "99 ordinary near 5"):
        assert _law(name, 0.1, "one weight, by signal energy")["w"][0] == 1.0
    wide, narrow = _law("aligned pair (2, 1) and (4, 2)", 0.2, "one weight, by signal energy"), _law("aligned pair (2, 1) and (4, 2)", 0.01, "one weight, by signal energy")
    assert narrow["w"][0] == 1.0 and math.log(narrow["hi"] / narrow["lo"]) < 0.1 * math.log(wide["hi"] / wide["lo"])


def test_a_zero_axis_inside_a_reading_still_counts():
    # T25's k-dependence: the 3D reading with an exactly zero third axis is discounted half as much as the 2D reading
    for body in ("one weight, by count", "one weight, by signal energy"):
        two = _law("one 2D reading F = (2, 0), a = (0, 1)", 0.01, body)["w"][0]
        three = _law("the same with an exactly zero third axis", 0.01, body)["w"][0]
        assert three == pytest.approx(2 * two, rel=0.01)


def test_unequal_weights_can_leave_the_law_more_concentrated_than_no_discount():
    # three scalar readings in noise units; down-weighting the dissenter lengthens the others' resultant
    ws = [(np.array([2.0]), np.array([2.0])), (np.array([2.2]), np.array([1.8])), (np.array([-1.6]), np.array([1.2]))]
    plain = lab.swing(ws, np.ones(3))
    per, _ = lab.weights(ws, "one weight per reading")
    one, _ = lab.weights(ws, "one weight, by count")
    assert plain == pytest.approx(12.382, abs=0.001)
    assert lab.swing(ws, per) == pytest.approx(15.032, abs=0.001) and per.max() == 1.0
    assert lab.swing(ws, one) == pytest.approx(6.434, abs=0.001)


def test_on_ordinary_weak_data_per_reading_weights_usually_sharpen_the_law():
    rng = np.random.default_rng(77)
    sharper = 0
    for _ in range(200):
        ws = [(rng.normal(size=3), rng.normal(size=3)) for _ in range(20)]        # 3D, no signal
        w, _ = lab.weights(ws, "one weight per reading")
        sharper += lab.swing(ws, np.where(np.isnan(w), 1.0, w)) > lab.swing(ws, np.ones(20))
    assert sharper > 140


def test_lumping_with_one_weight_the_glued_estimate_is_the_sum_of_the_parts_and_per_reading_it_is_not():
    F1, F2 = [5.0, 5.05, 6.0], [3.0, 3.9, 3.05]
    one, per = ss.glue(F1, F2, 0.1, 0.1, ss.one_weight), ss.glue(F1, F2, 0.1, 0.1, ss.per_reading)
    assert one["glued"] == pytest.approx(one["parts"], abs=1e-12)
    assert per["glued"] - per["parts"] == pytest.approx(0.623, abs=0.002)
    assert np.sqrt((ss.lumping_on_ordinary_data(ss.per_reading, reps=4000) ** 2).mean()) == pytest.approx(0.40, abs=0.03)


def test_per_reading_iterated_to_self_consistency_has_two_answers():
    x, s = np.array([5.0, 5.05, 6.0]), np.full(3, 0.1)
    assert ss.per_reading(x, s, steps=500)[0] == pytest.approx(5.028, abs=0.002)
    assert ss.per_reading(x, s, steps=500, start=6.0)[0] == pytest.approx(5.986, abs=0.002)


def test_by_energy_discounts_ordinary_weak_data_and_by_signal_energy_much_less():
    # 200 scalar readings at signal 1, the declared sigma true: weighing by the whole energy counts noise as signal
    share, median = lab.ordinary_share("one weight, by energy", 1, 200, 1.0, reps=60)
    assert share == 1.0 and median == pytest.approx(0.63, abs=0.05)
    share, median = lab.ordinary_share("one weight, by signal energy", 1, 200, 1.0, reps=60)
    assert share < 0.7 and median > 0.8


def test_by_signal_energy_with_thousands_of_noise_only_or_zero_readings():
    rng = np.random.default_rng(3)
    far = [([5.0], [1.0]), ([6.0], [1.0])]
    nulls = [([0.1 * rng.normal()], [0.002 * rng.normal()]) for _ in range(4000)]
    body = "one weight, by signal energy"
    base = lab.law(far, 0.1, 0.002, 5.0, body, grid=2001)["w"][0]
    assert base == pytest.approx(0.0202, abs=0.0002)
    assert lab.law(far + nulls, 0.1, 0.002, 5.0, body, grid=2001)["w"][0] == pytest.approx(base, rel=0.06)
    assert lab.law(far + [([0.0], [0.0])] * 4000, 0.1, 0.002, 5.0, body, grid=2001)["w"][0] == pytest.approx(base, rel=1e-9)


OWN = "one weight, energy less 2, own best line"


def test_energy_less_two_at_its_own_best_line_holds_the_scenes_and_keeps_ordinary_ones_at_face_value():
    for name in (C, C1, C2, "one 2D reading F = (2, 0), a = (0, 1)", "anti-parallel scalar reading F = 1, a = -1"):
        assert _bounded(name, OWN)
    for name in ("A: 5.5, 5.6 (each +- 0.1)", "99 ordinary near 5"):
        assert _law(name, 0.1, OWN)["w"][0] == 1.0
    assert _law("B: 5, 6", 0.1, OWN)["w"][0] == pytest.approx(0.0202, abs=0.0003)
    assert _law("B and 40 noise-only readings", 0.1, OWN)["w"][0] == pytest.approx(0.0202, abs=0.0004)
    assert _law("one 2D reading F = (2, 0), a = (0, 1)", 0.01, OWN)["w"][0] == pytest.approx(1e-4, rel=0.01)   # T30: (k - 1)/Q


def test_exactly_zero_vector_readings_change_the_signal_energy_weight_and_not_the_energy_less_two_one():
    # correction to item 5 of O14's first 18 Nov section: with weight C - 2k the expected side gains k (2 - 2k) per zero
    # reading when k > 1; with weight C - 2 it gains nothing, for every k
    far = [([5.0, 0.3, 0.0], [1.0, 0.0, 0.0]), ([6.0, -0.2, 0.1], [1.0, 0.0, 0.0])]
    zeros = [([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])] * 40

    def weight(readings, body):
        return lab.law(readings, 0.1, 0.002, 5.0, body, grid=2001)["w"][0]

    assert weight(far + zeros, OWN) == pytest.approx(weight(far, OWN), rel=1e-12)
    old = "one weight, by signal energy"
    assert weight(far + zeros, old) / weight(far, old) == pytest.approx(0.9996, abs=0.0002)
    assert weight(far + zeros, old) != pytest.approx(weight(far, old), rel=1e-6)
    assert weight(far + zeros, "one weight, by count") == 1.0
