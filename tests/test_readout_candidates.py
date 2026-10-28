"""O2's Gaussian candidate formulas, null behavior and representation symmetries."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from met.full import Reading, isotropic, ours29
from met.readout_candidates import CANDIDATE_IDS, candidate_weights, evaluate


def _reference(c):
    if c == 1:
        return 1.0
    if c < 1:
        return c * math.acosh(1 / c) / math.acos(c)
    return c * math.acos(1 / c) / math.acosh(c)


def _points(result):
    return {key: value["point"] for key, value in result["candidates"].items()}


@pytest.mark.parametrize("c", [.1, 1.0, 10.0])
def test_perpendicular_size_rules_share_a_shape_but_global_rms_is_distinct(c):
    reading = isotropic([1., 0.], [0., 1.], .3, .3)
    result = evaluate([reading], c)
    points = _points(result)
    expected = _reference(c)
    for name in ("cloud_mean", "mean_vector", "conditional_rms", "plain_sum_cloud"):
        assert points[name] == pytest.approx(expected, rel=2e-6)
    assert points["global_rms"] == pytest.approx(math.sqrt(c), rel=3e-6)
    assert points["law_median"] == pytest.approx(c, rel=2e-7)
    # A separate compact-angle integral checks the shared expression.
    numerator = quad(lambda phi: c * math.sin(phi) / (c * c * math.cos(phi) ** 2 + math.sin(phi) ** 2),
                     0, math.pi / 2, epsabs=1e-12)[0]
    denominator = quad(lambda phi: c * math.cos(phi) / (c * c * math.cos(phi) ** 2 + math.sin(phi) ** 2),
                       0, math.pi / 2, epsabs=1e-12)[0]
    assert expected == pytest.approx(numerator / denominator, rel=2e-12)


def test_perpendicular_weights_match_the_gaussian_conditional_moments():
    sigma = .3
    masses = np.array([.05, .2, 1., 5., 20.])
    result = candidate_weights([isotropic([1., 0.], [0., 1.], sigma, sigma)], masses)
    scaled = math.sqrt(1) / np.sqrt(1 + masses ** 2)
    np.testing.assert_allclose(result["weights"]["mean_vector"], scaled, rtol=1e-13)
    np.testing.assert_allclose(result["weights"]["conditional_rms"], math.sqrt(1 + 2 * sigma ** 2) * scaled,
                               rtol=1e-13)
    normalized_cloud = result["weights"]["cloud_mean"] / scaled
    np.testing.assert_allclose(normalized_cloud, np.full_like(masses, normalized_cloud[0]), rtol=2e-12)
    assert np.all(result["weights"]["mean_vector"] <= result["weights"]["cloud_mean"] + 1e-13)
    assert np.all(result["weights"]["cloud_mean"] <= result["weights"]["conditional_rms"] + 1e-13)


@pytest.mark.parametrize("c", [.1, 1., 10.])
def test_all_zero_mean_candidate_is_undefined_and_others_do_not_invent_c(c):
    result = evaluate([isotropic([0., 0.], [0., 0.], 1, 1)], c)
    assert result["candidates"]["mean_vector"] == {"status": "undefined", "point": None}
    for name in ("cloud_mean", "conditional_rms", "plain_sum_cloud"):
        assert result["candidates"][name]["point"] == pytest.approx(_reference(c), rel=2e-6)
    assert result["candidates"]["global_rms"]["point"] == pytest.approx(math.sqrt(c), rel=3e-6)
    assert result["candidates"]["law_median"]["point"] == pytest.approx(c, rel=2e-7)


def test_exact_null_padding_preserves_mean_candidate_and_law_but_moves_other_sizes():
    reading = isotropic([2., 0.], [1., 0.], 1, 1)
    zero = isotropic([0., 0.], [0., 0.], 1, 1)
    before, after = evaluate([reading], 1), evaluate([reading] + [zero] * 8, 1)
    np.testing.assert_allclose(before["law"]["log_mass"], after["law"]["log_mass"], atol=1e-12)
    np.testing.assert_allclose(before["law"]["probability"], after["law"]["probability"], atol=1e-12)
    assert before["law"]["eta"] == pytest.approx(after["law"]["eta"], abs=1e-12)
    assert before["candidates"]["mean_vector"]["point"] == pytest.approx(after["candidates"]["mean_vector"]["point"], rel=1e-12)
    assert before["candidates"]["law_median"]["point"] == pytest.approx(after["candidates"]["law_median"]["point"], rel=1e-12)
    for name in ("cloud_mean", "conditional_rms", "plain_sum_cloud", "global_rms"):
        ratio = after["candidates"][name]["point"] / before["candidates"][name]["point"]
        assert abs(math.log(ratio)) > 1e-4


def test_perpendicular_different_ratio_null_padding_has_exact_global_rms_shift():
    base = isotropic([1., 0.], [0., 1.], 1, 1)
    zero = isotropic([0., 0.], [0., 0.], 3, 1)
    before, after = evaluate([base], 1), evaluate([base, zero], 1)
    np.testing.assert_allclose(before["law"]["log_mass"], after["law"]["log_mass"], atol=1e-12)
    np.testing.assert_allclose(before["law"]["probability"], after["law"]["probability"], atol=1e-12)
    masses = np.exp(after["law"]["log_mass"])
    expected_second_moment = 3 / (1 + masses ** 2) + 18 / (9 + masses ** 2)
    np.testing.assert_allclose(after["second_moment"], expected_second_moment, rtol=2e-13)
    np.testing.assert_allclose(before["weights"]["mean_vector"], after["weights"]["mean_vector"],
                               rtol=2e-13)
    assert before["candidates"]["mean_vector"]["point"] == pytest.approx(1, rel=2e-6)
    assert after["candidates"]["mean_vector"]["point"] == pytest.approx(1, rel=2e-6)
    # Under the standard half-Cauchy law, E[1/(a²+m²)] = 1/[a(a+1)]
    # and E[m²/(a²+m²)] = 1/(a+1). The denominator is 3 and numerator 6.
    assert after["candidates"]["global_rms"]["point"] == pytest.approx(math.sqrt(2), rel=3e-6)
    # The added variance has a strictly increasing size ratio relative to the
    # perpendicular base, so these candidates reweight toward larger masses.
    for name in ("conditional_rms", "cloud_mean", "plain_sum_cloud"):
        assert after["candidates"][name]["point"] > before["candidates"][name]["point"]


def test_nonzero_observations_outside_the_declared_span_leave_mean_undefined():
    reading = Reading([0., 1.], [0., 1.], np.eye(4), subspace=[[1.], [0.]])
    result = evaluate([reading], 1)
    assert result["candidates"]["mean_vector"] == {"status": "undefined", "point": None}
    np.testing.assert_array_equal(result["weights"]["mean_vector"], 0)
    assert result["law"]["median"] == pytest.approx(1, rel=2e-7)


def test_cloud_row_is_exactly_the_production_value_and_law():
    readings = [isotropic([2., .2], [1., -.1], .6, .8)]
    result, production = evaluate(readings, 1.3), ours29(readings, 1.3)
    assert tuple(result["candidates"]) == CANDIDATE_IDS
    assert result["candidates"]["cloud_mean"]["point"] == production["readout"]
    np.testing.assert_array_equal(result["law"]["probability"], production["law"])
    np.testing.assert_array_equal(result["weights"]["cloud_mean"], production["A"])


def test_candidates_commute_with_channel_swap_and_units_with_full_covariance():
    factor = np.eye(6) + .05 * np.tril(np.ones((6, 6)), -1)
    covariance = .3 * factor @ factor.T
    original = Reading([2., .3, -.1], [1., -.2, .1], covariance)
    indices = [3, 4, 5, 0, 1, 2]
    swapped = Reading(original.acceleration, original.force, covariance[np.ix_(indices, indices)])
    scale = np.diag([3.] * 3 + [5.] * 3)
    units = Reading(3 * original.force, 5 * original.acceleration, scale @ covariance @ scale)
    base = _points(evaluate([original], 1.7))
    inverse, changed = _points(evaluate([swapped], 1 / 1.7)), _points(evaluate([units], 1.7 * 3 / 5))
    for name in CANDIDATE_IDS:
        assert base[name] * inverse[name] == pytest.approx(1, rel=3e-6)
        assert changed[name] / base[name] == pytest.approx(3 / 5, rel=3e-6)


def test_all_candidates_approach_the_compatible_strong_signal_mass():
    result = evaluate([isotropic([6., 0., 0.], [2., 0., 0.], .002, .002)], 1)
    for candidate in result["candidates"].values():
        assert candidate["status"] == "defined"
        assert candidate["point"] == pytest.approx(3, rel=4e-4)


def test_candidates_are_invariant_under_a_common_frame_change():
    factor = np.eye(6) + .05 * np.tril(np.ones((6, 6)), -1)
    covariance = .3 * factor @ factor.T
    original = Reading([2., .3, -.1], [1., -.2, .1], covariance,
                       subspace=[[1., 0.], [0., 1.], [0., 0.]])
    angle = .7
    rotation = np.array([[math.cos(angle), 0., math.sin(angle)],
                         [0., 1., 0.], [-math.sin(angle), 0., math.cos(angle)]])
    zero = np.zeros((3, 3))
    stacked_rotation = np.block([[rotation, zero], [zero, rotation]])
    rotated = Reading(rotation @ original.force, rotation @ original.acceleration,
                      stacked_rotation @ covariance @ stacked_rotation.T,
                      subspace=rotation @ original.subspace)
    base, changed = _points(evaluate([original], 1.7)), _points(evaluate([rotated], 1.7))
    for name in CANDIDATE_IDS:
        assert changed[name] == pytest.approx(base[name], rel=3e-6)


@pytest.mark.parametrize("masses", [[0, 1], [1, float("nan")], [], [[1, 2]]])
def test_weight_grid_rejects_invalid_masses(masses):
    with pytest.raises(ValueError):
        candidate_weights([isotropic([1.], [1.], 1, 1)], masses)


def test_empty_comparison_and_invalid_prior_are_not_silently_patched():
    with pytest.raises(ValueError):
        evaluate([], 1)
    with pytest.raises(ValueError):
        evaluate([isotropic([1.], [1.], 1, 1)], True)
