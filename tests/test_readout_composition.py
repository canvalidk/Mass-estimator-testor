"""O2 composition checks use propagated candidate laws, never fresh experiments."""
import json
import math

import numpy as np
import pytest

from met.readout_law_checks import (cancellation_counterexample, counterexamples,
                                   direction_counterexample, magnitude_readout, vector_readout)


@pytest.mark.parametrize("candidate,left,right", [
    ("cloud_mean", 4, 4), ("law_median", 5, 2),
    ("conditional_rms", 1 + math.sqrt(5), 13 / 4),
    ("global_rms", 4, 2 * math.sqrt(5)), ("mean_vector", 2, 10 / 3),
])
def test_exact_finite_law_composition_values(candidate, left, right):
    record = next(row for row in counterexamples() if row["candidate"] == candidate)
    assert record["L5_left"] == pytest.approx(left, abs=1e-13)
    assert record["L5_right"] == pytest.approx(right, abs=1e-13)
    assert record["L14_inverse_left"] == pytest.approx(left, abs=1e-13)
    assert record["L14_inverse_right"] == pytest.approx(right, abs=1e-13)
    if candidate == "cloud_mean":
        assert record["L5_absolute_error"] < 1e-13
        assert record["L14_absolute_error"] < 1e-13
    else:
        assert record["L5_absolute_error"] > .001
        assert record["L14_absolute_error"] > .001


@pytest.mark.parametrize("candidate", ["cloud_mean", "conditional_rms", "global_rms", "law_median"])
def test_magnitude_candidates_are_reciprocal_on_the_unique_median_case(candidate):
    x, y, p = [1, 3, 8], [1, 2, 2], [.2, .5, .3]
    mass = magnitude_readout(candidate, x, y, p)
    inverse = magnitude_readout(candidate, y, x, p)
    assert mass * inverse == pytest.approx(1, abs=1e-13)


def test_cloud_mean_composes_on_every_catalogue_common_channel_law():
    for record in counterexamples():
        x1, x2 = np.array(record["first_force"]), np.array(record["second_force"])
        y, p = np.array(record["common_acceleration"]), record["probabilities"]
        read = vector_readout if record["vector_atoms"] else magnitude_readout
        assert read("cloud_mean", x1 + x2, y, p) == pytest.approx(
            read("cloud_mean", x1, y, p) + read("cloud_mean", x2, y, p), abs=1e-13)
        assert 1 / read("cloud_mean", y, x1 + x2, p) == pytest.approx(
            1 / read("cloud_mean", y, x1, p) + 1 / read("cloud_mean", y, x2, p), abs=1e-13)


def _radial_law(record):
    force = np.linalg.norm(record["force"], axis=1)
    acceleration = np.linalg.norm(record["acceleration"], axis=1)
    result = {}
    for f, a, p in zip(force, acceleration, record["probabilities"]):
        result[(f, a)] = result.get((f, a), 0) + p
    return result


def test_conditional_mean_direction_counterexample_preserves_the_radial_law():
    record = direction_counterexample()
    assert _radial_law(record["first"]) == _radial_law(record["second"])
    assert record["first_readout"] == pytest.approx(1.5)
    assert record["second_readout"] == pytest.approx(2)


def test_balanced_mean_vectors_are_explicitly_undefined_without_a_fallback():
    record = cancellation_counterexample()
    assert record["status"] == "undefined" and record["value"] is None
    assert record["cloud_mean"] == pytest.approx(2)


def test_nonunique_discrete_median_does_not_silently_choose_a_branch():
    assert magnitude_readout("law_median", [1, 2], [1, 1], [.5, .5]) is None


def test_zero_weight_atoms_do_not_change_the_readout():
    for candidate in ("cloud_mean", "conditional_rms", "global_rms", "law_median"):
        assert magnitude_readout(candidate, [1, 2, 100], [1, 1, 1], [.7, .3, 0]) == pytest.approx(
            magnitude_readout(candidate, [1, 2], [1, 1], [.7, .3]))
    assert vector_readout("mean_vector", [[1], [2], [-100]], [[1], [1], [-1]], [.7, .3, 0]) == pytest.approx(
        vector_readout("mean_vector", [[1], [2]], [[1], [1]], [.7, .3]))


def test_counterexample_records_are_json_serializable_and_declare_their_scope():
    records = counterexamples() + [direction_counterexample(), cancellation_counterexample()]
    json.dumps(records, allow_nan=False)
    assert all("scope" in record for record in records)


@pytest.mark.parametrize("x,y,p", [([0, 1], [1, 1], None), ([1, 2], [1, 0], None),
                                  ([1, 2], [1], None), ([1, 2], [1, 1], [-1, 2]),
                                  ([1, 2], [1, 1], [0, 0])])
def test_invalid_finite_law_inputs_are_rejected(x, y, p):
    with pytest.raises(ValueError):
        magnitude_readout("cloud_mean", x, y, p)


def test_vector_atoms_must_be_positive_codirectional_pairs():
    with pytest.raises(ValueError, match="codirectional"):
        vector_readout("mean_vector", [[1, 0]], [[0, 1]], [1])
