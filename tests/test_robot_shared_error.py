"""Physical excitation and budget sensitivity on the maintained real-data loader."""
import hashlib

import numpy as np
import pytest

from met import robot
from met.robot_shared_error import ATTRIBUTION, _manifest, describe_readings


@pytest.fixture(scope="module")
def baseline():
    calibration = robot.calibrate()
    with_floor, _ = robot.windows("1-baseline", 0.2, calibration, floor=True)
    without_floor, _ = robot.windows("1-baseline", 0.2, calibration, floor=False)
    return calibration, describe_readings(with_floor, calibration.mass), describe_readings(without_floor, calibration.mass)


def test_robot_spectrum_uses_the_physical_common_frame(baseline):
    _, measured, _ = baseline
    acceleration = measured["acceleration"]
    np.testing.assert_allclose(measured["spectrum"]["singular_values"],
                               np.linalg.svd(acceleration, compute_uv=False), rtol=1e-13)
    np.testing.assert_allclose(measured["spectrum"]["centered_singular_values"],
                               np.linalg.svd(acceleration - acceleration.mean(0), compute_uv=False), rtol=1e-13)
    assert measured["readings"] == 31


def test_floor_changes_declared_resolution_without_changing_the_measured_excitation(baseline):
    _, with_floor, without_floor = baseline
    np.testing.assert_array_equal(with_floor["acceleration"], without_floor["acceleration"])
    np.testing.assert_array_equal(with_floor["force"], without_floor["force"])
    np.testing.assert_array_equal(with_floor["spectrum"]["singular_values"], without_floor["spectrum"]["singular_values"])
    assert with_floor["spectrum"]["noise_radius"] > without_floor["spectrum"]["noise_radius"]
    assert with_floor["summary"]["rank_lower"] <= without_floor["summary"]["rank_lower"]


def test_shared_fit_and_heldout_score_are_not_reported_as_pvalues(baseline):
    _, result, _ = baseline
    assert "pvalue" not in result["tilt"] and "pvalue" not in result["even_row_score"]
    assert result["tilt"]["fitted_cost"] <= result["tilt"]["mass_only_cost"] <= result["tilt"]["raw_cost"]
    assert any("noisy" in assumption for assumption in result["assumptions"])


def test_source_manifest_hashes_the_public_data_and_keeps_attribution():
    name = "0-calibration_fts-accel.csv"
    assert _manifest([name])[name] == hashlib.sha256((robot.RAW / name).read_bytes()).hexdigest()
    assert ATTRIBUTION["doi"] == "10.5281/zenodo.11096791"
    assert ATTRIBUTION["creator"] == "Aleksander Skrede" and ATTRIBUTION["license"] == "CC-BY-4.0"


@pytest.mark.xfail(strict=True, reason="Source-count discrepancy: met.robot.calibrate's existing docstring "
                   "says 23 static orientations, while 0-calibration_fts-accel.csv contains 24 data rows "
                   "and the source metadata specifies 24. The loader and earlier record are retained (O16).")
def test_calibration_count_matches_the_existing_23_orientation_description():
    assert len(robot._load("0-calibration_fts-accel")) == 23
