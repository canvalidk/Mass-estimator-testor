"""B10: direct ours30 checks for the four short proofs, independent of legacy readouts."""
import math

import numpy as np
import pytest

from met.full import Reading, _h_terms, _log_law, isotropic, ours30, reading_state


def _correlated_covariance():
    return np.array([[2.0, 0.3, 0.2, -0.1],
                     [0.3, 1.0, 0.1, 0.2],
                     [0.2, 0.1, 1.5, 0.4],
                     [-0.1, 0.2, 0.4, 1.2]])


@pytest.mark.parametrize("calibrated", [False, True], ids=["physical", "fixed_centre"])
def test_ours30_physical_and_calibrated_swaps_transform_the_median_and_law(calibrated):
    """Channel exchange and calibrated exchange transform data and covariance together."""
    centre = 2.3
    readings = [Reading(f, a, _correlated_covariance()) for f, a in [
        ([3.0, 0.2], [1.0, -0.4]),
        ([-3.0, 6.0], [3.9, -0.3]),
        ([0.4, 1.6], [0.2, 0.9]),
    ]]
    scale = centre if calibrated else 1.0
    transform = np.block([[np.zeros((2, 2)), scale * np.eye(2)],
                          [np.eye(2) / scale, np.zeros((2, 2))]])
    swapped = []
    for reading in readings:
        data = transform @ np.r_[reading.force, reading.acceleration]
        swapped.append(Reading(data[:2], data[2:], transform @ reading.cov @ transform.T,
                               reading.subspace))
    other_centre = centre if calibrated else 1 / centre
    numerator = centre**2 if calibrated else 1.0
    phis = np.array([-1.1, -0.2, 0.25, 1.2])
    for before, after in zip(readings, swapped):
        np.testing.assert_allclose(
            _h_terms(after.pre(), other_centre, math.pi / 2 - phis),
            _h_terms(before.pre(), centre, phis), rtol=2e-13, atol=1e-13)
        a, b = reading_state(before, centre), reading_state(after, other_centre)
        assert b["Q_plus"] == pytest.approx(a["Q_plus"], rel=1e-11, abs=1e-12)
        assert b["kappa"] == a["kappa"]

    before = ours30(readings, centre, readout=False, grid=(10.0, 2001))
    after = ours30(swapped, other_centre, readout=False, grid=(10.0, 2001))
    assert np.min(before["g"]) < 1  # the misfit factor is exercised, not disabled
    np.testing.assert_allclose(after["g"], before["g"], rtol=1e-11, atol=1e-12)
    assert after["eta"] == pytest.approx(before["eta"], rel=1e-10, abs=1e-12)
    np.testing.assert_allclose(after["lm"], math.log(numerator) - before["lm"][::-1],
                               rtol=0, atol=1e-14)
    np.testing.assert_allclose(after["law"], before["law"][::-1], rtol=1e-10, atol=1e-13)
    assert after["point"] == pytest.approx(numerator / before["point"], rel=1e-7)
    np.testing.assert_allclose(after["interval"], numerator / np.array(before["interval"])[::-1],
                               rtol=1e-7)


def test_ours30_union_recomputes_eta_in_two_reading_example():
    """Raw factors multiply; the common eta cannot be recovered by multiplying private laws."""
    readings = [isotropic([4.0], [1.0], 1.0, 1.0),
                isotropic([1.0], [4.0], 1.0, 1.0)]
    grams = [np.outer([4.0, 1.0], [4.0, 1.0]),
             np.outer([1.0, 4.0], [1.0, 4.0])]
    np.testing.assert_array_equal(sum(grams), [[17.0, 8.0], [8.0, 17.0]])
    union = ours30(readings, 1.0, readout=False)
    private = [ours30([r], 1.0, readout=False) for r in readings]
    assert union["g"] == [1.0, 1.0]
    assert [r["eta"] for r in private] == [1.0, 1.0]
    assert union["phi_hat"] == pytest.approx(math.pi / 4, abs=1e-12)
    assert union["curv"] == pytest.approx(16.0, abs=1e-11)
    assert union["J"] == pytest.approx(225 / 2, abs=1e-11)
    assert union["eta"] == pytest.approx(32 / 225, abs=1e-13)

    masses = np.array([1.0, 3.0])
    direction = np.array([masses, np.ones(2)]) / np.sqrt(1 + masses**2)
    summed_half_h = sum(np.einsum("in,ij,jn->n", direction, gram, direction)
                        for gram in grams) / 2
    np.testing.assert_allclose(summed_half_h, [12.5, 10.9], atol=1e-13)
    for reading, gram in zip(readings, grams):
        np.testing.assert_allclose(_h_terms(reading.pre(), 1.0, np.arctan(masses)),
                                   np.einsum("in,ij,jn->n", direction, gram, direction),
                                   atol=1e-13)
    lm = np.log(masses)
    # In log-mass density, divide by m to obtain the physical-mass density ratio.
    union_log_density = _log_law(readings, 1.0, union["g"], union["eta"], lm) - lm
    union_ratio = math.exp(union_log_density[0] - union_log_density[1])
    assert union_ratio == pytest.approx(5 * math.exp((32 / 225) * (8 / 5)), rel=1e-12)
    # Even after correcting for the duplicated prior, private tempering gives exp(sum ell).
    corrected_private_ratio = 5 * math.exp(8 / 5)
    private_log_densities = [_log_law([r], 1.0, p["g"], p["eta"], lm) - lm
                             for r, p in zip(readings, private)]
    corrected_log_density = sum(private_log_densities) + np.log1p(masses**2)
    assert math.exp(corrected_log_density[0] - corrected_log_density[1]) == pytest.approx(
        corrected_private_ratio, rel=1e-12)
    assert corrected_private_ratio == pytest.approx(24.765162121975575, rel=1e-12)
    assert union_ratio == pytest.approx(6.2776359466, abs=1e-10)
    assert corrected_private_ratio / union_ratio > 3


def test_discarding_rotated_covariances_changes_the_mass_factor():
    """A full-covariance rotation preserves H; dropping induced correlations changes its shape."""
    root_half = 1 / math.sqrt(2)
    rotation = root_half * np.array([[1.0, -1.0], [1.0, 1.0]])
    transform = np.kron(np.eye(2), rotation)
    reading = Reading([1.0, 0.0], [0.0, 1.0], np.diag([1.0, 4.0, 1.0, 1.0]))
    data = transform @ np.r_[reading.force, reading.acceleration]
    rotated_cov = transform @ reading.cov @ transform.T
    rotated = Reading(data[:2], data[2:], rotated_cov)
    discarded = Reading(data[:2], data[2:], np.diag(np.diag(rotated_cov)))
    np.testing.assert_allclose(rotated_cov[:2, :2], [[2.5, -1.5], [-1.5, 2.5]])
    masses = np.array([0.0, 1.0, 3.0])
    phis = np.arctan(masses)
    exact = masses**2 / (masses**2 + 1) + 4 / (masses**2 + 4)
    lost = ((4 / 25) * masses**2 + 1) / ((2 / 5) * masses**2 + 1)
    for full in [reading, rotated]:
        np.testing.assert_allclose(_h_terms(full.pre(), 1.0, phis), exact, atol=1e-13)
    np.testing.assert_allclose(_h_terms(discarded.pre(), 1.0, phis), lost, atol=1e-13)
    assert exact[0] == lost[0] == 1
    assert exact[1] == pytest.approx(13 / 10)
    assert lost[1] == pytest.approx(29 / 35)
    # These H differences are not a candidate-independent multiplicative constant.
    full, turned, diagonal = [ours30([r], 1.0, readout=False) for r in
                              [reading, rotated, discarded]]
    assert turned["point"] == pytest.approx(full["point"], rel=1e-11)
    np.testing.assert_allclose(turned["interval"], full["interval"], rtol=1e-11)
    assert abs(diagonal["point"] - full["point"]) > 0.1


@pytest.mark.parametrize("force,acceleration", [
    ([3.0, 0.0], [0.0, 1.5]),
    ([3.0, 0.0], [-1.5, 0.0]),
    ([3.0, 0.0], [0.0, 0.0]),
    ([0.0, 0.0], [1.5, 0.0]),
    ([0.0, 0.0], [0.0, 0.0]),
], ids=["perpendicular", "anti_parallel", "force_only", "acceleration_only", "all_zero"])
def test_ours30_projection_bound_and_totality(force, acceleration):
    reading = Reading(force, acceleration, _correlated_covariance())
    state = reading_state(reading, 2.3)
    phis = np.linspace(-math.pi / 2, math.pi / 2, 129)
    h = _h_terms(reading.pre(), 2.3, phis)
    assert h.min() >= -1e-12
    assert h.max() <= state["E"] + 1e-12
    result = ours30([reading], 2.3, readout=False)
    assert 0 < result["interval"][0] < result["point"] < result["interval"][1] < np.inf
    assert np.all(np.isfinite(result["law"])) and np.all(result["law"] >= 0)
    if state["E"] == 0:
        assert np.all(h == 0)
        assert result["point"] == pytest.approx(2.3, rel=1e-8)
