"""Independent geometric oracles for the one-reading misalignment investigation.

These checks establish what the diagnostics compute. They do not assert that a
particular mass law is the physically required response to incompatible data.
"""

import math

import numpy as np
import pytest
from scipy.special import i0e, i1e

from met.full import Reading, _stacked_mean_length, isotropic
from met.misalignment import isotropic_minima, residual_profile


def _correlated_reading(d=2):
    rng = np.random.default_rng(671 + d)
    factor = rng.normal(size=(2 * d, 2 * d))
    cov = factor @ factor.T + 0.4 * np.eye(2 * d)
    return Reading(rng.normal(size=d) * 2, rng.normal(size=d), cov)


@pytest.mark.parametrize("d", [2, 3])
def test_profile_matches_independent_whitened_least_squares(d):
    """Fit latent vectors directly, independently of a residual-covariance formula."""
    reading = _correlated_reading(d)
    masses = np.array([0.03, 0.4, 2.0, 30.0])
    got = residual_profile(reading, masses)
    observed = np.r_[reading.force, reading.acceleration]
    chol = np.linalg.cholesky(reading.cov)
    whitened_observed = np.linalg.solve(chol, observed)

    for j, mass in enumerate(masses):
        design = np.vstack([mass * np.eye(d), np.eye(d)])
        whitened_design = np.linalg.solve(chol, design)
        alpha, _, _, _ = np.linalg.lstsq(whitened_design, whitened_observed, rcond=None)
        candidate = design @ alpha
        corrections = observed - candidate
        whitened_corrections = np.linalg.solve(chol, corrections)
        residual_map = np.hstack([np.eye(d), -mass * np.eye(d)])

        np.testing.assert_allclose(got["candidate_force"][j], candidate[:d], rtol=2e-12, atol=2e-12)
        np.testing.assert_allclose(got["candidate_acceleration"][j], candidate[d:], rtol=2e-12, atol=2e-12)
        np.testing.assert_allclose(got["force_correction"][j], corrections[:d], rtol=2e-12, atol=2e-12)
        np.testing.assert_allclose(got["acceleration_correction"][j], corrections[d:], rtol=2e-12, atol=2e-12)
        np.testing.assert_allclose(got["residual"][j], residual_map @ observed, rtol=2e-12, atol=2e-12)
        np.testing.assert_allclose(
            got["residual_covariance"][j], residual_map @ reading.cov @ residual_map.T,
            rtol=2e-12, atol=2e-12,
        )
        assert got["correction_cost"][j] == pytest.approx(
            whitened_corrections @ whitened_corrections, rel=2e-12, abs=2e-12,
        )


def test_profile_cost_and_candidates_are_covariant_under_swap_and_units():
    reading = _correlated_reading()
    masses = np.array([0.05, 0.7, 3.0, 60.0])
    base = residual_profile(reading, masses)
    swap = np.block([[np.zeros((2, 2)), np.eye(2)], [np.eye(2), np.zeros((2, 2))]])
    swapped = Reading(reading.acceleration, reading.force, swap @ reading.cov @ swap.T)
    reverse = residual_profile(swapped, 1 / masses)
    np.testing.assert_allclose(reverse["correction_cost"], base["correction_cost"], rtol=2e-12)
    np.testing.assert_allclose(reverse["candidate_force"], base["candidate_acceleration"], rtol=2e-12, atol=2e-12)
    np.testing.assert_allclose(reverse["candidate_acceleration"], base["candidate_force"], rtol=2e-12, atol=2e-12)

    scale = 7.0
    transform = np.diag([scale, scale, 1.0, 1.0])
    converted = Reading(scale * reading.force, reading.acceleration, transform @ reading.cov @ transform.T)
    units = residual_profile(converted, scale * masses)
    np.testing.assert_allclose(units["correction_cost"], base["correction_cost"], rtol=2e-12)
    np.testing.assert_allclose(units["candidate_force"], scale * base["candidate_force"], rtol=2e-12, atol=2e-12)
    np.testing.assert_allclose(units["candidate_acceleration"], base["candidate_acceleration"], rtol=2e-12, atol=2e-12)


@pytest.mark.parametrize("force_sd,acceleration_sd", [(0.2, 0.2), (3.0, 0.4)])
@pytest.mark.parametrize("cosine", [1.0, math.sqrt(3) / 2])
def test_acute_minimum_matches_closed_form_with_unequal_instruments(force_sd, acceleration_sd, cosine):
    """Fixed whitened lengths 2 and 1 make the eigenvalue oracle independent of SDs."""
    force = force_sd * np.array([2.0, 0.0])
    acceleration = acceleration_sd * np.array([cosine, math.sqrt(max(0, 1 - cosine ** 2))])
    reading = isotropic(force, acceleration, force_sd, acceleration_sd)
    gap = math.sqrt(9 + 16 * cosine ** 2)
    expected_mass = (force_sd / acceleration_sd) * (3 + gap) / (4 * cosine)
    expected_cost = (5 - gap) / 2
    got = isotropic_minima(reading)
    assert got["projective_min_cost"] == pytest.approx(expected_cost, abs=2e-12)
    assert got["positive_min_cost"] == pytest.approx(expected_cost, abs=2e-12)
    assert got["positive_min_location"] == "interior"
    assert got["positive_min_mass"] == pytest.approx(expected_mass, rel=2e-12)
    assert residual_profile(reading, [expected_mass])["correction_cost"][0] == pytest.approx(expected_cost, abs=2e-12)


@pytest.mark.parametrize("sigma", [1.0, 0.1])
@pytest.mark.parametrize("acceleration,projective_cost", [([0.0, 1.0], 1.0), ([-1.0, 0.0], 0.0)])
def test_infinite_mass_infimum_keeps_positive_and_projective_domains_distinct(sigma, acceleration, projective_cost):
    reading = isotropic([2.0, 0.0], acceleration, sigma, sigma)
    got = isotropic_minima(reading)
    assert got["projective_min_cost"] == pytest.approx(projective_cost / sigma ** 2, abs=2e-10)
    assert got["positive_min_cost"] == pytest.approx(1 / sigma ** 2, rel=2e-12)
    assert got["positive_min_location"] == "mass_to_infinity"
    assert got["positive_min_mass"] is None
    # The positive-domain infimum is not attained by any of these finite candidates.
    costs = residual_profile(reading, [0.1, 2.0, 100.0])["correction_cost"]
    assert np.all(costs > got["positive_min_cost"])


def test_zero_endpoint_and_degenerate_minima_are_reported_explicitly():
    acceleration_only = isotropic([0.0, 0.0], [1.0, 0.0], 1.0, 1.0)
    zero_end = isotropic_minima(acceleration_only)
    assert zero_end["positive_min_location"] == "mass_to_zero"
    assert zero_end["positive_min_mass"] is None
    assert zero_end["positive_min_cost"] == pytest.approx(0.0, abs=2e-12)

    for reading, expected_cost in [
        (isotropic([0.0, 0.0], [0.0, 0.0], 1.0, 1.0), 0.0),
        (isotropic([2.0, 0.0], [0.0, 1.0], 2.0, 1.0), 1.0),
    ]:
        got = isotropic_minima(reading)
        assert got["positive_min_location"] == "all_masses"
        assert got["positive_min_mass"] is None
        assert got["positive_min_cost"] == pytest.approx(expected_cost, abs=2e-12)
        np.testing.assert_allclose(residual_profile(reading, [0.01, 1.0, 100.0])["correction_cost"], expected_cost, atol=2e-12)

    # Equal channel energies with opposition have two limiting minima, not a flat cost.
    opposed = isotropic([1.0, 0.0], [-1.0, 0.0], 1.0, 1.0)
    tied = isotropic_minima(opposed)
    assert tied["positive_min_location"] == "both_boundaries"
    assert tied["positive_min_mass"] is None
    assert tied["positive_min_cost"] == pytest.approx(1.0)
    assert tied["projective_min_cost"] == pytest.approx(0.0, abs=2e-12)
    assert residual_profile(opposed, [1.0])["correction_cost"][0] > tied["positive_min_cost"]


@pytest.mark.parametrize("masses", [[0.0], [-1.0], [np.nan], [np.inf], [-np.inf]])
def test_profile_rejects_nonpositive_or_nonfinite_masses(masses):
    with pytest.raises(ValueError):
        residual_profile(isotropic([2.0, 0.0], [0.0, 1.0], 1.0, 1.0), masses)


def test_unsupported_declarations_are_rejected_instead_of_silently_ignored():
    restricted = Reading([2.0, 0.0], [0.0, 1.0], np.eye(4), [[1.0], [0.0]])
    for function in (lambda r: residual_profile(r, [1.0]), isotropic_minima):
        with pytest.raises(ValueError):
            function(restricted)

    anisotropic = Reading([2.0, 0.0], [0.0, 1.0], np.diag([1.0, 2.0, 1.0, 1.0]))
    for reading in (anisotropic, _correlated_reading()):
        with pytest.raises(ValueError):
            isotropic_minima(reading)


@pytest.mark.parametrize("angle", [0, 30, 90, 180])
def test_readout_size_matches_independent_rice_mean(angle):
    """One 2D cloud admits a Bessel-function oracle for the diagnostic A(m)."""
    radians = math.radians(angle)
    force = np.array([2.0, 0.0])
    acceleration = np.array([math.cos(radians), math.sin(radians)])
    sigma = 0.1
    masses = np.unique(np.r_[0.5, np.exp(np.linspace(-5.0, 5.0, 21))])
    mean = (masses[:, None] * force + acceleration) / (1 + masses ** 2)[:, None]
    tau = sigma / np.sqrt(1 + masses ** 2)
    z = np.einsum("ij,ij->i", mean, mean) / (4 * tau ** 2)
    expected = tau * math.sqrt(math.pi / 2) * ((1 + 2 * z) * i0e(z) + 2 * z * i1e(z))
    got = _stacked_mean_length(masses, [isotropic(force, acceleration, sigma, sigma)])
    np.testing.assert_allclose(got, expected, rtol=2e-12, atol=0)
