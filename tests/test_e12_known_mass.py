"""Transport of the E12 witnesses using physical records, covariance and admissible masses."""
import math

import numpy as np
import pytest
from scipy.integrate import cumulative_trapezoid

from explore.e12_known_mass import known_mass_log_evidence
from explore.e12_residual_join import residual_log_evidence


ASYMMETRIC = np.array([[[3., 1.]], [[1., 2.]]])


def translated(matrix, pairs, covariance):
    return pairs @ matrix.T, matrix @ covariance @ matrix.T


@pytest.mark.parametrize("tau,median", [(0, 1.528611), (64, 1.510203)])
def test_direct_inference_on_the_totals_preserves_the_entire_law(tau, median):
    b = 2.
    total, cov = translated(np.array([[1, b], [0, 1]]), ASYMMETRIC, np.eye(2))
    angles = np.linspace(0, math.pi / 2, 60001)
    masses = np.tan(angles)
    base = residual_log_evidence(ASYMMETRIC.transpose(0, 2, 1) @ ASYMMETRIC, angles, tau=tau)
    joint = known_mass_log_evidence(total, cov, b + masses, bounds=(b, np.inf), tau=tau)
    np.testing.assert_allclose(joint, base, atol=5e-12)
    # Uniform base angular measure pulls forward to the translated half-Cauchy.
    cdf = cumulative_trapezoid(np.exp(joint - joint.max()), angles, initial=0)
    assert b + np.tan(np.interp(0.5, cdf / cdf[-1], angles)) == pytest.approx(b + median, abs=1e-6)


@pytest.mark.parametrize("rho", [-0.6, 0.0, 0.75])
@pytest.mark.parametrize("fit_domain", ["physical", "projective"])
def test_correlated_channel_reciprocity_subsequent_known_load_and_frame_rotation(rho, fit_domain):
    pairs = np.array([[[1., 2.], [3., -1.]], [[-1., 0.5], [2., 1.]]])
    covariance = np.array([[4., 2 * rho], [2 * rho, 1.]])
    masses = np.geomspace(0.25, 5, 31)
    kw = dict(tau=100, fit_domain=fit_domain)
    base = known_mass_log_evidence(pairs, covariance, masses, bounds=(0.25, 5), **kw)
    # Compose channel swap, unit conversion, and a known total mass offset.
    b, scale = 1.25, 1.7
    matrix = np.array([[b, scale], [1., 0.]])
    final, cov = translated(matrix, pairs, covariance)
    final[0] = np.array([[0., -1.], [1., 0.]]) @ final[0]
    final[1] = np.diag([-1, 1]) @ final[1]
    law = known_mass_log_evidence(final, cov, b + scale / masses,
                                  bounds=(b + scale / 5, b + scale / 0.25), **kw)
    np.testing.assert_allclose(law, base, atol=2e-9)


def test_only_projective_fitting_commutes_with_forgetting_the_mass_bounds():
    masses = np.linspace(.7, 1.6, 8)
    def contrast(bounds, fit_domain):
        law = known_mass_log_evidence(ASYMMETRIC, np.eye(2), masses, bounds=bounds, fit_domain=fit_domain)
        return law - law[0]
    np.testing.assert_allclose(contrast((.65, 1.7), "projective"),
                               contrast((0., np.inf), "projective"), atol=1e-10)
    assert np.max(abs(contrast((.65, 1.7), "physical") - contrast((0., np.inf), "physical"))) > 0.03


@pytest.mark.parametrize("tau", [0, 64])
def test_support_independent_extension_retains_split_but_loses_information_on_balanced_vector(tau):
    masses = np.array([.5, 1., 3.])
    logp = known_mass_log_evidence(np.eye(2)[None], np.eye(2), masses,
                                   fit_domain="projective", tau=tau)
    np.testing.assert_allclose(logp, 0.0)
    central = known_mass_log_evidence(ASYMMETRIC, np.eye(2), masses, fit_domain="projective", tau=0)
    alternative = known_mass_log_evidence(ASYMMETRIC, np.eye(2), masses, fit_domain="projective", tau=64)
    assert np.ptp(central - alternative) > .01


@pytest.mark.parametrize("rho", [-0.5, 0.5])
def test_fixed_compatible_conflict_and_tie_limits_on_genuine_correlated_instruments(rho):
    cov = np.array([[1., rho], [rho, 1.]])
    masses = np.array([0., 1., 2., np.inf])
    compatible = np.array([[[2., 1.]]])
    aligned = known_mass_log_evidence(compatible * 1e4, cov, masses)
    assert aligned[2] == pytest.approx(0., abs=1e-8)
    assert max(aligned[[0, 1, 3]]) < -1e6
    for pair in [ASYMMETRIC, np.linalg.cholesky(cov).T[None]]:
        reference = known_mass_log_evidence(pair * 1e3, cov, masses, tau=64)
        limit = known_mass_log_evidence(pair * 1e4, cov, masses, tau=0)
        assert .05 < np.ptp(reference) < 100
        np.testing.assert_allclose(reference, limit, atol=3e-5)


def test_obtuse_aperture_can_put_two_score_roots_inside_the_admissible_arc():
    # Physical channel correlation +1/2 gives a whitened aperture 2pi/3.
    pilot = math.pi / 12
    arc = 2 * math.pi / 3
    top = np.array([math.sin(pilot), math.cos(pilot)])
    B = 1000 * (np.eye(2) + 2 * np.outer(top, top))
    def standardized_score(theta):
        u = np.array([math.sin(theta), math.cos(theta)])
        v = np.array([math.cos(theta), -math.sin(theta)])
        return (v @ B @ u) / math.sqrt(u @ B @ u)
    assert 0 < pilot < pilot + math.pi / 2 < arc
    assert standardized_score(pilot) == pytest.approx(0., abs=1e-12)
    assert standardized_score(pilot + math.pi / 2) == pytest.approx(0., abs=1e-12)
    assert abs(standardized_score(pilot + math.pi / 4)) > 1.96
