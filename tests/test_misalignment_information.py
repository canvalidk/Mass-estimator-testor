"""Algebraic diagnostics of positivity, normalization and escaping mass laws.

These establish what the current operation does. They do not choose the physical
response to incompatible measurements. The geometric, angle-integral and symmetry
oracles are independent of the estimator's explained-energy construction.
"""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

from met.full import isotropic, misfit_weight, ours29
from met.misalignment import candidate_mean_length, investigate_cell, isotropic_minima, residual_profile
from met.preset import load


@pytest.mark.parametrize("force,acceleration", [
    ([2.0, 0.0], [1.0, 0.0]),
    ([2.0, 0.0], [math.sqrt(3) / 2, 0.5]),
    ([2.0, 0.0], [0.0, 1.0]),
    ([2.0, 0.0], [-math.sqrt(3) / 2, 0.5]),
    ([2.0, 0.0], [-1.0, 0.0]),
    ([1.0, 0.0], [-2.0, 0.0]),
    ([1.0, 0.0], [-1.0, 0.0]),
    ([1.0, 0.0], [0.0, 1.0]),
])
def test_positive_cost_bounds_against_direct_latent_vector_fits(force, acceleration):
    """Optimize independently over angle; each objective fits the latent vector by LS."""
    f, a = np.array(force), np.array(acceleration)
    sigma = 0.3
    observed = np.r_[f, a] / sigma

    def cost(phi):
        design = np.vstack([math.sin(phi) * np.eye(2), math.cos(phi) * np.eye(2)]) / sigma
        fitted, *_ = np.linalg.lstsq(design, observed, rcond=None)
        error = observed - design @ fitted
        return float(error @ error)

    low = minimize_scalar(cost, bounds=(0, math.pi / 2), method="bounded")
    high = minimize_scalar(lambda phi: -cost(phi), bounds=(0, math.pi / 2), method="bounded")
    lower = min(cost(0), cost(math.pi / 2), low.fun)
    upper = max(cost(0), cost(math.pi / 2), -high.fun)
    bounds = isotropic_minima(isotropic(f, a, sigma, sigma))
    assert bounds["positive_min_cost"] == pytest.approx(lower, abs=1e-8)
    assert bounds["positive_cost_supremum"] == pytest.approx(upper, abs=1e-8)
    assert bounds["positive_cost_span"] == pytest.approx(upper - lower, abs=1e-8)
    assert bounds["positive_min_cost"] == pytest.approx(
        bounds["projective_min_cost"] + bounds["positive_domain_penalty"], abs=1e-10)


@pytest.mark.parametrize("dimension", [1, 2, 3])
def test_positivity_cost_exists_even_without_collinearity_conditions(dimension):
    f, a = np.zeros(dimension), np.zeros(dimension)
    f[0], a[0] = 2.0, -1.0
    bounds = isotropic_minima(isotropic(f, a, 0.1, 0.1))
    assert bounds["projective_min_cost"] == pytest.approx(0, abs=1e-12)
    assert bounds["positive_domain_penalty"] == pytest.approx(100.0)
    assert bounds["positive_min_location"] == "mass_to_infinity"
    # nu counts equalities, not the sign inequality. In 1D the current rule
    # explicitly preserves weight 1, even if a positive-domain cost is supplied.
    if dimension == 1:
        assert misfit_weight(bounds["positive_min_cost"], 0) == 1.0


def test_cost_contrast_survives_a_large_common_baseline():
    reading = isotropic([1e10, 0.0], [-1e-7, 1e10], 1.0, 1.0)
    bounds = isotropic_minima(reading)
    # Each cost rounds to about 1e20; subtracting those float64 values loses
    # the 1000-nat-scale geometric contrast, which the dot product still carries.
    assert bounds["positive_cost_span"] == pytest.approx(1000.0)
    assert bounds["positive_domain_penalty"] == pytest.approx(1000.0)


@pytest.mark.parametrize("cosine", [math.sqrt(3) / 2, 0.0, -1.0])
@pytest.mark.parametrize("centre", [1.0, 3.0])
def test_normalized_law_uses_excess_cost_not_its_constant_baseline(cosine, centre):
    cell = dict(id="c000", case="cost identity", force=[2.0, 0.0],
                acceleration=[cosine, math.sqrt(max(0, 1 - cosine ** 2))],
                noise_scale=0.2, force_sd=0.2, acceleration_sd=0.2)
    row, curves = investigate_cell(cell, centre, dict(profile_log_half_span=12, profile_nodes=101))
    reading = isotropic(cell["force"], cell["acceleration"], 0.2, 0.2)
    np.testing.assert_allclose(curves["correction_cost"],
                               curves["positive_excess_cost"] + row["positive_min_cost"], rtol=2e-15)
    for label, weight in (("raw", 1.0), ("quoted", row["misfit_weight"] * row["bartlett_eta"])):
        lm = curves[f"{label}_lm"]
        mass = np.exp(lm)
        cost = residual_profile(reading, mass)["correction_cost"]
        log_prior = math.log(centre) + lm - np.logaddexp(2 * math.log(centre), 2 * lm)
        log_kernel = log_prior - weight * (cost - row["positive_min_cost"]) / 2
        density = np.exp(log_kernel - log_kernel.max())
        density /= np.trapezoid(density, lm)
        np.testing.assert_allclose(curves[f"{label}_density"], density, rtol=2e-12, atol=2e-14)
        assert row[label]["log_evidence_span"] == pytest.approx(weight * row["positive_cost_span"] / 2)


@pytest.mark.parametrize("cosine", [math.sqrt(3) / 2, 0.0, -math.sqrt(3) / 2])
def test_saturated_misfit_exactly_cancels_common_noise_scaling(cosine):
    """The frozen law follows from current g; this does not establish physical adequacy."""
    a = [cosine, math.sqrt(1 - cosine ** 2)]
    results = [ours29([isotropic([2, 0], a, sigma, sigma)], 1,
                     grid=(12, 1601), readout=False) for sigma in (0.2, 0.05)]
    assert all(result["Q"][0] > 1 for result in results)
    np.testing.assert_allclose(results[0]["law"], results[1]["law"], rtol=3e-12, atol=2e-15)
    assert results[0]["Q"][0] / results[1]["Q"][0] == pytest.approx(1 / 16)


@pytest.mark.parametrize("sigma", [0.5, 0.2, 0.1])
def test_equal_energy_antiparallel_centres_and_intervals_against_angle_integrals(sigma):
    """A stable point and median coexist with mass leaving every fixed finite band.

    At c=s=1 the prior is uniform in angle. The kernel exp(-sin(phi)cos(phi)/sigma²)
    gives an independent normalizer and probabilities at the reported quantiles.
    """
    reading = isotropic([1.0, 0.0], [-1.0, 0.0], sigma, sigma)
    result = ours29([reading], 1.0)
    assert result["readout"] == pytest.approx(1.0, rel=2e-10)
    assert result["interval"][0] * result["interval"][1] == pytest.approx(1.0, rel=2e-6)
    assert result["Q"][0] == pytest.approx(0, abs=2e-12)
    assert result["g"] == [1.0] and result["eta"] == 1.0

    def kernel(phi):
        return math.exp(-math.sin(phi) * math.cos(phi) / sigma ** 2)

    normalizer = quad(kernel, 0, math.pi / 2, epsabs=1e-13, epsrel=1e-12)[0]
    for mass, probability in zip((result["interval"][0], result["median"], result["interval"][1]),
                                 [0.025, 0.5, 0.975]):
        actual = quad(kernel, 0, math.atan(mass), epsabs=1e-13, epsrel=1e-12)[0] / normalizer
        assert actual == pytest.approx(probability, abs=8e-8)
    band = quad(kernel, math.atan(0.5), math.atan(2.0), epsabs=1e-100, epsrel=1e-12)[0] / normalizer
    assert band < {0.5: 0.3, 0.2: 0.001, 0.1: 1e-15}[sigma]

    masses = np.exp(np.linspace(-8, 8, 31))
    size = candidate_mean_length(reading, masses)
    reverse = candidate_mean_length(reading, 1 / masses)
    np.testing.assert_allclose(reverse, masses * size, rtol=3e-14)


# SD 0.1 was a strict xfail (2.341636397, the O13 median defect); fixed 27 Oct by
# met.full._plateau_centre. SD 0.05 and 0.02 sit deeper in the plateau.
@pytest.mark.parametrize("sigma", [0.5, 0.2, 0.1, 0.05, 0.02])
def test_equal_energy_antiparallel_mathematical_median_is_one(sigma):
    result = ours29([isotropic([1.0, 0.0], [-1.0, 0.0], sigma, sigma)], 1.0, readout=False)
    assert result["median"] == pytest.approx(1.0, rel=2e-7)


def test_positive_boundary_preset_keeps_the_swap_pair_and_equal_energy_case():
    spec = load("presets/misalignment_positive_boundaries.yaml")
    assert len(spec["cells"]) == 20
    assert spec["kind"] == "misalignment"
    assert {cell["case"] for cell in spec["cells"]} == {
        "aligned-equal-energy", "anti-parallel-force-dominant",
        "anti-parallel-acceleration-dominant", "anti-parallel-equal-energy"}
    for scale in (0.5, 0.2, 0.1, 0.05, 0.02):
        cells = {cell["case"]: cell for cell in spec["cells"] if cell["noise_scale"] == scale}
        force_dominant = cells["anti-parallel-force-dominant"]
        acceleration_dominant = cells["anti-parallel-acceleration-dominant"]
        np.testing.assert_array_equal(force_dominant["force"], -np.array(acceleration_dominant["acceleration"]))
        np.testing.assert_array_equal(force_dominant["acceleration"], -np.array(acceleration_dominant["force"]))
        assert force_dominant["force_sd"] == acceleration_dominant["acceleration_sd"]
        assert force_dominant["acceleration_sd"] == acceleration_dominant["force_sd"]


@pytest.mark.parametrize("c", [0.3, 3.0])
@pytest.mark.parametrize("sigma", [0.1, 0.05])
def test_plateau_quantiles_commute_with_the_swap(c, sigma):
    """The plateau rule off the symmetric case: swapping force and acceleration (with their
    SDs) and inverting c must invert every quantile, m -> 1/m (T5)."""
    import met.full as full
    fwd = ours29([isotropic([1.0, 0.0], [-1.0, 0.0], sigma, sigma)], c, readout=False)
    bwd = ours29([isotropic([-1.0, 0.0], [1.0, 0.0], sigma, sigma)], 1 / c, readout=False)
    assert fwd["median"] * bwd["median"] == pytest.approx(1, rel=1e-6)
    assert fwd["interval"][0] * bwd["interval"][1] == pytest.approx(1, rel=1e-6)


def test_plateau_rule_does_not_fire_on_an_ordinary_law(monkeypatch):
    import met.full as full
    reading = [isotropic([2.0, 0.3], [1.0, 0.1], 0.2, 0.1)]
    before = ours29(reading, 2.0, readout=False)
    monkeypatch.setattr(full, "PLATEAU_TOL", -1.0)
    after = ours29(reading, 2.0, readout=False)
    assert before["median"] == after["median"] and before["interval"] == after["interval"]
