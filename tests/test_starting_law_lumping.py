"""B26 checks for independent half-Cauchy composition and prior replacement."""

import math

import numpy as np
import pytest
from scipy import integrate, special

from met.composition_laws import LogMassLaw
from met.full import _log_prior
from tools.b26_starting_law_lumping import lumped_log_density, lumped_quantiles, reweight


PROBABILITIES = np.array([.025, .5, .975])


@pytest.mark.parametrize("operation", ["L5", "L14"])
def test_lumped_density_matches_independent_convolution(operation):
    centre = 2.3
    for mass in centre * np.array([.02, .6, 2.7, 31.]):
        # For the parallel sum, the reciprocals of the original variables
        # are half-Cauchy with scale 1/c. Take the reciprocal after summing.
        sum_argument = mass if operation == "L5" else 1 / mass
        sum_scale = centre if operation == "L5" else 1 / centre

        def half_cauchy(value):
            return 2 * sum_scale / (math.pi * (sum_scale**2 + value**2))

        density, _ = integrate.quad(
            lambda part: half_cauchy(part) * half_cauchy(sum_argument - part),
            0., sum_argument, epsabs=1e-13, epsrel=1e-12,
        )
        # The absolute Jacobian from log mass to the sum argument is that
        # argument for either operation.
        actual = math.exp(lumped_log_density(math.log(mass), centre, operation))
        assert actual == pytest.approx(sum_argument * density, rel=2e-11)


@pytest.mark.parametrize("operation", ["L5", "L14"])
def test_lumped_quantiles_match_log_density_quadrature(operation):
    centre = 1.8
    quantiles = lumped_quantiles(centre, operation)

    def log_mass_density(lm):
        return math.exp(lumped_log_density(lm, centre, operation))

    probabilities = [integrate.quad(log_mass_density, -40., math.log(mass),
                                    epsabs=1e-12, epsrel=1e-12)[0] for mass in quantiles]
    np.testing.assert_allclose(probabilities, PROBABILITIES, rtol=0, atol=2e-10)


def test_scaling_and_reciprocal_composition_identities():
    centre, scale = 2.2, 3.7
    lm = np.linspace(-8., 8., 49)
    summands = lumped_log_density(lm, centre, "L5")
    reciprocal = lumped_log_density(-lm, 1 / centre, "L14")
    np.testing.assert_allclose(reciprocal, summands, atol=2e-13, rtol=0)
    for operation in ("L5", "L14"):
        scaled = lumped_log_density(lm + math.log(scale), scale * centre, operation)
        np.testing.assert_allclose(scaled, lumped_log_density(lm, centre, operation),
                                   atol=2e-13, rtol=0)


@pytest.mark.parametrize("operation", ["L5", "L14"])
def test_lumped_half_cauchy_fails_reciprocal_symmetry_at_its_median(operation):
    lower, median, upper = lumped_quantiles(1., operation)
    # A continuous, strictly increasing CDF fixes the candidate centre at
    # its median. L4 would then force opposite quantiles to multiply to c².
    assert abs(math.log(lower * upper / median**2)) > .9


@pytest.mark.parametrize("operation", ["L5", "L14"])
def test_reweight_replaces_pure_prior_on_the_original_grid(operation):
    centre = 2.
    lm = math.log(centre) + np.linspace(-20., 20., 4001)
    original = LogMassLaw(lm, _log_prior(lm, centre))
    result = reweight(original, centre, operation)
    np.testing.assert_array_equal(result.lm, lm)
    log_quantiles = result.logquantile(np.log(PROBABILITIES), np.log1p(-PROBABILITIES))
    np.testing.assert_allclose(log_quantiles, np.log(lumped_quantiles(centre, operation)),
                               atol=2e-5, rtol=0)


def test_reweight_keeps_ours33_raw_law_then_applies_the_same_warp():
    centre, eta = 2., .36
    lm = math.log(centre) + np.linspace(-20., 20., 4001)
    # A bounded, nonconstant data factor, with the same raw/quoted setup as
    # an ours33 law. Its log factor tends to 0 on the left and 2 on the right.
    evidence = 2 * special.expit(lm - math.log(centre))
    original = LogMassLaw(lm, _log_prior(lm, centre) + evidence, eta=eta, warp=True,
                          left_log_constant=-math.log(centre),
                          right_log_constant=2 + math.log(centre))
    result = reweight(original, centre, "L5")
    assert result.warp and result.eta == eta
    expected = evidence + lumped_log_density(lm, centre, "L5")
    np.testing.assert_allclose(result.log_density - result.log_density[0], expected - expected[0],
                               atol=2e-13, rtol=0)

    raw_rank = math.exp(result.logcdf(math.log(centre), quoted=False)[0])
    quoted_rank = math.exp(result.logcdf(math.log(centre))[0])
    assert quoted_rank == pytest.approx(special.ndtr(math.sqrt(eta) * special.ndtri(raw_rank)),
                                        abs=2e-14)
    assert quoted_rank - raw_rank > .03
