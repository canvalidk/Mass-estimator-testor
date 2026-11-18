"""Checks for the conditional minimax test on the reciprocal tangent contrasts."""
import math

import numpy as np
import pytest
from scipy.special import hyp0f1
from scipy.stats import chi2, ncx2, norm

from explore.e12_tangent_pairings import tangent_pairings
from met.full import isotropic


def test_reciprocal_contrast_reports_conflict_at_mass_one():
    separate = tangent_pairings([isotropic([4], [1], 1, 1), isotropic([1], [4], 1, 1)], 1.)
    joint = tangent_pairings([isotropic([4, 1], [1, 4], 1, 1)], 1.)
    assert separate.between == pytest.approx(9.)
    assert joint.within == pytest.approx(9.)
    assert separate.degrees == (0, 1, 1)
    assert joint.degrees == (1, 0, 1)
    assert chi2.sf(separate.between, 1) == pytest.approx(2 * norm.sf(3.))
    assert separate.common == pytest.approx(0., abs=1e-12)


@pytest.mark.parametrize("rank", [.1, .5, .95])
def test_local_confidence_quantiles_are_uniquely_calibrated_at_each_rank(rank):
    state = tangent_pairings([isotropic([4], [1], 1, 1), isotropic([1], [4], 1, 1)], 1.)
    error_sd = 1 / math.sqrt(state.weight)
    for shift in [-.7, .2, 1.3]:
        # P(shift <= observed correction + offset) in the common tangent experiment.
        offset = norm.ppf(rank) * error_sd
        cutoff = shift - offset
        assert norm.sf(cutoff, loc=shift, scale=error_sd) == pytest.approx(rank)
        assert norm.sf(cutoff - .1 * error_sd, loc=shift, scale=error_sd) > rank


@pytest.mark.parametrize("q", [1, 2, 5])
def test_spherical_mixture_likelihood_ratio_orders_exactly_by_contrast_energy(q):
    size = .05
    delta = 1.4
    energies = np.array([.05, 1., 3., 9., 24.])
    averaged_ratio = math.exp(-delta**2 / 2) * hyp0f1(q/2, delta**2 * energies / 4)
    density_ratio = ncx2.pdf(energies, q, delta**2) / chi2.pdf(energies, q)
    np.testing.assert_allclose(averaged_ratio, density_ratio, rtol=1e-12)
    assert np.all(np.diff(averaged_ratio) > 0.)
    threshold = chi2.isf(size, q)
    assert ncx2.sf(threshold, q, delta**2) > size
    if q == 1:
        power = norm.sf(math.sqrt(threshold) - delta) + norm.cdf(-math.sqrt(threshold) - delta)
        assert ncx2.sf(threshold, q, delta**2) == pytest.approx(power)
