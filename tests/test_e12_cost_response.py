"""Finite implications of the reciprocal one-form and the two cost-response regimes."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from explore.e12_tangent_pairings import tangent_pairings
from met.full import isotropic


def candidate_cost(readings, angle):
    # Here both errors are one, so the mass-angle coordinate is atan(m).
    return sum(np.sum((row.force * math.cos(angle) - row.acceleration * math.sin(angle)) ** 2)
               for row in readings)


@pytest.mark.parametrize("angles", [(.13, .92), (.97, .23)])
def test_integrated_reciprocal_work_recovers_the_gaussian_contrast(angles):
    start, end = angles
    for observations in [[(4, 1), (1, 4)], [(2, 2), (2, 2)], [(3, -1), (.5, 4)]]:
        readings = [isotropic([f], [a], 1., 1.) for f, a in observations]
        work, _ = quad(lambda x: tangent_pairings(readings, math.tan(x)).score, start, end)
        cost_change = candidate_cost(readings, end) - candidate_cost(readings, start)
        assert work == pytest.approx(-cost_change / 2, abs=1e-11)


def test_rapid_cost_response_collapses_on_the_fixed_conflict():
    angles = np.linspace(0., math.pi / 2, 10001)
    delta = angles - math.pi/4
    k = 9 + 16 * np.sin(delta)**2
    raw = [isotropic([4.], [1.], 1., 1.), isotropic([1.], [4.], 1., 1.)]
    for i in [0, 1249, 4771, 5000, 7672, 9987]:
        assert k[i] == pytest.approx(candidate_cost(raw, angles[i]))
    # f(x)=exp(-sqrt(x)), with unbounded doubling ratio even though it
    # penalizes large costs far less than Gaussian evidence does.
    variances = []
    for scale in [1., 100., 10000.]:
        weights = np.exp(-math.sqrt(scale) * (np.sqrt(k) - 3))
        variances.append(np.average(delta**2, weights=weights))
    assert variances[-1] < .01
    assert np.all(np.diff(variances) < 0.)


def test_bounded_doubling_preserves_mass_away_from_the_true_angle():
    # A scatter equal to the repeated-push population limit times N.
    angle = np.linspace(.001, math.pi / 2 - .001, 2001)
    magnitude, truth = 2., .6
    shape = 1 + magnitude**2 * np.sin(angle-truth)**2
    target = angle > 1.
    limiting = np.mean(shape[target]**-1.5) / np.mean(shape**-1.5) * target.mean()
    for n in [100, 100000]:
        weights = (1 + n * shape)**-1.5
        assert np.max(weights)/np.min(weights) < (1+magnitude**2)**1.5
        assert np.sum(weights[target]) / np.sum(weights) == pytest.approx(limiting, abs=.002)
    assert limiting > .1


def test_a_count_dependent_power_changes_relative_odds_on_zero_padding():
    theta = np.array([.2, .8])
    cost = np.cos(theta)**2 + 4 * np.sin(theta)**2
    ratios = [((1+cost[0])/(1+cost[1]))**-n for n in (2, 3)]
    assert ratios[0] != pytest.approx(ratios[1])
