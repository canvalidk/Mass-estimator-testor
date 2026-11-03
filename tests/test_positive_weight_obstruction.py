"""T46: general identity, independent polar checks and strict weak-signal bias."""
import math

import numpy as np
import pytest
from scipy import special

from met.gaussian_information import gaussian_angle_cdf, gaussian_angle_quantile
from met.positive_weight_obstruction import (null_cross_moment,
                                             polar_weighted_scatter,
                                             weak_signal_certificate)


CONTROLS = dict(epsabs=2e-11, epsrel=2e-10, limit=150)
WEIGHTS = {
    "kink": lambda q: 1. if q <= .5 else .5/q,
    "inverse_gamma": lambda q: 10/(10+q),
    "exponential": lambda q: math.exp(-q),
    "constant": lambda q: 1.,
}


@pytest.mark.parametrize("name", list(WEIGHTS))
def test_general_null_identity_matches_closed_and_independent_polar_integrals(name):
    weight = WEIGHTS[name]
    bias, error = null_cross_moment(weight, zero_weight=1., **CONTROLS)
    expected = {
        "kink": (math.exp(-.5)-.5*special.exp1(.5))/math.pi,
        "inverse_gamma": (1-10*math.exp(10)*special.exp1(10))/math.pi,
        "exponential": 1/(2*math.pi),
        "constant": 0.,
    }[name]
    assert bias == pytest.approx(expected, abs=2e-12)
    assert error < 1e-10
    null = polar_weighted_scatter([0., 0.], weight, zero_weight=1., angular_order=40, **CONTROLS)
    assert null["scatter"][0, 1] == pytest.approx(bias, abs=2e-9)
    assert null["scatter"][0, 0] == pytest.approx(null["scatter"][1, 1], abs=2e-10)
    assert null["radial_error"] < 1e-8


@pytest.mark.parametrize("name", ["kink", "inverse_gamma", "exponential"])
def test_each_nonconstant_monotone_control_has_a_wrong_weak_gaussian_fit(name):
    weight, theta = WEIGHTS[name], math.atan(2)
    bias, _ = null_cross_moment(weight, zero_weight=1., **CONTROLS)
    # Exact shift bound = B/4, strictly less than |cos(2theta)|*B=3B/5.
    radius = math.sqrt(math.log1p((bias/4)**2))
    certificate = weak_signal_certificate(weight, zero_weight=1., angle=theta,
                                          radius=radius, **CONTROLS)
    assert certificate["score_shift_upper"] == pytest.approx(bias/4)
    assert certificate["score_upper"] < 0 and certificate["wrong_interior_fit"]
    u = np.array([math.sin(theta), math.cos(theta)])
    v = np.array([math.cos(theta), -math.sin(theta)])
    matrix = polar_weighted_scatter(radius*u, weight, zero_weight=1., angular_order=40, **CONTROLS)["scatter"]
    score = float(u@matrix@v)
    assert certificate["score_lower"] <= score <= certificate["score_upper"]
    assert score < 0 and matrix[0, 1] > 0
    eigenvalues, vectors = np.linalg.eigh(matrix)
    principal = vectors[:, 1]
    fitted_mass = abs(principal[0]/principal[1])
    assert eigenvalues[1] > eigenvalues[0] and abs(fitted_mass-2) > .5


def test_constant_control_is_not_declared_a_fisher_consistency_failure():
    result = weak_signal_certificate(WEIGHTS["constant"], zero_weight=1.,
                                      angle=math.atan(2), radius=.001, **CONTROLS)
    assert result["bias"] == 0 and not result["wrong_interior_fit"]


def test_reciprocal_symmetry_does_not_remove_the_bias():
    weight = WEIGHTS["inverse_gamma"]
    first = polar_weighted_scatter([.01, .02], weight, zero_weight=1., angular_order=40, **CONTROLS)["scatter"]
    second = polar_weighted_scatter([.02, .01], weight, zero_weight=1., angular_order=40, **CONTROLS)["scatter"]
    np.testing.assert_allclose(second, first[::-1, ::-1], atol=2e-10)
    assert first[0, 1] > 0


def test_uniform_angle_median_is_bracketed_for_every_tested_precision():
    # The theorem proves all nonnegative precisions; this deterministic grid
    # includes zero temperature and narrow laws on both sides of pi/4.
    truth, center = math.atan(2), math.pi/4
    for fitted in (.3, 1.):
        assert max(fitted, center) < truth
        for precision in (0., .1, 1., 10., 100., 1e4):
            median = gaussian_angle_quantile(.5, fitted_angle=fitted, precision=precision)
            assert min(fitted, center)-2e-10 <= median <= max(fitted, center)+2e-10
            lower = gaussian_angle_cdf(min(fitted, center), fitted_angle=fitted, precision=precision)[0]
            upper = gaussian_angle_cdf(max(fitted, center), fitted_angle=fitted, precision=precision)[0]
            assert lower <= .5+2e-10 and upper >= .5-2e-10
            assert median < truth


def test_explicit_bounded_weight_contract_and_controls():
    with pytest.raises(ValueError):
        null_cross_moment(lambda q: 2., zero_weight=1., **CONTROLS)
    with pytest.raises(ValueError):
        polar_weighted_scatter([0., 0.], WEIGHTS["constant"], zero_weight=1.,
                                angular_order=1, **CONTROLS)
