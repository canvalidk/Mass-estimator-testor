"""T47: independent probability, density, tail and concentration checks."""
import math

import numpy as np
import pytest
from scipy import integrate, special, stats

from met.cdf_calibration import (normal_cdf_warp, normal_log_cdf_warp,
                                 raw_quantile_log_levels)
from met.cut_consistency import (heavy_tail_limit_quantile, opposed_pair_angle_cdf,
                                  opposed_pair_eta)
from met.residual_mixture import pair_angle_cdf as student_pair_cdf, student_pair_parameters


def test_endpoints_identity_and_exact_median():
    ranks = np.array([0., .03, .5, .97, 1.])
    np.testing.assert_array_equal(normal_cdf_warp(ranks, 1.), ranks)
    for eta in (.01, .3, 1.):
        result = normal_cdf_warp(ranks, eta)
        assert result[0] == 0. and result[-1] == 1. and result[2] == .5
        assert np.all(np.diff(result) > 0)
        levels = raw_quantile_log_levels(.5, eta)
        assert levels["log_cdf"] == levels["log_sf"] == -math.log(2)


def test_semigroup_symmetry_and_outward_quantiles():
    ranks = np.array([.001, .03, .2, .5, .8, .97, .999])
    first = normal_cdf_warp(ranks, .2)
    np.testing.assert_allclose(normal_cdf_warp(first, .7), normal_cdf_warp(ranks, .14), atol=3e-15)
    np.testing.assert_allclose(normal_cdf_warp(1-ranks, .2), 1-first, atol=4e-15)
    assert np.all(first[:3] > ranks[:3]) and np.all(first[4:] < ranks[4:])
    levels = raw_quantile_log_levels(np.array([.025, .975]), .2)
    assert math.exp(levels["log_cdf"][0]) < .025
    assert math.exp(levels["log_sf"][1]) < .025


def test_gaussian_target_against_independent_normal_distribution():
    location, scale, eta = -.3, 1.7, .18
    x = np.linspace(-12., 12., 71)
    raw = stats.norm.cdf(x, loc=location, scale=scale)
    actual = normal_cdf_warp(raw, eta)
    expected = stats.norm.cdf(x, loc=location, scale=scale/math.sqrt(eta))
    np.testing.assert_allclose(actual, expected, atol=3e-9)


@pytest.mark.parametrize("eta", [.5, .8])
def test_non_gaussian_exponential_density_integrates_to_the_warped_cdf(eta):
    # Raw Exp(1). Independently differentiate H(F) and integrate the density,
    # retaining its upper tail through log survival rather than rounded F.
    def density(x):
        log_cdf = math.log(-math.expm1(-x))
        z = special.ndtri_exp(log_cdf) if x < math.log(2) else -special.ndtri_exp(-x)
        return math.exp(.5*math.log(eta)+.5*(1-eta)*z*z-x)

    low, error1 = integrate.quad(density, 0., 1., epsabs=2e-10, limit=200)
    high, error2 = integrate.quad(density, 1., np.inf, epsabs=2e-10, limit=200)
    assert low+high == pytest.approx(1., abs=2e-8)
    assert low == pytest.approx(normal_cdf_warp(-math.expm1(-1.), eta), abs=2e-8)
    assert error1+error2 < 1e-7


def test_units_increasing_coordinates_and_reciprocal_commute_with_the_law():
    raw = stats.lognorm(s=.8, scale=math.exp(.3))
    eta, masses = .3, np.array([.1, .5, 2., 8.])
    base = normal_cdf_warp(raw.cdf(masses), eta)
    changed_units = stats.lognorm(s=.8, scale=7*math.exp(.3))
    np.testing.assert_allclose(normal_cdf_warp(changed_units.cdf(7*masses), eta), base)
    # Re-expressing the same law in log mass uses identical CDF ranks.
    np.testing.assert_allclose(normal_cdf_warp(stats.norm.cdf(np.log(masses), .3, .8), eta), base)
    inverse = stats.lognorm(s=.8, scale=math.exp(-.3))
    np.testing.assert_allclose(normal_cdf_warp(inverse.cdf(1/masses), eta), 1-base, atol=2e-15)


def test_atomic_mass_laws_remain_proper_and_preserve_their_median():
    weights = np.array([.2, .6, .2])  # masses 1, 2, 4; inverse order reverses
    warped = np.diff(np.r_[0., normal_cdf_warp(np.cumsum(weights), .25)])
    assert np.all(warped > 0) and sum(warped) == pytest.approx(1.)
    assert warped[0] < .5 < sum(warped[:2])
    reversed_warp = np.diff(np.r_[0., normal_cdf_warp(np.cumsum(weights[::-1]), .25)])
    np.testing.assert_allclose(reversed_warp, warped[::-1])


def test_log_probability_semigroup_retains_extreme_and_near_one_tails():
    logs = np.array([-math.inf, -1e6, -1000., -10., -math.log(2), -1e-20, 0.])
    first = normal_log_cdf_warp(logs, .2)
    actual = normal_log_cdf_warp(first, .3)
    expected = normal_log_cdf_warp(logs, .06)
    # The extreme inverse-log-normal branch loses about 1e-12 relative
    # precision in log probabilities of magnitude 1e6 on the pinned SciPy.
    np.testing.assert_allclose(actual, expected, rtol=3e-12, atol=1e-13)
    assert np.all(np.isfinite(first[1:])) and first[0] == -math.inf
    np.testing.assert_array_equal(normal_log_cdf_warp(logs, 1.), logs)


def test_small_eta_quantiles_require_log_ranks_not_rounded_raw_cdf_levels():
    levels = raw_quantile_log_levels(.975, 1/90)
    assert levels["log_sf"] < -170.
    assert math.exp(levels["log_cdf"]) == 1.  # ordinary rank loses the tail
    assert normal_log_cdf_warp(levels["log_sf"], 1/90) == pytest.approx(math.log(.025), abs=2e-14)
    raw_gaussian_quantile = -special.ndtri_exp(levels["log_sf"])
    assert raw_gaussian_quantile == pytest.approx(stats.norm.ppf(.975)*math.sqrt(90), rel=2e-15)


def test_subnormal_log_survival_is_not_lost_before_a_second_warp():
    first = normal_log_cdf_warp(-1e-315, .99)
    assert first < 0. and first == pytest.approx(-1.356411194396e-312, rel=2e-10, abs=0.)
    composed = normal_log_cdf_warp(first, .0001)
    direct = normal_log_cdf_warp(-1e-315, .000099)
    assert composed == pytest.approx(direct, abs=2e-12)
    assert normal_cdf_warp(np.nextafter(0., 1.), .99) > 0.
    # Both returned fields retain representable subnormal tails.
    eta = (special.ndtri(.975)/38.)**2
    ranks = raw_quantile_log_levels(.975, eta)
    assert ranks["log_cdf"] < 0. and ranks["log_cdf"] > -1e-315
    assert ranks["log_cdf"] == pytest.approx(math.log1p(-math.exp(ranks["log_sf"])), rel=1e-10, abs=0.)


def test_warped_opposed_pair_tails_shrink_while_the_old_power_keeps_them():
    theta0, radius_error = math.atan(2), 3.
    outside, powered = [], []
    for radius in (100., 1000., 10000.):
        arguments = dict(push=radius, nu_eff=10., error_size=radius_error)
        eta = opposed_pair_eta(**arguments)
        left = opposed_pair_angle_cdf(theta0-.05, truth_angle=theta0, calibrate=False, **arguments)[0]
        right = opposed_pair_angle_cdf(math.pi/2-theta0-.05,
                                        truth_angle=math.pi/2-theta0, calibrate=False, **arguments)[0]
        outside.append(normal_cdf_warp(left, eta)+normal_cdf_warp(right, eta))
        pleft = opposed_pair_angle_cdf(theta0-.05, truth_angle=theta0, calibrate=True, **arguments)[0]
        pright = opposed_pair_angle_cdf(math.pi/2-theta0-.05,
                                         truth_angle=math.pi/2-theta0, calibrate=True, **arguments)[0]
        powered.append(pleft+pright)
    assert outside[0] > outside[1] > outside[2] > 0
    assert powered[-1] > .5
    assert outside[-1] < powered[-1]


def test_a_raw_median_failure_is_exactly_preserved():
    raw = heavy_tail_limit_quantile(.5, truth_angle=math.atan(2), tail_index=.5)
    for eta in (.01, .3, 1.):
        level = raw_quantile_log_levels(.5, eta)
        corrected = heavy_tail_limit_quantile(math.exp(level["log_cdf"]),
                                               truth_angle=math.atan(2), tail_index=.5)
        assert corrected == raw and math.tan(corrected) == pytest.approx(1.8219669037884)


def test_corrected_student_parent_and_cdf_warp_have_shrinking_pair_tails():
    # A combined compatibility check, not a production adapter or a rule
    # selected by its numbers. Both choices remain conditional in T47/T48.
    theta0, tails = math.atan(2), []
    controls = dict(epsabs=2e-12, epsrel=2e-11, limit=200)
    for radius in (100., 1000., 10000.):
        parameters = dict(radius=radius, nu_eff=10., error_size=3.)
        eta = student_pair_parameters(**parameters)["eta"]
        left = student_pair_cdf(theta0-.05, truth_angle=theta0, calibrate=False,
                                 **parameters, **controls)[0]
        right = student_pair_cdf(math.pi/2-theta0-.05, truth_angle=math.pi/2-theta0,
                                  calibrate=False, **parameters, **controls)[0]
        tails.append(normal_cdf_warp(left, eta)+normal_cdf_warp(right, eta))
    assert tails[0] > tails[1] > tails[2] > 0.


def test_vanishing_eta_can_amplify_exponentially_small_remote_mass():
    # P_N=(1-e^-N) N(0,1)+e^-N N(N^2,1) has TV error <=e^-N.
    # At N^2/2 the raw survival is asymptotically e^-N. After eta=1/N
    # it instead tends to sf(sqrt(2)), and escapes under sqrt(eta) scaling.
    tails = []
    for n in (20., 80., 320., 1280., 5120.):
        cutoff = n*n/2
        log_tail = np.logaddexp(math.log1p(-math.exp(-n))+stats.norm.logsf(cutoff),
                                -n+stats.norm.logsf(cutoff-n*n))
        tails.append(math.exp(normal_log_cdf_warp(log_tail, 1/n)))
    assert tails[-1] == pytest.approx(stats.norm.sf(math.sqrt(2)), abs=3e-4)
    assert tails[-1] > .07
    assert abs(tails[-1]-stats.norm.sf(math.sqrt(2))) < abs(tails[0]-stats.norm.sf(math.sqrt(2)))


@pytest.mark.parametrize("eta", [0., -1., 1.1, math.inf, math.nan])
def test_invalid_eta_is_rejected_without_adding_a_zero_eta_law(eta):
    with pytest.raises(ValueError):
        normal_cdf_warp(.3, eta)


def test_probability_and_log_probability_domains():
    with pytest.raises(ValueError):
        normal_cdf_warp(1.1, .5)
    with pytest.raises(ValueError):
        normal_log_cdf_warp(.1, .5)
    with pytest.raises(ValueError):
        raw_quantile_log_levels(1., .5)
