"""T44: component normalization, single-reading limits and the new score bias."""
import math

import numpy as np
import pytest
from scipy import special, stats

from met.cut_consistency import heavy_tail_limit_cdf
from met.cut_normalization import (discrete_mixture_cdf, discrete_mixture_quantile,
                                    gaussian_component, ig_normalized_cut_cdf,
                                    positive_q_bias)


CONTROLS = dict(epsabs=2e-10, epsrel=2e-9, limit=150)


def test_discrete_cut_keeps_q_while_factor_average_reweights_by_evidence():
    center = math.atan(2)
    arguments = dict(observed_angle=center, radius=10., variances=[.01, 100.],
                     weights=[.5, .5], component_temperature=1., **CONTROLS)
    cut = discrete_mixture_cdf(center, kind="normalized_cut", **arguments)
    factor = discrete_mixture_cdf(center, kind="factor_average", **arguments)
    np.testing.assert_allclose(cut["component_evidence"], [.01595849, .89099184], atol=5e-9)
    np.testing.assert_array_equal(cut["effective_weights"], [.5, .5])
    expected = cut["component_evidence"]/sum(cut["component_evidence"])
    np.testing.assert_allclose(factor["effective_weights"], expected)
    cut_median = math.tan(discrete_mixture_quantile(.5, kind="normalized_cut", **arguments))
    factor_median = math.tan(discrete_mixture_quantile(.5, kind="factor_average", **arguments))
    assert cut_median == pytest.approx(1.9773123883291086, abs=2e-9)
    assert factor_median == pytest.approx(1.1753003504694248, abs=2e-9)


def test_component_and_both_mixtures_return_prior_on_a_zero_reading():
    arguments = dict(observed_angle=.8, radius=0., component_temperature=.3, **CONTROLS)
    component = gaussian_component(.6, variance=2., **arguments)
    assert component["evidence"] == 1 and component["cdf"] == .6/(math.pi/2)
    for kind in ("normalized_cut", "factor_average"):
        result = discrete_mixture_cdf(.6, variances=[.1, 20.], weights=[.3, .7], kind=kind, **arguments)
        assert result["cdf"] == pytest.approx(.6/(math.pi/2))
        np.testing.assert_allclose(result["effective_weights"], [.3, .7])


@pytest.mark.parametrize("nu", [.2, .5, 2.])
def test_normalized_ig_cut_single_reading_band_mass_increases(nu):
    probabilities = []
    for radius in (10., 100., 1000.):
        args = dict(observed_angle=math.atan(2), radius=radius, nu_eff=nu,
                    component_temperature=1., **CONTROLS)
        low, le = ig_normalized_cut_cdf(math.atan(1.9), **args)
        high, he = ig_normalized_cut_cdf(math.atan(2.1), **args)
        probabilities.append(high-low)
        assert max(le, he) < 2e-6
    assert probabilities[0] < probabilities[1] < probabilities[2]
    if nu == .5:
        np.testing.assert_allclose(probabilities, [.1120329820589, .5564132508974, .8572065100464], atol=2e-7)
        limit = heavy_tail_limit_cdf(math.atan(2.1), truth_angle=math.atan(2), tail_index=.5) - \
            heavy_tail_limit_cdf(math.atan(1.9), truth_angle=math.atan(2), tail_index=.5)
        assert probabilities[-1] > .8 and limit == pytest.approx(.16095150216195464, abs=2e-12)


def _polar_weighted_scatter(mean, nu):
    # Independent Gamma-radial/Legendre-angular normal quadrature. Split at
    # axes and diagonals so the positive-Q formula is smooth in each sector.
    radial, rw = special.roots_laguerre(100)
    nodes, weights = np.polynomial.legendre.leggauss(64)
    result = np.zeros((2, 2))
    for sector in range(8):
        low, high = sector*math.pi/4, (sector+1)*math.pi/4
        angles = (low+high)/2+(high-low)*nodes/2
        aw = (high-low)*weights/2
        directions = np.stack([np.cos(angles), np.sin(angles)], axis=1)
        y = np.sqrt(2*radial)[:, None, None]*directions[None, :, :]
        q = np.where(y[:, :, 0]*y[:, :, 1] > 0, 0., np.minimum(y[:, :, 0]**2, y[:, :, 1]**2))
        w = nu/(nu+q)
        density_shift = np.exp(-.5*(mean@mean)+y@mean)/(2*math.pi)
        result += np.einsum("r,p,rp,rpi,rpj->ij", rw, aw, w*density_shift, y, y)
    return result


def test_positive_q_bias_and_weak_truth_score_by_independent_quadrature():
    nu, radius, theta = 10., .001, math.atan(2)
    bias, error = positive_q_bias(nu, **CONTROLS)
    closed = (1-nu*math.exp(nu)*special.exp1(nu))/math.pi
    assert bias == pytest.approx(.026854742133967015, abs=2e-13)
    assert bias == pytest.approx(closed, abs=2e-13) and error < 1e-10
    null = _polar_weighted_scatter(np.zeros(2), nu)
    assert null[0, 1] == pytest.approx(bias, abs=2e-12)
    assert null[0, 0] == pytest.approx(null[1, 1], abs=2e-12)
    u = np.array([math.sin(theta), math.cos(theta)])
    v = np.array([math.cos(theta), -math.sin(theta)])
    shifted = _polar_weighted_scatter(radius*u, nu)
    null_score, shifted_score = float(u@null@v), float(u@shifted@v)
    bound = math.sqrt(math.expm1(radius**2))
    assert null_score == pytest.approx(-.6*bias, abs=2e-12)
    assert abs(shifted_score-null_score) <= bound
    assert null_score+bound < 0 and shifted_score < 0


def test_component_temperature_is_inside_the_gaussian_factor():
    args = dict(observed_angle=.9, radius=5., **CONTROLS)
    first = gaussian_component(.8, variance=2., component_temperature=.5, **args)
    second = gaussian_component(.8, variance=4., component_temperature=1., **args)
    assert first["cdf"] == pytest.approx(second["cdf"], abs=1e-12)
    assert first["evidence"] == pytest.approx(second["evidence"], abs=1e-12)
    # Independent von-Mises probability integral and scaled Bessel normalizer.
    kappa = (.5*25/2)/4
    low = stats.vonmises.cdf(-1.8, kappa)
    total = stats.vonmises.cdf(math.pi-1.8, kappa)-low
    expected = (stats.vonmises.cdf(2*(.8-.9), kappa)-low)/total
    assert first["cdf"] == pytest.approx(expected, abs=2e-12)
    assert first["evidence"] == pytest.approx(2*special.i0e(kappa)*total, abs=2e-12)


def test_explicit_controls_and_proper_frozen_q_are_required():
    with pytest.raises(TypeError):
        positive_q_bias(10.)
    with pytest.raises(ValueError):
        discrete_mixture_cdf(.5, observed_angle=.7, radius=2., variances=[1., 2.],
                             weights=[1., 1.], kind="normalized_cut", component_temperature=1., **CONTROLS)
    with pytest.raises(ValueError):
        ig_normalized_cut_cdf(.5, observed_angle=.7, radius=2., nu_eff=0., component_temperature=1., **CONTROLS)
