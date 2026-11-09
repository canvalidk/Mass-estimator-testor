"""Independent numerical witnesses for the E10 null-limit exploration.

These finite quadratures check identities and counterexamples, not the
full-support obstruction's quantification over all probability measures.
"""
import math

import numpy as np
import pytest
from scipy.special import roots_legendre

from met.gaussian_information import flat_gaussian_fit
from met.null_limit import average_null_update


def _positive_prior(order=120, shape=0.):
    nodes, weights = roots_legendre(order)
    x = (nodes+1)*math.pi/2
    p = weights*(1+shape*np.cos(2*x))
    return x, p/p.sum()


def test_fixed_strength_matches_t42_sine_squared_kernel():
    x, p = _positive_prior(32, .3)
    phase_count, t = 128, 1.7
    phi = 2*math.pi*np.arange(phase_count)/phase_count
    # T42 has K=4*t and theta=x/2; the omitted exp(-t) cancels.
    factors = np.exp(-2*t*np.sin((x[:, None]-phi)/2)**2)
    joint = p[:, None]*factors
    expected = (joint/joint.sum(axis=0)).mean(axis=1)
    result = average_null_update(x, p, fixed_strength=t, phase_order=phase_count)
    np.testing.assert_allclose(result['weights'], expected, atol=2e-16)
    np.testing.assert_allclose(result['response']*p, expected, atol=2e-16)


@pytest.mark.parametrize('shape', [0., .6, -.6])
def test_endpoint_response_derivative_is_negative_for_full_support_shapes(shape):
    x, p = _positive_prior(120, shape)
    h = 1e-5
    a = average_null_update(x, p, queries=[-h, 0., h, math.pi],
                            radial_order=64, phase_order=384)
    assert a['derivative'][1] < -.5
    assert a['derivative'][3] > .5
    numerical = (a['response'][2]-a['response'][0])/(2*h)
    assert numerical == pytest.approx(a['derivative'][1], rel=2e-8)
    assert a['response'][1] == pytest.approx(a['response'][3], abs=2e-13)
    np.testing.assert_allclose(a['weights'], a['weights'][::-1], atol=3e-16)
    assert a['weights'].sum() == pytest.approx(1., abs=3e-15)


def test_half_cauchy_is_not_fixed_and_quadrature_refinement_preserves_drift():
    results = []
    for order, radial, phase in [(120, 48, 384), (180, 80, 512)]:
        x, p = _positive_prior(order)
        a = average_null_update(x, p, queries=[0., math.pi/2, math.pi],
                                radial_order=radial, phase_order=phase)
        results.append(a)
        tv = .5*np.abs(a['weights']-p).sum()
        assert .0954 < tv < .0957
        assert a['response'][0] > 1.7
        assert a['response'][1] < .8
    np.testing.assert_allclose(results[0]['response'], results[1]['response'], atol=2e-7)
    np.testing.assert_allclose(results[0]['derivative'], results[1]['derivative'], atol=2e-7)


def test_whole_circle_uniform_is_fixed_but_restriction_before_update_is_not():
    # Cyclic shifts of a uniform grid give an exact finite-quadrature control.
    x = 2*math.pi*np.arange(96)/96
    p = np.full(x.shape, 1/96)
    a = average_null_update(x, p, radial_order=32, phase_order=384)
    np.testing.assert_allclose(a['weights'], p, atol=4e-16)
    # Half-Cauchy is the restriction of the continuous uniform circle law.
    y, q = _positive_prior(100)
    b = average_null_update(y, q, radial_order=48, phase_order=384)
    assert np.abs(b['weights']-q).sum() > .18


@pytest.mark.parametrize('support', [[.7], [.4, 2.3]])
def test_atomic_fixed_points_do_not_contradict_full_support_obstruction(support):
    # A singleton is samplewise fixed; two equally weighted atoms exchange
    # under a reflection of the uniformly distributed phase.
    p = np.full(len(support), 1/len(support))
    a = average_null_update(support, p, radial_order=48, phase_order=2048)
    np.testing.assert_allclose(a['weights'], p, atol=2e-12)


def test_zero_strength_preserves_every_prior():
    x, p = _positive_prior(48, .5)
    a = average_null_update(x, p, beta=0., queries=[0., .7, math.pi])
    np.testing.assert_allclose(a['weights'], p, atol=2e-16)
    np.testing.assert_allclose(a['response'], 1., atol=2e-15)
    np.testing.assert_allclose(a['derivative'], 0., atol=0.)


def test_small_strength_matches_independent_second_order_expansion():
    x, p = _positive_prior(80)
    beta = .002
    a = average_null_update(x, p, beta=beta, queries=[0., math.pi/2],
                            radial_order=24, phase_order=128)
    expected = np.array([16/math.pi**2, 16/math.pi**2-8/math.pi])
    np.testing.assert_allclose((a['response']-1)/beta**2, expected,
                               rtol=2e-4, atol=2e-4)
    assert a['derivative'][0]/beta**2 == pytest.approx(-8/math.pi, rel=2e-4)


def test_gram_and_bartlett_strengths_have_the_claimed_denominators():
    # Algebra at an ordinary realized experiment. Their common null limit
    # uses separate LLNs, not this finite-data equality check.
    y = np.array([[2., .4], [-.3, 1.1], [.2, -.7], [1.2, .6]])
    fit = flat_gaussian_fit(y)
    h, j, largest = fit['curv'], fit['J'], fit['eigenvalues'][1]
    eta_g = 1-fit['eigenvalues'][0]/largest
    assert eta_g*h/4 == pytest.approx(h*h/(4*largest))
    assert fit['precision']/4 == pytest.approx(min(h/4, h*h/(4*j)))


@pytest.mark.parametrize('kwargs', [dict(beta=-1), dict(fixed_strength=-1),
                                   dict(phase_order=0)])
def test_invalid_quadrature_inputs_are_rejected(kwargs):
    with pytest.raises(ValueError):
        average_null_update([.2, .9], [.5, .5], **kwargs)
