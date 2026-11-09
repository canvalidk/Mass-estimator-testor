"""E9: independent Fisher, spectral, covariance and prior-shape checks.

These check conditional constructions, not a new mass equation or prior.
"""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

from met.mass_geometry import (functional_information, independent_channel_ratios,
                                mass_curve, scalar_channel_metric, structural_fisher)


def _experiment(dimension, rank):
    rng = np.random.default_rng(6900+10*dimension+rank)
    raw = rng.normal(size=(2*dimension, 2*dimension))
    cov = raw @ raw.T+np.eye(2*dimension)
    n, _ = np.linalg.qr(rng.normal(size=(dimension, rank)), mode="reduced")
    return cov, n


def _normalizer(ratios, correlations=None, weights=None):
    return quad(lambda m: np.sqrt(scalar_channel_metric(m, ratios, correlations, weights)),
                0, np.inf, epsabs=1e-11, epsrel=1e-11)[0]


def _median(ratios, weights=None):
    total = _normalizer(ratios, weights=weights)
    return brentq(lambda m: quad(lambda t: np.sqrt(scalar_channel_metric(t, ratios, weights=weights)),
                                0, m, epsabs=1e-11)[0]-total/2, 0, 100)


@pytest.mark.parametrize("dimension,rank", [(1, 1), (3, 1), (3, 2), (3, 3)])
@pytest.mark.parametrize("mass", [0.4, 2.0])
def test_full_covariance_structural_fisher_from_numerical_covariance_derivatives(dimension, rank, mass):
    cov, n = _experiment(dimension, rank)
    state = mass_curve(mass, cov, n)
    step = 1e-4
    p = lambda m: mass_curve(m, cov, n)["projector"]
    numeric_p = (p(mass-2*step)-8*p(mass-step)+8*p(mass+step)-p(mass+2*step))/(12*step)
    np.testing.assert_allclose(state["derivative"], numeric_p, atol=7e-12)
    assert state["metric"] == pytest.approx(np.trace(np.linalg.solve(state["gram"],
                                               state["transverse_gram"])), rel=3e-13)
    for omega in [0.2, 1.0, 3.0]:
        variance = lambda m, w: np.eye(2*dimension)+w*w*p(m)
        v = variance(mass, omega)
        derivatives = [(variance(mass+step, omega)-variance(mass-step, omega))/(2*step),
                       (variance(mass, omega+step)-variance(mass, omega-step))/(2*step)]
        solved = [np.linalg.solve(v, derivative) for derivative in derivatives]
        information = np.array([[np.trace(a @ b)/2 for b in solved] for a in solved])
        np.testing.assert_allclose(information, structural_fisher(mass, omega, cov, n),
                                   atol=3e-8, rtol=2e-7)
        # This is an actual added Gaussian population, in physical coordinates.
        b = np.linalg.solve(np.linalg.cholesky(cov), np.vstack([mass*n, n]))
        latent_cov = omega**2*np.linalg.inv(state["gram"])
        np.testing.assert_allclose(b @ latent_cov @ b.T, omega**2*state["projector"], atol=2e-14)


def test_functional_information_from_full_numerical_mean_jacobian():
    cov, n = _experiment(3, 2)
    mass, alpha = 1.7, np.array([0.4, -1.2])
    theta = np.r_[mass, alpha]
    mean = lambda t: np.r_[t[0]*(n @ t[1:]), n @ t[1:]]
    step = 1e-5
    jacobian = np.column_stack([(mean(theta+step*unit)-mean(theta-step*unit))/(2*step)
                               for unit in np.eye(3)])
    fisher = jacobian.T @ np.linalg.solve(cov, jacobian)
    schur = fisher[0, 0]-fisher[0, 1:] @ np.linalg.solve(fisher[1:, 1:], fisher[1:, 0])
    assert functional_information(mass, alpha, cov, n) == pytest.approx(schur, rel=3e-11)
    # Basis changes at the same physical point do not change efficient Fisher.
    change = np.array([[1.2, 0.3], [-0.2, 0.8]])
    assert functional_information(mass, np.linalg.solve(change, alpha), cov, n @ change) == pytest.approx(schur, rel=3e-11)


@pytest.mark.parametrize("correlation", [-0.7, 0.0, 0.6])
def test_correlated_scalar_closed_metric_matches_whitened_projectors(correlation):
    sf, sa = 2.4, 0.8
    covariance = np.array([[sf*sf, correlation*sf*sa], [correlation*sf*sa, sa*sa]])
    for mass in [0.1, 0.7, 3, 12]:
        actual = mass_curve(mass, covariance)["metric"]
        assert actual == pytest.approx(scalar_channel_metric(mass, [sf/sa], [correlation]), rel=2e-13)


@pytest.mark.parametrize("rank", [1, 2, 3])
def test_anisotropic_noncommuting_channels_reduce_to_generalized_noise_ratios(rank):
    rng = np.random.default_rng(6920+rank)
    raw1, raw2 = rng.normal(size=(2, 3, 3))
    force, acceleration = raw1 @ raw1.T+np.eye(3), raw2 @ raw2.T+0.5*np.eye(3)
    assert np.linalg.norm(force @ acceleration-acceleration @ force) > 0.1
    n, _ = np.linalg.qr(rng.normal(size=(3, rank)), mode="reduced")
    ratios = independent_channel_ratios(force, acceleration, n)
    covariance = np.block([[force, np.zeros((3, 3))], [np.zeros((3, 3)), acceleration]])
    for mass in [0.1, 1, 7]:
        assert mass_curve(mass, covariance, n)["metric"] == pytest.approx(scalar_channel_metric(mass, ratios), rel=3e-13)
    if rank < 3:
        # A projected marginal loses known-zero outside-mean noise information.
        marginal = np.block([[n.T @ force @ n, np.zeros((rank, rank))],
                             [np.zeros((rank, rank)), n.T @ acceleration @ n]])
        assert abs(mass_curve(1, covariance, n)["metric"]-mass_curve(1, marginal)["metric"]) > 1e-4


def test_graph_noise_scale_invariance_differs_from_functional_precision():
    cov, n = _experiment(3, 2)
    alpha = np.array([0.7, -0.4])
    for scale in [0.2, 3.0, 10.0]:
        assert mass_curve(1.3, scale*cov, n)["metric"] == pytest.approx(mass_curve(1.3, cov, n)["metric"], rel=3e-14)
        assert functional_information(1.3, alpha, scale*cov, n) == pytest.approx(
            functional_information(1.3, alpha, cov, n)/scale, rel=3e-14)


def test_units_and_physical_swap_transform_full_covariance_metric():
    cov, n = _experiment(3, 2)
    unit_change = np.diag([2.5]*3+[0.4]*3)
    ratio = 2.5/0.4
    swapped = np.block([[np.zeros((3, 3)), np.eye(3)], [np.eye(3), np.zeros((3, 3))]])
    for mass in [0.2, 1.0, 5.0]:
        metric = mass_curve(mass, cov, n)["metric"]
        assert mass_curve(ratio*mass, unit_change @ cov @ unit_change.T, n)["metric"]*ratio**2 == pytest.approx(metric, rel=4e-13)
        assert mass_curve(1/mass, swapped @ cov @ swapped.T, n)["metric"]/mass**4 == pytest.approx(metric, rel=4e-13)


def test_signal_path_choice_changes_volume_not_information_at_the_same_point():
    # One unit instrument: G=1+m², E=1/(1+m²).
    for mass in [0.1, 1, 10, 100]:
        assert functional_information(mass, [2.0], np.eye(2)) == pytest.approx(4/(1+mass**2))
        alpha_at_fixed_radius = 2/np.sqrt(1+mass**2)
        assert functional_information(mass, [alpha_at_fixed_radius], np.eye(2)) == pytest.approx(4/(1+mass**2)**2)
    # The former path has logarithmically divergent length; the latter is proper.
    partial = [quad(lambda m: np.sqrt(functional_information(m, [2.0], np.eye(2))), 0, limit)[0]
               for limit in [10, 100]]
    np.testing.assert_allclose(partial, 2*np.arcsinh([10, 100]), rtol=1e-13)
    assert quad(lambda m: 2/(1+m*m), 0, np.inf)[0] == pytest.approx(np.pi)


def test_common_structural_size_gives_product_geometry_even_when_instruments_differ():
    ratios = np.array([1.0, 4.0])
    for mass in [0.2, 1, 4, 20]:
        for omega in [0.3, 1, 3]:
            actual = sum(structural_fisher(mass, omega, np.diag([s*s, 1]))[0, 0] for s in ratios)
            expected = omega**4/(1+omega**2)*scalar_channel_metric(mass, ratios)
            assert actual == pytest.approx(expected, rel=3e-14)


def test_same_structural_family_has_different_fixed_nuisance_sections():
    # omega=tau*sqrt(1+m²) relates noise-radius size to physical acceleration size.
    tau, step = 0.8, 1e-5
    variance = lambda m, t: np.eye(2)+t*t*np.outer([m, 1.0], [m, 1.0])
    for mass in [0.3, 1.0, 4.0]:
        v = variance(mass, tau)
        derivatives = [(variance(mass+step, tau)-variance(mass-step, tau))/(2*step),
                       (variance(mass, tau+step)-variance(mass, tau-step))/(2*step)]
        solved = [np.linalg.solve(v, derivative) for derivative in derivatives]
        fisher = np.array([[np.trace(a @ b)/2 for b in solved] for a in solved])
        schur = fisher[0, 0]-fisher[0, 1]**2/fisher[1, 1]
        assert schur == pytest.approx(tau**4/(1+tau**2*(1+mass**2)), rel=2e-9)
        omega = tau*np.sqrt(1+mass**2)
        assert schur == pytest.approx(structural_fisher(mass, omega, np.eye(2))[0, 0], rel=2e-9)
        jacobian = np.array([[1, 0], [tau*mass/np.sqrt(1+mass**2), np.sqrt(1+mass**2)]])
        np.testing.assert_allclose(fisher, jacobian.T @ structural_fisher(mass, omega, np.eye(2)) @ jacobian,
                                   rtol=2e-9, atol=2e-10)
    # Fixed tau has a 1/m efficient-Fisher speed tail; fixed omega has 1/m².
    mass = 1e6
    assert mass*np.sqrt(tau**4/(1+tau**2*(1+mass**2))) == pytest.approx(tau, rel=2e-12)
    assert mass**2*np.sqrt(structural_fisher(mass, 1, np.eye(2))[0, 0]) == pytest.approx(1/np.sqrt(2), rel=1e-8)


def test_equal_metric_profiles_cancel_arbitrary_positive_size_allocations():
    masses = np.geomspace(0.01, 100, 51)
    base = scalar_channel_metric(masses, [2.0], [0.3])
    for weights in [[0.2, 8], [9, 1], [2, 2]]:
        mixture = scalar_channel_metric(masses, [2, 2], [0.3, 0.3], weights)
        np.testing.assert_allclose(mixture, sum(weights)*base, rtol=2e-15)


def test_orthogonal_size_nuisances_still_change_heterogeneous_mass_volume():
    # lambda(omega) can produce any positive weight: x=omega² solves x²/(1+x)=w.
    ratios, allocations = [1.0, 4.0], [[9, 1], [1, 9]]
    density_ratios, medians = [], []
    for weights in allocations:
        size2 = (np.array(weights)+np.sqrt(np.array(weights)**2+4*np.array(weights)))/2
        for mass in [1.0, 4.0]:
            fisher = [structural_fisher(mass, np.sqrt(x), np.diag([s*s, 1]))
                      for s, x in zip(ratios, size2)]
            assert all(info[0, 1] == 0 for info in fisher)
            assert sum(info[0, 0] for info in fisher) == pytest.approx(scalar_channel_metric(mass, ratios, weights=weights))
        density_ratios.append(np.sqrt(scalar_channel_metric(1, ratios, weights=weights)/
                                       scalar_channel_metric(4, ratios, weights=weights)))
        medians.append(_median(ratios, weights))
    assert abs(density_ratios[0]-density_ratios[1]) > 4
    assert medians[0] < medians[1]
    assert medians[0]*medians[1] == pytest.approx(4, rel=2e-12)


def test_distinct_scalar_noise_quadratics_have_uncancelled_double_poles():
    ratios, weights = [1, 2], [0.7, 3]
    delta = 1e-7
    for ratio, weight in zip(ratios, weights):
        pole = 1j*ratio
        value = sum(w*s*s/((pole+delta)**2+s*s)**2 for s, w in zip(ratios, weights))
        # A single SPD ray metric has just one upper-half-plane pole.
        assert delta**2*value == pytest.approx(-weight/4, abs=1e-7)


@pytest.mark.parametrize("ratios,correlations,weights,centre", [
    ([1, 4], [0.2, 0.2], [2, 2], 2),
    ([2, 2], [-0.5, 0.6], [1, 3], 2),
    ([1, 2, 4], [0, 0.3, 0], [1, 7, 1], 2),
])
def test_reciprocal_weighted_multiset_positive_controls(ratios, correlations, weights, centre):
    masses = np.geomspace(0.02, 100, 81)
    metric = scalar_channel_metric(masses, ratios, correlations, weights)
    reflected = scalar_channel_metric(centre**2/masses, ratios, correlations, weights)*centre**4/masses**4
    np.testing.assert_allclose(metric, reflected, rtol=2e-15)
    total = _normalizer(ratios, correlations, weights)
    below_centre = quad(lambda m: np.sqrt(scalar_channel_metric(m, ratios, correlations, weights)),
                        0, centre, epsabs=1e-11)[0]
    assert below_centre/total == pytest.approx(0.5, abs=2e-13)


def test_three_instrument_curve_has_no_l4_centre_despite_unit_swap_covariance():
    ratios = np.array([1.0, 2.0, 3.0])
    # Pole reflection must reverse ordered ratios; no constant c² can do so.
    products = ratios*ratios[::-1]
    np.testing.assert_array_equal(products, [3, 4, 3])
    total = _normalizer(ratios)
    assert total == pytest.approx(2.846513208821, abs=8e-12)
    assert _median(ratios) == pytest.approx(1.817989432560, abs=8e-12)
    # Endpoint density matching already fixes this NECESSARY candidate centre.
    centre = (np.sum(ratios**2)/np.sum(ratios**-2))**0.25
    direct = np.sqrt(scalar_channel_metric(1, ratios))/total
    image = np.sqrt(scalar_channel_metric(centre**2, ratios))*centre**2/total
    assert direct == pytest.approx(0.248411558041, abs=8e-12)
    assert image == pytest.approx(0.256055611015, abs=8e-12)
    assert abs(direct-image) > 0.007
    # Across transformed declarations the normalized measure is still covariant.
    reciprocal_total = _normalizer(1/ratios)
    assert reciprocal_total == pytest.approx(total, rel=2e-13)
    for mass in [0.2, 1.0, 7.0]:
        assert scalar_channel_metric(mass, ratios) == pytest.approx(
            scalar_channel_metric(1/mass, 1/ratios)/mass**4, rel=1e-14)


def test_full_covariance_curve_is_smooth_positive_with_integrable_end_tails():
    cov, n = _experiment(3, 2)
    swap = np.block([[np.zeros((3, 3)), np.eye(3)], [np.eye(3), np.zeros((3, 3))]])
    small = mass_curve(1e-6, cov, n)["metric"]
    large = mass_curve(1e6, cov, n)["metric"]*1e24
    expected_large = mass_curve(1e-6, swap @ cov @ swap.T, n)["metric"]
    assert small > 0 and large > 0
    assert large == pytest.approx(expected_large, rel=2e-8)
    total = quad(lambda m: np.sqrt(mass_curve(m, cov, n)["metric"]), 0, np.inf,
                 epsabs=1e-10, epsrel=1e-10)[0]
    assert 0 < total < 20
