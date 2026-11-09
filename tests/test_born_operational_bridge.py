"""E6/T70: counterexamples delimiting the Born reconstruction bridge.

These test written mathematical statements, not a new mass equation.
See explore/e6-born-bridge.md for operational assumptions and primary sources.
"""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import multivariate_normal

from met.declared_structure import (calibrated_discrete_joint, direction_log_factors,
                                    update_discrete_joint)
from met.full import _h_terms, isotropic


def _projector(phi):
    v = np.array([np.sin(phi), np.cos(phi)])
    return np.outer(v, v)


def _born(state, effect):
    return float(np.trace(state @ effect).real)


def _g(p, alpha):
    return p**alpha / (p**alpha + (1-p)**alpha)


def test_t70_real_tests_do_not_select_a_unique_complex_completion():
    gram = np.array([[2.0, 0.4], [0.4, 1.0]])
    energy = np.trace(gram)
    limit = np.sqrt(np.linalg.det(gram))
    states = [np.array([[gram[0, 0], gram[0, 1]+1j*t],
                        [gram[1, 0]-1j*t, gram[1, 1]]])/energy
              for t in [-limit, 0.0, limit]]
    for state in states:
        assert np.linalg.eigvalsh(state).min() >= -1e-15
        assert np.trace(state) == pytest.approx(1)
    for phi in np.linspace(0.01, np.pi/2-0.01, 21):
        probabilities = [_born(state, _projector(phi)) for state in states]
        np.testing.assert_allclose(probabilities, probabilities[0], atol=1e-15)
    assert np.trace(states[0] @ states[0]).real == pytest.approx(1)
    assert np.trace(states[1] @ states[1]).real < 0.7
    assert np.trace(states[2] @ states[2]).real == pytest.approx(1)


def test_perpendicular_reading_has_real_mixed_and_pure_imaginary_encodings():
    mixed = np.eye(2)/2
    psi = np.array([1.0, 1j])/np.sqrt(2)
    pure = np.outer(psi, psi.conj())
    for phi in np.linspace(-np.pi, np.pi, 31):
        assert _born(mixed, _projector(phi)) == pytest.approx(0.5)
        assert _born(pure, _projector(phi)) == pytest.approx(0.5)
    # An additional imaginary-basis test distinguishes them, but is not P_m.
    assert _born(mixed, pure) == pytest.approx(0.5)
    assert _born(pure, pure) == pytest.approx(1)


def test_positive_mass_rays_do_not_resolve_identity_with_uniform_angle_weight():
    # Even scaling the total trace correctly leaves a nonzero off-diagonal.
    integrated = np.array([[quad(lambda phi: 4/np.pi*_projector(phi)[i, j],
                                 0, np.pi/2)[0] for j in range(2)] for i in range(2)])
    np.testing.assert_allclose(integrated, [[1, 2/np.pi], [2/np.pi, 1]], atol=1e-14)
    rng = np.random.default_rng(6401)
    phis = rng.uniform(0.001, np.pi/2-0.001, 20)
    weights = rng.uniform(0.01, 1, 20)
    weights *= 2/weights.sum()  # necessary trace condition for a POVM
    summed = sum(w*_projector(phi) for w, phi in zip(weights, phis))
    assert summed[0, 1] > 0
    for phi in phis:
        np.testing.assert_allclose(np.eye(2)-_projector(phi), _projector(phi+np.pi/2), atol=1e-15)
        assert np.tan(phi+np.pi/2) < 0


@pytest.mark.parametrize("alpha", [0.5, 1.0, 2.0, 3.0])
def test_binary_normalization_and_covariance_allow_non_born_rules(alpha):
    rng = np.random.default_rng(6402)
    for _ in range(8):
        psi = rng.normal(size=2)+1j*rng.normal(size=2)
        psi /= np.linalg.norm(psi)
        state = np.outer(psi, psi.conj())
        phi = rng.uniform(0, np.pi)
        effect = _projector(phi)
        p = np.clip(_born(state, effect), 0, 1)
        assert _g(p, alpha)+_g(1-p, alpha) == pytest.approx(1)
        u, _ = np.linalg.qr(rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2)))
        transformed = _born(u @ state @ u.conj().T, u @ effect @ u.conj().T)
        assert _g(np.clip(transformed, 0, 1), alpha) == pytest.approx(_g(p, alpha))
    assert _g(0.0, alpha) == 0 and _g(1.0, alpha) == 1


def test_non_born_projector_assignment_cannot_have_a_born_density():
    state = np.diag([1.0, 0.0])
    # Certainty on this basis projector would fix any representing density.
    assert _g(_born(state, state), 2) == 1
    effect = _projector(np.pi/6)
    assert _born(state, effect) == pytest.approx(0.25)
    assert _g(_born(state, effect), 2) == pytest.approx(0.1)


def test_effect_refinement_excludes_the_displayed_scalar_extension():
    trine = [(2/3)*_projector(phi) for phi in [0, np.pi/3, 2*np.pi/3]]
    np.testing.assert_allclose(sum(trine), np.eye(2), atol=5e-16)
    mixed = np.eye(2)/2
    assert sum(_born(mixed, effect) for effect in trine) == pytest.approx(1)
    assert sum(_g(_born(mixed, effect), 2) for effect in trine) == pytest.approx(3/5)
    pure = np.diag([0.0, 1.0])
    assert sum(_g(_born(pure, effect), 2) for effect in trine) == pytest.approx(57/65)
    # This is the failure of F(E)=g(tr rho E); a different maximally-mixed
    # projector extension can of course be Born. See the written scope.


def test_affine_preparation_mixing_is_an_extra_constraint():
    first, second = np.diag([1.0, 0.0]), np.diag([0.0, 1.0])
    mixture = 0.25*first+0.75*second
    effect = first
    assert _g(_born(mixture, effect), 2) == pytest.approx(0.1)
    mixed_assignments = 0.25*_g(_born(first, effect), 2)+0.75*_g(_born(second, effect), 2)
    assert mixed_assignments == pytest.approx(0.25)


def test_normalized_qubit_discards_energy_needed_by_gaussian_mass_evidence():
    states, contrasts = [], []
    for radius in [1.0, 2.0, 10.0]:
        psi = radius*np.array([2.0, 1.0])
        states.append(np.outer(psi, psi)/(psi @ psi))
        reading = isotropic([2*radius], [radius], 1, 1)
        h = _h_terms(reading.pre(), 1, np.arctan(np.array([1.0, 2.0])))
        contrasts.append((h[1]-h[0])/2)
    np.testing.assert_allclose(states, np.tile(states[0], (3, 1, 1)), atol=1e-15)
    np.testing.assert_allclose(contrasts, [0.25, 1, 25], atol=1e-13)


def test_diagonal_born_effects_are_classical_decisions_and_refine_additively():
    joint = np.array([[0.10, 0.15], [0.25, 0.50]])
    rho = np.diag(joint.ravel())
    first = np.diag([1.0, 0, 0, 0])
    second = np.diag([0.0, 1, 0, 0])
    merged = first+second
    assert _born(rho, merged) == pytest.approx(_born(rho, first)+_born(rho, second))
    assert _born(rho, merged)+_born(rho, np.eye(4)-merged) == pytest.approx(1)
    # This representation is classical: every state/effect used here is diagonal.
    np.testing.assert_array_equal(first @ second, second @ first)
    marginal = np.trace(rho.reshape(2, 2, 2, 2), axis1=1, axis2=3)
    np.testing.assert_allclose(marginal, np.diag(joint.sum(axis=1)))
    assert not np.allclose(rho, np.kron(marginal, np.diag(joint.sum(axis=0))))


@pytest.mark.parametrize("power", [0.5, 1.0, 2.0])
def test_arbitrary_positive_kernels_have_the_same_diagonal_luders_representation(power):
    joint = np.array([[0.10, 0.15], [0.25, 0.50]])
    kernel = np.array([[0.8, 0.2], [0.3, 0.7]])**power
    rho, root_effect = np.diag(joint.ravel()), np.diag(np.sqrt(kernel.ravel()))
    updated = root_effect @ rho @ root_effect
    updated /= np.trace(updated)
    classical = update_discrete_joint(joint, np.log(kernel))
    np.testing.assert_allclose(np.diag(updated), classical.ravel(), atol=1e-16)
    assert np.linalg.eigvalsh(np.diag(kernel.ravel())).min() > 0
    assert np.linalg.eigvalsh(np.eye(4)-np.diag(kernel.ravel())).min() >= 0


def test_e6_partial_trace_preserves_mass_local_tests_but_not_fixture_tests():
    prior, qcal = np.array([1, 2, 1])/4, np.array([0.25, 0.75])
    masses, directions = [0.5, 1.0, 2.0], np.eye(2)
    empty = calibrated_discrete_joint(prior, qcal, np.zeros((0, 3, 2)))
    block = direction_log_factors([[1, 0], [1, 0]], [[1, 0], [-1, 0]], masses, directions)
    after = calibrated_discrete_joint(prior, qcal, block)
    rho0, rho1 = np.diag(empty.ravel()), np.diag(after.ravel())
    local_effect = np.kron(np.diag([0.0, 0.3, 1.0]), np.eye(2))
    assert _born(rho0, local_effect) == pytest.approx(_born(rho1, local_effect))
    future = np.exp(direction_log_factors([[2, 0]], [[1, 0]], masses, directions)[0])
    fixture_effect = np.diag(future.ravel())
    assert abs(_born(rho0, fixture_effect)-_born(rho1, fixture_effect)) > 0.05


def test_profile_factors_are_not_normalized_raw_data_densities():
    m = 2.0
    u = np.array([m, 1])/np.sqrt(1+m*m)
    v = np.array([1, -m])/np.sqrt(1+m*m)
    for parallel in [-10.0, 0.0, 10.0]:
        y = parallel*u+0.7*v
        residual = (y[0]-m*y[1])**2/(1+m*m)
        assert residual == pytest.approx(0.7**2)
    # Integrating in the parallel coordinate grows with the cutoff forever.
    normal_integral = quad(lambda r: np.exp(-r*r/2), -np.inf, np.inf)[0]
    truncated = [2*cutoff*normal_integral for cutoff in [1, 2, 5]]
    np.testing.assert_allclose(truncated, np.sqrt(2*np.pi)*np.array([2, 4, 10]))


def test_data_only_weight_can_normalize_profile_ratios_in_a_different_experiment():
    # Nonintegrability under Lebesgue measure does not exclude every other
    # normalized sampling family with those likelihood ratios.
    tau = 2.0
    y = np.array([2.0, -0.3])
    costs, log_densities, constants = [], [], []
    for mass in [0.5, 1.0, 2.0, 5.0]:
        projector = _projector(np.arctan(mass))
        residual = np.eye(2)-projector
        covariance = tau**2*projector+tau**2/(1+tau**2)*residual
        cost = y @ residual @ y
        log_density = multivariate_normal.logpdf(y, mean=np.zeros(2), cov=covariance)
        costs.append(cost)
        log_densities.append(log_density)
        constants.append(log_density+cost/2+(y @ y)/(2*tau**2))
        # This changes the declared unit-noise residual variance.
        v = np.array([1.0, -mass])/np.sqrt(1+mass**2)
        assert v @ covariance @ v == pytest.approx(4/5)
    np.testing.assert_allclose(constants, constants[0], atol=1e-14)
    np.testing.assert_allclose(np.diff(log_densities), -np.diff(costs)/2, atol=1e-14)
