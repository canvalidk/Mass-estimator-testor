"""Conditional-score candidate checks, including its false stationary branch."""
import math

import numpy as np
import pytest

from met.conditional_score import gaussian_score_test
from met.full import Reading
from met.joint_mass import (FiniteJointMassLaw, common_acceleration,
                            gaussian_joint_posterior, gaussian_profile_cost,
                            separate_accelerations)
from met.joint_score_law import (composite_quantile, gaussian_joint_score_posterior,
                                 gaussian_joint_score_statistic,
                                 scalar_common_score_from_scatter)
from met.operator_score import gaussian_operator_score_test


RANK = dict(rank_rtol=1e-11, rank_atol=0.)


def _covariance(size):
    factor = np.eye(size)+np.arange(size*size).reshape(size, size)/(7*size*size)
    return factor@factor.T


@pytest.mark.parametrize("mass", [.2, 2., 7.])
def test_scalar_restriction_matches_existing_exact_score_references(mass):
    forces = [[2., -.3], [.4, 1.2], [-.7, .2]]
    accelerations = [[.7, -.2], [.1, .4], [.3, .6]]
    covariance = _covariance(4)
    events = [common_acceleration([f], a, covariance) for f, a in zip(forces, accelerations)]
    legacy = [Reading(f, a, covariance) for f, a in zip(forces, accelerations)]
    actual = gaussian_joint_score_statistic(events, [mass], **RANK)
    scalar = gaussian_score_test(legacy, mass)
    operator = gaussian_operator_score_test(legacy, mass*np.eye(2),
                                            tangents=np.eye(2)[None], **RANK)
    assert actual["rank"] == 1
    assert actual["statistic"] == pytest.approx(scalar["statistic"], abs=3e-14)
    assert actual["statistic"] == pytest.approx(operator["statistic"], abs=3e-14)
    assert actual["score"][0]*(1+mass*mass) == pytest.approx(scalar["score"], abs=2e-14)


def test_fixed_projection_conditional_gaussian_integration_has_the_observed_rank():
    # Three events yield six residual coordinates but a two-dimensional
    # conditional score subspace. No signal population or true sizes are fitted.
    masses = np.array([1.2, .8])
    covariance = _covariance(3)
    chol = np.linalg.cholesky(covariance)
    declaration = common_acceleration([[0.], [0.]], [0.], covariance)
    rotation, _, _ = np.linalg.svd(np.linalg.solve(chol, declaration.design(masses)), full_matrices=True)
    nodes, weights = np.polynomial.hermite.hermgauss(3)
    expected_score, expected_second = np.zeros(2), np.zeros((2, 2))
    first_statistic, second_statistic = 0., 0.
    fixed_covariance = None
    for indices in np.ndindex(*(3,)*6):
        probability = math.prod(weights[list(indices)])/math.pi**3
        events = []
        for j, projection in enumerate((.7, -1.1, .4)):
            residual = math.sqrt(2)*rotation[:, 1:]@nodes[list(indices[2*j:2*j+2])]
            data = chol@(projection*rotation[:, 0]+residual)
            events.append(common_acceleration(data[:2, None], data[2:], covariance))
        result = gaussian_joint_score_statistic(events, masses, **RANK)
        assert result["rank"] == 2
        if fixed_covariance is None:
            fixed_covariance = result["covariance"]
        else:
            np.testing.assert_allclose(result["covariance"], fixed_covariance, atol=2e-14)
        expected_score += probability*result["score"]
        expected_second += probability*np.outer(result["score"], result["score"])
        first_statistic += probability*result["statistic"]
        second_statistic += probability*result["statistic"]**2
    np.testing.assert_allclose(expected_score, 0., atol=2e-15)
    np.testing.assert_allclose(expected_second, fixed_covariance, atol=2e-14)
    assert first_statistic == pytest.approx(2., abs=2e-13)
    assert second_statistic == pytest.approx(8., abs=2e-12)


def test_one_common_scalar_event_score_law_equals_profile_law_away_from_zero_projection():
    event = common_acceleration([[2.], [3.]], [1.], _covariance(3))
    prior = FiniteJointMassLaw([[1., 1.], [1., 2.], [2., 1.], [2., 2.]], [.1, .2, .3, .4])
    for masses in prior.atoms:
        score = gaussian_joint_score_statistic([event], masses, **RANK)
        assert score["rank"] == 2
        assert score["statistic"] == pytest.approx(gaussian_profile_cost(event, masses), abs=2e-14)
    score_law = gaussian_joint_score_posterior([event], prior, **RANK)
    profile_law = gaussian_joint_posterior([event], prior)
    np.testing.assert_allclose(score_law.probabilities, profile_law.probabilities, atol=2e-15)
    assert score_law.log_evidence <= 0
    assert score_law.probabilities.sum() == pytest.approx(1., abs=1e-15)


def test_exact_zero_padding_and_empty_series_return_the_supplied_joint_law():
    prior = FiniteJointMassLaw([[1., 1.], [1., 2.], [2., 2.]], [.2, .3, .5])
    events = [common_acceleration([[2.], [3.]], [1.], np.eye(3)),
              common_acceleration([[.4], [1.7]], [.8], _covariance(3))]
    zero = common_acceleration([[0.], [0.]], [0.], np.diag([2., 3., 4.]))
    before = gaussian_joint_score_posterior(events, prior, **RANK)
    after = gaussian_joint_score_posterior(events+[zero], prior, **RANK)
    np.testing.assert_array_equal(before.probabilities, after.probabilities)
    assert before.log_evidence == after.log_evidence
    for series in ([], [zero]):
        unchanged = gaussian_joint_score_posterior(series, prior, **RANK)
        np.testing.assert_allclose(unchanged.probabilities, prior.probabilities, atol=1e-16)
        result = gaussian_joint_score_statistic(series, [1., 2.], **RANK)
        assert result["rank"] == 0 and result["statistic"] == 0 and result["p_value"] == 1


def test_exceptional_zero_projection_convention_is_visible_to_an_atomic_prior():
    event = common_acceleration([[1.]], [-2.], np.eye(2))
    result = gaussian_joint_score_statistic([event], [2.], **RANK)
    assert result["rank"] == 0 and result["statistic"] == 0
    # This is the declared exceptional value, not a continuity assertion.
    nearby = gaussian_joint_score_statistic([event], [2.001], **RANK)
    assert nearby["rank"] == 1 and nearby["statistic"] > 4.9
    prior = FiniteJointMassLaw([[1.], [2.]], [.5, .5])
    posterior = gaussian_joint_score_posterior([event], prior, **RANK)
    assert posterior.probabilities[1] > .9


def test_full_source_covariance_units_and_frame_transform_the_same_joint_law():
    masses = np.array([1.3, .8])
    data = [np.array([2., -.3, .4, 1.2, .7, -.2]),
            np.array([.5, .4, -.2, .7, .1, .3])]
    covariance = _covariance(6)
    events = [common_acceleration(y[:4].reshape(2, 2), y[4:], covariance) for y in data]
    rotation = np.array([[0., -1.], [1., 0.]])
    units = [3., 7., .4]
    transform = np.zeros((6, 6))
    for j, unit in enumerate(units):
        transform[2*j:2*j+2, 2*j:2*j+2] = unit*rotation
    changed = []
    for y in data:
        output = transform@y
        changed.append(common_acceleration(output[:4].reshape(2, 2), output[4:],
                                            transform@covariance@transform.T))
    jacobian = np.diag(np.asarray(units[:2])/units[2])
    original = gaussian_joint_score_statistic(events, masses, **RANK)
    converted = gaussian_joint_score_statistic(changed, jacobian@masses, **RANK)
    inverse = np.linalg.inv(jacobian)
    assert converted["statistic"] == pytest.approx(original["statistic"], abs=2e-13)
    np.testing.assert_allclose(converted["score"], inverse@original["score"], atol=2e-14)
    np.testing.assert_allclose(converted["covariance"], inverse@original["covariance"]@inverse, atol=2e-14)
    prior = FiniteJointMassLaw([[.7, .6], [1., 2.], [2., 3.]], [.2, .3, .5])
    converted_prior = FiniteJointMassLaw(prior.atoms@jacobian, prior.probabilities)
    first = gaussian_joint_score_posterior(events, prior, **RANK)
    second = gaussian_joint_score_posterior(changed, converted_prior, **RANK)
    np.testing.assert_allclose(first.probabilities, second.probabilities, atol=3e-15)


def test_reciprocal_swap_for_the_declared_separate_acceleration_protocol():
    covariance = _covariance(4)
    permutation = [2, 3, 0, 1]
    events = [separate_accelerations([[2.5], [-.8]], [[.7], [1.2]], covariance),
              separate_accelerations([[.4], [1.7]], [[1.1], [-.3]], covariance)]
    swapped = [separate_accelerations(e.data[2:, None], e.data[:2, None],
                                      e.covariance[np.ix_(permutation, permutation)]) for e in events]
    masses = np.array([2., .7])
    assert gaussian_joint_score_statistic(events, masses, **RANK)["statistic"] == pytest.approx(
        gaussian_joint_score_statistic(swapped, 1/masses, **RANK)["statistic"], abs=2e-13)
    prior = FiniteJointMassLaw([[.5, .8], [1., 2.], [2., 3.]], [.2, .3, .5])
    inverse_prior = FiniteJointMassLaw(1/prior.atoms, prior.probabilities)
    original = gaussian_joint_score_posterior(events, prior, **RANK)
    inverse_law = gaussian_joint_score_posterior(swapped, inverse_prior, **RANK)
    lookup = {tuple(atom): p for atom, p in zip(inverse_law.atoms, inverse_law.probabilities)}
    for atom, probability in zip(original.atoms, original.probabilities):
        assert lookup[tuple(1/atom)] == pytest.approx(probability, abs=2e-15)


@pytest.mark.parametrize("covariance", [np.eye(2), np.array([[1., .6], [.6, 1.5]])])
def test_vectorized_scatter_statistic_matches_the_own_reading_svd(covariance):
    data = np.array([[2., .7], [.4, 1.1], [-.3, .5], [1.2, -.8]])
    events = [common_acceleration([[f]], [a], covariance) for f, a in data]
    masses = np.array([.1, .4, 1., 2., 7., 1e8])
    fast = scalar_common_score_from_scatter(data.T@data, masses, covariance=covariance, **RANK)
    for j, mass in enumerate(masses):
        reference = gaussian_joint_score_statistic(events, [mass], **RANK)
        assert fast["statistic"][j] == pytest.approx(reference["statistic"], rel=2e-13, abs=2e-14)
        assert fast["score"][j] == pytest.approx(reference["score"][0], rel=2e-13, abs=2e-14)
        assert fast["variance"][j] == pytest.approx(reference["covariance"][0, 0], rel=2e-13, abs=2e-14)
    scalar = scalar_common_score_from_scatter(data.T@data, 2., covariance=covariance, **RANK)
    assert isinstance(scalar["statistic"], float) and scalar["rank"] == 1


def test_positive_correlated_covariance_false_stationary_root_is_not_removed():
    covariance = np.array([[1., .9], [.9, 1.]])
    mean = np.array([6/5, 1.])
    scatter = covariance+np.outer(mean, mean)
    candidates = np.array([4/15, 6/5])
    result = scalar_common_score_from_scatter(scatter, candidates, covariance=covariance,
                                             rank_rtol=0., rank_atol=0.)
    assert np.all(result["rank"] == 1)
    assert np.all(result["statistic"] < 1e-28)
    # Both stationary branches remain equal-height maxima of exp(-T/2).
    np.testing.assert_allclose(np.exp(-result["statistic"]/2), 1., atol=0.)


def test_composite_quantiles_use_the_supplied_joint_pushforward():
    law = FiniteJointMassLaw([[1., 1.], [2., 2.]], [.4, .6])
    assert composite_quantile(law, .5, operation="addition") == 4.
    np.testing.assert_allclose(composite_quantile(law, [.2, .5, .9], operation="parallel"), [.5, 1., 1.])


def test_explicit_rank_controls_and_invalid_inputs():
    event = common_acceleration([[2.], [3.]], [1.], np.eye(3))
    result = gaussian_joint_score_statistic([event], [1., 2.], rank_rtol=0., rank_atol=1e6)
    assert result["rank"] == 0 and result["statistic"] == 0
    zero = scalar_common_score_from_scatter(np.zeros((2, 2)), [1., 2.], covariance=np.eye(2), **RANK)
    assert np.all(zero["rank"] == 0) and np.all(zero["statistic"] == 0)
    with pytest.raises(ValueError):
        gaussian_joint_score_statistic([event], [0., 1.], **RANK)
    with pytest.raises(ValueError):
        gaussian_joint_score_statistic([event], [1., 2.], rank_rtol=-1., rank_atol=0.)
    with pytest.raises(ValueError):
        scalar_common_score_from_scatter(np.diag([1., -1.]), 1., covariance=np.eye(2), **RANK)
    with pytest.raises(ValueError):
        composite_quantile(FiniteJointMassLaw([[1.]], [1.]), 0., operation="addition")
