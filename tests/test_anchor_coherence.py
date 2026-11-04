"""T62: posterior-only closure versus additive retained calibration state."""
import itertools
import math

import numpy as np
import pytest
from scipy import special

from met.joint_mass import (FiniteJointMassLaw, common_acceleration,
                            gaussian_joint_posterior, gaussian_profile_cost)
from met.joint_score_law import (gaussian_joint_score_posterior,
                                 gaussian_joint_score_statistic)
from met.joint_trace_anchor import joint_trace_anchor


RANK = dict(rank_rtol=1e-11, rank_atol=0.)


def _event(force, acceleration):
    return common_acceleration([[force]], [acceleration], np.eye(2))


def _log_odds(law):
    return math.log(law.probabilities[1]/law.probabilities[0])


@pytest.mark.parametrize("prefactor", [.001, .3, 1., 4.])
def test_every_positive_trace_prefactor_fails_posterior_only_sequential_closure(prefactor):
    # Equal observed numbers on two independently declared events. They are
    # not repeated use of one source observation.
    event = _event(1., 1.)
    prior = FiniteJointMassLaw([[1.], [2.]], [.5, .5])
    batch, _ = joint_trace_anchor([event, event], prior,
                                  anchor_prefactor=prefactor, **RANK)
    first, _ = joint_trace_anchor([event], prior,
                                  anchor_prefactor=prefactor, **RANK)
    sequential, _ = joint_trace_anchor([event], first,
                                       anchor_prefactor=prefactor, **RANK)
    expected_batch = -(1+prefactor/math.sqrt(5))/5
    expected_sequential = -(1+prefactor/math.sqrt(3))/5
    assert _log_odds(batch) == pytest.approx(expected_batch, abs=2e-15)
    assert _log_odds(sequential) == pytest.approx(expected_sequential, abs=2e-15)
    assert _log_odds(batch)-_log_odds(sequential) == pytest.approx(
        prefactor*(1/math.sqrt(3)-1/math.sqrt(5))/5, abs=3e-15)
    assert batch.probabilities[1] > sequential.probabilities[1]


def test_unanchored_score_pooling_also_fails_posterior_only_sequential_closure():
    first_event, second_event = _event(2., 1.), _event(1., 2.)
    prior = FiniteJointMassLaw([[1.], [2.]], [.5, .5])
    actual = np.array([[gaussian_joint_score_statistic(events, [mass], **RANK)["statistic"]
                        for mass in (1., 2.)]
                       for events in ([first_event], [second_event], [first_event, second_event])])
    np.testing.assert_allclose(actual, [[.5, 0.], [.5, 9/5], [0., 144/205]], atol=3e-15)
    batch = gaussian_joint_score_posterior([first_event, second_event], prior, **RANK)
    first = gaussian_joint_score_posterior([first_event], prior, **RANK)
    sequential = gaussian_joint_score_posterior([second_event], first, **RANK)
    assert _log_odds(batch) == pytest.approx(-72/205, abs=2e-15)
    assert _log_odds(sequential) == pytest.approx(-2/5, abs=2e-15)
    assert _log_odds(batch)-_log_odds(sequential) == pytest.approx(2/41, abs=2e-15)


def test_profile_kernel_is_a_posterior_only_sequential_closure_positive_control():
    # Prior weights are deliberately unequal. Exact profile costs add without
    # requiring the mass law to retain a calibration state.
    events = [_event(2., 1.), _event(1., 2.), _event(-.5, .8)]
    prior = FiniteJointMassLaw([[.2], [1.], [3.]], [.1, .3, .6])
    batch = gaussian_joint_posterior(events, prior)
    for order in itertools.permutations(range(3)):
        sequential = prior
        for i in order:
            sequential = gaussian_joint_posterior([events[i]], sequential)
        np.testing.assert_allclose(sequential.probabilities, batch.probabilities,
                                   rtol=2e-14, atol=1e-16)


@pytest.mark.parametrize("partition", [((0,), (1, 2)), ((0, 1), (2,)), ((2,), (1,), (0,))])
def test_additive_retained_state_reconstructs_the_union_law(partition):
    events = [_event(2., 1.), _event(1., 2.), _event(-.5, .8)]
    prior = FiniteJointMassLaw([[.2], [1.], [3.]], [.1, .3, .6])
    prefactor = 1.3
    batch, state = joint_trace_anchor(events, prior, anchor_prefactor=prefactor, **RANK)
    grouped_traces = [sum(events[i].data@events[i].data for i in group) for group in partition]
    trace = math.fsum(grouped_traces)
    assert trace == pytest.approx(state["observed_whitened_trace"], abs=2e-15)
    weight = prefactor/math.sqrt(1+trace)
    log_factors = []
    for mass in prior.atoms:
        groups = [[events[i] for i in group] for group in partition]
        fields = [gaussian_joint_score_statistic(group, mass, **RANK) for group in groups]
        score = sum(field["score"][0] for field in fields)
        covariance = sum(field["covariance"][0, 0] for field in fields)
        cost = math.fsum(gaussian_profile_cost(event, mass) for group in groups for event in group)
        statistic = score*score/covariance
        log_factors.append(-statistic/2-weight*cost/2)
    log_weights = np.log(prior.probabilities)+log_factors
    reconstructed = np.exp(log_weights-special.logsumexp(log_weights))
    np.testing.assert_allclose(reconstructed, batch.probabilities, rtol=2e-14, atol=1e-16)
