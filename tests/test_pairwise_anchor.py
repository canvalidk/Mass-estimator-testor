"""T67: independent physical GLS comparisons and finite-law identities."""
import numpy as np
import pytest

from met.joint_mass import FiniteJointMassLaw, common_acceleration
from met.pairwise_anchor import joint_pairwise_anchor, pairwise_profile_penalty

RANK = dict(rank_rtol=1e-12, rank_atol=0.)


def _case():
    prior = FiniteJointMassLaw([[.4, .8], [1.2, .8], [1.3, 2.], [3., 1.]], [.1, .2, .3, .4])
    events = [common_acceleration([[2.1], [.4]], [.7], np.diag([.8, 2., 1.1])),
              common_acceleration([[-.3], [1.4]], [-.2], np.array([[2., .2, .1], [.2, 1., -.1], [.1, -.1, .7]]))]
    return events, prior


def _fit(event, masses):
    c = np.r_[masses, 1.][:, None]
    precision = np.linalg.inv(event.covariance)
    mean = c @ np.linalg.solve(c.T @ precision @ c, c.T @ precision @ event.data)
    residual = event.data-mean
    return mean, residual @ precision @ residual


@pytest.mark.parametrize("prefactor", [.2, 1., 3.])
def test_pairwise_penalty_and_law_match_independent_physical_residuals(prefactor):
    events, prior = _case()
    state = pairwise_profile_penalty(events, prior)
    gain, difference = np.zeros((4, 4)), np.zeros((4, 4))
    score, covariance = np.zeros((4, 2)), np.zeros((4, 2, 2))
    for event in events:
        fits = [_fit(event, masses) for masses in prior.atoms]
        for i, (mean, cost) in enumerate(fits):
            for j, (other, other_cost) in enumerate(fits):
                delta = other-mean
                gain[i, j] += (cost-other_cost)/2
                difference[i, j] += delta @ np.linalg.solve(event.covariance, delta)
            masses = prior.atoms[i]
            constraint = np.column_stack([np.eye(2), -masses])
            residual = constraint @ event.data
            inverse = np.linalg.inv(constraint @ event.covariance @ constraint.T)
            latent = mean[-1]
            score[i] += latent*inverse @ residual
            covariance[i] += latent**2*inverse
    expected = (np.maximum(gain, 0.)/np.sqrt(1+difference)) @ prior.probabilities
    np.testing.assert_allclose(state["profile_contrast"], gain, atol=2e-14)
    np.testing.assert_allclose(state["projection_difference"], difference, atol=2e-14)
    np.testing.assert_allclose(state["penalty"], expected, atol=2e-14)
    exponent = np.array([-.5*s @ np.linalg.solve(c, s) for s, c in zip(score, covariance)])-prefactor*expected
    weights = prior.probabilities*np.exp(exponent-exponent.max())
    law, _ = joint_pairwise_anchor(events, prior, anchor_prefactor=prefactor, **RANK)
    np.testing.assert_allclose(law.probabilities, weights/weights.sum(), atol=2e-14)


def test_pairwise_fields_add_before_penalty_and_obey_energy_bounds():
    events, prior = _case()
    state = pairwise_profile_penalty(events, prior)
    parts = [pairwise_profile_penalty([event], prior) for event in events]
    for name in ("profile_contrast", "projection_difference"):
        np.testing.assert_allclose(state[name], parts[0][name]+parts[1][name], atol=2e-14)
    energy = sum(event.data @ np.linalg.solve(event.covariance, event.data) for event in events)
    assert np.max(state["projection_difference"]) <= energy
    assert np.max(state["penalty"]) <= energy/2
    np.testing.assert_allclose(state["profile_contrast"], -state["profile_contrast"].T, atol=2e-14)
    np.testing.assert_allclose(state["projection_difference"], state["projection_difference"].T, atol=2e-14)


def test_zero_events_prior_coupling_and_passive_units():
    events, prior = _case()
    law, state = joint_pairwise_anchor(events, prior, anchor_prefactor=1., **RANK)
    zero = common_acceleration([[0.], [0.]], [0.], np.eye(3))
    padded, _ = joint_pairwise_anchor([zero]+events, prior, anchor_prefactor=1., **RANK)
    np.testing.assert_array_equal(law.probabilities, padded.probabilities)
    empty, _ = joint_pairwise_anchor([], prior, anchor_prefactor=1., **RANK)
    np.testing.assert_allclose(empty.probabilities, prior.probabilities, atol=2e-16)
    # Force units x3 and acceleration units x2 imply masses x1.5.
    transform = np.diag([3., 3., 2.])
    changed = [common_acceleration((transform@e.data)[:2, None], (transform@e.data)[2:],
                                   transform@e.covariance@transform.T) for e in events]
    scaled_prior = FiniteJointMassLaw(prior.atoms*1.5, prior.probabilities)
    scaled, other = joint_pairwise_anchor(changed, scaled_prior, anchor_prefactor=1., **RANK)
    np.testing.assert_allclose(scaled.probabilities, law.probabilities, atol=2e-14)
    np.testing.assert_allclose(other["penalty"], state["penalty"], atol=2e-14)
    weights = prior.probabilities[::-1]
    other = pairwise_profile_penalty(events, FiniteJointMassLaw(prior.atoms, weights))
    assert not np.allclose(other["penalty"], state["penalty"])


@pytest.mark.parametrize("bad", [0., -1., np.inf, np.nan])
def test_prefactor_validation(bad):
    events, prior = _case()
    with pytest.raises(ValueError, match="prefactor"):
        joint_pairwise_anchor(events, prior, anchor_prefactor=bad, **RANK)


def test_overflowing_finite_prefactor_is_an_explicit_numerical_error():
    events, prior = _case()
    exaggerated = [common_acceleration([[100.], [5.]], [1.], np.eye(3))]+events
    with pytest.raises(ValueError, match="floating-point"):
        joint_pairwise_anchor(exaggerated, prior, anchor_prefactor=1e308, **RANK)


def test_independent_event_frames_signs_and_scalar_channel_swap():
    prior = FiniteJointMassLaw([[.4], [1.2], [3.]], [.2, .5, .3])
    event = common_acceleration([[2., -.3]], [.8, .2], np.diag([1., 2., .7, 1.1]))
    other = common_acceleration([[-.4, 1.]], [.3, .7], np.diag([.8, 1., 2., .9]))
    expected, _ = joint_pairwise_anchor([event, other], prior, anchor_prefactor=1., **RANK)
    rotation = np.array([[.6, -.8], [.8, .6]])
    transform = np.kron(np.eye(2), rotation)
    y = transform @ other.data
    rotated = common_acceleration([-y[:2]], -y[2:], transform@other.covariance@transform.T)
    actual, _ = joint_pairwise_anchor([event, rotated], prior, anchor_prefactor=1., **RANK)
    np.testing.assert_allclose(actual.probabilities, expected.probabilities, atol=2e-14)
    scalar = common_acceleration([[2.]], [.7], np.array([[2., .4], [.4, 1.]]))
    swapped = common_acceleration([[.7]], [2.], scalar.covariance[::-1, ::-1])
    original, _ = joint_pairwise_anchor([scalar], prior, anchor_prefactor=1., **RANK)
    reciprocal, _ = joint_pairwise_anchor([swapped],
        FiniteJointMassLaw(1/prior.atoms, prior.probabilities), anchor_prefactor=1., **RANK)
    # FiniteJointMassLaw sorts its atoms; reciprocal exchange reverses them.
    np.testing.assert_allclose(reciprocal.atoms[:, 0], 1/original.atoms[::-1, 0])
    np.testing.assert_allclose(reciprocal.probabilities, original.probabilities[::-1], atol=2e-14)


def test_common_large_score_cost_does_not_erase_pairwise_ordering(monkeypatch):
    import met.pairwise_anchor as module
    events, prior = _case()
    expected = pairwise_profile_penalty(events, prior)["penalty"]
    monkeypatch.setattr(module, "gaussian_joint_score_statistic",
                        lambda *args, **kwargs: {"statistic": 1e200})
    law, _ = joint_pairwise_anchor(events, prior, anchor_prefactor=1., **RANK)
    weights = prior.probabilities*np.exp(-expected)
    np.testing.assert_allclose(law.probabilities, weights/weights.sum(), atol=2e-14)
