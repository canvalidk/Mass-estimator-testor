"""Proper finite trace-anchor laws, event frames and matrix-local references."""
import math

import numpy as np
import pytest

from met.joint_calibration import gaussian_joint_score_moments, joint_covariance_benchmark
from met.joint_mass import (FiniteJointMassLaw, common_acceleration,
                            gaussian_profile_cost, scalar_composite_law,
                            separate_accelerations)
from met.joint_score_law import gaussian_joint_score_statistic
from met.joint_trace_anchor import anchored_local_precision, joint_trace_anchor
from met.score_law_obstruction import ScalarScoreQuadrature, witness_constants


RANK = dict(rank_rtol=1e-11, rank_atol=0.)


def _covariance(size):
    factor = np.eye(size)+np.arange(size*size).reshape(size, size)/(7*size*size)
    return factor@factor.T


def _events():
    return [common_acceleration([[2., .3], [4., -.1]], [1., .2], _covariance(6)),
            common_acceleration([[-.2, 1.], [.4, .8]], [.3, .5], np.diag([1., 3., 2., 5., .7, 1.1]))]


def _prior():
    return FiniteJointMassLaw([[1., 1.], [1., 3.], [2., 2.], [3., 1.]], [.1, .2, .3, .4])


@pytest.mark.parametrize("prefactor", [.3, 1., 4.])
def test_k2_finite_law_is_proper_with_declared_covariances_and_explicit_prefactor(prefactor):
    events, prior = _events(), _prior()
    law, state = joint_trace_anchor(events, prior, anchor_prefactor=prefactor, **RANK)
    expected_trace = sum(e.data@np.linalg.solve(e.covariance, e.data) for e in events)
    assert state["observed_whitened_trace"] == pytest.approx(expected_trace, rel=2e-14)
    assert state["anchor_weight"] == pytest.approx(prefactor/math.sqrt(1+expected_trace), rel=2e-14)
    assert state["body_count"] == 2 and state["spatial_dimension"] == 2
    assert np.all(law.probabilities > 0)
    assert law.probabilities.sum() == pytest.approx(1., abs=1e-15)
    assert math.isfinite(law.log_evidence) and law.log_evidence <= 1e-14
    np.testing.assert_array_equal(law.atoms, prior.atoms)


def test_k1_point_factor_ratios_match_the_continuous_scalar_reference():
    case = witness_constants()
    scatter = 31*case["population_scatter"]
    data = np.linalg.cholesky(scatter).T
    events = [common_acceleration([[f]], [a], case["covariance"]) for f, a in data]
    prior = FiniteJointMassLaw([[.2], [.4], [1.], [1.5], [3.]], [.1, .2, .25, .3, .15])
    law, state = joint_trace_anchor(events, prior, anchor_prefactor=.7, **RANK)
    q = ScalarScoreQuadrature(scatter, case["covariance"], prior_center=1.,
        anchor_exponent=.5, anchor_schedule="spectral", anchor_weight_override=state["anchor_weight"],
        epsabs=1e-11, epsrel=1e-10, limit=500)
    factors = np.array([q.kernel(math.atan(m), "anchored_score") for m in prior.atoms[:, 0]])
    expected = prior.probabilities*factors
    expected /= expected.sum()
    # This checks likelihood ratios at supplied atoms, not continuous-prior
    # quadrature accuracy or an atomic-prior TV Gaussian limit.
    np.testing.assert_allclose(law.probabilities, expected, atol=3e-14)


def test_independent_event_frame_rotation_preserves_trace_and_law_while_pooled_eigenvalue_changes():
    first = common_acceleration([[2., 0.]], [1., 0.], np.eye(4))
    rotated = common_acceleration([[0., 2.]], [0., 1.], np.eye(4))
    prior = FiniteJointMassLaw([[1.], [2.], [3.]], [.2, .5, .3])
    same_frame = [first, first]
    independent_frame = [first, rotated]
    original, state = joint_trace_anchor(same_frame, prior, anchor_prefactor=1., **RANK)
    changed, transformed = joint_trace_anchor(independent_frame, prior, anchor_prefactor=1., **RANK)
    assert state["observed_whitened_trace"] == transformed["observed_whitened_trace"] == 10.
    assert state["anchor_weight"] == transformed["anchor_weight"]
    np.testing.assert_allclose(original.probabilities, changed.probabilities, atol=1e-16)
    pooled_first = sum(np.outer(e.data, e.data) for e in same_frame)
    pooled_changed = sum(np.outer(e.data, e.data) for e in independent_frame)
    assert np.linalg.eigvalsh(pooled_first)[-1] == pytest.approx(10.)
    assert np.linalg.eigvalsh(pooled_changed)[-1] == pytest.approx(5.)


def test_per_event_frames_and_units_preserve_the_finite_joint_law_and_state():
    events, prior = _events(), _prior()
    original, state = joint_trace_anchor(events, prior, anchor_prefactor=1.3, **RANK)
    factors = np.array([3., 7., .4])
    changed = []
    for event, angle in zip(events, (.3, -.7)):
        rotation = np.array([[math.cos(angle), -math.sin(angle)],
                             [math.sin(angle), math.cos(angle)]])
        transform = np.zeros((6, 6))
        for j, factor in enumerate(factors):
            transform[2*j:2*j+2, 2*j:2*j+2] = factor*rotation
        data = transform@event.data
        changed.append(common_acceleration(data[:4].reshape(2, 2), data[4:],
                                            transform@event.covariance@transform.T))
    jacobian = factors[:2]/factors[2]
    transformed_prior = FiniteJointMassLaw(prior.atoms*jacobian, prior.probabilities)
    converted, converted_state = joint_trace_anchor(changed, transformed_prior,
                                                    anchor_prefactor=1.3, **RANK)
    assert converted_state["observed_whitened_trace"] == pytest.approx(state["observed_whitened_trace"], rel=2e-14)
    assert converted_state["anchor_weight"] == pytest.approx(state["anchor_weight"], rel=2e-14)
    np.testing.assert_allclose(converted.probabilities, original.probabilities, atol=2e-14)


def test_exact_zero_padding_empty_law_and_supplied_prior_coupling():
    events, prior = _events(), _prior()
    original, state = joint_trace_anchor(events, prior, anchor_prefactor=.9, **RANK)
    zero = common_acceleration(np.zeros((2, 2)), np.zeros(2), 3*np.eye(6))
    padded, padded_state = joint_trace_anchor(events+[zero], prior, anchor_prefactor=.9, **RANK)
    np.testing.assert_array_equal(padded.probabilities, original.probabilities)
    assert padded.log_evidence == original.log_evidence
    assert padded_state == state
    for series in ([], [zero]):
        unchanged, observed = joint_trace_anchor(series, prior, anchor_prefactor=.9, **RANK)
        np.testing.assert_allclose(unchanged.probabilities, prior.probabilities, atol=1e-16)
        assert observed["observed_whitened_trace"] == 0 and observed["anchor_weight"] == .9
    coupled = FiniteJointMassLaw([[1., 1.], [2., 2.]], [.5, .5])
    independent = FiniteJointMassLaw([[1., 1.], [1., 2.], [2., 1.], [2., 2.]], [.25]*4)
    first, _ = joint_trace_anchor([], coupled, anchor_prefactor=1., **RANK)
    second, _ = joint_trace_anchor([], independent, anchor_prefactor=1., **RANK)
    assert len(scalar_composite_law(first, operation="addition").support) == 2
    assert len(scalar_composite_law(second, operation="addition").support) == 3


def test_separate_large_common_contrasts_do_not_erase_the_joint_prior():
    event = common_acceleration([[1e100], [-1e100]], [0.], np.eye(3))
    prior = FiniteJointMassLaw([[1., 1.], [2., 2.], [3., 3.]], [.1, .2, .7])
    law, state = joint_trace_anchor([event], prior, anchor_prefactor=1., **RANK)
    np.testing.assert_allclose(law.probabilities, prior.probabilities, atol=1e-16)
    assert state["minimum_profile_cost"] == 2e200
    assert state["minimum_score_statistic"] == 0.
    assert math.isfinite(law.log_evidence) and law.log_evidence < -1e99


def test_large_common_combined_contrast_is_shifted_before_adding_prior_weights():
    covariance = np.array([[1., .9], [.9, 1.]])
    chol = np.linalg.cholesky(covariance)
    true_line = np.linalg.solve(chol, [1.2, 1.])
    true_line /= np.linalg.norm(true_line)
    false_line = np.linalg.solve(chol, [4/15, 1.])
    false_line /= np.linalg.norm(false_line)
    halfway = chol@((true_line+false_line)/math.sqrt(2))
    atoms = np.array([[halfway[0]/halfway[1]], [4/15]])
    radius = 1e10
    data = np.column_stack([math.sqrt(2)*radius*(chol@true_line),
                            radius*(chol@false_line)]).T
    events = [common_acceleration([[y[0]]], [y[1]], covariance) for y in data]
    statistics = np.array([gaussian_joint_score_statistic(events, atom, **RANK)["statistic"]
                           for atom in atoms])
    costs = np.array([sum(gaussian_profile_cost(e, atom) for e in events) for atom in atoms])
    trace = math.fsum(float(np.linalg.solve(chol, e.data)@np.linalg.solve(chol, e.data))
                      for e in events)
    weight = (statistics[0]-statistics[1])/(costs[1]-costs[0])
    prefactor = weight*math.sqrt(1+trace)
    for candidate in (prefactor, np.nextafter(prefactor, 0.), np.nextafter(prefactor, math.inf)):
        effective = candidate/math.sqrt(1+trace)
        combined = (.5*(statistics-statistics.min())
                    +.5*effective*(costs-costs.min()))
        if combined[0] == combined[1]:
            prefactor = candidate
            break
    else:
        pytest.fail("physical regression did not produce its tied represented contrasts")
    assert combined[0] > 1e18
    # This fixture checks normalization of equal represented contrasts, not
    # arbitrary-precision geometry or a proposed data-dependent kappa rule.
    prior = FiniteJointMassLaw(atoms, [.1, .9])
    law, state = joint_trace_anchor(events, prior, anchor_prefactor=prefactor, **RANK)
    np.testing.assert_allclose(law.probabilities, prior.probabilities, atol=1e-16)
    assert state["minimum_combined_contrast"] == combined[0]
    expected_log_evidence = -(combined[0]+statistics.min()/2+effective*costs.min()/2)
    assert law.log_evidence == pytest.approx(expected_log_evidence, rel=2e-16)


def test_genuine_matrix_benchmark_and_weak_path_anchor_distortion():
    masses = np.array([1., 2.])
    covariance = np.diag([1., 3., 2., 5., .7, 1.1])
    event = common_acceleration(np.zeros((2, 2)), np.zeros(2), covariance)
    moments = gaussian_joint_score_moments([event], masses, [[1., 0.]])
    a1, geometry = moments["sensitivity"], moments["noise_geometry"]
    first = np.linalg.inv(np.diag([1., 2.])+.7*np.outer(masses, masses))
    second = np.linalg.inv(np.diag([3., 5.])+1.1*np.outer(masses, masses))
    gram_first = 1+4/2+1/.7
    gram_second = 1/3+4/5+1/1.1
    np.testing.assert_allclose(a1, first, atol=2e-15)
    np.testing.assert_allclose(geometry, first/gram_first+second/gram_second, atol=2e-15)
    distortions = []
    for count in (1e4, 1e8, 1e12):
        energy = count**.6
        a = energy*a1
        b = a+count*geometry
        # This trace is a stated Gaussian-experiment expectation used only
        # in the theoretical A/B comparison, never a law-update input.
        weight = 1/math.sqrt(1+6*count+gram_first*energy)
        reference = joint_covariance_benchmark(a, b)
        result = anchored_local_precision(a, b, anchor_weight=weight)
        np.testing.assert_allclose(result["score_precision"],
            np.linalg.inv(reference["sampling_covariance"]), rtol=3e-14)
        np.testing.assert_allclose(result["total_precision"],
            result["score_precision"]+weight*a, rtol=3e-14)
        eigenvalues = result["relative_anchor_eigenvalues"]
        assert eigenvalues[1]-eigenvalues[0] > .001*eigenvalues[0]
        distortions.append(eigenvalues[-1])
    assert distortions[0] > distortions[1] > distortions[2]
    assert distortions[-1] < .03
    unanchored = anchored_local_precision(a1, a1+geometry, anchor_weight=0.)
    expected = joint_covariance_benchmark(a1, a1+geometry)["sampling_covariance"]
    np.testing.assert_allclose(unanchored["covariance"], expected, rtol=3e-14)


def test_bounded_protocol_and_explicit_control_validation():
    prior = _prior()
    separate = separate_accelerations([[1.], [2.]], [[.5], [.7]], np.eye(4))
    permuted = common_acceleration([[1.], [2.]], [.5], np.eye(3), body_indices=[1, 0])
    for event in (separate, permuted):
        with pytest.raises(ValueError, match="canonical"):
            joint_trace_anchor([event], prior, anchor_prefactor=1., **RANK)
    for prefactor in (0., -1., math.nan, math.inf):
        with pytest.raises(ValueError):
            joint_trace_anchor([], prior, anchor_prefactor=prefactor, **RANK)
    with pytest.raises(ValueError):
        joint_trace_anchor([], prior, anchor_prefactor=1., rank_rtol=-1., rank_atol=0.)
    with pytest.raises(ValueError):
        anchored_local_precision(np.zeros((2, 2)), np.eye(2), anchor_weight=1.)
