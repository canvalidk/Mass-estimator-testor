"""Prior-projector variance checked with independent physical GLS geometry."""
import math

import numpy as np
import pytest
from scipy import special

from met.joint_mass import (FiniteJointMassLaw, common_acceleration,
                            gaussian_profile_cost, separate_accelerations)
from met.joint_score_law import gaussian_joint_score_statistic
from met.mass_contrast_anchor import (joint_mass_contrast_anchor,
                                      observed_mass_contrast_energy,
                                      scalar_prior_projector_contrast)


RANK = dict(rank_rtol=1e-11, rank_atol=0.)


def _covariance(size):
    factor = np.eye(size)+np.arange(size*size).reshape(size, size)/(5*size*size)
    return factor@factor.T


def _prior():
    return FiniteJointMassLaw([[.4, .8], [1.2, .8], [1.3, 2.], [3., 1.]], [.1, .2, .3, .4])


def _events(dimension=1):
    if dimension == 1:
        return [common_acceleration([[2.1], [.4]], [.7], _covariance(3)),
                common_acceleration([[-.3], [1.4]], [-.2], np.diag([.7, 2.3, 1.1]))]
    return [common_acceleration([[2., .3], [4., -.1]], [1., .2], _covariance(6)),
            common_acceleration([[-.2, 1.], [.4, .8]], [.3, .5], np.diag([1., 3., 2., 5., .7, 1.1]))]


def _physical_fit(event, masses):
    design = np.vstack([mass*np.eye(event.dimension) for mass in masses]+[np.eye(event.dimension)])
    precision_design = np.linalg.solve(event.covariance, design)
    fitted = np.linalg.solve(design.T@precision_design, precision_design.T@event.data)
    return design@fitted


def _physical_energy(events, prior):
    total = 0.
    for event in events:
        fits = np.array([_physical_fit(event, masses) for masses in prior.atoms])
        differences = fits-prior.probabilities@fits
        total += sum(p*(difference@np.linalg.solve(event.covariance, difference))
                     for p, difference in zip(prior.probabilities, differences))
    return total


def _physical_scalar_kernel(events, masses, weight):
    # K scalar body masses, d=1: physical residual F-m*A has covariance V.
    score, covariance, cost = np.zeros(len(masses)), np.zeros((len(masses), len(masses))), 0.
    for event in events:
        channel = np.r_[masses, 1.]
        constraint = np.column_stack([np.eye(len(masses)), -masses])
        residual_covariance = constraint@event.covariance@constraint.T
        precision = np.linalg.inv(residual_covariance)
        residual = constraint@event.data
        precision_channel = np.linalg.solve(event.covariance, channel)
        fitted = precision_channel@event.data/(channel@precision_channel)
        score += fitted*precision@residual
        covariance += fitted*fitted*precision
        cost += residual@precision@residual
    return -.5*(score@np.linalg.solve(covariance, score)+weight*cost)


@pytest.mark.parametrize("prefactor", [.3, 1., 4.])
def test_finite_joint_kernel_matches_physical_residuals_and_prior_projection_variance(prefactor):
    events, prior = _events(), _prior()
    law, state = joint_mass_contrast_anchor(events, prior, anchor_prefactor=prefactor, **RANK)
    expected_u = _physical_energy(events, prior)
    assert state["observed_mass_contrast_energy"] == pytest.approx(expected_u, rel=3e-14)
    weight = prefactor/math.sqrt(1+expected_u)
    assert state["anchor_weight"] == pytest.approx(weight, rel=2e-14)
    log_factors = np.array([_physical_scalar_kernel(events, masses, weight) for masses in prior.atoms])
    log_weights = np.log(prior.probabilities)+log_factors
    np.testing.assert_allclose(law.probabilities, np.exp(log_weights-special.logsumexp(log_weights)),
                               rtol=4e-14, atol=3e-15)
    assert law.log_evidence == pytest.approx(special.logsumexp(log_weights), abs=3e-14)
    assert np.all(law.probabilities > 0) and law.probabilities.sum() == pytest.approx(1.)
    total_energy = sum(event.data@np.linalg.solve(event.covariance, event.data) for event in events)
    assert 0 < expected_u < total_energy/4


def test_centered_projection_variance_equals_the_operator_average_formula():
    prior = _prior()
    for event in _events(2):
        chol = np.linalg.cholesky(event.covariance)
        whitened = np.linalg.solve(chol, event.data)
        average = np.zeros((6, 6))
        for mass, probability in zip(prior.atoms, prior.probabilities):
            design = np.vstack([m*np.eye(2) for m in mass]+[np.eye(2)])
            basis = np.linalg.solve(chol, design)
            projector = basis@np.linalg.solve(basis.T@basis, basis.T)
            average += probability*projector
        operator = average-average@average
        np.testing.assert_allclose(operator, operator.T, atol=3e-16)
        assert np.linalg.eigvalsh(operator)[0] > -2e-16
        assert np.linalg.eigvalsh(operator)[-1] < .25+2e-16
        assert observed_mass_contrast_energy([event], prior) == pytest.approx(
            whitened@operator@whitened, rel=4e-14, abs=3e-15)


def test_independent_event_frames_units_signs_and_prior_relabeling_preserve_the_law():
    events, prior = _events(2), _prior()
    law, state = joint_mass_contrast_anchor(events, prior, anchor_prefactor=.8, **RANK)
    converted_events = []
    factors = np.array([3., 7., .4])
    for event, angle, sign in zip(events, (.3, -.7), (1., -1.)):
        rotation = sign*np.array([[math.cos(angle), -math.sin(angle)],
                                  [math.sin(angle), math.cos(angle)]])
        transform = np.kron(np.diag(factors), rotation)
        point = transform@event.data
        converted_events.append(common_acceleration(point[:4].reshape(2, 2), point[4:],
                                                    transform@event.covariance@transform.T))
    converted_prior = FiniteJointMassLaw(prior.atoms*factors[:2]/factors[2], prior.probabilities)
    converted, changed = joint_mass_contrast_anchor(converted_events, converted_prior,
                                                   anchor_prefactor=.8, **RANK)
    assert changed["observed_mass_contrast_energy"] == pytest.approx(
        state["observed_mass_contrast_energy"], rel=5e-14)
    assert changed["anchor_weight"] == pytest.approx(state["anchor_weight"], rel=3e-14)
    np.testing.assert_allclose(converted.probabilities, law.probabilities, rtol=4e-14, atol=4e-15)
    reversed_prior = FiniteJointMassLaw(prior.atoms[::-1], prior.probabilities[::-1])
    reordered, reordered_state = joint_mass_contrast_anchor(events[::-1], reversed_prior,
                                                           anchor_prefactor=.8, **RANK)
    np.testing.assert_allclose(reordered.probabilities, law.probabilities, atol=4e-15)
    assert reordered_state["observed_mass_contrast_energy"] == pytest.approx(
        state["observed_mass_contrast_energy"], rel=3e-15)


def test_exact_zero_padding_single_atom_and_huge_equal_projections_preserve_the_supplied_prior():
    events, prior = _events(), _prior()
    law, state = joint_mass_contrast_anchor(events, prior, anchor_prefactor=.8, **RANK)
    zero = common_acceleration([[0.], [0.]], [0.], 3*np.eye(3))
    padded, padded_state = joint_mass_contrast_anchor(events+[zero], prior, anchor_prefactor=.8, **RANK)
    np.testing.assert_array_equal(padded.probabilities, law.probabilities)
    assert padded_state == state
    for series in ([], [zero]):
        unchanged, observed = joint_mass_contrast_anchor(series, prior, anchor_prefactor=.8, **RANK)
        np.testing.assert_allclose(unchanged.probabilities, prior.probabilities, atol=1e-16)
        assert observed["observed_mass_contrast_energy"] == 0.
    singleton = FiniteJointMassLaw([[1., 2.]], [1.])
    assert observed_mass_contrast_energy(events, singleton) == 0.
    huge = common_acceleration([[1e100], [-1e100]], [0.], np.eye(3))
    equal_fits = FiniteJointMassLaw([[1., 1.], [2., 2.], [3., 3.]], [.1, .2, .7])
    unchanged, observed = joint_mass_contrast_anchor([huge], equal_fits, anchor_prefactor=1., **RANK)
    assert observed["observed_mass_contrast_energy"] == 0.
    np.testing.assert_allclose(unchanged.probabilities, equal_fits.probabilities, atol=1e-16)
    assert math.isfinite(unchanged.log_evidence)


def test_joint_prior_coupling_changes_the_projector_variance_without_inventing_independence():
    event = common_acceleration([[2.], [.3]], [1.], np.eye(3))
    coupled = FiniteJointMassLaw([[1., 1.], [3., 3.]], [.5, .5])
    anticoupled = FiniteJointMassLaw([[1., 3.], [3., 1.]], [.5, .5])
    # Both body marginals are the same equal law on {1,3}.
    for column in (0, 1):
        np.testing.assert_array_equal(np.sort(coupled.atoms[:, column]),
                                      np.sort(anticoupled.atoms[:, column]))
    u1 = observed_mass_contrast_energy([event], coupled)
    u2 = observed_mass_contrast_energy([event], anticoupled)
    assert abs(u1-u2) > .01
    for prior, expected in ((coupled, u1), (anticoupled, u2)):
        law, state = joint_mass_contrast_anchor([event], prior, anchor_prefactor=.7, **RANK)
        np.testing.assert_array_equal(law.atoms, prior.atoms)
        assert state["observed_mass_contrast_energy"] == expected


def test_weak_prior_limit_is_continuous_at_a_single_graph_without_density_or_fit():
    event = common_acceleration([[2.]], [1.], np.eye(2))
    prior = FiniteJointMassLaw([[1.], [3.]], [.5, .5])
    distance_squared = 4*observed_mass_contrast_energy([event], prior)
    previous = math.inf
    for epsilon in (.1, .01, .001, 1e-5):
        approaching = FiniteJointMassLaw([[1.], [3.]], [1-epsilon, epsilon])
        law, state = joint_mass_contrast_anchor([event], approaching, anchor_prefactor=.7, **RANK)
        assert state["observed_mass_contrast_energy"] == pytest.approx(
            epsilon*(1-epsilon)*distance_squared, rel=3e-14)
        assert state["observed_mass_contrast_energy"] < previous
        previous = state["observed_mass_contrast_energy"]
        assert law.probabilities[1] < 3*epsilon
    assert state["anchor_weight"] == pytest.approx(.7, abs=1e-5)


def test_combined_large_contrasts_are_shifted_before_adding_unequal_prior_weights():
    covariance = np.array([[1., .9], [.9, 1.]])
    chol = np.linalg.cholesky(covariance)
    true_line = np.linalg.solve(chol, [1.2, 1.])
    true_line /= np.linalg.norm(true_line)
    false_line = np.linalg.solve(chol, [4/15, 1.])
    false_line /= np.linalg.norm(false_line)
    halfway = chol@((true_line+false_line)/math.sqrt(2))
    atoms = np.array([[halfway[0]/halfway[1]], [4/15]])
    radius = 1e10
    data = np.column_stack([math.sqrt(2)*radius*(chol@true_line), radius*(chol@false_line)]).T
    events = [common_acceleration([[point[0]]], [point[1]], covariance) for point in data]
    prior = FiniteJointMassLaw(atoms, [.1, .9])
    # Use the prior's canonical order throughout this floating-point fixture.
    statistics = np.array([gaussian_joint_score_statistic(events, atom, **RANK)["statistic"]
                           for atom in prior.atoms])
    costs = np.array([sum(gaussian_profile_cost(event, atom) for event in events)
                      for atom in prior.atoms])
    weight = (statistics[0]-statistics[1])/(costs[1]-costs[0])
    u = observed_mass_contrast_energy(events, prior)
    prefactor = weight*math.sqrt(1+u)
    for candidate in (prefactor, np.nextafter(prefactor, 0.), np.nextafter(prefactor, math.inf)):
        effective = candidate/math.sqrt(1+u)
        combined = .5*(statistics-statistics.min())+.5*effective*(costs-costs.min())
        if combined[0] == combined[1]:
            prefactor = candidate
            break
    else:
        pytest.fail("physical regression did not produce tied represented contrasts")
    assert combined[0] > 1e18
    law, state = joint_mass_contrast_anchor(events, prior, anchor_prefactor=prefactor, **RANK)
    np.testing.assert_allclose(law.probabilities, prior.probabilities, atol=1e-16)
    assert state["minimum_combined_contrast"] == combined[0]


def test_explicit_protocol_and_numerical_controls_are_checked():
    prior = _prior()
    separate = separate_accelerations([[1.], [2.]], [[.5], [.7]], np.eye(4))
    for control in (0., -1., math.nan, math.inf):
        with pytest.raises(ValueError):
            joint_mass_contrast_anchor([], prior, anchor_prefactor=control, **RANK)
    with pytest.raises(ValueError, match="canonical"):
        joint_mass_contrast_anchor([separate], prior, anchor_prefactor=1., **RANK)
    with pytest.raises(ValueError):
        joint_mass_contrast_anchor([], prior, anchor_prefactor=1., rank_rtol=-1., rank_atol=0.)
    with pytest.raises(ValueError, match="FiniteJointMassLaw"):
        observed_mass_contrast_energy([], None)


@pytest.mark.parametrize("ratio", [1., 3., 1000.])
def test_continuous_scalar_coefficient_matches_independent_half_cauchy_integrals(ratio):
    result = scalar_prior_projector_contrast(np.diag([ratio*ratio, 1.]),
        prior_center=1., mass_bounds=(0., math.inf), epsabs=1e-12, epsrel=1e-11, limit=500)
    off_diagonal = (1/math.pi if ratio == 1
                    else 2*ratio*math.log(ratio)/(math.pi*(ratio*ratio-1)))
    expected = ratio/(ratio+1)**2-off_diagonal**2
    assert result["coefficient"] == pytest.approx(expected, rel=2e-10, abs=2e-13)
    assert 0 < result["coefficient"] <= .25


def test_compact_prior_coefficient_retains_tiny_mass_contrast_in_degenerate_geometry():
    from scipy import integrate
    lower, upper = .1, 10.
    normalization = math.atan(upper)-math.atan(lower)
    mean = integrate.quad(lambda m: m/(1+m*m), lower, upper)[0]/normalization
    second = integrate.quad(lambda m: m*m/(1+m*m), lower, upper)[0]/normalization
    variance = second-mean*mean
    for ratio in (1e4, 1e10, 1e30, 1e40):
        result = scalar_prior_projector_contrast(np.diag([ratio*ratio, 1.]),
            prior_center=1., mass_bounds=(lower, upper), epsabs=1e-12, epsrel=1e-11, limit=500)
        assert result["coefficient"] > 0
        tolerance = 2e-6 if ratio == 1e4 else 2e-12
        assert result["coefficient"]*ratio*ratio == pytest.approx(variance, rel=tolerance)


def test_continuous_scalar_coefficient_respects_units_and_reciprocal_prior_pushforward():
    covariance = np.array([[2., .4], [.4, .7]])
    controls = dict(epsabs=1e-12, epsrel=1e-11, limit=500)
    original = scalar_prior_projector_contrast(covariance, prior_center=2.,
        mass_bounds=(.1, 8.), **controls)["coefficient"]
    transform = np.diag([3., .4])
    converted = scalar_prior_projector_contrast(transform@covariance@transform,
        prior_center=15., mass_bounds=(.75, 60.), **controls)["coefficient"]
    swapped = scalar_prior_projector_contrast(covariance[::-1, ::-1],
        prior_center=.5, mass_bounds=(1/8, 10.), **controls)["coefficient"]
    assert converted == pytest.approx(original, rel=8e-14)
    assert swapped == pytest.approx(original, rel=8e-14)
    for scale in (1e-200, 1e200):
        scaled = scalar_prior_projector_contrast(scale*covariance, prior_center=2.,
            mass_bounds=(.1, 8.), **controls)["coefficient"]
        assert scaled == pytest.approx(original, rel=8e-14)
