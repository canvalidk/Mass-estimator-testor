"""Independent geometry checks for O16's fixed relative sensor tilt."""
import numpy as np
import pytest

from met.full import Reading, bartlett_eta, isotropic, reading_state
from met.misalignment import isotropic_minima
from met.shared_error_geometry import excitation_spectrum, rotation, score_tilt, shared_tilt_diagnostic, skew


def _readings(accelerations, omega, mass=2.0):
    return [isotropic(mass * a + mass * np.cross(omega, a), a, 0.1, 0.2)
            for a in np.asarray(accelerations)]


def test_skew_and_rodrigues_match_geometric_rotation():
    a, b = np.array([1., 2., 3.]), np.array([-2., 0.5, 1.])
    np.testing.assert_allclose(skew(a) @ b, np.cross(a, b))
    result = rotation([0, 0, 2], np.pi / 2)
    np.testing.assert_allclose(result @ [1, 0, 0], [0, 1, 0], atol=1e-15)
    np.testing.assert_allclose(result.T @ result, np.eye(3), atol=1e-15)
    assert np.linalg.det(result) == pytest.approx(1)


def test_local_fit_recovers_a_common_rotation_and_mass_perturbation():
    acceleration = np.vstack((np.eye(3), -np.eye(3), [[1, 2, -1]]))
    omega, mass, delta = np.array([0.002, -0.004, 0.003]), 2.0, 0.01
    readings = [isotropic((mass + delta) * a + mass * np.cross(omega, a), a, .1, .2)
                for a in acceleration]
    result = shared_tilt_diagnostic(readings, mass)
    assert result["delta_mass"] == pytest.approx(delta, abs=1e-14)
    np.testing.assert_allclose(result["rotation_vector"], omega, atol=1e-14)
    assert result["rotation_rank"] == 3 and result["fit_rank"] == 4
    assert result["fitted_cost"] < 1e-24


def test_exact_small_rotation_has_second_order_linearization_error():
    axis = np.array([1., 2., -1.]) / np.sqrt(6)
    errors = []
    for angle in (0.02, 0.01):
        matrix = rotation(axis, angle)
        readings = [isotropic(2 * matrix @ a, a, .1, .2) for a in np.eye(3)]
        result = shared_tilt_diagnostic(readings, 2)
        errors.append(np.linalg.norm(result["residual"] - result["predicted_residual"]))
        np.testing.assert_allclose(result["rotation_vector"], angle * axis, atol=angle ** 3)
    assert errors[0] / errors[1] == pytest.approx(4, rel=2e-3)


def test_rank_one_excitation_detects_off_axis_tilt_but_not_axis_rotation():
    accelerations = np.array([[0., 0., 1.], [0., 0., 2.], [0., 0., -1.]])
    omega = np.array([0.02, -0.01, 0.3])
    result = shared_tilt_diagnostic(_readings(accelerations, omega), 2)
    assert result["rotation_rank"] == 2
    np.testing.assert_allclose(result["rotation_vector"], omega * [1, 1, 0], atol=1e-14)
    np.testing.assert_allclose(result["rotation_projector"], np.diag([1, 1, 0]), atol=1e-14)
    matrix = rotation([0, 0, 1], 0.7)
    np.testing.assert_allclose(accelerations @ matrix.T, accelerations, atol=1e-15)


def test_two_excitation_directions_already_identify_all_rotation_components():
    result = shared_tilt_diagnostic(_readings([[1, 0, 0], [0, 1, 0]], [.02, -.01, .03]), 2)
    assert result["rotation_rank"] == 3
    np.testing.assert_allclose(result["rotation_projector"], np.eye(3), atol=1e-14)


def test_rotational_information_eigenvalues_follow_the_excitation_scatter():
    acceleration = np.array([[3., 0., 0.], [0., 2., 0.], [0., 0., 1.]])
    mass, sf, sa = 2.0, .1, .2
    readings = [isotropic(mass * a, a, sf, sa) for a in acceleration]
    observed = shared_tilt_diagnostic(readings, mass)["rotation_singular_values"] ** 2
    scatter = acceleration.T @ acceleration / (sf ** 2 + mass ** 2 * sa ** 2)
    s = np.linalg.eigvalsh(scatter)
    expected = np.sort(mass ** 2 * (s.sum() - s))[::-1]
    np.testing.assert_allclose(observed, expected, rtol=1e-13)


def test_shared_axis_gain_can_create_mass_disagreement_and_tempering():
    """A common non-scalar sensor map can disagree about mass despite Q_i=0.

    True m=1, force-channel map diag(1,2,3), excitations e1 and e2. This
    refutes reading O16's 'shared error produces none' as a universal statement.
    """
    powers = []
    for sd in (.1, .03):
        readings = [isotropic([1, 0, 0], [1, 0, 0], sd, sd),
                    isotropic([0, 2, 0], [0, 1, 0], sd, sd)]
        assert max(reading_state(r, 1.0)["Q"] for r in readings) < 1e-9
        powers.append(bartlett_eta(readings, 1.0, [1.0, 1.0])[0])
    assert 0 < powers[1] < powers[0] < 1
    assert powers[1] / powers[0] == pytest.approx((.03 / .1) ** 2, rel=1e-8)


def test_fixed_tilt_does_not_make_every_q_large_or_residual_directions_parallel():
    matrix = rotation([0, 0, 1], 0.1)
    acceleration = np.eye(3)
    readings = [isotropic(2 * matrix @ a, a, .01, .01) for a in acceleration]
    result = shared_tilt_diagnostic(readings, 2)
    assert np.linalg.norm(np.cross(result["residual"][0], result["residual"][1])) > 0.01
    # A push on the fixed rotation axis is unchanged, regardless of tilt strength.
    assert isotropic_minima(readings[2])["projective_min_cost"] == pytest.approx(0, abs=1e-12)


def test_excitation_bounds_keep_mean_direction_and_signed_noise_correction():
    acceleration = np.tile([0., 0., 2.], (5, 1))
    cov = np.broadcast_to(.01 * np.eye(3), (5, 3, 3))
    result = excitation_spectrum(acceleration, cov)
    assert result["certified_uncentered_rank_lower_bound"] == 1
    assert result["certified_centered_rank_lower_bound"] == 0
    assert result["centered_singular_values"].max() == 0
    assert result["noise_subtracted_gram_eigenvalues"][-1] < 0
    assert result["noise_frobenius_expected_squared"] == pytest.approx(.15)


def test_weyl_intervals_cover_a_fixed_bounded_perturbation():
    latent = np.diag([3., 2., 1.])
    error = np.array([[.1, -.2, .1], [.1, 0, .1], [0, .1, -.1]])
    cov = np.broadcast_to(.02 * np.eye(3), (3, 3, 3))
    result = excitation_spectrum(latent + error, cov, alpha=.2)
    assert np.linalg.norm(error, "fro") < result["noise_radius"]
    singular = np.linalg.svd(latent, compute_uv=False)
    assert np.all(singular >= result["singular_value_lower"])
    assert np.all(singular <= result["singular_value_upper"])


def test_fit_cost_and_reference_frame_covariance_match_direct_whitening():
    acceleration = np.array([[1., 0., 0.], [0., 2., 0.], [1., 1., 2.]])
    omega = np.array([.01, -.03, .02])
    factor = np.eye(6) + .1 * np.tril(np.ones((6, 6)), -1)
    covariance = factor @ factor.T
    readings = [Reading(2 * a + 2 * np.cross(omega, a) + [.01, -.02, .03], a, covariance)
                for a in acceleration]
    result = shared_tilt_diagnostic(readings, 2)
    constraint = np.c_[np.eye(3), -2 * np.eye(3)]
    residual = np.array([r.force - 2 * r.acceleration for r in readings])
    residual_cov = constraint @ covariance @ constraint.T
    direct_cost = sum(r @ np.linalg.solve(residual_cov, r) for r in residual)
    assert result["raw_cost"] == pytest.approx(direct_cost, rel=1e-13)
    change = rotation([1, -1, 2], .4)
    block = np.block([[change, np.zeros((3, 3))], [np.zeros((3, 3)), change]])
    transformed = [Reading(change @ r.force, change @ r.acceleration, block @ r.cov @ block.T)
                   for r in readings]
    other = shared_tilt_diagnostic(transformed, 2)
    np.testing.assert_allclose(other["rotation_vector"], change @ result["rotation_vector"], atol=1e-13)
    np.testing.assert_allclose(other["rotation_singular_values"], result["rotation_singular_values"], atol=1e-13)
    assert other["raw_cost"] == pytest.approx(result["raw_cost"], rel=1e-13)
    assert other["fitted_cost"] == pytest.approx(result["fitted_cost"], abs=1e-15)


def test_heldout_scoring_uses_supplied_parameters_and_can_get_worse():
    omega = np.array([.02, -.01, .03])
    fit = shared_tilt_diagnostic(_readings(np.eye(3), omega), 2)
    heldout = _readings([[1, 2, 3], [-2, 1, 0]], omega)
    correct = score_tilt(heldout, 2, fit["delta_mass"], fit["rotation_vector"])
    wrong = score_tilt(heldout, 2, 0.0, -omega)
    assert correct["cost_after"] < 1e-25
    assert correct["improvement"] > 0
    assert wrong["cost_after"] > wrong["cost_before"]
    assert wrong["improvement"] < 0


@pytest.mark.parametrize("value", [0, [1, 2], [1, np.nan, 2]])
def test_bad_vectors_are_rejected(value):
    with pytest.raises(ValueError):
        skew(value)


def test_bad_covariance_and_probability_are_rejected():
    acceleration = np.eye(3)
    for covariance in (np.zeros((3, 3, 3)), np.full((3, 3, 3), np.nan)):
        with pytest.raises(ValueError):
            excitation_spectrum(acceleration, covariance)
    with pytest.raises(ValueError):
        excitation_spectrum(acceleration, np.broadcast_to(np.eye(3), (3, 3, 3)), alpha=True)
