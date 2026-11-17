"""T71 commutants and single-reading fits in isotropic noise units."""

import numpy as np
import pytest

from met.full import isotropic, reading_state


J = np.array([[0., -1.], [1., 0.]])
SF, SA = .3, .1


def commutant_dimension(*generators):
    """Count solutions of RM=MR by stacking the column-vectorized equations."""
    dimension = generators[0].shape[0]
    identity = np.eye(dimension)
    constraints = np.vstack([
        np.kron(identity, rotation) - np.kron(rotation.T, identity)
        for rotation in generators
    ])
    return dimension**2 - np.linalg.matrix_rank(constraints)


def cross_matrix(axis):
    return np.cross(axis, np.eye(3)).T


def axis_map(axis, longitudinal, planar_scale, planar_turn):
    projection = np.outer(axis, axis)
    return (longitudinal * projection + planar_scale * (np.eye(3) - projection)
            + planar_turn * cross_matrix(axis))


def quarter_turn(axis):
    return np.outer(axis, axis) + cross_matrix(axis)


def state(force, acceleration):
    return reading_state(isotropic(SF * force, SA * acceleration, SF, SA), SF / SA)


@pytest.mark.parametrize("k, rotation_dimension", [(1, 1), (2, 2), (3, 1)])
def test_t71_rotation_and_orthogonal_commutants(k, rotation_dimension):
    rotations = {
        1: (np.eye(1),),
        2: (J,),
        3: (quarter_turn(np.eye(3)[0]), quarter_turn(np.eye(3)[2])),
    }[k]
    reflection = np.eye(k)
    reflection[0, 0] = -1
    assert commutant_dimension(*rotations) == rotation_dimension
    assert commutant_dimension(*rotations, reflection) == 1


def test_t71_axial_commutant_and_planar_reflection():
    axis = np.array([1., 2., -3.]) / np.sqrt(14.)
    rotation = quarter_turn(axis)
    operator = axis_map(axis, 2., .7, -.4)
    assert commutant_dimension(rotation) == 3
    np.testing.assert_allclose(operator @ rotation, rotation @ operator, atol=1e-14)
    np.testing.assert_allclose(operator @ axis, 2. * axis, atol=1e-14)
    reflection = np.diag([1., -1.])
    np.testing.assert_array_equal(reflection @ J @ reflection.T, -J)
    np.testing.assert_array_equal((-np.eye(2)) @ J @ (-np.eye(2)), J)


def test_t71_complex_fit_vanishes_but_real_misfit_matches_reading_state():
    rng = np.random.default_rng(2026112302)
    for force, acceleration in rng.normal(size=(4, 2, 2)):
        ratio = complex(*force) / complex(*acceleration)
        complex_map = ratio.real * np.eye(2) + ratio.imag * J
        np.testing.assert_allclose(complex_map @ acceleration, force, atol=1e-13)
        singular_values = np.linalg.svd(np.column_stack((force, acceleration)),
                                        compute_uv=False)
        real_misfit = singular_values[-1]**2
        assert real_misfit > 0.
        assert state(force, acceleration)["Q"] == pytest.approx(real_misfit, abs=1e-9)


def test_t71_antiparallel_pairs_satisfy_reality_but_fail_positivity():
    reading = state(np.array([-2., 0.]), np.array([1., 0.]))
    assert reading["Q"] == pytest.approx(0., abs=1e-9)
    assert reading["Q_plus"] == pytest.approx(1., abs=1e-9)


def test_t71_push_axis_has_scalar_misfit_and_an_independent_axis_fits():
    rng = np.random.default_rng(2026112303)
    for force, acceleration in rng.normal(size=(4, 2, 3)):
        pair = np.column_stack((force, acceleration))
        left, singular_values, right = np.linalg.svd(pair, full_matrices=False)
        # Nearest rank-one pair over a latent push axis selected with the fit.
        tied_axis = left[:, 0]
        force_along, push_along = singular_values[0] * right[0]
        fitted_acceleration = push_along * tied_axis
        push_map = axis_map(tied_axis, force_along / push_along, .7, -.4)
        fitted_pair = np.column_stack((push_map @ fitted_acceleration,
                                        fitted_acceleration))
        tied_misfit = np.linalg.norm(pair - fitted_pair, "fro")**2
        assert tied_misfit > 0.
        assert tied_misfit == pytest.approx(singular_values[-1]**2, abs=1e-12)
        assert state(force, acceleration)["Q"] == pytest.approx(tied_misfit, abs=1e-9)

        # A separately chosen axis leaves both channels in its transverse plane.
        free_axis = np.cross(acceleration, force)
        free_axis /= np.linalg.norm(free_axis)
        size_squared = acceleration @ acceleration
        scale = (acceleration @ force) / size_squared
        turn = np.cross(free_axis, acceleration) @ force / size_squared
        free_map = axis_map(free_axis, 1., scale, turn)
        rotation = quarter_turn(free_axis)
        np.testing.assert_allclose(free_map @ rotation, rotation @ free_map, atol=1e-13)
        np.testing.assert_allclose(free_map @ acceleration, force, atol=1e-13)


def test_t71_perpendicular_pair_separates_push_stabilizer_from_free_axis():
    force, free_axis, push_axis = np.eye(3)
    assert state(force, push_axis)["Q"] == pytest.approx(1., abs=1e-9)
    # With the axis declared, both fitted channels lie on it. Their
    # transverse residuals have total squared length 1 in this example.
    transverse = np.eye(3) - np.outer(push_axis, push_axis)
    pair = np.column_stack((force, push_axis))
    assert np.linalg.norm(transverse @ pair, "fro")**2 == pytest.approx(1.)
    # A scalar-tied fit can attain the same residual with nonzero force and
    # acceleration by choosing their common line along the sum.
    fitted = (force + push_axis) / 2
    assert np.linalg.norm(pair - np.column_stack((fitted, fitted)), "fro")**2 == pytest.approx(1.)
    untied_map = axis_map(free_axis, 1., 0., 1.)
    np.testing.assert_array_equal(untied_map @ push_axis, force)
