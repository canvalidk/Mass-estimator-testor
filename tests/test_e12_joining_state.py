"""Exact-state implications of appending readings to the E12 completed law."""
import math

import numpy as np
import pytest

from explore.e12_scalar_completion import completed_log_evidence
from explore.e12_scalar_witness import pooled_witness
from explore.e12_tangent_pairings import tangent_pairings
from met.conditional_score import gaussian_score_accepts
from met.full import isotropic


ANGLES = np.linspace(.001, math.pi / 2 - .001, 301)


def rows_from_lines(energies, angles):
    return np.sqrt(energies)[:, None] * np.column_stack((np.sin(angles), np.cos(angles)))


def contrasts(rows, alpha=1/8):
    value = completed_log_evidence(rows, ANGLES, alpha=alpha)
    return value - value[0]


@pytest.mark.parametrize("rows", [np.array([[1., 2.], [3., 4.]]),
                                  np.array([[2., 1.], [-4., -2.]]), np.zeros((2, 2))])
def test_three_positive_extensions_identify_the_preappend_scatter(rows):
    scatter = rows.T @ rows
    baseline = np.linalg.det(scatter)
    times = np.array([1., 2., 3.])
    extensions = np.column_stack((np.ones(3), times))
    changes = []
    for probe in extensions:
        joined = np.vstack((rows, probe))
        changes.append(np.linalg.det(joined.T @ joined) - baseline)
    coeff = np.linalg.solve(np.vander(times, 3, increasing=True), changes)
    inferred = np.array([[coeff[2], -coeff[1] / 2], [-coeff[1] / 2, coeff[0]]])
    np.testing.assert_allclose(inferred, scatter, atol=2e-12)


@pytest.mark.parametrize("alpha", [1/16, 1/8, 1.0])
def test_current_law_point_and_score_hulls_agree_but_a_positive_append_changes_law(alpha):
    u = np.full(2, 1/math.sqrt(2))
    left = np.array([u, 2*u])
    right = np.array([math.sqrt(2.5)*u, math.sqrt(2.5)*u])
    np.testing.assert_allclose(left.T @ left, right.T @ right, atol=2e-15)
    np.testing.assert_allclose(contrasts(left, alpha), contrasts(right, alpha), atol=2e-14)
    # The score-hull's other input is the entire T45 score set, fixed by the
    # common Gram; spot-check the repo implementation at several levels.
    left_readings = [isotropic([f], [a], 1., 1.) for f, a in left]
    right_readings = [isotropic([f], [a], 1., 1.) for f, a in right]
    candidates = np.tan(ANGLES[::40])
    for p in [.65, .8, .99]:
        left_scores = gaussian_score_accepts(left_readings, candidates, probability=p)
        right_scores = gaussian_score_accepts(right_readings, candidates, probability=p)
        np.testing.assert_array_equal(left_scores, right_scores)
    probe = np.array([1., 2.])
    left_joined = np.vstack((left, probe))
    right_joined = np.vstack((right, probe))
    np.testing.assert_allclose(left_joined.T @ left_joined,
                               right_joined.T @ right_joined, atol=3e-15)
    assert np.ptp(contrasts(left_joined, alpha) - contrasts(right_joined, alpha)) > 1.5e-7


def log_factor_derivative(rows, z):
    """Derivative in z using the displayed real-angle factor after complex substitution."""
    fit = pooled_witness(rows, alpha=1/8)
    energy = np.sum(rows**2, axis=1)
    c = energy / (1 + energy)
    phase = np.exp(2j * np.arctan2(rows[:, 0], rows[:, 1]))
    minimum_diagonal = np.min(np.sum(rows**2, axis=0))
    cap = min(fit.precision, 1 + minimum_diagonal**2)
    pilot_phase = np.exp(2j * fit.angle)
    value = (cap / 8) * (1/pilot_phase - pilot_phase/z**2)
    gap_precision = fit.gap**2 / fit.largest
    gate = math.exp(-gap_precision) * fit.positive_residual / (1 + fit.positive_residual)
    for scale, q in zip(c, phase):
        numerator = -.25 * scale * (1/q - q/z**2)
        denominator = 1 + .5*scale - .25*scale*(z/q + q/z)
        value -= gate/(1 + np.sum(c)) * numerator/denominator
    return value


def test_complex_curve_derivative_agrees_with_completed_law_and_recovers_rows():
    energies = np.array([1., 1., 2., .6, 1., 3.])
    angles = np.array([.23, .23, .23, -.7, .85, 1.3])
    rows = rows_from_lines(energies, angles)
    rows[1] *= -1
    # Check the rational derivative against numerical differentiation of the actual law.
    angles_to_test = np.array([.22, .71, 1.1])
    z = np.exp(2j * angles_to_test)
    step = 1e-6
    finite_diff = (completed_log_evidence(rows, angles_to_test+step, alpha=1/8)
                   - completed_log_evidence(rows, angles_to_test-step, alpha=1/8)) / (2*step)
    np.testing.assert_allclose(2j*z*log_factor_derivative(rows, z), finite_diff, atol=2e-9)

    # Integrate the derivative around isolated roots; the marker at 1.3 is known once.
    c = energies/(1+energies)
    radius = (np.sqrt(1+c)-1)/(np.sqrt(1+c)+1)
    roots = radius*np.exp(2j*angles)
    distinct = roots[[0, 2, 3, 4, 5]]
    residues = []
    for i, pole in enumerate(distinct):
        others = np.delete(distinct, i)
        contour_radius = min(abs(pole), np.min(np.abs(others-pole))) / 8
        offsets = contour_radius*np.exp(2j*math.pi*np.arange(128)/128)
        residues.append(np.mean(offsets*log_factor_derivative(rows, pole + offsets)))
    np.testing.assert_allclose(np.array(residues)/residues[-1], [2, 1, 1, 1, 1], atol=2e-10)
    radial = np.abs(distinct)
    recovered_c = ((1+radial)/(1-radial))**2 - 1
    np.testing.assert_allclose(recovered_c/(1-recovered_c), energies[[0, 2, 3, 4, 5]], atol=2e-13)


def test_completed_updates_ignore_padding_signs_and_order_with_common_extensions():
    rows = np.array([[2., 1.], [1., 3.], [-1., .5]])
    changed = np.vstack((-rows[2], rows[0], np.zeros(2), -rows[1]))
    for append in [np.empty((0, 2)), np.array([[1., 2.]]), np.array([[2., -1.], [.2, .3]])]:
        np.testing.assert_allclose(contrasts(np.vstack((rows, append))),
                                   contrasts(np.vstack((changed, append))), atol=1e-14)


def vector_from_shape(energy, centre, width):
    force = np.array([math.sqrt(energy), 0.])
    acceleration = math.sqrt(energy) * np.array([-centre, width])
    return isotropic(force, acceleration, 1., 1.)


def test_local_tangent_within_curve_and_merged_shapes_agree_for_all_sampled_masses():
    parts = [vector_from_shape(1., .3, 1.5), vector_from_shape(3., .3, 1.5),
             vector_from_shape(2., -.5, .4)]
    merged = [vector_from_shape(4., .3, 1.5), vector_from_shape(2., -.5, .4)]
    for mass in [.001, .5, 1., 3., 11.]:
        left = tangent_pairings(parts, mass)
        right = tangent_pairings(merged, mass)
        expected = (1 + mass**2) * (4 * 1.5**2 / ((mass - .3)**2 + 1.5**2)
                                    + 2 * .4**2 / ((mass + .5)**2 + .4**2))
        assert left.within == pytest.approx(expected)
        np.testing.assert_allclose([left.within, left.between, left.common],
                                   [right.within, right.between, right.common])
        assert left.degrees != right.degrees  # This classification is of energies.


def test_same_declared_count_scatter_and_exterior_split_can_hide_tangent_disagreement():
    pairs = [
        [np.diag([3., 1.]), np.diag([1., 3.])],
        [np.array([[2., 1.], [1., 2.]]), np.array([[2., -1.], [-1., 2.]])],
    ]
    for matrices, expected in zip(pairs, [(3., 1.), (4., 0.)]):
        np.testing.assert_allclose(sum(matrices), 4*np.eye(2))
        assert sum(np.linalg.det(g) for g in matrices) == pytest.approx(6.)
        measurements = []
        for gram in matrices:
            force, accel = np.linalg.cholesky(gram)
            measurements.append(isotropic(force, accel, 1., 1.))
        result = tangent_pairings(measurements, 1.)
        np.testing.assert_allclose([result.within, result.between], expected, atol=1e-13)
        assert result.common == pytest.approx(0., abs=1e-13)


def test_antiparallel_scalar_reading_creates_one_positive_within_energy_spike():
    diagonal = [isotropic([1.], [1.], 1., 1.), isotropic([1.], [-1.], 1., 1.)]
    root_two = math.sqrt(2.)
    axes = [isotropic([root_two], [0.], 1., 1.), isotropic([0.], [root_two], 1., 1.)]
    for mass in [.25, 1., 4.]:
        diag_split = tangent_pairings(diagonal, mass)
        axis_split = tangent_pairings(axes, mass)
        expected = [2., 0., 0.] if mass == 1. else [0., 2., 0.]
        np.testing.assert_allclose([diag_split.within, diag_split.between, diag_split.common],
                                   expected, atol=2e-13)
        np.testing.assert_allclose([axis_split.within, axis_split.between, axis_split.common],
                                   [0., 2., 0.], atol=2e-13)


def test_regular_polygon_annihilates_all_retained_nonconstant_harmonics():
    angles = 2 * math.pi * np.arange(17)/17
    frequency = np.arange(1, 16)
    np.testing.assert_allclose(np.sum(np.exp(1j*angles[:, None]*frequency), axis=0), 0., atol=4e-14)
    rows = rows_from_lines(np.ones(17), angles)
    np.testing.assert_allclose(rows.T @ rows, 8.5*np.eye(2), atol=4e-15)
