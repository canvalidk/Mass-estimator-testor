"""Proof checks for the shared-drive tangent and positive coupling counterexample."""
from collections import defaultdict
from fractions import Fraction
from itertools import permutations, product
from math import prod

import numpy as np
import pytest

from met.conditional_score import gaussian_score_test
from met.full import Reading


def test_shared_drive_projection_with_correlated_instruments():
    masses = np.array([1.25, 3.0, 0.5])
    data = np.array([0.3, 2.5, 1.0, 0.8])
    lower = np.array([[2, 0, 0, 0], [1, 2, 0, 0], [-1, 0, 1, 0], [.3, .4, .2, 1]])
    covariance = lower @ lower.T
    defect = np.column_stack([np.eye(3), -masses])
    noise = defect @ covariance @ defect.T
    drive_row = np.eye(4)[-1] - np.linalg.solve(noise, defect @ covariance[:, -1]) @ defect
    drive = drive_row @ data
    np.testing.assert_allclose(drive_row @ covariance @ defect.T, 0, atol=1e-14)
    graph = np.r_[masses, 1]
    weighted_graph = np.linalg.solve(covariance, graph)
    assert drive == pytest.approx(weighted_graph @ data / (weighted_graph @ graph))
    correction = covariance @ defect.T @ np.linalg.solve(noise, defect @ data)
    np.testing.assert_allclose(data - correction, graph * drive)


def test_single_raw_lump_refits_drive_but_preserves_scalar_fieller_score():
    masses = np.array([1., 2.])
    acceleration = 1.
    forces = np.array([2., 2.])
    lump = np.ones(2)
    residual = forces - masses * acceleration
    noise = np.eye(2) + np.outer(masses, masses)
    fine_drive = acceleration + masses @ np.linalg.solve(noise, residual)
    lump_variance = lump @ noise @ lump
    coarse_drive = acceleration + masses.sum() * (lump @ residual) / lump_variance
    assert fine_drive == pytest.approx(7/6)
    assert coarse_drive == pytest.approx(14/11)
    fine_center = lump @ residual / fine_drive
    coarse_center = lump @ residual / coarse_drive
    assert fine_center != pytest.approx(coarse_center)
    fine_pivot = fine_center**2 / (lump_variance / fine_drive**2)
    coarse = Reading([forces.sum()], [acceleration], np.diag([2., 1.]))
    assert fine_pivot == pytest.approx(1/11)
    assert gaussian_score_test([coarse], masses.sum())["statistic"] == pytest.approx(fine_pivot)


@pytest.mark.parametrize("complementary", [True, False])
def test_conditional_lump_information_gap(complementary):
    lump = np.ones(2)
    force_covariances = [np.diag([1., 9.]), np.diag([9., 1.])]
    if not complementary:
        # Different variances on the body contrast, equal coupling to the lump.
        force_covariances = [np.eye(2), np.array([[2., -1.], [-1., 2.]])]
    noises = [force + np.outer(lump, lump)/10 for force in force_covariances]
    drives = [1., 1.]
    precision = sum(x*x*np.linalg.inv(noise) for x, noise in zip(drives, noises))
    reciprocal = np.linalg.solve(precision, lump)
    fine_variance = lump @ reciprocal
    coarse_variance = 1/sum(x*x/(lump @ noise @ lump) for x, noise in zip(drives, noises))
    fine_weights = [x*np.linalg.solve(noise, reciprocal) for x, noise in zip(drives, noises)]
    coarse_weights = [x*coarse_variance*lump/(lump @ noise @ lump)
                      for x, noise in zip(drives, noises)]
    np.testing.assert_allclose(sum(x*w for x, w in zip(drives, fine_weights)), lump)
    np.testing.assert_allclose(sum(x*w for x, w in zip(drives, coarse_weights)), lump)
    gap = sum((z-w) @ noise @ (z-w) for z, w, noise in zip(coarse_weights, fine_weights, noises))
    assert coarse_variance - fine_variance == pytest.approx(gap)
    normalized = [noise @ lump / (lump @ noise @ lump) for noise in noises]
    if complementary:
        assert fine_variance == pytest.approx(2.)
        assert coarse_variance == pytest.approx(5.2)
        assert not np.allclose(normalized[0], normalized[1])
    else:
        np.testing.assert_allclose(normalized[0], normalized[1])
        np.testing.assert_allclose(fine_weights, coarse_weights)


def test_rotating_anisotropic_single_spatial_push_preserves_tensor_gap():
    eye = np.eye(2)
    lump = np.ones(2)
    defect = np.column_stack([np.eye(4), -np.kron(lump[:, None], eye)])
    noise = defect @ np.diag([1., 9., 9., 1., .1, .1]) @ defect.T
    sum_forces = np.kron(lump, eye)
    observed = np.array([1., 0., 0., 0.])
    rotation = np.array([[.8, -.6], [.6, .8]])
    corrections = []
    for frame in (eye, rotation):
        drive = frame @ np.ones(2)
        frame_stack = np.kron(eye, frame)
        tangent = np.kron(eye, drive[:, None])
        frame_noise = frame_stack @ noise @ frame_stack.T
        reciprocal_normal = np.linalg.solve(frame_noise, tangent).T
        precision = reciprocal_normal @ tangent
        correction = np.linalg.solve(precision, reciprocal_normal @ frame_stack @ observed)
        corrections.append(correction)
        fine_variance = lump @ np.linalg.solve(precision, lump)
        lump_noise = sum_forces @ frame_noise @ sum_forces.T
        coarse_variance = 1/(drive @ np.linalg.solve(lump_noise, drive))
        assert fine_variance == pytest.approx(2.)
        assert coarse_variance == pytest.approx(26/5)
        # The returned spatial corrections on each body cannot both lie along the drive.
        allocated = frame_noise @ sum_forces.T @ np.linalg.solve(lump_noise, drive)
        projection = tangent @ np.linalg.lstsq(tangent, allocated, rcond=None)[0]
        assert np.linalg.norm(allocated - projection) > .5
    np.testing.assert_allclose(corrections[0], corrections[1], atol=1e-14)


def _exchangeable_pair_perturbation():
    """Cofactor construction in the ring of signed multiplicative-scale measures."""
    half = Fraction(1, 2)
    shapes = [Fraction(i) for i in (1, 2, 3, 4)]
    factors = [
        [[(r, half), (1/r, half)] for r in shapes],
        [[(r + 1/r, 1)] for r in shapes],
        [[(1/(r + 1/r), 1)] for r in shapes],
    ]
    joint = defaultdict(Fraction)
    for j, shape in enumerate(shapes):
        other = [k for k in range(4) if k != j]
        for columns in permutations(other):
            inversions = sum(a > b for i, a in enumerate(columns) for b in columns[i+1:])
            signed_cofactor = (-1)**(j + inversions)
            row_terms = [row[column] for row, column in zip(factors, columns)]
            for terms in product(*row_terms):
                scale = prod(value for value, _ in terms)
                weight = signed_cofactor * prod(probability for _, probability in terms)
                x, y = scale*shape, scale/shape
                joint[x, y] += weight*half
                joint[y, x] += weight*half
    return {pair: weight for pair, weight in joint.items() if weight}


def _push(joint, observable):
    weights = defaultdict(Fraction)
    for pair, weight in joint.items():
        weights[observable(*pair)] += weight
    return {value: weight for value, weight in weights.items() if weight}


def test_exchangeable_joint_laws_agree_on_both_part_laws_lump_and_reduced_mass():
    signed = _exchangeable_pair_perturbation()
    positive_mass = sum(weight for weight in signed.values() if weight > 0)
    assert positive_mass == -sum(weight for weight in signed.values() if weight < 0) == 12
    assert len(signed) == 72
    assert all(min(x, y) > 0 and signed[y, x] == weight for (x, y), weight in signed.items())
    for observable in (lambda x, y: x, lambda x, y: y, lambda x, y: x+y,
                       lambda x, y: x*y/(x+y)):
        assert _push(signed, observable) == {}
    gap = _push(signed, lambda x, y: 2*x+y)[Fraction(80, 17)]
    assert abs(gap)/positive_mass == Fraction(1, 24)
