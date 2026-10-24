"""Exact/numerical examples for robot_payload_literature_bridge.md.

These verify the proposed data reduction and moment readout, not robot
performance or any published solver. Requires Python 3 and NumPy.
Prints JSON; writes nothing unless --output is supplied.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def h(x: float, sigma: float) -> float:
    """Positive-truncated Gaussian mean, sufficient for the nonnegative examples."""
    t = x / sigma
    assert t >= 0 and sigma > 0
    phi = math.exp(-0.5 * t * t) / math.sqrt(2 * math.pi)
    cdf = 0.5 * (1 + math.erf(t / math.sqrt(2)))
    return x + sigma * phi / cdf


def skew(v: np.ndarray) -> np.ndarray:
    x, y, z = v
    return np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])


def pseudo(points: np.ndarray, masses: np.ndarray) -> np.ndarray:
    q = np.column_stack((points, np.ones(len(points))))
    return np.einsum("i,ij,ik->jk", masses, q, q)


def checks() -> dict:
    result = {}
    mass = 2.
    gravity = np.array([0., 0., -9.81])
    force = -mass * gravity
    wrong = h(float(np.linalg.norm(force)), .05) / h(0., .05)
    right = h(float(np.linalg.norm(force)), .05) / h(9.81, .05)
    assert np.isclose(right, mass, atol=1e-12)
    result["stationary_known_direction"] = {
        "true_mass_kg": mass, "contact_force_N": force.tolist(),
        "raw_kinematic_acceleration_m_s2": [0., 0., 0.],
        "effective_acceleration_m_s2": (-gravity).tolist(),
        "wrong_pair_readout_kg": wrong, "correct_pair_readout_kg": right,
        "force_sigma_N": .05, "acceleration_sigma_m_s2": .05,
    }

    a_sensor = np.zeros(3)
    omega = np.array([0., 0., 3.])
    angular_acceleration = np.array([0., 2., 0.])
    center = np.array([.2, 0., 0.])
    rotational = (np.cross(angular_acceleration, center)
                  + np.cross(omega, np.cross(omega, center)))
    a_center = a_sensor + rotational
    contact = mass * (a_center - gravity)
    matrix_form = (mass * (a_sensor - gravity)
                   + (skew(angular_acceleration) + skew(omega) @ skew(omega))
                   @ (mass * center))
    assert np.allclose(contact, matrix_form)
    assert np.linalg.norm(np.cross(contact, a_sensor - gravity)) > 1.
    result["rotating_sensor_origin"] = {
        "center_offset_m": center.tolist(), "omega_rad_s": omega.tolist(),
        "angular_acceleration_rad_s2": angular_acceleration.tolist(),
        "center_acceleration_m_s2": a_center.tolist(),
        "contact_force_N": contact.tolist(),
        "naive_force_over_sensor_effective_acceleration_kg":
            float(np.linalg.norm(contact) / np.linalg.norm(a_sensor - gravity)),
        "correct_force_over_center_effective_acceleration_kg":
            float(np.linalg.norm(contact) / np.linalg.norm(a_center - gravity)),
    }

    # Unknown constant sensor bias and known effective acceleration b_i.
    # Force is m*b_i + bias. Rank is 4 iff at least two b_i differ.
    known_b = np.tile(-gravity, (5, 1))
    varying_b = known_b.copy()
    varying_b[-1, 0] += .2
    def design(b):
        return np.vstack([np.column_stack((row, np.eye(3))) for row in b])
    assert np.linalg.matrix_rank(design(known_b)) == 3
    assert np.linalg.matrix_rank(design(varying_b)) == 4
    force_sigma = .05
    D = design(varying_b)
    fisher = D.T @ D / force_sigma**2
    schur = float(fisher[0, 0] - fisher[0, 1:] @
                  np.linalg.solve(fisher[1:, 1:], fisher[1:, 0]))
    centered_information = float(np.sum((varying_b - varying_b.mean(0))**2)
                                 / force_sigma**2)
    assert np.isclose(schur, centered_information, rtol=1e-9)
    result["constant_bias_identifiability"] = {
        "constant_b_design_rank": 3, "varying_b_design_rank": 4,
        "mass_information_after_unknown_bias": centered_information,
        "schur_complement_check": schur,
    }

    # Transform independent raw (force,acceleration) errors to differences
    # with one shared reference. Use nonzero channel cross-covariance.
    n = 20
    B = np.eye(6)
    B[0, 3] = .3
    C = B @ B.T
    L = np.column_stack((-np.ones(n), np.eye(n)))
    transform = np.kron(L, np.eye(6))
    covariance = transform @ np.kron(np.eye(n + 1), C) @ transform.T
    expected = np.kron(np.eye(n) + np.ones((n, n)), C)
    assert np.allclose(covariance, expected)
    averaging = np.kron(np.ones((1, n)) / n, np.eye(6))
    mean_covariance = averaging @ covariance @ averaging.T
    assert np.allclose(mean_covariance, (1 + 1/n) * C)
    result["shared_tare_covariance"] = {
        "subsequent_samples": n,
        "true_mean_variance_multiplier_of_raw_variance": 1 + 1/n,
        "independence_assumption_multiplier": 2/n,
        "variance_understatement_factor": (1 + 1/n) / (2/n),
    }

    # Two mass distributions can produce the same gravitational wrench
    # while having different inertias. Uniform centered spheres suffice.
    radii = [.05, .10]
    inertias = [.4 * mass * radius**2 for radius in radii]
    result["gravity_only_inertia_nonidentifiability"] = {
        "same_mass_kg": mass, "same_center_m": [0., 0., 0.],
        "sphere_radii_m": radii,
        "isotropic_inertias_kg_m2": inertias,
        "same_static_force_N": force.tolist(), "same_static_torque_Nm": [0., 0., 0.],
    }
    assert np.isclose(inertias[1] / inertias[0], 4.)

    # Positive rotational inertia does not imply physical realizability.
    I = np.diag([1., 1., 3.])
    S = .5 * np.trace(I) * np.eye(3) - I
    assert np.linalg.eigvalsh(I).min() > 0
    assert np.linalg.eigvalsh(S).min() < 0
    result["positive_but_impossible_inertia"] = {
        "inertia_diagonal": np.diag(I).tolist(),
        "spatial_second_moment_diagonal": np.diag(S).tolist(),
    }

    # An acceleration-weighted pseudo-inertia readout extends the mass
    # ratio algebraically while preserving the cone; it selects no law.
    points = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [-1., -1., -1.]])
    J1 = pseudo(points, np.array([.2, .3, .1, .4]))
    J2 = pseudo(.5 * points + [.2, 0., 0.], np.array([.4, .2, .8, .6]))
    laws = np.array([.25, .75])
    acceleration = np.array([.2, 1.5])
    weights = laws * acceleration
    weights /= weights.sum()
    Jhat = weights[0] * J1 + weights[1] * J2
    force_magnitudes = acceleration * np.array([J1[3, 3], J2[3, 3]])
    ratio = float(laws @ force_magnitudes / (laws @ acceleration))
    assert np.isclose(Jhat[3, 3], ratio)
    assert np.linalg.eigvalsh(Jhat).min() > 0
    # Rigid change of spatial frame acts by congruence.
    angle = .7
    H = np.eye(4)
    H[:3, :3] = [[math.cos(angle), -math.sin(angle), 0.],
                 [math.sin(angle), math.cos(angle), 0.], [0., 0., 1.]]
    H[:3, 3] = [.1, .2, -.3]
    averaged_transform = weights[0]*(H @ J1 @ H.T) + weights[1]*(H @ J2 @ H.T)
    assert np.allclose(averaged_transform, H @ Jhat @ H.T)
    result["weighted_pseudo_inertia_readout"] = {
        "readout_mass_kg": ratio,
        "mass_matrix_entry_kg": float(Jhat[3, 3]),
        "minimum_pseudo_inertia_eigenvalue_in_chosen_units":
            float(np.linalg.eigvalsh(Jhat).min()),
        "frame_congruence_residual": float(np.max(np.abs(
            averaged_transform - H @ Jhat @ H.T))),
    }
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"purpose": "mechanics and algebra checks; no empirical robot benchmark",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "checks": checks()}
    rendered = json.dumps(report, indent=2)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
