"""Verify study algebra and reproducibility; no probabilistic pass threshold."""
import numpy as np
import pytest
from scipy import special
import yaml

from met.joint_calibration import gaussian_joint_observed_moments
from met.joint_mass import common_acceleration, gaussian_profile_cost
from met.joint_sanity import JointQuadrature, profile_derivatives, run


def test_scatter_acceleration_matches_event_derivatives_and_independent_differences():
    data = np.array([[2., .7, 1.2], [1.1, 2.3, -.2], [.7, 1.2, .4]])
    masses, variance = np.array([1.2, .8]), np.array([1., 9., 1.])
    events = [common_acceleration(y[:2, None], y[2:], np.diag(variance)) for y in data]
    observed = gaussian_joint_observed_moments(events, masses)
    _, score, sensitivity = profile_derivatives(masses, data.T@data, variance)
    np.testing.assert_allclose(score, observed["score"], rtol=1e-13, atol=1e-14)
    np.testing.assert_allclose(sensitivity, observed["sensitivity"], rtol=1e-13, atol=1e-14)
    step = 1e-5
    numeric = np.column_stack([(profile_derivatives(masses+step*d, data.T@data, variance)[1]
                              -profile_derivatives(masses-step*d, data.T@data, variance)[1])/(2*step)
                              for d in np.eye(2)])
    np.testing.assert_allclose(sensitivity, -numeric, rtol=1e-8, atol=1e-9)


def test_quadrature_matches_direct_joint_profile_weights():
    quad = JointQuadrature(12, .95)
    data, scale = np.array([1.2, 1.7, .9]), .3
    event = common_acceleration(data[:2, None], data[2:], scale**2*np.eye(3))
    costs = np.array([gaussian_profile_cost(event, point) for point in quad.masses.reshape(-1, 2)])
    logweights = np.log(quad.prior_weights.ravel())-.5*costs
    weights = np.exp(logweights-special.logsumexp(logweights))
    expected = []
    for j, order in enumerate(quad.orders):
        indices = np.searchsorted(np.cumsum(weights[order]), quad.ranks)
        expected.append(quad.values[order[indices], j])
    np.testing.assert_array_equal(quad.intervals(data, scale)[0], expected)


def test_small_study_replays_identically_and_records_numerical_controls(tmp_path):
    config = yaml.safe_load(open("presets/joint_sanity.yaml", encoding="utf-8"))
    config["replicates"] = 3
    config["small_noise"].update(scales=[.15], angle_nodes=24,
        refinement_nodes=32, refinement_replicates=2)
    config["repeated"]["blocks"] = [50]
    first, meta = run(config, tmp_path/"first")
    second, _ = run(config, tmp_path/"second")
    assert first == second
    assert meta["preset"] == config
    assert len(meta["refinement"]) == 1
    with np.load(tmp_path/"first"/"values.npz") as a, np.load(tmp_path/"second"/"values.npz") as b:
        assert a.files == b.files
        for key in a.files:
            np.testing.assert_array_equal(a[key], b[key])


def test_standard_preset_runner_routes_the_study_and_refuses_parallel_reinterpretation(tmp_path):
    from met import preset, runner
    config = yaml.safe_load(open("presets/joint_sanity.yaml", encoding="utf-8"))
    config["replicates"] = 2
    config["small_noise"].update(scales=[.15], angle_nodes=12,
        refinement_nodes=16, refinement_replicates=2)
    config["repeated"]["blocks"] = [20]
    path = tmp_path/"recipe.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    assert preset.load(path)["kind"] == "joint_sanity"
    target, record = runner.run(path, tmp_path/"results", workers=1, progress=None)
    assert record["preset"] == config
    assert (target/"summary.csv").is_file()
    with pytest.raises(ValueError, match="serial"):
        runner.run(path, tmp_path/"other", workers=2, progress=None)
    config["small_noise"]["refinement_nodes"] = 8
    with pytest.raises(ValueError, match="refinement"):
        preset.validate(config)
