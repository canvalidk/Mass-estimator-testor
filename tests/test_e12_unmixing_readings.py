"""Population identification in the additive Gaussian force–acceleration experiment."""
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import erf, roots_hermitenorm
from scipy.stats import multivariate_normal, norm


def test_deconvolution_recovers_joint_atoms_with_correlated_noise():
    # Labelled scalar readings (F1,F2,a) live on a known finite integer grid
    # in this fixture. Characteristic values on its dual grid suffice.
    means = np.array([[2, 3, 1], [2, 4, 2], [-1, -3, -1], [0, 0, 0]])
    weights = np.array([.2, .25, .4, .15])
    covariance = np.array([[.17, .05, -.02], [.05, .26, .04], [-.02, .04, .22]])
    side = 9
    freq = 2 * np.pi * np.fft.fftfreq(side)
    points = np.stack(np.meshgrid(freq, freq, freq, indexing="ij"), axis=-1)
    projected_var = np.einsum("...i,ij,...j->...", points, covariance, points)

    # Integrate exp(i t.Y) under the observed Gaussian mixture. For each t,
    # the one-dimensional projection t.Y has known variance t.Sigma.t.
    nodes, quadrature = roots_hermitenorm(32)
    quadrature /= np.sqrt(2 * np.pi)
    observed = np.zeros(points.shape[:-1], dtype=complex)
    for weight, mean in zip(weights, means):
        locations = (points @ mean)[..., None] + np.sqrt(projected_var)[..., None] * nodes
        observed += weight * np.sum(np.exp(1j * locations) * quadrature, axis=-1)
    latent_transform = observed / np.exp(-.5 * projected_var)
    recovered = np.fft.fftn(latent_transform) / side**3

    expected = np.zeros_like(recovered, dtype=float)
    for mean, weight in zip(means, weights):
        expected[tuple(mean % side)] = weight
    np.testing.assert_allclose(recovered, expected, atol=2e-13)
    identifiable = {}
    for first, second, push in means:
        if push:
            masses = (first / push, second / push)
            identifiable[masses] = recovered[first % side, second % side, push % side].real
    assert identifiable == pytest.approx({(2, 3): .2, (1, 2): .25, (1, 3): .4})
    assert recovered[0, 0, 0].real == pytest.approx(.15)


def test_real_density_shifts_recover_the_conditional_mean_transform():
    covariance = np.array([[1.2, -.3], [-.3, .7]])
    means = np.array([[1.5, .5], [-2., -1.], [0., 0.]])
    weights = np.array([.3, .25, .45])
    kernels = [multivariate_normal(mean=mu, cov=covariance) for mu in means]

    def joint_weights(point):
        return weights * np.array([kernel.pdf(point) for kernel in kernels])

    y = np.array([1.2, -.6])
    predictive = sum(joint_weights(y))
    posterior = joint_weights(y) / predictive
    for t in [np.array([.2, -.1]), np.array([1.1, .4]), np.array([-.5, -.6])]:
        ratio = sum(joint_weights(y + covariance @ t)) / predictive
        from_readings = np.exp(t @ y + .5 * t @ covariance @ t) * ratio
        assert from_readings == pytest.approx(posterior @ np.exp(means @ t), abs=2e-14)


def test_vector_reciprocal_recovers_each_mass_and_lumped_mass_after_rotation():
    rotation = np.array([[1., 2., 2.], [2., 1., -2.], [-2., 2., -1.]]) / 3
    np.testing.assert_allclose(rotation.T @ rotation, np.eye(3), atol=1e-15)
    assert np.linalg.det(rotation) == pytest.approx(1.)
    for push, masses in [(np.array([2., -.5, -1.]), np.array([.2, 3., 1.1])),
                         (np.array([-.4, -3., 1.5]), np.array([4., 1.6, 2.]))]:
        forces = masses[:, None] * push
        for basis in [np.eye(3), rotation]:
            accel = basis @ push
            reciprocal = accel / (accel @ accel)
            rotated_forces = forces @ basis.T
            np.testing.assert_allclose(rotated_forces @ reciprocal, masses, atol=2e-15)
            assert rotated_forces.sum(axis=0) @ reciprocal == pytest.approx(masses.sum())


@pytest.mark.parametrize("y", [[0., 0.], [2., .5], [-2.5, 1.]])
def test_zero_push_mass_changes_posterior_by_exactly_the_unresolved_share(y):
    weights = np.array([.3, .4, .3])
    pushes = np.array([1.25, -.8, 0.])
    covariance = np.array([[.8, -.3], [-.3, 1.3]])
    known_masses = np.array([1.5, 3.])
    outcomes, high_mass = [], []
    for hidden in [.5, 8.]:
        masses = np.r_[known_masses, hidden]
        locations = np.column_stack((masses * pushes, pushes))
        joint = weights * np.array([multivariate_normal.pdf(y, mean=mu, cov=covariance)
                                    for mu in locations])
        outcomes.append(joint.sum())
        high_mass.append(joint[masses > 2].sum() / joint.sum())
    assert outcomes[0] == outcomes[1]
    unresolved = weights[-1] * multivariate_normal.pdf(y, cov=covariance) / outcomes[0]
    assert 0 < unresolved < 1
    assert high_mass[1] - high_mass[0] == pytest.approx(unresolved)
    assert weights[:2][known_masses > 2].sum() == pytest.approx(.4)
    assert high_mass[0] > 0
    assert high_mass[1] < 1


def test_same_ensemble_mass_marginal_identifies_the_full_posterior():
    # The true five-component ensemble. Two masses occur under both active and
    # zero pushes. The last mass occurs exclusively at zero push.
    labels = np.array([0, 1, 0, 1, 2])
    masses = np.array([1., 3., 4.])
    pushes = np.array([2., -1., 0., 0., 0.])
    weights = np.array([.25, .35, .1, .05, .25])
    means = np.column_stack((masses[labels] * pushes, pushes))
    covariance = np.array([[1.2, -.3], [-.3, .7]])
    prior = np.bincount(labels, weights=weights)
    active_mass = np.array([.25, .35, 0.])
    zero_mass = prior - active_mass
    np.testing.assert_allclose(zero_mass, [.1, .05, .25])
    assert zero_mass.sum() == pytest.approx(.4)
    incompatible_prior = np.array([.24, .51, .25])
    assert np.min(incompatible_prior - active_mass) < 0
    for y in [[0., 0.], [-1., -.5], [2., 1.]]:
        density = np.array([multivariate_normal.pdf(y, mean=mu, cov=covariance)
                            for mu in means])
        joint_components = weights * density
        direct = np.bincount(labels, weights=joint_components) / joint_components.sum()
        active_joint = np.r_[active_mass[:2] * density[:2], 0.]
        recovered = active_joint + zero_mass * multivariate_normal.pdf(y, cov=covariance)
        np.testing.assert_allclose(recovered / recovered.sum(), direct, atol=2e-15)


def test_prior_correction_is_continuous_when_alternating_active_atoms_reach_zero():
    prior = np.array([.3, .3, .4])
    masses = np.array([1., 4., 11.])
    y = np.array([1., -.3])
    noise = multivariate_normal(cov=[[.8, .15], [.15, 1.]])
    baseline = noise.pdf(y)
    limit_joint = prior * [noise.pdf(y - [1., 1.]), baseline, baseline]
    limit = limit_joint / limit_joint.sum()

    for i in range(1, 8):
        scale = .1**i
        pushes = np.array([1., scale if i % 2 else 0., 0. if i % 2 else scale])
        means = np.column_stack((masses * pushes, pushes))
        active = pushes != 0
        # The null mass oscillates between .3 and .4 and jumps to .7 in the
        # weak limit. Only the complete posterior should converge to its limit.
        assert prior[~active].sum() == pytest.approx(.4 if i % 2 else .3)
        joint = prior * noise.pdf(y - means)
        corrected = baseline * prior
        for mu, weight in zip(means[active], prior[active]):
            mass = mu[0] / mu[1]
            index = np.argmin(np.abs(masses - mass))
            assert masses[index] == pytest.approx(mass)
            corrected[index] += weight * (noise.pdf(y - mu) - baseline)
        np.testing.assert_allclose(corrected, joint, atol=1e-15)
        posterior = corrected / joint.sum()
        assert np.max(np.abs(posterior - limit)) < 10 * scale


@pytest.mark.parametrize("scale", [1., .1, .001])
def test_distinct_nonzero_mass_experiments_get_arbitrarily_close_at_weak_drive(scale):
    covariance = np.array([[1.1, -.4], [-.4, .7]])
    means = scale * np.array([[1., 1.], [4., 1.]])
    difference = means[1] - means[0]
    # The set favouring the second mean is a half-space in this coordinate.
    direction = np.linalg.solve(covariance, difference)
    distance = np.sqrt(difference @ direction)
    midpoint = means.mean(axis=0)
    centers = (means - midpoint) @ direction / distance
    exact_tv, _ = quad(lambda z: norm.pdf(z - centers[1]) - norm.pdf(z - centers[0]),
                       0., np.inf)
    assert exact_tv == pytest.approx(erf(distance / (2 * np.sqrt(2))), abs=1e-13)
    assert 0 < exact_tv <= distance / np.sqrt(2 * np.pi)
    if scale == .001:
        assert exact_tv < .002


def test_total_only_force_loses_labels_even_when_pushes_are_all_nonzero():
    pushes = np.array([1., -.5])
    before = np.array([[1., 4.], [2., 1.]])
    after = np.array([[2., 3.], [.5, 2.5]])
    forces = [m * pushes[:, None] for m in (before, after)]
    labelled = [np.column_stack((f, pushes)) for f in forces]
    total = [np.column_stack((f.sum(axis=1), pushes)) for f in forces]
    np.testing.assert_array_equal(total[0], total[1])
    assert not np.array_equal(labelled[0], labelled[1])
    np.testing.assert_array_equal(total[0][:, 0] / pushes, [5., 3.])
