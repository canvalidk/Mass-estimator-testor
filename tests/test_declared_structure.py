"""E6: independent checks of declared-constraint Gaussian profile statements."""

import numpy as np
import pytest
from scipy.linalg import block_diag
from scipy.optimize import minimize_scalar
from scipy.special import hyp1f1

from met.declared_structure import gaussian_profile, isotropic_profile, linear_mass_profile
from met.full import Reading, _h_terms


@pytest.mark.parametrize("dimension", [1, 2, 3])
@pytest.mark.parametrize("mass", [0.2, 1.0, 4.0])
def test_closed_linear_profiles_match_independent_gls(dimension, mass):
    rng = np.random.default_rng(6100 + dimension)
    f, a = rng.normal(size=(2, 5, dimension))
    n = rng.normal(size=dimension)
    n /= np.linalg.norm(n)
    bases = {
        "free": np.eye(f.size),
        "known_line": block_diag(*([n[:, None]] * len(f))),
        "equal": np.tile(np.eye(dimension), (len(f), 1)),
    }
    for structure, basis in bases.items():
        expected = linear_mass_profile(f, a, mass, basis)
        actual = isotropic_profile(f, a, mass, structure, n)
        assert actual == pytest.approx(expected, abs=2e-12)


@pytest.mark.parametrize("mass", [0.1, 0.6, 1.5, 5.0])
def test_unknown_line_matches_direct_direction_search(mass):
    rng = np.random.default_rng(6106)
    f, a = rng.normal(size=(2, 6, 2))

    def cost(theta):
        n = np.array([np.cos(theta), np.sin(theta)])
        return isotropic_profile(f, a, mass, "known_line", n)

    # Independently fit the physical direction, not a singular value formula.
    grid = np.linspace(0, np.pi, 401)
    j = np.argmin([cost(t) for t in grid])
    fitted = minimize_scalar(cost, bounds=(grid[j] - 0.02, grid[j] + 0.02), method="bounded")
    assert isotropic_profile(f, a, mass, "unknown_line") == pytest.approx(fitted.fun, abs=1e-9)


def test_same_reading_grams_do_not_determine_common_unknown_line_evidence():
    f = np.array([[2.0, 0.0], [0.5, 0.0]])
    a = np.array([[1.0, 0.0], [1.0, 0.0]])
    fr, ar = f.copy(), a.copy()
    fr[1], ar[1] = [0.0, 0.5], [0.0, 1.0]
    for i in range(2):
        np.testing.assert_allclose(np.stack((f[i], a[i])) @ np.stack((f[i], a[i])).T,
                                   np.stack((fr[i], ar[i])) @ np.stack((fr[i], ar[i])).T)
    for m in [0.25, 1.0, 2.0, 4.0]:
        assert isotropic_profile(f, a, m) == pytest.approx(isotropic_profile(fr, ar, m))
    fits = [minimize_scalar(lambda m: isotropic_profile(x, y, m, "unknown_line"),
                            bounds=(0.1, 5), method="bounded").x
            for x, y in [(f, a), (fr, ar)]]
    assert fits[0] == pytest.approx((9 + np.sqrt(481)) / 20, abs=2e-6)
    assert fits[1] == pytest.approx(2.0, abs=2e-6)


def test_unknown_common_line_is_zero_safe_despite_sharing_a_nuisance():
    rng = np.random.default_rng(6107)
    f, a = rng.normal(size=(2, 4, 3))
    for m in [0.1, 0.8, 2.0, 10.0]:
        original = isotropic_profile(f, a, m, "unknown_line")
        padded = isotropic_profile(np.vstack((f, np.zeros(3))), np.vstack((a, np.zeros(3))),
                                   m, "unknown_line")
        assert padded == pytest.approx(original, abs=1e-12)


def test_zero_independent_observation_of_a_shared_vector_changes_mass_contrast():
    f, a = np.array([[2.0]]), np.array([[1.0]])
    fp, ap = np.vstack((f, [0.0])), np.vstack((a, [0.0]))
    before = isotropic_profile(f, a, 1.0, "equal") - isotropic_profile(f, a, 2.0, "equal")
    after = isotropic_profile(fp, ap, 1.0, "equal") - isotropic_profile(fp, ap, 2.0, "equal")
    assert before == pytest.approx(0.5)
    assert after == pytest.approx(before / 2)
    # It is not an m-independent normalizer difference.
    assert after != before


@pytest.mark.parametrize("mass", [0.2, 1.0, 3.0])
def test_equal_vector_reduces_to_mean_with_covariance_and_ancillary_contrasts(mass):
    rng = np.random.default_rng(6108)
    f, a = rng.normal(size=(2, 6, 3))
    mean_pair = np.r_[f.mean(axis=0), a.mean(axis=0)]
    design = np.vstack((mass * np.eye(3), np.eye(3)))
    contrasts = np.sum((f - f.mean(axis=0))**2) + np.sum((a - a.mean(axis=0))**2)
    assert isotropic_profile(f, a, mass, "equal") == pytest.approx(
        contrasts + gaussian_profile(mean_pair, np.eye(6) / len(f), design))


def test_averaging_free_readings_loses_mass_evidence():
    # The contrast is itself a mass reading when signals are free.
    f, a = np.array([[2.0], [-2.0]]), np.array([[1.0], [-1.0]])
    assert np.all(f.mean(axis=0) == 0) and np.all(a.mean(axis=0) == 0)
    assert isotropic_profile(f, a, 2.0) == pytest.approx(0)
    assert isotropic_profile(f, a, 1.0) == pytest.approx(1.0)


def test_full_experiment_recoding_preserves_profile_under_spd_correlated_noise():
    rng = np.random.default_rng(6109)
    f, a = rng.normal(size=(2, 3, 2))
    b = rng.normal(size=(6, 3))
    z = rng.normal(size=(12, 12))
    cov = z @ z.T + np.eye(12)
    transform = rng.normal(size=(12, 12)) + 4 * np.eye(12)
    y, design = np.r_[f.ravel(), a.ravel()], np.vstack((1.7 * b, b))
    base = gaussian_profile(y, cov, design)
    assert gaussian_profile(transform @ y, transform @ cov @ transform.T,
                            transform @ design) == pytest.approx(base, rel=1e-11)
    assert linear_mass_profile(f, a, 1.7, b @ np.diag([1, 2, 3]), cov) == pytest.approx(base)
    perm = rng.permutation(12)
    assert gaussian_profile(y[perm], cov[np.ix_(perm, perm)], design[perm]) == pytest.approx(base)


def test_independent_frame_changes_require_transporting_the_declared_relation():
    rng = np.random.default_rng(6110)
    f, a = rng.normal(size=(2, 3, 2))
    rotations = [np.linalg.qr(rng.normal(size=(2, 2)))[0] for _ in range(3)]
    r = block_diag(*rotations)
    b = np.tile(np.eye(2), (3, 1))
    base = linear_mass_profile(f, a, 1.4, b)
    fr, ar = (r @ f.ravel()).reshape(3, 2), (r @ a.ravel()).reshape(3, 2)
    assert linear_mass_profile(fr, ar, 1.4, r @ b) == pytest.approx(base)
    assert abs(linear_mass_profile(fr, ar, 1.4, b) - base) > 0.01


def test_product_constraint_profiles_factor_with_block_independent_noise():
    rng = np.random.default_rng(6111)
    readings = []
    expected = 0.0
    for dimension in [2, 3]:
        z = rng.normal(size=(2 * dimension, 2 * dimension))
        cov = z @ z.T + np.eye(2 * dimension)
        y = rng.normal(size=2 * dimension)
        b = rng.normal(size=(dimension, 1))
        x = np.vstack((1.3 * b, b))
        expected += gaussian_profile(y, cov, x)
        readings.append((y, cov, x))
    assert gaussian_profile(np.concatenate([r[0] for r in readings]),
                            block_diag(*[r[1] for r in readings]),
                            block_diag(*[r[2] for r in readings])) == pytest.approx(expected)


def test_free_profile_matches_current_full_type_geometry():
    rng = np.random.default_rng(6112)
    f, a = rng.normal(size=(2, 4, 2))
    m = 1.6
    readings = [Reading(x, y, np.eye(4)) for x, y in zip(f, a)]
    explained = sum(_h_terms(r.pre(), 1.0, np.array([np.arctan(m)]))[0] for r in readings)
    assert isotropic_profile(f, a, m) == pytest.approx(np.sum(f*f + a*a) - explained)


def test_single_acceleration_acquisition_is_not_independent_copies():
    f, observed_a = np.array([1.0, 3.0]), 1.0
    for m in [0.5, 1.0, 2.0, 4.0]:
        correct = gaussian_profile(np.r_[f, observed_a], np.eye(3), np.array([[m], [m], [1.0]]))
        reduced = 2.0 + (2.0 - m)**2 / (0.5 + m*m)
        assert correct == pytest.approx(reduced)
        naive = isotropic_profile(f[:, None], np.full((2, 1), observed_a), m, "equal")
        assert naive == pytest.approx(2.0 + (2.0 - m)**2 / (0.5 + 0.5*m*m))
    assert gaussian_profile([1, 3, 1], np.eye(3), [[1], [1], [1]]) == pytest.approx(8/3)


def test_duplicate_reports_with_singular_covariance_retain_one_observation():
    y = np.array([1.0, 3.0, 1.0, 1.0])
    cov = block_diag(np.eye(2), np.ones((2, 2)))
    assert y[2] == y[3]  # deterministic support must accompany the pseudoinverse
    precision = np.linalg.pinv(cov)
    for m in [0.5, 1.0, 3.0]:
        x = np.array([m, m, 1.0, 1.0])
        u = (x @ precision @ y) / (x @ precision @ x)
        residual = y - x*u
        supported_cost = residual @ precision @ residual
        assert supported_cost == pytest.approx(gaussian_profile(y[:3], np.eye(3), x[:3, None]))
    with pytest.raises(np.linalg.LinAlgError):
        gaussian_profile(y, cov, np.ones((4, 1)))


def test_relative_factorization_requires_product_signal_constraints():
    def mixed_difference(basis):
        def contrast(y):
            f, a = y[:2, None], y[2:, None]
            return (linear_mass_profile(f, a, 1.0, basis) -
                    linear_mass_profile(f, a, 2.0, basis))
        u, v = np.eye(4)[0], np.eye(4)[1]
        return contrast(u+v) - contrast(u) - contrast(v) + contrast(np.zeros(4))

    assert mixed_difference(np.eye(2)) == pytest.approx(0, abs=1e-14)
    # Nonzero mixed finite difference cannot be removed by any data-only constant.
    assert mixed_difference(np.ones((2, 1))) == pytest.approx(0.3)


def test_ours28b_count_changes_relative_evidence_without_any_vector_coherent_sum():
    # Exact 1F1 factor in the historical equation: M counts one new zero.
    # This checks the formula independently, without running old readout code.
    h = lambda m: (2*m + 1)**2 / (1 + m*m)
    ratios = [hyp1f1(0.5, 1 + count/2, h(1)/2) /
              hyp1f1(0.5, 1 + count/2, h(2)/2) for count in [1, 2]]
    assert abs(ratios[0] - ratios[1]) > 1e-3


@pytest.mark.parametrize("structure", ["free", "known_line", "unknown_line", "equal"])
def test_profiles_bounded_by_observed_energy(structure):
    rng = np.random.default_rng(6113)
    f, a = rng.normal(size=(2, 4, 3))
    for m in np.geomspace(1e-6, 1e6, 11):
        k = isotropic_profile(f, a, m, structure, [1, 2, 3])
        assert 0 <= k <= np.sum(f*f + a*a) + 1e-12
