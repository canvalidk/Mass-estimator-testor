"""B21 / T69: maximal invariant likelihood, its scope and its premise failures.

Research comparisons only. Strict xfails are the requested Z and T30 failures;
L10 and L19's permission for a shared Delta is not decided by a test. The
separately seeded coverage preset reports property N, without selecting a law.
"""
import math

import numpy as np
import pytest
from scipy.integrate import quad, cumulative_trapezoid
from scipy.optimize import minimize_scalar
from scipy.special import hyp0f1, i0

from tools.maximal_invariant import (angular_law, gram_energy,
                                     integrated_log_kernel, invariant_log_kernel,
                                     profile_curvature_ratio, profile_delta,
                                     profile_log_kernel)


@pytest.mark.parametrize("d", [1, 2, 3, 12])
def test_wishart_tilt_matches_independent_sphere_integrals(d):
    h, delta = 3.2, 1.7
    x = math.sqrt(h*delta)
    if d == 1:
        sphere = math.cosh(x)
    elif d == 2:
        sphere = quad(lambda theta: math.exp(x*math.cos(theta)), 0, 2*math.pi,
                      epsabs=1e-12)[0]/(2*math.pi)
    else:
        numerator = quad(lambda t: math.exp(x*t)*(1-t*t)**((d-3)/2), -1, 1,
                         epsabs=1e-12)[0]
        denominator = quad(lambda t: (1-t*t)**((d-3)/2), -1, 1, epsabs=1e-12)[0]
        sphere = numerator/denominator
    assert math.exp(float(invariant_log_kernel(h, d, delta))) == pytest.approx(
        math.exp(-delta/2)*sphere, rel=2e-11)


@pytest.mark.parametrize("d", [1, 2, 3, 15, 600])
def test_scaled_kernel_matches_direct_hypergeometric_in_overlap(d):
    h = np.array([0., .1, 4., 120., 300.])
    delta = np.array([3., 0., .4, 120., 500.])
    expected = -.5*delta+np.log(hyp0f1(d/2, delta*h/4))
    np.testing.assert_allclose(invariant_log_kernel(h, d, delta), expected,
                               atol=1e-12, rtol=2e-12)


@pytest.mark.parametrize("d", [1, 2, 3, 15, 600])
def test_profile_matches_independent_scalar_optimization(d):
    h = np.array([0., d/2, d, d+0.1, 2*d, 6*d+10])
    fitted = profile_delta(h, d)
    assert np.array_equal(fitted[:3], np.zeros(3))
    for hh, dd in zip(h[3:], fitted[3:]):
        fit = minimize_scalar(lambda energy: -float(invariant_log_kernel(hh, d, energy)),
                              bounds=(0, hh), method="bounded", options={"xatol": 1e-8})
        assert dd == pytest.approx(fit.x, rel=2e-5, abs=2e-5)
        assert float(profile_log_kernel(hh, d)) == pytest.approx(-fit.fun, abs=2e-10)


@pytest.mark.parametrize("d", [1, 2, 3, 15])
def test_integrated_kernel_matches_delta_quadrature(d):
    for h in [0., .3, 4., 20.]:
        # Integrate the original nuisance likelihood, not the 1F1 expression.
        answer = quad(lambda delta: math.exp(float(invariant_log_kernel(h, d, delta))),
                      0, np.inf, epsabs=1e-9, epsrel=1e-10)[0]
        assert math.log(answer/2) == pytest.approx(float(integrated_log_kernel(h, d)), abs=2e-10)


def test_integrated_special_cases_and_large_values():
    h = np.array([0., .01, 3., 500., 5000., 100000.])
    np.testing.assert_array_equal(integrated_log_kernel(h, 2), h/2)
    for d in [1, 3, 15, 600]:
        values = integrated_log_kernel(h, d)
        assert np.all(np.isfinite(values)) and np.all(np.diff(values) > 0)
    # In d=3 the half-line integral has its independent error-function form.
    from scipy.special import erf
    x = 40.
    exact = x+math.log(math.sqrt(math.pi)/(2*math.sqrt(x))*erf(math.sqrt(x)))
    assert float(integrated_log_kernel(2*x, 3)) == pytest.approx(exact, abs=1e-12)


@pytest.mark.parametrize("d", [1, 3, 15, 600])
def test_profile_curvature_is_the_envelope_derivative(d):
    h = float(3*d+7)
    step = h*1e-4
    derivative = (profile_log_kernel(h+step, d)-profile_log_kernel(h-step, d))/(2*step)
    assert float(2*derivative) == pytest.approx(float(profile_curvature_ratio(h, d)), rel=1e-7)
    assert float(profile_curvature_ratio(d, d)) == 0.


@pytest.mark.parametrize("signal_per_dimension", [.1, 1., 4.])
def test_profile_tempering_reaches_the_many_reading_target(signal_per_dimension):
    ds = np.array([30, 300, 3000])
    r = signal_per_dimension
    ratios = np.array([profile_curvature_ratio(d*(1+r), int(d)) for d in ds])
    target = r/(1+r)
    assert abs(ratios[-1]-target) < 1e-3
    assert abs(ratios[-1]-target) < abs(ratios[0]-target)/20
    # eta_G with lambda_min=d has the asymptotic target exactly; profile differs at finite d.
    assert abs(ratios[0]-target) > 1e-4


@pytest.mark.parametrize("d", [1, 3, 15])
def test_fixed_row_count_large_signal_recovers_flat_curvature(d):
    ratios = profile_curvature_ratio(np.array([100., 1000., 10000.]), d)
    assert ratios[-1] > .998
    assert 1-ratios[-1] <= (1-ratios[0])*.02+1e-12


@pytest.mark.parametrize("k,n", [(1, 7), (2, 3), (3, 4)])
def test_row_group_preserves_gram_and_model(k, n):
    rng = np.random.default_rng(2100+k)
    d = k*n
    z = rng.normal(size=(d, 2))
    rotation = np.linalg.qr(rng.normal(size=(d, d)))[0]
    np.testing.assert_allclose((rotation@z).T@(rotation@z), z.T@z, atol=2e-14)
    alpha = rng.normal(size=d)
    v = np.array([2., 1.])/math.sqrt(5)
    np.testing.assert_allclose(rotation@np.outer(alpha, v), np.outer(rotation@alpha, v), atol=1e-15)
    assert np.linalg.norm(rotation@alpha) == pytest.approx(np.linalg.norm(alpha), rel=1e-14)
    # Independent force/acceleration rotations do not in general preserve scalar-mass candidates.
    wrong = np.column_stack([rotation@alpha*v[0], alpha*v[1]])
    assert np.linalg.norm(wrong[:, 0]-2*wrong[:, 1]) > .1


def test_gram_maximality_reconstructs_the_row_orbit():
    rng = np.random.default_rng(2169)
    z = rng.normal(size=(8, 2))
    q = np.linalg.qr(rng.normal(size=(8, 8)))[0]
    other = q@z
    # Orthogonal Procrustes constructs a row transformation from equal Grams.
    u, _, vt = np.linalg.svd(other@z.T)
    recovered = u@vt
    np.testing.assert_allclose(recovered@z, other, atol=5e-15)
    np.testing.assert_allclose(recovered.T@recovered, np.eye(8), atol=2e-15)


def test_row_group_classes_follow_noise_ratios_not_absolute_scales():
    truth = 3.
    sf, sa = np.array([2., 6., 4., 12.]), np.array([1., 3., 1., 3.])
    acceleration = np.array([.7, -1.2, .9, 1.1])
    whitened = np.column_stack([truth*acceleration/sf, acceleration/sa])
    angle = .63
    rotation = np.array([[math.cos(angle), -math.sin(angle)],
                         [math.sin(angle), math.cos(angle)]])
    blocks = np.zeros((4, 4))
    blocks[:2, :2] = rotation
    blocks[2:, 2:] = rotation
    mixed = blocks@whitened
    # Absolute covariances differ within each block, but c=sf/sa does not.
    np.testing.assert_allclose(mixed[:, 0]*sf, truth*mixed[:, 1]*sa, atol=1e-15)
    np.testing.assert_allclose(mixed[:2].T@mixed[:2], whitened[:2].T@whitened[:2], atol=1e-15)
    cross = np.eye(4)
    cross[np.ix_([0, 2], [0, 2])] = rotation
    crossed = cross@whitened
    # Mixing distinct ratios changes the candidate line of the addressed rows.
    assert np.linalg.norm(crossed[:, 0]*sf-truth*crossed[:, 1]*sa) > .1


def test_wishart_sampling_moments_use_one_total_energy():
    rng = np.random.default_rng(2190)
    d, delta, batches = 6, 9., 30000
    v = np.array([2., 1.])/math.sqrt(5)
    z = rng.normal(size=(batches, d, 2))
    z[:, 0, :] += math.sqrt(delta)*v
    grams = np.einsum("nki,nkj->nij", z, z)
    np.testing.assert_allclose(grams.mean(0), d*np.eye(2)+delta*np.outer(v, v), atol=.08)
    assert grams[:, 0, 1].var() == pytest.approx(d+delta, rel=.03)
    assert grams[:, 0, 0].var() == pytest.approx(2*d+4*delta*v[0]**2, rel=.03)


@pytest.mark.parametrize("kind", ["profile", "integrated"])
def test_l2_and_all_zero_measured_inputs(kind):
    rng = np.random.default_rng(216)
    cases = [np.zeros((3, 2)), np.array([[2., 1.]]), np.array([[0., 1.]]),
             np.array([[1., 0.]]), np.array([[-2., 1.]]), rng.normal(size=(12, 2))]
    for z in cases:
        theta, density, log_kernel = angular_law(z.T@z, len(z), kind, 1001)
        assert np.all(np.isfinite(log_kernel)) and np.all(density > 0)
        assert np.trapezoid(density, theta) == pytest.approx(1., abs=2e-15)
        cdf = cumulative_trapezoid(density, theta, initial=0)
        median_theta = np.interp(.5, cdf, theta)
        assert 0 < median_theta < math.pi/2
    np.testing.assert_allclose(angular_law(np.zeros((2, 2)), 3, kind, 101)[1],
                               2/math.pi, atol=2e-15)


@pytest.mark.parametrize("kind", ["profile", "integrated"])
def test_l3_units_and_l4_swap_of_complete_laws(kind):
    rng = np.random.default_rng(2193)
    force, acceleration = rng.normal(size=(2, 12))
    sf, sa = .3, .7
    z = np.column_stack([force/sf, acceleration/sa])
    changed = np.column_stack([5*force/(5*sf), .2*acceleration/(.2*sa)])
    np.testing.assert_allclose(changed, z, atol=1e-15)
    theta, law, kernel = angular_law(z.T@z, len(z), kind, 1001)
    _, scaled, _ = angular_law(changed.T@changed, len(z), kind, 1001)
    _, swapped, swapped_kernel = angular_law((z[:, ::-1]).T@z[:, ::-1], len(z), kind, 1001)
    np.testing.assert_allclose(scaled, law, rtol=1e-11, atol=1e-13)
    np.testing.assert_allclose(swapped[::-1], law, rtol=1e-11, atol=1e-13)
    np.testing.assert_allclose(swapped_kernel[::-1], kernel, rtol=1e-11, atol=1e-12)
    # c scales by lambda/mu, so identical angular laws give correspondingly scaled masses.
    c = sf/sa
    ms = c*np.tan(theta[1:-1])
    np.testing.assert_allclose((sa/sf)*np.tan(theta[-2:0:-1]), 1/ms, rtol=1e-12)


@pytest.mark.parametrize("kind", ["profile", "integrated"])
@pytest.mark.xfail(strict=True, raises=AssertionError, reason="B21 / T69: exact-zero padding changes d=Nk while preserving G; both invariant laws fail Z.")
def test_z_requires_unchanged_law_after_appending_zero_readings(kind):
    gram = np.array([[14., 6.], [6., 5.]])
    first = angular_law(gram, 3, kind, 801)[1]
    padded = angular_law(gram, 18, kind, 801)[1]
    np.testing.assert_allclose(first, padded, rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize("kind", ["profile", "integrated"])
@pytest.mark.parametrize("pair", ["misaligned", "anti_1d", "anti_3d"])
@pytest.mark.xfail(strict=True, raises=AssertionError, reason="B21 / T69: invariant evidence spans grow as SD^-2 on fixed incompatible pairs; T30 fails.")
def test_t30_requires_bounded_tilt_as_declared_noise_shrinks(kind, pair):
    f, a = {"misaligned": ([math.sqrt(3), 1., 0.], [1., 0., 0.]),
            "anti_1d": ([-2.], [1.]), "anti_3d": ([-2., 0., 0.], [1., 0., 0.])}[pair]
    z = np.column_stack([f, a])
    theta = np.linspace(0, math.pi/2, 501)
    fn = profile_log_kernel if kind == "profile" else integrated_log_kernel
    spans = []
    for sd in [.1, .01]:
        kernel = fn(gram_energy(z.T@z/sd**2, theta), len(z))
        spans.append(float(np.ptp(kernel)))
    assert spans[1] <= 2*spans[0]


@pytest.mark.parametrize("k,signal", [(1, 1.), (1, 2.), (3, 1.), (3, 2.)])
def test_l11_fixed_signal_population_optimum_is_the_truth(k, signal):
    truth, n = 2., 180000
    v = np.array([truth, 1.])/math.sqrt(1+truth**2)
    expected = n*k*np.eye(2)+n*signal**2*np.outer(v, v)
    top = np.linalg.eigh(expected)[1][:, -1]
    assert top[0]/top[1] == pytest.approx(truth, abs=1e-14)
    # Seeded finite-world check beside the exact population identity.
    rng = np.random.default_rng(21000+k+int(signal))
    z = rng.normal(size=(n, k, 2))
    z[:, 0] += signal*v
    gram = np.einsum("nki,nkj->ij", z, z)
    fitted = np.linalg.eigh(gram)[1][:, -1]
    assert fitted[0]/fitted[1] == pytest.approx(truth, abs=.06)
    # For every positive Delta both laws are increasing in the same h.
    angles = np.arctan([1.8, truth, 2.2])
    h = gram_energy(expected/n, angles)
    assert np.argmax(profile_log_kernel(h, k)) == 1
    assert np.argmax(integrated_log_kernel(h, k)) == 1


def test_mile_exceptions_are_flat_or_boundary_not_an_interior_tls_point():
    theta = np.linspace(0., math.pi/2, 501)
    # Every angle maximizes the profile below the fitted-energy threshold.
    np.testing.assert_array_equal(profile_log_kernel(gram_energy(np.diag([1., 2.]), theta), 3),
                                  np.zeros_like(theta))
    # Anti-aligned Gram: the positive-domain supremum is at its excluded endpoint.
    gram = np.array([[16., -8.], [-8., 4.]])
    kernel = profile_log_kernel(gram_energy(gram, theta), 3)
    assert np.argmax(kernel) == len(theta)-1
    assert np.all(kernel[:-1] < kernel[-1])
