"""B16: independent checks of the evidence-only normalization obstruction.

These are proof checks, not a sampled coverage study. A finite integration
window cannot establish propriety on the entire positive mass line.
"""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

from met.full import Reading, isotropic
from met.ours31 import _per_reading
from met.ours33 import _l_terms


def _cases():
    yield []
    yield [isotropic([0.], [0.], 1., 1.)] * 3
    yield [isotropic([2.], [1.], 1., 1.)]
    yield [isotropic([2.], [1.], 1., 1.), isotropic([1.], [4.], 1., 1.)]
    angle = math.atan(2.)
    line = np.array([math.sin(angle), math.cos(angle)])
    normal = np.array([math.cos(angle), -math.sin(angle)])
    for radius in (100., 1000.):
        points = [radius * line + sign * 3. * normal for sign in (-1, 1)]
        yield [isotropic([p[0]], [p[1]], 1., 1.) for p in points]
    rng = np.random.default_rng(4)
    signal = np.array([1., 0., 0.])
    yield [isotropic(signal + rng.normal(size=3), signal + rng.normal(size=3), 1., 1.)
           for _ in range(6)]
    # Full cross-channel covariance and a proper declared spatial subspace.
    chol = np.array([[1., 0., 0., 0.], [.2, 1.4, 0., 0.],
                     [-.3, .1, .8, 0.], [.1, -.2, .3, 1.1]])
    yield [Reading([2., -.4], [.3, 1.], chol @ chol.T,
                   subspace=np.array([[.6], [.8]]))]


def _independent_projection(reading, mass):
    """QR/lstsq reference, independent of the estimator's cached H/K code.

    Endpoints 0 and infinity are the vertical and horizontal graph limits.
    Only the declared two-channel span enters the kernel under test.
    """
    n = reading.subspace
    zero = np.zeros_like(n)
    chol = np.linalg.cholesky(reading.cov)
    span = np.block([[n, zero], [zero, n]])
    basis, _ = np.linalg.qr(np.linalg.solve(chol, span), mode="reduced")
    y = basis.T @ np.linalg.solve(chol, np.r_[reading.force, reading.acceleration])
    candidate = np.vstack((n, zero)) if math.isinf(mass) else np.vstack((mass * n, n))
    design = basis.T @ np.linalg.solve(chol, candidate)
    fitted = design @ np.linalg.lstsq(design, y, rcond=None)[0]
    residual = y - fitted
    return float(y @ y), float(residual @ residual)


@pytest.mark.parametrize("readings", list(_cases()))
def test_every_finite_series_has_a_positive_uniform_evidence_floor(readings):
    c = .7
    nu = [.5 + 2. * index for index in range(len(readings))]
    states = _per_reading(readings, c, nu)
    energies = [_independent_projection(r, 0.)[0] for r in readings]
    floor = -sum((v + r.k) / 2 * math.log1p(e / v)
                 for r, v, e in zip(readings, nu, energies))
    masses = np.exp(np.linspace(-12., 12., 101))
    terms = _l_terms(states, c, np.arctan(masses / c), 0)
    actual = sum((t[0] for t in terms), np.zeros_like(masses))
    reference = np.array([-sum((v + r.k) / 2 * math.log1p(
        _independent_projection(r, m)[1] / v) for r, v in zip(readings, nu))
        for m in masses])
    np.testing.assert_allclose(actual, reference, rtol=2e-10, atol=2e-10)
    assert np.all(actual >= floor - 2e-10)
    assert np.all(actual <= 2e-10)
    assert math.exp(floor) > 0.  # These fixtures do not underflow the analytic bound.


@pytest.mark.parametrize("readings", list(_cases()))
def test_both_endpoint_graphs_leave_strictly_positive_evidence(readings):
    nu = [10.] * len(readings)
    states = _per_reading(readings, 1., nu)
    expected = []
    for endpoint in (0., math.inf):
        expected.append(-sum((v + r.k) / 2 * math.log1p(
            _independent_projection(r, endpoint)[1] / v)
            for r, v in zip(readings, nu)))
    terms = _l_terms(states, 1., np.array([0., math.pi / 2]), 0)
    observed = sum((t[0] for t in terms), np.zeros(2))
    np.testing.assert_allclose(observed, expected, rtol=1e-12, atol=1e-10)
    assert all(math.isfinite(x) and math.exp(x) > 0. for x in expected)


def _aligned_evidence(mass, power=1.):
    # One actual compatible reading (F,a)=(2,1), nu=10, k=1.
    # Independent analytic residual; no prior factor or estimator grid.
    residual = ((2. - mass) / math.hypot(1., mass)) ** 2
    return math.exp(-power * 5.5 * math.log1p(residual / 10.))


@pytest.mark.parametrize("power", [0., .1, 1., 3.])
def test_truncated_integrals_grow_in_mass_and_log_mass_under_every_named_power(power):
    floor = math.exp(-power * 5.5 * math.log1p(5. / 10.))
    tail_low, tail_high = _aligned_evidence(0., power), math.exp(-power * 5.5 * math.log1p(.1))
    mass_ratios, log_ratios = [], []
    for bound in (64., 256., 1024.):
        integral = quad(lambda m: _aligned_evidence(m, power), 0., bound,
                        points=[2.], epsabs=1e-9)[0]
        assert integral >= floor * bound - 1e-9
        mass_ratios.append(integral / bound)
    for bound in (16., 64., 128.):
        integral = quad(lambda u: _aligned_evidence(math.exp(u), power),
                        -bound, bound, points=[0., math.log(2.)], epsabs=1e-9)[0]
        assert integral >= floor * 2 * bound - 1e-9
        log_ratios.append(integral / (2 * bound))
    assert abs(mass_ratios[-1] - tail_high) < abs(mass_ratios[0] - tail_high) + 1e-12
    assert abs(log_ratios[-1] - (tail_low + tail_high) / 2) < \
        abs(log_ratios[0] - (tail_low + tail_high) / 2) + 1e-12


def test_finite_window_quantiles_follow_the_cutoff_instead_of_defining_evidence_width():
    quantiles = []
    for cutoff in (1000., 10000.):
        def integral(upper):
            return quad(_aligned_evidence, 0., upper, points=[2.] if upper > 2. else None)[0]
        normalizer = integral(cutoff)
        q = [brentq(lambda m: integral(m) / normalizer - p, 0., cutoff)
             for p in (.025, .5, .975)]
        np.testing.assert_allclose(np.array(q) / cutoff, [.025, .5, .975], atol=.006)
        quantiles.append(q)
    # With the same reading, the artificial 95% width grows with the window.
    assert (quantiles[1][2] - quantiles[1][0]) / (quantiles[0][2] - quantiles[0][0]) > 9.
