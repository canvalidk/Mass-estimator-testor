"""B11: the formal integrated-scale kernel, including its propriety failures.

These are test-local calculations, not a production law. A proper v integral
at each fixed m does not ensure a proper marginal law of m.
"""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

import met.full as full


def _residual(readings, mass):
    phi = np.array([math.atan(mass)])
    total = 0.0
    for reading in readings:
        y = np.r_[reading.force, reading.acceleration]
        energy = y @ np.linalg.solve(reading.cov, y)
        total += energy - full._h_terms(reading.pre(), 1.0, phi)[0]
    return total


def _count(readings):
    return sum(2 * r.dimension - r.k for r in readings)


def _scalar(force, acceleration):
    return full.Reading([force], [acceleration], np.eye(2))


@pytest.mark.parametrize("d,k", [(1, 1), (2, 2), (3, 3), (2, 1), (3, 1), (3, 2)])
def test_correlated_geometry_and_scale_integral_include_outside_subspace(d, k):
    rng = np.random.default_rng(1100 + 10 * d + k)
    a = rng.normal(size=(2 * d, 2 * d))
    cov = a @ a.T + np.eye(2 * d)
    subspace = np.linalg.qr(rng.normal(size=(d, d)))[0][:, :k]
    y = rng.normal(size=2 * d)
    reading = full.Reading(y[:d], y[d:], cov, subspace)
    precision = np.linalg.inv(cov)
    assert _count([reading]) == 2 * d - k
    for mass in (0.6, 2.3):
        candidate = np.vstack([mass * subspace, subspace])
        gram = candidate.T @ precision @ candidate
        projector = candidate @ np.linalg.solve(gram, candidate.T @ precision)
        residual = y - projector @ y
        direct = residual @ precision @ residual
        actual = _residual([reading], mass)
        assert actual == pytest.approx(direct, abs=2e-14)
        assert np.trace(np.eye(2 * d) - projector) == pytest.approx(2 * d - k)
        if k < d:
            # pre.energy is only the span energy; using it loses pure noise.
            span_residual = reading.pre()["energy"] - full._h_terms(
                reading.pre(), 1.0, np.array([math.atan(mass)]))[0]
            assert actual > span_residual + 1e-5
        beta = (2 * d - k) / 2  # a=1, the branch specified in B11.
        closed = math.gamma(beta) * (actual / 2) ** -beta
        numerical = quad(lambda v: v ** (-beta - 1) * math.exp(-actual / (2 * v)),
                         0, np.inf, epsabs=1e-10, epsrel=2e-10, limit=250)[0]
        assert numerical == pytest.approx(closed, rel=2e-9)


@pytest.mark.parametrize("copies", [1, 5])
def test_aligned_positive_readings_have_a_nonintegrable_mass_pole(copies):
    readings = [_scalar(2, 1) for _ in range(copies)]
    beta = copies / 2
    assert _count(readings) == copies
    for mass in (0.5, 1.5, 2.0, 3.0):
        exact = copies * (2 - mass) ** 2 / (1 + mass * mass)
        assert _residual(readings, mass) == pytest.approx(exact, abs=2e-14)
    # On disjoint shells m=2+[epsilon,2epsilon], the unnormalised density
    # is copies^-beta (1+m²)^(beta-1) |m-2|^-copies (pi's 2/pi omitted).
    # Each shell has mass bounded away from zero for copies=1, or growing
    # as epsilon^-4 for copies=5. No improper-normaliser quadrature is used.
    shells = []
    for epsilon in (1e-2, 1e-3, 1e-4):
        shell = copies ** -beta * epsilon ** (1 - copies) * quad(
            lambda t: (1 + (2 + epsilon * t) ** 2) ** (beta - 1) * t ** -copies,
            1, 2, epsabs=1e-12)[0]
        shells.append(shell)
        if copies == 1:
            assert shell >= math.log(2) / math.sqrt(1 + 2.02 ** 2)
        else:
            assert shell >= (15 / 320) * epsilon ** -4
    assert min(shells) > 0.3
    if copies == 5:
        assert shells[-1] > shells[0] * 9e7
    ours = full.ours30(readings, 1.0, readout=False)
    assert np.isfinite(ours["law"]).all()
    assert np.sum(ours["law"]) == pytest.approx(1, abs=3e-14)
    assert 0 < ours["interval"][0] < ours["median"] < ours["interval"][1] < np.inf


@pytest.mark.parametrize("force,acceleration,boundary", [(0, 1, 1e-6), (1, 0, 1e6)])
def test_finite_pointwise_scale_integral_has_an_improper_mass_boundary(force, acceleration, boundary):
    reading = _scalar(force, acceleration)
    exact = lambda mass: (force - mass * acceleration) ** 2 / (1 + mass * mass)
    for mass in (0.2, 1, 4):
        assert _residual([reading], mass) == pytest.approx(exact(mass), abs=2e-14)
    for mass in (1e-6, 1, 1e6):
        residual = exact(mass)  # Avoid cancellation in energy-H near the boundary.
        assert residual > 0
        assert math.isfinite(math.gamma(0.5) * (residual / 2) ** -0.5)
    # pi*K^-1/2 behaves as 2/(pi*m): near zero for (0,1), and at infinity
    # for (1,0). Either logarithmic divergence prevents normalisation,
    # despite the v integral being finite at every positive mass.
    density = (2 / math.pi) / (1 + boundary ** 2) / math.sqrt(exact(boundary))
    assert boundary * density == pytest.approx(2 / math.pi, abs=1e-11)


@pytest.mark.parametrize("copies", [1, 2, 5])
def test_all_zero_scale_integral_diverges_while_ours30_returns_prior(copies):
    readings = [_scalar(0, 0) for _ in range(copies)]
    assert [_residual(readings, mass) for mass in (0.2, 1, 4)] == [0, 0, 0]
    beta = _count(readings) / 2
    # K=0 leaves v^(-beta-1): its [epsilon,2epsilon] mass diverges.
    shell = lambda epsilon: epsilon ** -beta * (1 - 2 ** -beta) / beta
    assert shell(1e-4) / shell(1e-3) == pytest.approx(10 ** beta)
    ours = full.ours30(readings, 1.0, readout=False)
    assert np.isfinite(ours["law"]).all()
    assert np.sum(ours["law"]) == pytest.approx(1, abs=3e-14)
    assert ours["median"] == pytest.approx(1, abs=2e-13)


@pytest.mark.parametrize("ratio", [0.5, 2.0])
def test_distinct_proper_paths_to_all_zero_have_distinct_limiting_laws(ratio):
    # Two readings (epsilon*s,0), (0,epsilon*t) have beta=1 and
    # pi*K^-1 proportional to 1/(s²+t²m²), for every epsilon>0.
    s, t = ratio, 1.0
    for epsilon in (1e-1, 1e-3):
        readings = [_scalar(epsilon * s, 0), _scalar(0, epsilon * t)]
        assert _count(readings) == 2
        for mass in (0.2, 1, 4):
            expected = epsilon ** 2 * (s * s + t * t * mass * mass) / (1 + mass * mass)
            assert _residual(readings, mass) == pytest.approx(expected, rel=2e-14)
        # The normalised density is half-Cauchy with centre s/t. Thus a
        # path-independent continuation at all-zero is not supplied by K.
        density = lambda mass: 2 * s * t / (math.pi * (s * s + t * t * mass * mass))
        assert quad(density, 0, np.inf)[0] == pytest.approx(1, abs=2e-13)
        assert quad(density, 0, ratio)[0] == pytest.approx(0.5, abs=2e-13)


@pytest.mark.parametrize("zeros,expected", [
    (0, 1.5811388300841898), (1, 1.6099334758590365),
    (5, 1.6180126120954565), (20, 1.6180339887498763),
])
def test_null_padding_changes_the_median_of_a_proper_series(zeros, expected):
    base = [_scalar(2, 1), _scalar(1, 1)]
    padded = base + [_scalar(0, 0) for _ in range(zeros)]
    assert _count(padded) == 2 + zeros
    for mass in (0.2, 1, 2, 4):
        assert _residual(padded, mass) == _residual(base, mass)
    beta = _count(padded) / 2
    minimum = (7 - math.sqrt(45)) / 2
    assert minimum > 0

    def angular(phi):
        residual = (2 * math.cos(phi) - math.sin(phi)) ** 2
        residual += (math.cos(phi) - math.sin(phi)) ** 2
        return (residual / minimum) ** -beta

    normalizer = quad(angular, 0, math.pi / 2, epsabs=1e-12, epsrel=1e-12)[0]
    phi = brentq(lambda x: quad(angular, 0, x, epsabs=1e-12)[0] - normalizer / 2,
                 0, math.pi / 2, xtol=1e-14)
    assert math.tan(phi) == pytest.approx(expected, abs=2e-10)

    def direct(mass):
        residual = ((2 - mass) ** 2 + (1 - mass) ** 2) / (1 + mass * mass)
        return (residual / minimum) ** -beta / (1 + mass * mass)

    mass_normalizer = quad(direct, 0, np.inf, epsabs=1e-12, epsrel=1e-12)[0]
    assert mass_normalizer == pytest.approx(normalizer, abs=2e-12)
    assert quad(direct, 0, expected, epsabs=1e-12)[0] / mass_normalizer == pytest.approx(0.5, abs=2e-11)
