"""The sign numerator (O13, 29 Oct; Claude's sort, provisional): g = min(1, kappa_k / Q_+) with
kappa_k = k - 1/2, as `ours29(..., positivity=True)`. The default (ours29) is unchanged.

Q_+ is the misfit at the best positive mass (L21: "misfit at the best mass is evidence about the
reading"; L1: m > 0). kappa_k is the most Q_+ the declared noise produces in expectation over every
allowed signal size (the ruled meaning of the numerator, O12, read size-free under L10)."""
import math

import numpy as np
import pytest

import met.full as full

LM = np.linspace(-12.0, 12.0, 4001)


def _closed_q_plus(x, y):
    """T32: Q_+ = lambda_- when C > 0, min(A, B) otherwise (whitened, isotropic channels)."""
    a, b, c = x @ x, y @ y, x @ y
    if c > 0:
        return (a + b - math.sqrt((a - b) ** 2 + 4 * c * c)) / 2
    return min(a, b)


@pytest.mark.parametrize("k", [1, 2, 3])
def test_q_plus_is_t32s_closed_form(k):
    rng = np.random.default_rng(290 + k)
    for _ in range(40):
        f, a = rng.normal(size=k) * rng.uniform(0.1, 5), rng.normal(size=k) * rng.uniform(0.1, 5)
        sf, sa = rng.uniform(0.2, 2), rng.uniform(0.2, 2)
        for c in (0.3, sf / sa, 4.0):
            st = full.reading_state(full.isotropic(f, a, sf, sa), c)
            assert st["Q_plus"] == pytest.approx(_closed_q_plus(f / sf, a / sa), rel=1e-9, abs=1e-10)
            assert st["Q_plus"] >= st["Q"] - 1e-12


def test_q_plus_under_a_general_covariance_against_a_dense_grid():
    """Independent of the closed form: E - max of H over a dense grid of positive masses."""
    rng = np.random.default_rng(2910)
    for d in (2, 3):
        for _ in range(15):
            m = rng.normal(size=(2 * d, 2 * d))
            cov = m @ m.T + 0.3 * np.eye(2 * d)
            r = full.Reading(rng.normal(size=d) * 3, rng.normal(size=d), cov)
            st = full.reading_state(r, 1.0)
            phi = np.linspace(0.0, math.pi / 2, 200001)
            brute = r.pre()["energy"] - full._h_terms(r.pre(), 1.0, phi).max()
            assert st["Q_plus"] == pytest.approx(max(brute, st["Q"]), rel=1e-7, abs=1e-8)


def _span(force, acceleration, sd, positivity, c=1.0):
    readings = [full.isotropic(force, acceleration, sd, sd)]
    st = full.reading_state(readings[0], c)
    g = [full.misfit_weight(st["Q_plus"], st["kappa"]) if positivity else full.misfit_weight(st["Q"], st["nu"])]
    eta = full.bartlett_eta(readings, c, g)[0]
    ev = full._log_law(readings, c, g, eta, LM) - full._log_prior(LM, c)
    return float(ev.max() - ev.min())


@pytest.mark.parametrize("force,acc", [([1.0, 0.0], [-1.0, 0.0]), ([1.0], [-1.0]), ([2.0, 0.0, 0.0], [-1.0, 0.0, 0.0])])
def test_t30_bounded_tilt_on_anti_parallel_pairs_under_the_candidate(force, acc):
    """T30's requirement, which ours29 fails (tests/test_sign_fork.py, strict xfail): bounded
    evidence on a fixed pair no positive mass fits, in 1D too (kappa_1 = 1/2, not 0)."""
    spans = [_span(force, acc, sd, True) for sd in (0.1, 0.03, 0.01)]
    assert max(spans) - min(spans) < 1e-6 * max(spans)
    assert 0 < spans[0] < 10
    assert _span(force, acc, 0.01, False) > 100 * spans[0]          # ours29: unbounded


def test_no_gate_and_aligned_1d_readings_keep_full_weight():
    """L21: g > 0 everywhere; a 1D reading is weighted only when its sign contradicts m > 0."""
    rng = np.random.default_rng(2911)
    for _ in range(200):
        k = int(rng.integers(1, 4))
        r = full.isotropic(rng.normal(size=k) * 3, rng.normal(size=k), 0.5, 0.5)
        g = full.ours29([r], 2.0, positivity=True, readout=False)["g"][0]
        assert 0 < g <= 1
    assert full.ours29([full.isotropic([2.0], [1.0], 0.1, 0.1)], 1.0, positivity=True, readout=False)["g"] == [1.0]
    assert full.ours29([full.isotropic([2.0], [-1.0], 0.1, 0.1)], 1.0, positivity=True, readout=False)["g"][0] == pytest.approx(0.5 / 100)


def test_exact_zero_readings_change_nothing():
    """Z on exact zeros (T24, T25): Q_+ = 0, g = 1, law and eta unchanged."""
    sig = [full.isotropic([2.0, 0.3, 0.0], [1.0, 0.1, 0.05], 0.3, 0.2), full.isotropic([-1.0, 0.2, 0.1], [0.9, 0.0, 0.0], 0.3, 0.2)]
    base = full.ours29(sig, 1.5, positivity=True, readout=False)
    pad = full.ours29(sig + [full.isotropic([0.0] * 3, [0.0] * 3, 0.3, 0.2)] * 5, 1.5, positivity=True, readout=False)
    assert pad["g"][2:] == [1.0] * 5
    assert pad["eta"] == pytest.approx(base["eta"], rel=1e-12)
    assert pad["median"] == pytest.approx(base["median"], rel=1e-9)


def test_swap_and_units():
    """The swap (force <-> acceleration, c -> 1/c) and a change of units leave every g alone."""
    rng = np.random.default_rng(2912)
    for _ in range(30):
        k = int(rng.integers(1, 4))
        f, a, sf, sa = rng.normal(size=k) * 2, rng.normal(size=k), 0.4, 0.3
        g = full.ours29([full.isotropic(f, a, sf, sa)], 1.7, positivity=True, readout=False)["g"][0]
        swap = full.ours29([full.isotropic(a, f, sa, sf)], 1 / 1.7, positivity=True, readout=False)["g"][0]
        units = full.ours29([full.isotropic(5 * f, a, 5 * sf, sa)], 5 * 1.7, positivity=True, readout=False)["g"][0]
        assert swap == pytest.approx(g, rel=1e-9) and units == pytest.approx(g, rel=1e-9)


@pytest.mark.parametrize("k", [2, 3])
def test_acute_pairs_are_discounted_less_by_the_factor_k_minus_half_over_k_minus_1(k):
    """The cost of counting the sign: on an acute pair Q_+ = Q, so g moves from (k-1)/Q to
    (k-1/2)/Q wherever both are below 1. This changes the ruled nu = k - 1 (R-2026-09-29-9)."""
    f = np.zeros(k); a = np.zeros(k); f[0] = 2.0; a[0] = math.cos(0.4); a[1] = math.sin(0.4)
    r = full.isotropic(f, a, 0.1, 0.1)
    new = full.ours29([r], 1.0, positivity=True, readout=False)
    old = full.ours29([r], 1.0, readout=False)
    assert new["Q_plus"] == pytest.approx(old["Q"], rel=1e-12)
    assert new["g"][0] / old["g"][0] == pytest.approx((k - 0.5) / (k - 1), rel=1e-12)


@pytest.mark.parametrize("k", [1, 2, 3])
def test_the_numerator_is_the_supremum_over_sizes(k):
    """Seeded Monte Carlo over whitened channel strengths x = R e1 + xi, y = b e1 + zeta:
    sup E[Q] = k - 1 (both channels strong: the ruled nu is this supremum for the sign-blind
    misfit), sup E[Q_+] = k - 1/2 (one channel strong, the other without signal). Proof owed: B02."""
    rng = np.random.default_rng(2913 + k)
    n = 100_000
    xi, ze = rng.standard_normal((n, k)), rng.standard_normal((n, k))

    def moments(R, b):
        x, y = xi.copy(), ze.copy()
        x[:, 0] += R
        y[:, 0] += b
        A, B, C = (x * x).sum(1), (y * y).sum(1), (x * y).sum(1)
        q = (A + B - np.sqrt((A - B) ** 2 + 4 * C * C)) / 2
        return q.mean(), np.where(C > 0, q, np.minimum(A, B)).mean()

    grid = (0.0, 0.3, 1.0, 3.0, 10.0, 100.0)
    table = {(R, b): moments(R, b) for R in grid for b in grid}
    tol = 0.03 + 0.01 * k
    assert table[(100.0, 100.0)][0] == pytest.approx(k - 1, abs=tol)
    assert table[(100.0, 0.0)][1] == pytest.approx(k - 0.5, abs=tol)
    assert max(v[0] for v in table.values()) <= k - 1 + tol
    assert max(v[1] for v in table.values()) <= k - 0.5 + tol
