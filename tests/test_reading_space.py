"""T36: the reading space, its metric, and one residual (definition/theorems/T36.md).

(1)-(3) tie the geometry to met/full.py: H_i is the projection onto the mass subspace in the
noise metric, K_i(m) = E_i - H_i(m) is the residual, Q_i its minimum. (4) is the dimension
count: k_i per reading at the true mass, and the per-reading plug-in (26 Oct, O14) against one
residual at a common mass. Seeded; whitened isotropic channels for (4).
"""
import math

import numpy as np
import pytest

from met.full import Reading, _h_terms, reading_state


def _random_reading(rng, k, correlated):
    a = rng.normal(size=(2 * k, 2 * k))
    cov = a @ a.T + 0.5 * np.eye(2 * k) if correlated else np.diag(rng.uniform(0.2, 2.0, 2 * k))
    return Reading(rng.normal(size=k) * 3, rng.normal(size=k) * 2, cov)


def _closed_form_residual(reading, m):
    d = reading.dimension
    p = np.hstack([np.eye(d), -m * np.eye(d)])
    r = reading.force - m * reading.acceleration
    return float(r @ np.linalg.solve(p @ reading.cov @ p.T, r))


@pytest.mark.parametrize("k", [1, 2, 3])
@pytest.mark.parametrize("correlated", [False, True])
def test_t36_misfit_at_m_is_the_residual_in_the_noise_metric(k, correlated):
    rng = np.random.default_rng(36 + k + 10 * correlated)
    c = 1.3
    for _ in range(20):
        rd = _random_reading(rng, k, correlated)
        pre = rd.pre()
        phi = rng.uniform(-math.pi / 2 + 0.05, math.pi / 2 - 0.05, size=7)
        k_met = pre["energy"] - _h_terms(pre, c, phi)
        k_closed = np.array([_closed_form_residual(rd, c * math.tan(p)) for p in phi])
        assert np.allclose(k_met, k_closed, rtol=1e-9, atol=1e-9)


@pytest.mark.parametrize("k", [1, 2, 3])
def test_t36_own_line_misfit_is_the_least_residual(k):
    rng = np.random.default_rng(360 + k)
    c = 0.8
    for _ in range(20):
        rd = _random_reading(rng, k, correlated=True)
        q = reading_state(rd, c)["Q"]
        phi = np.linspace(-math.pi / 2, math.pi / 2, 20001)[1:]
        k_grid = rd.pre()["energy"] - _h_terms(rd.pre(), c, phi)
        assert q <= k_grid.min() + 1e-9
        assert q >= k_grid.min() - 1e-6 * max(1.0, k_grid.min())


def _series(rng, n_series, n, k, r, m=1.0):
    """Whitened isotropic readings, unit noise: Y[s, i] is the 2 x k data matrix (force; acceleration)."""
    u = rng.normal(size=(n_series, n, k))
    u /= np.linalg.norm(u, axis=-1, keepdims=True)
    a = r * u / math.sqrt(1 + m * m)
    f = m * a + rng.normal(size=a.shape)
    acc = a + rng.normal(size=a.shape)
    return np.stack([f, acc], axis=2), m


@pytest.mark.parametrize("k", [1, 2, 3])
def test_t36_residual_has_k_dimensions_per_reading_at_the_true_mass(k):
    rng = np.random.default_rng(3600 + k)
    for r in (0.0, 1.0, 5.0):
        y, m = _series(rng, 1, 20000, k, r)
        resid = ((y[0, :, 0] - m * y[0, :, 1]) ** 2).sum(-1) / (1 + m * m)   # K_i(m0), unit noise
        assert abs(resid.mean() / k - 1) < 0.03, (k, r, resid.mean() / k)


def _least_eig(g):
    t = g[..., 0, 0] + g[..., 1, 1]
    det = g[..., 0, 0] * g[..., 1, 1] - g[..., 0, 1] ** 2
    return (t - np.sqrt(np.maximum(t * t - 4 * det, 0))) / 2


def test_t36_one_residual_against_the_per_reading_plug_in():
    """3D, no signal, true scale 1: the plug-in sum Q / sum(k - 1) sits at 1/2 for every N
    (Codex's exact null value, O14); one residual at the series' best line over (Nk - 1)
    is far closer and approaches 1 as N grows. In 1D the plug-in is undefined."""
    rng = np.random.default_rng(29102026)
    k = 3
    out = {}
    for n in (5, 100):
        y, _ = _series(rng, 400, n, k, 0.0)
        gram = np.einsum("sipk,siqk->sipq", y, y)                 # per reading, 2 x 2
        plug = _least_eig(gram).sum(1) / (n * (k - 1))
        one = _least_eig(gram.sum(1)) / (n * k - 1)
        out[n] = (plug.mean(), one.mean())
    for n, (plug, one) in out.items():
        assert abs(plug - 0.5) < 0.05, (n, plug)
        assert one > plug + 0.15, (n, plug, one)
    assert out[100][1] > out[5][1] + 0.1
    assert out[100][1] > 0.88
    y, _ = _series(rng, 400, 20, 1, 2.0)
    one_1d = _least_eig(np.einsum("sipk,siqk->spq", y, y)) / (20 - 1)
    assert 0.9 < one_1d.mean() < 1.05
