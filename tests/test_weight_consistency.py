"""G22 (31 Oct, Claude): which per-reading weights keep the law's population optimum at the truth.

The tempered log-law is sum_i g_i H_i(theta) with g_i fixed per reading, so as N grows at fixed
per-reading signal its maximiser tends to the maximiser of E[g H(theta)]
= sin^2 E[gA] + cos^2 E[gB] + 2 sin cos E[gC] (whitened, c = 1). A weight that is a function of the
Gram eigenvalues alone (ours29's g on the whole-line misfit) keeps that maximiser at the true angle:
given the eigenvalues the eigenvector angle's density is even about it (T38's proof, step 1).
ours30's g on Q_+ depends on which quadrant the eigenvector falls in, is correlated with the score,
and moves the maximiser (B05, T23: the literal L11 counterexample). Seeded population checks.
"""
import os

import numpy as np
import pytest

from met.full import isotropic, ours29, ours30


def _population_mass(k, r, weight, m=2.0, n=400_000, seed=2210):
    rng = np.random.default_rng(seed)
    d = r * np.array([m, 1.0]) / np.sqrt(1 + m * m)
    x = rng.normal(size=(n, k)); x[:, 0] += d[0]
    y = rng.normal(size=(n, k)); y[:, 0] += d[1]
    A, B, C = (x * x).sum(1), (y * y).sum(1), (x * y).sum(1)
    lam = (A + B) / 2 - np.sqrt(((A - B) / 2) ** 2 + C ** 2)
    if weight == "ours30":
        q = np.where(C > 0, lam, np.minimum(A, B))
        g = np.minimum(1.0, (k - 0.5) / np.maximum(q, 1e-300))
    else:
        g = np.ones_like(lam) if k == 1 else np.minimum(1.0, (k - 1) / np.maximum(lam, 1e-300))
    s = np.array([[(g * A).mean(), (g * C).mean()], [(g * C).mean(), (g * B).mean()]])
    u = np.linalg.eigh(s)[1][:, 1]
    return u[0] / u[1]


@pytest.mark.parametrize("k", [2, 3])
def test_g22_eigenvalue_only_weights_keep_the_population_optimum_at_the_truth(k):
    # at r = 0.5 the eigen-gap is small and 4e5 draws scatter by about 0.05 across seeds, so the
    # check starts at r = 1 (scatter about 0.015); with 2e6 draws r = 0.5 gives 2.005 (2D), 2.002 (3D)
    for r in (1.0, 2.0):
        assert _population_mass(k, r, "ours29") == pytest.approx(2.0, abs=0.05)


@pytest.mark.parametrize("k, r, upper", [(1, 1.0, 1.80), (2, 1.0, 1.85), (3, 1.0, 1.86), (3, 2.0, 1.98)])
def test_g22_ours30_weights_move_the_population_optimum_below_the_truth(k, r, upper):
    # recorded 31 Oct with 2e6 draws: 1D 1.735, 2D 1.790, 3D 1.813 at r = 1; 3D 1.957 at r = 2
    assert _population_mass(k, r, "ours30") < upper


@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason="3000 readings through met.full (about 70 s)")
def test_g22_production_ours30_concentrates_at_the_population_optimum_not_the_truth():
    # recorded 31 Oct: ours30 median 1.7638, [1.5970, 1.9586]; ours29 2.0019, [1.7685, 2.2887]
    # 1D, per-reading signal 1, true mass 2, c = s = 1: B05's world with a larger signal.
    rng = np.random.default_rng(3110)
    mu = np.array([2.0, 1.0]) / np.sqrt(5)
    z = rng.normal(size=(3000, 2)) + mu
    readings = [isotropic([f], [a], 1.0, 1.0) for f, a in z]
    p30 = ours30(readings, c=1.0)
    p29 = ours29(readings, c=1.0)
    lo30, hi30 = p30["interval"]
    assert hi30 < 2.0, "ours30's 95% interval should exclude the true mass here"
    lo29, hi29 = p29["interval"]
    assert lo29 <= 2.0 <= hi29
