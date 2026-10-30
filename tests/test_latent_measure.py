"""The latent measure and calibration (problems/O15.md, 30 Oct, later): a declared size law on the
pushes restores calibration; the flat measure does not; eta recovers part of it.

Bayesian coverage: m ~ pi (half-Cauchy, c = 1), pushes alpha_i ~ N(0, w^2 I) in whitened units, unit
noise, 3D, N = 20. The structural law integrates the pushes under the true size law. Seeded; MET_SLOW=1
gives the 300-series numbers quoted in O15.
"""
import math
import os

import numpy as np
import pytest

from met.full import _h_terms, _log_prior, _quantiles, isotropic, ours30

LM = np.linspace(math.log(1e-3), math.log(1e3), 6001)
M = np.exp(LM)
PHI = np.arctan(M)


def _interval(logf):
    (lo, _, hi), _ = _quantiles(LM, _log_prior(LM, 1.0) + logf, [0.025, 0.5, 0.975])
    return lo, hi


def structural_log_factor(readings, w):
    """Each push integrated under alpha ~ N(0, w^2 I): y ~ N(0, I + w^2 C_m C_m^T), whitened, c = 1."""
    lf = np.zeros_like(LM)
    for r in readings:
        s = w * w * (1 + M * M)
        lf += -(r.k / 2) * np.log1p(s) + (s / (2 * (1 + s))) * _h_terms(r.pre(), 1.0, PHI)
    return lf


def test_structural_law_is_the_flat_law_with_lindleys_pull_as_w_grows():
    rng = np.random.default_rng(5)
    rd = [isotropic(rng.normal(size=3) + [2, 0, 0], rng.normal(size=3) + [1, 0, 0], 1.0, 1.0) for _ in range(4)]
    big = structural_log_factor(rd, 1e4)
    flat = sum(_h_terms(r.pre(), 1.0, PHI) / 2 for r in rd)
    pull = -(3 * len(rd) / 2) * np.log1p(M * M)          # (1 + m^2)^(-Nk/2): Lindley 1966's factor
    diff = big - flat - pull
    assert np.ptp(diff) < 1e-3                            # equal up to a constant


def _coverage(w, n_series, seed):
    rng = np.random.default_rng([3113, seed])
    out = {"flat": [], "ours30": [], "structural": []}
    for _ in range(n_series):
        m = abs(rng.standard_cauchy())
        rd = []
        for _ in range(20):
            a = w * rng.normal(size=3)
            rd.append(isotropic(m * a + rng.normal(size=3), a + rng.normal(size=3), 1.0, 1.0))
        lo, hi = _interval(sum(_h_terms(r.pre(), 1.0, PHI) / 2 for r in rd))
        out["flat"].append(lo <= m <= hi)
        lo, hi = ours30(rd, 1.0, readout=False)["interval"]
        out["ours30"].append(lo <= m <= hi)
        lo, hi = _interval(structural_log_factor(rd, w))
        out["structural"].append(lo <= m <= hi)
    return {k: float(np.mean(v)) for k, v in out.items()}


def test_a_declared_size_law_restores_calibration_and_the_flat_measure_does_not():
    n = 300 if os.environ.get("MET_SLOW") else 40
    c = _coverage(0.3, n, 0)
    assert c["flat"] < 0.80
    assert c["structural"] > 0.88
    assert c["flat"] < c["ours30"] < c["structural"] + 0.05
    if os.environ.get("MET_SLOW"):
        print(f"\nw = 0.3: flat {c['flat']:.3f} | ours30 {c['ours30']:.3f} | structural {c['structural']:.3f}")
