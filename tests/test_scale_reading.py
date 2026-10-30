"""What g_i and eta are: noise scales read on the floor branch (problems/O14.md, 30 Oct sort).

1. Tempering a reading's evidence by w is declaring its covariance Sigma / w (H is homogeneous
   of degree -1 in Sigma): g_i and eta are noise scales, 1/g_i and 1/eta, never below 1.
2. T30's bounded confidence is a consequence of reading the scale from the misfit: once g < 1,
   g * K_i(m) does not change when the declared Sigma is scaled, so sigma -> 0 cannot sharpen it.
3. On worlds where the declaration is off, g and eta do different jobs: eta finds a common
   inflation (alone in 1D, where g is inert); g finds one bad reading in 3D; on a common
   inflation in 3D they over-correct together. Seeded; MET_SLOW=1 gives the 150-series numbers
   quoted in O14.
"""
import math
import os

import numpy as np
import pytest

from met.full import Reading, _h_terms, isotropic, ours30, reading_state


def _reading(rng, k):
    a = rng.normal(size=(2 * k, 2 * k))
    return Reading(rng.normal(size=k) * 3, rng.normal(size=k) * 2, a @ a.T + 0.5 * np.eye(2 * k))


@pytest.mark.parametrize("k", [1, 2, 3])
def test_tempering_a_reading_is_declaring_its_noise_larger(k):
    rng = np.random.default_rng(1030 + k)
    phi = np.linspace(-1.4, 1.4, 11)
    for _ in range(10):
        r = _reading(rng, k)
        w = rng.uniform(0.05, 0.95)
        wider = Reading(r.force, r.acceleration, r.cov / w)
        assert np.allclose(w * _h_terms(r.pre(), 1.0, phi), _h_terms(wider.pre(), 1.0, phi), rtol=1e-10)
        assert w * r.pre()["energy"] == pytest.approx(wider.pre()["energy"], rel=1e-10)


@pytest.mark.parametrize("k", [2, 3])
def test_t30_bounded_confidence_follows_from_reading_the_scale_from_the_misfit(k):
    """A fixed misaligned pair (30 degrees), declared noise shrinking: once g < 1 the tempered
    residual g * (E - H) is the same function of m at every declared scale."""
    th = math.radians(30)
    f = np.zeros(k); f[0], f[1] = 2 * math.cos(th), 2 * math.sin(th)
    a = np.zeros(k); a[0] = 1.0
    phi = np.linspace(0.2, 1.3, 9)
    curves = []
    for sd in (1e-2, 1e-4, 1e-6):
        r = isotropic(f, a, sd, sd)
        st = reading_state(r, 1.0)
        g = min(1.0, st["kappa"] / st["Q_plus"])
        assert g < 1
        curves.append(g * (r.pre()["energy"] - _h_terms(r.pre(), 1.0, phi)))
    for c in curves[1:]:
        assert np.allclose(c, curves[0], rtol=1e-6)


def _series(rng, n, k, r, v0, outlier):
    out = []
    for i in range(n):
        u = rng.normal(size=k)
        u /= np.linalg.norm(u)
        a = r * u / math.sqrt(2)
        sd = math.sqrt(v0 * (25.0 if i == outlier else 1.0))
        out.append(isotropic(a + sd * rng.normal(size=k), a + sd * rng.normal(size=k), 1.0, 1.0))
    return out


def _cell(seed, k, v0, outlier, n_series):
    rng = np.random.default_rng(seed)
    cover, ge = [], []
    for _ in range(n_series):
        res = ours30(_series(rng, 40, k, 20.0, v0, outlier), 1.0, readout=False)
        lo, hi = res["interval"]
        cover.append(lo <= 1.0 <= hi)
        ge.append(float(np.median(res["g"]) * res["eta"]))
    return float(np.mean(cover)), float(np.median(ge))


def test_g_and_eta_on_worlds_where_the_declaration_is_off():
    n = 150 if os.environ.get("MET_SLOW") else 20
    cov1, ge1 = _cell(301, 1, 4.0, None, n)      # common inflation, 1D: eta alone
    cov3, ge3 = _cell(303, 3, 4.0, None, n)      # common inflation, 3D: g and eta together
    covo, geo = _cell(305, 3, 1.0, 0, n)         # one bad reading, 3D: g catches it
    assert 0.2 < ge1 < 0.34                      # 1/v0 = 0.25
    assert ge3 < 0.22                            # over-corrected: below 1/v0
    assert geo > 0.9                             # eta stays near 1
    if os.environ.get("MET_SLOW"):
        print(f"\n1D x4: cover {cov1:.2f} g*eta {ge1:.3f} | 3D x4: cover {cov3:.2f} g*eta {ge3:.3f}"
              f" | 3D one bad: cover {covo:.2f} g*eta {geo:.3f}")
