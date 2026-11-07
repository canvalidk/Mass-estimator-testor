"""T70: the reading as a qubit (Claude, 7 Nov; Can's question of the same day).

In noise units on one isotropic instrument (c = s), a reading is the qubit psi = (|f|, |alpha| e^{i theta}),
theta the angle between force and acceleration; a mass hypothesis is the real state (sin phi, cos phi).
"""
import math

import numpy as np
import pytest

from met.full import _h_terms, bartlett_eta, gram_eta, isotropic, reading_state

SF, SA = 0.3, 0.1
C = SF / SA
PHIS = np.linspace(-1.5, 1.5, 13)


def _qubit(reading):
    f, al = reading.force / SF, reading.acceleration / SA
    nf, na = float(np.linalg.norm(f)), float(np.linalg.norm(al))
    if nf == 0 or na == 0:
        return np.array([nf, na], dtype=complex)
    th = math.acos(float(np.clip(f @ al / (nf * na), -1.0, 1.0)))
    return np.array([nf, na * np.exp(1j * th)])


def _gram(readings):
    g = np.zeros((2, 2))
    for r in readings:
        f, al = r.force / SF, r.acceleration / SA
        g += np.array([[f @ f, f @ al], [f @ al, al @ al]])
    return g


def _readings(seed, n, k):
    rng = np.random.default_rng(seed)
    return [isotropic(rng.normal(size=k) * 2.0, rng.normal(size=k), SF, SA) for _ in range(n)]


@pytest.mark.parametrize("k", [1, 2, 3])
def test_t70_explained_energy_is_the_born_rule(k):
    for r in _readings(70 + k, 4, k):
        psi = _qubit(r)
        born = np.array([abs(math.sin(p) * psi[0] + math.cos(p) * psi[1]) ** 2 for p in PHIS])
        np.testing.assert_allclose(_h_terms(r.pre(), C, PHIS), born, rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("k", [1, 2, 3])
def test_t70_whole_line_misfit_is_the_distance_to_the_nearest_real_state(k):
    """Q = E lambda_min(Re rho) = (E - |sum psi_j^2|) / 2."""
    for r in _readings(80 + k, 4, k):
        psi = _qubit(r)
        e = float(np.sum(np.abs(psi) ** 2))
        lam_min = float(np.linalg.eigvalsh(np.real(np.outer(psi, psi.conj())))[0])
        q = reading_state(r, C)["Q"]
        assert q == pytest.approx(lam_min, rel=1e-8, abs=1e-9)
        assert q == pytest.approx((e - abs(np.sum(psi ** 2))) / 2, rel=1e-8, abs=1e-9)


def test_t70_the_anti_parallel_pair_is_real_so_the_phase_cannot_see_the_sign():
    """theta = pi gives a real state with a negative amplitude: zero misfit on the whole line.
    Positivity (L1) is a sign condition, not a phase condition (Q_plus sees it, Q does not)."""
    r = isotropic([1.0, 0.0, 0.0], [-0.5, 0.0, 0.0], SF, SA)
    st = reading_state(r, C)
    assert st["Q"] == pytest.approx(0.0, abs=1e-9)
    assert st["Q_plus"] > 1.0


def test_t70_the_default_starting_law_is_uniform_over_real_states():
    """Half-Cauchy at c = s, pushed to phi = arctan(m / c): density 2 / pi on (0, pi/2)."""
    m = np.geomspace(1e-3, 1e3, 9) * C
    phi_density = (2 / (math.pi * C * (1 + (m / C) ** 2))) * C / np.cos(np.arctan(m / C)) ** 2
    np.testing.assert_allclose(phi_density, 2 / math.pi, rtol=1e-12)


def test_t70_the_flat_series_law_is_an_axial_von_mises():
    """sum_i H_i(phi) = E/2 + (R/2) cos(2 phi - 2 phi_hat): the doubled angle, concentration R/4 in exp(H/2)."""
    rs = _readings(90, 5, 3)
    g = _gram(rs)
    e = float(np.trace(g))
    bz, bx = g[1, 1] - g[0, 0], 2 * g[0, 1]
    r_len, ang = math.hypot(bx, bz), math.atan2(bx, bz)
    total = sum(_h_terms(r.pre(), C, PHIS) for r in rs)
    np.testing.assert_allclose(total, e / 2 + (r_len / 2) * np.cos(2 * PHIS - ang), rtol=1e-10)


def test_t70_one_misaligned_reading_and_two_disagreeing_aligned_readings_give_the_same_flat_law():
    """Every real 2x2 Gram is the real part of one pure qubit, so the flat law cannot tell them apart."""
    pair = [isotropic([2.0 * SF * 3], [SA * 3], SF, SA), isotropic([0.5 * SF * 2], [SA * 2], SF, SA)]
    g = _gram(pair)
    nf, na = math.sqrt(g[0, 0]), math.sqrt(g[1, 1])
    th = math.acos(g[0, 1] / (nf * na))
    one = isotropic([nf * SF, 0.0], [na * SA * math.cos(th), na * SA * math.sin(th)], SF, SA)
    np.testing.assert_allclose(_gram([one]), g, rtol=1e-12)
    np.testing.assert_allclose(_h_terms(one.pre(), C, PHIS),
                               sum(_h_terms(r.pre(), C, PHIS) for r in pair), rtol=1e-10)
    assert reading_state(one, C)["Q"] > 0                       # the single reading is misaligned
    assert all(reading_state(r, C)["Q"] < 1e-9 for r in pair)   # each 1D reading is not


def test_t70_a_gram_only_tempering_cannot_bound_a_misaligned_reading():
    """(6): along sigma -> 0 a fixed misaligned pair's Gram scales with fixed shape, so a tempering that
    reads only the Gram's shape (as eta_G does) leaves its evidence growing like sigma^-2; T30 needs a
    count of what noise alone produces, which ours30's per-reading g supplies."""
    f, a = [2.0, 0.0, 0.0], [math.cos(0.5), math.sin(0.5), 0.0]
    spans, bounded = [], []
    for sd in (0.1, 0.03, 0.01):
        r = isotropic(f, a, sd, sd)
        eta = gram_eta([r], 1.0, [1.0])[0]
        h = _h_terms(r.pre(), 1.0, np.linspace(0.05, 1.5, 400))
        spans.append(eta * (h.max() - h.min()) / 2)
        g = min(1.0, reading_state(r, 1.0)["kappa"] / reading_state(r, 1.0)["Q_plus"])
        bounded.append(g * (h.max() - h.min()) / 2)
    assert max(bounded) - min(bounded) < 1e-3 * max(bounded)    # with g: frozen (T30)
    assert spans[2] > 50 * spans[0]                              # Gram shape only: grows as sd^-2


def test_t70_gram_eta_is_the_eigenvalue_ratio_on_one_instrument():
    rs = _readings(100, 6, 3)
    lam = np.linalg.eigvalsh(_gram(rs))
    assert gram_eta(rs, C, [1.0] * len(rs))[0] == pytest.approx(1 - lam[0] / lam[1], rel=1e-8)


def test_t70_gram_eta_leaves_a_zero_reading_alone_where_the_count_does_not():
    rs = _readings(101, 5, 3)
    padded = rs + [isotropic([0.0] * 3, [0.0] * 3, SF, SA)]
    assert gram_eta(padded, C, [1.0] * 6)[0] == pytest.approx(gram_eta(rs, C, [1.0] * 5)[0], rel=1e-10)
    m = sum(r.k for r in rs)
    h = gram_eta(rs, C, [1.0] * 5)[2]
    assert 1 - (m + 3) / h < 1 - m / h                           # the size law's 1 - M/H^2 moves (O14)


def test_t70_gram_eta_does_not_depend_on_units_swap_or_frame():
    rs = _readings(102, 4, 3)
    base = gram_eta(rs, C, [1.0] * 4)[0]
    lam, mu = 3.0, 0.5
    scaled = [isotropic(lam * r.force, mu * r.acceleration, lam * SF, mu * SA) for r in rs]
    assert gram_eta(scaled, C * lam / mu, [1.0] * 4)[0] == pytest.approx(base, rel=1e-8)
    swapped = [isotropic(r.acceleration, r.force, SA, SF) for r in rs]
    assert gram_eta(swapped, 1 / C, [1.0] * 4)[0] == pytest.approx(base, rel=1e-8)
    q, _ = np.linalg.qr(np.random.default_rng(5).normal(size=(3, 3)))
    rotated = [isotropic(q @ r.force, q @ r.acceleration, SF, SA) for r in rs]
    assert gram_eta(rotated, C, [1.0] * 4)[0] == pytest.approx(base, rel=1e-8)
    assert gram_eta(rs, 7.0 * C, [1.0] * 4)[0] == pytest.approx(base, rel=1e-6)   # not on c


def test_t70_bartlett_eta_and_gram_eta_agree_when_nothing_disagrees():
    """An exactly aligned noiseless series: no misfit, so eta_G = 1; Bartlett's J = 0 gives 1 too."""
    rs = [isotropic([2.0 * s * SF], [s * SA], SF, SA) for s in (1.0, 2.0, 3.0)]
    assert gram_eta(rs, C, [1.0] * 3)[0] == pytest.approx(1.0, abs=1e-9)
    assert bartlett_eta(rs, C, [1.0] * 3)[0] == pytest.approx(1.0, abs=1e-9)
