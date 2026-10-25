"""met.full: ours29 in the full input type (any covariance per reading, a declared subspace).

Three kinds of check:
  - agreement with the reference implementation (reference/noise_vs_signal/ours29_reference.py)
    at the reference's own grid, and with met's batch pipeline on isotropic noise;
  - the definition's exact properties: rotation with covariance (T22), the swap (T5), units (L3),
    Z for the law (T24), zero components (T25), the misfit as an eigenvalue (T32), totality (T13);
  - the exact branch (T1).
"""
import math
import sys
import pathlib

import numpy as np
import pytest

from met.full import Reading, isotropic, instrument_ratio, ours29, reading_state, exact_mass
from met.analyst import estimate

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "reference" / "noise_vs_signal"))
from ours29_reference import ours29 as reference_ours29  # noqa: E402

LAW = {"reference": "cartesian2", "nuisance": "radius"}
OURS29 = {"rule": "likelihood_product", "prior": [{"uniform_angle": {"center": "first_reading"}}],
          "misfit": "kink", "calibrate": "bartlett"}


def _rng(seed):
    return np.random.default_rng(seed)


def _full_cov(rng, d, scale=1.0):
    a = rng.normal(size=(2 * d, 2 * d))
    return scale * (a @ a.T / (2 * d) + 0.3 * np.eye(2 * d))


def _series_full(seed, d, n, mass=2.0, signal=1.5):
    """n readings in d dimensions, each with its own full covariance (cross-channel correlation
    included), a true mass and a random true push; errors drawn from the declared covariance."""
    rng = _rng(seed)
    out = []
    for _ in range(n):
        cov = _full_cov(rng, d)
        push = rng.normal(size=d) * signal
        truth = np.r_[mass * push, push]
        y = truth + np.linalg.cholesky(cov) @ rng.normal(size=2 * d)
        out.append(Reading(y[:d], y[d:], cov))
    return out


def _as_reference(readings):
    out = []
    for r in readings:
        y = np.r_[r.force, r.acceleration]
        if r.k < r.dimension:
            out.append((y, r.cov, r.subspace))
        else:
            out.append((y, r.cov))
    return out


def _rotate(reading, rot):
    block = np.kron(np.eye(2), rot)
    sub = None if reading.k == reading.dimension else rot @ reading.subspace
    return Reading(rot @ reading.force, rot @ reading.acceleration, block @ reading.cov @ block.T, sub)


def _swap(reading):
    d = reading.dimension
    perm = np.block([[np.zeros((d, d)), np.eye(d)], [np.eye(d), np.zeros((d, d))]])
    sub = None if reading.k == d else reading.subspace
    return Reading(reading.acceleration, reading.force, perm @ reading.cov @ perm.T, sub)


def _units(reading, lam):
    d = reading.dimension
    scale = np.diag(np.r_[np.full(d, lam), np.ones(d)])
    sub = None if reading.k == d else reading.subspace
    return Reading(lam * reading.force, reading.acceleration, scale @ reading.cov @ scale, sub)


def _random_rotation(rng, d):
    q, r = np.linalg.qr(rng.normal(size=(d, d)))
    return q * np.sign(np.diag(r))


# ------------------------------------------------------------------ agreement with the reference

@pytest.mark.parametrize("d,n,seed", [(1, 5, 11), (2, 4, 12), (3, 3, 13), (3, 8, 14)])
def test_agrees_with_the_reference_at_its_grid(d, n, seed):
    readings = _series_full(seed, d, n)
    mine = ours29(readings, 2.0, grid=(9.0, 721))
    ref = reference_ours29(_as_reference(readings), c=2.0)
    # the misfit weights and eta: the reference's best lines are golden-section and its
    # derivatives finite differences (h = 1e-4), so it is good to about 1e-7
    assert np.allclose(mine["g"], ref["g"], rtol=1e-6, atol=1e-9)
    assert abs(mine["eta"] - ref["eta"]) <= 1e-6
    # the readout: the reference truncates its E|X| integral at log t = +-20 (about 1e-8)
    assert mine["readout"] == pytest.approx(ref["readout"], rel=1e-6)
    # the median and interval: same grid, so the only difference is the quantile rule
    # (the reference's midpoint CDF against met's fourth-order one): about a grid step
    step = 18.0 / 720
    assert abs(math.log(mine["median"] / ref["median"])) <= step
    assert abs(math.log(mine["interval"][0] / ref["lo"])) <= step
    assert abs(math.log(mine["interval"][1] / ref["hi"])) <= step


def test_declared_subspace_agrees_with_the_reference():
    rng = _rng(21)
    cov = _full_cov(rng, 3)
    line = np.array([[0.2], [-0.35], [0.915]]) / np.linalg.norm([0.2, -0.35, 0.915])
    plane = np.linalg.qr(rng.normal(size=(3, 2)))[0]
    readings = [Reading([2.0, 1.0, 0.5], [0.9, 0.6, 0.1], cov),
                Reading([-1.0, 2.0, 1.2], [-0.3, 1.1, 0.7], cov, plane),
                Reading([0.4, -0.8, 1.9], [0.3, -0.2, 1.1], cov, line)]
    mine = ours29(readings, 2.0, grid=(9.0, 721))
    ref = reference_ours29(_as_reference(readings), c=2.0)
    assert np.allclose(mine["g"], ref["g"], rtol=1e-6)
    assert mine["readout"] == pytest.approx(ref["readout"], rel=1e-6)
    assert abs(math.log(mine["median"] / ref["median"])) <= 18.0 / 720


def test_eta_below_one_agrees_with_the_reference():
    force = [[4.0, 0.3, 0.0], [6.1, -0.2, 0.1], [1.0, 0.8, 0.0], [3.0, 0.1, -0.2]]
    acc = [[2.0, 0.0, 0.1], [2.0, 0.1, 0.0], [0.6, -0.4, 0.0], [1.0, 0.0, 0.0]]
    readings = [isotropic(f, a, 0.2, 0.1) for f, a in zip(force, acc)]
    mine = ours29(readings, 2.0)
    ref = reference_ours29(_as_reference(readings), c=2.0)
    assert mine["eta"] < 0.2
    assert mine["eta"] == pytest.approx(ref["eta"], abs=1e-7)
    assert mine["readout"] == pytest.approx(ref["readout"], rel=1e-8)


# ------------------------------------------------------------------ agreement with met's batch pipeline

@pytest.mark.parametrize("d,seed", [(1, 1), (2, 2), (3, 3)])
def test_isotropic_noise_agrees_with_the_batch_pipeline(d, seed):
    rng = _rng(seed)
    force = rng.normal(size=(6, d)) * 2 + 1
    acc = rng.normal(size=(6, d)) + 0.4
    mine = ours29([isotropic(f, a, 0.5, 0.25) for f, a in zip(force, acc)], 2.0)
    batch = estimate(force, acc, 0.5, 0.25, ["ratio_of_means_rss", "median", "interval_95"],
                     law=LAW, combine=OURS29)
    assert mine["readout"] == pytest.approx(float(batch["ratio_of_means_rss"]), rel=1e-9)
    assert mine["median"] == pytest.approx(float(batch["median"]), rel=1e-6)
    assert np.allclose(mine["interval"], batch["interval_95"], rtol=1e-6)


def test_both_steps_off_is_ours27():
    rng = _rng(4)
    force = rng.normal(size=(5, 3)) * 2 + 1
    acc = rng.normal(size=(5, 3)) + 0.4
    mine = ours29([isotropic(f, a, 0.5, 0.25) for f, a in zip(force, acc)], 2.0,
                  misfit=False, calibrate=False)
    batch = estimate(force, acc, 0.5, 0.25, ["ratio_of_means_rss", "median"], law=LAW,
                     combine={"rule": "likelihood_product",
                              "prior": [{"uniform_angle": {"center": "first_reading"}}]})
    assert mine["eta"] == 1.0 and mine["g"] == [1.0] * 5
    assert mine["readout"] == pytest.approx(float(batch["ratio_of_means_rss"]), rel=1e-9)
    assert mine["median"] == pytest.approx(float(batch["median"]), rel=1e-6)


# ------------------------------------------------------------------ exact properties

@pytest.mark.parametrize("d,seed", [(2, 31), (3, 32)])
def test_rotation_with_covariance_changes_nothing(d, seed):
    """T22: each reading rotated with its covariance (and its subspace)."""
    readings = _series_full(seed, d, 5, signal=0.6)
    rng = _rng(seed + 100)
    rotated = [_rotate(r, _random_rotation(rng, d)) for r in readings]
    a, b = ours29(readings, 2.0), ours29(rotated, 2.0)
    assert b["eta"] == pytest.approx(a["eta"], abs=1e-11)
    assert np.allclose(b["g"], a["g"], rtol=1e-10, atol=1e-12)
    assert b["readout"] == pytest.approx(a["readout"], rel=1e-10)
    assert b["median"] == pytest.approx(a["median"], rel=1e-10)
    assert np.allclose(b["interval"], a["interval"], rtol=1e-10)


def test_swap_gives_the_reciprocal():
    """T5 and L4: force and acceleration swapped, centre c -> 1/c: readout and median invert."""
    readings = _series_full(41, 3, 5, signal=0.8)
    a = ours29(readings, 2.0)
    b = ours29([_swap(r) for r in readings], 0.5)
    assert b["eta"] == pytest.approx(a["eta"], abs=1e-11)
    assert a["readout"] * b["readout"] == pytest.approx(1.0, rel=1e-9)
    assert a["median"] * b["median"] == pytest.approx(1.0, rel=1e-8)
    assert a["interval"][0] * b["interval"][1] == pytest.approx(1.0, rel=1e-7)


def test_units_scale_as_a_mass():
    """L3: force units times lambda, with the centre: every mass output times lambda."""
    readings = _series_full(42, 2, 4)
    a = ours29(readings, 2.0)
    b = ours29([_units(r, 7.0) for r in readings], 14.0)
    assert b["readout"] == pytest.approx(7.0 * a["readout"], rel=1e-10)
    assert b["median"] == pytest.approx(7.0 * a["median"], rel=1e-9)
    assert b["eta"] == pytest.approx(a["eta"], abs=1e-12)


def test_a_zero_reading_leaves_the_law_alone_and_moves_the_readout():
    """T24 (Z holds for the law, for any covariance) and T29 (the cloud readout fails Z)."""
    readings = _series_full(43, 3, 4, signal=0.7)
    rng = _rng(143)
    zero = Reading(np.zeros(3), np.zeros(3), _full_cov(rng, 3))
    a, b = ours29(readings, 2.0, grid=(12.0, 3001)), ours29(readings + [zero], 2.0, grid=(12.0, 3001))
    assert b["eta"] == pytest.approx(a["eta"], abs=1e-12)
    assert b["median"] == pytest.approx(a["median"], rel=1e-12)
    assert np.allclose(b["law"], a["law"], rtol=1e-10, atol=1e-300)
    assert abs(b["readout"] / a["readout"] - 1) > 1e-6


@pytest.mark.parametrize("cov_kind", ["isotropic", "full"])
def test_the_all_zero_experiment_returns_the_input_centre_as_its_median(cov_kind):
    """T24, TASK.md section 2: on null input the law is pi, so its median is c, for any covariance."""
    rng = _rng(44)
    cov = np.diag([1.0] * 3 + [0.25] * 3) if cov_kind == "isotropic" else _full_cov(rng, 3)
    zeros = [Reading(np.zeros(3), np.zeros(3), cov) for _ in range(3)]
    for c in (0.3, 2.0, 6.0):
        r = ours29(zeros, c)
        assert r["median"] == pytest.approx(c, rel=1e-6)
        assert r["eta"] == 1.0


def test_the_all_zero_readout_is_the_centre_only_when_the_centre_is_s():
    """T14 holds only at c = s; at c = 3s the cloud readout gives 2.095 s (ours26 sweep, E1)."""
    zeros = [isotropic(np.zeros(3), np.zeros(3), 1.0, 0.5) for _ in range(2)]
    s = instrument_ratio(zeros[0])
    assert s == pytest.approx(2.0)
    assert ours29(zeros, s)["readout"] == pytest.approx(s, rel=1e-7)
    assert ours29(zeros, 3 * s)["readout"] / s == pytest.approx(2.095, abs=5e-4)


def test_zero_components_carry_nothing_and_a_k1_reading_is_the_1d_reading():
    """T25: exact zeros on some axes leave the law as without those axes; k = 1 equals 1D."""
    one_d = [isotropic([4.0], [2.0], 1.0, 1.0)]
    three_d = [isotropic([4.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0, 1.0)]
    k1 = [Reading([4.0, 0.0, 0.0], [2.0, 0.0, 0.0], np.eye(6), np.array([[1.0], [0.0], [0.0]]))]
    a, b, c = (ours29(r, 1.0, grid=(9.0, 1801), readout=False) for r in (one_d, three_d, k1))
    assert np.allclose(b["law"], a["law"], rtol=1e-10)
    assert np.allclose(c["law"], a["law"], rtol=1e-10)


@pytest.mark.parametrize("seed", [51, 52, 53])
def test_the_misfit_is_the_smaller_eigenvalue_of_the_whitened_gram_matrix(seed):
    """T32: with diagonal channel noise, Q is lambda_- of the reading's whitened 2 x 2 Gram matrix."""
    rng = _rng(seed)
    for d in (2, 3):
        f, a = rng.normal(size=d) * 3, rng.normal(size=d)
        sf, sa = 0.7, 0.3
        x, y = f / sf, a / sa
        gram = np.array([[x @ x, x @ y], [x @ y, y @ y]])
        lam_minus = np.linalg.eigvalsh(gram)[0]
        for c in (0.5, sf / sa, 5.0):
            q = reading_state(isotropic(f, a, sf, sa), c)["Q"]
            assert q == pytest.approx(lam_minus, rel=1e-9, abs=1e-12)


@pytest.mark.parametrize("force,acc", [
    ([3.0, 0.0], [0.0, 1.5]),       # perpendicular
    ([3.0, 0.0], [-1.5, 0.0]),      # anti-parallel
    ([3.0, 0.0], [0.0, 0.0]),       # force only
    ([0.0, 0.0], [1.5, 0.0]),       # acceleration only
])
def test_every_measured_input_gives_a_mass(force, acc):
    """T13 and the 29 Sept ruling: never undefined on measured input."""
    rng = _rng(61)
    for cov in (np.diag([0.04, 0.04, 0.01, 0.01]), _full_cov(rng, 2, 0.05)):
        r = ours29([Reading(force, acc, cov)], 2.0)
        for key in ("readout", "median"):
            assert np.isfinite(r[key]) and r[key] > 0
        assert 0 < r["interval"][0] < r["interval"][1] < np.inf


# ------------------------------------------------------------------ the exact branch and the type

@pytest.mark.parametrize("force,acc,kind,value", [
    ([4.0, 0.0], [2.0, 0.0], "one", 2.0),
    ([1.0, 2.0, 2.0], [0.5, 1.0, 1.0], "one", 2.0),
    ([0.0, 0.0], [0.0, 0.0], "every", None),
    ([4.0, 0.0], [0.0, 0.0], "none", None),
    ([0.0, 0.0], [2.0, 0.0], "none", None),
    ([4.0, 0.0], [0.0, 2.0], "none", None),
    ([4.0, 0.0], [-2.0, 0.0], "none", None),
])
def test_exact_branch(force, acc, kind, value):
    """T1 in pairing form: one mass iff codirectional with positive sizes; every mass at 0/0."""
    got = exact_mass(force, acc)
    assert got[0] == kind
    if value is not None:
        assert got[1] == pytest.approx(value)


def test_the_type_refuses_what_it_should():
    with pytest.raises(ValueError, match="positive definite"):
        Reading([1.0], [1.0], np.array([[1.0, 2.0], [2.0, 1.0]]))
    with pytest.raises(ValueError, match="orthonormal"):
        Reading([1.0, 0.0], [1.0, 0.0], np.eye(4), np.array([[1.0], [1.0]]))
    with pytest.raises(ValueError, match="2 x 2|4 x 4"):
        Reading([1.0, 0.0], [1.0, 0.0], np.eye(2))
    with pytest.raises(ValueError, match="centre"):
        ours29([isotropic([1.0], [1.0], 1.0, 1.0)], 0.0)
    assert instrument_ratio(Reading([1.0, 0.0], [1.0, 0.0], np.diag([1.0, 2.0, 1.0, 1.0]))) is None


@pytest.mark.parametrize("d", [1, 2, 3])
@pytest.mark.parametrize("n", [1, 5, 40])
def test_stacked_mean_length_closed_form(d, n):
    """A(m) for all-zero readings: the stacked candidates are N(0, v(m) I_{nd}), so
    E|X| = sqrt(2 v) Gamma((k+1)/2) / Gamma(k/2), k = nd. The integral's window and step
    (T_HALF, T_NODES, set 25 Oct) give this to 3e-14 (rounding over up to 120 components);
    the old +-40 window gave 6e-13."""
    from scipy.special import gammaln
    from met.full import _stacked_mean_length
    readings = [isotropic(np.zeros(d), np.zeros(d), 0.2, 0.1) for _ in range(n)]
    m = np.exp(np.linspace(-3.0, 3.0, 13))
    v = 1.0 / (m ** 2 / 0.2 ** 2 + 1.0 / 0.1 ** 2)
    k = n * d
    exact = np.sqrt(2 * v) * math.exp(gammaln((k + 1) / 2) - gammaln(k / 2))
    np.testing.assert_allclose(_stacked_mean_length(m, readings), exact, rtol=1e-13, atol=0)
