"""ours30 (R-2026-10-29-1): the definition since 29 Oct. The misfit weight reads each reading's misfit
at its best positive mass, g_i = min(1, (k_i - 1/2) / Q_plus_i); the rest is ours29's.

Here: met/full's ours30 is ours29 with positivity; the batch pipeline's `positivity: true` agrees with
it (isotropic noise end to end, per-axis noise for the weights); 1D readings; the preset validator.
The properties of the weight (T30's bounded tilt, no gate, Z, swap, units, the supremum) are in
tests/test_sign_numerator.py."""
import numpy as np
import pytest

from met.analyst import Estimator, estimate, misfit_weights
from met.full import Reading, isotropic, ours29, ours30, reading_state
from met.law import ReadingSet

LAW = {"reference": "cartesian2", "nuisance": "radius"}
OURS30 = {"rule": "likelihood_product", "prior": [{"uniform_angle": {"center": "first_reading"}}],
          "misfit": "kink", "positivity": True, "calibrate": "bartlett"}


def _series(seed, d, n=7):
    """Readings with some force/acceleration pairs pointing apart, so Q_plus differs from Q."""
    rng = np.random.default_rng(seed)
    force = rng.normal(size=(n, d)) * 2 + 0.6
    acc = rng.normal(size=(n, d)) + 0.3
    return force, acc


def test_ours30_is_ours29_with_positivity():
    force, acc = _series(300, 3)
    readings = [isotropic(f, a, 0.5, 0.25) for f, a in zip(force, acc)]
    a, b = ours30(readings, 2.0), ours29(readings, 2.0, positivity=True)
    assert a["point"] == b["point"] and a["interval"] == b["interval"] and a["g"] == b["g"]
    assert a["point"] == a["median"]


@pytest.mark.parametrize("d,seed", [(1, 301), (2, 302), (3, 303)])
def test_isotropic_noise_agrees_with_the_batch_pipeline(d, seed):
    force, acc = _series(seed, d)
    readings = [isotropic(f, a, 0.5, 0.25) for f, a in zip(force, acc)]
    mine = ours30(readings, 2.0)
    assert any(np.dot(f, a) < 0 for f, a in zip(force, acc))          # the sign is exercised
    batch = estimate(force, acc, 0.5, 0.25, ["ratio_of_means_rss", "median", "interval_95"],
                     law=LAW, combine=OURS30)
    assert mine["median"] == pytest.approx(float(batch["median"]), rel=1e-6)
    assert np.allclose(mine["interval"], batch["interval_95"], rtol=1e-6)
    assert mine["readout"] == pytest.approx(float(batch["ratio_of_means_rss"]), rel=1e-9)
    rs = ReadingSet(force[None], acc[None], np.full((1, len(force)), 0.5), np.full((1, len(force)), 0.25))
    assert np.allclose(misfit_weights(rs, OURS30)[0], mine["g"], rtol=1e-9, atol=0)


@pytest.mark.parametrize("d,seed", [(2, 311), (3, 312)])
def test_per_axis_noise_weights_agree_with_met_full(d, seed):
    """The batch pipeline's positive-slope search against met/full's Q_plus, per-axis SDs."""
    force, acc = _series(seed, d)
    rng = np.random.default_rng(seed + 1)
    sf, sa = rng.uniform(0.3, 1.0, size=d), rng.uniform(0.1, 0.5, size=d)
    rs = ReadingSet(force[None], acc[None], np.broadcast_to(sf, force.shape)[None],
                    np.broadcast_to(sa, force.shape)[None])
    batch = misfit_weights(rs, OURS30)[0]
    full = []
    for f, a in zip(force, acc):
        st = reading_state(Reading(f, a, np.diag(np.r_[sf ** 2, sa ** 2])), 1.0)
        full.append(min(1.0, st["kappa"] / st["Q_plus"]) if st["Q_plus"] > 0 else 1.0)
    assert np.allclose(batch, full, rtol=1e-7)
    assert np.any(np.array(full) < 1)


def test_one_dimensional_readings_are_weighted_only_when_the_sign_contradicts_m_positive():
    force = np.array([[2.0], [3.0], [-1.5]])
    acc = np.array([[1.0], [1.4], [1.2]])
    rs = ReadingSet(force[None], acc[None], np.full((1, 3), 0.1), np.full((1, 3), 0.1))
    g = misfit_weights(rs, OURS30)[0]
    assert g[0] == g[1] == 1.0
    assert g[2] == pytest.approx(0.5 / min(150.0 ** 2 / 1e2, 120.0 ** 2 / 1e2), rel=1e-12)
    old = misfit_weights(rs, {**OURS30, "positivity": False})[0]
    assert np.all(old == 1.0)                                          # ours29: inert in 1D


def test_the_preset_validator_takes_positivity_only_with_misfit():
    base = {"name": "x", "reduce": "none", "law": LAW, "readouts": ["median"]}
    Estimator(dict(base, combine=OURS30))
    with pytest.raises(ValueError, match="positivity"):
        Estimator(dict(base, combine={k: v for k, v in OURS30.items() if k != "misfit"}))
    with pytest.raises(ValueError, match="positivity"):
        Estimator(dict(base, combine=dict(OURS30, positivity="yes")))


# ------------------------------------------------------------------ the robot table under ours30

@pytest.fixture(scope="module")
def robot_baseline():
    from met import robot
    cal = robot.calibrate()
    rd, _ = robot.windows("1-baseline", 0.2, cal)
    return cal, ours30(rd, robot.series_ratio(rd))


def test_robot_baseline_under_ours30(robot_baseline):
    """tools/robot_run.py's default row, 29 Oct: the median as ours29's to 1e-4 (the robot's pairs
    are acute, so only the discount changes), the interval slightly narrower."""
    _, r = robot_baseline
    assert r["point"] == pytest.approx(0.9213, abs=1e-4)
    assert r["interval"] == pytest.approx((0.9145, 0.9282), abs=1e-4)
    assert min(r["g"]) == pytest.approx(0.18, abs=0.005)


@pytest.mark.xfail(strict=True, reason=(
    "Experiment vs ours30 (29 Oct): 1-baseline, 0.2 s windows, floor budget gives [0.9145, 0.9282]; "
    "the static calibration mass is 0.9323. As under ours29: the robot's declared budget is not "
    "yet right (problems/O16.md, 28 Oct; brief B01)."))
def test_robot_baseline_interval_covers_the_calibration_mass_under_ours30(robot_baseline):
    cal, r = robot_baseline
    lo, hi = r["interval"]
    assert lo <= cal.mass <= hi


# ------------------------------------------------------------------ the equation note d3, section 3

def test_the_closed_forms_in_equation_note_d3():
    """outbox/VD-docs/the_mass_estimator_equation_d3.md §3: with isotropic noise and c = s,
    H(t) = (t^2 A + 2tC + B)/(1 + t^2), Q = (A + B - Delta)/2, Q_plus = Q if C > 0 else min(A, B),
    and eta from Z = (B - A) + 2iC, written independently of met/full and checked against it."""
    import math
    import met.full as full
    rng = np.random.default_rng(2915)
    for _ in range(60):
        d = int(rng.integers(1, 4))
        sf, sa = rng.uniform(0.2, 2), rng.uniform(0.2, 2)
        f, a = rng.normal(size=d) * rng.uniform(0.1, 4), rng.normal(size=d) * rng.uniform(0.1, 4)
        x, y = f / sf, a / sa
        A, B, C = x @ x, y @ y, x @ y
        r = isotropic(f, a, sf, sa)
        for t in (0.0, 0.3, 1.0, 7.0):
            h = full._h_terms(r.pre(), sf / sa, np.array([math.atan(t)]))[0]
            assert h == pytest.approx((t * t * A + 2 * t * C + B) / (1 + t * t), rel=1e-11, abs=1e-12)
        st = reading_state(r, sf / sa)
        q = 0.5 * (A + B - math.sqrt((A - B) ** 2 + 4 * C * C))
        assert st["Q"] == pytest.approx(q, rel=1e-9, abs=1e-10)
        assert st["Q_plus"] == pytest.approx(q if C > 0 else min(A, B), rel=1e-9, abs=1e-10)
    for _ in range(3):                      # series whose readings disagree, so eta < 1
        acc = rng.normal(size=(8, 3))
        m = rng.uniform(0.5, 4, size=8)
        force = m[:, None] * acc + 0.05 * rng.normal(size=(8, 3))
        acc = acc + 0.025 * rng.normal(size=(8, 3))
        out = ours30([isotropic(fi, ai, 0.05, 0.025) for fi, ai in zip(force, acc)], 2.0, readout=False)
        g = np.array(out["g"])
        x, y = force / 0.05, acc / 0.025
        z = ((y * y).sum(1) - (x * x).sum(1)) + 2j * (x * y).sum(1)
        tot = (g * z).sum()
        score = -0.5 * g * np.imag(np.conj(z) * np.exp(1j * np.angle(tot)))
        assert out["eta"] < 0.5
        assert out["eta"] == pytest.approx(min(1.0, abs(tot) / (score ** 2).sum()), rel=1e-9)
