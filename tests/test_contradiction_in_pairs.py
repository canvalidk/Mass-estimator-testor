"""What no single mass explains, read in pairs with no mass (problems/O14.md, 18 Nov, night;
tools/contradiction_in_pairs.py). Three identities that hold for every configuration of signals with no true mass, no
fitted mass and no count; what Z selects in them; and what "expected over observed" gives on the scenes."""
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import contradiction_in_pairs as cp                             # noqa: E402
import discount_lab as lab                                      # noqa: E402

SCENES = lab.scenes()


def _ws(name, sigma):
    readings, _, _, ratio = SCENES[name]
    return lab.whiten(readings, sigma, sigma * ratio)


@pytest.mark.parametrize("k", [1, 3])
@pytest.mark.parametrize("pushes", [[0.0, 0.0], [5.0, 5.0], [0.0, 2.0, 7.0], [0.5] * 5])
def test_the_three_identities_hold_at_every_signal_size_with_no_mass(k, pushes):
    r = cp.identity_cell(k, pushes, 200000)
    assert r["between"] == pytest.approx(1.0, abs=max(5 * r["between_se"], 0.01))
    if k > 1:
        assert r["inside"] == pytest.approx(1.0, abs=0.02)
    assert r["signal"] == pytest.approx(r["signal_exact"], abs=max(5 * r["signal_se"], 0.02))


def test_the_series_determinant_is_the_sum_over_pairs_and_insides():
    rng = np.random.default_rng(3)
    ws = [(rng.normal(size=k) + 1.0, rng.normal(size=k) + 0.5) for k in (1, 3, 2, 1)]
    a, b, c, _ = cp.parts(ws)
    pairs = sum(a[i] * b[j] + b[i] * a[j] - 2 * c[i] * c[j] for i in range(4) for j in range(i + 1, 4))
    G = sum(lab.gram(f, al) for f, al in ws)
    assert np.linalg.det(G) == pytest.approx(pairs + (a * b - c * c).sum(), rel=1e-10)


def test_no_mass_is_in_them_and_no_zero_reading():
    """Turning every reading by one angle in the (alpha, f) plane moves every line and leaves all four sums alone; so
    does adding exactly zero readings of any size."""
    rng = np.random.default_rng(4)
    ws = [(rng.normal(size=k) * 3, rng.normal(size=k) * 3 + 2) for k in (1, 3, 2, 1, 1)]
    base = cp.pair_sums(ws)
    th = 0.83
    turned = [(math.cos(th) * f + math.sin(th) * al, -math.sin(th) * f + math.cos(th) * al) for f, al in ws]
    zeros = ws + [(np.zeros(1), np.zeros(1)), (np.zeros(3), np.zeros(3))]
    for other in (cp.pair_sums(turned), cp.pair_sums(zeros)):
        for key, value in base.items():
            assert other[key] == pytest.approx(value, rel=1e-9), key


def test_z_selects_the_two_constants():
    """Weights C - c each have an exact identity; its expected side vanishes at a zero reading only for c = 2 between
    readings and c = 4 inside one."""
    C = np.array([0.5, 3.0, 40.0])
    assert np.allclose(cp.expected_between(0.0, C, 2.0), 0) and np.allclose(cp.expected_inside(0.0, 3, 4.0), 0)
    assert not np.allclose(cp.expected_between(0.0, C, 3.0), 0) and cp.expected_inside(0.0, 3, 2.0) != 0
    assert np.allclose(cp.expected_between(C, C[::-1], 2.0), C * C[::-1] * (C + C[::-1] - 6))
    assert np.allclose(cp.expected_inside(C, 3, 4.0), 2 * C * (C - 5))
    # the identity with another constant is as exact, and is not zero-safe
    rng = np.random.default_rng(5)
    f, al = cp.draw(rng, 1, [1.5, 0.0], 400000)
    Ci, Cj = (f[:, 0, 0] ** 2 + al[:, 0, 0] ** 2), (f[:, 1, 0] ** 2 + al[:, 1, 0] ** 2)
    X = (f[:, 0, 0] * al[:, 1, 0] - al[:, 0, 0] * f[:, 1, 0]) ** 2
    obs = (Ci - 3) * (Cj - 3) * X
    exact = (1.5 ** 4 + 3 * 1.5 ** 2) * 1.0 + 2 * (1.5 ** 2 + 1) * 1.0          # a_i b_j + b_i a_j + 2 b_i b_j at c = 3
    assert obs.mean() == pytest.approx(exact, abs=5 * obs.std() / math.sqrt(obs.size))
    assert cp.expected_between(Ci, Cj, 3.0).mean() == pytest.approx(exact, abs=0.2)
    f, al = cp.draw(rng, 3, [2.0], 400000)
    a, b, c = (f * f).sum(2)[:, 0], (al * al).sum(2)[:, 0], (f * al).sum(2)[:, 0]
    assert ((a + b - 2) * (a * b - c * c)).mean() == pytest.approx(cp.expected_inside(a + b, 3, 2.0).mean(), rel=0.02)


def test_where_the_fit_is_harmless_the_pairs_give_the_fitted_weight():
    for name in ("B: 5, 6", "5, 5.05, 6", "B and 40 noise-only readings", "B and 40 exactly zero readings",
                 "99 ordinary and one at 9"):
        ws = _ws(name, 0.1)
        assert cp.pair_weight(ws) == pytest.approx(float(lab.weights(ws, cp.FITTED)[0][0]), rel=2e-3), name
    assert cp.pair_weight(_ws("B: 5, 6", 0.1)) == pytest.approx(0.02024, rel=1e-3)
    for name in ("A: 5.5, 5.6 (each +- 0.1)", "99 ordinary near 5"):
        assert cp.pair_weight(_ws(name, 0.1)) == 1.0
    assert cp.pair_weight(_ws("aligned pair (2, 1) and (4, 2)", 0.01)) == 1.0


def test_right_angles_are_seen_and_regrouping_a_reading_into_its_axes_keeps_the_weight():
    """Scene C as two scalar readings and as one 2D reading: 1.25 sigma^2 either way (the fitted body: 1.18 and 1.00)."""
    two = cp.pair_weight(_ws("C: force 2, no acceleration; acceleration 1, no force", 0.01)) / 1e-4
    one = cp.pair_weight(_ws("one 2D reading F = (2, 0), a = (0, 1)", 0.01)) / 1e-4
    assert two == pytest.approx(1.25, abs=0.005) and one == pytest.approx(1.25, abs=0.005)
    fitted_two = float(lab.weights(_ws("C: force 2, no acceleration; acceleration 1, no force", 0.01), cp.FITTED)[0][0]) / 1e-4
    fitted_one = float(lab.weights(_ws("one 2D reading F = (2, 0), a = (0, 1)", 0.01), cp.FITTED)[0][0]) / 1e-4
    assert fitted_two == pytest.approx(1.18, abs=0.01) and fitted_one == pytest.approx(1.00, abs=0.01)
    # turning the scene (C'' is C' turned by 60 degrees) changes nothing
    a, b = (cp.pair_weight(_ws(n, 0.05)) for n in ("C': two readings (2, 0) and one (0, 1)",
                                                    "C'': two readings on m = 1.73 and one at right angles"))
    assert a == pytest.approx(b, rel=1e-9)


def test_the_sign_condition_is_not_in_the_pairs():
    """An anti-parallel reading lies on a line of negative mass: no pair is contradicted, and the weight is 1. The
    fitted body, read over positive masses only, discounts it."""
    ws = _ws("anti-parallel scalar reading F = 1, a = -1", 0.05)
    assert cp.pair_weight(ws) == 1.0
    assert float(lab.weights(ws, cp.FITTED)[0][0]) < 0.01
