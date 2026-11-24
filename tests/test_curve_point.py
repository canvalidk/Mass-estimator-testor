"""The point read off the evidence curve: where the curve is least (T81; problems/O20.md).

Harvested 23 Nov 2026 from the exploration branch `explore/e14-ours34-walls` (Antigravity with Can,
22 Nov; its Theorems E14-1 and E14-2) and checked here on met.evidence, which shares no code with
that branch. Nothing in this file adopts a point: `definition/equation.md` and `met/ours33.py` are unchanged.
"""
import math

import numpy as np
import pytest
from scipy import stats

from met.evidence import Curve, _least_on, common_force, leftover, lump
from met.full import Reading, isotropic


def _any_reading(rng, m, d, k, r):
    """A true pair of size r along a random direction of a random k-dimensional declared subspace,
    read through a full correlated covariance."""
    q, _ = np.linalg.qr(rng.normal(size=(d, d)))
    n = q[:, :k]
    a = n @ rng.normal(size=k)
    a *= r / np.linalg.norm(a)
    mix = rng.normal(size=(2 * d, 2 * d))
    cov = 0.3 * (mix @ mix.T) / (2 * d) + np.diag([0.04] * d + [0.01] * d)
    y = np.concatenate([m * a, a]) + np.linalg.cholesky(cov) @ rng.normal(size=2 * d)
    return Reading(y[:d], y[d:], 0.5 * (cov + cov.T), subspace=n)


def _least_mass(curve):
    least, u = curve.least()
    return curve.mass_at(u), least


def _body(rng, m, n, r, sf=0.2, sa=0.1, d=3):
    out = []
    for _ in range(n):
        a = rng.normal(size=d)
        a *= r * sa / np.linalg.norm(a)
        out.append(isotropic(m * a + sf * rng.normal(size=d), a + sa * rng.normal(size=d), sf, sa))
    return out


# --- (1) the least composes, under any way of composing ------------------------------------------

def test_the_least_of_a_composite_curve_is_the_composite_of_the_parts_leasts():
    """Functional invariance of maximum likelihood (Zehna 1966) on T74's composition: full
    covariances, declared subspaces of every dimension, two bodies. Lumping, common force, and a
    composition with no physics in it (the product of the two masses) all bottom where the parts
    do, so nothing here is special to L5 and L14."""
    rng = np.random.default_rng(20261123)
    c1 = Curve([_any_reading(rng, 2.0, 3, k, 6.0) for k in (3, 2, 3, 1, 2)], scale=2.0)
    c2 = Curve([_any_reading(rng, 3.0, 3, k, 6.0) for k in (2, 3, 3, 2)], scale=2.0)
    (m1, k1), (m2, k2) = _least_mass(c1), _least_mass(c2)

    totals = np.linspace(m1 + m2 - 0.4, m1 + m2 + 0.4, 161)
    joined = lump(c1, c2, totals)
    assert totals[np.argmin(joined)] == pytest.approx(m1 + m2, abs=6e-3)
    assert joined.min() == pytest.approx(k1 + k2, rel=1e-6)

    reduced = 1.0 / (1.0 / m1 + 1.0 / m2)
    grid = np.linspace(reduced - 0.15, reduced + 0.15, 121)
    assert grid[np.argmin(common_force(c1, c2, grid))] == pytest.approx(reduced, abs=3e-3)

    def product_curve(p):
        return _least_on(lambda t: c1(np.exp(t)) + c2(p / np.exp(t)), math.log(0.2), math.log(20.0))[0]
    products = np.linspace(m1 * m2 - 0.8, m1 * m2 + 0.8, 81)
    assert products[np.argmin([product_curve(p) for p in products])] == pytest.approx(m1 * m2, abs=2.5e-2)


def test_level_set_midpoints_do_not_compose():
    """Two quadratic curves, centres 4 and 6, widths 1 and 2 (T74 (2): the lump has centre 10 and
    width sqrt 5). The geometric midpoint of the set K <= 1 fails lumping; the arithmetic midpoint
    passes it on these curves and fails the swap."""
    quad = lambda mu, s: (lambda m: (np.asarray(m, dtype=float) - mu) ** 2 / s ** 2)
    totals = np.linspace(5.0, 15.0, 2001)
    inside = totals[lump(quad(4.0, 1.0), quad(6.0, 2.0), totals) <= 1.0]
    assert inside[0] == pytest.approx(10 - math.sqrt(5), abs=6e-3) and inside[-1] == pytest.approx(10 + math.sqrt(5), abs=6e-3)
    geometric = lambda lo, hi: math.sqrt(lo * hi)
    parts = geometric(3.0, 5.0) + geometric(4.0, 8.0)
    assert parts == pytest.approx(9.530, abs=1e-3) and geometric(inside[0], inside[-1]) == pytest.approx(9.747, abs=6e-3)
    # the swap sends the set [3, 5] to [1/5, 1/3]; its arithmetic midpoint is 4/15, and 1/4 is asked
    assert 0.5 * (1 / 5 + 1 / 3) == pytest.approx(4 / 15) and abs(4 / 15 - 1 / 4) > 0.01


# --- (2) the uniqueness construction ---------------------------------------------------------------

def _doubling(curve, grid):
    """On an angle grid: n -> 2n for the scaled lumping, G(m) = inf over (x + y)/2 = m of g(x) + g(y),
    and for the scaled common force, the same on inverse masses."""
    mass = curve.scale * np.tan(grid)
    mass[-1] = math.inf

    def at(vals, x):
        x = np.asarray(x, dtype=float)
        u = np.where(np.isinf(x), 0.5 * math.pi, np.arctan(np.where(np.isinf(x), 0.0, x) / curve.scale))
        return np.interp(u, grid, vals)

    def inverse(t):
        return np.where(t > 0, 1.0 / np.where(t > 0, t, 1.0), math.inf)

    def mean(vals):
        out = np.empty_like(vals)
        for i, m in enumerate(mass):
            if math.isinf(m):
                out[i] = vals[-1] + vals.min()
                continue
            x = np.linspace(0.0, 2.0 * m, 1201)
            out[i] = np.min(at(vals, x) + at(vals, 2.0 * m - x))
        return out

    def harmonic(vals):
        out = np.empty_like(vals)
        for i, m in enumerate(mass):
            if m == 0:
                out[i] = vals[0] + vals.min()
                continue
            w = 0.0 if math.isinf(m) else 1.0 / m
            p = np.linspace(0.0, 2.0 * w, 1201)
            out[i] = np.min(at(vals, inverse(p)) + at(vals, inverse(2.0 * w - p)))
        return out

    return mass, mean, harmonic


def test_lumping_alone_leaves_the_heavy_side_bounded_and_the_two_together_pinch_both_sides():
    """E14-1's construction on the curve of a real series (three 3D readings). Composing n bodies
    with this curve at a common acceleration, rescaled, keeps the least where it was and raises the
    wall on the light side without bound; on the heavy side it can never exceed the curve's own
    value at infinite mass. The same under a common force, the sides exchanged. One after the other
    raises both walls, above the lower bound in E14-1's proof, and the least has not moved."""
    curve = Curve(_body(np.random.default_rng(7), 2.0, 3, 2.5))
    least, u_star = curve.least()
    m_star = curve.mass_at(u_star)
    grid = np.linspace(0.0, 0.5 * math.pi, 361)
    k = curve.at_angle(grid) - least
    mass, mean, harmonic = _doubling(curve, grid)
    eps = 0.5
    left, right = mass <= m_star - eps, mass >= m_star + eps
    depth = min(k[mass <= m_star - eps / 2].min(), k[mass >= m_star + eps / 2].min())   # c(eps / 2)

    lumped = [k]
    for _ in range(3):
        lumped.append(mean(lumped[-1]))
    for n_power, g in enumerate(lumped):
        assert g[right].min() <= k[-1] + 1e-9                    # never above K(inf) - K_min
    assert lumped[3][left].min() > 6 * lumped[0][left].min()     # n = 8: the light side keeps rising

    walls = []
    for n_power in (1, 2):
        q = lumped[n_power]
        for _ in range(n_power):
            q = harmonic(q)
        q = q - q.min()
        n = 2 ** n_power
        bound = eps * depth * n / (2.0 * (m_star + eps))
        walls.append(min(q[left].min(), q[right].min()))
        assert walls[-1] >= bound
        assert mass[np.argmin(q)] == pytest.approx(m_star, abs=0.02)
    assert walls[1] > 3 * walls[0] and walls[1] > k[-1] * 0.4    # both sides, and still growing


# --- (3) on the data themselves the least does not add ---------------------------------------------

def _same_data_defects(size, n_series, seed, sf=0.2, sa=0.1):
    """The 28 Oct scene (tests/test_point_composition.py): bodies of mass 1 and 2 share each
    acceleration reading and the lump's force reading is F1 + F2, its noise added."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_series):
        a = rng.standard_normal((4, 3))
        a *= size / np.linalg.norm(a, axis=1, keepdims=True)
        at = a + sa * rng.standard_normal((4, 3))
        f1 = 1.0 * a + sf * rng.standard_normal((4, 3))
        f2 = 2.0 * a + sf * rng.standard_normal((4, 3))
        p = [leftover([isotropic(f, x, s, sa) for f, x in zip(force, at)])["mass"]
             for force, s in ((f1, sf), (f2, sf), (f1 + f2, sf * math.sqrt(2)))]
        if all(math.isfinite(v) and v > 0 for v in p):
            out.append(abs(p[2] - p[0] - p[1]) / (p[0] + p[1]))
    return np.array(out)


def test_on_the_same_data_the_least_fails_lumping_like_every_other_point():
    """Measured 23 Nov, 30 series, seed 7: 0.122% (strong), 3.35% (mid), 27.9% (weak), where the
    median of ours29's law had 0.273%, 2.62%, 23.25% on 28 Oct. So the 28 Oct finding stands for the
    least as well: fitting the summed data is not the sum of the fits at finite noise."""
    strong, mid, weak = (float(np.median(_same_data_defects(s, 12, 7))) for s in (2.0, 0.4, 0.1))
    assert strong < 0.01 < mid < weak and weak > 0.08


# --- (4) where there is no single interior least ----------------------------------------------------

def test_on_one_instrument_the_least_is_interior_exactly_when_the_summed_dot_product_is_positive():
    rng = np.random.default_rng(11)
    for _ in range(150):
        n, d = int(rng.integers(1, 5)), int(rng.integers(1, 4))
        readings = [isotropic(rng.normal(size=d), rng.normal(size=d), 0.3, 0.7) for _ in range(n)]
        dot = sum(float(r.force @ r.acceleration) for r in readings)
        assert (dot > 0) == (not leftover(readings)["at_end"])


def test_at_no_signal_half_of_all_series_have_their_least_at_zero_or_infinite_mass():
    """One scalar reading whose true pair is r noise units in each channel: the least is at an end
    with probability 2 Phi(r) Phi(-r): 0.50, 0.27, 0.044 and 0.003 at r = 0, 1, 2, 3. Measured on
    the curve itself at r = 0 and r = 1 (400 readings each)."""
    rng = np.random.default_rng(5)
    for r in (0.0, 1.0):
        ends = sum(leftover([isotropic([2.0 * 0.1 * r + 0.2 * rng.normal()], [0.1 * r + 0.1 * rng.normal()], 0.2, 0.1)])["at_end"]
                   for _ in range(400))
        expected = 2 * stats.norm.cdf(r) * stats.norm.cdf(-r)
        assert ends / 400 == pytest.approx(expected, abs=0.07)
    assert 2 * stats.norm.cdf(3.0) * stats.norm.cdf(-3.0) < 0.003


@pytest.mark.parametrize("n,r,recorded,slack", [(5, 1.0, 0.149, 0.012), (5, 2.0, 0.001, 0.001), (20, 1.0, 0.021, 0.005)])
def test_how_often_a_series_of_3d_readings_has_its_least_at_an_end(n, r, recorded, slack):
    """By the rule of the test above (the summed dot product is not positive), 20 000 series each;
    r is the true pair's size in noise units in each channel."""
    rng = np.random.default_rng(1000 * n + int(10 * r))
    a = rng.normal(size=(20000, n, 3))
    a *= r * 0.1 / np.linalg.norm(a, axis=2, keepdims=True)
    force = 2.0 * a + 0.2 * rng.normal(size=a.shape)
    seen = a + 0.1 * rng.normal(size=a.shape)
    share = float(np.mean(np.einsum("sni,sni->s", force, seen) <= 0))
    assert share == pytest.approx(recorded, abs=slack)


def test_the_kinds_of_curve_with_no_single_interior_least():
    # flat with a leftover: a force-only and an acceleration-only reading of equal size
    flat = leftover([isotropic([1.0], [0.0], 1.0, 1.0), isotropic([0.0], [1.0], 1.0, 1.0)])
    assert flat["flat"] and flat["least"] == pytest.approx(1.0) and flat["mass"] is None
    # least at one end only: the body did not move, or was not pushed
    assert leftover([isotropic([10.0], [0.0], 1.0, 1.0)])["mass"] == math.inf
    assert leftover([isotropic([0.0], [10.0], 1.0, 1.0)])["mass"] == 0.0
    # force against acceleration: the curve is highest inside and equally low at both ends
    against = Curve([isotropic([1.0], [-1.0], 1.0, 1.0)])
    assert against(np.array([0.0, 1.0, math.inf])) == pytest.approx([1.0, 2.0, 1.0])
    # two equal interior lows need two noise ratios (on one the curve has a single turn, T73 (4))
    two = Curve([isotropic([0.1], [1.0], 0.1, 1.0), isotropic([1.0], [0.1], 1.0, 0.1)], scale=1.0)
    low_a, low_b, between, end_a, end_b = two(np.array([0.1, 10.0, 1.0, 0.0, math.inf]))
    assert low_a == pytest.approx(low_b, rel=1e-9) and low_a == pytest.approx(0.980, abs=1e-3)
    assert between == pytest.approx(1.604, abs=1e-3) and end_a == pytest.approx(2.0) and end_b == pytest.approx(2.0)


def test_the_least_jumps_from_one_end_to_the_other_as_a_reading_against_itself_tips():
    # force against acceleration, the two almost equal in noise units: one part in a million decides
    heavy = leftover([isotropic([1.000001], [-1.0], 1.0, 1.0)])
    light = leftover([isotropic([1.0], [-1.000001], 1.0, 1.0)])
    assert heavy["mass"] == math.inf and light["mass"] == 0.0


def test_no_positive_value_at_a_one_sided_curve_keeps_composition_at_the_point():
    """E14-2. A body at 3 (F = 3, a = 1, unit noise) under a common force with a body that did not
    move (F = 10, a = 0): the composite's curve has one interior least, at 3. If the composite's
    point is to be the reduced mass of the parts' points, 1/3 = 1/3 + 1/R, the second body's point
    is infinite. The same with lumping and a body that was not pushed gives zero."""
    three = Curve([isotropic([3.0], [1.0], 1.0, 1.0)], scale=1.0)
    still = Curve([isotropic([10.0], [0.0], 1.0, 1.0)], scale=1.0)
    unpushed = Curve([isotropic([0.0], [10.0], 1.0, 1.0)], scale=1.0)
    grid = np.linspace(0.5, 12.0, 461)
    under_force = common_force(three, still, grid)
    assert grid[np.argmin(under_force)] == pytest.approx(3.0, abs=0.03) and under_force.min() == pytest.approx(0.0, abs=1e-6)
    assert np.sum(under_force <= under_force.min() + 1e-4) <= 16          # one low, not a plateau
    lumped = lump(three, unpushed, grid)
    assert grid[np.argmin(lumped)] == pytest.approx(3.0, abs=0.03) and lumped.min() == pytest.approx(0.0, abs=1e-6)


def test_the_two_endpoints_of_padding_with_exact_zeros():
    """The 22 Nov section of problems/O20.md (Antigravity with Can). The old cloud readout's limit
    under zero readings is s E[sin t] / E[cos t], t = atan(m / s): on equal-weight two-point laws it
    fails lumping (1.618 + 1.618 against 3.075) and reads its anchor (2.000 at c = 2, 1.766 at c = 1).
    The unknown-scale law's limit is where the curve is least: the golden ratio on (2, 1), (1, 1)."""
    def circular(masses, s):
        t = np.arctan(np.array(masses) / s)
        return s * np.sin(t).mean() / np.cos(t).mean()
    assert 2 * circular([1.0, 3.0], 1.0) == pytest.approx(3.236, abs=1e-3)
    assert circular([2.0, 6.0], 1.0) == pytest.approx(3.075, abs=1e-3)
    assert circular([1.0, 4.0], 2.0) == pytest.approx(2.000, abs=1e-3)
    assert circular([1.0, 4.0], 1.0) == pytest.approx(1.766, abs=1e-3)
    golden = leftover([isotropic([2.0], [1.0], 1.0, 1.0), isotropic([1.0], [1.0], 1.0, 1.0)])["mass"]
    assert golden == pytest.approx(0.5 * (1 + math.sqrt(5)), rel=1e-7)


# --- (5) against the other rows ---------------------------------------------------------------------

def test_the_least_ignores_exact_zeros_and_follows_units_the_swap_and_the_coordinate():
    rng = np.random.default_rng(9)
    series = _body(rng, 2.0, 4, 5.0)
    base = leftover(series)["mass"]
    zeros = [isotropic(np.zeros(3), np.zeros(3), 5.0, 0.01), isotropic(np.zeros(1), np.zeros(1), 0.01, 7.0)]
    assert leftover(series + zeros, scale=Curve(series).scale)["mass"] == base        # another instrument's zeros
    lam, mu = 7.0, 0.3                                                                # units: m -> m lam / mu
    scaled = [isotropic(lam * r.force, mu * r.acceleration, lam * 0.2, mu * 0.1) for r in series]
    assert leftover(scaled)["mass"] == pytest.approx(base * lam / mu, rel=1e-6)
    swapped = [isotropic(r.acceleration, r.force, 0.1, 0.2) for r in series]          # the swap: m -> 1 / m
    assert leftover(swapped)["mass"] == pytest.approx(1.0 / base, rel=1e-6)
    # the coordinate the curve is read in does not matter: a different declared scale re-maps the arc
    assert leftover(series, scale=37.0)["mass"] == pytest.approx(base, rel=1e-6)
