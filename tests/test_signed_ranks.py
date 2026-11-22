"""How much of a law T45's signed ranks supply without a starting law (T79; exploration E12,
`e12-conditional-ranks.md`, Codex, 21 Nov; checked by Claude). Scalar readings of one instrument
in noise units; G is the series' 2x2 matrix of dot products; Z(theta) = u'Gv / sqrt(u'Gu) is T45's
signed score and Phi(-Z) its rank, exactly uniform at the true mass."""
import math

import numpy as np
import pytest
from scipy import stats

from met.conditional_score import gaussian_score_test
from met.full import isotropic


def _rank(gram, theta):
    u = np.stack([np.sin(theta), np.cos(theta)])
    v = np.stack([np.cos(theta), -np.sin(theta)])
    along = np.einsum("ia,ij,ja->a", u, gram, u)
    return stats.norm.cdf(-np.einsum("ia,ij,ja->a", u, gram, v) / np.sqrt(along)), along


def _rise(values):
    return float(np.sum(np.clip(np.diff(values), 0.0, None)))


def _share(rows):
    """W over the whole line of slopes, and the rising share M on the positive masses."""
    rows = np.asarray(rows, dtype=float)
    gram = rows.T @ rows
    low, high = np.linalg.eigvalsh(gram)
    if high <= 0:
        return 0.0, 0.0
    whole = 2 * stats.norm.cdf(math.sqrt(high) - math.sqrt(max(low, 0.0))) - 1
    return whole, _rise(_rank(gram, np.linspace(0, 0.5 * math.pi, 200001)[1:-1])[0])


def test_the_ranks_rise_by_w_and_fall_by_w_over_a_full_turn():
    rng = np.random.default_rng(79)
    theta = np.linspace(0.0, math.pi, 400001)
    for _ in range(8):
        rows = rng.normal(size=(rng.integers(2, 7), 2)) * rng.uniform(0.3, 3) + rng.normal(size=2) * rng.uniform(0, 3)
        gram = rows.T @ rows
        low, high = np.linalg.eigvalsh(gram)
        rank, along = _rank(gram, theta)
        expected = 2 * stats.norm.cdf(math.sqrt(high) - math.sqrt(low)) - 1
        assert _rise(rank) == pytest.approx(expected, abs=1e-9)
        assert _rise(-rank) == pytest.approx(expected, abs=1e-9)
        # it rises exactly where the candidate line holds more than the geometric mean of the eigenvalues
        step = np.diff(rank)
        moving = np.abs(step) > 1e-13
        mid = 0.5 * (along[1:] + along[:-1])
        assert np.mean((step > 0)[moving] == (mid > math.sqrt(high * low))[moving]) > 0.9999


def test_the_rank_is_t45s_and_it_is_not_a_cdf_on_the_positive_masses():
    # forces (3, 0), accelerations (0, 1) in noise units: G = diag(9, 1). The rank is 1/2 at both
    # ends of the positive masses and below 1/2 between them: no CDF has these ordinates.
    gram = np.diag([9.0, 1.0])
    rank, _ = _rank(gram, np.linspace(0.0, 0.5 * math.pi, 9))
    assert rank[0] == pytest.approx(0.5) and rank[-1] == pytest.approx(0.5) and np.all(rank[1:-1] < 0.5)
    series = [isotropic([3.0], [0.0], 1.0, 1.0), isotropic([0.0], [1.0], 1.0, 1.0)]
    for mass in (0.3, 1.0, 4.0):
        rank_here = float(_rank(gram, np.array([math.atan(mass)]))[0][0])
        assert abs(gaussian_score_test(series, mass)["normal_score"]) == pytest.approx(abs(stats.norm.ppf(rank_here)), rel=1e-9)


@pytest.mark.parametrize("name,rows,whole,positive", [
    ("the all-zero experiment", [[0.0, 0.0]], 0.0, 0.0),
    ("right angles, equal sizes (scene C)", [[20.0, 0.0], [0.0, 20.0]], 0.0, 0.0),
    ("one reading at half a noise unit a channel", [[0.5, 0.5]], 0.5205, 0.3829),
    ("one reading at (3, 2)", [[3.0, 2.0]], 0.9997, 0.9759),
    ("the near pair", [[55.0, 500.0], [56.0, 500.0]], 1.0, 1.0),
    ("the far pair", [[50.0, 500.0], [60.0, 500.0]], 1.0, 1.0),
])
def test_how_much_of_a_law_the_readings_alone_supply(name, rows, whole, positive):
    got = _share(rows)
    assert got[0] == pytest.approx(whole, abs=1e-4) and got[1] == pytest.approx(positive, abs=1e-4), name
    assert got[1] <= got[0] + 1e-9
