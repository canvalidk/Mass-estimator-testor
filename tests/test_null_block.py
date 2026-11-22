"""A block of noise-only readings beside a fixed informative reading (T77; exploration E2, Codex,
17 Nov). E2's own checks (tests/test_null_combination.py, tests/test_null_combination_score.py) are
exact identities with no simulation; this file is Claude's simulation of the limit they lead to,
on the scalar one-instrument slice in noise coordinates, where T45's statistic is
(u'Gv)^2 / (u'Gu) with G the series' 2x2 matrix of dot products."""
import math

import numpy as np
import pytest

ANGLES = np.linspace(0.0, 0.5 * math.pi, 1441)                    # the positive masses
U = np.stack([np.sin(ANGLES), np.cos(ANGLES)])                    # along each candidate line
V = np.stack([np.cos(ANGLES), -np.sin(ANGLES)])                   # across it


def _null_gram(rng, n, reps):
    """The matrix of dot products of n noise-only scalar readings, by Bartlett's decomposition."""
    a = np.sqrt(rng.chisquare(n, reps))
    b = rng.normal(size=reps)
    c = np.sqrt(rng.chisquare(n - 1, reps))
    return np.stack([np.stack([a * a, a * b], axis=-1), np.stack([a * b, b * b + c * c], axis=-1)], axis=-2)


def _largest_statistic(gram):
    score = np.einsum("ia,rij,ja->ra", U, gram, V)
    variance = np.einsum("ia,rij,ja->ra", U, gram, U)
    return np.max(score ** 2 / variance, axis=1)


def _kept_width(gram):
    """The median, over series, of the total angle of positive masses the 95% score set keeps."""
    score = np.einsum("ia,rij,ja->ra", U, gram, V)
    variance = np.einsum("ia,rij,ja->ra", U, gram, U)
    return float(np.median(np.sum(score ** 2 / variance <= 3.8415, axis=1)) * (ANGLES[1] - ANGLES[0]))


def _informative(rng, reps, signal=20.0, angle=math.atan(2.0)):
    y = signal * np.array([math.sin(angle), math.cos(angle)]) + rng.normal(size=(reps, 2))
    return y[:, :, None] * y[:, None, :]


def test_enough_noise_only_readings_erase_a_fixed_informative_reading_from_the_score_set():
    # one reading at signal 20 (true mass 2 in noise units), alone and beside n noise-only readings.
    # Alone, never is every positive mass kept. Beside 1e8 of them the 95% set is every positive
    # mass in 0.8535 of series, 1 - exp(-3.8415 / 2): the block's own limit, with the reading gone.
    rng = np.random.default_rng(77)
    reps = 4000
    block = _informative(rng, reps)
    assert np.mean(_largest_statistic(block) <= 3.8415) == 0.0
    beside = block + _null_gram(rng, 10 ** 8, reps)
    assert np.mean(_largest_statistic(beside) <= 3.8415) == pytest.approx(1 - math.exp(-3.8415 / 2), abs=0.02)
    # on the way there (measured 22 Nov): at 400, T76 (4)'s count for the constant weight, the kept
    # set is sqrt(2) times as wide and never everything; at 4000 it is everything in 0.10 of series,
    # at 40000 in 0.72, at 1.6e6 in 0.854
    at_400 = block + _null_gram(rng, 400, reps)
    assert np.mean(_largest_statistic(at_400) <= 3.8415) == 0.0
    assert _kept_width(at_400) / _kept_width(block) == pytest.approx(math.sqrt(2.0), abs=0.08)
    assert np.mean(_largest_statistic(block + _null_gram(rng, 40000, reps)) <= 3.8415) > 0.6


def test_the_true_mass_stays_covered_whatever_the_block_does():
    # coverage is not what is lost: the test of the true mass is exact with or without the block
    rng = np.random.default_rng(78)
    reps = 20000
    u = np.array([math.sin(math.atan(2.0)), math.cos(math.atan(2.0))])
    v = np.array([u[1], -u[0]])
    for n in (2, 400, 10 ** 6):
        gram = _informative(rng, reps) + _null_gram(rng, n, reps)
        statistic = np.einsum("i,rij,j->r", u, gram, v) ** 2 / np.einsum("i,rij,j->r", u, gram, u)
        assert np.mean(statistic <= 3.8415) == pytest.approx(0.95, abs=0.006)
