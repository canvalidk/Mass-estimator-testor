"""Finite couplings, equality cases, and conditioning in mass-ray coordinates."""
from collections import Counter
from fractions import Fraction
from itertools import permutations

import numpy as np
import pytest


def reduced(x, y):
    return x * y / (x + y)


def laplace(sample, t):
    return np.mean(np.exp(-t * np.asarray(sample, dtype=float)))


def test_common_shape_one_laplace_value_detects_every_off_ray_coupling_on_fixture():
    s = (Fraction(1, 2), Fraction(2), Fraction(2), Fraction(11))
    a, b = Fraction(3), Fraction(5)
    inv_s = [1 / u for u in s]
    for other in set(permutations(s)):
        pairs = [(a * u, b * v) for u, v in zip(s, other)]
        totals = [(x + y) / (a + b) for x, y in pairs]
        inverse_reduced = [reduced(a, b) / reduced(x, y) for x, y in pairs]
        on_ray = s == other
        assert (Counter(totals) == Counter(s)) == on_ray
        assert (Counter(inverse_reduced) == Counter(inv_s)) == on_ray
        for normalized, reference in [(totals, s), (inverse_reduced, inv_s)]:
            deficit = laplace(reference, .61) - laplace(normalized, .61)
            if on_ray:
                assert abs(deficit) < 1e-14
            else:
                assert deficit > 1e-6


def test_all_quantiles_compose_exactly_only_for_the_comonotone_pairing_on_fixture():
    # These are different marginal shapes. Repeated atoms allow tied ranks.
    xs = (Fraction(1), Fraction(1), Fraction(4), Fraction(7))
    ys = (Fraction(2), Fraction(3), Fraction(3), Fraction(8))
    assert len({y / x for x, y in zip(xs, ys)}) > 1
    target_sum = [x + y for x, y in zip(xs, ys)]
    target_h = [reduced(x, y) for x, y in zip(xs, ys)]
    for other in set(permutations(ys)):
        comonotone = all((x - u) * (y - v) >= 0
                         for x, y in zip(xs, other) for u, v in zip(xs, other))
        sums = sorted(x + y for x, y in zip(xs, other))
        hs = sorted(reduced(x, y) for x, y in zip(xs, other))
        assert (sums == target_sum) == comonotone
        assert (hs == target_h) == comonotone
        for actual, target in [(sums, target_sum),
                               ([1 / h for h in hs], [1 / h for h in target_h])]:
            deficit = laplace(target, .51) - laplace(actual, .51)
            if comonotone:
                assert abs(deficit) < 1e-14
            else:
                assert deficit > 1e-6


def test_mixture_of_three_ray_points_breaks_both_median_compositions():
    pairs = [(Fraction(1), Fraction(1)), (Fraction(1), Fraction(4)),
             (Fraction(4), Fraction(1))]
    xs, ys = zip(*pairs)
    median_x, median_y = sorted(xs)[1], sorted(ys)[1]
    totals = [x + y for x, y in pairs]
    hs = [reduced(x, y) for x, y in pairs]
    assert (median_x + median_y, sorted(totals)[1]) == (2, 5)
    assert (reduced(median_x, median_y), sorted(hs)[1]) == (Fraction(1, 2), Fraction(4, 5))
    for (x, y), total, h in zip(pairs, totals, hs):
        r = x / total
        assert (r * total, (1 - r) * total, r * (1 - r) * total) == (x, y, h)


def gap_ratio(sqrt_left, sqrt_right, w):
    # Arguments are exp(-t U/2) and exp(-t V/2). Their expectations use equal weights.
    gap = np.mean(w * sqrt_left**2 + (1 - w) * sqrt_right**2
                  - sqrt_left**(2 * w) * sqrt_right**(2 * (1 - w)))
    return gap / np.mean((sqrt_left - sqrt_right)**2)


@pytest.mark.parametrize("w", [.07, .32, .5, .73, .96])
def test_quantitative_jensen_deficit_for_equal_margins(w):
    rng = np.random.default_rng(782)
    sqrt_left = rng.uniform(.005, .97, size=100)
    for _ in range(12):
        sqrt_right = rng.permutation(sqrt_left)
        ratio = gap_ratio(sqrt_left, sqrt_right, w)
        assert ratio >= min(w, 1 - w) - 1e-13
        if w == .5:
            assert ratio == pytest.approx(.5, abs=1e-13)


@pytest.mark.parametrize("w", [.27, .73])
def test_cyclic_couplings_approach_sharp_unequal_weight_coefficient(w):
    # Equal marginals on a mesh have many small steps plus one large return step.
    mesh = np.linspace(1e-4, .9, 4097)
    successor = np.roll(mesh, -1)
    left, right = (mesh, successor) if w < .5 else (successor, mesh)
    coefficient = min(w, 1 - w)
    assert coefficient <= gap_ratio(left, right, w) < coefficient + .002


def condition(prior, likelihood):
    weights = {point: mass * likelihood(point) for point, mass in prior.items()
               if likelihood(point) > 0}
    evidence = sum(weights.values())
    return {point: mass / evidence for point, mass in weights.items()}


def marginal(prior, readout):
    weights = Counter()
    for point, mass in prior.items():
        weights[readout(point)] += mass
    return weights


def test_sharp_and_noisy_allocation_readouts_leave_the_within_ray_radial_law():
    low, high = Fraction(1, 4), Fraction(3, 4)
    prior = {(4, low): Fraction(1, 8), (8, low): Fraction(3, 8), (4, high): Fraction(1, 2)}
    exact = condition(prior, lambda p: p[1] == low)
    assert marginal(exact, lambda p: p[0]) == {4: low, 8: high}
    noisy = condition(prior, lambda p: high if p[1] == low else low)
    assert marginal(noisy, lambda p: p[1]) == {low: high, high: low}
    for r in (low, high):
        assert condition(noisy, lambda p: p[1] == r) == condition(prior, lambda p: p[1] == r)
    exact_total = condition(prior, lambda p: p[0] == 4)
    assert marginal(exact_total, lambda p: p[1]) == {low: Fraction(1, 5), high: Fraction(4, 5)}


def test_identical_total_and_allocation_marginals_can_give_different_conditional_masses():
    half, low, high = Fraction(1, 2), Fraction(1, 4), Fraction(3, 4)
    direct = {(4, low): half, (8, high): half}
    crossed = {(8, low): half, (4, high): half}
    for coordinate in [0, 1]:
        assert marginal(direct, lambda p: p[coordinate]) == marginal(crossed, lambda p: p[coordinate])
    sharp = lambda p: p[0] == 4
    mass_x = lambda p: p[0] * p[1]
    assert marginal(condition(direct, sharp), mass_x) == {1: 1}
    assert marginal(condition(crossed, sharp), mass_x) == {3: 1}


def median_interval(sample):
    ordered = sorted(sample)
    return ordered[(len(ordered) - 1) // 2], ordered[len(ordered) // 2]


@pytest.mark.parametrize("xs", [[1, 2, 5, 7], [1, 3, 4, 6, 9]])
def test_median_sets_reflect_on_a_sharp_total_and_dualize_on_a_sharp_reduced_mass(xs):
    total = 10
    ys = [total - x for x in xs]
    lo, hi = median_interval(xs)
    assert median_interval(ys) == (total - hi, total - lo)
    if lo == hi:
        med_y, med_y_upper = median_interval(ys)
        assert lo + med_y == total and med_y == med_y_upper
        inverse_x = [Fraction(1, x) for x in xs]
        inverse_y = [Fraction(1, y) for y in ys]
        med_inv_x, med_inv_y = median_interval(inverse_x)[0], median_interval(inverse_y)[0]
        assert set(reduced(x, y) for x, y in zip(inverse_x, inverse_y)) == {Fraction(1, total)}
        assert reduced(med_inv_x, med_inv_y) == Fraction(1, total)
    else:
        # Independently choosing both lower medians need not recover the measured sum.
        assert lo + median_interval(ys)[0] < total


@pytest.mark.parametrize("middle_sign", [-1, 0, 1])
def test_interval_total_median_error_can_attain_either_endpoint(middle_sign):
    total, epsilon = 10, Fraction(1, 4)
    xs = [1, 3, 5, 7, 9]
    ys = [total - x + sign * epsilon for x, sign in zip(xs, [-1, 1, middle_sign, -1, 1])]
    assert {min(x + y for x, y in zip(xs, ys)), max(x + y for x, y in zip(xs, ys))} == {
        total - epsilon, total + epsilon}
    median_x, median_y = median_interval(xs)[0], median_interval(ys)[0]
    assert median_x + median_y - total == middle_sign * epsilon


@pytest.mark.parametrize("pairs,one_ray,shared_ranks", [
    ([(1, 2), (2, 4), (3, 6)], True, True),
    ([(1, 2), (2, 4), (3, 5)], False, True),
    ([(1, 2), (2, 4), (3, 1)], False, False),
])
def test_strictly_positive_likelihood_preserves_ray_and_comonotone_support(pairs, one_ray, shared_ranks):
    prior = {pair: Fraction(1, len(pairs)) for pair in pairs}
    post = condition(prior, lambda p: 1 if p[1] == 2 * p[0] else Fraction(1, 10**9))
    assert set(prior) == set(post)
    for law in (prior, post):
        assert (len({Fraction(x, y) for x, y in law}) == 1) == one_ray
        discordance = sum(prob * other_prob for (x, y), prob in law.items()
                          for (a, b), other_prob in law.items() if (x - a) * (y - b) < 0)
        assert (discordance == 0) == shared_ranks
    if not one_ray:
        on_ray = sum(weight for (x, y), weight in post.items() if y == 2 * x)
        assert Fraction(999999, 1000000) < on_ray < 1
