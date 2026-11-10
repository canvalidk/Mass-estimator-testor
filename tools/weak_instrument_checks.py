"""Deterministic B24 checks: conditional critical values and set geometry.

Run: python tools/weak_instrument_checks.py
No random draws, coverage simulation, production estimator, or Wishart integral.
The orientation checks partition [0, pi) at every boundary crossing of the
positive-angle window, then classify each intervening open interval exactly
up to scalar root/quadrature accuracy. Isolated accepted points are retained.
"""

from functools import lru_cache
from math import acos, asin, atan, cos, exp, pi, sin, sqrt

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import chi2


P = 0.95
ALPHA = 1 - P
L = pi / 2
Q = float(chi2.ppf(P, 1))


def close(actual, expected, label, atol=2e-8):
    assert abs(actual - expected) <= atol, (label, actual, expected)


def clr_tail(c, a, k):
    """P(U/c + V/(a+c) > 1); substitute U=c*t^2 in the integral."""
    if c <= 0:
        return 1.0
    if k == 1:
        return float(chi2.sf(c, 1))
    integral = quad(
        lambda t: exp(-c * t * t / 2)
        * chi2.sf((a + c) * (1 - t * t), k - 1),
        0, 1, epsabs=2e-11, epsrel=2e-11,
    )[0]
    return float(chi2.sf(c, 1) + sqrt(2 * c / pi) * integral)


@lru_cache(maxsize=None)
def clr_critical(a, k):
    upper = float(chi2.ppf(P, k))
    if k == 1 or a == 0:
        return upper
    return brentq(lambda c: clr_tail(c, a, k) - ALPHA,
                  Q, upper, xtol=2e-11)


def unique(values, tolerance=1e-9):
    answer = []
    for value in sorted(values):
        if not answer or value - answer[-1] > tolerance:
            answer.append(value)
    return answer


def boundaries(f):
    """Independent sign-change roots; include exact grid tangencies."""
    grid = np.linspace(0, pi, 129)
    values = [f(t) for t in grid]
    roots = []
    for t, value in zip(grid[:-1], values[:-1]):
        if abs(value) < 1e-11:
            roots.append(float(t))
    for left, right, fl, fr in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if fl * fr < 0:
            roots.append(brentq(f, left, right, xtol=1e-13))
    return unique(roots)


def orientation_probabilities(f):
    """Uniform rotation, positive window (0,L), including zero-length pieces."""
    roots = boundaries(f)
    cuts = unique([0.0, pi] + [(-t) % pi for t in roots]
                  + [(L - t) % pi for t in roots])
    result = dict(whole=0., empty=0., disconnected=0., triples=0.,
                  endpoints_low=0., endpoints_high=0.)
    for left, right in zip(cuts[:-1], cuts[1:]):
        phi = (left + right) / 2
        weight = (right - left) / pi
        inner = unique([(phi + t) % pi for t in roots
                        if 1e-9 < (phi + t) % pi < L - 1e-9])
        edges = [0.] + inner + [L]
        # Alternate open segments and their shared boundary points. This
        # counts the minor-axis isolated point in rank-one score inversion.
        states = []
        for index, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
            states.append(f((a + b) / 2 - phi) >= 0)
            if index < len(inner):
                states.append(f(b - phi) >= -1e-8)
        components = sum(value and (i == 0 or not states[i - 1])
                         for i, value in enumerate(states))
        result['whole'] += weight * all(states)
        result['empty'] += weight * (components == 0)
        result['disconnected'] += weight * (components >= 2)
        result['triples'] += weight * (components >= 3)
        endpoints = [f(-phi) >= 0, f(L - phi) >= 0]
        result['endpoints_low'] += weight * all(endpoints)
        result['endpoints_high'] += weight * (not any(endpoints))
    return result


def compare_shapes(actual, expected, label):
    for name in ('whole', 'empty', 'disconnected', 'triples'):
        close(actual[name], expected.get(name, 0.), f'{label}: {name}')


def single_arc(delta):
    wide = max(2 * delta / pi - .5, 0.)
    return dict(whole=wide, empty=max(.5 - 2 * delta / pi, 0.),
                disconnected=wide)


def invariance_checks():
    data = np.array([[1., -2.], [3., .5], [-.7, 1.2], [2., 4.]])

    def statistics(y, theta):
        x = y @ np.array([sin(theta), cos(theta)])
        w = y @ np.array([cos(theta), -sin(theta)])
        a, b, d = float(x @ x), float(w @ w), float(x @ w)
        lr = (b - a + sqrt((a - b) ** 2 + 4 * d * d)) / 2
        return dict(score=d * d / a if a else 0., a=a, ar=b, d_squared=d * d,
                    clr=lr, clr_tail=clr_tail(lr, a, len(y)))

    def same_statistics(actual, expected, label):
        for name in actual:
            close(actual[name], expected[name], f'{label}: {name}')

    for theta in (.23, .71, 1.1):
        reference = statistics(data, theta)
        close(statistics(np.vstack([data, [0., 0.]]), theta)['score'],
              reference['score'], 'score Z')
        same_statistics(statistics(data[[2, 0, 3, 1]], theta), reference, 'order')
        same_statistics(statistics(data[:, ::-1], L - theta), reference, 'swap')
    close(statistics(np.zeros((3, 2)), .71)['score'], 0., 'score zero convention')
    sigmas, mass = np.array([.8, 1.7]), .9
    physical = data * sigmas
    scales = np.array([2.3, .6])
    theta = atan(mass / (sigmas[0] / sigmas[1]))
    changed_sigmas = sigmas * scales
    changed_mass = mass * scales[0] / scales[1]
    changed_theta = atan(changed_mass / (changed_sigmas[0] / changed_sigmas[1]))
    close(changed_theta, theta, 'unit conversion angle')
    same_statistics(statistics(physical * scales / changed_sigmas, changed_theta),
                    statistics(physical / sigmas, theta), 'physical units')
    print('AR/CLR/score: order, channel-swap and physical-units checks pass; '
          'score zero-padding and all-zero convention pass.')


def zero_padding_checks():
    # Gram diag(A,B), unchanged by a zero pair; K nevertheless increases.
    a, b, d = 1., 7., 0.
    lr = (b - a + sqrt((a - b) ** 2 + 4 * d * d)) / 2
    critical = {k: clr_critical(a, k) for k in (2, 3)}
    assert b > chi2.ppf(P, 2) and b < chi2.ppf(P, 3), 'AR Z fixture'
    assert lr > critical[2] and lr < critical[3], 'CLR Z fixture'
    for k, c in critical.items():
        close(clr_tail(c, a, k), ALPHA, f'CLR critical K={k}', 2e-10)
        close(clr_tail(chi2.ppf(P, k), 0., k), ALPHA,
              f'CLR chi-square limit K={k}', 2e-10)
    print(f'Zero-padding counterexamples: AR B={b:g}, q2={chi2.ppf(P, 2):.9f}, '
          f'q3={chi2.ppf(P, 3):.9f}; CLR LR={lr:g}, '
          f'c2(1)={critical[2]:.9f}, c3(1)={critical[3]:.9f}.')


def lm_shapes():
    for u, v in [(2., 1.), (10., 1.), (40., .5), (5., 0.), (12., 0.)]:
        def f(t):
            a = u * cos(t) ** 2 + v * sin(t) ** 2
            d = (u - v) * sin(t) * cos(t)
            return Q * a - d * d

        actual = orientation_probabilities(f)
        if (sqrt(u) - sqrt(v)) ** 2 <= Q:
            expected = dict(whole=1.)
        else:
            discriminant = (u + v - Q) ** 2 - 4 * u * v
            rminus = (u + v - Q - sqrt(discriminant)) / 2
            rplus = (u + v - Q + sqrt(discriminant)) / 2
            a = acos(sqrt(np.clip((rplus - v) / (u - v), 0., 1.)))
            b = asin(sqrt(np.clip((rminus - v) / (u - v), 0., 1.)))
            whole = (max(2 * a - L, 0.) + max(2 * b - L, 0.)) / pi
            expected = dict(whole=whole, disconnected=2 * (a + b) / pi - whole,
                            triples=whole)
        compare_shapes(actual, expected, f'LM eigenvalues {(u, v)}')
        if v == 0 and u == 5:
            assert actual['triples'] > .1, 'isolated minor point was dropped'
    print('LM: five orientation kernels pass, including rank-one triples.')


def clr_shapes():
    for k, u, v in [(2, 2., 1.), (2, 12., 1.), (2, 20., 5.), (3, 12., 1.)]:
        def f(t):
            a = u * cos(t) ** 2 + v * sin(t) ** 2
            return clr_tail(max(u - a, 0.), a, k) - ALPHA

        actual = orientation_probabilities(f)
        if u <= v + clr_critical(v, k):
            expected = dict(whole=1.)
        else:
            threshold = brentq(lambda a: a + clr_critical(a, k) - u, v, u)
            delta = acos(sqrt((threshold - v) / (u - v)))
            expected = single_arc(delta)
        compare_shapes(actual, expected, f'CLR K={k}, eigenvalues {(u, v)}')
    print('CLR: four orientation kernels pass against direct conditional-tail roots.')


def ar_shapes():
    for k, u, v in [(2, 3., 1.), (2, 10., 7.), (2, 10., 1.),
                    (3, 8., 1.), (3, 20., 1.)]:
        q = float(chi2.ppf(P, k))
        f = lambda t: q - (u * sin(t) ** 2 + v * cos(t) ** 2)
        actual = orientation_probabilities(f)
        expected = (dict(whole=1.) if u <= q else dict(empty=1.) if v > q
                    else single_arc(asin(sqrt((q - v) / (u - v)))))
        compare_shapes(actual, expected, f'AR K={k}, eigenvalues {(u, v)}')
        # Averaging these identities over a central Wishart Gram gives the
        # stated null formulas: independent diagonal chi-squares replace
        # endpoints_low by P^2 and endpoints_high by (1-P)^2.
        close(actual['whole'], (actual['endpoints_low'] + (u <= q)) / 2, 'AR whole')
        close(actual['empty'], (actual['endpoints_high'] + (v > q)) / 2, 'AR empty')
        close(actual['disconnected'], (actual['endpoints_low'] - (u <= q)) / 2,
              'AR disconnected')
    print('AR: five orientation kernels and positive-null probability identities pass.')


if __name__ == '__main__':
    zero_padding_checks()
    invariance_checks()
    lm_shapes()
    clr_shapes()
    ar_shapes()
    print('All deterministic B24 checks passed; no coverage study was run.')
