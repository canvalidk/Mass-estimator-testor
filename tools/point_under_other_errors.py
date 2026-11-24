"""The curve's least under errors that are not Gaussian, and on the whole arc of masses (T84; problems/O19.md).
The simulations behind definition/theorems/T84.md; tests/test_point_under_other_errors.py imports them.
Measurements, not a rule. Run: python tools/point_under_other_errors.py  (about five minutes)"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from met.evidence import Curve, leftover, span_of        # noqa: E402
from met.full import Reading                             # noqa: E402


# ---- the full type: mixed dimensions, declared subspaces, correlated covariances; an error is the
# ---- covariance's Cholesky factor times independent unit-variance draws of the named shape

def unit(rng, n, shape):
    if shape == "skewed":                                   # exponential, centred
        return rng.exponential(size=n) - 1.0
    if shape == "two-point":                                # -0.5 with probability 0.8, 2 otherwise
        return np.where(rng.random(n) < 0.8, -0.5, 2.0)
    if shape == "t5":
        return rng.standard_t(5, size=n) / math.sqrt(5.0 / 3.0)
    return rng.normal(size=n)


def declare(rng, n, mass):
    out = []
    for i in range(n):
        d = (3, 2, 1)[i % 3]
        k = int(rng.integers(1, d + 1))
        q, _ = np.linalg.qr(rng.normal(size=(d, d)))
        mix = rng.normal(size=(2 * d, 2 * d))
        cov = 0.6 * (mix @ mix.T) / (2 * d) + 0.4 * np.eye(2 * d)
        cov = 0.5 * (cov + cov.T)
        a = q[:, :k] @ rng.normal(size=k)
        out.append((d, q[:, :k], cov, np.concatenate([mass * a, a]), np.linalg.cholesky(cov)))
    return out


def true_readings(declared, factor=1.0):
    return [Reading(factor * y[:d], factor * y[d:], cov, subspace=sub) for d, sub, cov, y, _ in declared]


def observed(rng, declared, shape, factor=1.0):
    out = []
    for d, sub, cov, y, chol in declared:
        e = chol @ unit(rng, 2 * d, shape)
        out.append(Reading(factor * y[:d] + e[:d], factor * y[d:] + e[d:], cov, subspace=sub))
    return out


def expected_rise(shape, masses, reps, seed=84, n=30, mass=2.0, factor=0.8):
    """The noiseless rise of the curve from the true mass, and the mean and standard error of the observed rise."""
    rng = np.random.default_rng(seed)
    declared = declare(rng, n, mass)
    masses = np.asarray(masses, dtype=float)
    at_truth = int(np.argmin(np.abs(masses - mass)))
    noiseless = Curve(true_readings(declared, factor), scale=1.0)(masses)
    rises = []
    for _ in range(reps):
        values = Curve(observed(rng, declared, shape, factor), scale=1.0)(masses)
        rises.append(values - values[at_truth])
    rises = np.array(rises)
    return noiseless, rises.mean(0), rises.std(0) / math.sqrt(reps)


def normalised_errors(shape, n, reps, seed):
    """Median and 90% point of |angle of the least - true angle| times Delta / sqrt(D + Delta),
    with the signal energy set to Delta = 6 D^0.75 (so Delta / sqrt(D) grows slowly)."""
    rng = np.random.default_rng(seed)
    scaled = []
    for _ in range(reps):
        declared = declare(rng, n, 2.0)
        energy = sum(float(span_of(r).z @ span_of(r).z) for r in true_readings(declared))
        count = sum(sub.shape[1] for _, sub, *_ in declared)
        target = 6.0 * count ** 0.75
        factor = math.sqrt(target / energy)
        angle = leftover(observed(rng, declared, shape, factor), scale=1.0)["angle"]
        scaled.append(abs(angle - math.atan(2.0)) * target / math.sqrt(count + target))
    return float(np.median(scaled)), float(np.quantile(scaled, 0.9))


# ---- one reading on an instrument of ratio 1 with signal energy 8 log N at mass 1, and N noise-only
# ---- readings on an instrument of ratio N (sigma_F = N, sigma_a = 1); x, y are those in noise units

def far_instrument_curve(f1, a1, x, y, n, masses):
    m = np.asarray(masses, dtype=float)
    finite = np.isfinite(m)
    mm = np.where(finite, m, 1.0)
    first = np.where(finite, (f1 - mm * a1) ** 2 / (1.0 + mm * mm), a1 * a1)
    sxx, sxy, syy = x @ x, x @ y, y @ y
    rest = np.where(finite, (n * n * sxx - 2.0 * n * mm * sxy + mm * mm * syy) / (n * n + mm * mm), syy)
    return first + rest


def far_instrument_shares(n, trials, seed):
    """Share of series whose least on the whole arc is at a mass of n/2 or more, and share whose least on
    the interval 0.2 to 5 is within 0.5 of the true mass 1."""
    rng = np.random.default_rng(seed)
    size = math.sqrt(8.0 * math.log(n))
    masses = np.concatenate([[0.0], np.geomspace(1e-3, 1e3 * n, 6001), [math.inf]])
    inside = (masses >= 0.2) & (masses <= 5.0)
    gone = near = 0
    for _ in range(trials):
        f1, a1 = size / math.sqrt(2.0) + rng.normal(), size / math.sqrt(2.0) + rng.normal()
        curve = far_instrument_curve(f1, a1, rng.normal(size=n), rng.normal(size=n), n, masses)
        gone += masses[curve.argmin()] >= 0.5 * n
        near += abs(masses[inside][curve[inside].argmin()] - 1.0) < 0.5
    return gone / trials, near / trials


def far_instrument_conditions(n):
    """B27's weighted condition J^2 / (W + J), and the unweighted Delta / sqrt(D), on that series."""
    weighted_count = 0.25 + n * (n / (1.0 + n * n)) ** 2
    weighted_energy = 0.25 * 8.0 * math.log(n)
    return weighted_energy ** 2 / (weighted_count + weighted_energy), 8.0 * math.log(n) / math.sqrt(n + 1.0)


# ---- one instrument, scalar readings in noise units, true mass 2; errors Gaussian or a spherical
# ---- Student law standardised to the identity covariance (2.5 degrees: no fourth moment; 5: one)

MASSES = np.geomspace(0.05, 80.0, 801)
COS, SIN = np.cos(np.arctan(MASSES)), np.sin(np.arctan(MASSES))


def one_instrument_curve(f, a):
    return COS * COS * (f @ f) - 2.0 * SIN * COS * (f @ a) + SIN * SIN * (a @ a)


def errors(rng, n, degrees):
    if degrees is None:
        return rng.normal(size=(n, 2))
    return rng.normal(size=(n, 2)) * np.sqrt((degrees - 2.0) / rng.chisquare(degrees, size=n))[:, None]


def misses(degrees, n, power, coefficient, trials, seed, log_curve=False):
    """Share of series whose point is 0.5 or more from the true mass 2, at Delta = coefficient N^power.
    The point is the curve's least, or with log_curve the least of sum log(1 + K_i / 4)."""
    rng = np.random.default_rng(seed)
    along = np.array([2.0, 1.0]) / math.hypot(2.0, 1.0)
    size = math.sqrt(coefficient * n ** power / n)
    missed = 0
    for _ in range(trials):
        e = errors(rng, n, degrees)
        f, a = size * along[0] + e[:, 0], size * along[1] + e[:, 1]
        if log_curve:
            values = np.log1p((np.outer(f, COS) - np.outer(a, SIN)) ** 2 / 4.0).sum(0)
        else:
            values = one_instrument_curve(f, a)
        missed += not abs(MASSES[values.argmin()] - 2.0) < 0.5
    return missed / trials


if __name__ == "__main__":
    print("(1) the curve's expected rise from the true mass (30 full-type readings, 2000 series), at masses 0.3, 1, 4, 12")
    for shape in ("skewed", "two-point", "t5"):
        noiseless, mean, se = expected_rise(shape, [0.3, 1.0, 2.0, 4.0, 12.0], 2000)
        keep = [0, 1, 3, 4]
        print(f"   {shape:9s} noiseless {np.round(noiseless[keep], 2)}  observed {np.round(mean[keep], 2)}  +- {np.round(se[keep], 2)}")
    print("(2) |angle error| times Delta / sqrt(D + Delta), median (90% point), Delta = 6 D^0.75, 40 series a cell")
    for shape in ("gaussian", "skewed", "two-point"):
        cells = [normalised_errors(shape, n, 40, 1000 + n) for n in (60, 240, 720)]
        print(f"   {shape:9s} " + "   ".join(f"n={n}: {m:.2f} ({h:.2f})" for n, (m, h) in zip((60, 240, 720), cells)))
    print("(5) one reading at ratio 1, energy 8 log N, beside N noise-only readings at ratio N; 600 series")
    for n in (100, 500, 2500):
        away, near = far_instrument_shares(n, 600, 1720 + n)
        weighted, unweighted = far_instrument_conditions(n)
        print(f"   N={n:5d}  least on the whole arc at N/2 or more: {away:.3f}   on [0.2, 5] within 0.5 of the truth: {near:.3f}"
              f"   J^2/(W+J) = {weighted:.2f}   Delta/sqrt(D) = {unweighted:.2f}")
    print("(6) share of series whose point is 0.5 or more from the true mass 2; 300 series a cell")
    for label, power, coefficient in (("Delta = 8 N^0.65", 0.65, 8.0), ("Delta = 6 N^0.95", 0.95, 6.0)):
        for name, degrees in (("Gaussian", None), ("Student 5", 5.0), ("Student 2.5", 2.5)):
            least = [misses(degrees, n, power, coefficient, 300, 1730 + n + int(100 * power)) for n in (100, 400, 1600, 6400)]
            logs = [misses(degrees, n, power, coefficient, 300, 1790 + n + int(100 * power), log_curve=True) for n in (100, 400, 1600)]
            print(f"   {label}  {name:11s} the least: {np.round(least, 3)}   least of sum log(1 + K_i/4): {np.round(logs, 3)}")
