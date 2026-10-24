"""Reproduce the Gaussian-conditioning comparison (3D, independent channels).

Research calculation, not a replacement estimator. Python + NumPy only.
Default prints checked results; --output optionally writes machine-readable JSON.
"""
import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

LAB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB))
from estimator import estimate, isotropic_covariance

SOURCE_COMMIT = '8782a293297330f7a354d87ebe62792018fa8866'
MEASURES = ('flat', 'tube', 'angle')


@lru_cache(None)
def gauss(order, upper):
    x, w = np.polynomial.legendre.leggauss(order)
    return (x+1)*upper/2, w*upper/2


def sphere_average(x, y, force, acceleration):
    """Average Gaussian likelihood over a uniform shared direction on S^2.

    x,y are standardized magnitudes. The angular exponential integrates to
    sinh(k)/k, k=|x F/sf + y a/sa|. Include the data constant for stability.
    """
    p, q = float(force @ force), float(acceleration @ acceleration)
    dot = float(force @ acceleration)
    k = np.sqrt(np.maximum(0., p*x*x + q*y*y + 2*dot*x*y))
    factor = np.ones_like(k)
    np.divide(-np.expm1(-2*k), 2*k, out=factor, where=k > 0)
    return np.exp(-.5*(x*x+y*y+p+q)+k)*factor


def summary(force, acceleration, measure, order=128, coordinates='polar', padding=14.):
    """Unit sigmas; callers multiply mass by sf/sa after standardizing data.

    Quadrature has caller-set finite cutoff and order, checked below on the
    printed cases by refinement and an independent Cartesian parameterization.
    """
    force, acceleration = np.array(force, dtype=float), np.array(acceleration, dtype=float)
    if force.shape != (3,) or acceleration.shape != (3,):
        raise ValueError('This comparison is specifically three dimensional')
    if measure not in MEASURES:
        raise ValueError('Unknown reference measure')
    upper = math.sqrt(float(force @ force + acceleration @ acceleration)) + padding
    if coordinates == 'polar':
        r, wr = gauss(order, upper)
        theta, wt = gauss(order, math.pi/2)
        x, y = r[:, None]*np.sin(theta), r[:, None]*np.cos(theta)
        weights = wr[:, None]*wt*r[:, None]
    elif coordinates == 'cartesian':
        t, w = gauss(order, upper)
        x, y = t[:, None], t[None, :]
        weights = w[:, None]*w[None, :]
    else:
        raise ValueError('Unknown integration coordinates')
    rho = 1. if measure == 'flat' else x*x+y*y if measure == 'tube' else (x*y)**2
    weighted = weights*rho*sphere_average(x, y, force, acceleration)
    z = float(weighted.sum())
    mean_x, mean_y = float((weighted*x).sum()/z), float((weighted*y).sum()/z)
    return {'normalization': z, 'mean_x': mean_x, 'mean_y': mean_y, 'mass_over_s': mean_x/mean_y}


def null_results():
    # Angular density: p(theta)=4/pi sin^2(2 theta).
    # Its CDF is 2 theta/pi - sin(4 theta)/(2 pi).
    lo, hi = 0., math.pi/2
    for _ in range(65):
        mid = (lo+hi)/2
        if 2*mid/math.pi-math.sin(4*mid)/(2*math.pi) < .025:
            lo = mid
        else:
            hi = mid
    low = math.tan((lo+hi)/2)
    return {
        'flat': {'interval95': [math.tan(math.pi/80), 1/math.tan(math.pi/80)],
                 'log_sd': math.pi/2, 'mean_x': math.sqrt(2/math.pi)},
        'tube': {'interval95': [math.tan(math.pi/80), 1/math.tan(math.pi/80)],
                 'log_sd': math.pi/2, 'mean_x': 1.5*math.sqrt(2/math.pi)},
        'angle': {'interval95': [low, 1/low], 'log_sd': math.sqrt(math.pi**2/4-2),
                  'mean_x': 2*math.sqrt(2/math.pi), 'mean_ratio': 4/math.pi,
                  'second_ratio_moment': 3.}}


def correlated_null_results():
    """Separate sensitivity: retain flat magnitudes, vary joint covariance."""
    ell, weights = gauss(256, 40.)
    results = {}
    for correlation in (-.8, 0., .8):
        beta = math.pi/2+math.asin(correlation)
        a = math.sqrt(1-correlation**2)/beta
        density = a/(2*(np.cosh(ell)-correlation))
        expected = (math.pi**2-beta**2)/3
        assert abs(float(2*weights @ density)-1) < 2e-11
        assert abs(float(2*weights @ (ell**2*density))-expected) < 2e-11
        # For l>=40 and |c|<=.8, density <= 2a exp(-l).
        tail_bound = 4*a*math.exp(-40)*(40**2+2*40+2)
        results[str(correlation)] = {'log_variance': expected,
                                     'omitted_log_variance_bound': tail_bound}
    return results


def run():
    cases = {
        'zero_zero': ([0, 0, 0], [0, 0, 0]),
        'force_only': ([3, 0, 0], [0, 0, 0]),
        'oblique': ([4, 0, 0], [1, 1, 0]),
        'anti_aligned': ([4, 0, 0], [-1, 0, 0]),
        'strong_aligned': ([8, 0, 0], [4, 0, 0]),
    }
    results, worst = {}, 0.
    null = null_results()
    for name, (force, acceleration) in cases.items():
        results[name] = {}
        for measure in MEASURES:
            low = summary(force, acceleration, measure, 128)
            high = summary(force, acceleration, measure, 256)
            independent = summary(force, acceleration, measure, 256, 'cartesian')
            extended = summary(force, acceleration, measure, 256, padding=18.)
            for alternate in (low, independent, extended):
                for key in high:
                    error = abs(alternate[key]/high[key]-1)
                    worst = max(worst, error)
                    assert error < 2e-10, (name, measure, key, error)
            swapped = summary(acceleration, force, measure, 128)
            assert abs(high['mass_over_s']*swapped['mass_over_s']-1) < 2e-10
            if name == 'zero_zero':
                assert abs(high['mean_x']/null[measure]['mean_x']-1) < 2e-11
                assert abs(high['mass_over_s']-1) < 2e-11
            results[name][measure] = high
        base = estimate(force, acceleration, isotropic_covariance(3, 1., 1.),
                        direction_order=32, ratio_order=512)
        assert abs(results[name]['flat']['mass_over_s']/base.mass-1) < 2e-7

    # Independently integrate the angular limit's known density and moments.
    theta, w = gauss(512, math.pi/2)
    density = 4/math.pi*np.sin(2*theta)**2
    for power, expected in ((0, 1.), (1, 4/math.pi), (2, 3.)):
        assert abs(float(w @ (density*np.tan(theta)**power))-expected) < 2e-11
    # Log-density integration avoids logarithmic endpoint convergence issues.
    ell, wl = gauss(256, 30.)
    log_density = 2/math.pi/np.cosh(ell)**3
    assert abs(float(2*wl @ log_density)-1) < 2e-11
    assert abs(float(2*wl @ (ell**2*log_density))-(math.pi**2/4-2)) < 2e-11
    # Explicit two-sided omitted second moment bound from sech(l)^3 <= 8e^-3l.
    log_tail_bound = 32/math.pi*math.exp(-90)*(30**2/3+2*30/9+2/27)
    return {'source_commit': SOURCE_COMMIT, 'cases': results, 'null': null,
            'correlated_flat_null': correlated_null_results(),
            'worst_relative_quadrature_difference': worst,
            'angular_log_variance_tail_bound': log_tail_bound,
            'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'baseline_sha256': hashlib.sha256((LAB/'estimator.py').read_bytes()).hexdigest()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    print('case                 flat          tube         angle     (mass / s)')
    for name, row in result['cases'].items():
        print(f"{name:20s} " + '  '.join(f"{row[m]['mass_over_s']:12.8f}" for m in MEASURES))
    print('Null summaries:', json.dumps(result['null']))
    print('Worst relative quadrature difference:', result['worst_relative_quadrature_difference'])
    print('All assertions passed.')
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
