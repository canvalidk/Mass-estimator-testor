"""Focused repeated-sampling calibration for the three 3D compatible measures."""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'03_trace_and_evaluator/mass_estimator/investigation'))
from gaussian_conditioning import summary, null_results

MEASURES = ('flat', 'tube', 'angle')


@lru_cache(None)
def nodes(n):
    t, w = np.polynomial.legendre.leggauss(n)
    return (t+1)/2, w/2


def infer_batch(force, acceleration, true_ratio, nr=96, nt=64, padding=12.):
    """Point and CDF at the true mass, unit channel sigmas, shared direction.

    A central 95% credible interval covers the fixed true mass precisely when
    this CDF lies between .025 and .975. Splitting the angular quadrature at
    atan(true_ratio) avoids estimating quantiles from a fixed interpolant.
    """
    force, acceleration = np.asarray(force, float), np.asarray(acceleration, float)
    if force.shape != acceleration.shape or force.ndim != 2 or force.shape[1] != 3:
        raise ValueError('Expected two matching arrays of three-dimensional readings')
    if not 0 < true_ratio < math.inf:
        raise ValueError('Positive finite mass ratio required')
    p, q = (force**2).sum(1), (acceleration**2).sum(1)
    dot = (force*acceleration).sum(1)
    z, wz = nodes(nr)
    t, wt = nodes(nt)
    split = math.atan(true_ratio)
    theta = np.r_[t*split, split+t*(math.pi/2-split)]
    wtheta = np.r_[wt*split, wt*(math.pi/2-split)]
    upper = np.sqrt(p+q)+padding
    r = upper[:, None, None]*z[None, :, None]
    x, y = r*np.sin(theta), r*np.cos(theta)
    k = np.sqrt(np.maximum(0., p[:, None, None]*x*x+q[:, None, None]*y*y
                           +2*dot[:, None, None]*x*y))
    factor = np.ones_like(k)
    np.divide(-np.expm1(-2*k), 2*k, out=factor, where=k > 0)
    base = (np.exp(-.5*(r*r+p[:, None, None]+q[:, None, None])+k)*factor
            *r*upper[:, None, None]*wz[None, :, None]*wtheta)
    points, cdfs = [], []
    for rho in (1., r*r, (x*y)**2):
        w = base*rho
        norm = w.sum((1, 2))
        points.append((w*x).sum((1, 2))/(w*y).sum((1, 2)))
        cdfs.append(w[:, :, :nt].sum((1, 2))/norm)
    return np.stack(points, 1), np.stack(cdfs, 1)


def infer(force, acceleration, true_ratio, nr=96, nt=64):
    points, cdfs = [], []
    for start in range(0, len(force), 32):
        point, cdf = infer_batch(force[start:start+32], acceleration[start:start+32],
                                true_ratio, nr, nt)
        points.append(point)
        cdfs.append(cdf)
    return np.concatenate(points), np.concatenate(cdfs)


def wilson(successes, n):
    p, z = successes/n, 1.959963984540054
    center = (p+z*z/(2*n))/(1+z*z/n)
    radius = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return [center-radius, center+radius]


def checks():
    for force, acceleration, true_ratio in (([0,0,0],[0,0,0],4.),
                                             ([4,0,0],[1,1,0],2.),
                                             ([64,0,0],[4,0,0],16.)):
        point, cdf = infer_batch([force], [acceleration], true_ratio)
        point2, cdf2 = infer_batch([force], [acceleration], true_ratio, 128, 96)
        assert np.max(np.abs(point/point2-1)) < 1e-7, (force, acceleration, point, point2)
        assert np.max(np.abs(cdf-cdf2)) < 1e-7, (force, acceleration, cdf, cdf2)
        inverse, icdf = infer_batch([acceleration], [force], 1/true_ratio, 128, 96)
        assert np.max(np.abs(point2*inverse-1)) < 1e-10
        assert np.max(np.abs(cdf2+icdf-1)) < 1e-10
        for j, measure in enumerate(MEASURES):
            reference = summary(force, acceleration, measure, 256)
            assert abs(point2[0,j]/reference['mass_over_s']-1) < 1e-10
        if force == [0,0,0]:
            angle = math.atan(true_ratio)
            expected = [2*angle/math.pi, 2*angle/math.pi,
                        2*angle/math.pi-math.sin(4*angle)/(2*math.pi)]
            assert np.max(np.abs(cdf2[0]-expected)) < 1e-11


def run(n=2048, seed=20260915):
    checks()
    rng = np.random.default_rng(seed)
    # Common random numbers across scenarios and measures; iid within a scenario.
    noise = rng.normal(size=(n, 2, 3))
    results = []
    max_point_error, max_cdf_error = 0., 0.
    for mass in (1., 4., 16.):
        for excitation in (0., .25, 1., 4.):
            force, acceleration = noise[:,0].copy(), noise[:,1].copy()
            force[:,0] += mass*excitation
            acceleration[:,0] += excitation
            points, cdfs = infer(force, acceleration, mass)
            # Check a fixed subset plus observations nearest either coverage boundary.
            proximity = np.min(np.minimum(abs(cdfs-.025), abs(cdfs-.975)), axis=1)
            selected = np.unique(np.r_[np.arange(min(32,n)), np.argsort(proximity)[:32]])
            refined_points, refined_cdfs = infer(force[selected], acceleration[selected], mass, 128, 96)
            ep = float(np.max(abs(points[selected]/refined_points-1)))
            ec = float(np.max(abs(cdfs[selected]-refined_cdfs)))
            max_point_error, max_cdf_error = max(max_point_error,ep), max(max_cdf_error,ec)
            assert ep < 1e-6 and ec < 1e-6, (mass,excitation,ep,ec)
            flags = (cdfs[selected] >= .025)&(cdfs[selected] <= .975)
            refined_flags = (refined_cdfs >= .025)&(refined_cdfs <= .975)
            assert np.array_equal(flags, refined_flags), (mass,excitation,'classification')
            measures = {}
            for j, measure in enumerate(MEASURES):
                covered = (cdfs[:,j] >= .025)&(cdfs[:,j] <= .975)
                count = int(covered.sum())
                log_error = np.log(points[:,j]/mass)
                measures[measure] = {
                    'coverage': count/n, 'covered': count, 'wilson95': wilson(count,n),
                    'median_point_over_truth': float(np.median(points[:,j]/mass)),
                    'rms_log_point_error': float(np.sqrt(np.mean(log_error**2))),
                    'truth_below_interval': float(np.mean(cdfs[:,j] < .025)),
                    'truth_above_interval': float(np.mean(cdfs[:,j] > .975))}
            row = {'mass_over_s':mass, 'acceleration_snr':excitation, 'measures':measures}
            results.append(row)
            print(f"r={mass:4g}, alpha/sa={excitation:4g}: " + ', '.join(
                f"{m} coverage={measures[m]['coverage']:.4f}, median/truth={measures[m]['median_point_over_truth']:.4f}"
                for m in MEASURES), flush=True)
    return {'seed':seed, 'replicates_per_scenario':n, 'dimension':3,
            'sigma_force':1., 'sigma_acceleration':1., 'shared_direction_known_to_estimator':False,
            'max_relative_point_refinement_error':max_point_error,
            'max_absolute_cdf_refinement_error':max_cdf_error,
            'results':results}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replicates', type=int, default=2048)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run(args.replicates)
    print('All assertions passed.', flush=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
