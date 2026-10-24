"""Reproduce weighting, pooling, and repeated-sampling calibration findings.

Python + NumPy. Research checks, not a replacement estimator or confidence rule.
Default: analytic checks and 2048 experiments per scenario. --checks-only skips
Monte Carlo. --output writes results only when explicitly provided.
"""
import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

LAB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB))
from gaussian_conditioning import summary
from estimator import positive_normal_mean
from practical_formulas import aligned_local_uncertainty

MEASURES = ('flat', 'tube', 'angle')
SOURCE_COMMIT = '8782a293297330f7a354d87ebe62792018fa8866'


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


def pooling_checks():
    expected = []
    for order in (64, 128, 256):
        z, w = np.polynomial.hermite.hermgauss(order)
        means = [float(w @ np.array([positive_normal_mean(t+math.sqrt(2)*v, 1.)
                                    for v in z])/math.sqrt(math.pi)) for t in (2., 1.)]
        expected.append(means)
    assert np.max(np.abs(np.array(expected)/expected[-1]-1)) < 1e-12
    f, a = expected[-1]
    assert abs(f/a-1.514700842103585) < 1e-12
    assert 1 < f/a < 2
    # Earlier audit's universal bracketing assertion has this exact counterexample.
    rom, regression, inverse_regression = Fraction(10,3), Fraction(19,5), Fraction(82,19)
    assert rom < regression < inverse_regression
    return {'true_force':2., 'true_acceleration':1., 'sigma_force':1.,
            'sigma_acceleration':1., 'true_mass':2., 'expected_fitted_force':f,
            'expected_fitted_acceleration':a, 'separate_fit_pooling_limit':f/a,
            'bracketing_counterexample':{'ratio_of_means':str(rom),
                                       'regression':str(regression),
                                       'inverse_regression':str(inverse_regression)}}


def fisher_checks():
    rng = np.random.default_rng(9152026)
    worst = 0.
    for d in (1, 2, 3):
        for _ in range(20):
            frame, _ = np.linalg.qr(rng.normal(size=(d,d)))
            u, tangent = frame[:,0], frame[:,1:]
            f, alpha, sf, sa = np.exp(rng.uniform(-1,1,size=4))
            m, a = f/alpha, alpha*u
            jac = np.zeros((2*d,d+1))
            jac[:d,0], jac[d:,1] = u, u
            jac[:d,2:], jac[d:,2:] = f*tangent, alpha*tangent
            precision = np.diag([sf**-2]*d+[sa**-2]*d)
            observed = np.sqrt(np.linalg.det(jac.T @ precision @ jac))
            expected = (f*f/sf**2+alpha*alpha/sa**2)**((d-1)/2)/(sf*sa)
            worst = max(worst, abs(observed/expected-1))
            assert abs(observed/expected-1) < 2e-12

            b = rng.normal(size=(2*d,2*d))
            covariance = b @ b.T+np.eye(2*d)
            mass_jac = np.zeros((2*d,d+1))
            mass_jac[:d,0] = a
            mass_jac[:d,1:] = m*np.eye(d)
            mass_jac[d:,1:] = np.eye(d)
            information = mass_jac.T @ np.linalg.solve(covariance,mass_jac)
            effective = information[0,0]-information[0,1:] @ np.linalg.solve(
                information[1:,1:],information[1:,0])
            contrast = np.c_[np.eye(d),-m*np.eye(d)]
            formula = float(a @ np.linalg.solve(contrast @ covariance @ contrast.T,a))
            local = aligned_local_uncertainty(m*a,a,covariance)
            errors = [abs(formula/effective-1), abs(local*math.sqrt(formula)-1)]
            worst = max(worst,*errors)
            assert max(errors) < 2e-12
    return {'cases':60, 'worst_relative_difference':worst}


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
    return {'pooling':pooling_checks(), 'fisher':fisher_checks()}


def run(n=2048, seed=20260915):
    analytic_checks = checks()
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
    return {'analytic_checks':analytic_checks, 'seed':seed, 'replicates_per_scenario':n, 'dimension':3,
            'sigma_force':1., 'sigma_acceleration':1., 'shared_direction_known_to_estimator':False,
            'max_relative_point_refinement_error':max_point_error,
            'max_absolute_cdf_refinement_error':max_cdf_error,
            'results':results}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replicates', type=int, default=2048)
    parser.add_argument('--checks-only', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.replicates < 64:
        parser.error('Use at least 64 replicates')
    result = {'analytic_checks':checks()} if args.checks_only else run(args.replicates)
    result['source_commit'] = SOURCE_COMMIT
    result['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['baseline_sha256'] = hashlib.sha256((LAB/'estimator.py').read_bytes()).hexdigest()
    result['integration_reference_sha256'] = hashlib.sha256(
        Path(__file__).with_name('gaussian_conditioning.py').read_bytes()).hexdigest()
    print('All assertions passed.', flush=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
