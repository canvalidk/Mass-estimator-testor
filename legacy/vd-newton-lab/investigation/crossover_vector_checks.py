"""Reproduce the isotropic 3D transition study (Python and NumPy only).

Exact radial and direction integrals leave one angular quadrature. This is
research code for a declared model, not a replacement production estimator.
All comparisons use the same simulated readings; no trial is rejected.
"""
import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np

LAB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB))
from estimator import estimate, isotropic_covariance, positive_normal_mean

SOURCE_COMMIT = '8782a293297330f7a354d87ebe62792018fa8866'


@lru_cache(None)
def nodes(order):
    x, w = np.polynomial.legendre.leggauss(order)
    return (x+1)/2, w/2


def infer(force, acceleration, ratio=1., order=64):
    """Flat and tube points, plus CDF at ratio, in standardized 3D units.

    Tube uses (f^2+alpha^2) df d(alpha) dOmega. Channel sigmas are one.
    This measure is a declared sensitivity comparator, not the user's law.
    """
    force, acceleration = np.atleast_2d(force), np.atleast_2d(acceleration)
    assert force.shape == acceleration.shape and force.shape[1] == 3
    p, q = np.sum(force**2, axis=1), np.sum(acceleration**2, axis=1)
    dot = np.sum(force*acceleration, axis=1)
    x, w = nodes(order)
    cut = math.atan(ratio)
    theta = np.r_[cut*x, cut+(math.pi/2-cut)*x]
    weight = np.r_[cut*w, (math.pi/2-cut)*w]
    sn, cs = np.sin(theta), np.cos(theta)
    h2 = np.maximum(0., p[:, None]*sn**2+q[:, None]*cs**2
                    +2*dot[:, None]*sn*cs)
    # Maximum on the positive quadrant, for stable common scaling.
    peak = np.where(dot >= 0, (p+q+np.hypot(p-q, 2*dot))/2, np.maximum(p, q))
    exp = np.exp((h2-peak[:, None])/2)
    h = np.sqrt(h2)
    erf = np.fromiter((math.erf(float(z)/math.sqrt(2)) for z in h.flat),
                      float, count=h.size).reshape(h.shape)
    factor = np.ones_like(h)
    np.divide(math.sqrt(math.pi/2)*erf, h, out=factor, where=h > 1e-12)
    j1 = exp*factor
    j3 = (1+h2)*j1+np.exp(-peak[:, None]/2)
    first = (exp, (h2+3)*exp)
    density = (j1, j3)
    points = np.column_stack([(v@(weight*sn))/(v@(weight*cs)) for v in first])
    cdf = np.column_stack([(v[:, :order]@weight[:order])/(v@weight) for v in density])
    return points, cdf


def profile_point(force, acceleration):
    """Positive-slope profile optimum; 0/inf mean unattained boundaries.

    At p=q,dot=0 every positive mass minimizes the profile: return NaN to
    represent nonidentification. At p=q,dot<0 either boundary is optimal;
    infinity is a reporting convention, and still counted as a point failure.
    """
    p, q = np.sum(force**2, axis=1), np.sum(acceleration**2, axis=1)
    dot = np.sum(force*acceleration, axis=1)
    disc = np.hypot(p-q, 2*dot)
    result = np.where(p >= q, np.inf, 0.)
    good = dot > 0
    high = good & (p >= q)
    low = good & (p < q)
    result[high] = (p[high]-q[high]+disc[high])/(2*dot[high])
    result[low] = 2*dot[low]/(q[low]-p[low]+disc[low])
    result[(p == q) & (dot == 0)] = np.nan
    return result


def metrics(readout, truth):
    valid = (readout > 0) & np.isfinite(readout)
    success = (readout >= truth/2) & (readout <= 2*truth)
    result = {'finite_positive_fraction': float(np.mean(valid)),
              'within_factor_two': float(np.mean(success)),
              'factor_two_mcse': float(np.std(success, ddof=1)/math.sqrt(len(success)))}
    if valid.all():
        error = np.log(readout/truth)
        result.update(log_mse=float(np.mean(error**2)),
                      log_mse_mcse=float(np.std(error**2, ddof=1)/math.sqrt(len(error))),
                      mean_log_error=float(np.mean(error)),
                      median_ratio=float(np.median(readout/truth)),
                      q90_absolute_log_error=float(np.quantile(abs(error), .9)))
    else:
        result['log_mse'] = None  # Do not report risk conditional on success.
    return result, success


def crossover_identities():
    c = math.sqrt(2/math.pi)
    def correction(t):
        return math.exp(-t*t/2)/(math.sqrt(2*math.pi)*t*.5*math.erfc(-t/math.sqrt(2)))
    thresholds = []
    for tolerance in (.1,.05,.01,.001):
        lo, hi = .0001, 10.
        for _ in range(65):
            mid = (lo+hi)/2
            if correction(mid) > tolerance:
                lo = mid
            else:
                hi = mid
        thresholds.append({'tolerance':tolerance,'minimum_observed_snr':(lo+hi)/2})
    for tf in (.1,.3,1.,2.,5.,10.):
        for ta in (.1,.3,1.,2.,5.,10.):
            m = positive_normal_mean(tf,1)/positive_normal_mean(ta,1)
            assert min(1,tf/ta)-1e-13 <= m <= max(1,tf/ta)+1e-13
            assert abs(m/(tf/ta)-1) <= correction(min(tf,ta))+1e-12
    for t in np.linspace(-10,12,1001):
        h = positive_normal_mean(float(t),1)
        lam = h-t
        slope = (1-lam*h)/h
        assert 0 < slope <= 1
    curves = []
    x,w = nodes(192)
    for t in (0.,.25,.5,1.,2.,3.,5.,10.):
        # Scaled Dawson integral, avoiding large exponentials.
        dawson = t*float(w @ np.exp(((t*x)**2-t*t)/2))
        r3 = math.erf(t/math.sqrt(2))/dawson if t else c
        point,_ = infer([0,0,0],[t,0,0],order=128)
        assert abs(point[0,0]-c/r3) < 2e-12
        curves.append({'snr':t,'H':positive_normal_mean(t,1),'R3':r3})
    return {'thresholds':thresholds,'axis_curves':curves,
            'zero_log_sensitivity':(1-c*c)/c}


def checks():
    cases = [([0,0,0],[0,0,0]), ([8,0,0],[0,0,0]),
             ([3,0,0],[1,.5,0]), ([2,0,0],[0,1,0]),
             ([3,0,0],[-2,0,0]), ([1,2,0],[2,-1,1])]
    worst = 0.
    for force, acceleration in cases:
        point, _ = infer(force, acceleration, order=128)
        baseline = estimate(force, acceleration, np.eye(6), direction_order=40, ratio_order=1024)
        error = abs(point[0,0]/baseline.mass-1)
        worst = max(worst, error)
        assert error < 2e-7, (force, acceleration, error)
        reverse, _ = infer(acceleration, force, order=128)
        assert np.max(abs(point*reverse-1)) < 2e-12
    rows = []
    for f in (0., .5, 1., 2., 3., 5., 8.):
        for a in (0., .5, 1., 2., 3., 5., 8.):
            point, _ = infer([f,0,0], [a,0,0], order=128)
            scalar = positive_normal_mean(f,1)/positive_normal_mean(a,1)
            rows.append({'observed_force_snr':f, 'observed_acceleration_snr':a,
                         'flat_point':float(point[0,0]), 'tube_point':float(point[0,1]),
                         'known_direction_point':scalar})
    angle_rows = []
    for degrees in (0, 30, 60, 90, 120, 150, 180):
        phi = math.radians(degrees)
        point, _ = infer([4,0,0], [math.cos(phi),math.sin(phi),0], order=128)
        angle_rows.append({'degrees':degrees, 'flat_point':float(point[0,0]),
                           'tube_point':float(point[0,1]), 'norm_ratio':4.})
    assert np.isnan(profile_point(np.array([[1.,0,0]]),np.array([[0.,1,0]]))[0])
    return {'baseline_relative_error':worst, 'observed_aligned_path':rows,
            'crossover_identities':crossover_identities(),
            'observed_angle_path':angle_rows}


def run(n):
    rng = np.random.default_rng(2026091601)
    errors = rng.standard_normal((n,6))
    scenarios = [(8.,x) for x in (.25,.5,1.,2.,3.,5.,8.)]
    scenarios += [(4.,x) for x in (.5,1.,2.,3.)]
    scenarios += [(x,x) for x in (.25,.5,1.,2.,3.)]
    results = []
    refinement = 0.
    for f, a in scenarios:
        force, acceleration = errors[:,:3].copy(), errors[:,3:].copy()
        force[:,0] += f
        acceleration[:,0] += a
        truth = f/a
        points = np.empty((n,2))
        cdfs = np.empty((n,2))
        for start in range(0,n,512):
            end = min(n,start+512)
            points[start:end], cdfs[start:end] = infer(force[start:end], acceleration[start:end], truth)
        refined, rcdf = infer(force[:128], acceleration[:128], truth, order=128)
        refinement = max(refinement,float(np.max(abs(refined/points[:128]-1))),
                         float(np.max(abs(rcdf-cdfs[:128]))))
        nf, na = np.linalg.norm(force,axis=1), np.linalg.norm(acceleration,axis=1)
        floor = math.sqrt(2/math.pi)
        methods = {'flat_joint':points[:,0], 'tube_joint':points[:,1],
                   'norm_ratio':nf/na, 'norm_floor':np.maximum(nf,floor)/np.maximum(na,floor),
                   'positive_profile':profile_point(force,acceleration)}
        row = {'true_force_snr':f,'true_acceleration_snr':a,'true_mass':truth,'methods':{}}
        indicators = {}
        for name, values in methods.items():
            row['methods'][name], indicators[name] = metrics(values,truth)
            if name == 'positive_profile':
                # Gaussian noise has full support: boundary outcomes have
                # positive probability even if this finite sample misses them.
                row['methods'][name]['log_mse'] = None
                row['methods'][name]['population_log_risk'] = 'infinite for boundary-valued rule'
        for name in ('flat_joint','tube_joint'):
            idx = 0 if name == 'flat_joint' else 1
            covered = (cdfs[:,idx] >= .025) & (cdfs[:,idx] <= .975)
            row['methods'][name]['posterior_central_95_coverage'] = float(np.mean(covered))
            row['methods'][name]['coverage_mcse'] = float(np.std(covered,ddof=1)/math.sqrt(n))
        row['paired_factor_two_gains'] = {}
        for name in methods:
            if name == 'flat_joint':
                continue
            diff = indicators['flat_joint'].astype(float)-indicators[name].astype(float)
            row['paired_factor_two_gains'][name] = {'gain':float(np.mean(diff)),
                                                  'mcse':float(np.std(diff,ddof=1)/math.sqrt(n))}
        results.append(row)
    assert refinement < 2e-8, refinement
    return {'n_per_scenario':n,'seed':2026091601,'scenarios':results,
            'max_refinement_difference':refinement}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples',type=int,default=32768)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = {'managed_commit':SOURCE_COMMIT,'python':platform.python_version(),'numpy':np.__version__,
              'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'estimator_sha256':hashlib.sha256((LAB/'estimator.py').read_bytes()).hexdigest(),
              'checks':checks(), 'sampling':run(args.samples)}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'assertions':'passed','baseline_relative_error':result['checks']['baseline_relative_error'],
                      'refinement_difference':result['sampling']['max_refinement_difference'],
                      'samples_per_scenario':args.samples,'scenario_count':len(result['sampling']['scenarios'])},indent=2))
