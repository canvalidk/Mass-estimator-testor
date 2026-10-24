"""Exact 1D reduction for three radial references in the isotropic 3D model.

Research checks only; does not alter the baseline estimator. Uses NumPy and
the standard library. --output optionally saves the reproducible results.
"""
import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from gaussian_conditioning import SOURCE_COMMIT, summary
from weighting_checks import infer_batch

MEASURES = ('flat', 'tube', 'candidate_covariance')


@lru_cache(None)
def nodes(order):
    t, w = np.polynomial.legendre.leggauss(order)
    return (t+1)/2, w/2


def radial_integrals(theta, force, acceleration):
    """Return scaled J1,J3,J3+2J1 and their radial first-moment integrals.

    Scaling is common to all angles and both sets of integrals. J1 denotes
    integral r exp(-r²/2) sinh(rh)/(rh) dr; J3 has radial power three.
    """
    force, acceleration = np.asarray(force, float), np.asarray(acceleration, float)
    p, q, c = float(force@force), float(acceleration@acceleration), float(force@acceleration)
    sn, cs = np.sin(theta), np.cos(theta)
    h2 = np.maximum(0., p*sn**2+q*cs**2+2*c*sn*cs)
    top = .5*(p+q+math.sqrt((p-q)**2+4*c*c)) if c >= 0 else max(p, q)
    h = np.sqrt(h2)
    erf = np.fromiter((math.erf(float(x)/math.sqrt(2)) for x in h.flat),
                      dtype=float, count=h.size).reshape(h.shape)
    factor = np.ones_like(h)
    np.divide(math.sqrt(math.pi/2)*erf, h, out=factor, where=h > 1e-10)
    exp = np.exp((h2-top)/2)
    j1 = exp*factor
    j3 = (1+h2)*j1+math.exp(-top/2)
    j2 = math.sqrt(math.pi/2)*exp
    j4 = (h2+3)*j2
    return np.stack((j1, j3, j3+2*j1)), np.stack((j2, j4, j4+2*j2))


def infer(force, acceleration, ratio=1., order=192):
    t, wt = nodes(order)
    cut = math.atan(ratio)
    theta = np.r_[t*cut, cut+t*(math.pi/2-cut)]
    weight = np.r_[wt*cut, wt*(math.pi/2-cut)]
    dens, first = radial_integrals(theta, force, acceleration)
    z = dens@weight
    cdf = (dens[:, :order]@weight[:order])/z
    point = (first@(weight*np.sin(theta)))/(first@(weight*np.cos(theta)))
    return point, cdf


def quantile(force, acceleration, probability, measure_index, order=192):
    low, high = -36., 36.
    _, left = infer(force, acceleration, math.exp(low), order)
    _, right = infer(force, acceleration, math.exp(high), order)
    if not left[measure_index] < probability < right[measure_index]:
        raise ValueError('Quantile lies outside the research calculation bracket')
    for _ in range(55):
        mid = (low+high)/2
        _, cdf = infer(force, acceleration, math.exp(mid), order)
        if cdf[measure_index] < probability:
            low = mid
        else:
            high = mid
    return math.exp((low+high)/2)


def log_moments(force, acceleration, order=768, cutoff=36.):
    t, wt = nodes(order)
    ell, weight = cutoff*(2*t-1), 2*cutoff*wt
    theta = np.arctan(np.exp(ell))
    dens, _ = radial_integrals(theta, force, acceleration)
    mass_weight = dens*(weight/(2*np.cosh(ell)))
    z = mass_weight.sum(1)
    means = (mass_weight@ell)/z
    variances = (mass_weight@(ell*ell))/z-means*means
    return means, variances


def run():
    rng = np.random.default_rng(20260919)
    cases = [([0,0,0], [0,0,0]), ([4,0,0], [-1,0,0]),
             ([4,0,0], [1,1,0]), ([8,0,0], [4,0,0])]
    cases += [(rng.normal(size=3)*3, rng.normal(size=3)*3) for _ in range(8)]
    worst_point, worst_cdf = 0., 0.
    for force, acceleration in cases:
        point, cdf = infer(force, acceleration, 2.3)
        point2, cdf2 = infer(force, acceleration, 2.3, 320)
        independent_point, independent_cdf = infer_batch([force], [acceleration], 2.3, 160, 128)
        worst_point = max(worst_point, float(np.max(abs(point/point2-1))),
                          float(np.max(abs(point[:2]/independent_point[0,:2]-1))))
        worst_cdf = max(worst_cdf, float(np.max(abs(cdf-cdf2))),
                        float(np.max(abs(cdf[:2]-independent_cdf[0,:2]))))
        flat = summary(force, acceleration, 'flat', 256)
        tube = summary(force, acceleration, 'tube', 256)
        expected = ((tube['normalization']*tube['mean_x']+2*flat['normalization']*flat['mean_x']) /
                    (tube['normalization']*tube['mean_y']+2*flat['normalization']*flat['mean_y']))
        worst_point = max(worst_point, abs(point[2]/expected-1))
        expected_cdf = ((tube['normalization']*independent_cdf[0,1]+2*flat['normalization']*independent_cdf[0,0]) /
                        (tube['normalization']+2*flat['normalization']))
        worst_cdf = max(worst_cdf, abs(cdf[2]-expected_cdf))
        inverse, icdf = infer(acceleration, force, 1/2.3)
        assert np.max(abs(point*inverse-1)) < 2e-10
        assert np.max(abs(cdf+icdf-1)) < 2e-10
    assert worst_point < 2e-9, worst_point
    assert worst_cdf < 2e-9, worst_cdf

    matched = {}
    for name, cosine in [('aligned',1.), ('perpendicular',0.), ('opposed',-1.)]:
        force, acceleration = [3.,0,0], [3*cosine,3*math.sqrt(1-cosine*cosine),0]
        points, cdfs = infer(force, acceleration)
        means, variances = log_moments(force, acceleration)
        _, refined_variances = log_moments(force, acceleration, 1024)
        _, extended_variances = log_moments(force, acceleration, 1280, 44.)
        assert np.max(abs(variances-refined_variances)) < 2e-10
        assert np.max(abs(variances-extended_variances)) < 2e-10
        assert np.max(abs(points-1)) < 2e-12
        assert np.max(abs(cdfs-.5)) < 2e-12
        assert np.max(abs(means)) < 2e-12
        row = {}
        for j, measure in enumerate(MEASURES):
            lo = quantile(force, acceleration, .025, j)
            hi = quantile(force, acceleration, .975, j)
            assert abs(lo*hi-1) < 2e-10
            row[measure] = {'point':float(points[j]), 'interval95':[lo,hi],
                            'log_variance':float(variances[j])}
        if cosine > 0:
            assert variances[1] < variances[2] < variances[0]
        elif cosine < 0:
            assert variances[0] < variances[2] < variances[1]
        else:
            assert np.max(abs(variances-math.pi**2/4)) < 2e-10
            assert max(abs(row[m]['interval95'][0]/math.tan(math.pi/80)-1) for m in MEASURES) < 2e-10
        matched[name] = row

    weak = []
    for delta in (.04, .02, .01):
        force, acceleration = np.array([2.,0,0])*delta, np.array([1.,.5,0])*delta
        means, variances = log_moments(force, acceleration)
        point, _ = infer(force, acceleration)
        expected_means = delta**2*2.75*np.array([1/6,1/3,1/4])
        expected_variances = math.pi**2/4-delta**2*2*math.pi*np.array([1/9,2/9,1/6])
        expected_points = 1+delta**2*2.75*np.array([1/6,5/18,7/30])
        errors = np.r_[means-expected_means, variances-expected_variances, point-expected_points]
        assert np.max(abs(errors)) < 3*delta**4
        weak.append({'scale':delta, 'max_absolute_second_order_remainder':float(np.max(abs(errors))),
                     'log_means':means.tolist(), 'log_variances':variances.tolist(), 'points':point.tolist()})

    return {'source_commit':SOURCE_COMMIT,
            'independent_cases':len(cases), 'worst_relative_point_difference':worst_point,
            'worst_absolute_cdf_difference':worst_cdf, 'matched_observed_norm_3':matched,
            'weak_data_expansion':weak,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'surface_integrator_sha256':hashlib.sha256(Path(__file__).with_name('gaussian_conditioning.py').read_bytes()).hexdigest(),
            'cdf_integrator_sha256':hashlib.sha256(Path(__file__).with_name('weighting_checks.py').read_bytes()).hexdigest()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    print(json.dumps(result, indent=2))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print('All assertions passed.')
