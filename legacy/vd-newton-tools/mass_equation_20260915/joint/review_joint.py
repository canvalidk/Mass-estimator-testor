"""Direct-vector audit of the high-N score simulation, without posterior rerun."""
from pathlib import Path
import math
from statistics import NormalDist
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'03_trace_and_evaluator/mass_estimator/investigation'))
from joint_checks import coverage_flags, posterior


def run():
    rng = np.random.default_rng(20260916)
    nrep = 8192
    scenarios = ((1,16.,.25),(1,2.,1.),(8,2.,1.),(64,2.,1.),(512,2.,1.))
    zcut = NormalDist().inv_cdf(.975)
    zvalues, centers, quadratic_covered, grams, standardized_centers = [], [], [], [], []
    for trials, r, excitation in scenarios:
        k = 3*trials
        for start in range(0, nrep, 32):
            count = min(32,nrep-start)
            x,y = rng.normal(size=(2,count,trials,3))
            if trials != 512:
                continue
            for i in range(trials):
                x[:,i,i%3] += r*excitation
                y[:,i,i%3] += excitation
            xv,yv = x.reshape(count,k),y.reshape(count,k)
            residual = (xv-r*yv)/math.sqrt(1+r*r)
            fitted = (r*xv+yv)/math.sqrt(1+r*r)
            z = np.sum(residual*fitted,axis=1)/np.linalg.norm(fitted,axis=1)
            zvalues.extend(z.tolist())
            a,b,c = np.sum(xv*xv,axis=1),np.sum(yv*yv,axis=1),np.sum(xv*yv,axis=1)
            flags = coverage_flags(a,b,c,k,r)['score']
            assert np.array_equal(flags,abs(z)<=zcut)
            root = (a-b+np.sqrt((a-b)**2+4*c*c))/(2*c)
            theta_error = np.arctan(root)-math.atan(r)
            curvature = np.sqrt((a-b)**2+4*c*c)
            grams.extend(np.stack([a,b,c],axis=1).tolist())
            standardized_centers.extend((theta_error*np.sqrt(curvature)).tolist())
            centers.extend((math.sqrt(trials)*theta_error).tolist())
            quadratic_covered.extend((abs(theta_error)*np.sqrt(curvature)<=zcut).tolist())
    zvalues,centers = np.array(zvalues),np.array(centers)
    grams,standardized_centers = np.array(grams),np.array(standardized_centers)
    selected = np.unique(np.r_[np.argsort(abs(abs(standardized_centers)-zcut))[:32],
                                np.argsort(abs(standardized_centers))[-16:]])
    abc=grams[selected].T
    point,cdf=posterior(*abc,1536,2.)
    point_hi,cdf_hi=posterior(*abc,1536,2.,128,192)
    assert np.max(abs(cdf-cdf_hi)) < 1e-8
    assert np.max(abs(point/point_hi-1)) < 1e-8
    assert np.array_equal((cdf>=.025)&(cdf<=.975),(cdf_hi>=.025)&(cdf_hi<=.975))
    return {
        'direct_vector_score_coverage':float(np.mean(abs(zvalues)<=zcut)),
        'score_mean':float(zvalues.mean()),
        'score_variance':float(zvalues.var(ddof=1)),
        'profile_quadratic_coverage':float(np.mean(quadratic_covered)),
        'N_times_variance_theta_center':float(centers.var(ddof=1)),
        'asymptotic_N_times_sampling_variance_theta':8/25,
        'asymptotic_N_times_posterior_variance_theta':1/5,
        'asymptotic_credible_coverage':2*NormalDist().cdf(zcut/math.sqrt(1.6))-1,
        'additional_posterior_checks':len(selected),
        'max_added_absolute_cdf_discrepancy':float(np.max(abs(cdf-cdf_hi))),
        'max_added_relative_point_discrepancy':float(np.max(abs(point/point_hi-1))),
    }


if __name__=='__main__':
    print(run())
    print('All direct-vector score flags agree with Gram implementation.')
