"""Checks for weak-acceleration mass information, with NumPy only.

Research companion to weak_acceleration_information.md. The existing estimator
is unchanged. --output optionally writes the complete results to a chosen path.
Known-direction sampling comparisons use the original flat magnitude measure.
"""
import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist
import sys

import numpy as np

LAB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB))
from estimator import estimate, positive_normal_mean, radial_logs
from operational_checks import infer as vector_infer, quantile as vector_quantile

SOURCE_COMMIT = '8782a293297330f7a354d87ebe62792018fa8866'
NORMAL = NormalDist()
Z95 = NORMAL.inv_cdf(.95)


def phi(x):
    return np.exp(-np.asarray(x)**2/2)/math.sqrt(2*math.pi)


def Phi(x):
    x = np.asarray(x, dtype=float)
    return np.fromiter((.5*math.erfc(-float(t)/math.sqrt(2)) for t in x.flat),
                       float, count=x.size).reshape(x.shape)


@lru_cache(None)
def nodes(n):
    t, w = np.polynomial.legendre.leggauss(n)
    return (t+1)/2, w/2


def wilson(rate,n):
    z=NORMAL.inv_cdf(.975)
    center=(rate+z*z/(2*n))/(1+z*z/n)
    radius=z*math.sqrt(rate*(1-rate)/n+z*z/(4*n*n))/(1+z*z/n)
    return [center-radius,center+radius]


def bisect(function, target, lo, hi, steps=70):
    """Increasing continuous function, with a verified bracket."""
    assert function(lo) <= target <= function(hi), (lo, hi, target)
    for _ in range(steps):
        mid = (lo+hi)/2
        if function(mid) < target:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def alpha_quantile(t, p=.95):
    # Survival form avoids subtracting a tiny normal-tail probability from one.
    if t>-8:
        return t-NORMAL.inv_cdf((1-p)*float(Phi(t)))
    def logcdf(x):
        if x>-36:
            return math.log(.5*math.erfc(-x/math.sqrt(2)))
        # Mills expansion in a region where its first terms converge rapidly.
        term=total=1.
        for k in range(1,40):
            term*=-(2*k-1)/(x*x)
            total+=term
            if abs(term)<1e-17:
                break
        return -.5*x*x-math.log(-x)-.5*math.log(2*math.pi)+math.log(total)
    base=logcdf(t)
    return bisect(lambda q:base-logcdf(t-q),-math.log1p(-p),0.,12.)


def known_force_results():
    rows = []
    for theta in (.05, .25, 1., 2., 4., 8.):
        threshold = bisect(alpha_quantile, theta, -100., max(10., theta+5.))
        coverage = float(Phi(theta-threshold))
        # h(t) is increasing: exact central sampling quantiles of mhat/mtrue.
        ratios = [theta/positive_normal_mean(theta-NORMAL.inv_cdf(p), 1.)
                  for p in (.05, .5, .95)]
        assert coverage >= .95-2e-14
        rows.append({'true_acceleration_over_sigma':theta,
                     'lower_bound_coverage':coverage,
                     'point_over_true_mass_sampling_quantiles_05_50_95':ratios})
    zero = {'point_over_F_over_sigma':math.sqrt(math.pi/2),
            'posterior_lower95_over_F_over_sigma':1/alpha_quantile(0.),
            'gaussian_lower95_over_F_over_sigma':1/Z95}
    # Independent positive-axis quadrature checks the truncated mean/quantile.
    t, w = nodes(512)
    worst = 0.
    for reading in (-4., -1., 0., 1., 4.):
        upper = max(0., reading)+14
        a = upper*t
        density = phi(a-reading)/float(Phi(reading))
        mean = float(upper*w@(a*density))
        worst = max(worst, abs(mean/positive_normal_mean(reading,1.)-1))
        q = alpha_quantile(reading)
        prob = float(q*w@(phi(q*t-reading)/float(Phi(reading))))
        worst = max(worst, abs(prob-.95))
    assert worst < 3e-11, worst
    return {'zero_reading':zero, 'coverage_and_point':rows,
            'worst_independent_integral_error':worst}


def censor_cdf(theta, k):
    """P(alpha/sigma <= theta | only the bin |Y/sigma| <= k is reported)."""
    def J(v):
        return v*float(Phi(v))+float(phi(v))
    return (k-J(k-theta)+J(-k-theta))/k


def coarsened_coverage(theta, k, p=.95):
    qbin = bisect(lambda v:censor_cdf(v,k), p, 0., k+12.)
    tstar = bisect(lambda t:alpha_quantile(t,p), theta, -100., max(10.,theta+5.))
    # Outside the bin, exact readings are reported and used.
    above = float(Phi(theta-max(k,tstar)))
    below = max(0., float(Phi(-k-theta)-Phi(tstar-theta))) if tstar < -k else 0.
    inbin = float(Phi(k-theta)-Phi(-k-theta)) if theta <= qbin else 0.
    return above+below+inbin


def censoring_results(rng, n=300000):
    rows = []
    t,w = nodes(512)
    for k in (.1, 1., 2.):
        mean = ((k*k+1)*(float(Phi(k))-.5)+k*float(phi(k)))/k
        q = bisect(lambda a:censor_cdf(a,k), .95, 0., k+12.)
        a = (k+14)*t
        likelihood = Phi(k-a)-Phi(-k-a)
        assert abs(float((k+14)*w@likelihood)/k-1) < 3e-12
        assert abs(float((k+14)*w@(a*likelihood))/k/mean-1) < 3e-12
        just_above = q+1e-5
        exact_coverage = coarsened_coverage(just_above,k)
        readings = just_above+rng.normal(size=n)
        inbin = abs(readings)<=k
        # Exact-reading posterior CDF evaluated at true acceleration; coverage
        # means that true acceleration is below its posterior upper quantile.
        cdf_at_truth = 1-Phi(readings-just_above)/Phi(readings)
        covers = np.where(inbin, just_above<=q, cdf_at_truth<=.95)
        simulation = float(covers.mean())
        se = math.sqrt(exact_coverage*(1-exact_coverage)/n)
        assert abs(simulation-exact_coverage) < 5*se+1/n
        rows.append({'nondetection_halfwidth_over_sigma':k,
                     'mean_acceleration_over_sigma':mean,
                     'point_over_F_over_sigma':1/mean,
                     'posterior_upper95_acceleration_over_sigma':q,
                     'posterior_lower95_mass_over_F_over_sigma':1/q,
                     'coverage_right_limit':float(Phi(q-k)),
                     'tested_theta':just_above,'exact_coverage':exact_coverage,
                     'simulation_coverage':simulation, 'simulation_n':n,
                     'safe_gaussian_lower95_over_F_over_sigma':1/(k+Z95)})
    return rows


def scalar_batch(force, acceleration, mass, correlation=0., order=64):
    """Known positive direction; unit marginal errors, full 2x2 covariance.

    Returns original readout and posterior CDF at the supplied mass. Analytic
    radial integration plus Gauss quadrature split exactly at the CDF boundary.
    """
    x,y = np.broadcast_arrays(np.asarray(force,float),np.asarray(acceleration,float))
    x,y = x.ravel(),y.ravel()
    if not (mass>0 and abs(correlation)<1):
        raise ValueError('Positive mass and positive definite covariance required')
    t,w = nodes(order)
    cut=math.atan(mass)
    theta=np.r_[cut*t,cut+(math.pi/2-cut)*t]
    weight=np.r_[cut*w,(math.pi/2-cut)*w]
    sn,cs=np.sin(theta),np.cos(theta)
    c=correlation
    aa=(1-2*c*sn*cs)/(1-c*c)
    bb=((x[:,None]-c*y[:,None])*sn+(y[:,None]-c*x[:,None])*cs)/(1-c*c)
    j1,j2=radial_logs((bb/np.sqrt(aa)).ravel())
    lz=j1.reshape(bb.shape)-np.log(aa)
    lr=j2.reshape(bb.shape)-1.5*np.log(aa)
    peak=np.maximum(lz.max(1),lr.max(1))[:,None]
    z,r=np.exp(lz-peak)*weight,np.exp(lr-peak)*weight
    return (r@sn)/(r@cs),z[:,:order].sum(1)/z.sum(1)


def confidence_lower(force, acceleration, correlation=0.):
    """Infimum of {m>0: (x-my)/sqrt(1+m²-2cm)<=z95}.

    0 means no positive lower restriction; inf marks an empty candidate set,
    never a physical infinite-mass estimate. Complete set may be disconnected.
    """
    x,y=np.broadcast_arrays(np.asarray(force,float),np.asarray(acceleration,float))
    result=np.zeros(x.shape)
    initial_rejection=(x>Z95)|((x==Z95)&(y<Z95*correlation))
    empty=initial_rejection&(y<=-Z95)
    result[empty]=math.inf
    active=initial_rejection&~empty
    xx,yy=x[active],y[active]
    low=np.zeros_like(xx)
    high=np.ones_like(xx)
    def g(m):
        return xx-m*yy-Z95*np.sqrt(1+m*m-2*correlation*m)
    for _ in range(80):
        grow=g(high)>0
        if not grow.any():
            break
        high[grow]*=2
    else:
        raise ArithmeticError('Confidence lower endpoint was not bracketed')
    for _ in range(70):
        mid=(low+high)/2
        reject=g(mid)>0
        low=np.where(reject,mid,low)
        high=np.where(reject,high,mid)
    result[active]=(low+high)/2
    return result


def scalar_checks():
    worst_point=worst_cdf=0.
    for c in (-.8,0.,.8):
        for x,y,m in ((0.,0.,1.),(4.,0.,4.),(3.,-.5,6.),(4.,1.,4.)):
            point,cdf=scalar_batch([x],[y],m,c,96)
            point2,cdf2=scalar_batch([x],[y],m,c,192)
            baseline=estimate([x],[y],[[1,c],[c,1]],known_direction=[1.],ratio_order=2048)
            worst_point=max(worst_point,abs(float(point[0]/point2[0]-1)),
                            abs(float(point[0]/baseline.mass-1)))
            worst_cdf=max(worst_cdf,abs(float(cdf[0]-cdf2[0])))
            inv,icdf=scalar_batch([y],[x],1/m,c,96)
            assert abs(point[0]*inv[0]-1)<2e-10
            assert abs(cdf[0]+icdf[0]-1)<2e-10
            if c==0:
                exact=positive_normal_mean(x,1.)/positive_normal_mean(y,1.)
                worst_point=max(worst_point,abs(float(point[0]/exact-1)))
                # Independent Cartesian integral over f: P(alpha>=f/m|y).
                t,w=nodes(512)
                upper=max(0.,x)+14
                f=upper*t
                cart=float(upper*w@(phi(f-x)/float(Phi(x))*Phi(y-f/m)/float(Phi(y))))
                worst_cdf=max(worst_cdf,abs(float(cdf[0]-cart)))
            lower=float(confidence_lower(np.array([x]),np.array([0.]),c)[0])
            expected=(c+math.sqrt(c*c+x*x/Z95**2-1)) if x>Z95 else 0.
            assert abs(lower-expected)<2e-12
    assert worst_point<2e-8, worst_point
    assert worst_cdf<2e-10, worst_cdf
    # Boundary x=z95 and disconnected-set example deserve explicit checks.
    assert confidence_lower([Z95],[-.5],0.)[0]>0
    assert confidence_lower([Z95],[0.],0.)[0]==0
    assert math.isinf(confidence_lower([3.],[-2.],0.)[0])
    return {'worst_relative_point_difference':worst_point,'worst_absolute_cdf_difference':worst_cdf}


def sampling_results(rng,n):
    noise=rng.normal(size=(n,2))
    scenarios=[(8.,a,0.) for a in (.05,.25,1.,4.)]
    scenarios += [(3.,.5,c) for c in (-.8,0.,.8)]
    scenarios += [(1.,1.,0.),(.25,1.,0.)]
    rows=[]
    worst_point=worst_cdf=0.
    for f,a,c in scenarios:
        x=f+noise[:,0]
        y=a+c*noise[:,0]+math.sqrt(1-c*c)*noise[:,1]
        mass=f/a
        points,cdfs=[],[]
        for start in range(0,n,256):
            p,v=scalar_batch(x[start:start+256],y[start:start+256],mass,c)
            points.extend(p);cdfs.extend(v)
        points,cdfs=np.asarray(points),np.asarray(cdfs)
        p2,v2=scalar_batch(x[:128],y[:128],mass,c,128)
        worst_point=max(worst_point,float(np.max(abs(points[:128]/p2-1))))
        worst_cdf=max(worst_cdf,float(np.max(abs(cdfs[:128]-v2))))
        lower=confidence_lower(x,y,c)
        pivot=(x-mass*y)/math.sqrt(1+mass*mass-2*c*mass)
        exact_set=pivot<=Z95
        hull=lower<=mass
        assert np.all(hull|~exact_set)
        set_rate=float(exact_set.mean())
        assert abs(set_rate-.95)<5*math.sqrt(.95*.05/n)
        row={'true_force_over_sigma':f,'true_acceleration_over_sigma':a,
             'correlation':c,'true_mass_in_instrument_units':mass,
             'posterior_lower95_coverage':float((cdfs>=.05).mean()),
             'one_sided_pivot_set_coverage':set_rate,
             'pivot_lower_hull_coverage':float(hull.mean()),
             'pivot_lower_hull_is_zero_rate':float((lower==0).mean()),
             'pivot_empty_set_rate':float(np.isinf(lower).mean()),
             'point_over_true_mass_median':float(np.median(points/mass)),
             'point_within_factor_two_rate':float(((points>=mass/2)&(points<=mass*2)).mean())}
        row['posterior_coverage_mc_wilson95']=wilson(row['posterior_lower95_coverage'],n)
        rows.append(row)
    assert worst_point<2e-7, worst_point
    assert worst_cdf<2e-7, worst_cdf
    return {'n_per_scenario':n,'rows':rows,'refinement_subset_per_scenario':128,
            'worst_relative_point_refinement':worst_point,'worst_absolute_cdf_refinement':worst_cdf}


def vector_results():
    # chi-square_3 .95, independently obtained from its exact elementary CDF.
    cdf=lambda q:math.erf(math.sqrt(q/2))-math.sqrt(2*q/math.pi)*math.exp(-q/2)
    cutoff=bisect(cdf,.95,0.,30.)
    rows=[]
    for f in (1.,3.,8.,20.):
        point,_=vector_infer([f,0,0],[0,0,0])
        lower=vector_quantile([f,0,0],[0,0,0],.05,0)
        refined=vector_quantile([f,0,0],[0,0,0],.05,0,320)
        assert abs(lower/refined-1)<2e-9
        residual_lower=math.sqrt(max(0.,f*f-cutoff)/cutoff)
        rows.append({'observed_force_over_sigma':f,
                     'flat_point_in_instrument_units':float(point[0]),
                     'flat_posterior_lower95':lower,
                     'vector_residual_set_lower_endpoint':residual_lower})
    return {'chi_square_3_95':cutoff,'observed_zero_acceleration_rows':rows}


def run(n=16384,seed=20260925):
    rng=np.random.default_rng(seed)
    results={'source_commit':SOURCE_COMMIT,'seed':seed,'known_force':known_force_results(),
             'censoring':censoring_results(rng),'scalar_checks':scalar_checks(),
             'sampling':sampling_results(rng,n),'unknown_direction_3d':vector_results()}
    results['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    results['baseline_sha256']=hashlib.sha256((LAB/'estimator.py').read_bytes()).hexdigest()
    results['runtime']={'python':sys.version,'numpy':np.__version__}
    dependencies=[LAB/'practical_formulas.py']+[Path(__file__).with_name(name) for name in
                   ('operational_checks.py','gaussian_conditioning.py','weighting_checks.py')]
    results['helper_sha256']={str(p.relative_to(LAB)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in dependencies}
    return results


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples',type=int,default=16384)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.samples<1024:
        parser.error('--samples must be at least 1024 for the sampling checks')
    result=run(args.samples)
    print(json.dumps(result,indent=2,allow_nan=False))
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('All assertions passed.')
