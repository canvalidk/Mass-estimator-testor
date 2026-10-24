"""Reproduce the first equation-investigation contribution.

Python + NumPy only. This is a sensitivity study, not a replacement estimator.
Default: run checks and print a compact table. Optional --output writes JSON.
All quoted probabilities are conditional on a chosen latent measure.
"""
from functools import lru_cache
from fractions import Fraction
import argparse
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
LAMBDAS = (-4., 0., 4., 16., 64.)


@lru_cache(None)
def nodes(order):
    return np.polynomial.legendre.leggauss(order)


def integrate(function, lower, upper, order=192):
    x, w = nodes(order)
    t = lower + (upper-lower)*(x+1)/2
    return float((upper-lower)/2 * (w @ function(t)))


def null_summary(lam, order=192):
    """Fixed-experiment null law rho=exp(lam*sin(2 theta)); tested -4..256.

    theta=atan((f/sigma_F)/(alpha/sigma_a)). Values of mass below are divided
    by s=sigma_F/sigma_a. Smooth log-density quadrature extends to |log r|=50;
    its omitted second moment is bounded explicitly in the returned result.
    Resolution is caller-controlled, not an adaptive convergence certificate.
    """
    if not math.isfinite(lam):
        raise ValueError('Finite lambda required')
    shift = max(lam, 0.)

    def weight(theta):
        return np.exp(lam*np.sin(2*theta)-shift)

    z = integrate(weight, 0., math.pi/2, order)
    if z <= 0:
        raise ArithmeticError('Quadrature failed to resolve the density')

    def quantile(probability):
        lo, hi = 0., math.pi/2
        for _ in range(64):
            mid = (lo+hi)/2
            if integrate(weight, 0., mid, order)/z < probability:
                lo = mid
            else:
                hi = mid
        return math.tan((lo+hi)/2)

    def log_density(ell):
        sech = 1/np.cosh(ell)
        return np.exp(lam*sech-shift)*sech/(2*z)

    limit = 50
    variance = 2*sum(integrate(lambda ell: ell*ell*log_density(ell), i, i+1, 32)
                     for i in range(limit))
    # For |ell|>=L, sech ell<=2 exp(-|ell|) and
    # exp(lam sech ell)<=exp(max(lam,0)*2 exp(-L)).
    log_z = math.log(z)+shift
    omitted_variance_bound = (2*math.exp(max(lam,0)*2*math.exp(-limit)-log_z-limit)
                              * (limit*limit+2*limit+2))
    point = (integrate(lambda th: np.sin(th)*weight(th), 0., math.pi/2, order)
             / integrate(lambda th: np.cos(th)*weight(th), 0., math.pi/2, order))
    lower, median, upper = [quantile(p) for p in (.025,.5,.975)]
    probability_check = 2*sum(integrate(log_density, i, i+1, 32) for i in range(limit))
    return dict(lambda_value=lam, point_over_scale=point, lower_over_scale=lower,
                median_over_scale=median, upper_over_scale=upper,
                log_standard_deviation=math.sqrt(variance), log_normalizer=log_z,
                log_tail_coefficient=-log_z,
                omitted_log_variance_bound=omitted_variance_bound,
                log_density_normalization=probability_check)


def cartesian_null_check(lam, order=256):
    """Independent original-coordinate integration over x,y in [0,12].

    Used only for moderate lambdas. This integral does not use polar/ratio
    density formulas. The Gaussian truncation is negligible for these cases;
    the origin's directional limit is tested by quadrature refinement.
    """
    t,w = nodes(order)
    x,weights = 6*(t+1),6*w
    xx,yy = x[:,None],x[None,:]
    density = np.exp(-.5*(xx*xx+yy*yy)+lam*2*xx*yy/(xx*xx+yy*yy)-max(lam,0))
    joint = weights[:,None]*weights[None,:]*density
    z = joint.sum()
    # E[radius]=sqrt(pi/2) for the separate Rayleigh radial law.
    return dict(normalizer_scaled=float(z), point=float((joint*xx).sum()/(joint*yy).sum()),
                mean_x=float((joint*xx).sum()/z))


def exact_examples():
    p = Fraction(1,2)
    laws = (((1,1),(9,9)), ((1,9),(9,1)))
    output = []
    for law in laws:
        point = sum(p*f for f,a in law)/sum(p*a for f,a in law)
        output.append(dict(pairs=law, point=str(point),
                           masses=[str(Fraction(f,a)) for f,a in law]))
        assert point == 1
    assert output[0]['masses'] == ['1','1']
    assert output[1]['masses'] == ['1/9','9']

    weights = (Fraction(3,4),Fraction(1,4))
    features = (Fraction(9,16),Fraction(1,16))
    force = (1,2)
    direction_results=[]
    for lam in (0,1,4,16):
        numerator = sum(p*(1+lam*t)*f for p,t,f in zip(weights,features,force))
        denominator = sum(p*(1+lam*t) for p,t in zip(weights,features))
        ratio = numerator/denominator
        assert ratio == Fraction(80+29*lam,64+28*lam)
        # Here alpha=1, hence Q_alpha=P and Cov(feature,M)=-3/32.
        assert ratio-Fraction(5,4) == lam*Fraction(-3,32)/(1+lam*Fraction(7,16))
        direction_results.append(dict(lambda_value=lam,point=str(ratio)))
    return dict(coupling_examples=output,direction_family=direction_results)


def run():
    rows=[null_summary(lam) for lam in LAMBDAS]
    refined=[null_summary(lam,384) for lam in LAMBDAS]
    keys=('point_over_scale','lower_over_scale','upper_over_scale','log_standard_deviation')
    refinement=max(abs(a[k]/b[k]-1) for a,b in zip(rows,refined) for k in keys)
    assert refinement<1e-10
    for row in rows:
        assert abs(row['point_over_scale']-1)<1e-12
        assert abs(row['median_over_scale']-1)<1e-12
        assert abs(row['lower_over_scale']*row['upper_over_scale']-1)<1e-10
        assert abs(row['log_density_normalization']-1)<1e-11
        assert row['omitted_log_variance_bound']<1e-16
        # Numerical illustration of the derived tail coefficient, not a proof
        # that an ordinary moment diverges (that conclusion is analytic).
        r=1e8
        coefficient_ratio=math.exp(row['lambda_value']*2*r/(1+r*r))*r*r/(1+r*r)
        assert abs(coefficient_ratio-1)<2e-6
    assert all(a['log_standard_deviation']>b['log_standard_deviation'] for a,b in zip(rows,rows[1:]))
    baseline=rows[1]
    assert abs(baseline['lower_over_scale']/math.tan(math.pi*.025/2)-1)<1e-12
    assert abs(baseline['log_standard_deviation']-math.pi/2)<1e-12

    existing=estimate([0,0],[0,0],isotropic_covariance(2,2,1),ratio_order=4096)
    assert abs(existing.mass-2)<1e-12
    assert abs(existing.quantile(.025)/(2*baseline['lower_over_scale'])-1)<1e-10

    cartesian=[]
    for lam in (-4.,0.,4.,16.):
        direct=cartesian_null_check(lam,512)
        doubled=cartesian_null_check(lam,1024)
        exact=null_summary(lam)
        z=math.exp(exact['log_normalizer']-max(lam,0))
        error=abs(doubled['normalizer_scaled']/z-1)
        assert error<1e-7
        assert abs(direct['normalizer_scaled']/doubled['normalizer_scaled']-1)<1e-6
        assert abs(doubled['point']-1)<1e-12
        cartesian.append(dict(lambda_value=lam,normalizer_relative_error=error))

    asymptotic=[]
    for lam in (64.,128.,256.):
        r=null_summary(lam,384)
        asymptotic.append(lam*r['log_standard_deviation']**2)
    assert all(a>b>1 for a,b in zip(asymptotic,asymptotic[1:]))
    assert abs(asymptotic[-1]-1)<.007
    files=[Path(__file__),LAB/'estimator.py']
    return dict(source_commit=SOURCE_COMMIT,
                scope='Measured zero vectors, independent isotropic Gaussian errors; conditional latent-measure sensitivity.',
                status='All assertions passed.',rows=rows,
                maximum_relative_quadrature_change=refinement,
                cartesian_checks=cartesian,
                lambda_times_log_variance=asymptotic,
                exact_examples=exact_examples(),
                sha256={str(p.relative_to(LAB)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='Optional JSON output path; parent must exist.')
    args=parser.parse_args()
    result=run()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('lambda   point/s   95% lower/s   95% upper/s   SD(log(M/s))')
    for r in result['rows']:
        print(f"{r['lambda_value']:6g} {r['point_over_scale']:9.6f} {r['lower_over_scale']:13.6f} "
              f"{r['upper_over_scale']:13.6f} {r['log_standard_deviation']:15.6f}")
    print(result['status'])
    print(f"Largest quadrature refinement change: {result['maximum_relative_quadrature_change']:.3g}")


if __name__=='__main__':
    main()
