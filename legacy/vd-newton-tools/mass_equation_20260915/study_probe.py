"""First investigation checks; development copy kept out of the document space."""
from functools import lru_cache
import json
import math
from pathlib import Path
import sys

import numpy as np

LAB = Path(__file__).resolve().parents[2] / '03_trace_and_evaluator' / 'mass_estimator'
sys.path.insert(0, str(LAB))
from estimator import estimate, isotropic_covariance


@lru_cache(None)
def nodes(order):
    return np.polynomial.legendre.leggauss(order)


def integrate(function, lower, upper, order=192):
    x, w = nodes(order)
    t = lower + (upper-lower)*(x+1)/2
    return float((upper-lower)/2 * (w @ function(t)))


def null_summary(lam, order=192):
    shift = max(lam, 0)
    def weight(theta):
        return np.exp(lam*np.sin(2*theta)-shift)
    z = integrate(weight, 0, math.pi/2, order)
    def quantile(p):
        lo, hi = 0., math.pi/2
        for _ in range(64):
            mid = (lo+hi)/2
            if integrate(weight, 0, mid, order)/z < p:
                lo = mid
            else:
                hi = mid
        return math.tan((lo+hi)/2)
    def log_density(ell):
        sech = 1/np.cosh(ell)
        return np.exp(lam*sech-shift)*sech/(2*z)
    # Smooth integration split into unit intervals; remaining tails analytic bound.
    var = 2*sum(integrate(lambda ell: ell*ell*log_density(ell), i, i+1, 32)
                for i in range(50))
    point = (integrate(lambda th: np.sin(th)*weight(th), 0, math.pi/2, order)
             / integrate(lambda th: np.cos(th)*weight(th), 0, math.pi/2, order))
    qlo, median, qhi = [quantile(p) for p in (.025,.5,.975)]
    return dict(lam=lam, point=point, lower=qlo, median=median, upper=qhi,
                log_sd=math.sqrt(var), log_normalizer=math.log(z)+shift,
                log_tail_coefficient=-math.log(z)-shift)


def main():
    rows = [null_summary(x) for x in (-4,0,4,16,64)]
    refined = [null_summary(x,384) for x in (-4,0,4,16,64)]
    assert max(abs(r['point']-1) for r in rows)<1e-12
    assert max(abs(r['lower']*r['upper']-1) for r in rows)<1e-10
    assert max(abs(r['median']-1) for r in rows)<1e-12
    assert max(abs(r[k]/s[k]-1) for r,s in zip(rows,refined)
               for k in ('point','lower','upper','log_sd'))<1e-10
    base = rows[1]
    assert abs(base['lower']/math.tan(math.pi*.025/2)-1)<1e-12
    assert abs(base['log_sd']-math.pi/2)<1e-12
    actual = estimate([0,0],[0,0],isotropic_covariance(2,2,1),ratio_order=4096)
    assert abs(actual.mass-2)<1e-12
    assert abs(actual.quantile(.025)/(2*base['lower'])-1)<1e-10
    for row in rows:
        lam,zlog=row['lam'],row['log_normalizer']
        r=1e8
        scaled_tail=math.exp(lam*2*r/(1+r*r))*r*r/(1+r*r)
        assert abs(scaled_tail-1)<2e-6
    print(json.dumps(dict(rows=rows,checks='passed'),indent=2))


if __name__=='__main__':
    main()
