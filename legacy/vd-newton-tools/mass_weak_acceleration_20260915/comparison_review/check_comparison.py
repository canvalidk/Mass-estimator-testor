"""Independent checks of the one-sided scalar benchmark; NumPy only."""
import json
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np

ETA = .05
Z = NormalDist().inv_cdf(1-ETA)
CDF = np.vectorize(lambda x: .5*math.erfc(-float(x)/math.sqrt(2)))


def lower(x, y, sf, sa, c):
    """Upward hull infimum, retaining empty-set infimum as infinity."""
    if x < Z*sf:
        return 0.
    if x == Z*sf and y >= Z*c/sf:
        return 0.
    if y <= -Z*sa:
        return math.inf
    def g(m):
        return x-m*y-Z*math.sqrt(sf*sf+sa*sa*m*m-2*c*m)
    hi = sf/sa
    while g(hi) > 0:
        hi *= 2
    lo = 0.
    for _ in range(100):
        mid = (lo+hi)/2
        if g(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def hull_coverage_formula(m, a, sf, sa, c):
    s = math.sqrt(sf*sf+m*m*sa*sa-2*m*c)
    r = (sf*sf-m*c)/(sf*s)
    nodes, weights = np.polynomial.legendre.leggauss(200)
    t = Z+(nodes+1)*(12-Z)/2
    values = np.exp(-t*t/2)/math.sqrt(2*math.pi)
    values *= CDF((Z-m*a/sf-r*t)/math.sqrt(1-r*r))
    return 1-ETA+float(weights@values)*(12-Z)/2


def flat_known_force_coverage(delta):
    """Exact numerical coverage of flat posterior lower mass bound."""
    # U_flat(t) = delta iff Phi(t-delta)/Phi(t) = ETA.
    def h(t):
        return math.log(.5*math.erfc((delta-t)/math.sqrt(2)))-math.log(.5*math.erfc(-t/math.sqrt(2)))-math.log(ETA)
    lo, hi = -30., 30.
    assert h(lo) < 0 < h(hi)
    for _ in range(100):
        mid=(lo+hi)/2
        if h(mid) < 0:
            lo=mid
        else:
            hi=mid
    threshold=(lo+hi)/2
    return {'acceleration_snr':delta,'threshold_y_over_sigma':threshold,
            'coverage':float(CDF(delta-threshold))}


def run():
    rng=np.random.default_rng(20260915)
    errors=[]
    for _ in range(1000):
        sf,sa=np.exp(rng.uniform(-2,2,2))
        c=rng.uniform(-.98,.98)*sf*sa
        x=(Z+rng.uniform(.001,10))*sf
        expected=(c+math.sqrt(c*c+sa*sa*(x*x/(Z*Z)-sf*sf)))/(sa*sa)
        observed=lower(x,0.,sf,sa,c)
        errors.append(abs(expected-observed)/expected)
    assert max(errors)<1e-11
    # Correlated boundary: x=z*sf and y=0 need not have lower endpoint zero.
    assert abs(lower(Z,0.,1.,1.,.5)-1.)<1e-12
    assert lower(Z,0.,1.,1.,-.5)==0
    assert math.isinf(lower(Z, -2.,1.,1.,.5))
    samples=600000
    coverage=[]
    for m,a,c in [(1.,.1,0.),(1.,1.,0.),(1.,4.,0.),(.2,1.,.7),(5.,.1,-.7),(1.,.1,.7)]:
        data=rng.multivariate_normal([m*a,a],[[1,c],[c,1]],samples)
        x,y=data.T
        t=(x-m*y)/math.sqrt(1+m*m-2*m*c)
        exactset=float(np.mean(t<=Z))
        hull=float(np.mean(~((x>Z)&(t>Z))))
        formula=hull_coverage_formula(m,a,1,1,c)
        assert abs(exactset-.95)<.0015
        assert abs(hull-formula)<.0015
        coverage.append({'m':m,'a':a,'c':c,'set_mc':exactset,
                         'hull_mc':hull,'hull_formula':formula})
    output={'zero_y_root_max_relative_error':max(errors),
            'mc_trials_per_case':samples,'coverage':coverage,
            'known_force_flat_coverage':[flat_known_force_coverage(d) for d in [.25,.5,1.,2.,5.]],
            'zero_y_known_force_unit_scale':{'confidence_lower':1/Z,
             'flat_lower':1/NormalDist().inv_cdf(.975),
             'ratio_of_means_point':math.sqrt(math.pi/2)}}
    destination=Path(__file__).with_name('checks.json')
    destination.write_text(json.dumps(output,indent=2)+'\n',encoding='utf8')
    print(json.dumps(output,indent=2))


if __name__=='__main__':
    run()
