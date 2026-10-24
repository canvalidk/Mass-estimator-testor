"""Check precise scalar links to published normal-mean and ratio constructions.

Liseo (2003), equations (10), (12), (13): condition the flat-original-means
ratio posterior on positive ratio and change the dimensional scale. This
conditioning is OUR specialization, not a claimed formula printed by Liseo.
Katz formula: as reproduced in Chang et al.'s 2017 manuscript, equation (1.4).
No claim that either source selects the full VD vector ratio-of-means readout.
"""
import json
import math
from pathlib import Path

import numpy as np

from estimator import estimate, isotropic_covariance


def cdf(x):
    return .5 * math.erfc(-x/math.sqrt(2))


def katz(x, sigma):
    z=x/sigma
    return x + sigma*math.exp(-z*z/2)/math.sqrt(2*math.pi)/cdf(z)


def liseo_positive_angle(theta, zf, za):
    # Density in angle theta=atan(m/(sf/sa)), not density in mass m.
    # The ratio Jacobian sec(theta)^2 cancels Liseo's 1/(1+r^2).
    s=za*np.cos(theta)+zf*np.sin(theta)
    delta=np.array([.5*math.erf(x/math.sqrt(2)) for x in s])
    phi=np.exp(-s*s/2)/math.sqrt(2*math.pi)
    same_sign=cdf(zf)*cdf(za)+cdf(-zf)*cdf(-za)
    return math.exp(-.5*(zf*zf+za*za))/math.pi*(1+s*delta/phi)/same_sign


def run():
    rows=[]
    # Positive, negative, mixed-sign, one-zero, and both-zero readings; two
    # dimensional uncertainty ratios. Values avoid cancellation in closed forms.
    for zf,za in [(0.,0.),(0.,2.),(2.,0.),(3.,2.),(-2.,-3.),(1.,-1.),(-2.,.5)]:
        for sf,sa in [(1.,1.),(.7,1.3)]:
            f,a=sf*zf,sa*za
            unknown=estimate([f],[a],isotropic_covariance(1,sf,sa),ratio_order=4096)
            published=liseo_positive_angle(unknown.ratio_angle,zf,za)
            rel=float(np.max(abs(unknown.angle_density/published-1)))
            assert rel<1e-9
            known=estimate([f],[a],isotropic_covariance(1,sf,sa),
                           known_direction=[1.],ratio_order=4096)
            expected=katz(f,sf)/katz(a,sa)
            point_rel=abs(known.mass/expected-1)
            assert point_rel<1e-9
            rows.append(dict(z_force=zf,z_acceleration=za,sigma_force=sf,sigma_acceleration=sa,
                             max_density_relative_error=rel,
                             known_direction_point_relative_error=point_rel))
    out=dict(scope='Independent Gaussian scalar channels; full S^0 law and known-positive-direction point.',
             liseo_source='https://www.researchgate.net/publication/5182203_Bayesian_and_conditional_frequentist_analyses_of_the_Fieller%27s_problem_A_critical_review',
             katz_formula_source='https://m.sc.niigata-u.ac.jp/~hirukawa/seminar/niigata2017_program/Chang.pdf',
             cases=rows,
             max_density_relative_error=max(r['max_density_relative_error'] for r in rows),
             max_known_direction_point_relative_error=max(r['known_direction_point_relative_error'] for r in rows))
    target=Path(__file__).resolve().parent/'results/literature_equivalence.json'
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='cases'},indent=2))


if __name__=='__main__': run()
