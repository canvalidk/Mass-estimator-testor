"""Independent Cartesian audit of correlated scalar posterior integrals."""
from pathlib import Path
import json
import math
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'03_trace_and_evaluator/mass_estimator/investigation'))
import weak_acceleration_checks as w

def Phi(x):
    x=np.asarray(x)
    return np.array([.5*math.erfc(-float(v)/math.sqrt(2)) for v in x.ravel()]).reshape(x.shape)

def phi(x):
    return np.exp(-x*x/2)/math.sqrt(2*math.pi)

def independent_cartesian(x,y,m,c,order=1024):
    t,weights=np.polynomial.legendre.leggauss(order)
    upper=max(0.,y)+16
    sd=math.sqrt(1-c*c)
    crossing=(x-c*y)/(m-c) if m != c else -1.
    cuts=sorted(set([0.,upper]+[v for v in (crossing-8*sd/abs(m-c),crossing,crossing+8*sd/abs(m-c)) if 0<v<upper])) if m != c else [0.,upper]
    a=np.concatenate([lo+(t+1)*(hi-lo)/2 for lo,hi in zip(cuts[:-1],cuts[1:])])
    weights=np.concatenate([weights*(hi-lo)/2 for lo,hi in zip(cuts[:-1],cuts[1:])])
    marginal=phi(a-y)
    conditional_mean=x+c*(a-y)
    positivity=Phi(conditional_mean/sd)
    normalizer=weights@(marginal*positivity)
    ratio_probability=weights@(marginal*(Phi((m*a-conditional_mean)/sd)-Phi(-conditional_mean/sd)))/normalizer
    first_force=weights@(marginal*(conditional_mean*positivity+sd*phi(conditional_mean/sd)))
    first_acceleration=weights@(marginal*a*positivity)
    return first_force/first_acceleration,ratio_probability

rows=[]
for c in (-.98,-.8,0.,.8,.98):
    for x,y,m in ((0.,0.,1.),(4.,0.,4.),(3.,-.5,6.),(4.,1.,4.),(-2.,2.,.3),(8.,.05,160.)):
        p,v=w.scalar_batch([x],[y],m,c,192)
        pp,vv=independent_cartesian(x,y,m,c)
        rows.append({'c':c,'x':x,'y':y,'m':m,'point_relative_difference':float(abs(p[0]/pp-1)),
                     'cdf_absolute_difference':float(abs(v[0]-vv))})
print(sorted(rows,key=lambda r:r['cdf_absolute_difference'],reverse=True)[:3])
assert max(r['point_relative_difference'] for r in rows)<2e-8
assert max(r['cdf_absolute_difference'] for r in rows)<2e-9

bins=[]
for k in (.1,1.,2.):
    b=w.bisect(lambda a:w.censor_cdf(a,k),.95,0.,k+12.)
    at=w.coarsened_coverage(b,k)
    right=w.coarsened_coverage(b+1e-8,k)
    assert abs(at-float(Phi(b+k)))<1e-12
    assert abs(right-float(Phi(b-k)))<4e-9
    bins.append({'k':k,'B':b,'coverage_at_B':at,'coverage_right_of_B':right})

result={'cartesian_cases':len(rows),'max_point_relative_difference':max(r['point_relative_difference'] for r in rows),
        'max_cdf_absolute_difference':max(r['cdf_absolute_difference'] for r in rows),'bins':bins}
Path(__file__).with_name('audit_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
