"""Independent bounded audit of conditional_checks; scratch only."""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'03_trace_and_evaluator/mass_estimator/investigation'))
import conditional_checks as cc

xg,wg=np.polynomial.legendre.leggauss(192)
ug=(xg+1)/2
wg=wg/2

def closed3(c,s):
    if s<=200:
        j=float(np.dot(wg,np.exp(-s*(1-ug*ug)/2)))
    else:
        v=100*ug
        # Omitted tail <= 2 exp(-50)/s, since v<=s and exponent<=-v/2.
        j=float(np.dot(wg,np.exp(-v+v*v/(2*s))))*100/s
    return math.erf(math.sqrt(c/2))-math.sqrt(2*c/math.pi)*math.exp(-c/2)*j

def main():
    cases=[]
    for c in np.geomspace(1e-8,1e12,41):
        for s in np.r_[0.,np.geomspace(1e-8,1e12,41)]:
            ref=closed3(float(c),float(s))
            got=float(cc.conditional_cdf(c,s,3))
            cases.append((abs(got-ref),float(c),float(s),got,ref))
    ordinary=[q for q in cases if q[1]<=1e5 and q[2]<=1e5]
    extremes=sorted(cases,reverse=True)[:12]
    roots=[]
    for k in (1,3,24,1536,12288):
        for L in (1e3,1e4,1e5,1e6,1e8,1e10,1e12):
            if L<=cc.chi_quantile(k,.95): continue
            b=float(cc.fitted_threshold(L,k))
            c=L-b
            ref=closed3(c,b) if k==3 else None
            roots.append(dict(k=k,L=L,b=b,c=c,closed3=ref,
                              refined=float(cc.conditional_cdf(c,b,k,256))))
    # Independent conditional simulation: no Gram summaries or production LR.
    rng=np.random.default_rng(18375921)
    n=1_000_000
    mc=[]
    for k in (1,3,24,1536,12288):
        z=rng.standard_normal(n)
        U=z*z
        W=rng.chisquare(k-1,n) if k>1 else np.zeros(n)
        for s in (0.,float(k),float(k)*10):
            lo,hi=0.,cc.chi_quantile(k,.95)
            for _ in range(48):
                mid=(lo+hi)/2
                if float(cc.conditional_cdf(mid,s,k))<.95: lo=mid
                else: hi=mid
            c=(lo+hi)/2
            coverage=float(np.mean(U/c+W/(c+s)<=1))
            mc.append(dict(k=k,s=s,c=c,coverage=coverage,
                           standard_errors=(coverage-.95)/math.sqrt(.95*.05/n)))
    report=dict(k3_grid_count=len(cases),
                ordinary_max_error=max(q[0] for q in ordinary),
                extremes=extremes,roots=roots,conditional_mc=mc)
    out=Path(__file__).with_name('review_results.json')
    out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
