"""Scoped boundary and asymptotic calculations; all errors retained."""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'03_trace_and_evaluator/mass_estimator/investigation'))
from crossover_vector_checks import infer
ANTITHETIC=False

@lru_cache(None)
def quadrature(order):
    x,w=np.polynomial.legendre.leggauss(order)
    return (x+1)/2,w/2

def point(F,A,order=96):
    p,q=np.sum(F*F,axis=1),np.sum(A*A,axis=1)
    d=np.sum(F*A,axis=1)
    cut=np.where(d>0,.5*np.arctan2(2*d,q-p),np.where(p>=q,math.pi/2,0.))
    x,w=quadrature(order)
    t=np.concatenate((cut[:,None]*x,cut[:,None]+(math.pi/2-cut[:,None])*x),axis=1)
    weights=np.concatenate((cut[:,None]*w,(math.pi/2-cut[:,None])*w),axis=1)
    sn,cs=np.sin(t),np.cos(t)
    h2=p[:,None]*sn*sn+q[:,None]*cs*cs+2*d[:,None]*sn*cs
    peak=np.where(d>=0,(p+q+np.hypot(p-q,2*d))/2,np.maximum(p,q))
    density=np.exp((h2-peak[:,None])/2)*weights
    return np.sum(sn*density,axis=1)/np.sum(cs*density,axis=1)

def stat(v):
    if ANTITHETIC:
        v=(v[:len(v)//2]+v[len(v)//2:])/2
    return dict(mean=float(np.mean(v)),mcse=float(np.std(v,ddof=1)/math.sqrt(len(v))))

def compare(f,a,eps,order=96):
    F=eps[:,:3]+[f,0,0]
    A=eps[:,3:]+[a,0,0]
    m=np.concatenate([point(F[i:i+2048],A[i:i+2048],order) for i in range(0,len(eps),2048)])
    nf,na=np.linalg.norm(F,axis=1),np.linalg.norm(A,axis=1)
    norm=nf/na
    floor=math.sqrt(2/math.pi)
    floored=np.maximum(nf,floor)/np.maximum(na,floor)
    ours=np.log(m)-math.log(f/a)
    theirs=np.log(norm)-math.log(f/a)
    delta={
        'squared_log':stat(ours**2-theirs**2),
        'abs_log':stat(abs(ours)-abs(theirs)),
        'root':stat(4*np.sinh(ours/4)**2-4*np.sinh(theirs/4)**2),
        'factor_two_failure':stat((abs(ours)>math.log(2)).astype(float)-(abs(theirs)>math.log(2)).astype(float)),
        'squared_log_vs_floor':stat(ours**2-(np.log(floored)-math.log(f/a))**2)}
    refined=point(F[:64],A[:64],order*2)
    refinement=float(np.max(abs(refined/m[:64]-1)))
    assert refinement<1e-8,(f,a,refinement)
    if max(f,a)<=16:
        ref,_=infer(F[:32],A[:32],f/a,order=192)
        assert np.max(abs(ref[:,0]/m[:32]-1))<1e-8
    return dict(f=f,a=a,truth=f/a,delta=delta,
                ours=dict(squared_log=stat(ours**2),bias=float(np.mean(ours)),variance=float(np.var(ours)),success2=float(np.mean(abs(ours)<=math.log(2)))),
                norm=dict(squared_log=stat(theirs**2),bias=float(np.mean(theirs)),variance=float(np.var(theirs)),success2=float(np.mean(abs(theirs)<=math.log(2)))),
                high_signal_gap_prediction=3*(f**-4+a**-4)+2/(f*f*a*a),
                refinement=refinement)

def main():
    global ANTITHETIC
    p=argparse.ArgumentParser()
    p.add_argument('--n',type=int,default=65536)
    p.add_argument('--seed',type=int,default=2026091801)
    p.add_argument('--mode',choices=['coarse','high','custom'],default='coarse')
    p.add_argument('--f',type=float,default=8.)
    p.add_argument('--a',type=float,nargs='*')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--antithetic',action='store_true')
    args=p.parse_args()
    ANTITHETIC=args.antithetic
    eps=np.random.default_rng(args.seed).standard_normal((args.n//2 if ANTITHETIC else args.n,6))
    if ANTITHETIC:
        assert args.n%2==0
        eps=np.concatenate((eps,-eps))
    if args.mode=='coarse':
        pairs=[(8.,a) for a in [.025,.05,.075,.1,.125,.15,.175,.2,.25,.5,1.,1.5,1.75,2.,2.125,2.25,2.5,3.,4.,5.,8.,12.,16.,32.]]
    elif args.mode=='high':
        pairs=[(t,t) for t in (8.,12.,16.,24.,32.,48.,64.)]
        pairs +=[(2*t,t) for t in (8.,16.,32.)]
    else:
        pairs=[(args.f,a) for a in args.a]
    rows=[]
    start=time.perf_counter()
    for f,a in pairs:
        row=compare(f,a,eps)
        rows.append(row)
        print(json.dumps(row),flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(n=args.n,seed=args.seed,antithetic=ANTITHETIC,rows=rows,seconds=time.perf_counter()-start),indent=2)+'\n')
if __name__=='__main__':
    main()
