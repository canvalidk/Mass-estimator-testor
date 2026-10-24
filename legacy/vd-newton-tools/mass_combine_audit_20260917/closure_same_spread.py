"""Scratch witness: same ratio-of-means point/log SD, different mass curves."""
import argparse
import json
import math
import numpy as np
from math import erf

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--order', type=int, default=32)
parser.add_argument('--cutoff', type=float, default=48.)
args = parser.parse_args()
CUTOFF = args.cutoff
gx,gw=np.polynomial.legendre.leggauss(args.order)
def quad(fn,a,b,**unused):
    count=max(1,math.ceil((b-a)/4))
    edges=np.linspace(a,b,count+1)
    result=0.
    for left,right in zip(edges[:-1],edges[1:]):
        mid,half=(left+right)/2,(right-left)/2
        result+=half*sum(w*fn(mid+half*x) for x,w in zip(gx,gw))
    return result,None

def brentq(fn,a,b,xtol=1e-13):
    fa,fb=fn(a),fn(b)
    assert fa*fb<0
    for _ in range(64):
        mid=(a+b)/2
        fm=fn(mid)
        if fm==0 or b-a<xtol:
            return mid
        if fa*fm<0:
            b=mid
        else:
            a,fa=mid,fm
    return (a+b)/2

def raw_log_density(ell, A, B, C):
    sech = 1 / np.cosh(ell)
    h2 = (A+B)/2 + (A-B)/2*np.tanh(ell) + C*sech
    h = np.sqrt(max(0., h2))
    J = 1. if h < 1e-8 else math.sqrt(math.pi/2)*np.exp(h2/2)*erf(h/math.sqrt(2))/h
    return sech*J/2

def law(A, B, C):
    args = (A, B, C)
    z = quad(lambda l: raw_log_density(l,*args), -CUTOFF,CUTOFF,epsabs=2e-11,epsrel=2e-13,limit=200)[0]
    return lambda l: raw_log_density(l,*args)/z

def stats(q):
    mean = quad(lambda l:l*q(l),-CUTOFF,CUTOFF,epsabs=2e-12,epsrel=2e-13)[0]
    var = quad(lambda l:(l-mean)**2*q(l),-CUTOFF,CUTOFF,epsabs=2e-12,epsrel=2e-13)[0]
    def cdf(x):
        return quad(q,-CUTOFF,x,epsabs=2e-12,epsrel=2e-13)[0]
    quantiles = [math.exp(brentq(lambda l:cdf(l)-p,-40,40,xtol=1e-13)) for p in (.025,.5,.975)]
    return {'log_mean':mean,'log_sd':math.sqrt(var),'mass_quantiles_025_50_975':quantiles,'cdf_at_half':cdf(math.log(.5)),'log_density_at_0_1_3':[q(l) for l in (0,1,3)]}

q1=law(4,4,2)
target=quad(lambda l:l*l*q1(l),-CUTOFF,CUTOFF,epsabs=2e-12,epsrel=2e-13)[0]
def discrepancy(c):
    q=law(8,8,c)
    return quad(lambda l:l*l*q(l),-CUTOFF,CUTOFF,epsabs=2e-12,epsrel=2e-13)[0]-target
c2=brentq(discrepancy,0,8,xtol=1e-13)
q2=law(8,8,c2)
q3=law(9,1,3)
def combine(q):
    raw=lambda l:q(l)*q3(l)*math.exp(-l)
    z=quad(raw,-CUTOFF,CUTOFF,epsabs=2e-13,epsrel=2e-13)[0]
    return lambda l:raw(l)/z

result={'parameters1':[4,4,2,1], 'parameters2':[8,8,c2,1],
        'point_both_exact':1,'curve1':stats(q1),'curve2':stats(q2),
        'density_l1_distance':quad(lambda l:abs(q1(l)-q2(l)),-CUTOFF,CUTOFF,epsabs=2e-12,limit=200)[0],
        'third_parameters':[9,1,3,1],'combined1':stats(combine(q1)),'combined2':stats(combine(q2)),
        'cutoff':CUTOFF, 'order':args.order}
print(json.dumps(result,indent=2))
