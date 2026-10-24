"""Independent NumPy-only derivation check, not production inference."""
import math
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / '03_trace_and_evaluator/mass_estimator/investigation'))
from comparison_posterior import infer_batch

def beta_nodes(a,b,n):
    # Jacobi alpha=b-1, beta=a-1; weights normalized to expectation.
    alpha,beta=b-1,a-1
    k=np.arange(n,dtype=float)
    diag=(beta*beta-alpha*alpha)/((2*k+alpha+beta)*(2*k+alpha+beta+2))
    j=k[1:]
    off=2/(2*j+alpha+beta)*np.sqrt(j*(j+alpha)*(j+beta)*(j+alpha+beta)/((2*j+alpha+beta-1)*(2*j+alpha+beta+1)))
    eig,vec=np.linalg.eigh(np.diag(diag)+np.diag(off,1)+np.diag(off,-1))
    return (eig+1)/2,vec[0]**2

def gamma_nodes(k,n):
    j=np.arange(n,dtype=float)
    off=np.sqrt(j[1:]*(j[1:]+k-1))
    eig,vec=np.linalg.eigh(np.diag(2*j+k)+np.diag(off,1)+np.diag(off,-1))
    return eig,vec[0]**2

def radial(h2):
    h=np.sqrt(np.maximum(h2,0))
    erf=np.fromiter((math.erf(float(x)/math.sqrt(2)) for x in h.flat),float,h.size).reshape(h.shape)
    factor=np.divide(math.sqrt(math.pi/2)*erf,h,out=np.ones_like(h),where=h>1e-12)
    return np.exp(h2/2)*factor

def collapsed(F,A,sf,sa,nuf,nua,z,order=80):
    kf,ka=nuf/2,nua/2
    bf,ba=kf*sf**2,ka*sa**2
    u,w=beta_nodes(kf+1,ka+1,order)
    P=u*np.dot(F,F)/bf
    Q=(1-u)*np.dot(A,A)/ba
    D=np.sqrt(u*(1-u)/(bf*ba))*np.dot(F,A)
    massscale=np.sqrt((1-u)*bf/(u*ba))
    zz=np.asarray(z)[None,:]-np.log(massscale[:,None])
    theta=np.arctan(np.exp(zz))
    sn,cs=np.sin(theta),np.cos(theta)
    H=P[:,None]*sn**2+Q[:,None]*cs**2+2*D[:,None]*sn*cs
    C=1+(P[:,None]+Q[:,None]-H)/2
    v,vw=np.polynomial.legendre.leggauss(128)
    v,vw=(v+1)/2,vw/2
    R=np.zeros_like(H)
    for vv,ww in zip(v,vw):
        R+=ww*np.exp(-(kf+ka+2)*np.log(C+H*vv*vv/2))
    density=kf*ka/(bf*ba)*np.sum(w[:,None]*R/(2*np.cosh(zz)),axis=0)
    tx,tw=np.polynomial.legendre.leggauss(256)
    theta=(tx+1)*math.pi/4
    sn,cs=np.sin(theta),np.cos(theta)
    H=P[:,None]*sn**2+Q[:,None]*cs**2+2*D[:,None]*sn*cs
    C=1+(P[:,None]+Q[:,None]-H)/2
    first=np.exp(-(kf+ka+1.5)*np.log(C))
    nf=np.sum(w*np.sqrt(bf/u)*np.sum(tw*sn*first,axis=1))
    na=np.sum(w*np.sqrt(ba/(1-u))*np.sum(tw*cs*first,axis=1))
    return density,nf/na

def mixture(F,A,sf,sa,nuf,nua,z,order=48):
    kf,ka=nuf/2,nua/2
    bf,ba=kf*sf**2,ka*sa**2
    lf,wf=gamma_nodes(kf,order)
    la,wa=gamma_nodes(ka,order)
    lf,la=np.meshgrid(lf/bf,la/ba,indexing='ij')
    weights=np.outer(wf,wa).ravel()
    lf,la=lf.ravel(),la.ravel()
    P=lf*np.dot(F,F)
    Q=la*np.dot(A,A)
    D=np.sqrt(lf*la)*np.dot(F,A)
    zz=np.asarray(z)[None,:]-.5*np.log(la[:,None]/lf[:,None])
    theta=np.arctan(np.exp(zz))
    sn,cs=np.sin(theta),np.cos(theta)
    H=P[:,None]*sn**2+Q[:,None]*cs**2+2*D[:,None]*sn*cs
    den=weights*lf*la*np.exp(-(P+Q)/2)
    density=np.sum(den[:,None]*radial(H)/(2*np.cosh(zz)),axis=0)
    tx,tw=np.polynomial.legendre.leggauss(256)
    theta=(tx+1)*math.pi/4
    sn,cs=np.sin(theta),np.cos(theta)
    H=P[:,None]*sn**2+Q[:,None]*cs**2+2*D[:,None]*sn*cs
    first=np.exp((H-P[:,None]-Q[:,None])/2)
    nf=np.sum(weights*np.sqrt(lf)*la*np.sum(tw*sn*first,axis=1))
    na=np.sum(weights*lf*np.sqrt(la)*np.sum(tw*cs*first,axis=1))
    return density,nf/na

if __name__=='__main__':
    z=np.array([-4.,-1.,0.,1.,4.])
    for F,A,sf,sa,nuf,nua in [
        (np.array([0.,0.,0.]),np.array([0.,0.,0.]),1.,1.,6,6),
        (np.array([2.,.5,0.]),np.array([1.,-.2,.1]),1.2,.7,6,10),
        (np.array([3.,0.,0.]),np.array([-2.,1.,.3]),.7,1.4,10,6),
    ]:
        c,cp=collapsed(F,A,sf,sa,nuf,nua,z)
        m,mp=mixture(F,A,sf,sa,nuf,nua,z)
        print(dict(nuf=nuf,nua=nua,collapsed_density=c.tolist(),mixture_density=m.tolist(),max_relative_density_difference=float(np.max(np.abs(c-m)/c)),collapsed_point=cp,mixture_point=mp))
    F,A=np.array([2.,.5,0.]),np.array([1.,-.2,.1])
    known=infer_batch(F/1.2,A/.7)['points']['flat_joint'][0]*1.2/.7
    for nu in [100.,1000.,10000.,100000.]:
        _,point=collapsed(F,A,1.2,.7,nu,nu,z,32)
        print(dict(nu=nu,point=point,known=known,relative_difference=point/known-1))
