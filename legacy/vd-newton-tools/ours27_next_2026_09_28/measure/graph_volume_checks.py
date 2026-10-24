"""Independent pullback-metric and known-empty limit checks.

Whiten candidate Z=U_m V (2d x N), U_m orthonormal, then graph [I_N;Z].
The canonical real Grassmann metric is
 tr[(I+ZZ')^-1 dZ (I+Z'Z)^-1 dZ'].
Its pullback volume in (m,V) coordinates is
 det(I+VV')^(-(N+d)/2) sqrt(tr[T VV'(I+VV')^-1]),
 T=U_m'^T(I-U_m U_m')U_m'. The symmetric K^1/2 frame need not be horizontal;
the Schur complement removes its tangent rotation automatically.

This geometric formula is a mathematical candidate, NOT an adopted equation.
The test also demonstrates its known-empty discontinuity: a second measured
zero axis with noise scaled epsilon is unchanged in standardized coordinates,
so its mass-law effect is independent of epsilon even as epsilon -> 0.
"""

import json
import numpy as np

def symroot(A,power):
    lam,Q=np.linalg.eigh(A)
    return (Q*lam**power)@Q.T

def geometry(Sigma,m):
    d=Sigma.shape[0]//2
    W=symroot(Sigma,-.5)
    def u(x):
        H=W@np.vstack((x*np.eye(d),np.eye(d)))
        K=np.linalg.inv(H.T@H)
        return H@symroot(K,.5)
    U=u(m)
    step=1e-5*max(m,1)
    Up=(u(m+step)-u(m-step))/(2*step)
    T=Up.T@(np.eye(2*d)-U@U.T)@Up
    return U,Up,T

def density(Sigma,m,V):
    d,N=V.shape
    U,Up,T=geometry(Sigma,m)
    S=V@V.T
    B=S@np.linalg.inv(np.eye(d)+S)
    return float(np.sqrt(np.trace(T@B))*np.linalg.det(np.eye(d)+S)**(-(N+d)/2))

def metric_test():
    rng=np.random.default_rng(92827)
    d,N,m=2,3,1.7
    C=rng.normal(size=(2*d,2*d))
    Sigma=C@C.T+np.diag([.2,.3,.4,.5])
    V=rng.normal(size=(d,N))
    U,Up,T=geometry(Sigma,m)
    Z=U@V
    AL=np.linalg.inv(np.eye(2*d)+Z@Z.T)
    AR=np.linalg.inv(np.eye(N)+Z.T@Z)
    dz=[Up@V]
    for j in range(d*N):
        E=np.zeros((d,N));E.flat[j]=1
        dz.append(U@E)
    metric=np.array([[np.trace(AL@P@AR@Q.T) for Q in dz] for P in dz])
    direct=float(np.sqrt(np.linalg.det(metric)))
    closed=density(Sigma,m,V)
    schur=metric[0,0]-metric[0,1:]@np.linalg.solve(metric[1:,1:],metric[1:,0])
    phi=.63
    O=np.array([[np.cos(phi),-np.sin(phi)],[np.sin(phi),np.cos(phi)]])
    O2=np.zeros((4,4));O2[:2,:2]=O;O2[2:,2:]=O
    rot=density(O2@Sigma@O2.T,m,O@V)
    S=np.diag([3.,3.,1.,1.])
    units=3*density(S@Sigma@S.T,3*m,V)
    P=np.block([[np.zeros((d,d)),np.eye(d)],[np.eye(d),np.zeros((d,d))]])
    swap=density(P@Sigma@P.T,1/m,V)/(m*m)
    return {'direct_volume':direct,'closed_volume':closed,
      'relative_volume_error':direct/closed-1,
      'Schur_mass_metric':float(schur),
      'closed_mass_metric':float(np.trace(T@(V@V.T)@np.linalg.inv(np.eye(d)+V@V.T))),
      'nonhorizontal_frame_cross_norm':float(np.linalg.norm(metric[0,1:])),
      'rotation_relative_error':rot/closed-1,
      'force_units_with_jacobian_relative_error':units/closed-1,
      'swap_with_jacobian_relative_error':swap/closed-1}

def null_axis_test(n=240):
    x,w=np.polynomial.legendre.leggauss(n)
    r=(x+1)*6
    wr=w*6
    th=(x+1)*np.pi/4
    wt=w*np.pi/4
    h=2*np.sin(th)
    base1=2*r*(1+r*r)**-1.5*np.exp(-r*r/2)*wr
    base2=2*np.pi*r*r*(1+r*r)**-2*np.exp(-r*r/2)*wr
    L1=base1@np.cosh(r[:,None]*h[None,:])
    L2=base2@np.i0(r[:,None]*h[None,:])
    odds1=float(base1@np.cosh(2*r)/np.sum(base1))
    odds2=float(base2@np.i0(2*r)/np.sum(base2))
    prob_light1=float(np.sum(wt[th<np.pi/4]*L1[th<np.pi/4])/np.sum(wt*L1))
    prob_light2=float(np.sum(wt[th<np.pi/4]*L2[th<np.pi/4])/np.sum(wt*L2))
    return {'reading_force':[2,0],'reading_acceleration':[0,0],
      'one_axis_endpoint_law_odds':odds1,
      'with_null_axis_endpoint_law_odds_for_every_epsilon_positive':odds2,
      'one_axis_P_mass_below1':prob_light1,
      'with_null_axis_P_mass_below1_for_every_epsilon_positive':prob_light2,
      'statement':'Second axis noise SD=(epsilon,epsilon); epsilon cancels from law exactly.'}

def null_axis_readout(n=120):
    x,w=np.polynomial.legendre.leggauss(n)
    r=(x+1)*6;wr=w*6
    th=(x+1)*np.pi/4;wt=w*np.pi/4
    phi=(x+1)*np.pi/4;wp=w*np.pi/4
    h=2*np.sin(th)
    A1=(2*r*r*(1+r*r)**-1.5*np.exp(-r*r/2)*wr)@np.cosh(r[:,None]*h[None,:])
    point1=float((wt*np.sin(th))@A1/((wt*np.cos(th))@A1))
    cosh=np.cosh(r[:,None,None]*h[None,:,None]*np.cos(phi)[None,None,:])
    base=4*r**3*(1+r*r)**-2*np.exp(-r*r/2)*wr
    points={}
    for eps in [1.,.1,.001,0.]:
        angular=cosh@(wp*np.sqrt(np.cos(phi)**2+eps*eps*np.sin(phi)**2))
        A2=base@angular
        points[str(eps)]=float((wt*np.sin(th))@A2/((wt*np.cos(th))@A2))
    return {'one_axis_point':point1,'with_null_axis_point_by_noise_epsilon':points,
      'epsilon0':'Limit of positive-noise augmented model; not the model with an axis removed.'}

if __name__=='__main__':
    print(json.dumps({'metric':metric_test(),'known_empty':null_axis_test(),
      'known_empty_readout':null_axis_readout()},indent=2))
