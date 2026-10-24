"""Check the high-SNR fourth-order risk coefficient by exact Gaussian moments.

Both true channel SNRs tend to infinity with fixed positive ratio. Put
x=1/f, y=1/a, k=2 transverse dimensions, z=log(||F||/||A||), and
T=tanh(log(f/a))=(y^2-x^2)/(x^2+y^2). For positive observed dot product,
log(m_TLS)=asinh(sinh(z)/cos(gamma)). On resolved interior data the flat
ratio-of-means differs from this by exponentially small boundary terms.

With independent longitudinal normals u,v and transverse U,V, the first
two orders in norm-ratio log error are e1,e2 below. The first two orders
of 1-cos(gamma) are L2,L3. Expand the risk difference through order four:
  2 E[e1*(T*L3+(1-T^2)*e1*L2)] + 2 E[e2*T*L2] + E[T^2*L2^2].
It simplifies to 3(x^4+y^4)+2*x^2*y^2 in 3D.

This does NOT assign a finite population risk to the boundary-valued TLS
rule. It uses TLS only on the high-probability interior event. For the
original flat estimator, 0<=||F sin(theta)+A cos(theta)||^2<=p+q gives
|log(m_flat)|<=(p+q)/2. Its squared-log loss thus has a polynomial envelope;
Gaussian probability outside the resolved event is exponentially small.
This permits the local expansion to describe its actual expected risk.
"""
import itertools
import json
import numpy as np

# Three-point Gaussian quadrature integrates each coordinate's polynomials
# through degree five, so all fourth-order moments below are exact.
v=np.array([-np.sqrt(3.),0.,np.sqrt(3.)])
w=np.array([1/6,2/3,1/6])
indices=np.array(list(itertools.product(range(3),repeat=6)))
data=v[indices]
weight=np.prod(w[indices],axis=1)
u,v=data[:,0],data[:,3]
U,V=data[:,1:3],data[:,4:6]
U2,V2=np.sum(U*U,axis=1),np.sum(V*V,axis=1)
UV=np.sum(U*V,axis=1)
rows=[]
for x,y in [(1.,1.),(.5,1.),(1.,.5),(.25,1.),(2.,3.)]:
    T=(y*y-x*x)/(x*x+y*y)
    e1=x*u-y*v
    e2=.5*x*x*(U2-u*u)-.5*y*y*(V2-v*v)
    e3=x**3*(u**3/3-u*U2)-y**3*(v**3/3-v*V2)
    L2=.5*(x*x*U2+y*y*V2-2*x*y*UV)
    L3=-x**3*u*U2-y**3*v*V2+x*y*(x*u+y*v)*UV
    delta4=2*e1*(T*L3+(1-T*T)*e1*L2)+2*e2*T*L2+T*T*L2*L2
    actual=float(weight@delta4)
    formula=3*(x**4+y**4)+2*x*x*y*y
    norm4=float(weight@(e2*e2+2*e1*e3))
    expected_norm4=-.25*(x**4+y**4)-.5*x*x*y*y
    assert abs(actual-formula)<1e-10,(actual,formula)
    assert abs(norm4-expected_norm4)<1e-10,(norm4,expected_norm4)
    rows.append(dict(x=x,y=y,exact_moment_gap=actual,formula=formula,norm_fourth_order=norm4))
print(json.dumps(dict(assertions='passed',rows=rows),indent=2))
