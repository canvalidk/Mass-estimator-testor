"""Scratch audit of candidate-graph Grassmann joint volume; no adopted equation.

Single active scalar reading (F,a)=(2,1), unit noise. Compare d=1 with d=2
having an additional measured zero axis and channel noises (epsilon,epsilon).
In whitened coordinates the added axis and its graph prior do not depend on
epsilon. Therefore the d=2 mass law below is its epsilon->0 limit as well.

For N=1 with equal instrument ratios, theta=atan(m), joint graph volume is
  r/(1+r*r)**((d+2)/2) dtheta d^d v.
After angular integration, likelihood factors independent of theta are omitted.
Quadrature is radial 1D plus the mass angle. This avoids stochastic error.
"""
import json
import argparse
import numpy as np


F,a=2.0,1.0
parser=argparse.ArgumentParser()
parser.add_argument('--angles',type=int,default=1201)
parser.add_argument('--radial',type=int,default=300)
args=parser.parse_args()
theta=np.linspace(0.0,np.pi/2,args.angles)
proj=F*np.sin(theta)+a*np.cos(theta)


nodes,weights=np.polynomial.legendre.leggauss(args.radial)
radial=7*(nodes+1)
weights=7*weights


def evidence(t,d):
    if d==1:
        # exp(-r^2/2)*cosh(r*t), kept bounded on integration interval.
        val=np.exp(-radial*radial/2)*np.cosh(radial*t)*radial/(1+radial*radial)**1.5
    else:
        val=np.exp(-radial*radial/2)*np.i0(radial*t)*radial*radial/(1+radial*radial)**2
    return np.dot(weights,val)


results=[]
densities=[]
for d in (1,2):
    raw=np.array([evidence(t,d) for t in proj])
    Z=np.trapezoid(raw,theta)
    p=raw/Z
    cdf=np.r_[0.0,np.cumsum((p[1:]+p[:-1])*np.diff(theta)/2)]
    qs=np.tan(np.interp([.025,.5,.975],cdf,theta))
    densities.append(p)
    results.append({"dimension":d,"mass_quantiles_025_50_975":qs.tolist(),
                    "density_theta_at_0_pi4_pi2":p[[0,(args.angles-1)//2,args.angles-1]].tolist()})


def mu_matrix(V):
    d,n=V.shape
    S=V@V.T
    return np.sqrt(np.trace(S@np.linalg.inv(np.eye(d)+S)))*np.linalg.det(np.eye(d)+S)**(-(n+d)/2)


matrices=[np.eye(2),np.array([[np.sqrt(2),0.],[0.,0.]])]
flat_mu=np.sqrt(2)/(1+2)**3
print(json.dumps({
    "empty_axis_law_comparison":results,
    "total_variation_between_mass_laws":float(.5*np.trapezoid(np.abs(densities[0]-densities[1]),theta)),
    "axis_flattening_same_squared_norm_2": [{"matrix":V.tolist(),"matrix_graph_mu":float(mu_matrix(V)),"flattened_graph_mu":float(flat_mu)} for V in matrices],
},indent=2))
