"""Known-empty continuity test for noise-metric coarea candidate law.

Base: one scalar reading F=2,a=0,sigma_F=sigma_a=1. Under ambient Gaussian
conditioning on matching signs, mhat=[2(Phi(2)-1/2)+phi(2)]/phi(0).

Append one observed-zero axis with sigma_F=sigma_a=epsilon>0. Force and
acceleration candidates on this axis have physical sizes proportional to
epsilon. In standardized candidate coordinates, however, the coarea
Jacobian remains sqrt(v1^2+v2^2), independent of epsilon. At angle theta,
v1~N(2 sin(theta),1),v2~N(0,1), and Q=4 cos(theta)^2.

The physical stacked acceleration size is
cos(theta)*sqrt(v1^2+epsilon^2 v2^2). Therefore the ratio-of-means point uses
 W_epsilon(theta)=E[sqrt(v1^2+epsilon^2 v2^2)*sqrt(v1^2+v2^2)].
As epsilon->0, this tends to E[|v1|sqrt(v1^2+v2^2)], whereas the base
one-coordinate law uses E[v1^2]=1+4 sin(theta)^2.

The resulting point DOES NOT converge to the known-empty base point.
The angle law likewise uses E sqrt(v1^2+v2^2) for every epsilon>0, so is
independent of epsilon and differs from the one-coordinate law's E|v1|.
Hence noise-metric coarea/GFI fails a known-empty-axis continuity property,
despite being proper and invariant under regrouping a fixed observation
space. This is a new rejection test, not a proposed fix or production edit.

We integrate folded v1,v2 on [0,infinity)^2 by Gauss-Legendre quadrature;
this avoids a moving absolute-value kink in Gauss-Hermite coordinates.
Truncations at 14 and 12 are >10 noise SDs above relevant means.
"""
from pathlib import Path
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss


def phi(x):
    return np.exp(-np.asarray(x)**2/2)/np.sqrt(2*np.pi)


def grid(upper,n):
    z,w=leggauss(n)
    return upper*(z+1)/2,w*upper/2


def run(nodes):
    theta,wt=grid(np.pi/2,nodes)
    sn,cs=np.sin(theta),np.cos(theta)
    mean=2*sn
    t,w1=grid(14.,nodes)
    u,w2=grid(12.,nodes)
    folded1=phi(t[None,:]-mean[:,None])+phi(t[None,:]+mean[:,None])
    half2=2*phi(u)*w2
    norm=np.sqrt(t[:,None]**2+u[None,:]**2)
    likelihood=np.exp(-2*cs*cs)
    def point(weight):
        return float(np.sum(wt*likelihood*sn*weight)/np.sum(wt*likelihood*cs*weight))
    base=point(1+mean*mean)
    base_mean_abs=np.sqrt(2/np.pi)*np.exp(-mean*mean/2)+mean*np.array([math.erf(a/np.sqrt(2)) for a in mean])
    augmented_mean_norm=folded1@(w1*(norm@half2))
    base_law=wt*likelihood*base_mean_abs;base_law/=base_law.sum()
    augmented_law=wt*likelihood*augmented_mean_norm;augmented_law/=augmented_law.sum()
    rows=[]
    for epsilon in (1.,.3,.1,.03,.01,.001,.00001,0.):
        weight_t=(np.sqrt(t[:,None]**2+epsilon*epsilon*u[None,:]**2)*norm)@half2
        weight=folded1@(w1*weight_t)
        value=point(weight)
        rows.append({"epsilon":epsilon,"coarea_augmented_point":value,
                     "relative_difference_from_base":value/base-1})
    return {"nodes_per_integral":nodes,"base_point_quadrature":base,
            "augmented_points":rows,
            "base_angle_law_probability_theta_above_pi_over_4":float(base_law[theta>np.pi/4].sum()),
            "augmented_angle_law_probability_theta_above_pi_over_4":float(augmented_law[theta>np.pi/4].sum()),
            "angle_law_total_variation":float(.5*np.sum(abs(base_law-augmented_law)))}


def main():
    base_exact=(2*((1+math.erf(np.sqrt(2)))/2-.5)+float(phi(2)))/float(phi(0))
    coarse=run(160);fine=run(320)
    error=max(abs(a["coarea_augmented_point"]-b["coarea_augmented_point"])
              for a,b in zip(coarse["augmented_points"],fine["augmented_points"]))
    assert error<2e-7
    assert abs(fine["base_point_quadrature"]-base_exact)<1e-12
    assert abs(fine["augmented_points"][-1]["relative_difference_from_base"])>.001
    out={"source_commit":"b4bfc2040031ef010416334d0c5f4aa2ad05cbd8",
         "base_point_closed_form":base_exact,"coarse_to_fine_point_max_difference":error,
         "fine":fine,"status":"Coarea bridge rejected: fails continuity as a measured-zero axis becomes known-empty."}
    Path(__file__).with_name("coarea_known_empty_results.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
