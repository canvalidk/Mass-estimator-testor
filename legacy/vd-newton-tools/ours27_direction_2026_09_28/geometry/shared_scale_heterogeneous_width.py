"""Shared dt/t scale leaves a regular heterogeneous-covariance width defect.

This is a deterministic asymptotic calculation, not a tuned score or a Monte
Carlo selection. It checks the candidate construction, not a final equation.

REPEATED EXPERIMENT
-------------------
Each independent reading has two component axes. Force SDs are s=(1,10),
acceleration SDs are (1,1), and tested/generated m=1. Use u=log(m). In each
axis's whitened two-channel plane the candidate line has angle
theta_i(u)=atan(exp(u)/s_i), with slope at u=0
  w_i = dtheta_i/du = s_i/(1+s_i^2).
Let the squared candidate length in each whitened line be r_i^2=(0.1,9.9).
Thus one possible physical candidate acceleration is
  a_i = r_i / sqrt(1+1/s_i^2), F_i=a_i.
This merely specifies a Gaussian candidate experiment; it does not assert an
ontic mass behind the operational definition.

At u=0, resolve a noisy observation into parallel and perpendicular channels:
  A_ij=r_i+xi_ij, B_ij=eta_ij,
with independent standard normals. The profiled log likelihood's score is
  S_N=sum_ij w_i A_ij B_ij.
The score variance and expected curvature per reading are
  Var(S_1)=J+G, J=sum r_i^2*w_i^2, G=sum w_i^2.
Indeed ell'' has terms w_i' A_i B_i+w_i^2(B_i^2-A_i^2), whose expectations
are -r_i^2*w_i^2. Under the regular independent-repetition limit,
  Var(u_hat) ~ (J+G)/(N*J^2).
The expected contrast has a unique positive minimum at m=1, so the estimator
is consistent and the usual local expansion is applicable. This is not a
weak-identification or boundary regime: N*J tends to infinity.

SHARED OBSERVABLE-VARIANCE MIXTURE
--------------------------------
The scalar t>=1 candidate law has starting measure dt/t and conditional
candidate covariance (t-1) times the whitened subspace identity. For M=2N,
the marginal evidence is exp(-T/2)*F_M(H(u)), with F_M strictly increasing.
Consequently it preserves all profile likelihood modes exactly, including
their ordering, without choosing or averaging a local branch. A fixed smooth
positive mass prior contributes only O(1) log curvature and cannot change the
leading N-order result here.

Its exact conditional shrinkage is kappa(H)=E[beta|H], beta=1-1/t. Since
  H(u_hat)/N -> r_1^2+r_2^2+2 = 12,
  M/N = 2,
the truncated-gamma scale posterior gives kappa -> 1-M/H = 5/6. The marginal
log evidence curvature at a stationary point is EXACTLY kappa times profile
curvature because H'=0 there. Hence its local law variance is
  V_law ~ 1/(N*kappa*J),
and
  V_sampling/V_law -> kappa*(J+G)/J.

The scale inferred from total standardized energy treats the two component
directions alike, but their ability to change the mass angle differs by w_i.
Signal lies mostly on the less sensitive axis. Therefore the shared scale
does not determine the required local uncertainty in this permitted diagonal
subcase of full covariance. It corrects the equal-angle-slope case, not the
general geometry. Integrating the scale rather than plugging it in does not
remove this leading error because its posterior concentrates as N grows.
"""
from fractions import Fraction
from pathlib import Path
import json
import math


def main():
    w=[Fraction(1,2),Fraction(10,101)]
    r2=[Fraction(1,10),Fraction(99,10)]
    J=sum(rr*ww*ww for rr,ww in zip(r2,w))
    G=sum(ww*ww for ww in w)
    beta=Fraction(5,6)
    required=J/(J+G)
    ratio=beta/required
    output={
        'source_commit':'b4bfc2040031ef010416334d0c5f4aa2ad05cbd8',
        'force_sd':[1,10], 'acceleration_sd':[1,1], 'mass':1,
        'squared_whitened_candidate_lengths':[float(x) for x in r2],
        'log_mass_angle_slopes':[float(x) for x in w],
        'J_per_reading':float(J),'J_exact':str(J),
        'G_per_reading':float(G),'G_exact':str(G),
        'sampling_score_variance_per_reading':float(J+G),
        'shared_scale_beta_limit':float(beta),
        'required_local_beta_for_variance_match':float(required),
        'required_local_beta_exact':str(required),
        'sampling_variance_over_law_variance_limit':float(ratio),
        'variance_ratio_exact':str(ratio),
        'sampling_sd_over_law_sd_limit':math.sqrt(float(ratio)),
        'sampling_variance_times_N_limit':float((J+G)/(J*J)),
        'law_variance_times_N_limit':float(1/(beta*J)),
        'status':'Regular asymptotic under-width for shared scalar scale in heterogeneous geometry.'
    }
    assert 2.6<float(ratio)<2.61
    assert .319<float(required)<.320
    Path(__file__).with_name('shared_scale_heterogeneous_width_results.json').write_text(json.dumps(output,indent=2),encoding='utf-8')
    print(json.dumps(output,indent=2))


if __name__=='__main__':
    main()
