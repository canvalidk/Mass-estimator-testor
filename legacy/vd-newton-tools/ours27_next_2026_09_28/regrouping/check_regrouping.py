"""Scratch counterexamples for the full candidate-covariance proposal.

Not a project record or proposed production estimator. Analysis uses the managed
snapshot b4bfc2040031ef010416334d0c5f4aa2ad05cbd8 and the scratch tensor derivation.

1. Reinterpreting axes as repeated scalar draws changes the covariance family.
2. Keeping one unknown candidate fixed across repetitions is different from
   separate candidates sharing only mass. Instruments/data cannot identify that
   experimental relation by themselves.
3. Under q_i|Omega iid N(0,Omega), finite fourth moments, and arbitrary component
   regrouping, E[q11^2 q12^2]-E[q11^2 q22^2]=2 E[Omega12^2]. Regrouping invariance
   forces Omega12=0 almost surely. Requiring the same after every physical
   rotation forces Omega=lambda I almost surely. Tensor anisotropy requires
   retaining replication/coupling metadata through representation changes.
"""

import json
import math


def gain(x):
    return x - 1.0 - math.log(x) if x > 1.0 else 0.0


# One 2D reading: F=(3,0), a=(0,0), both error covariance blocks I.
# U_m^T z=(3*m/sqrt(1+m*m),0). Compare full 2x2 covariance (one row)
# versus common scalar covariance after calling the axes two separate rows.
profiles = []
for m in [0.1, 1/math.sqrt(8), 0.5, 1.0, 2.0, 10.0]:
    energy = 9*m*m/(1+m*m)
    profiles.append({
        "m": m,
        "projected_energy": energy,
        "one_2d_log_profile": 0.5*gain(energy),
        "two_1d_log_profile": gain(energy/2),
    })

# Equal instruments, two readings (2,1) and (-2,-1). Shared stable candidate
# averages to (0,0): its residual log factor is constant. Separate candidates
# sharing m yield likelihood exp[-(2-m)^2/(1+m^2)].
replication = []
for m in [0.5, 1.0, 2.0, 4.0]:
    replication.append({
        "m": m,
        "stable_pair_mass_factor": 1.0,
        "varying_pairs_mass_factor": math.exp(-(2-m)**2/(1+m*m)),
    })

# Uniformly oriented 2D covariance with eigenvalues 4 and 1/4. Rotation
# averaging is exact: E[Omega12^2]=(lambda1-lambda2)^2/8.
l1, l2 = 4.0, 0.25
offdiag_second = (l1-l2)**2/8
cross_rows = l1*l2 + offdiag_second
same_row = cross_rows + 2*offdiag_second

print(json.dumps({
    "regrouping_profile_counterexample": profiles,
    "stable_vs_varying_candidates": replication,
    "rotation_invariant_tensor_prior_moments": {
        "E_Omega12_squared": offdiag_second,
        "E_q11_squared_q12_squared_same_row": same_row,
        "E_q11_squared_q22_squared_cross_rows": cross_rows,
        "difference": same_row-cross_rows,
    },
}, indent=2))
