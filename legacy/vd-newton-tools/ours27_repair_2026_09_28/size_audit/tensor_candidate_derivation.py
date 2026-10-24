"""Constructive tensor-size route, not an adopted estimator or project record.

SCOPE
-----
One group of N independent readings with the same positive-definite 2d x 2d
error covariance Sigma. The argument extends to a fixed number of such groups,
each with N_g growing. It does not solve entirely unreplicated heterogeneous
covariances, a universal starting measure, or finite-N totality.

GEOMETRY
--------
W = Sigma**(-1/2); H(m)=W B_m; K=(H.T H)**(-1);
U=H K**(1/2), P=U U.T. U is an orthonormal frame for the candidate subspace.
In a locally horizontal frame U.T U'=0, put D=(I-P)U'. Any original frame can
be put in this gauge locally; a full freely varying R absorbs its tangent
rotation. D.T D is orthogonally equivalent to K**(1/2) V**(-1) K**(1/2).

WORKING CANDIDATE MEASURE
-------------------------
At candidate m, v_i ~ N(0,R), R a full symmetric positive-definite tensor,
and z_i | v_i,m ~ N(U v_i,I). The marginal is

    C(m,R) = I + U R U.T,
    log L = const - N/2 log|I+R|
            + 1/2 sum t_i.T R(I+R)**(-1) t_i,   t_i=U.T z_i.

The omitted const is -1/2 sum ||z_i||**2 and noise normalizers, independent
of both m and R. The law of candidates given data is

    v_i | m,R,z_i ~ N(Beta t_i, Beta), Beta=R(I+R)**(-1).

Physical acceleration candidates are K**(1/2) v_i in the displayed frame.
Force candidates remain m times acceleration. A ratio-of-mean stacked sizes
can therefore be computed from this joint law exactly as in ours27. No new
point-summary rule is proposed here.

FIXED CANDIDATES, NO ONTIC GAUSSIAN PUSH ASSUMPTION
-------------------------------------------------
For the calibration argument only, let an admissible fixed candidate sequence
be z_i=U_0 v_i + epsilon_i, epsilon_i~N(0,I). Let

    R_N = (1/N) sum v_i v_i.T -> R_0 > 0.

R_N is an uncentered second moment, so a nonzero average candidate is allowed.
At (m_0,R_N), the summed expected log working likelihood depends only on
sum E[z_i z_i.T]=N(I+U_0 R_N U_0.T). Thus its expected Hessian is exactly the
Gaussian covariance Fisher Hessian, despite fixed, arbitrary v_i.

In the horizontal frame write z_parallel=v_i+e_i and z_perp=f_i, with e_i
and f_i independent standard Gaussian. The mass score is

    s_m = sum z_perp_i.T D R(I+R)**(-1) z_parallel_i.

Its mean is zero. At R=R_N, its variance is exactly

    sum tr[D R(I+R)**(-1) (I+v_i v_i.T)
                (I+R)**(-1) R D.T]
      = N tr[D R(I+R)**(-1) R D.T] = I_mass.

The expected negative mass curvature is the same I_mass. Every nuisance
R-score is a centered quadratic function of z_parallel only. Its covariance
with s_m is zero, since each s_m contains the independent centered f_i.
Expected cross curvature also vanishes (orthogonal covariance blocks).
Consequently regular nuisance estimation does not change this first-order
mass variance, and the working posterior curvature is first-order calibrated
even against fixed non-Gaussian candidate sequences.

Regularity still matters: fixed d, enough repetitions, identifiable m,
R_N -> R_0 positive-definite, and a local CLT/no dominating reading condition.
Singular R_0 requires boundary analysis; no theorem about that case is claimed.

WHY THIS EXCEEDS A SCALAR TEMPERATURE
-----------------------------------
For R=omega**2 I,

    I_mass = N omega**4/(1+omega**2) tr(D.T D) = beta J.

For general R, the matrix Beta weighs different tangent directions differently.
Its information need not equal any isotropic mean-size correction. Those
directions reflect physical candidate variation as well as instrument geometry.

At fixed m, profiling R >= 0 is explicit. If S_m=(1/N)sum t_i t_i.T has
eigenvalues lambda_j, then R_hat=(S_m-I)_+, and

    log L_profile = const + N/2 sum_{lambda_j>1}
                                      (lambda_j - 1 - log lambda_j).

That is useful for checking but is NOT a replacement for integrating R when
defining a law over candidates. A nuisance measure is still required.

OPEN MEASURE / TOTALITY ISSUE
-----------------------------
The natural covariance Jeffreys volume |I+R|**(-(d+1)/2)dR on R>0 (noise floor
I) is invariant under orthogonal tangent-frame changes. Combined with N
readings it has an integrable tail only if N>d-1. For one eigenvalue lambda
going to infinity with others bounded, the eigenvalue volume contributes
lambda**(d-1), the prior lambda**(-(d+1)/2), and likelihood lambda**(-N/2),
giving tail lambda**((d-3-N)/2). This is nonintegrable for N<=d-1.
It is therefore not an all-input replacement (e.g. N=1,d=3 fails).
This route also learns from appended nulls, since their likelihood includes
|I+R|**(-1/2); exact appended-null invariance would need changing.

BOUNDED MATRIX FORM
-------------------
Set Beta=R(I+R)**(-1), so 0<Beta<I. The symmetric-matrix Jacobian is
dR=|I-Beta|**(-(d+1)) dBeta. The conditional nuisance Fisher volume times
likelihood becomes

    p(Beta | m,z) proportional
       |I-Beta|**((N-d-1)/2) exp[tr(Beta S_sum(m))/2] dBeta.

This is an exponential tilt of a matrix-Beta law with parameters
alpha=(d+1)/2, b=N/2. Its integrability condition b>(d-1)/2 is N>d-1.
For d=1 it is the scalar tilted Beta(1,N/2), precisely the ordered reference
common-size construction. This gives a useful untuned expression on the
regular replicated domain, not a solution to the missing small-N domain.

Also, the mass prior is not automatically sqrt(G) in this tensor model.
Its efficient information tr[(D.T D) R**2 (I+R)**(-1)] does not separate
between m and R unless D.T D is scalar. A full ordered-reference construction
must specify the compact exhaustion and perform its expectation explicitly.
"""

import json
import numpy as np

# Monte Carlo verification of the cross-block score and information using a
# deliberately non-Gaussian, nonzero-mean, fixed candidate list.
rng = np.random.default_rng(280927)
v = np.array([[3., .2], [1., -1.], [.1, 2.], [-.4, .7], [2., -1.3]])
N,d=v.shape
R=v.T@v/N
D=np.array([[.35,.1],[-.2,.7]])
Beta=R@np.linalg.inv(np.eye(d)+R)
L=D@Beta
predicted=N*np.trace(D@R@np.linalg.inv(np.eye(d)+R)@R@D.T)
trials=400000
e=rng.normal(size=(trials,N,d))
f=rng.normal(size=(trials,N,d))
t=v[None,:,:]+e
score=np.einsum('tni,ij,tnj->t',f,L,t)
print(json.dumps({
    'fixed_candidate_tensor':R.tolist(),
    'predicted_score_variance_and_curvature':float(predicted),
    'simulated_score_mean':float(score.mean()),
    'simulated_score_variance':float(score.var()),
},indent=2))
