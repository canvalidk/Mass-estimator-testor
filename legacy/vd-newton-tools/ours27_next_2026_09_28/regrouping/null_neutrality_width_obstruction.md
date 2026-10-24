# Restricted null-neutrality / width obstruction

Scratch derivation, 2026-09-28. This is a conditional theorem about a specified
class of equations, not a ruling that the project's accepted premises conflict.
No managed or production files were changed.

## Assumptions

Use one force/acceleration noise pair for all M scalar axes. Whiten those
channels, write the positive mass angle as theta in (0, pi/2), and let
U_theta v be an aligned candidate, where v is in R^M and U_theta is an
isometric embedding into R^(2M). The estimator's mass law is obtained by
ordinary Gaussian likelihood integration against a positive candidate measure
mu_M, with one fixed mass-angle prior pi(theta):

    p_M(theta | z) proportional to
        pi(theta) exp(-T/2) F_M(H_theta),
    T = ||z||^2,
    H_theta = ||U_theta' z||^2,
    F_M(||t||^2) = integral exp(-||v||^2/2 + t'v) mu_M(dv).

Assume:

1. mu_M is independent of theta and invariant under all O(M) rotations of
   the standardized, stacked candidate vector. This is stronger than only
   simultaneous physical laboratory rotation.
2. The construction is defined for every finite M. F_M(H) is finite and
   strictly positive for every H >= 0. mu_M need not be proper.
3. The same pi is used across M; it is positive throughout the angle interior.
4. Appending a measured 0/0 axis with the same noise pair leaves the whole
   normalized mass law exactly unchanged, for every original data vector.

Candidate-prior projectivity is NOT assumed. Neither is prior properness.
The calibration requirement considered below is a further requirement,
not a consequence of assumptions 1-4 or an established VD ruling.

## Null neutrality forces one generating function

An appended zero axis leaves T and H_theta unchanged. Equality of the
normalized mass laws therefore says F_(M+1)(H_theta)/F_M(H_theta) is constant
in theta for each data vector. Take nonzero force, zero acceleration, with
force norm r: H_theta runs over (0, r^2). These intervals overlap and exhaust
(0, infinity). Continuity at zero fixes the constant as
F_(M+1)(0)/F_M(0). Thus

    F_M(H) / F_M(0) = F(H) for every M and H >= 0,
    F(0) = 1.

Now define the probability measure

    rho_M(dv) = exp(-||v||^2/2) mu_M(dv) / F_M(0).

Its moment generating function is F(||t||^2). Setting trailing coordinates
of t to zero proves that the marginals of rho_(M+k) and rho_M have identical
moment generating functions, hence are identical. The rho_M family is
therefore projectively consistent and fully orthogonally invariant.

The infinite orthogonal-invariance representation theorem gives one mixing
probability nu on [0, infinity) such that

    rho_M = integral N(0, beta I_M) nu(d beta),
    F(H) = integral exp(beta H/2) nu(d beta).

The theorem requires no moments; here finiteness of F already supplies all
positive exponential moments of nu. A primary reference is Diaconis and
Freedman (1987), *A dozen de Finetti-style results in search of a theory*,
Theorem (3), pages 399-400:
https://www.numdam.org/item/AIHPB_1987__23_S2_397_0.pdf

Consequently

    gamma(H) = 2 (log F)'(H) = E_H[beta],
    gamma'(H) = Var_H(beta)/2 >= 0,
    lim_(H->infinity) gamma(H) = ess sup_nu beta = b,

where b may be infinite. The last limit follows by exponential tilting:
every interval below the support supremum eventually loses to any interval
with positive measure strictly above it. Flat mu_M is included: rho_M is
standard Gaussian, nu is a point mass at beta=1, and gamma=1.

## Fixed-candidate sampling width

Let the data really be generated at an interior angle theta_0, with a fixed
candidate vector q_M, known independent standard Gaussian errors, and
R_M^2 = ||q_M||^2. In coordinates aligned with theta_0,

    x = q_M + epsilon,    r = eta,
    epsilon, eta independent N(0, I_M).

The maximum of H_theta is the largest eigenvalue Hhat of

    [ ||x||^2   x'r   ]
    [ x'r      ||r||^2].

Write Delta for this matrix's eigenvalue gap. If R_M^2/M -> lambda in
(0, infinity), then Hhat/M -> 1+lambda and Delta/R_M^2 -> 1 in probability.
The fitted angle satisfies

    thetahat - theta_0 = (x'r)/R_M^2 + o_p(M^(-1/2)),
    Var(x'r) = R_M^2 + M,
    sqrt(M) (thetahat-theta_0)
        -> N(0, (lambda+1)/lambda^2).

This result uses fixed candidates, not Gaussian-random candidates. O(M)
symmetry of the errors allows q_M to be rotated to (R_M,0,...,0); the CLT
then follows from the leading Gaussian term and the remaining independent
noise products. At the peak H'' = -2 Delta, so the negative curvature of
log F(H_theta) is exactly gamma(Hhat) Delta.

For finite b>0, the limit gamma(H)->b gives the regular local Gaussian law
directly: on an angle displacement h/sqrt(M),

    H_theta = Hhat - lambda h^2 + o_p(1),
    log F(H_theta) - log F(Hhat) -> -b lambda h^2/2.

A fixed smooth positive pi contributes o_p(1) locally. Thus its first-order
posterior angular variance is 1/(b R_M^2). Correct fixed-sequence sampling
width instead is (R_M^2+M)/R_M^4. Matching requires

    b = lambda/(1+lambda).

One common b cannot satisfy this at two distinct positive lambda values.
For example lambda=1 requires b=1/2, whereas lambda=3 requires b=3/4.
If b=0, nu is concentrated at zero and the mass law gets no data information.
If b=infinity, every fixed nonzero M^(-1/2) displacement loses unbounded log
mass, so the law is sharper than the required regular sampling scale.

Therefore no equation in assumptions 1-4 has universally calibrated regular
many-axis mass-law widths at all positive fixed per-axis signal strengths.
This includes improper candidate measures with finite Gaussian evidence.
It says nothing by itself about ratio-of-mean-sizes point readout calibration.

## Scope of the managed sources

All managed sources below were read with the installed GitHub connector at
canvalidk/VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.

* `Newton-analysis/_dscn_mass_estimator/ours27/ours27_candidate_size_and_null_readings_2026-09-28.md`,
  section 5, lines 144-146, says a known-empty axis does not change the answer;
  a measured-empty axis does not change the law, and its point effect vanishes
  as measurement becomes exact. The text calls this the property as it holds.
  Its clean version is expressly restricted to one noise pair (lines 160-161).
  The exact mass-law neutrality premise above follows the described property,
  but requiring every future equation to preserve it is a substantive choice.
* The same source, section 7, lines 223-226, explicitly leaves undecided the
  trade between calibrated widths and the weak-signal point, and calls the
  planted-mass yardstick unsettled under mass-is-its-estimator. Therefore the
  theorem MUST NOT be presented as a contradiction of already accepted aims.
* `Newton-analysis/_dscn_mass_estimator/ours25b_the_many_reading_readout_2026-09-25.md`,
  sections 5.3-5.4, lines 265-282, explicitly derives arbitrary O(M) stacked
  rotation invariance from the flat law and the statistics P,Q,D. Its setting
  uses one noise pair throughout. This is a documented property, not an
  independent theorem that every instrument-only candidate measure is O(M).
* `Newton-analysis/_dscn_mass_estimator/ours27/ours27_full_covariance_2026-09-27.md`,
  section 4, lines 224-230, says the stated axis-splitting premise applies only
  when covariance is diagonal in the chosen axes, and L8 owes restatement.
  Physical stacked candidate size survives; arbitrary standardized rotations
  under heterogeneous/correlated instruments are not licensed by that remark.
* `Newton-analysis/_dscn_mass_estimator/mass_estimator_lemmas_2026-09-27.md`,
  L8-L13, lines 58-63: axis splitting, isotropy, consistency, one joining rule,
  and Gaussian noise are premises; L11 is point consistency, not calibrated
  interval width. T18 reports the width defect; it does not outlaw it.
* `Newton-analysis/_dscn_mass_estimator/the_combination_equation_is_exact_2026-09-26.md`,
  sections 1-2, records exact sequential/batch equivalence by adding data state,
  with one prior applied once. That operation is distinct from prior
  projectivity. The new proof derives the latter only for rho_M, using the
  additional exact whole-law null-neutrality assumption.

## What can escape the theorem

An estimator can fall outside the result by changing any load-bearing
assumption: null axes may affect the mass law; the candidate measure may
depend on mass/covariance geometry rather than be one radial noise-unit law;
candidate covariance can retain declared experiment/repetition structure;
the mass uncertainty law can be defined by a different inferential operation;
or the universal frequentist width yardstick can be rejected or reformulated.
For example a direct score using both fitted and residual energy need not
be an ordinary radial integrated likelihood F(H). The theorem neither chooses
among those moves nor supplies the missing physical premise that would.
