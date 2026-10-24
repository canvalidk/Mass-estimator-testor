# Internal continuation state — candidate-law bridge

Scratch only, 2026-09-28. No next equation has been adopted. Sources remain
pinned to VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8. This continues
../ours27_next_2026_09_28/WORKING_STATE.md. No production or curated edits.

## Concrete positive diagnostic

For a common noise pair, let h>=q be the eigenvalues of the standardized
2x2 data Gram matrix and Delta=h-q. At the unrestricted projective TLS
direction the exact conditional score has derivative -Delta/sqrt(h).
Its local information is Delta^2/h, whereas ours27's profile curvature is
Delta. Within scalar powers of that profile this uniquely selects

    beta_obs=Delta/h,
    p(theta) proportional pi(theta) exp[-beta_obs Q(theta)/2].

All-zero data use beta=0. Same-coordinate rotation, units, channel swap,
regrouping and literal measured-null neutrality pass. Regular asymptotic
calibration holds. Twelve numerical cells (4,000 draws each) are in
local_information_probe.py and local_information_results.json. Example:
M=63,R^2=21,m=s, nominal95% coverage .7005 -> .94575; M400,R^2=400,m=s,
.83575 -> .95625. These are diagnostics, not tuning or a uniqueness proof.

The exact Gram construction extends to separable covariance
C_channel tensor S_component by separate channel/component whitening. It
does not automatically extend to arbitrary differing instrument ratios.

## Why this is not an adopted replacement

1. Local curvature fixes only a second derivative. Many global laws share
   it, the same symmetries and the same regular limit. For example with
   I=Delta^2/h and u=sin^2(theta-thetahat), both logL=-Iu/2 and
   logL=-I(u+u^2)/2 work locally but differ on weak data. See reference/
   local_information_audit.py. Root also checked a quartic-in-profile-loss
   counterfamily in local_information_probe.py.
2. Keeping ours27's conditional Gaussian candidates under the powered law
   gives a mathematical joint distribution, but no fixed candidate measure
   plus the declared Gaussian likelihood can produce it. The required
   measure depends on the observed data. A separate isotropic marginal
   odds argument also excludes it. See composition/data_power_bayes_fork.md.
   Do not misstate this as violating before-data pi alone, or claim that
   mass-marginal batch products must work for every hierarchical model.
3. Arbitrary-covariance beta=j/Vc at a global profile maximum fails. Exact
   two-instrument example: s1=sqrt(2)+1,s2=sqrt(2)-1, s_a=1, F_i=k*s_i,a_i=k.
   In u=log m, ell has a unique global maximum at0, zero quadratic curvature,
   and fourth derivative -3k^2/sqrt(2). beta=0 would erase a quartic peak
   which narrows as k^-1/2. Root independently proved the stable global-gap
   identity and verified it in full_covariance_quartic_counterexample.py.
   exact/profile_power_ties.py additionally enumerates stationary roots for
   two globally tied modes having different betas, proving a jump under
   arbitrarily small perturbations. beta can also exceed one.

## Two principled joint-law constructions investigated

- Noise-metric coarea / l2 generalized fiducial geometry gives Jacobian
  sqrt(det(B'Sigma^-1 B))*sqrt(a' H a), with H the mass tangent after removing
  the candidate-vector tangent. Its mass law is exp(-Q/2) times a conditional
  Gaussian expectation of sqrt(a'H a). This is a specified joint rule, but
  retains ours27's leading width defect and fails known-empty continuity:
  F2,a0 base point2.52791131 vs zero-axis limit2.47806104. See exact/
  bridge_candidates.py and coarea_known_empty.py.
- Reference analysis of the invariant rank-one experiment gives uniform
  angle times sqrt(I_RR(R))dR. This uses statistical distinguishability to
  select a radial measure after choosing the reduced experiment. It obtains
  correct regular information R^4/(R^2+M), but fails known-empty continuity:
  point2.32379859 ->1.91213162. Reference-before-reduction gives a different
  radial measure, so the order/reduction is a substantive premise. See
  reference/radial_reference.py.

## Exact global conditional benchmark

In whitened data z, let P_m project onto the candidate subspace. Condition
on P_m z and replace the orthogonal component by independent standard
Gaussian noise: z*=P_m z+W_m xi. Calibrate the WHOLE profile statistic
Lambda_m(z)=sup_t z'P_t z-z'P_m z against that conditional experiment.
This avoids fitted-mode selection and second-order-only approximations.
Root implementation full_covariance_conditional_lr.py uses a finite grid
including the tested mass. Its conditional calibration is valid for that
grid statistic; Monte Carlo approximates the critical values. Three cells
with independent20,000-draw calibration/assessment gave .94725-.95045
coverage at nominal.95. The quartic example's informative contour survives.

But exact conditional LR ALSO fails known-empty continuity: F2,a0,testedm1,
unitnoise gives p=.1572992071; adding a zero axis with bothSD epsilon gives
p=.2391605842 for EVERY epsilon>0. Exact integration and algebra checks are
in exact/conditional_lr_known_empty.py. This remains a benchmark only.

## Composition bridge

The composition source explicitly assumes one propagated joint candidate
law; it does not equate that law with a fresh fit after information is lost.
Sum/contrast factorization of two force measurements with common acceleration
shows that the contrast generally informs total mass. The correct coarsening
identity is the posterior tower identity, with acceleration-weighted averaging
for the ratio-of-means point. A law-level composition rule must push measures
forward: df1 df2 -> f df under f=f1+f2, not reset df. This constrains coherent
families but leaves many possibilities. See composition/composition_bridge.md.

## Still unresolved

A full, unique finite-data joint candidate law satisfying the project's
operational requirements. Neither benchmark calibration, local information,
canonical surface volume, reference priors, nor shared-law composition alone
has selected it. The meaning of reported probability and the operation for
updating/conditioning it remain load-bearing. Do not rename a good simulation
correction as the sought definition or infer a universal impossibility.

## User's Fieller correction: coverage is not an everywhere requirement

The user correctly challenged any implicit demand that finite probability
intervals maintain coverage at arbitrarily weak signal. Fieller confidence
sets can retain coverage by being unbounded; the impossibility concerns
always-bounded sets. Primary sources checked: von Luxburg and Franz,
https://arxiv.org/pdf/0711.0198, and Dufour,
https://jeanmariedufour.github.io/Dufour_1994_ImpossibilityTheorems_W.pdf.

For our Gaussian candidate experiment, as total standardized signal R tends
to zero, the data law tends in total variation to the same noise law P0 for
every generating candidate ratio. A proper law on positive finite mass has
finite interior upper quantile U(Z), so P0(U >= m) tends to zero as m tends
to infinity. Taking positive R sufficiently small for each m proves infimal
coverage zero for its finite central interval. Likewise a positive lower
quantile gives failure toward zero mass. This need not be undercoverage at
every weak-signal mass: default-centre overcoverage is possible.

This is a property of specified candidate experiments, not a premise that
an observer possesses a hidden exact mass. Also distinguish observed exact
zero data from a zero or almost-zero candidate measured with random noise.

weak_identification_check.py uses 10,000 draws per cell, M=63, R^2=.01,
the same uniform-angle instrument prior, and unchanged fixed noise. At
candidate m/s=.001,.01,1,100,1000, finite95% interval coverage for the
local-information diagnostic is respectively 0%,9%,92.49%,8.87%,0%.
Ours27 is .02%,25.08%,56.68%,24.84%,.02%. Refining quadrature from96 to160
nodes changed tested CDFs by less than2e-14. Results are a boundary check,
NOT a ranking or a target for further tuning.

The earlier regime R^2=lambda*M with fixed lambda>0 has regular information
R^4/(R^2+M)=lambda^2*M/(1+lambda) tending to infinity. Its persistent width
mismatch is separate from this unavoidable weak-identification obstruction.
For regular interior comparisons also keep scaled angular distance from the
positive-mass boundaries large. Repeated-experiment calibration remains an
unratified definitional criterion; do not silently promote it into an axiom.

Research consequence: weak-signal undercoverage alone cannot reject an
otherwise valid total mass law, nor select the default prior. The restricted
null-neutrality/regular-width obstruction stays conditional on adopting its
width requirement. Quartic-peak erasure, known-empty discontinuities, and
unproved uniqueness are distinct issues, not explained away by Fieller.

## User hypothesis: apparent flatness supplies unwarranted information

The user asked whether the flat measure itself makes the result too certain.
This has an exact mechanism in the common-noise case, and is already the
diagnostic in the managed ours26_learned_omega_is_a_hint record (ours26/,
not ours27/), reread at the pinned commit above together with the learned
push-scale derivation. Its prior rejection as an adopted solution stands.

For standardized x,y and theta=atan(m/s), set A=x sin(theta)+y cos(theta),
H=||A||^2, Q=||x cos(theta)-y sin(theta)||^2 and T=H+Q. Giving the candidate
stack v the proper N(0,omega^2 I_M) law yields marginal evidence
(1+omega^2)^(-M/2) exp[-T/2 + beta H/2], beta=omega^2/(1+omega^2).
At fixed omega this has theta dependence exp[-beta Q/2]. Conditional
v is N(beta A,beta I), so changing the measure also changes the joint
candidate law and readout weight. Flat dv is the omega-to-infinity limit
of normalized angular laws. It is not averaging over unknown omega.
The correction multiplier 1/beta concerns local VARIANCE; standard
deviation/interval width scales as 1/sqrt(beta) in the local normal regime.

Flat Cartesian dv carries radial factor r^(M-1) dr dOmega, so it is not
flat in candidate magnitude. Isotropy permits arbitrary radial weighting;
it alone cannot select Lebesgue measure. This supports the user's causal
diagnosis but does not select a replacement radial law or make fixed finite
omega a sanctioned assumption. Do not call the improper measure literally
a probability assertion that all pushes are infinite. The proven statement
is its infinite-scale limiting equivalence within this Gaussian family.
