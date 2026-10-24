# Continuation after the new September 28 sources

Scratch only. No new human-facing records requested or created in this turn.
No production or managed-document changes.

## Source snapshot and scope

Resolved VD-docs main once through the GitHub connector to
592d144531fa4a30d1c2229a4c3b528e84043d77. Its sole parent is the prior snapshot
b4bfc2040031ef010416334d0c5f4aa2ad05cbd8. The commit adds eight records, modifies
the ours27 README and catalog, and changes no other managed files. Four of the
eight are our previous four local records, now received into ours27; leave
the local copies alone. The other four were read in full and cached under
sources/ here:

- ours28_candidate_push_scale_as_nuisance_2026-09-28.md
- ours28_size_measure_jeffreys_check_2026-09-28.md
- ours28_known_empty_vs_measured_empty_2026-09-28.md
- ours27_c3_from_betweenness_and_the_scale_fork_2026-09-28.md

The catalog's relevant entries were also read. The older canonical sources
are unchanged between the snapshots. The readout agent fetched the September
6 composition source at the new snapshot; its blob remains
1bb77de5c16b1d1102b2e74866e5a4219ed10c53.
Local scope includes our prior four records and relevant .tools investigations.
This is a scoped continuation of the mass-estimator derivation, not an
exhaustive synthesis of all VD documents.

User steering: read the new documents and work on solving the equation;
do not spend the turn correcting or coordinating other user-owned sessions.
No messages were sent to other user-owned chats.

## Root result 1: fixed per-reading shrinkage cannot close regular widths

For common-noise experiment at a tested candidate angle,
 A ~ N(0,t I_d), B ~ N(0,I_d), independently, t=1+omega^2.
For a fixed per-reading radial evidence F(H), H=|A|^2,
 kappa=2 dlogF/dH, score=kappa(H) A.B.
Exact identities:
 V=E[H kappa(H)^2],
 J=-E ell''=(1-1/t)E[H kappa(H)].
Gaussian integration by parts gives the second identity from
E[kappa(H)(H-d)-2 H kappa'(H)].

Under regular independent replication, the mass-law variance is 1/(NJ);
the sampling variance of its local center is V/(NJ^2). These concern local
law/root widths, with the earlier smooth-readout qualification still in force.

If V=J for EVERY t>1, put q=1/(2t). Uniqueness of Laplace transforms yields
 2 kappa' + (d/H-1)kappa + kappa^2 = 0,
 equivalently 2 H F''+(d-H)F'=0.
The boundary term is zero because H^(d/2)kappa(H)->0 at H=0.
A positive finite Gaussian scale-mixture evidence is analytic at zero, with
kappa(0)>0 for any nonzero beta mixture; this ODE requires kappa(0)=0.
Thus no nontrivial such FIXED PER-READING law can satisfy all-scale regular
width matching. This does not reject a shared scale/distribution learned from
the data, or calibration averaged under the law's own prior experiment.

For the arcsine mixture, d3, omega=.5:
 V=.6705896324340224, J=.3006954206277883,
 V/J=2.230129181994103, asymptotic central-95 coverage=.8106319735532941.
The signal per reading is fixed positive and information grows with N;
this is not the Fieller nonidentification limit. Finite N central-mass cells
can cover more because the posterior remains broad on the bounded angle.

Verified by two independent gamma-weighted quadratures of J, and 240->400
refinement. See local_scale_width.py and results. This calculation directs
the constructive search toward LEARNING shared states, not tuning a new
fixed Beta shape or applying a sandwich multiplier.

## Root result 2: an actual learned-group law repairs the heterogeneous example

For independent channel-pair blocks, group identical declared noise laws.
Group g has M_g candidate coordinates and projected squared data H_g(m).
Use a proper joint starting measure on beta_g in (0,1), independent of m:

 p(m,d beta|D) proportional pi(m) P_Sigma(d beta)
   product_g (1-beta_g)^(M_g/2) exp[beta_g H_g(m)/2].

For initially independent arcsine betas this is
 pi(m) product_g beta_g^(-1/2)
   (1-beta_g)^((M_g-1)/2) exp[beta_g H_g(m)/2].
Keep beta states across readings, not a fresh draw per reading. Conditional
candidate coordinates v_g are N(beta_g A_g,beta_g I).
Use the existing ratio E||F*|| / E||a*|| on the resulting joint law.

Let u=log m, w_g=dtheta_g/du at the tested mass, and let fixed candidate
coefficients have mean square C_g>0 in the regular replicated limit.
Then H_g/M_g ->1+C_g and posterior beta_g -> C_g/(1+C_g).
The mass and beta score blocks are orthogonal at the test candidate.

 J=sum_g M_g beta_g C_g w_g^2,
 V=sum_g M_g beta_g^2(C_g+1)w_g^2=J.

This holds for fixed non-Gaussian candidate sequences under appropriate
concentration/no-dominating-reading conditions; only their second moments
enter. Finite fixed number of groups and nonsingular positive prior density
are part of the regime.

On the previous heterogeneous two-axis case (force SD1,10; acceleration SD1,1;
m1; candidate squared sizes.1,9.9), beta limits are .09090909,.90825688.
J=V=.09041842984371105 per repetition; variance*N ->11.059692163738164.
The single shared-scalar law had variance ratio2.60722676; this law's ratio
is exactly1 in this regular limit. This also changes the efficient center,
so do not describe it as merely widening the old center's law.

grouped_learning.py verifies exact rational arithmetic, integration,
sequential combination and independent-null factors. Hard grouping is
discontinuous at equal noise; do NOT present that hard rule as the final
full-covariance definition.

## Root result 3: continuous coupling from noise-law overlap

A simpler continuous alternative to hard grouping has now been constructed
for independent paired-channel blocks. Let S_g be the SPD 2x2 Gaussian noise
covariance of a block, including force/acceleration correlation if supplied.
The centered Gaussian density is n_g.

 rho_gh = integral sqrt(n_g n_h)
        = |S_g|^(1/4)|S_h|^(1/4) / |(S_g+S_h)/2|^(1/2).
 Z ~ N(0,rho), beta_g=sin^2[(pi/2)Phi(Z_g)].

Use this P_Sigma in the learned-group equation above. This is an explicit
Gaussian coupling principle on square-root noise laws; it is not derived
uniquely from the current mass axioms. It avoids a free correlation-length
parameter and an arbitrary grouping threshold.

Properties with proofs:

- rho is a PSD Gram matrix of square-root densities; distinct Gaussian noise
  densities are linearly independent, so finite distinct types give positive
  definite rho. The field has consistent marginal laws as instrument types
  are added. Identical noise laws have rho=1 and share one Z/beta.
- Each beta is exactly arcsine, since Phi(Z) is uniform.
- Common channel unit changes and channel swapping are changes of variable
  in the overlap integral, so leave rho unchanged.
- Covariances approaching equality make rho->1 continuously. No discrete
  merge is used. For commuting spatial blocks the common degenerate noise
  eigenspace gets a single scalar beta, so arbitrary basis changes there
  cannot change the equation.
- For S_active=I_2 and S_null=epsilon^2 I_2,
  rho=2epsilon/(1+epsilon^2)->0. Noisy zero directions gradually decouple
  as their absolute precision improves. Their likelihood then becomes a
  mass-independent independent factor. Their physical candidate covariance
  is O(epsilon^2), so they vanish in the ratio-of-mean-norms point too.
- All-zero data leave p(m)=pi. An added finite-noise zero may still inform
  related scale states. Literal arbitrary-zero neutrality is NOT claimed.
- A finite fixed collection of distinct instrument types has a positive
  smooth prior density on the beta interior; its prior correlations wash
  out in regular replicated inference, so the J=V calculation survives.

For a general 2x2 noise block S_g, use
 u_g(m)=S_g^-1/2(m,1)' / sqrt((m,1)S_g^-1(m,1)'),
 z_gi=S_g^-1/2(F_gi,a_gi)', H_g=sum_i(z_gi.u_g)^2.
The displayed likelihood is exact. Across-axis correlations not reducible
to declared independent blocks still require a full matrix construction.

The law is proper at every finite input for proper pi: H_g(m)<=sum||z_gi||^2,
so all data factors are uniformly bounded; they are positive for interior
beta. The physical conditional candidate means/covariances follow from
beta times the flat conditional Gaussian, giving a finite positive norm
ratio under the same endpoint bounds as ours27.

Numerical aligned example F3,a1, two zero axes:
max angular-density difference from the 1D law is
 .18293738 at epsilon1,
 .01765442 at .1,
 .00166510 at .01,
 .00016552 at .001,
 0 at the independent limit.
The epsilon1 endpoint equals the scalar 3D arcsine law; the zero endpoint
equals the 1D law, each checked to <1.4e-11. Quadrature160->200 changes the
intermediate density by4.3e-12. Units/swap errors5.6e-16.
See noise_overlap_candidate.py and results.

This is a candidate equation for independent instrument blocks, not an
announcement of the project's final equation or a general uniqueness claim.

## Full-covariance construction developed in parallel

geometry/instrument_covariance_field.py defines an instrument-only Gaussian
matrix field. Its exact joint likelihood is

 p(m,dR|y) proportional pi(m) P_Sigma(dR)
   product_i |I-Beta_i|^(1/2)
     exp[-T_i/2+a_bar_i' Beta_i a_bar_i/2].

Here K_m=(B_m'Sigma^-1B_m)^-1, b_m=K_m B_m'Sigma^-1 y,
a_bar=K_m^-1/2 b_m, T=y'Sigma^-1 y, Beta=R(I+R)^-1.
The conditional physical acceleration is
 N(K_m^1/2 Beta a_bar, K_m^1/2 Beta K_m^1/2).

The shared field's instrument descriptors and stationary covariance are fully
specified in the script; full covariance, independent frames, units and swap
checks are below9e-16. Eigenvalue-coincidence continuity and an independent
precise-null limit have been verified. Repeated diagonal groups recover
the same J=V result.

This proves a compatible smooth full-matrix construction EXISTS, but its
kernel, log-width, reference prescription and sharing are extra choices.
Arbitrary full-covariance width calibration is not proved. Do not conflate
it with the simpler noise-overlap construction; they illustrate different
ways to meet the geometric conditions, hence those conditions alone have
not selected a unique prior.

## Readout advance

readout/finite_slice_and_covariant_score_proof.md closes the finite-support
flat-slice gap directly: strict betweenness + continuity on finite simplices
give a unique (up to positive multiple) fixed score at each output level r.
Readout covariance under an involution forces its score multiplier to be1.
Reflections generate O(M), so full stacked orthogonal recoding forces the
score to depend only on (||F*||,||a*||), without requiring the candidate law
or its generating measure to be isotropic.

Wasserstein-1 continuity extends this to integrable continuous candidate laws;
the existing composition theorem then gives ratio of mean norms.
The current own-frame/sign rulings alone give block O(d), leaving the
per-reading size-profile survivor. A sufficient extra recoding operation is
invertible orthogonal recombination of scalar force channels across readings
and the identical recombination of acceleration channels, carrying the full
joint law/covariance. This need not reset scale sharing.

The extra W1 regularity is real: old continuity at constants alone admits an
atomic-direction weight1 / nonatomic-direction weight2 counterexample.
Do not claim unconditional C3 uniqueness from the existing listed rulings.

## Tensor alternative and limits

tensor/proper_tensor.py contains:
- proper matrix-Beta(d/2,d/2) / Grassmann construction;
- proper projective R=T*IW_d(d,I), T~chi-square1, which recovers half-Cauchy
  scalar marginals but loses matrix inversion symmetry;
- exact collapse of the one-reading tensor to a scalar radial mixture;
- shared tensor calibration for common instrument geometry, and the failure
  of automatic transfer to arbitrary heterogeneous designs;
- the need for declared frame transport when actual tensor orientations
  are shared across readings.

Scalar instrument-defined sharing is the cleaner immediate constructive path.

## Next work

1. Derive or replace the noise-overlap joint measure from an explicit
   operational principle; do not call a convenient Gaussian coupling unique.
2. Generalize the simple block construction to arbitrary full covariance
   while keeping its regular width theorem, or prove the matrix-field law's
   calibration with its precise support and sharing conditions.
3. Resolve whether readout-only stacked recoding and candidate-cloud
   continuity follow from the project's existing operational rules.
4. Preserve the no-meta-mass reading: repeated candidate experiments are
   diagnostic tests, not a replacement definition of mass.

All generated work is scratch. The user's previous four human-facing records
were left unchanged, as was unrelated active code.
