# A composition uniqueness proof for the mass readout

**Date:** 2026-09-06  
**Status:** new local mathematical handoff; unsubmitted to VD-docs. The theorem is proved under the explicit assumptions below. Its applicability is not a new ratified estimator policy.  
**Purpose:** replace the earlier restricted-family uniqueness sketch with a direct argument, and identify exactly what remains outside the result.

## Result and scope

For a readout determined by the **joint probability law of the positive force and acceleration magnitudes**, the two composition identities, point calibration, and a mild continuity requirement force

$$
R(f,\alpha)=\frac{E[f]}{E[\alpha]}.
$$

The proof does not choose a loss function or assume a weighted-mean family. It does not require Gaussian errors or independent magnitudes. The algebraic step applies to arbitrarily dependent magnitudes.

This is a theorem about the readout of an already supplied magnitude law. It does not uniquely construct that law from measurements. In particular, magnitude-law dependence must be stated explicitly; dependence on the full vector/direction law is a weaker restriction and admits the counterexample below.

## 1. Domain and assumptions

Choose reference units F0 and a0 and the corresponding mass unit m0 = F0/a0. Work with dimensionless positive magnitudes X = f/F0 and Y = alpha/a0; write R for the dimensionless mass readout. The two appearances of 1 below mean separate force and acceleration reference magnitudes, not an identification of their physical units.

The domain contains every joint law of strictly positive random variables X,Y with finite means. R is finite and strictly positive. Domain closure includes deterministic inputs, addition of positive constants, arbitrary common-channel couplings, and iid extensions. These closure assumptions matter: a restricted Gaussian-posterior family alone is not the theorem's domain.

**C1. Force composition, on one propagated joint law:**

$$
R(X_1+X_2,Y)=R(X_1,Y)+R(X_2,Y).
$$

**C2. Acceleration composition for inverse mass, on one propagated joint law:**

$$
\frac1{R(X,Y_1+Y_2)}=\frac1{R(X,Y_1)}+\frac1{R(X,Y_2)}.
$$

These requirements apply to every allowed coupling with the indicated common channel. They do not equate a propagated law with a fresh fit using new evidence or a newly assigned measure.

**C3. Magnitude-law dependence:** R depends only on the joint law of (X,Y). It cannot use probability-space labels, an extra latent variable, or further direction information absent from that law.

**C4. Point calibration:** R(x,y)=x/y for deterministic positive x,y.

**C5. Constant-limit continuity:** define T(X)=R(X,1) and S(Y)=1/R(1,Y). If positive variables X_n converge in mean absolute error to a positive constant c, then T(X_n) tends to c; likewise S(Y_n) tends to c. Full L1 continuity of T and S is sufficient but stronger than needed.

C5 is the required meaning of continuity; the earlier sketch's phrase 'in the usual sense' did not specify a topology. This proof does not need a separate monotonicity axiom or assume a particular regular family.

## 2. Composition alone forces separation

For any joint law of (X,Y), put

$$
r=R(X,Y),\qquad x=R(X,1),\qquad y=\frac1{R(1,Y)}.
$$

All three numbers are finite and positive. Compute R(X+1,Y+1) in two orders.

First split the force side, then use inverse-mass composition:

$$
R(X+1,Y+1)
=R(X,Y+1)+R(1,Y+1)
=\frac{xr}{x+r}+\frac1{y+1}.
$$

Alternatively split the acceleration side first, then use force composition:

$$
R(X+1,Y+1)
=\left[\frac1{R(X+1,Y)}+\frac1{R(X+1,1)}\right]^{-1}
=\frac{(r+1/y)(x+1)}{r+1/y+x+1}.
$$

Here R(1,1)=1 by calibration. The difference between the two expressions is identically

$$
-\frac{(ry-x)^2}{(x+r)(y+1)(yr+1+xy+y)}.
$$

They describe the same readout, so their difference is zero. The denominator is strictly positive, hence ry=x. Therefore

$$
\boxed{R(X,Y)=\frac{T(X)}{S(Y)}.}
$$

This is the step missing from the previous sketch. It handles every coupling directly, rather than checking product laws and then restricting to an assumed coefficient family. No distributional independence is used.

## 3. Law dependence and continuity force expectation

By C1, T is additive over arbitrarily coupled positive variables. By C3 it depends only on the marginal law of its argument, and by C4 it fixes positive constants.

Take iid copies X_1,...,X_n of an arbitrary positive integrable X, and let Xbar_n be their arithmetic mean. Additivity supplies integer and rational homogeneity, so

$$
T(\overline X_n)=\frac1n\sum_{i=1}^n T(X_i)=T(X).
$$

The sample mean converges in L1 to E[X]. One elementary verification truncates X at K: the bounded sample mean has mean absolute deviation at most sqrt(Var(min(X,K))/n), while the two tail corrections total at most 2 E[(X-K)_+]. First take n to infinity, then K to infinity.

By C5 and calibration,

$$
T(X)=\lim_n T(\overline X_n)=T(E[X])=E[X].
$$

C2 makes S additive. The identical argument gives S(Y)=E[Y]. Substitution into the separation identity proves

$$
\boxed{R(X,Y)=\frac{E[X]}{E[Y]}.}
$$

Restoring units gives R_physical(f,alpha)=E[f]/E[alpha]. Conversely this readout satisfies C1-C5 on the stated domain, so existence and uniqueness both hold.

## 4. What is stronger now

The previous audit passed from arbitrary joint laws into a proposed regular family midway through its proof. This argument removes that restriction: the square identity forces separation directly. It also gives an explicit, sufficient continuity condition and a precise integrable domain.

Reciprocity, exactness on any law supported on f=m alpha, and the readout's independence of the coupling when the two magnitude marginals are fixed follow from the derived formula. They do not have to be separately selected through a loss function.

The result answers the user's 'eleventh solution' objection at the readout level: under C1-C5 there is no alternative readout. It still matters whether the premises are independently required by the intended mass operation. In particular, C1/C2 are operational composition commitments, not consequences of writing F=ma alone.

## 5. A necessary caution about common direction

If 'law dependence' allows the readout to use the full law of (f,alpha,u), rather than only its induced magnitude-pair law, the preceding theorem cannot be asserted without C3.

Here is a concrete family showing why. Let u be a common unit direction, v=E[u], and h(u)=1+(u dot v)^2. Define

$$
R_h(f,\alpha,u)=\frac{E[h(u)f]}{E[h(u)\alpha]}.
$$

The same h is used whenever the common direction variable is retained during either composition. This readout satisfies both composition identities, reciprocal symmetry, units covariance, point calibration, and exactness when all states share one mass. It is invariant under simultaneous rotations, since u dot v is. It is determined by the full vector-pair law. It also satisfies constant-limit continuity for integrable magnitudes: since 1 <= h <= 2, the deviation of E[h X]/E[h] from a constant c is at most 2 E[abs(X-c)]. What it violates is C3: the readout can change when the magnitude-pair law stays fixed but its association with direction changes.

For example, let (u,f,alpha) be (e1,1,1) with probability 3/4 and (e2,2,1) with probability 1/4. Ordinary magnitude means give 5/4. Here v=(3/4,1/4), h(e1)=25/16, h(e2)=17/16, and the tilted readout is 109/92, which differs from 5/4.

This is not a recommended estimator. It exposes the exact premise needed to exclude a second direction-based weighting after the measurement model has already assigned probabilities. The paper should justify that once P has been constructed, the scalar readout uses its magnitude law without additional ancillary-direction weighting. Calling the full vector-law version of the earlier proposition proved would overstate this result.

## 6. What still needs finishing in the full estimator

1. **Justify the theorem's scope.** Establish the two common-channel composition requirements, domain closure, and magnitude-law sufficiency as requirements of the intended operation. The mathematics above then fixes the readout without choosing a loss.
2. **Specify the measure and measurement model.** The theorem accepts P as supplied. It does not choose the flat positive-magnitude measure over cone volume or another measure, nor Gaussian errors over an experiment-specific model. Separate genuinely experimental inputs from residual discretionary choices.
3. **Specify uncertainty and downstream standing.** Induced-ratio probabilities, fixed-mass compatibility, and repeated-sampling uncertainty answer different questions. Measured zero-zero can retain a positive point and unrestricted compatibility, but propagation and admission must preserve that distinction. Infinite ratio moments are not a unique null detector.
4. **Finish repeated-trial scope.** Repeated observations of one latent pair, different latent trials sharing one mass, and regrouping one existing law are distinct operations. The common-mass measure problem is not solved by the single-law readout theorem or by pooling independently fitted moments.
5. **Finish numerical and gate validation.** A reference implementation needs explicit integration error control, tail handling, and calibrated compatibility behavior, including boundary cases. More compute should improve numerical accuracy without changing the declared estimator or creating physical information.

These items do not by themselves call for a new point formula. Some are specification work, while the measure and the independent justification of the readout axioms remain substantive modeling questions.

## 7. Verification, provenance, and handoff

The two-order identity was checked by exact polynomial expansion with integer coefficients in a JavaScript calculation during this session. Its numerator is 2*x*y*r - y^2*r^2 - x^2, exactly -(ry-x)^2; the difference from that expression had every coefficient zero. The uniqueness proof itself is the argument in sections 2-3, not a numerical test or a restricted family search.

**Managed snapshot:** canvalidk/VD-docs@743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6. main was resolved again on 2026-09-06 and matched the catch-up snapshot. Source inventory: catalog.yaml. Sources used from that snapshot:

- Newton-analysis/_dscn_mass_estimator/what_forces_the_ratio_of_means_readout.md — section 4's composition commitments and incomplete uniqueness sketch.
- Newton-analysis/_dscn_mass_estimator/what_forces_the_chain_to_the_estimator.md — distinction between framework premises, mathematical consequences, and rulings.
- Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md — current construction, direction integration, output contract, and remaining policy list.
- Newton-analysis/_dscn_mass_estimator/mass_estimator_originality_research_2026-09-05.md — the proof gap and measure/uncertainty distinctions.
- Newton-analysis/_dscn_mass_estimator/mass_estimator_additional_derivations_2026-09-05.md — conditional-on-P scope, finite-loss derivation, and aggregation limits.

**Local/session context:** the user's latest direction in the prior task was to eliminate arbitrary choice rather than add a preferred eleventh solution; the present task asks whether the function is essentially finished and can be made more forced. Existing E10 production files and the three previously deleted local staging copies were not modified or restored.

**Coverage:** new internal derivation, not a fresh literature/originality search or an exhaustive context audit. No historical-priority claim is made. No managed source has been edited; no source correction or policy has been ratified by saving this file.

**Submission:** saved under the workspace-root pass to VD-docs/ policy. Awaiting transfer; not submitted to VD-docs.
