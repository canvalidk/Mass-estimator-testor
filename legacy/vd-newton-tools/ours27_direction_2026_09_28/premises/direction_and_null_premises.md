# Direction, measured zeros, and known-empty limits

Scratch premise audit, 2026-09-28. No managed or production changes.
No replacement equation is adopted here.

## 1. What the managed sources actually require

All managed sources were read through the installed GitHub connector at
canvalidk/VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8, using its catalog.

**Ratified shared quantities.**
`Newton-analysis/_dscn_mass_estimator/why_reconciliation_cannot_be_upstream_of_E10.md`
ratifies on August 1 independent obtainability of the force and acceleration
inputs, and places their reconciliation inside mass determination. The
ratification does not prescribe a candidate measure, independence of candidate
components, or neutrality of a noisy zero component. In fact its closing
scope explicitly leaves cross-trial persistence open. Learning candidate
direction or size inside E10 leaves independently obtained original inputs
intact and does not reopen the ratified argument.

The catalog now places the September 6 shared-quantity expansion under
`unsanctioned/`. Its old root path returns 404. It was not promoted into this
argument; the ratified August source above is sufficient.

**Rotation.**
`Newton-analysis/_dscn_mass_estimator/ours27/ours27_full_covariance_2026-09-27.md`
requires simultaneous transformation of observations and covariance. Lines
47-52 permit correlated component and channel errors; input obtainability
therefore cannot be read as probabilistic independence of all these noises.
Rotation equivalence preserves the dimension and the full experiment. It
does not say that deleting zero-valued observed components is harmless.

**Axis splitting and its extra combination assumption.**
`Newton-analysis/_dscn_mass_estimator/axis_splitting_and_the_24_september_change_2026-09-25.md`
states the elementary bookkeeping requirement in section 1. Immediately after
it, however, "combined" is defined as separately integrating each reading's
latent pair and multiplying the resulting evidence. Section 4(b) then infers
a product candidate weighting over axes and uses isotropy plus no size scale
to reach flatness. The operational bookkeeping claim does not itself force
that particular nuisance-state elimination order. A shared nuisance state can
be retained while axes are split and regrouped. The same record, section 4(a),
explicitly says dropping y and z readings is less data, not regrouping.

The later full-covariance record, lines 224-230, says L8 as stated only applies
when covariance is diagonal in the chosen axes and owes restatement. The
noise-free stacked candidate-size requirement survives separately.

**Measured-null neutrality.**
`Newton-analysis/_dscn_mass_estimator/ours27/ours27_candidate_size_and_null_readings_2026-09-28.md`
section 5, lines 144-146, describes the property as it holds for ours27:
known-empty axes do not change the answer; measured-empty axes leave the mass
law unchanged and change the point by an amount vanishing with noise.
The document says at line 9 that it edits nothing and binds nothing. Its
clean version is restricted to one noise pair (lines 160-161). The text records
a resolution within the existing equation; it does not prove that rotation,
shared quantities, or instrument-only defaults force this property for every
future candidate measure.

**L10 is already open.**
`Newton-analysis/_dscn_mass_estimator/ours26/ours26_learned_omega_is_a_hint_2026-09-27.md`
rejects the learned-omega proposal as adoption, keeps it as a hint, and
expressly locates the width issue in what L10 should mean. The correction
does not require all finite-scale alternatives to be rejected forever.

Thus changing noisy-null feedback is not a violation of the ratified
shared-quantity argument. It does change a documented ours27 property and
the stronger product-after-nuisance-integration reading of axis splitting.
Those changes must be explicit, not disguised as an algebra correction.

## 2. What the rotation story exposes

In one common noise pair, write standardized observations x,y and

    A_theta = x sin(theta)+y cos(theta),
    H_theta = ||A_theta||^2,
    Q_theta = ||x cos(theta)-y sin(theta)||^2.

Under ours27, conditional candidates are v|theta,z ~ N(A_theta,I). When the
observed vectors lie on the first axis, transverse candidate coordinates are
still independent N(0,1). Candidate directional uncertainty was not omitted.
Flat integration makes their contribution to the mass evidence constant:
the cancellation removes mass-law feedback, not the candidate uncertainty.

A different joint size/direction measure can couple the transverse
uncertainty to the active candidate size. Then zero-valued noisy components
can alter the relative weights of candidate masses. Rotation covariance
permits this. No hidden true mass is needed: the update acts on hypothetical
aligned candidates using the independently obtained noisy observations.

There is also a sampling caution. If axes are selected from the observed
vectors, the resulting zero components were obtained by that selection;
conditional on the selected frame they are not new independent Gaussian
measurements. Pointwise rotation of data, covariance, and candidates is valid.
Deleting dimensions afterward and resetting an estimator is a different
experiment. A rotation-invariant likelihood still retains the original
dimension and its directional integration measure.

## 3. A precise weakening that permits learning

Keep these requirements:

* coordinate changes and axis splitting retain the SAME full joint candidate
  model, covariance, and shared nuisance state;
* the candidate measure is specified without knowledge of the object's mass;
* likelihoods of independent measurements multiply conditional on every shared
  state, with the original observations and residuals retained;
* physically absent candidate coordinates are removed as structural constraints;
* all-zero readings introduce no preferred mass angle when the declared model
  has the corresponding mass symmetry.

Replace exact measured-null mass-law neutrality with:

    A noisy zero can inform shared candidate size/direction. Any influence on
    mass must arise through that declared joint model, not an extra mass tilt
    attached by hand to the reading.

Consequently an all-zero dataset can have p(m)=pi while still updating size.
The full joint state must be carried forward. Carrying only p(m)=pi would
forget information that may matter once a nonzero reading arrives. The
earlier inference "mass law unchanged, therefore no information was supplied"
is specific to the flat model and is not a general consequence of a zero value.

This weakening preserves elementary axis bookkeeping. It abandons independent
axis nuisance integration before combination, not axis invariance itself.

## 4. The dt/t candidate has a mathematical reference justification

In the declared Gaussian covariance family

    z | theta,t ~ N(0, I+(t-1)P_theta),    t>=1,
    v | t ~ N(0,(t-1)I_M),

the information for total projected variance t is I_tt=M/(2t^2). Conditional
Fisher/reference weighting therefore gives dt/t, up to a constant. Together
with a uniform angle in the common-noise geometry it produces

    p(theta|z) proportional to
      pi(theta) exp(-Q_theta/2)
      gamma(M/2,H_theta/2)/(H_theta/2)^(M/2),

with the continuous H=0 limit used. The t posterior is proper for M>0:
its large-t integrand decays as t^(-M/2-1). At zero data the theta law is pi.
Appending a measured-zero coordinate adds a factor t^(-1/2), which learns
size and can alter the mass law when other readings are active.

This removes an arbitrary Beta shape parameter inside this chosen family.
It does not derive the family itself. A common Gaussian variance, its role
as the nuisance state, and reference-information weighting remain substantive
premises. Calling dt/t a Haar measure needs qualification: t>=1 is not closed
under the full multiplicative scale group. The clean statement is that dt/t
is the Fisher length measure of this conditional variance experiment.

The same-law composition theorem remains available for any normalized joint
candidate law produced this way. It does not select this reference principle.

## 5. Known-empty is not automatically the epsilon-noise limit

There are three distinct operations:

1. **Structural absence:** a component is declared exactly zero in every
   candidate, independently of all size states.
2. **Increasing precision:** a component remains an unknown candidate variable,
   but a zero observation is supplied with decreasing noise.
3. **Changing the default experiment:** the instrument covariance changes and
   the candidate measure or mass starting law is reconstructed from the new
   design, possibly changing retained components too.

Their equality is an additional compatibility requirement. It is not implied
by rotation or by independent input obtainability.

An exact one-coordinate model makes the distinction visible. Let an uncertain
transverse candidate coordinate be b=delta sqrt(t-1) xi and measure

    Y = b + epsilon eta,   xi,eta independent N(0,1).

For Y=0 the likelihood factor for t, after removing a t-independent constant,
is

    [1 + (delta/epsilon)^2 (t-1)]^(-1/2).

* If delta/epsilon -> 0, the factor tends to 1. The candidate coordinate is
  known absent more tightly than the instrument can distinguish. It can be
  removed continuously without changing the active joint law.
* If delta=epsilon, the factor is t^(-1/2) for EVERY epsilon>0. This is the
  noise-relative shared-variance model. Both physical candidate spread and
  measurement error vanish, but their ratio and information about t do not.
* If epsilon -> 0 with delta fixed, the observation increasingly constrains a
  physical unknown coordinate. It generally continues to inform shared size;
  exact zero is then conditioning information, not deletion of a dummy input.

Therefore the earlier epsilon tests correctly demonstrated failure of the
specific identity "unconstrained measured zero at arbitrarily small noise
equals a lower-dimensional fresh default fit". It was too strong to infer
from that alone a violation of the shared-quantity or rotation premises.

This is a singular statistical limit: weak collapse of b to zero is not enough
for posterior irrelevance when observation precision improves at the same
rate. The standardized experiment still distinguishes t.

## 6. Two coherent versions of known-empty continuity

**Structural restriction continuity.** Declare a geometrical constraint on the
candidate space and strengthen it toward exact absence, holding the retained
experiment and its measure fixed. The condition delta/epsilon -> 0 above
gives the desired continuous reduction. At exact structural absence the
independent added noise factor is parameter-free and cancels. This allows
measured-null learning for unconstrained axes while preserving the irrelevance
of truly absent axes. The constraint is information about the admissible
motion/geometry, not prior knowledge of mass.

**Conditioning continuity.** If exact emptiness is instead learned by precise
measurement, the lower-dimensional reference must be the conditional/coarea
restriction of the SAME joint physical candidate measure. It need not equal
the lower-dimensional default obtained by forgetting the constrained
coordinate and constructing a new measure. Comparing against that conditional
restriction can preserve continuity without asserting that the observation
was uninformative about the remaining shared state.

If the user specifically requires the old unconstrained epsilon-to-zero
experiment to equal the fresh lower-dimensional default, dt/t with noise-relative
candidate spread still fails it. Keeping that stronger condition would need
another measure construction or a justified design-dependent scaling such as
delta=o(epsilon); the latter is a new physical premise, not a free consequence
of instrument-only weighting. The earlier global-design escape constructions
exhibit mathematical possibilities but do not select one uniquely.

## 7. Does a precise added instrument necessarily change pi?

No. With one unchanged channel-noise ratio, a uniform angle starting law can
remain exactly the same as the extra axis's absolute noise tends to zero.
The dt/t conditional information factor sqrt(M/2) also disappears on
normalization; it does not change the measure's shape. Hence the shared-size
effect above cannot be dismissed as merely changing the explicit pi.

More general instrument-geometry rules can change pi, or the active candidate
measure, when an instrument is added. That is allowed only as an explicit
change of declared specification. `ours26_starting_law_2026-09-27.md`, lines
94-95, distinguishes changing that specification from updating with readings.
It does not prove that all design changes must alter pi. A continuity test
must state which design, starting measure, and structural knowledge are held
fixed, and must compare the induced retained model rather than silently reset
it. Changing the active prior and then demanding identical output confounds
different definitions.

An additional full-covariance caveat: even a structurally zero candidate axis
can provide information about other errors if its observation noise is
correlated with theirs. For covariance blocks (S,epsilon C;epsilon C',epsilon^2 V),
conditioning the zero-axis reading changes the retained noise covariance to
S-C V^(-1) C', for every positive epsilon. The covariance entries can converge
to a singular block matrix while the conditional experiment retains this
information. Dummy-axis irrelevance therefore also requires independent added
noise, or preservation of the appropriate conditional experiment. The old
independent-axis counterexamples do not use this caveat, but a blanket claim
for arbitrary Sigma would be too strong.

## Conclusion of the audit

Measured-null neutrality is a revisable modeling commitment, not a theorem of
the project's ratified shared-quantity structure. Rotation and elementary
regrouping are compatible with size/direction learning. A defensible revision
preserves structural known-empty restriction and full joint-state recoding,
while allowing noisy zeros to update shared nuisance information. The exact
epsilon-to-fresh-default continuity requirement is a separate, stronger
choice. The reference dt/t candidate makes that fork concrete but does not
settle which physical uncertainty law the definition must use.
