# Measured Zero–Zero: Two Cases for the Mass Operation

**Date:** 2026-09-05  
**Status:** competing argument drafts; neither option is ratified by this document  
**Purpose:** develop each option in its own committed voice, so that its premises, explanatory strengths, and consequences can be examined.  
**Local status:** active, unsubmitted VD-Newton analysis and handoff; not a replacement for the managed source notes.

## The question

When an empirical mass determination receives zero central readings for both net force and inertial acceleration, with nonzero uncertainties, should it return:

1. **mass undefined from this trial**, or
2. **a positive scalar mass estimate**, which under the current symmetric policy is the apparatus scale $\sigma_F/\sigma_a$?

Throughout, $0/5$, $5/0$, and $0/0$ abbreviate **force–acceleration measurement inputs**, with appropriate units. They are not literal vector-division operations. The one-zero formulas use a one-dimensional, known-direction limit with the nonzero channel well resolved.

Both cases preserve the distinction between a measurement result and an explicitly exact premise. Neither proposal licenses correcting an exact contradiction. The disagreement concerns the empirical zero–zero input and what counts as a successful mass determination.

**Uncertainty development, added 2026-09-05:** the user proposed that the estimator may already have the desired behaviour: its point value remains positive while the accompanying mass uncertainty exposes the absence of a resolved mass. The two cases below now address that proposal. The key distinctions are a finite point estimate, divergent moments of an induced mass distribution, and a compatibility set containing every positive mass.

---

# Case I — A Null Trial Must Be Allowed to Leave Mass Undefined

## 1. The mass operation must preserve an information boundary

I defend a mass-determination operation that can conclude:

> Newton II is compatible with this trial, but this trial does not determine a mass.

That conclusion is a successful analysis. The operation has located the boundary of what its evidence supplies. Replacing that boundary with an apparatus-selected number would make the output look more complete while weakening what it means.

Mass remains a positive scalar wherever a mass value is established. The determination operation is partial as a function into positive mass, or total as a function into an explicit result type:

```text
MassDeterminationResult :=
    MassEstimate(positive_scalar, standing)
  | MassNotIdentified(reason)
  | PoorModelCompatibility(reason)
```

This distinction matters. We need not claim that the object has an undefined physical property. We claim that this route has not supplied a value for it. Any later operation requiring that value must obtain it through another licensed route or halt.

## 2. The one-zero cases contain a constraint that zero–zero lacks

Consider measured $0/5$. Acceleration is resolved away from zero. A sufficiently large mass would require a force that the force instrument should have resolved. Its non-detection therefore constrains the mass toward the small-mass regime.

For measured $5/0$, force is resolved away from zero. A sufficiently small mass would produce an acceleration that the acceleration instrument should have resolved. Its non-detection therefore constrains the mass toward the large-mass regime.

Under a declared positive-magnitude averaging policy, representative values in those regimes can be

$$
\widehat m_{0/5}\approx\frac{c\sigma_F}{5},
\qquad
\widehat m_{5/0}\approx\frac{5}{c\sigma_a},
\qquad c=\sqrt{\frac2\pi}.
$$

These point values are policy-dependent summaries. The evidence itself usually supports a one-sided restriction more directly than it supports one precise mass. Nevertheless, there is a physical restriction to summarize: some masses predict a signal incompatible with what was observed.

At zero–zero, that particular mechanism disappears. Every positive mass admits the underlying state

$$
\mathbf F^*=\mathbf0,
\qquad
\mathbf a^*=\mathbf0.
$$

No mass needs a large force or acceleration merely to account for this pair of readings.

## 3. A residual calculation makes the distinction explicit

For isotropic independent Gaussian errors, the fixed-mass profile discrepancy is

$$
Q(m)=
\frac{
|\widetilde{\mathbf F}|^2
+m^2|\widetilde{\mathbf a}|^2
-2m\widetilde{\mathbf F}\cdot\widetilde{\mathbf a}
}{\sigma_F^2+m^2\sigma_a^2}.
$$

At measured zero–zero,

$$
Q(m)=0\qquad\text{for every }m>0.
$$

At $0/5$,

$$
Q(m)=\frac{25m^2}{\sigma_F^2+m^2\sigma_a^2},
$$

which penalizes large masses when the acceleration is well resolved. At $5/0$,

$$
Q(m)=\frac{25}{\sigma_F^2+m^2\sigma_a^2},
$$

which penalizes small masses when the force is well resolved.

A compatibility allowance therefore distinguishes a proper restricted region of positive masses from the entire positive ray. This supplies a principled foundation for accepting a qualified summary in the one-zero cases while declining to select a mass in the null case.

The profile calculation does not itself produce the finite point estimates above: its one-zero optima approach zero or infinity. I use it to establish the evidential restriction, and a declared summary policy to report a representative positive value where such a restriction exists. The two operations have different jobs.

## 4. Positive averaging introduces a choice that the null data do not require

Under the current flat latent-magnitude policy, zero readings produce half-normal distributions over positive magnitudes:

$$
E[f]=c\sigma_F,
\qquad
E[\alpha]=c\sigma_a.
$$

Their ratio is

$$
\widehat m=\frac{\sigma_F}{\sigma_a}.
$$

The calculation is coherent. Its coherence does not oblige us to admit the number as the trial's mass determination.

Integration weights the neighbouring compatible states according to a selected latent measure. That weighting supplies information not contained in the statement that every mass has a perfectly fitting zero–zero state. A different experimental or prior model can assign different relative support to masses.

If the experiment establishes a truly unexcited state, and its instrument noise is independent of mass, observations in that state have no mass-dependent likelihood. An apparatus scale cannot acquire object-specific meaning merely because the output field is named `mass`.

## 5. Removing free particles does not remove balanced interacting objects

The proposed free-particle separation is useful. Established absence of interactions can answer the relevant motion question without attempting this mass inference.

But an object with interactions can also have a zero resultant. Equal-and-opposite contributions can balance. The object remains interacting even though its net force and acceleration vanish.

Uncertain measurements of those contributions do not establish exact cancellation. They also do not guarantee that the trial identifies a mass. If the resultant and acceleration are both unresolved, a large range of masses can remain compatible with the evidence.

The fact that force-sum produced the input establishes its provenance. It does not automatically establish its identifying power. Interaction evidence and mass evidence must retain their separate roles.

## 6. Continuity of a formula does not require continuity of admission

The positive-magnitude formula varies continuously as the measurements approach zero. That is a useful numerical property. It is not a reason to forbid an explicit non-identification result.

A determination can lose identifying power while its numerical summary continues to exist. The trace must be allowed to record that loss.

I therefore do not propose a rule triggered solely by whether rounded displays contain the character `0`. The admission rule concerns the mass information in the full measurement result. Exact central zero–zero with a flat profile is its clearest example; sufficiently weak nearby cases may also fail a declared identification requirement.

Such a rule requires an explicit policy for sufficient identification. That is a proper burden for this design. It places the decision where it belongs: at the point where a numerical candidate is accepted as an answer that other entries may use.

## 7. The type must describe what the operation has earned

The strongest objection to my position is that an undefined result prevents mass from being a total positive-scalar function. I accept that consequence.

If the evidence has not supplied a mass, making a scalar available for multiplication does not repair the missing evidence. It permits unsupported calculations to proceed. A trace that stops with `MassNotIdentified` preserves more meaning than one that supplies an arbitrary finite value with an easily neglected warning.

The hypothetical force-free Sun makes the danger visible. A calculation may return $1\,\mathrm{kg}$ because that is the apparatus scale. Reusing it as the Sun's mass to predict a later acceleration would exceed what the null trial established. The determination should expose that missing licence at its own output.

## 8. Infinite uncertainty strengthens the need to distinguish a candidate from a determination

The strongest new reply to my position is that the estimator can return $1\,\mathrm{kg}$ while its uncertainty calculation already says that no finite precision has been achieved. Perhaps the desired information boundary is present in the complete result without making its mass field undefined.

I accept the mathematical finding. Under independent Gaussian uncertainties and a flat measure over positive latent magnitudes, the induced ratio $M=f/\alpha$ at zero–zero has a half-Cauchy distribution with scale $s=\sigma_F/\sigma_a$. The ratio-of-means estimate and the induced median both equal $s$, but $E[M]$ diverges and there is no finite standard deviation. The expected squared deviation from the reported point is infinite.

That finding must not be read as equal support for all positive masses. With $s=1\,\mathrm{kg}$, the same distribution puts 95% of its probability between approximately $0.039$ and $25.45\,\mathrm{kg}$. It still assigns strong relative preferences through the chosen latent weighting. Infinite squared-error uncertainty does not erase those preferences or establish that the Sun's mass has been adequately represented.

The profile-based compatibility result is different and more decisive for my argument:

$$
\{m>0:Q(m)\leq q\}=(0,\infty)
\qquad\text{for any }q\geq0
$$

at central zero–zero. Every positive mass passes. This result states the non-identification directly, without assigning relative probability by a latent measure.

I welcome returning that set beside a diagnostic candidate. But if the complete output says that the trial excludes no positive mass, the trace should state `MassNotIdentified` explicitly. The availability of an apparatus-scale point should not weaken that conclusion.

Nor can infinite standard deviation serve as a special zero–zero admission test. In the simple $5/0$ limit, positive acceleration support extends arbitrarily close to zero, so $M=5/\alpha$ also has a divergent high-mass tail. The useful distinction remains which masses the evidence excludes, rather than whether a second moment exists.

## Conclusion of Case I

**An empirical mass operation must be allowed to return no mass value when its evidence does not identify one.**

The one-zero cases can constrain mass toward small or large values because one channel supplies a resolved physical reference. Zero–zero can lack that restriction altogether. A policy-generated apparatus scale is available for diagnostics, but its existence does not require admitting it as a successful mass determination.

I retain positive mass as the value type and make failure to determine such a value explicit. An unrestricted compatibility set is valuable evidence for that result; a finite candidate accompanied by divergent moments does not overrule it. The surrounding entries must respect that information boundary.

---

# Case II — Uncertain Zero–Zero Must Remain Inside the Positive-Scalar Estimator

## 1. The empirical operation must be defined on its actual inputs

I defend one positive-scalar estimation operation for all admissible empirical force–acceleration results, including measured zero–zero.

The input is never merely the pair of central readings. It includes the uncertainty model and the estimation policy:

$$
(\widetilde{\mathbf F},\widetilde{\mathbf a},\Sigma,\Pi).
$$

A zero central reading does not discard the rest of this input. The estimator must continue to use it.

Declaring measured zero–zero undefined because the exact equation $\mathbf0=m\mathbf0$ admits every mass would apply an exact-input answer to a richer empirical question. The empirical operation asks which compatible possibilities are supported under the declared policy. It is entitled to use the uncertainty information that makes this question different.

## 2. Establish freeness before asking the estimator to represent it

The decisive distinction begins before arithmetic.

If a particle is established to have no interactions, the free-particle account can supply its motion without determining mass from zero net force and zero acceleration. That case need not be represented as an ordinary null measurement inside this mass-inference route.

For an interacting object, in the law-house considered here, the route to net force that does not already consume the demanded mass is force-sum:

$$
\widetilde{\mathbf F}_{\mathrm{net}}
=\sum_i\widetilde{\mathbf F}_i.
$$

Each contribution has a measurement standing. The resultant inherits propagated uncertainty, including correlations:

$$
\Sigma_F=\sum_{i,j}\operatorname{Cov}(\boldsymbol\epsilon_i,\boldsymbol\epsilon_j).
$$

Consequently a zero central sum can coexist with unresolved positive and negative resultant components. The sum says what the measured contributions add to. It does not by itself establish exact physical cancellation.

This is the zero that the empirical estimator must handle: **an unresolved resultant of an interaction account**.

## 3. The mechanism that repairs one unresolved channel also handles two

In measured $0/5$, a zero force reading permits positive latent force magnitudes. Under the current Gaussian policy their mean is approximately $c\sigma_F$, so the mass readout is

$$
\widehat m\approx\frac{c\sigma_F}{5}.
$$

In measured $5/0$, a zero acceleration reading permits positive latent acceleration magnitudes. Their mean is approximately $c\sigma_a$, giving

$$
\widehat m\approx\frac{5}{c\sigma_a}.
$$

The same operation at measured zero–zero gives

$$
\widehat m
=\frac{E[f]}{E[\alpha]}
=\frac{c\sigma_F}{c\sigma_a}
=\frac{\sigma_F}{\sigma_a}.
$$

No new arithmetic has been introduced. Both channels receive the treatment already accepted when either one is unresolved.

To accept the positive latent mean in the numerator of $0/5$ and the denominator of $5/0$, but forbid their joint use at $0/0$, requires an additional admission rule. The positive-scalar operation needs no such interruption.

## 4. A perfectly fitting point does not exhaust an uncertain input

The objection that every mass can fit the latent origin is correct. It establishes that the best achievable residual alone does not distinguish the masses.

But our estimator is not defined as the selection of the best latent point. It weights the full family of plausible compatible pairs. The existence of one perfect point for every mass does not make all such weighted families identical.

The uncertainty model and latent measure are declared inputs. Using them at low signal is exactly what the estimator was built to do. They must not suddenly become inadmissible because the central readings reach zero.

The resulting apparatus dependence is explicit. When the measurements carry little resolving power, the policy has more influence on the returned value. The estimate's standing records that influence.

## 5. Positive mass can remain the numerical output type

For finite, strictly positive Gaussian uncertainty scales and the current flat positive-magnitude measure, the compatible-pair distribution has finite, strictly positive first moments:

$$
0<E[f]<\infty,
\qquad
0<E[\alpha]<\infty.
$$

Therefore

$$
0<\widehat m=\frac{E[f]}{E[\alpha]}<\infty.
$$

This is the relevant totality claim. It is a property of this estimator on its declared empirical domain. Nonzero uncertainty by itself would not prove it for every possible distribution or policy.

The design contract is clear: an admissible uncertainty model and estimation policy must support the required finite positive readout. Explicitly exact inputs belong to their own checked operation. They do not undermine totality on the empirical domain.

## 6. The null result is continuous with the weak-signal result

Hold the uncertainty model fixed and let the force and acceleration readings approach zero. The current estimator approaches $\sigma_F/\sigma_a$ continuously.

That behaviour matches the gradual loss of resolving power. Object-specific signal gives way to apparatus and policy influence without a numerical singularity.

An undefined result introduced only at central zero–zero would puncture this continuous operation. A broader identification threshold would avoid dependence on rounding, but would still add a separate decision about when an existing estimate may cease to count as a mass answer.

I reject that additional interruption. The estimate should remain available, and its uncertainty and provenance should state how much it can support.

## 7. The Sun objection must respect the input and the output standing

The objection asks whether this design calls the mass of a nonaccelerating Sun $1\,\mathrm{kg}$. It combines a stronger description of the physical state with a weaker empirical input.

If absence of interactions has been established, the free-particle account handles the relevant motion without this null mass inference. If the Sun is interacting and its resultant is unresolved, the estimator receives an uncertain force-sum result. We cannot replace that input with knowledge that its underlying resultant is exactly zero.

If independent information already establishes the Sun's mass, that information must be retained in the larger account. A fresh, weak trial cannot silently supersede it. It may instead test compatibility with the established mass.

If no mass information is available beyond the null measurements, the current policy can indeed return $1\,\mathrm{kg}$ with apparatus-dominated standing. That number is a conditional estimate whose evidential limits must travel with it. Its availability does not license confident downstream prediction or unqualified mass reuse.

I accept the obligation this creates: the VD must preserve the standing of the positive scalar when it is consumed elsewhere. Hiding that standing would invalidate the design. Preserving it lets one numerical type cover strong, weak, and apparatus-dominated determinations without treating low information as a new kind of mass value.

## 8. Compatibility and numerical totality have different obligations

The estimator must still expose poor compatibility with the Newton-II account. A normalized compatible-pair calculation can return a positive number even when every plausible fit is bad.

I therefore preserve both a total numerical candidate and a separate compatibility verdict. A rejected candidate is not admitted for physical use. This proposal resolves the zero–zero issue for trials that pass the compatibility gate; it does not claim that every conceivable empirical input earns an accepted mass result.

At zero–zero the best-fit compatibility problem is absent: the origin fits. The contested question is whether weak identification should suppress the positive candidate. My answer is no. It should qualify the candidate's standing.

## 9. The uncertainty calculation supplies the behaviour we wanted

The objection to a $1\,\mathrm{kg}$ null estimate considered the point in isolation. I defend the complete operation: a positive scalar accompanied by an explicit account of what the trial resolves. The new uncertainty calculation shows how that operation can retain both numerical continuity and an honest information boundary.

At zero central readings, under independent isotropic Gaussian errors, a flat positive-magnitude measure, and no additional bounds or excitation weighting, write

$$
f=\sigma_F|Z_1|,
\qquad
\alpha=\sigma_a|Z_2|,
$$

where $Z_1$ and $Z_2$ are independent standard normal variables. Then

$$
M=\frac f\alpha
=s\frac{|Z_1|}{|Z_2|},
\qquad s=\frac{\sigma_F}{\sigma_a}.
$$

The induced mass density is

$$
p(m)=\frac{2s}{\pi(s^2+m^2)},
\qquad m>0.
$$

The positive scalar readout remains

$$
\widehat m=\frac{E[f]}{E[\alpha]}=s,
$$

and the median of $M$ is also $s$. Yet

$$
E[M]=\infty,
\qquad
E[(M-s)^2]=\infty.
$$

There is no finite mass standard deviation. Strictly, a conventional variance about a finite mean is unavailable because that mean does not exist. These statements concern the induced latent mass distribution; they do not assert that the sampling variance of the ratio-of-means estimator across repeated experiments is the same object.

This is a substantive answer to the fear that a finite readout must pretend to finite precision. The same declared calculation yields a finite positive summary and exposes the failure of a finite squared-error uncertainty description. The uncertainty is computed from the model, rather than added as an unsupported warning after the point has been chosen.

The fixed-mass compatibility calculation goes further. Because $Q(m)=0$ for every positive $m$, its allowed mass set is the entire positive ray. A complete null result can therefore retain:

```text
point estimate:                 sigma_F / sigma_a
fixed-mass compatibility set:   all positive masses
induced mass distribution:      half-Cauchy under the declared policy
mass standard deviation:        no finite value
standing:                      apparatus-dominated; no mass excluded by fit
```

That is the behaviour I want. The positive scalar remains available, while the result explicitly denies that this trial has narrowed the compatible mass range. The ratio-of-means readout and the compatibility set are different functionals of the declared model; the result names both rather than presenting one as the other.

### What the Sun example now tests

For an interacting Sun whose resultant and acceleration are unresolved, this complete result does not exclude a large solar mass on the basis of the fit. Its compatibility set includes that mass. The apparatus-scale point does not become a precise object fact merely by being stored in the mass field.

The half-Cauchy distribution nevertheless retains policy-dependent relative preferences. Its quantiles are

$$
m_p=s\tan\left(\frac{\pi p}{2}\right).
$$

For $s=1\,\mathrm{kg}$:

| Summary | Result |
|---|---|
| Ratio-of-means point and induced median | $1\,\mathrm{kg}$ |
| Central 68% induced interval | approximately $[0.257,3.895]\,\mathrm{kg}$ |
| Central 95% induced interval | approximately $[0.039,25.45]\,\mathrm{kg}$ |
| Fixed-mass compatibility set | $(0,\infty)$ |

I therefore do not defend the notation "$1\,\mathrm{kg}\pm\infty$" as a substitute for the full result. Infinite standard deviation is not complete ignorance, and the finite induced intervals are not fit-based exclusions. If those intervals were reported as though the null trial itself had excluded a large mass, the Sun objection would remain unanswered. The output must preserve the origin and meaning of each uncertainty statement.

### A general uncertainty rule, not a special null patch

The one-zero cases fit the same architecture. In the resolved-channel limits, $0/5$ can exclude large masses and $5/0$ can exclude small masses, while zero–zero excludes neither. The $5/0$ induced ratio can also lack finite moments because its latent denominator approaches zero. Thus a divergent moment is not a unique signature of null data, and need not become a new branching trigger.

I keep the positive readout throughout. The compatibility set, distribution, and provenance express the changing information content. A richer uncertainty output makes an additional undefined numerical branch unnecessary for the purpose of representing that change.

## Conclusion of Case II

**Measured zero–zero should return the positive scalar selected by the declared uncertainty-aware estimator.**

Established freeness is handled before this attempted mass inference. Within the interacting empirical account, a zero force sum is a measurement result with its own uncertainty. It does not assert the exact null state whose algebraic inversion is undefined.

The same operation that yields a small positive estimate at $0/5$ and a large positive estimate at $5/0$ yields $\sigma_F/\sigma_a$ at $0/0$. At the null input, the full calculation can also return an unrestricted mass compatibility set and an induced distribution with no finite standard deviation. That supplies the information boundary without deleting the positive scalar.

I retain that continuity, preserve the positive-scalar numerical type, and require uncertainty, compatibility, policy, and provenance to govern how the result may be used. The success of this design belongs to the complete output. A bare point, or a finite probability interval whose policy dependence has been hidden, does not satisfy its contract.

---

## Source basis and handoff notes

**Managed snapshot:** `canvalidk/VD-docs@58dba2d5c3d5d77802f3d57e0d1e2ba2f42b7369`, the `main` commit used for the preceding review. Inventory and draft status were read from `catalog.yaml` at that commit.

The arguments draw on the previously reviewed note workbench, particularly:

- `Newton-analysis/_note_why_mass_estimation/why_inertial_mass_must_be_estimated.md` — the managed current full draft; exact/empirical distinction and scope limits.
- `Newton-analysis/_note_why_mass_estimation/disagreement_rule_mass_estimator_trace-informing.md` — compatibility versus identification.
- `Newton-analysis/_note_why_mass_estimation/median_based_estimator_inertial_mass-informing.md` — latent-state weighting and the dependence of null inference on experimental assumptions.
- `Newton-analysis/_note_why_mass_estimation/inertial_mass_estimator_full_working_spec_d1-informing.md` — positive-magnitude mean readout, null cases, compatibility gate, and output standing.
- `Newton-analysis/_note_why_mass_estimation/what_forces_the_chain_to_the_estimator-informing.md` — exact branch, identification, and corrections to the stationary-root account.
- `Newton-analysis/vd_design_3_newton_i_correction.md` — question-dependent use of Newton I and Newton II; corrects the older universal free-particle routing proposal.

The informing copies identify their canonical homes in `Newton-analysis/_dscn_mass_estimator/`. They were used as the workbench inputs already reviewed, without creating a second maintained source collection.

**Active local authority consulted:** `08_nm_lesson_drafts/nml3/NML3_FORCE_SUM_TRIPLET_NAMING_DECISION_2026-08-26.md`, which supplies the current `net-force / impressed-force / interaction-set` names and force-sum relation.

**New contribution from this conversation:** the user proposed separating established freeness from the interacting mass-inference domain, then treating measured zero resultant force as an uncertain force-sum result. Case II develops that proposal in a committed voice. It is not recorded here as an already ratified replacement for the managed Newton I correction.

**Further contribution, 2026-09-05:** the user proposed that the mass uncertainty at zero–zero should be infinite, and then identified the combination of a positive point and unbounded uncertainty as potentially the desired estimator behaviour. Case II §9 develops that argument; Case I §8 answers it. The half-Cauchy density and quantiles are derived here from the stated Gaussian/flat-measure model. The general Cauchy distribution and its absence of finite mean and standard deviation are documented by [NIST](https://itl.nist.gov/div898/handbook/eda/section3/eda3663.htm). The positive half has divergent first and second raw moments. The unrestricted compatibility set follows independently from $Q(m)=0$ and must not be conflated with the half-Cauchy probability intervals.

**Boundary to preserve on ingestion:** removing a free null case from an attempted mass determination does not prohibit Newton II from calculating a free massive particle's acceleration using an already supplied mass. Neither case in this document establishes a blanket particle-routing rule.

**Outstanding decision exposed by the two cases:** whether insufficient identification suppresses an otherwise compatible positive candidate, or remains a qualification attached to an admitted estimate. Numerical totality and admission for downstream use must be stated separately whichever design is chosen.

**Uncertainty consequence for a later estimator specification:** declare separately the induced-ratio probability intervals, the fixed-mass compatibility set, and any sampling uncertainty of the point estimator. A generic field called `mass uncertainty` cannot silently stand for all three. Divergent moments must be reported as such rather than replaced by finite values caused by numerical integration cutoffs.
