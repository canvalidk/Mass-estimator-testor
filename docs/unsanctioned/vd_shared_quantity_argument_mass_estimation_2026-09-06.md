# Why the mass operation must accept potentially unaligned vectors

## The VD argument and its immediate consequences

**Date:** 2026-09-06  
**Status:** local explanatory draft and handoff, awaiting transfer to VD-docs. Develops the argument clarified by the user in this session and the previously ratified shared-quantity argument. The exposition and assembly below are proposed wording, not a new estimator-policy ruling.  
**Purpose:** provide a robust, reusable explanation for the introduction and the extended problem discussion. Establish the central dependency argument first, then connect the pieces that follow from it without introducing the proposed estimator prematurely.

## 1. The claim in one paragraph

The mass entry must accept potentially unaligned empirical force and acceleration results because those quantities must also be usable in the other two Newton-II operations. To determine force from mass and acceleration, acceleration must be obtainable without first knowing that force. To determine acceleration from force and mass, force must be obtainable without first knowing that acceleration. Aligning either vector to the other requires information from the other vector. Such alignment therefore cannot be a mandatory condition for obtaining the shared input quantity: it would make at least one operation require its own output before its input was available. Independently obtainable measurements can consequently reach the mass operation without agreeing in direction. The operation must address their compatibility where both are available; it cannot assume that their determination has already done this work for it.

## 2. The three expressions must be usable operations over the same quantities

For one material object A, over one relevant interval and in one inertial-frame-qualified account, Newton II relates

$$
\mathbf F[A]=m[A]\mathbf a[A],\qquad m[A]>0.
$$

The VD reads the three directions as operations with different demands:

$$
\begin{aligned}
\operatorname{Force}(m,\mathbf a)&=m\mathbf a,\\
\operatorname{Acceleration}(\mathbf F,m)&=\mathbf F/m,\\
\operatorname{Mass}(\mathbf F,\mathbf a)&=\text{the operation under investigation}.
\end{aligned}
$$

The first two are useful precisely because the demanded quantity may not already be available. We can have an independently supplied mass and acceleration and demand force. We can instead have an independently supplied mass and force and demand acceleration.

These are three operations over the same force, acceleration, and mass quantities. The acceleration supplied to the force operation must not acquire a new compulsory force-dependent meaning merely because acceleration is also an input to the mass operation. The corresponding requirement holds for force.

Here, 'the same quantity' does not mean that every later estimate must retain an earlier numerical value. Measurements can be refined and model-constrained estimates can be constructed. It means that such a change must not silently replace an independently available input with a differently conditioned object while claiming that the shared input definition has stayed the same.

The operative requirement is therefore:

> **Each direction must remain usable to determine its output from the other two quantities, including cases where that output has not otherwise been obtained.**

This is the VD premise doing the work. It is not merely a preference for detailed documentation.

## 3. Why acceleration cannot require prior knowledge of force

Consider a demand for net force. Mass is already supplied, and an acceleration-side procedure is available. The intended route is:

```text
acceleration-side evidence -> acceleration ----+
                                              +-> F = m a -> force
independently supplied mass ------------------+
```

Suppose we change what is required before the acceleration can be accepted. Its direction must first be reconciled with the object's net-force direction. Schematically,

$$
\mathbf a_{\mathrm{admitted}}
=\mathcal A(D_a,\mathbf F),
$$

where D_a is the acceleration-side evidence and the additional force input is needed to enforce alignment.

The demand now has the dependency:

```text
to determine force
    -> obtain an admitted acceleration
        -> obtain the force direction
            -> obtain the force being demanded
```

We no longer possess the input needed to run F = ma in the case for which we wanted that operation: force was not already available.

The defect is not that multiplying two known numbers would be impossible. It is that the proposed qualification of the input prevents those inputs from being available in a legitimate force-determination task.

Therefore, acceleration must have an admissible route that does not require the force being demanded. Prior alignment to that force cannot be a universal condition for qualifying as the acceleration input.

## 4. The symmetric requirement on force

Now demand acceleration from an independently supplied mass and a force-side result:

```text
force-side evidence -> net force -------------+
                                              +-> a = F / m -> acceleration
independently supplied mass ------------------+
```

If admitting the force requires aligning it with the object's acceleration, then supplying this operation's force input requires the acceleration it is supposed to produce.

Therefore, force must also have an admissible route that does not require that acceleration. Prior alignment to acceleration cannot be a universal condition for qualifying as the force input.

Notice that a reconciliation need not change both vectors to create the problem. If only acceleration is made force-dependent, the force-determination route is affected. If only force is made acceleration-dependent, the acceleration-determination route is affected. If both are made mutually dependent, both routes are affected. Disabling one of the required directions is already enough to defeat the proposed general repair.

## 5. What independent obtainability means

The conclusion is a requirement about necessary dependencies:

> **There must be a route to acceleration without the force whose determination it will support, and a route to force without the acceleration whose determination it will support.**

It cannot mean that force and acceleration may never be calculated from one another. That would prohibit the very Newton-II operations being defended. It means that dependence on the other vector cannot be compulsory for every admissible instance of the shared input quantity.

Nor does it mean that the channel errors must be statistically independent. Shared calibration, frame information, or instrument effects may create correlations. Such correlations can be represented in the uncertainty model without requiring that a procedure first know the output it is meant to help determine.

Independent routes also need not exist for every object in every possible experimental circumstance. A particular task may lack the required evidence. The claim is that the definitions must preserve the routes when their evidence is available; they cannot make those routes unavailable by imposing prior mutual alignment on the input quantities.

Finally, an external geometric restriction may supply a direction without consulting the unknown other vector. That does not create the dependency above. A specially constrained experiment is therefore compatible with the argument, but it does not establish a general guarantee of alignment for all empirical mass inputs.

## 6. What this forces the mass entry to receive

The mass operation consumes the same independently obtainable quantities. For empirical work these arrive as measurement results, schematically

$$
M_F=(\widetilde{\mathbf F},U_F,P_F),
\qquad
M_a=(\widetilde{\mathbf a},U_a,P_a),
$$

where U denotes uncertainty information and P denotes provenance. A joint uncertainty description may additionally record cross-channel correlations.

The procedures producing these results need not preserve exact codirectionality of their central vectors. Neither one's admission as a measurement result can generally require inspecting and aligning to the other vector first. The mass operation must therefore be able to receive an unaligned pair as a pair of legitimate empirical results.

This is the crucial connection:

> **The other two corners determine what the mass corner is allowed to assume about its inputs. They require independently obtainable quantities; those quantities do not carry a general certificate of mutual alignment.**

'Accept' here means accept for assessment as empirical inputs to the determination. It does not mean that every pair must yield an admitted positive mass, that poor provenance must be ignored, or that an incompatible exact pair must be silently repaired. The operation can discover a failed model fit or insufficient information. It must handle the situation as part of the demand rather than pretend the input quantities could only have existed after reconciliation.

Under a continuous error model with a full-dimensional joint density, exact nonzero codirectionality in two or three dimensions occupies a lower-dimensional set and has probability zero. This explains why the problem is ordinary in unconstrained vector measurement. The dependency argument does not need that probability-zero claim to be universal: the possibility of valid unaligned empirical inputs is enough to defeat a general admission requirement of prior alignment. Rounded displayed alignment also does not change a measurement result into an exact premise.

## 7. Why 'first align the inputs' does not preserve the proposed definition

There are three different proposals that this phrase can conceal.

### 7.1 Redefine the shared inputs as already aligned

The force and acceleration delivered as the shared quantities have been adjusted toward one another, and that reconciliation is made a condition of obtaining them. This is the proposal ruled out by sections 3-4: at least one input now requires the output of another Newton-II operation.

The response 'metrology will do it' does not alter that dependency. An alignment performed by metrology needs the same information as alignment performed by any other agent.

### 7.2 Keep the independently obtained quantities and require the mass entry to receive only compatible ones

A restricted exact-recovery helper is mathematically legitimate. On a certified compatible, nonzero pair, it can return the magnitude quotient.

But declaring that helper's domain does not discharge the empirical mass demand. The independently obtainable measurements established above can still be unaligned. A complete account must explain how that demand proceeds from them, or why it terminates with a specified non-mass outcome. The restricted helper alone is not the complete empirical definition.

### 7.3 Preserve the inputs and construct a separate compatible latent pair

This is coherent. Both measurement results are available, and a procedure constructs or weights distinct latent quantities satisfying

$$
\mathbf F^*=m\mathbf a^*,\qquad m>0.
$$

The starred quantities are outputs of joint inference. They do not retroactively become what either independent measurement had to provide before entering the triplet.

This proposal establishes the required conclusion: the mass determination includes the joint construction. A final scalar quotient may be evaluated afterward, but it is only a component of the composed operation.

The numerical routine can be external, run first within a software pipeline, or return an intermediate object. None of those implementation choices removes it from the operation that starts with the independent measurements and answers the mass demand. 'Internal to mass determination' concerns that dependency and scope; it does not require the dictionary to contain the estimator's algorithm.

## 8. The immediate consequence for magnitude division

The opening hook is that vector division is not a general operation. The proposed repair is to use magnitudes. On a compatible nonzero exact pair, that repair works:

$$
\mathbf F=m\mathbf a,\quad m>0,\quad\mathbf a\ne\mathbf0
\quad\Longrightarrow\quad
m=\frac{\|\mathbf F\|}{\|\mathbf a\|}.
$$

The quotient reads the scalar after the relation has been established. It does not establish the relation.

For example, take

$$
\widetilde{\mathbf F}=(10,0)\ \mathrm N,
\qquad
\widetilde{\mathbf a}=(0,5)\ \mathrm{m\,s^{-2}}.
$$

The magnitude quotient is 2 kg, but

$$
(2\ \mathrm{kg})\widetilde{\mathbf a}=(0,10)\ \mathrm N
\ne\widetilde{\mathbf F}.
$$

The number has the correct dimensions and is positive, yet it does not reconcile the supplied vectors. If the inputs are exact, no positive mass satisfies them. If they are measured, the uncertainty account is needed to assess whether a compatible underlying state is plausible.

Taken alone, this example could be dismissed as use of the quotient outside its exact domain. The VD argument supplies the missing force: independently obtained empirical inputs of this kind cannot be excluded by defining the shared force and acceleration quantities as already mutually reconciled. The complete mass operation has work to do that the quotient does not perform.

This does not show that a declared estimator can never return the numerical magnitude ratio in a particular case. It shows that the bare quotient cannot, by itself, constitute the general empirical mass determination. Equal numerical outputs in a special case do not make two operations identical.

## 9. The nearby exact-input result: none, one, or every positive mass

For exact supplied vectors, define the compatibility set

$$
C(\mathbf F,\mathbf a)
=\{m>0:\mathbf F=m\mathbf a\}.
$$

The cases are exhaustive:

| Exact pair | Compatible positive masses | What the mass demand establishes |
|---|---|---|
| Nonzero and codirectional | Exactly one | Recover the unique scalar by a magnitude quotient. |
| Both zero | Every positive mass | The pair agrees with the law but does not identify mass. |
| Exactly one zero, or nonzero and not codirectional | None | The supplied exact premises do not satisfy finite positive-mass Newton II. |

The algebraic asymmetry follows immediately. Multiplying any acceleration vector by a positive mass gives a force vector. Dividing any force vector by a positive mass gives an acceleration vector. Two arbitrary vectors, however, need not identify one positive scalar relating them.

This is a statement about exact algebraic input types. It does not claim that experimental force and acceleration determinations are free of estimation or uncertainty propagation. The additional issue in the mass corner is the relational compatibility and identification condition on its two vector inputs.

The phrase 'undefined' alone hides the difference between no compatible mass and every mass being compatible. The trace needs to preserve which situation has occurred.

## 10. Why the zero cases strengthen the problem

The zero cases reveal a difficulty that remains even after all directions have been reduced to a common scalar axis.

For an exact zero-zero pair, multiplication by zero removes all dependence on m:

$$
0=m0\qquad\text{for every }m>0.
$$

An inverse procedure cannot identify a distinction its inputs do not contain. Moreover, exact pairs (F,a)=(kt,t) approach the origin with any chosen positive mass k. No one number assigned at exact zero-zero continuously completes all those recoveries.

The expectation that mass should provide a positive scalar adds operational pressure. A zero quotient violates that promised value range, and an undefined result does not directly supply the input needed by a later Newton-II calculation. But the physical mass type need not contain `undefined`: a determination can instead return a non-identification or failure result. The unresolved task is to specify the determination's output and what can legitimately be done with it.

For empirical results, a zero central reading is not automatically an exact physical zero. A zero-force reading and a resolved nonzero acceleration can carry different information from an exact force of zero. The same applies to a zero-acceleration reading and to two unresolved channels. Their uncertainty models are part of the input that the mass operation must assess.

This establishes why the zero readings require a measurement-aware operation. It does not yet select a finite sigma/sigma readout, require an undefined empirical-null result, or settle downstream admission. Those stronger choices belong to the later estimator discussion. Exact non-identification and empirical inference must not be merged by an arithmetic shorthand such as 0/0.

## 11. Estimation is located, but its policy is not selected by this argument

Once both independent measurement results are in scope, a mass determination can infer a compatible latent relation while preserving the independent inputs as its evidence. That inference can construct a fitted pair or a distribution over compatible pairs; the dependency argument does not select between those representations.

Newton II alone does not assign the relative plausibility of different changes to the reported vectors. For nonzero measured vectors, several scalar candidate expressions can agree on exact codirectional pairs and disagree elsewhere. For example, using

$$
\widetilde{\mathbf F}=(10,0)\ \mathrm N,
\qquad
\widetilde{\mathbf a}=(4,3)\ \mathrm{m\,s^{-2}},
$$

gives

$$
\frac{\|\widetilde{\mathbf F}\|}{\|\widetilde{\mathbf a}\|}=2\ \mathrm{kg},
\qquad
\frac{\widetilde{\mathbf a}\cdot\widetilde{\mathbf F}}
     {\widetilde{\mathbf a}\cdot\widetilde{\mathbf a}}=1.6\ \mathrm{kg},
\qquad
\frac{\widetilde{\mathbf F}\cdot\widetilde{\mathbf F}}
     {\widetilde{\mathbf F}\cdot\widetilde{\mathbf a}}=2.5\ \mathrm{kg}.
$$

Each expression reduces to the exact scalar on the positive nonzero law relation. Agreement there does not select their behavior away from it. These are illustrative competing readouts, not a recommendation to use any of them without a measurement policy.

The immediate conclusion is that the complete empirical operation must supply and justify the reconciliation or inference rule, including the uncertainty assumptions relevant to the task. Stronger composition requirements may constrain the eventual readout, as the separate readout work investigates. That is a further argument; it is not contained in the dependency proof.

For the dictionary, the required operation can be named without teaching or owning its numerical implementation. The source ruling places estimator competence outside the dictionary's authority while requiring its invocation within mass determination.

## 12. The adjacent evidential requirement: preserve what was changed

A latent pair constructed subject to F* = m a* satisfies the imposed relation by construction. That equality alone cannot tell us how well the independent measurements supported it.

If the purpose of the returned trace is to preserve the grounds of the result, the output must retain enough of the original measurements, uncertainty, provenance, and discrepancy to explain its standing. A distribution normalized over compatible states similarly expresses relative support conditional on that construction; normalization alone does not establish absolute compatibility with the observations.

This motivates a separate fit assessment and the distinction between a conditional numerical estimate and a result licensed for use.

The extra premise here is evidential accountability. It is closely adjacent to the VD dependency argument, but it is not the same proof. The dependency argument locates the joint operation and prevents compulsory cross-dependent input definitions. The evidential argument explains why the discrepancy must survive the operation. Keeping them distinct prevents every output-field requirement from being mistaken for part of the core definition.

## 13. A compact passage for the introduction

The following passage assembles the central argument with the magnitude-division transition. It is reusable draft prose, not a replacement for the full introduction.

> Newton II must allow us to determine force when force is not already known. The acceleration used in that operation must therefore be obtainable without first consulting the force. Otherwise, obtaining the input would require the output we were trying to calculate. The converse direction imposes the symmetric requirement: force must be obtainable without the acceleration that it will be used to determine.
>
> This matters when we turn to mass. The force and acceleration used there are the same quantities used in the other two operations. We cannot require them to have been aligned with each other before they qualify as inputs, because alignment of either to the other would introduce precisely the dependency that one of the other operations cannot require. Their independent measurement results may consequently arrive pointing in different directions.
>
> A quotient of magnitudes does not resolve that situation. It produces a positive scalar even for perpendicular vectors, although no scalar multiplication turns one of those vectors into the other. We also cannot dismiss the difficulty by asking for the 'properly aligned' quantities instead: we must explain the joint operation that constructs them from the independent results. That construction is already part of determining mass.
>
> The quotient is therefore available after a suitable relation has been established, but it does not describe the full operation that must establish or infer that relation. The zero cases expose further questions about whether a scalar is identified at all. What, then, does determining inertial mass consist in?

## 14. How the pieces fit, and where this draft stops

The assembled chain is:

```text
three usable Newton-II directions over the same quantities
    -> independent input routes must remain available
    -> prior alignment to the other vector cannot be compulsory
    -> empirical mass inputs may be unaligned
    -> a bare magnitude quotient does not perform the required assessment
    -> joint inference belongs within the mass determination

exact scalar recovery and the zero cases
    -> distinguish incompatibility from non-identification
    -> preserve the exact/empirical input distinction

an accountable empirical result
    -> retain evidence and discrepancy
    -> distinguish a computed point from what its evidence warrants
```

The first line is a VD framework commitment. The input-dependency consequences follow from it. Potential empirical misalignment additionally uses the character of the measurement model. The exact-case table is algebra. Discrepancy retention uses the additional evidential premise. Selection of the full estimator, its latent measure, the empirical-null policy, and its numerical realization remain subsequent work.

This document deliberately makes no claim that practical mass measurement universally uses a particular sigma formula, and no claim of historical originality. Those require their own evidence. The argument here explains why the quotient cannot remain the complete empirical definition under the stated shared-quantity commitments.

## 15. Provenance, relation to earlier records, and submission

**Managed snapshot:** canvalidk/VD-docs@743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6. main was resolved again for this task on 2026-09-06 and matched the preceding discussion. Inventory/status came from catalog.yaml at that snapshot.

Task sources at that commit:

- [Why reconciliation cannot be upstream of E10](https://github.com/canvalidk/VD-docs/blob/743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6/Newton-analysis/_dscn_mass_estimator/why_reconciliation_cannot_be_upstream_of_E10.md), especially the three operations, the input dependencies, and ratification status. The shared-quantity argument is recorded there as ratified on 2026-08-01.
- [What forces the chain to the estimator](https://github.com/canvalidk/VD-docs/blob/743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6/Newton-analysis/_dscn_mass_estimator/what_forces_the_chain_to_the_estimator.md), especially section 4's clarification that inoperability concerns a demand whose output is unavailable, and its treatment of a distinct latent pair.
- [Why inertial mass must be estimated](https://github.com/canvalidk/VD-docs/blob/743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6/Newton-analysis/_note_why_mass_estimation/why_inertial_mass_must_be_estimated.md), catalog-current d3, for exact recovery, measured versus latent quantities, the reconciliation argument, and its stated scope.
- [Full working estimator specification](https://github.com/canvalidk/VD-docs/blob/743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6/Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md), for the current separation of compatible-state inference, absolute compatibility, output standing, and exact/empirical branches.

**Active local context:** `mass_estimation_introduction_development_plan_2026-09-06.md` in this same handoff folder, especially stages 2-5. The user asked for a more robust standalone explanation of the VD argument and its immediate adjacent consequences before drafting the full introduction.

**What this exposition sharpens:** mandatory versus optional dependencies; the role of an unavailable demanded output; the fact that changing only one vector can already disable one corner; the distinction between independent obtainability and statistical independence; the compatibility of a separately labeled latent construction with preserving the original inputs; and the difference between accepting an empirical input for assessment and promising a successful mass output.

These sharpen the explanation of the ratified argument. They are not a new ruling to edit managed sources, change estimator policy, or revise active entries. In particular, the older source's broad wording about neither quantity being defined through the other is read as a prohibition on compulsory input dependencies, not as a prohibition on Newton II itself providing one of the available determination routes.

**Coverage:** focused source reading and new explanatory synthesis, not an exhaustive historical or inbox audit and not a new external literature search. The quantitative examples and exact-case classification are elementary calculations stated in the document.

**Submission:** saved in the workspace-root `pass to VD-docs/` folder. Awaiting transfer; not submitted to VD-docs. No managed document, existing introduction, entry packet, engine, or evaluator file was modified.
