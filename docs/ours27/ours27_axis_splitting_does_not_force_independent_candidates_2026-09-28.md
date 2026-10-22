# Axis splitting does not force independent candidate weighting

**Record — 28 September 2026.** Written by Codex at the user's request.

**Finding.** Representing a vector reading as a list of component readings does not, by itself, require independent candidate weighting on those components. The argument from axis splitting to flat candidate measure needs an additional combination rule: eliminate each component's candidate state independently, then multiply the resulting mass evidence. A shared candidate-size or covariance state can instead be retained through the change of representation.

This is a correction to the scope of a premise argument. It does not adopt a new estimator, prove a unique replacement measure, or discard the September 25 theorem under its stronger assumptions.

## 1. The passage that needs qualification

The [September 25 axis-splitting record](https://github.com/canvalidk/VD-docs/blob/b4bfc2040031ef010416334d0c5f4aa2ad05cbd8/Newton-analysis/_dscn_mass_estimator/axis_splitting_and_the_24_september_change_2026-09-25.md) begins with an elementary bookkeeping requirement: a vector reading and the corresponding component readings should give the same mass law. It then specifies combination by separately integrating each reading's latent pair and multiplying the resulting evidence, with the mass starting law applied once.

Section 4(b) uses that specification to require a product candidate weighting across axes. Isotropy then sharply restricts product laws; the further treatment of size scale leads to the flat weighting used in ours27.

Two requirements must therefore be distinguished:

1. **Representation invariance:** the same joint experiment gives the same answer when its components are relisted or regrouped.
2. **Independent nuisance elimination:** each component has a candidate state that can be integrated out independently before the component evidence is combined.

The second entails a particular factorization. It does not follow from the first.

The ratified [shared-quantity argument](https://github.com/canvalidk/VD-docs/blob/b4bfc2040031ef010416334d0c5f4aa2ad05cbd8/Newton-analysis/_dscn_mass_estimator/why_reconciliation_cannot_be_upstream_of_E10.md) requires force and acceleration to be independently obtainable and places their reconciliation inside mass determination. It does not assert probabilistic independence of candidate components. Newton II holding component by component also does not establish that independence.

## 2. A proper counterexample

Let a candidate vector $v\in\mathbb R^d$, with $d\geq2$, have a shared variance state $\tau$:

$$
\Pr(\tau=1)=\Pr(\tau=4)=\frac12,
\qquad
v\mid\tau\sim N(0,\tau I_d).
$$

This is a proper, isotropic candidate law. Conditional on $\tau$, the components are independent. Marginally, they are not:

$$
\operatorname{Cov}(v_1^2,v_2^2)
=\operatorname{Var}(\tau)
=\frac94>0.
$$

Suppose the measurement likelihood factors over axes conditional on the candidate, with component likelihoods $L_k(m,v_k)$. Define

$$
g_k(m,\tau)
=\int L_k(m,v_k)\,\phi(v_k;0,\tau)\,dv_k,
$$

where $\phi$ denotes a Gaussian density. The vector calculation gives

$$
L_{\mathrm{joint}}(m)
=\frac12\prod_{k=1}^d g_k(m,1)
 +\frac12\prod_{k=1}^d g_k(m,4).
$$

Splitting the axes and retaining the same $\tau$ gives exactly the same expression. Thus the calculation is invariant under this representation change.

In contrast, integrating a fresh $\tau$ separately on every axis gives

$$
L_{\mathrm{reset}}(m)
=\prod_{k=1}^d
\left\{\frac12g_k(m,1)+\frac12g_k(m,4)\right\}.
$$

For two axes, the difference is

$$
L_{\mathrm{joint}}-L_{\mathrm{reset}}
=\frac14
\bigl[g_1(m,1)-g_1(m,4)\bigr]
\bigl[g_2(m,1)-g_2(m,4)\bigr].
$$

It is generally nonzero and mass dependent. The reset calculation describes independent variance states $\tau_1,\tau_2$, rather than a different notation for one shared state.

This counterexample uses a proper finite mixture so that the issue cannot be attributed to an improper integral or an interchange of divergent operations. It refutes the implication from representation invariance alone to product candidate weighting. Its chosen scales are illustrative; it is not offered as an instrument-only default or as a solution of the no-size-scale requirement.

## 3. What must travel with a coordinate change

For independent measurement errors, multiplication of component likelihoods is legitimate **conditional on all shared states**:

$$
L(m)
=\int \nu(d\tau)\prod_k g_k(m,\tau).
$$

Integrating before or after regrouping this same joint integral changes nothing. Replacing it by a product of separately averaged factors changes the dependence structure.

With correlated errors, the full joint measurement likelihood must also travel with the recoding. The [full-covariance ours27 record](https://github.com/canvalidk/VD-docs/blob/b4bfc2040031ef010416334d0c5f4aa2ad05cbd8/Newton-analysis/_dscn_mass_estimator/ours27/ours27_full_covariance_2026-09-27.md) already allows spatial and force–acceleration error correlations and notes the need to restate the diagonal-axis version of L8.

For example, suppose $N$ candidate vectors share a covariance tensor $\Omega$. Their stacked conditional covariance is $I_N\otimes\Omega$. Under an invertible linear recoding $T$, it becomes

$$
T(I_N\otimes\Omega)T^\top.
$$

The data, measurement covariance, candidate measure and sharing rule must all be pushed forward together. A transformed model need not retain the same simple block notation. Fitting new independent scalar models to its components is an additional model choice.

This also distinguishes two claims about repeated readings. A shared mass does not by itself imply that readings repeat one candidate vector, or that their different candidate vectors are independent, or that they share one size tensor. Those are separate specifications. The present result does not close the project's persistence question.

## 4. Consequence for the search

A proposed clarification of axis splitting is:

> Relisting or regrouping components preserves the mass law when the full joint experiment, candidate measure, covariance and shared states are preserved.

Under this clarification, candidate-size or candidate-direction learning is compatible with elementary axis bookkeeping. The stronger rule requiring independent nuisance elimination remains available, but it must be stated and justified independently. Its consequences cannot be attributed to representation invariance alone.

The conditional isotropy/factorization argument in the September 25 record therefore survives within its stated combination model. What remains open is whether that model is forced by the project's requirements.

The [zero-reading record](ours27_zero_readings_and_shared_candidate_information_2026-09-28.md) gives a concrete reason to retain shared state. The [repair-limits record](ours27_limits_of_scalar_and_tensor_size_repairs_2026-09-28.md) explains why reopening that possibility still does not supply a final equation.

## 5. Provenance and status

Managed sources were selected through the VD-docs catalog and read at the single snapshot **canvalidk/VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8**. The links above identify the relevant canonical sources. The associated [September 27 lemmas](https://github.com/canvalidk/VD-docs/blob/b4bfc2040031ef010416334d0c5f4aa2ad05cbd8/Newton-analysis/_dscn_mass_estimator/mass_estimator_lemmas_2026-09-27.md) provide the L8–L10 context.

The counterexample and distinction in this record are direct calculations from the September 28 investigation. Its local premise audit is in .tools/ours27_direction_2026_09_28/premises/direction_and_null_premises.md. This record contains the substantive argument and does not require that scratch file.

This is a local record for the user to carry. It has not been submitted to VD-docs, and no managed premise or production equation has been changed.
