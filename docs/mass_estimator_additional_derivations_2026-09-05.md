# Additional derivations and constraints for the mass readout

**Date:** 2026-09-05  
**Status:** local active, unsubmitted mathematical handoff; proposed interpretations and diagnostics, not ratified estimator policy.  
**Purpose:** identify useful properties and derivations beyond the preceding originality assessment. “Additional” means additional to the current inspected argument, not established historical originality.

## Scope and notation

Let $P$ be the already constructed normalized distribution over Newton-compatible positive magnitudes $f,\alpha$ and their common direction. Write

$$
B=E_P[f],\qquad A=E_P[\alpha],\qquad C=E_P[\sqrt{f\alpha}],
\qquad \widehat m=B/A.
$$

Assume $f,\alpha>0$ almost surely and $0<A,B<\infty$. The default Gaussian/flat-positive-magnitude construction satisfies these assumptions. Most results below do not otherwise require Gaussian errors, independent magnitudes, or a fixed latent direction. They concern the readout applied to the resulting magnitude distribution, not a replacement for reconciliation or the absolute compatibility gate.

## 1. A finite, reciprocal-symmetric loss selects exactly this readout

For a candidate mass $m>0$, define

$$
\ell(m;f,\alpha)
=\left(\frac{\sqrt f}{m^{1/4}}-\sqrt\alpha\,m^{1/4}\right)^2
=\frac f{\sqrt m}+\alpha\sqrt m-2\sqrt{f\alpha}.
$$

This is a squared mismatch between positive square-root magnitudes, with the mass conversion split symmetrically between the two channels.

It has several exact properties:

- It is nonnegative, and vanishes precisely when $f=m\alpha$.
- Exchanging $f$ with $\alpha$ and $m$ with $1/m$ leaves it unchanged.
- Under independent unit rescalings $f'=c_Ff$, $\alpha'=c_a\alpha$, and $m'=(c_F/c_a)m$, it multiplies by the candidate-independent factor $\sqrt{c_Fc_a}$. The selected physical mass therefore does not depend on unit conventions.
- Its expectation is finite for every finite positive candidate mass. Indeed, $C\le\sqrt{AB}$ by Cauchy–Schwarz, and $0\le\ell\le f/\sqrt m+\alpha\sqrt m$.

The expected loss is

$$
\mathcal R(m)=\frac B{\sqrt m}+A\sqrt m-2C,
\qquad
\mathcal R'(m)=\frac{Am-B}{2m^{3/2}}.
$$

It decreases for $m<B/A$ and increases for $m>B/A$. Therefore its unique global minimizer is

$$
\boxed{\arg\min_{m>0}E_P[\ell(m;f,\alpha)]=\frac{E_P[f]}{E_P[\alpha]}.}
$$

Equivalently, the risk is strictly convex in $t=\log(m/m_0)$ for any fixed positive reference mass $m_0$. Relative to the optimum,

$$
\mathcal R(\widehat m e^t)
=2\sqrt{AB}\cosh(t/2)-2C.
$$

Equal multiplicative errors above and below the chosen point have equal excess loss. This is a symmetry in relative error, not symmetry under equal additive kilogram errors.

**Why it helps:** we now have a direct decision-theoretic derivation that simultaneously respects the exact positive mass relation, channel reciprocity, units, and finite expected loss at the zeros. In particular, it avoids the integrability failure of trying to minimize $E[\alpha(M-m)^2]$ when $M=f/\alpha$ has a boundary-induced divergent mean. It does not require $E_P[M]$ to exist.

**What has and has not been established:** the optimizer theorem is exact under the stated assumptions. The proposed loss is an additional defensible criterion; Newton II and reciprocity have not been proved to force this loss uniquely. This derives the readout conditional on $P$, not the Gaussian model or the latent measure used to construct $P$.

The square-root construction is related to established Hellinger geometry. A separate literature follow-up also located loss-function characterizations of ratios of expected outcomes in another field; therefore neither ratio-selecting losses in general nor square-root geometry should be described as newly invented here. The precise mass/inverse-mass application above is our deduction in this note. [Gorczyca and Kang, 2024, section 2.3](https://onlinelibrary.wiley.com/doi/10.1002/sta4.680), [Hellinger distance and affinity definitions](https://pmc.ncbi.nlm.nih.gov/articles/PMC10137612/).

## 2. Improving absolute precision cannot manufacture mass precision at null

Under the default independent isotropic Gaussian model and flat positive-magnitude measure, zero–zero gives

$$
f=\sigma_F U,\qquad \alpha=\sigma_a V,
\qquad U,V\ \text{independent standard half-normal variables}.
$$

If both uncertainty scales shrink by the same factor $\lambda>0$, the new magnitude law is the old law with $(f,\alpha)$ replaced by $(\lambda f,\lambda\alpha)$. Consequently,

$$
\widehat m_\lambda=\widehat m,
\qquad
M_\lambda\overset d=M,
\qquad M=f/\alpha.
$$

The entire induced mass distribution is unchanged, including every quantile and every interval obtained from those quantiles. Meanwhile, the absolute latent magnitudes become smaller. This gives a useful behavioral constraint: **increasing precision around an unexcited null state must not be mistaken for identifying its mass ratio.** The default estimator passes this test.

One repeated-measurement specialization is exact. Suppose independent Gaussian measurements all concern the same latent force–acceleration pair, with known unchanged channel variances, and their sufficient sample means are both zero. After $n$ measurements, the effective scales are $\sigma_F/\sqrt n$ and $\sigma_a/\sqrt n$. The induced null mass distribution is the same as for one such observation. There is no artificial narrowing merely from repeating that null determination.

This is not a statement that all repeated experiments leave mass uncertainty unchanged. New nonzero excitation can identify mass. Trials with different latent accelerations but a common mass require a different joint model. Unequal improvements to the two channel scales change the apparatus ratio and can change the policy-induced mass distribution. Finally, the limit as $\lambda\to0$ does not turn an exact physical $0/0$ premise into an algebraically identified ratio; the boundary remains singular.

## 3. Grouping the same uncertainty into cases cannot change the answer

Partition the same law $P$ into cases $K=k$, with probabilities $p_k$. Define

$$
A_k=E[\alpha\mid K=k],\quad B_k=E[f\mid K=k],
\quad m_k=B_k/A_k.
$$

By total expectation,

$$
\widehat m
=\frac{\sum_k p_kB_k}{\sum_k p_kA_k}
=\sum_k w_km_k,
\qquad
w_k=\frac{p_kA_k}{\sum_j p_jA_j}.
$$

Thus the combined readout is an acceleration-weighted average of the conditional readouts. Inverse-mass readouts combine analogously with force weights. This is the conditional version of the size-bias interpretation already found in the research. [Coordinate size biasing, Arratia, Goldstein and Kochman, section 2.3](https://arxiv.org/pdf/1308.2729).

Consequences include:

- The combined point lies between the smallest and largest conditional points.
- If every case gives the same mass, forgetting which case occurred preserves that mass.
- Splitting a case into subcases and recombining them preserves the answer, provided probability and first-moment information are retained.
- Aggregation can be implemented associatively by adding the unnormalized triples $(p_k,p_kA_k,p_kB_k)$.

For example, two equally probable exact pairs $(f,\alpha)=(1,1)$ and $(100,10)$ have masses 1 and 10 in consistent units. Their combined readout is $101/11\approx9.18$, with mass weights $1/11$ and $10/11$. Averaging the two masses with equal weights would lose the required magnitude information.

**Why it helps:** the output is insensitive to arbitrary bookkeeping choices about how the same uncertainty is grouped. It also tells us exactly what information must be preserved when a calculation is summarized. Retaining only conditional mass points and their probabilities is insufficient.

This is an aggregation theorem about one probability law. Combining new experiments by averaging separately fitted mass estimates is not automatically licensed: shared parameters and measurement dependence must first be modeled.

## 4. A finite, symmetric diagnostic emerges from the same derivation

Define

$$
\kappa=\frac{E_P[\sqrt{f\alpha}]}{\sqrt{E_P[f]E_P[\alpha]}},
\qquad D=1-\kappa.
$$

Cauchy–Schwarz gives $0<\kappa\le1$, so $0\le D<1$. Equality $D=0$ occurs exactly when $f$ is a fixed positive multiple of $\alpha$ almost surely. In other words, every supported latent pair agrees on one mass. This can happen even when the individual force and acceleration magnitudes remain uncertain.

There is also a familiar probability interpretation. Define the two size-biased distributions

$$
dQ_f=\frac fB\,dP,\qquad dQ_\alpha=\frac\alpha A\,dP.
$$

Their Hellinger affinity is $\kappa$, and their squared Hellinger distance, under the convention $H^2=1-\int\sqrt{dQ_fdQ_\alpha}$, is $D$. Thus $D$ measures how differently the force and acceleration weight the same compatible states. The loss minimum satisfies

$$
\min_m\mathcal R(m)=2\sqrt{AB}\,D.
$$

This is useful because it remains finite when the ordinary mass-ratio variance does not. It is invariant under units and channel exchange. At the default null law,

$$
\kappa_0=\frac{\Gamma(3/4)^2}{\sqrt\pi}\approx0.847213,
\qquad D_0\approx0.152787.
$$

At support concentrated on one mass, $D=0$. It therefore separates the exact null law from a law concentrated on one mass without appealing to divergent moments.

**Diagnostic limit:** this is a discrepancy between two size-biased laws, not the complete spread of $M$ under $P$, not an absolute model-fit test, and not a calibrated admission threshold. States with very small magnitudes can be downweighted in both biased laws. A concentrated compatible-state distribution can also coexist with poor absolute compatibility with the observations. Keep the original compatibility gate and the probability/compatibility distinctions.

## 5. Logarithmic mass uncertainty is well behaved at null

With $s=\sigma_F/\sigma_a$, the null induced ratio is half-Cauchy. Set $L=\log(M/s)$, which is dimensionless. Changing variables from $M$ gives

$$
p_L(\ell)=\frac1{\pi\cosh\ell},\qquad
E[L]=0,\qquad \operatorname{Var}(L)=\frac{\pi^2}{4}.
$$

For example, the moment-generating integral is $E[e^{tL}]=\sec(\pi t/2)$ for $|t|<1$; its first two derivatives at zero give the displayed moments. The log-scale standard deviation is therefore $\pi/2$. It is a finite, reciprocal-symmetric measure of multiplicative dispersion even though the ordinary positive mass mean diverges.

Inverting mass changes the sign of the centered logarithm and preserves its variance. Under the default finite-data Gaussian magnitude model, every fixed-order absolute log-ratio moment is finite: logarithmic singularities at the positive boundaries are integrable, and Gaussian tails control large magnitudes.

This offers an additional descriptive field without changing the declared point estimator. A log-scale standard deviation is not a confidence interval, and it does not replace the unrestricted fixed-mass compatibility set at zero–zero. Away from null, the geometric mean of the ratio need not equal the ratio-of-means point.

## What these results add

The square-root loss is the strongest additional derivation: the exact readout is the unique optimum of a finite criterion respecting reciprocal and dimensional symmetry. The null scaling theorem supplies a direct test against false precision. Conditional aggregation supplies a representation-consistency test. Hellinger and log-scale summaries offer finite diagnostics where ordinary ratio moments fail.

Together they strengthen the mathematical case for this readout. They do not yet prove uniqueness among all physically acceptable criteria or settle the latent-measure and downstream-admission policies.

## Verification and source boundary

The algebra above was checked directly. Numerical checks over 1,000 randomly generated positive discrete joint laws confirmed the loss identity and minimum, channel exchange, independent unit rescaling, mixture aggregation, and Hellinger identity. Independent quadrature checked the null affinity and log-ratio variance. These checks supplement the analytical proofs; they do not calibrate sampling coverage or any admission rule.

**Managed snapshot:** `canvalidk/VD-docs@58dba2d5c3d5d77802f3d57e0d1e2ba2f42b7369`. `main` was resolved again for this turn and matched the preceding research snapshot. `catalog.yaml` and the already-read canonical sources at this commit supplied the current comparison target, especially:

- `Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`
- `Newton-analysis/_dscn_mass_estimator/what_forces_the_ratio_of_means_readout.md`

The inspected specification and readout audit did not contain the square-root-loss derivation, Hellinger diagnostic, or full null scaling formulation above. This is a scoped comparison, not an exhaustive claim about every managed or historical note. The local `mass_estimator_originality_research_2026-09-05.md` and `mass_estimation_zero_zero_two_cases.md` in this handoff folder supplied the active conversation context. No separate inbox or historical collection was exhaustively searched.

**External follow-up:** the Hellinger definitions were checked in *Empirical Squared Hellinger Distance Estimator and Generalizations to a Family of α-Divergence Estimators*; size biasing was checked against Arratia et al. Gorczyca and Kang's 2024 paper supplies an additional neighboring example of losses selecting ratios of expectations, not this exact Newton construction. Access date: 2026-09-05. No historical originality claim is made for the deductions in this note.

**Submission:** saving this new record under `pass to VD-docs/` does not itself submit it to the managed collection. Existing managed and local argument documents were not overwritten.
