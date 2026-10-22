# Is the VD mass estimator original?

**Deep Research assessment — 5 September 2026**  
**Status:** local, unsubmitted research handoff for VD-docs; findings and proposed corrections, not ratified estimator policy.  
**Comparison target:** the positive compatible-pair estimator at managed VD-docs commit `58dba2d5c3d5d77802f3d57e0d1e2ba2f42b7369`, together with the current local zero–zero debate.

## Finding

**Most of the striking individual properties have clear antecedents.** Positive estimates after a zero reading, the noise-scale ratio at zero–zero, a half-Cauchy latent ratio, and unbounded ratio inference belong to established statistical constructions. Their appearance here is a coherent rediscovery and application of that mathematics.

**I did not find the complete VD estimator presented with its exact combination of a shared latent direction, positive magnitudes, ratio-of-means readout, two composition requirements, and operational admission rules.** That combination remains a candidate contribution. An unsuccessful exact-match search does not establish its originality or uniqueness.

The research also changes the argument we were making: **divergent latent-ratio moments are not confined to zero–zero.** Under the default Gaussian/flat-magnitude policy, they occur at every finite empirical input. Therefore, “the uncertainty becomes infinite, so the estimator automatically recognizes the null trial” is too strong. The unrestricted compatibility set at zero–zero is the stronger identification result.

## 1. What was compared

The estimator reconciles measured vectors with compatible latent pairs

$$
\mathbf F^*=f\mathbf u,\qquad \mathbf a^*=\alpha\mathbf u,
\qquad f>0,\ \alpha>0,\ \|\mathbf u\|=1.
$$

For the default independent isotropic Gaussian model, its normalized support is proportional to

$$
\exp\!\left[-\frac{\|\widetilde{\mathbf F}-f\mathbf u\|^2}{2\sigma_F^2}
-\frac{\|\widetilde{\mathbf a}-\alpha\mathbf u\|^2}{2\sigma_a^2}\right]
\,df\,d\alpha\,d\Omega.
$$

The declared point is

$$
\widehat m=\frac{E[f]}{E[\alpha]},\qquad
\widehat\mu=\frac{E[\alpha]}{E[f]}.
$$

This is a ratio of positive latent-magnitude expectations. It is neither the mean of the latent ratio $M=f/\alpha$ nor an ordinary least-squares slope. The flat magnitude measure is a substantive modeling choice. Calling it a support policy rather than a prior does not remove its influence on normalized probabilities.

Three different outputs must be distinguished throughout:

| Output | Question answered |
|---|---|
| Induced conditional distribution of $M=f/\alpha$ | How does the chosen normalized latent-state weight distribute mass values? |
| Fixed-mass compatibility set | Which masses can fit these observations within the chosen discrepancy threshold? |
| Repeated-sampling distribution of $\widehat m$ | How would the complete point-estimation procedure vary across repeated experiments? |

## 2. Property-by-property assessment

| Property | Closest established antecedent | Assessment |
|---|---|---|
| Positive magnitude estimate from a zero Gaussian reading | Katz's positive normal-mean estimator; Roe–Woodroofe's positive posterior | Exact scalar building block; high confidence |
| Small positive $0/A$, large positive $F/0$ | The same positive normal-mean formula, applied to numerator and denominator | Immediate limiting consequences; no exact Newton application located |
| $\widehat m=\sigma_F/\sigma_a$ at zero–zero | Half-normal expectations and exchange symmetry | Standard consequence; not a defensible standalone novelty claim |
| Half-Cauchy induced mass distribution | Normal-ratio theory; Liseo's null ratio posterior; Kim's folded-normal ratios | Very close or exact distributional antecedents |
| Divergent first/second ratio moments | Classical normal-ratio tail pathology | Established; not exclusive to null data |
| Every positive mass compatible at zero–zero | Fieller geometry; flat ratio profile likelihood | Established identification geometry |
| Ratio of posterior means | Explicitly used in other statistical applications, including Little | General device already known; purposes differ |
| Reciprocal mass/inverse-mass points | Algebraic symmetry of ratios; equivariant medians also qualify | Does not uniquely select this readout |
| Acceleration-weighted mean interpretation | Coordinate size biasing | Exact specialization of known change-of-measure mathematics |
| Both Newtonian composition requirements selecting the readout | General expectation linearity is established | Exact characterization not located; current proof remains a sketch |
| Finite point accompanied by uncertainty and provenance | Standard metrological reporting | Established architecture; VD admission rules may add a particular implementation |
| Complete VD construction and constitutive interpretation | No complete match found | Unresolved novelty; not established by the component results |

The sources and the differences behind these judgments follow. “Not located” means absent from this bounded search, not absent from the literature.

## 3. The zero-channel mechanism dates back at least to Katz

For a scalar Gaussian reading $x$ of a nonnegative quantity, a flat positive prior produces the truncated-normal posterior. Its mean is

$$
h_\sigma(x)=x+\sigma\frac{\phi(x/\sigma)}{\Phi(x/\sigma)}.
$$

Chang, Shinozaki and Strawderman's 2017 research manuscript explicitly reproduces this as Katz's 1961 estimator, with the uniform-positive-prior interpretation, on page 2, equation (1.4). Katz's original article is *Admissible and Minimax Estimates of Parameters in Truncated Spaces*. The original full text was not accessible in this search; the displayed formula was independently checked in the later manuscript. [Chang et al., equation (1.4)](https://m.sc.niigata-u.ac.jp/~hirukawa/seminar/niigata2017_program/Chang.pdf), [Katz, 1961](https://doi.org/10.1214/aoms/1177705146).

Roe and Woodroofe also give the positive Gaussian posterior explicitly in a physical-measurement setting: their equations (1)–(2) normalize the likelihood over the nonnegative parameter domain. This is an especially transparent antecedent for treating a measured zero as uncertain evidence about a positive latent quantity. [Roe and Woodroofe, *Setting Confidence Belts*, 2000 preprint](https://arxiv.org/pdf/hep-ex/0007048).

**Our specialization to the VD readout:** put $c=\sqrt{2/\pi}$. Then $h_\sigma(0)=c\sigma$, while $h_\sigma(x)\approx x$ for a well-resolved positive reading. Consequently:

| Force and acceleration readings | Ratio-of-means result |
|---|---|
| $0,A$, with $A\gg\sigma_a$ | $\widehat m\approx c\sigma_F/A$ |
| $F,0$, with $F\gg\sigma_F$ | $\widehat m\approx F/(c\sigma_a)$ |
| $0,0$ | $\widehat m=\sigma_F/\sigma_a$ |

The first two rows are scalar/known-direction limits. The full vector estimator also integrates over direction, so those approximations are not universal closed forms. In particular, the two one-zero cases do **not** both give “sigma divided by five”: exchanging the channels reverses the ratio and its units.

These formulas are compelling evidence of mathematical ancestry. They do not show that Katz proposed the full Newton estimator.

## 4. The sigma/sigma and half-Cauchy results have very close precedents

At zero–zero, the default model factorizes in the positive magnitudes:

$$
f=\sigma_F|Z_1|,\qquad \alpha=\sigma_a|Z_2|,
\qquad Z_1,Z_2\stackrel{\mathrm{iid}}\sim N(0,1).
$$

Writing $s=\sigma_F/\sigma_a$, direct transformation gives

$$
p_M(m)=\int_0^\infty \alpha\,p_f(m\alpha)p_\alpha(\alpha)\,d\alpha
=\frac{2s}{\pi(s^2+m^2)},\qquad m>0.
$$

Thus the induced ratio is half-Cauchy, its median is $s$, and the ratio-of-means point also happens to be $s$. A scale of $1\,\mathrm{kg}$ occurs only when the **dimensional uncertainty ratio** equals $1\,\mathrm{kg}$.

Marsaglia's normal-ratio work explicitly derives the Cauchy case and standardization by the ratio of standard deviations. Kim's 2006 paper on ratios of folded normal variables explicitly includes half-standard Cauchy in its family. The folded and positive-truncated normal constructions coincide at zero location; they generally differ away from zero. Kim's abstract was inspected, not the full theorem. [Marsaglia, 2006](https://www.jstatsoft.org/article/view/v016i04), [Kim, 2006](https://doi.org/10.1080/03610920600672229).

**The closest Bayesian comparison is Liseo (2003).** His equations (10), (12) and (13), pages 137–138, use a prior flat in the original normal means, expressed with its Jacobian in ratio coordinates. Setting both observed means to zero reduces the marginal ratio posterior to Cauchy. Restricting to the positive quadrant and rescaling gives our half-Cauchy specialization. His profile likelihood, equation (5), is instead flat at null. Both sides of our probability-versus-compatibility distinction are already present in that analysis. [Liseo, *Bayesian and conditional frequentist analyses of the Fieller's problem*, 2003](https://www.researchgate.net/publication/5182203_Bayesian_and_conditional_frequentist_analyses_of_the_Fieller%27s_problem_A_critical_review).

An additional elementary observation limits uniqueness claims: if $f/\sigma_F$ and $\alpha/\sigma_a$ are exchangeable with finite positive means, their means agree. The point is then $s$ even without the particular Gaussian derivation. The null point alone does not identify a unique measure or readout.

## 5. Infinite moments do not finish the argument

Marsaglia emphasizes that normal ratios can have nonexistent moments even when their central distributions look well behaved at high signal. This is not just a pathological zero–zero phenomenon. [Marsaglia, section 4](https://pdfs.semanticscholar.org/9ead/12ca2b04f00533587ecfe90df7ea334d2f58.pdf).

**We can establish the relevant result directly for the actual VD model.** Integrate out direction and call the normalized positive-magnitude density $b(f,\alpha)$. With finite observations, nondegenerate Gaussian covariance, and the default flat magnitude measure, this density remains continuous and strictly positive as $\alpha$ approaches zero while $f$ stays in any compact interval $[f_0,f_1]$ with $f_0>0$.

On a sufficiently small rectangle, it therefore has a positive lower bound. Its contribution to the ratio mean is bounded below by a positive constant times

$$
\int_{f_0}^{f_1} f\,df\int_0^\varepsilon\frac{d\alpha}{\alpha}=\infty.
$$

The second raw moment also diverges. This proof does not require the two magnitudes to remain independent after direction integration. It applies to **every finite empirical observation under that default policy**, including strong codirectional data. At high signal, the boundary weight can be exceptionally tiny; that changes practical tail probabilities, not mathematical divergence.

Strict wording matters: the induced ratio has infinite positive mean and second raw moment, so it has **no conventional finite standard deviation about a finite mean**. Its mean-square deviation about the reported finite point is infinite. Writing “$1\,\mathrm{kg}\pm\infty$” conceals these distinctions.

Nor does this establish infinite sampling variance of $\widehat m$. In the independent fixed-direction model, $\widehat m=h_{\sigma_F}(X)/h_{\sigma_a}(Y)$. For large negative $y$, $h_\sigma(y)\sim\sigma^2/|y|$; for large positive $y$, it grows like $y$. Its reciprocal therefore grows at most linearly in the negative tail. Gaussian readings give this ratio-of-means rule finite sampling moments. This is a deduction for that simplified model, not a completed sampling analysis of the full vector algorithm.

**The measure also matters.** As a counterexample, replace the flat magnitude weight at null by a weight proportional to $(f/\sigma_F)^2(\alpha/\sigma_a)^2$. The standardized magnitudes become independent chi variables with three degrees of freedom. Exchange symmetry still gives $\widehat m=s$, but

$$
E[M]=\frac4\pi s,\qquad E[M^2]=3s^2,
\qquad \operatorname{Var}(M)=\left(3-\frac{16}{\pi^2}\right)s^2<\infty.
$$

This is a sensitivity counterexample, not a recommended replacement policy. A sigma/sigma point does not force infinite ratio uncertainty.

The revised conclusion is that a positive point can coexist with an honest record of non-identification. **Infinite moments alone are insufficient to make that record honest or diagnostic.**

## 6. Fieller geometry supplies the clearer distinction between the zero cases

Fieller-type ratio confidence sets can be bounded, unbounded, or the whole real line. Von Luxburg and Franz give a geometric construction: when the uncertainty ellipse contains the origin, their ratio set is unrestricted. Restricting the parameter domain to positive mass gives the corresponding all-positive set. This is a direct antecedent for allowing a trial to constrain neither small nor large mass. [Von Luxburg and Franz, construction 1 and theorem 1](https://arxiv.org/pdf/0711.0198).

**Our scalar specialization:** with independent Gaussian errors, profiling over latent acceleration gives

$$
Q(m)=\frac{(F-mA)^2}{\sigma_F^2+m^2\sigma_a^2}.
$$

For the nonnegative readings considered here, the minimizing latent state lies in the allowed region or its boundary limit. At zero–zero, $Q(m)=0$ for every $m>0$. No mass fits better than another by this criterion.

For an absolute threshold $k>0$, the set $Q(m)\le k$ has useful one-zero forms:

| Readings | Positive masses admitted, when the nonzero channel is resolved |
|---|---|
| $0,A$, $A^2>k\sigma_a^2$ | $0<m\le\sqrt{k}\sigma_F/\sqrt{A^2-k\sigma_a^2}$ |
| $F,0$, $F^2>k\sigma_F^2$ | $m\ge\sqrt{F^2-k\sigma_F^2}/(\sqrt{k}\sigma_a)$ |
| $0,0$ | Every $m>0$ |

If the nonzero channel does not clear the displayed condition, even that one-zero case admits all positive masses at this threshold. A nonzero numeral alone does not guarantee identification. These are discrepancy sets; their frequentist coverage requires appropriate calibration, especially for boundaries and the full vector model.

This also resolves what the Sun objection tests. The null half-Cauchy has central 95% probability interval approximately $[0.0393s,25.45s]$, despite divergent moments. With $s=1\,\mathrm{kg}$, it assigns little probability to enormously larger masses. The flat compatibility profile does not exclude those masses. A policy-generated probability interval must not be presented as if the null trial had identified the object's mass scale.

## 7. The mass-estimation literature already discusses zero force and zero acceleration

Hans Butler's 2004 published technical disclosure *Controlling a position of a mass, especially in a lithographic apparatus* describes estimating stage mass from force and acceleration. Section 2.4 explicitly discusses a constant-velocity phase with nominal controller force and acceleration both zero: noise dominates, and adaptation is switched off or bounded. Section 2.3 discusses nonzero force offsets at zero acceleration being interpreted as infinite mass. This is a striking direct application precedent, although its recursive least-squares estimator differs from ours. [EP1455231A2, sections 2.3–2.4, published 8 September 2004](https://patents.google.com/patent/EP1455231A2/en).

It therefore would be inaccurate to claim that prior mass-estimation work simply ignored the zeros. Our potential difference is the particular positive readout and the way VD records and admits its result.

Uncertainty in both measurement channels is likewise established methodology. Bartel, Stoudt and Possolo present errors-in-variables regression and Monte Carlo uncertainty evaluation for force calibration. That paper supports the general methodological ancestry, not the exact VD readout. [NIST/Metrologia, 2016](https://www.nist.gov/publications/force-calibrations-using-errors-variables-regression-and-monte-carlo-uncertainty).

## 8. What remains distinctive about the readout argument?

Ratios of posterior expectations already occur in other fields. Little explicitly replaces a posterior expectation of a ratio by a ratio of posterior expectations in survey inference, as an approximation ignoring terms of order $1/n$. Our readout instead selects that functional exactly, including low-signal cases where the ratio mean diverges. The precedent establishes that the device is known; it does not justify our choice. [Little, 2012, page 319](https://www.scb.se/contentassets/ca21efb41fee47d293bbee5bf7be7fb3/calibrated-bayes-an-alternative-inferential-paradigm-for-official-statistics.pdf).

The weighted-mean identity also has a standard interpretation. Under coordinate size biasing,

$$
dQ_\alpha=\frac{\alpha}{E_P[\alpha]}dP,
\qquad \widehat m=E_{Q_\alpha}[M].
$$

The inverse point uses a different tilt, $dQ_f=f\,dP/E_P[f]$, to give $\widehat\mu=E_{Q_f}[1/M]$. These are applications of the coordinate-biasing definition in Arratia, Goldstein and Kochman, section 2.3, equations (19)–(20). They do not turn the point into $E_P[M]$. [*Size bias for one and all*, first posted 2013](https://arxiv.org/pdf/1308.2729).

Reciprocity alone is weak selection evidence. Any separable rule $T(f)/T(\alpha)$ reciprocates when its inputs are exchanged. For a continuous positive ratio, quantiles obey $q_p(1/M)=1/q_{1-p}(M)$; the median is self-reciprocal, but a general same-level quantile is not. The managed audit's blanket quantile wording needs that qualification.

The more interesting proposed result is the pair of composition identities:

$$
R(f_1+f_2,\alpha)=R(f_1,\alpha)+R(f_2,\alpha),
$$

$$
\frac1{R(f,\alpha_1+\alpha_2)}
=\frac1{R(f,\alpha_1)}+\frac1{R(f,\alpha_2)}.
$$

The ratio of expectations satisfies both on the same propagated joint law. I found no exact published characterization using this particular pair as Newtonian requirements. However, the managed document explicitly offers a sketch and narrows to a regular family during the argument. It has not yet established the full advertised uniqueness theorem.

Composition alone also holds for $R_L(f,\alpha)=L(f)/L(\alpha)$ with a suitable common positive linear functional $L$. Selecting specifically expectation under the declared law requires the law-dependence and regularity assumptions to do real work. Moreover, composing already propagated latent laws is different from re-estimating a combined body using a new likelihood and measure.

**This characterization is the strongest candidate mathematical contribution exposed by the search.** Its originality and correctness should be investigated together through a complete theorem or counterexample analysis.

## 9. Implications for “the estimator is the mass”

Metrology already treats a measurement result as an attributed quantity value accompanied by relevant information, generally including uncertainty. Thus “the usable empirical result includes more than a bare scalar” has an established home. It does not follow that any particular estimator is uniquely constitutive of the physical quantity. [JCGM VIM, entry 2.9](https://jcgm.bipm.org/vim/en/2.9.html).

The interacting-object argument also remains useful: cancellation of measured force contributions need not establish an exact physical zero. Propagating a sum requires the contribution covariances,

$$
\operatorname{Cov}\!\left(\sum_i\epsilon_i\right)
=\sum_{i,j}\operatorname{Cov}(\epsilon_i,\epsilon_j).
$$

This is standard uncertainty propagation, including correlated inputs. Nonzero uncertainty for every contribution does not by itself guarantee a nonzero uncertainty for the sum if correlations impose exact cancellation. [JCGM GUM, section 5.2](https://www.iso.org/sites/JCGM/GUM/JCGM100/C045315e-html/C045315e_FILES/MAIN_C045315e/05_e.html).

Removing established freeness from this attempted inference is a proposed VD routing decision. It must preserve the managed correction that an already supplied mass can still be used with Newton II to compute a free massive particle's acceleration. Neither the statistical literature nor our calculations establish that all physical mass must be defined through this single trial procedure.

My position after the research is more precise than earlier enthusiasm: **the construction is a defensible candidate empirical mass operation, but the infinity observation does not settle its identification policy or prove that it uniquely “is” mass.** The choice between an admitted scalar with standing and an explicit not-identified result remains a specification decision.

## 10. What can responsibly be claimed now

A defensible description is:

> We develop a positive compatible-pair estimator for inertial mass, using established positive-parameter and ratio-inference mathematics. Its proposed contribution is the particular Newtonian composition justification and the explicit separation of point readout, compatibility, uncertainty, provenance, and admission for downstream use.

Before making a stronger originality claim:

1. Complete the two-composition characterization, including its domain, continuity, law dependence, and coupling assumptions.
2. Specify which of the three uncertainty objects is reported and which determines admission. Do not use divergent moments as the null detector.
3. State the measure dependence and examine repeated trials, where a common-mass latent measure requires additional choices.
4. Compare the full procedure against constrained errors-in-variables and Bayesian ratio methods in informative, weak, null, and incompatible experiments.

These are focused ways to establish a contribution. The old components do not make the synthesis worthless; they determine what the synthesis must add.

## 11. Sources, search coverage, and handoff boundaries

The search covered statistical ratio inference, positive normal-mean estimation, reciprocal and composition characterizations, physical measurement, inertial-mass identification, and official metrology guidance. Discovery used three research lanes, followed by primary-source reading and independent checks of the decisive equations. Representative query families included “mass estimator ratio of means,” “force acceleration posterior means,” “mass estimation half-Cauchy,” positive normal means/Katz, Fieller zero denominator, ratio of posterior expectations, coordinate size bias, and mass estimation/persistent excitation. Exact-phrase follow-ups did not locate the full VD construction. Searches were conducted on 5 September 2026 without a lower date restriction.

**Evidence limitations:** this is a bounded literature assessment, not a proof of first publication. It is primarily English-language and web-accessible; it does not exhaust subscription databases, dissertations, historical operationalist literature, or unpublished work. Katz's original full article and Kim's full theorem were not accessible. Katz's formula was verified in an original follow-on research manuscript; Kim's claim is limited to its publisher abstract. Liseo's full article text was read through ResearchGate because the publisher-hosted PDF failed to open. Marsaglia's article was read in a PDF mirror with its journal identity intact. Patent text was used as a published technical disclosure, not to assess patent rights.

**Managed/current sources:** `catalog.yaml` supplied inventory and status at `canvalidk/VD-docs@58dba2d5c3d5d77802f3d57e0d1e2ba2f42b7369`. The following canonical sources were read at that same commit for this assessment:

- `Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`
- `Newton-analysis/_dscn_mass_estimator/what_forces_the_ratio_of_means_readout.md`
- `Newton-analysis/_dscn_mass_estimator/what_forces_the_chain_to_the_estimator.md`
- `Newton-analysis/_note_why_mass_estimation/why_inertial_mass_must_be_estimated.md`
- `Newton-analysis/vd_design_3_newton_i_correction.md`

The current full specification defines the comparison target. Older stationary-root and median proposals are context, not interchangeable versions of this estimator. Their dedicated historical files were not exhaustively reread for this research.

**Local active/unsubmitted sources:** `pass to VD-docs/mass_estimation_zero_zero_two_cases.md` supplies the conversation's current competing arguments and infinity proposal; `08_nm_lesson_drafts/nml3/NML3_FORCE_SUM_TRIPLET_NAMING_DECISION_2026-08-26.md` supplies local ratified terminology. A scoped filename search also located older entry packets and handoffs; these were not substituted for the managed estimator or audited exhaustively. No separate inbox/unclassified collection was exhaustively searched. Accordingly, the report assesses the specified current construction rather than claiming complete coverage of every VD note.

**Handoff action:** preserve this report as a new research record alongside the two-cases argument. Its principal correction strengthens that argument's existing warning about divergent moments: under the default policy, the issue extends to all finite empirical inputs. Managed sources and the competing advocacy draft were not edited by this research. Saving this file in `pass to VD-docs/` does not itself submit it to the managed collection.
