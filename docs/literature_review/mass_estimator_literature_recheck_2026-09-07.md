# Mass estimator literature recheck — 7 September 2026

Status: local research record awaiting transfer to VD-docs. No managed source changed or submission made.

## Finding

The claim that nobody has discussed this mathematics is not supported. Several ingredients have exact antecedents. This search still did not locate a publication specifying the complete three-dimensional compatible-pair integral, its flat positive-magnitude measure, and its particular ratio-of-magnitude-means point together. That is an unresolved priority question, not evidence of first discovery.

This pass strengthens the September 5 assessment with an explicit scalar equivalence check and examines additional work on symmetric regression and Bayesian shape alignment. It concerns the equation, not an empirical rejection gate: the current local steering conditions the investigation on Newton II being correct.

## Exact comparison target

With measured vectors y=(F,a), latent x=(f u,alpha u), positive magnitudes and unit common direction, the default construction is

$$
Q=(y-x)^T\Sigma^{-1}(y-x),\qquad
\widehat m=
\frac{\int f\exp(-Q/2)\,df\,d\alpha\,d\Omega}
     {\int \alpha\exp(-Q/2)\,df\,d\alpha\,d\Omega}.
$$

The magnitude measure is flat and direction uniform. The point is a ratio of expectations under the same normalized compatible-pair law. It is not the posterior mean of f/alpha or the minimum-discrepancy slope. These differences determine whether a candidate is actually the same estimator.

## 1. Known-direction scalar point: exact older building blocks

The positive normal-mean formula attributed to Katz (1961) is explicitly reproduced in equation (1.4) of [Chang, Shinozaki and Strawderman's 2017 manuscript](https://m.sc.niigata-u.ac.jp/~hirukawa/seminar/niigata2017_program/Chang.pdf):

$$
h_\sigma(x)=x+\sigma\frac{\phi(x/\sigma)}{\Phi(x/\sigma)}.
$$

It is the mean obtained from a Gaussian observation and a flat prior over the nonnegative parameter. The original Katz article was not accessible in this pass; attribution rests on that original follow-on research manuscript. The corresponding published 2017 article lists Chang and Strawderman, while the linked manuscript also lists Shinozaki; do not interchange those bibliographic records.

**Our mathematical specialization:** with known common direction and independent Gaussian channels, the magnitude likelihood factorizes. Therefore our point is exactly

$$
\widehat m=\frac{h_{\sigma_F}(F)}{h_{\sigma_a}(a)}.
$$

At two zero readings, both means are their respective sigma times sqrt(2/pi), giving sigma_F/sigma_a. This establishes the ancestry of the zero-channel mechanism; it does not establish that Katz proposed this mass readout.

## 2. Unknown-sign scalar mass law: exact specialization of ratio inference

[Liseo (2003), section 3, equations (10), (12), (13)](https://www.researchgate.net/publication/5182203_Bayesian_and_conditional_frequentist_analyses_of_the_Fieller%27s_problem_A_critical_review) gives a normal-means ratio posterior using a prior flat in the original two means. In ratio/nuisance coordinates the prior contains the Jacobian |lambda|. The author-uploaded article was inspected; the publisher PDF timed out.

**Our equivalence derivation:** let omega_1 and omega_2 be the two underlying signed scalar means. Restrict that flat-means model to omega_1*omega_2>0. Write

$$
\omega_1=u f,\quad\omega_2=u\alpha,
\quad u\in\{-1,+1\},\quad f,\alpha>0.
$$

Each quadrant has unit Jacobian and the two signs have equal base weight. This is exactly our one-dimensional compatible-pair law. Unequal known channel uncertainties are handled by standardizing each channel and restoring the mass scale sigma_F/sigma_a. The equivalence applies away from zero as well as at zero. Positivity conditioning here is our specialization; it is not claimed as Liseo's stated physical application.

This identifies an exact scalar distributional precedent. It does not identify our E[f]/E[alpha] point as Liseo's proposed estimator, or reproduce the unknown-direction sphere integral in three spatial dimensions.

## 3. Ratio-of-posterior-expectations is an established operation

[Little (2012), page 319](https://www.scb.se/contentassets/ca21efb41fee47d293bbee5bf7be7fb3/calibrated-bayes-an-alternative-inferential-paradigm-for-official-statistics.pdf) explicitly uses a ratio of posterior expectations in survey inference. He uses it as an approximation to a posterior expectation of a ratio, neglecting terms of order 1/n. Our equation selects that functional directly, including cases where the ratio mean diverges. This is a precedent for the operation, with a different justification and model.

## 4. Symmetric inference with error in both channels is established

[Leonard (2011), Estimating a bivariate linear relationship](https://arxiv.org/pdf/1202.0957), equations (1)–(6) and (13)–(21), develops errors-in-variables inference invariant under exchanging and rescaling the two observed coordinates. It supplies a slope distribution and uses a Cauchy prior for a standardized slope. Its model has multiple observations, unspecified error variances and normally distributed true values. Those choices differ from our specified-covariance compatible-pair measure and point readout. It is a close methodological precedent, not an exact equation match.

An additional vector-related candidate, [Mardia et al. (2013), Bayesian alignment of similarity shapes](https://arxiv.org/pdf/1312.1840), models noisy configurations and positive scale with latent alignment. Equations (3)–(5) give a conditional scale distribution with a gamma prior. This differs from our flat magnitude-pair model and ratio of means; no full match was located there.

## 5. Force/acceleration mass estimation and zeros were discussed directly

[Hans Butler's 2004 technical disclosure, EP1455231A2](https://patents.google.com/patent/EP1455231A2/en), sections 2.3–2.4, explicitly discusses nonzero force offsets at zero acceleration and noise-dominated phases with nominal force and acceleration both zero. Its response includes filtering and limiting or stopping recursive least-squares adaptation. This is a direct physical precedent for the problem, using a different estimation rule. The disclosure is cited as technical literature, without any assessment of patent rights.

## Independent equivalence checks

New script: `03_trace_and_evaluator/mass_estimator/literature_equivalence_check.py`. Results: `results/literature_equivalence.json` in the same directory. The existing estimator was not changed.

Fourteen scenarios cover positive, negative, opposite-sign, one-zero and both-zero readings with two uncertainty-scale ratios. The script compares:

- The whole scalar angle density against Liseo's equation (12), conditioned on positive ratio. Maximum relative difference: **7.99e-15**.
- The known-positive-direction point against the quotient of the two Katz formulas. Maximum relative difference: **6.89e-14**.

These calculations check the formula translation and numerical implementation. The change-of-variables argument establishes the model equivalence; a finite test grid alone would not prove it.

## What can be said

We can say: “We formulate a particular positive, symmetric mass estimator using established constrained normal-mean and ratio-inference mathematics. Its scalar uncertainty model has a direct classical specialization. We have not located the complete vector construction with this exact readout in the literature searched.”

We cannot yet say that the extended equation is new, that the zero/zero result is a new statistical discovery, or that no previous researcher discussed mass inference at zeros. Conversely, finding its components does not establish that a source contains the exact full estimator. The composition characterization is a separate possible contribution; its historical priority was not established here.

## Search scope and custody

Managed snapshot: `canvalidk/VD-docs@adfe14bbaafb29856f0a30da088495b7c9d483d2`, resolved from main for this pass. Consulted the relevant catalog entries and these canonical files at that commit:

- `Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`
- `Newton-analysis/_dscn_mass_estimator/mass_estimator_originality_research_2026-09-05.md`

Local active/unsubmitted scope: current estimator implementation; `../Mass estimator/mass_estimator_law_assumed_correct_steering_2026-09-07.md`; the opening statement of `../Mass estimator/mass_estimator_practical_limits_and_uncertainty_2026-09-07.md`; the composition theorem's stated assumptions and conclusion in `../Mass estimator/mass_readout_composition_uniqueness_2026-09-06.md`. All named local records are under `pass to VD-docs/`. The local theorem is later than the proof-sketch assessment in the September 5 research; that old statement is not treated as the present proof status. No exhaustive historical, inbox, or local-document audit was performed.

Fresh query families included exact phrases “ratio of posterior means,” “ratio of posterior expectations,” “force acceleration ratio of means,” “force acceleration half-Cauchy,” “positive collinear posterior magnitudes,” Bayesian scale-factor/vector estimation, Bayesian total least squares, reciprocal ratio-of-means estimation, positive normal-mean/Katz estimation, and Liseo/Fieller inference. Reference following led to the primary materials above. Roe–Woodroofe's positive-parameter paper and other applications were also screened, without a closer full-equation match.

This is a bounded English-language web search and equation comparison, not an exhaustive subscription-database search, citation-network audit, or proof of first publication. An exact construction may be published under other notation or terminology. No contact with authors was made.
