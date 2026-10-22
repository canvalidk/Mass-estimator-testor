# Mass-equation comparison: institutional access and five new source groups

Date: 2026-09-08. Status: local active research record awaiting transfer to
VD-docs; NOT submitted. The user authorized reading Katz's original and
examining five fresh candidates, then stopping to discuss findings and
diminishing returns. This round is complete. Mardia was not reopened.

## Main finding

The broader search was useful. Huard–Mailhot (2006) and Mantz (2016) provide
closer Bayesian frameworks than the earlier least-squares and scalar-boundary
comparisons. Their published likelihood and latent-input machinery can support
the entire stipulated observation model. We explicitly instantiated that
machinery with fixed, declared priors and proved that the resulting point
estimate and quantiles cover all finite observations with positive-definite
joint Gaussian covariance. This is a successful functional competitor
construction within earlier published frameworks, not another failure case.

The extra priors and point-summary choices are ours in this audit. We have not
found an inspected publication that prints the complete target measure and
point prescription. It would be inaccurate either to attribute our whole
construction to those authors or to conclude that the broad functional
capability is unavailable in previous Bayesian methods.

Consequently, finite positive output, directional signs, full covariance and
zero coverage alone do not strongly distinguish the target mathematically.
The remaining specific comparison concerns its amplitude–direction measure,
its ratio-of-magnitude-means decision, the resulting mass law, and any
demonstrated advantages. Historical priority for that combination remains
unresolved; no exact predecessor was established in the inspected sources.

## Completed source groups and verdicts

| Source | Verified published contribution | Verdict for this comparison |
|---|---|---|
| [Katz 1961, original](https://www.jstor.org/stable/2237613) | Restricted scalar-parameter point estimation; p.139 directly gives the positive normal-mean rule | Earlier access gap closed; scalar coverage confirmed; original scope does not supply the two-vector mass construction |
| [Marchand–Strawderman 2004](https://doi.org/10.1214/lnms/1196285377) | Multivariate restricted-normal mean estimation and a two-vector difference-constraint result | Demonstrated limitations of direct substitution: required convex ambient region and difference restriction do not describe positive proportionality; exact predecessor unresolved |
| [Huard–Mailhot 2006](https://espace.inrs.ca/9481/1/P1330.pdf) | General vector input/output error likelihood, latent-input integration, and a positive-slope calibration example | Full coverage proved for our explicit instantiation of the published framework; their worked example alone does not establish the full contract |
| [Kelly 2007](https://arxiv.org/pdf/0705.2774) | Bayesian regression with noisy correlated channels and a latent Gaussian-mixture input model | Native defaults do not meet positive-scale contract; a particular flat-positive-slope specialization is improper at zero; proper-prior adaptation succeeds for scalar pairs |
| [Mantz 2016; preprint 2015](https://arxiv.org/pdf/1509.00908) | Multivariate extension with arbitrary covariance among all measurements of one object | Full required observation likelihood is present; our common-positive-scale/proper-prior specialization has full coverage; unchanged default sampler and target-equation identity are not established |
| [Zhou–Kou–Li–Fang 2016](https://uwaterloo.ca/geospatial-intelligence/sites/default/files/uploads/files/2016-zhou-kou-li-fang-jse.pdf) | Structured errors-in-variables fitting with full cross-covariance and local uncertainty approximation | Demonstrated limitation: positive scale optimizer need not exist even with strictly positive-definite correlated errors |

Five groups are new. Katz is an update to an existing group. Across all three
rounds there are fourteen audited source groups, including background and
framework sources; this does not mean fourteen independent complete competitors.
It is also different from the original supplied fourteen-item list.

Detailed evidence and primary equation anchors:

- [Updated Katz note](mass_candidate_katz_evidence_2026-09-08.md).
- [Marchand–Strawderman note](mass_candidate_marchand_strawderman_evidence_2026-09-08.md).
- [Huard–Mailhot note](mass_candidate_huard_mailhot_evidence_2026-09-08.md).
- [Kelly note](mass_candidate_kelly_evidence_2026-09-08.md).
- [Mantz note](mass_candidate_mantz_evidence_2026-09-08.md).
- [Zhou et al. note](mass_candidate_zhou_etal_evidence_2026-09-08.md).

## A concrete full-coverage alternative

The following choices and proof were made in this audit, using the general
latent-input Bayesian construction in Huard–Mailhot Eq.(3) and the full
measurement covariance of Mantz Eq.(1). They are not quoted as those papers'
printed positive-mass estimator.

Let `v=a*` be the true acceleration vector and `m>0`. In fixed reference units,
choose once `v~N_d(0,I)` and `log(m)~N(0,1)`, independently. These are declared
model assumptions, held fixed across every input and covariance, rather than
additional measurements or per-case adjustments. Positive scales different
from one can instead be specified explicitly. Let

$$
z=(F,a),\qquad B_m=\begin{bmatrix}mI\\I\end{bmatrix},\qquad
z\mid(m,v)\sim N_{2d}(B_mv,\Sigma),\quad\Sigma\succ0.
$$

The true force is exactly `mv`: no fitted rotation, offset, or intrinsic
scatter is needed. Integration over the true vector gives

$$
p(m\mid z,\Sigma)\propto
\pi(m)\,\phi_{2d}(z;0,\Sigma+B_mB_m^T),\qquad m>0.
$$

Use its mean as the estimate and its quantiles as uncertainty. The priors have
support at every finite true vector and every positive mass; no upper mass
cutoff is imposed. This construction generally differs from the target's
prior, point, tails and exchange properties. Its legitimacy under the broad
contract does not establish that it is preferable for the intended use.

For any fixed finite observation and SPD covariance, the unintegrated
Gaussian likelihood is strictly positive and at most
`C_Sigma=(2*pi)^(-d)|Sigma|^(-1/2)`. Integrating against the proper vector
prior preserves this bound. The posterior normalizer therefore satisfies
`0<Z<=C_Sigma`. Its mean satisfies

$$
0<E[m\mid z,\Sigma]
\le C_\Sigma E_\pi[m]/Z<\infty.
$$

The mass density is positive throughout `(0,infinity)`, so every interior
quantile is finite and strictly positive. Local dominated convergence over
finite observations and the open SPD domain proves continuity of the
normalizer, mean and CDF; strict increase of the CDF gives continuous
interior quantiles. These arguments establish coverage everywhere, including
measured zeros; they do not assert statistical identification at zero.

Directional signs survive integration. With unit error covariance the
exponent contains a nonzero multiple of `m F dot a`. Reversing one vector
changes the mass likelihood. Arbitrary covariance uses the whole quadratic
form. Thus this construction does not meet positivity merely by ignoring
direction or always returning a fixed prior mean.

At both measured vectors zero, its answer depends on the declared priors
and the covariance. For example, with 3D unit latent prior and measurement
standard deviation 0.25 in each channel, its mean is approximately 0.6193,
median 0.4919, and central 95% interval `[0.0996,1.8894]`. The target's
equal-error zero–zero point is 1 with a half-Cauchy uncertainty law.
Different weak-data answers expose the different assumptions; neither is
an assumption-free identification of a physical mass.

## Where the target sits in that same likelihood family

This is our change-of-variables calculation, independently checked in two
lanes. It does not attribute a previously unpublished choice to a source.

The target uses compatible vectors `(fu,alpha u)`, positive amplitudes and
common unit direction, with base measure `df d(alpha) dOmega`. Set
`m=f/alpha` and `v=alpha u`. Since `df d(alpha)=alpha dm d(alpha)` and
`d^dv=alpha^(d-1)d(alpha)dOmega`,

$$
df\,d\alpha\,d\Omega
=\|v\|^{2-d}\,dm\,d^dv.
$$

The target therefore has the same errors-in-variables likelihood, with
improper prior measure `1_{m>0} ||v||^(2-d) dm d^dv`, and point

$$
\widehat m_{\rm target}
=\frac{E[m\|v\|\mid z]}{E[\|v\|\mid z]}.
$$

For one dimension the prior factor is `|v|`; for two it is constant; for
three it is `1/||v||`. The target's induced uncertainty is the posterior
law of `m`. In two dimensions the law is exactly that obtained by
marginalizing the Gaussian EIV likelihood with a flat positive slope and
flat latent-vector measure. This is an exact model identity, not merely
similar outputs on selected examples.

This result identifies what remains to search for: the prior measure and
the point functional. In particular, ordinary posterior mean of mass is
not the target point. The earlier heavy-tail analysis still applies.

Marchand–Strawderman's direct ambient-uniform formula does not automatically
choose this measure. The compatible-pair set is nonconvex and has zero
ambient measure for `d>=2`. Even uniform induced surface measure is
different: the embedding metric gives
`dH_(d+1)=(f^2+alpha^2)^((d-1)/2) df d(alpha) dOmega`.
These distinctions must be explicit before claiming equation equivalence.

## Bounded numerical validation

New script:
`03_trace_and_evaluator/mass_estimator/coverage_third_pass.py`.
Output:
`03_trace_and_evaluator/mass_estimator/results/coverage_third_pass_2026-09-08.json`.

The Gaussian latent vector is integrated analytically. Numerical integration
is one-dimensional in log mass, with stable log densities. No MCMC mixing
claim is needed. All assertions passed using the existing NumPy runtime.

- Eighteen 3D fixtures: both zero, either vector zero, aligned, perpendicular
  and anti-aligned, each under small isotropic, large isotropic and fully
  correlated anisotropic error covariances.
- Twenty-three path points: changing relative angle and passing either
  measured vector continuously through zero.
- Refinement from 6,001 points on log mass `[-10,10]` to 14,401 points on
  `[-12,12]` changed means or quantiles by at most `1.200e-5` relatively.
  This is a numerical refinement estimate, not a certified uniform bound.
- The aligned ratio-two fixture with standard errors 0.03 gave mean
  `1.9979901`. With errors 0.25 the aligned mean was `1.8627259`, while
  the same-length anti-aligned fixture gave `8.7039336`: directional signs
  materially affect the inference.
- Six independent 2D checks compared the target's direction/radial
  integration with analytic Cartesian Gaussian marginalization under the
  reparameterized flat measure. The largest absolute angle-density
  difference was `9.77e-14`, including zero and adverse-angle cases with
  correlated covariance. This verifies the implemented 2D identity; the
  Jacobian argument establishes it mathematically.

For Zhou et al., a separate numeric sweep would add little: the exact
correlated counterexample in its note proves nonexistence of the positive
optimizer. Similarly, Kelly's scalar flat-prior counterexample has a
`1/m` likelihood tail and a divergent normalizer, so the issue is exact
propriety, not a failing numerical sampler. Proper priors can alter it.

## Originality and diminishing returns: review position

The strongest defensible update is:

> The construction belongs to an established Bayesian errors-in-variables
> family. In the publications inspected so far, we have not identified the
> complete combination of its flat amplitude–direction measure, induced
> positive-scale law and ratio-of-posterior-magnitude-means decision. Its
> global finite-input coverage is established, but that capability is also
> attainable within earlier published frameworks under other declared priors.

This is a scoped description of findings, not a declaration of originality.
New combinations or derivations can still be meaningful contributions. An
account of novelty should state the exact combination and any proved or
demonstrated benefit, instead of treating every prior framework as excluded
because its printed example has a different prior, sampler, or application.

The broadening produced material information, so the previous supplied-list
search had not exhausted useful territory. More generic OLS/TLS papers now
have low expected return: the boundary mechanism repeats even after full
cross-covariance is included. Additional generic Bayesian existence papers
would also increasingly repeat the point already established here.

If work resumes after discussion, the more discriminating next step is an
exact equation comparison against the Huard–Mailhot/Mantz model family,
including which priors and point choices appear in older cited sources and
what the target's choice buys. That work has not been started as another
unbounded branch. This round stops here as requested.

## Search coverage, provenance and custody

Managed `canvalidk/VD-docs` main was resolved again via the installed GitHub
connector and is unchanged:
`adfe14bbaafb29856f0a30da088495b7c9d483d2`.
The already-read coherent sources at that commit remain:

- `catalog.yaml`, relevant mass-estimator inventory selection;
- `Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`;
- `Newton-analysis/_dscn_mass_estimator/mass_estimator_originality_research_2026-09-05.md`.

The active local September 7 law-assumed-correct steering remains the task's
authority for applying no rejection gate. The source relationships and
local developments absent from that catalog selection are recorded in the
first and second reviews. No new managed-context reconciliation occurred.

Katz pp.136–142 were read through the user's Imperial-authenticated JSTOR
viewer; the printed normal equation was checked visually on p.139. Browser
access thus closed a real original-text gap. It did not cause all subsequent
discoveries: the five new source groups were accessible through primary
author/institutional copies or publisher text without institutional access.
No claim is made that Scopus or Web of Science was searched in this round.

Selection broadened through the Katz related-literature trail into restricted
parameter estimation; through Bayesian measurement-error regression in
astronomy from Kelly to Mantz; and through positive-slope input-uncertainty
calibration and structured full-covariance EIV. Search concepts included
positive scale/gain, noisy inputs and outputs, Bayesian calibration,
multivariate regression, correlated measurement errors, and structured or
weighted total least squares. Exact titles, DOI/publisher records and
author-hosted texts were used for verification. Readable manuscript vs
publication pagination is explicitly distinguished in the candidate notes.

Local inspection covered the prior two reviews, relevant existing candidate
notes, `estimator.py`, `coverage_first_pass.py`, `coverage_second_pass.py`,
the five new candidate notes, the updated Katz note, and the new script and
results. This is not an exhaustive scan of managed historical collections,
inbox/unclassified material, local notes, literature databases or citation
networks. No author was contacted and no external account content was edited.

Only the local Katz research note was updated among earlier research records;
the other files for this round are new. Earlier reviews remain historical
records; their Katz-access statement is superseded by this update. Existing
unrelated deletions and untracked files were left alone. No VDfirst or
managed VD-docs file was changed. No commit or push was made. Local saving
does not mean these documents have been submitted to VD-docs.
