# Mass-equation comparison: six more candidates, second review

Date: 2026-09-08. Status: local research record awaiting transfer to VD-docs;
not submitted. The user authorized approximately five or six more candidates
and then a stop to discuss the results and diminishing returns. **Six source
groups were examined; work is paused at that boundary.** Katz and the Chang
corroborating manuscript count as one group, not two independent candidates.

## Result of this pass

The new sources establish strong scalar precedents and clarify the limits of
least-squares approaches. None has been established as a published replacement
meeting the complete two-vector contract. This is not evidence of first
discovery. Katz's original full text remains unavailable, so its complete
scope is unresolved; the formula attributed to it is independently corroborated.

| Candidate | Strongest supported result | Verdict against the complete requested problem | What this adds |
|---|---|---|---|
| Butler (2004) | Direct force–acceleration mass estimation, with explicit treatment of zero-excitation phases | Demonstrated limitation: least-squares point can be negative/zero or unidentified; native zero handling retains history or stops adaptation | Distinguishes maintaining a controller estimate from inferring mass and uncertainty from the current measured pair |
| Bartel–Stoudt–Possolo (2016) | Errors-in-variables fitting and Monte Carlo uncertainty, including specified correlations | Demonstrated limitation: the explicit linear specialization can have no finite positive optimum; zero–zero leaves slopes tied | Adding uncertainty propagation does not itself solve an optimizer's existence problem |
| Liseo (2003) | Positive-conditioned equation (12) exactly matches our unknown-sign scalar mass distribution under its stated Gaussian/flat-means model | Full scalar coverage established for the specialization; full vector prescription and our particular point not established by this source | Firm scalar distributional precedent; correlated scalar extension derived and checked separately |
| Katz (1961), corroborated in Chang et al. | All finite scalar readings give a finite positive normal-mean estimate | Full coverage of the corroborated scalar rule; full two-vector answer not supplied by that rule; original-paper scope unresolved | Earlier positive estimate at zero; the known-direction quotient is our explicit composition of two such rules |
| Roe–Woodroofe (2000) | Explicit restricted-normal posterior and boundary credible intervals | Full scalar uncertainty coverage; native model does not specify the two-vector relationship | Published uncertainty companion to the same positivity mechanism; interval branches need not be ad hoc repairs |
| von Luxburg–Franz (2007) | Covariance-aware Fieller confidence sets, including unbounded sets | Demonstrated limitation of the native point-plus-set package: raw ratio fails zero/positivity; positive-restricted set can be empty | Strong scalar uncertainty method with a different probability meaning from a posterior |

Together with Aach, Leonard, and Mardia in the
[first review](mass_equation_comparison_review_2026-09-08.md), this is nine
candidate groups. It is not nine independent failed versions of our equation:
several are successful building blocks or uncertainty constructions.

## Contract and comparison discipline retained

The input is two finite measured vectors in the same fixed coordinate frame,
with a known finite symmetric positive-definite joint **Gaussian** error
covariance. The true vectors satisfy F*=m a*, m>0, by premise. The desired
output is a finite positive point plus stated uncertainty throughout this
domain, retaining direction and relative signs, using one mathematical
prescription, including measured zeros and adverse angles.

Coverage of that input domain, repeated-sampling confidence coverage, and
posterior credibility are different properties. Finite interior quantiles
are permitted even where a mass mean or standard deviation does not exist.
Finite endpoints are not required of every confidence set. A finite point
does not by itself establish identification.

Our target remains the Gaussian likelihood on (fu,alpha u), f,alpha>0,
under flat df d(alpha) and uniform common direction, with point E[f]/E[alpha]
and the induced distribution of M=f/alpha. No rejection gate is applied.
Different priors or point functionals remain eligible competitors if declared.
Published prescriptions, our specializations, and exact equation identity
are assessed separately throughout.

## 1. Butler: physical precedent, with operational zero handling

[Butler's published disclosure](https://patents.google.com/patent/EP1455231A2/en)
uses recursive least squares for stage mass. Section 3.1 gives

$$
A^TA\widehat m=A^TF.
$$

Our mapping of sample indices to vector coordinates yields
m-hat=(a dot F)/(a dot a) when a is nonzero. With joint covariance I4,
a=(1,0) and F=(-1,0) gives -1; F=(0,1) gives zero. At a=0 the normal equation
does not select a unique mass. The path a=(t,0), F=(1,0) gives 1/t.

Section 2.4 explicitly disables adaptation at zero setpoint acceleration or
limits its gain. Retaining a previously initialized mass can remain finite;
that uses history rather than solving the specified current-pair problem.
The audited prescription does not provide the required joint-covariance
mass uncertainty. This does not rule out an augmented controller.

**Verdict:** demonstrated limitation. A direct physical precedent has been
verified, and the zero issue was explicitly discussed. Further tests of the
same least-squares singularity would add little. See the
[Butler evidence note](mass_candidate_butler_evidence_2026-09-08.md).

## 2. Bartel: EIV uncertainty does not guarantee a positive interior optimum

[Bartel, Stoudt and Possolo's NIST-hosted manuscript](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=919758),
equation (3), minimizes weighted force and transducer-response residuals over
calibration coefficients and latent force values. Section 5.3 propagates
uncertainty through repeated fits, including specified systematic dependence.
It would be inaccurate to say it ignores all correlation.

Our explicit specialization sets polynomial degree one and intercept zero,
maps pairs to vector coordinates, and uses unit error covariance. Eliminating
the latent vector gives

$$
Q(m)=\frac{\|F-ma\|^2}{1+m^2}.
$$

At a=(1,0), F=(-1,0), its unconstrained optimum is negative. Adding m>0 gives
Q(m)=1+2m/(1+m²)>1, with unattained infimum 1 at zero or infinity. At both
measured vectors zero, every slope is tied. Resampling the optimizer does
not supply a unique positive solution to either problem. Bounds or prior
regularization would require a separate specification.

**Verdict:** demonstrated limitation for the requested contract, with the
specialization stated explicitly. See the
[Bartel evidence note](mass_candidate_bartel_evidence_2026-09-08.md).

## 3. Liseo: a precise scalar distributional identity

[Liseo (2003)](https://www.researchgate.net/publication/5182203_Bayesian_and_conditional_frequentist_analyses_of_the_Fieller%27s_problem_A_critical_review),
equations (10), (12)–(13), gives the ratio posterior associated with a prior
flat in two signed Gaussian means. The paper's inspected model uses equal
independent errors. In ratio/denominator coordinates (m,t), its Jacobian is |t|.

Our positive-conditioning specialization is

$$
p_+(m)\propto\mathbf1_{m>0}\int_{-\infty}^{\infty}
|t|\,\phi_2((mt,t);(F,a),C)\,dt.
$$

Writing the same-sign means as (u f,u alpha), u=±1, gives exactly our scalar
base measure and ratio distribution. This proves an identity throughout
that scalar model, not only at zero. The positive conditioning is ours;
the displayed publication does not claim our ratio-of-magnitude-means point.

Replacing the paper's covariance by any known positive-definite 2×2 C gives
a direct correlated scalar extension. We derived and tested it; we do not
attribute that extension to the paper.

**Verdict:** exact scalar distributional precedent after the stated
conditioning; full vector coverage remains unresolved from this source.
See the [Liseo evidence note](mass_candidate_liseo_evidence_2026-09-08.md).

## 4. Katz: the positive scalar boundary mechanism survives

The [Chang–Shinozaki–Strawderman manuscript](https://m.sc.niigata-u.ac.jp/~hirukawa/seminar/niigata2017_program/Chang.pdf),
equation (1.4), explicitly attributes to Katz the positive-normal-mean rule

$$
h_\sigma(x)=x+\sigma\frac{\phi(x/\sigma)}{\Phi(x/\sigma)}.
$$

It is the mean of the Gaussian likelihood normalized over a nonnegative
parameter under a flat prior. Hence it is finite and strictly positive for
every finite reading, with h_sigma(0)=sigma sqrt(2/pi). The proper posterior
also supplies finite interior quantiles.

Our known-direction, independent-channel construction is
h_sigmaF(F)/h_sigmaa(a). It equals our restricted point and stays finite and
positive through both zeros. This quotient is our composition of the rules,
not a mass formula attributed to Katz. A known-direction correlated model
requires joint truncation; multiplying the marginal rules discards correlation.

**Verdict:** full scalar coverage for the corroborated rule; unknown-direction
two-vector coverage not supplied by that rule. The original
[Katz article](https://doi.org/10.1214/aoms/1177705146) remained unreadable, so
its complete scope is unresolved. See the
[Katz evidence note](mass_candidate_katz_evidence_2026-09-08.md).

## 5. Roe–Woodroofe: boundary uncertainty is already explicit

[Roe and Woodroofe](https://arxiv.org/pdf/hep-ex/0007048), equations (1)–(6),
publishes, for unit error scale,

$$
p(\theta\mid x)=\frac{\phi(x-\theta)}{\Phi(x)},\quad\theta\ge0,
$$

and shortest posterior credible intervals. The posterior remains proper for
every finite x. Its mean is the Katz rule; at x=0 the posterior is half-normal.
The paper's interval changes from a boundary interval to an interior interval
as the reading grows. Both are solutions of the same density-threshold
prescription, not repairs to a failed zero calculation.

A zero lower interval endpoint does not violate the requirement of a positive
point estimate. Nor should its posterior probability be represented as an
unqualified guarantee of nominal repeated-sampling coverage; the paper
distinguishes these properties.

**Verdict:** full scalar boundary uncertainty coverage under the stated
Bayesian interpretation. Its native scalar model does not specify a noisy
denominator, unknown spatial direction, or full joint vector covariance.
See the [Roe–Woodroofe evidence note](mass_candidate_roe_woodroofe_evidence_2026-09-08.md).

## 6. Fieller geometry: successful uncertainty, different output

[von Luxburg and Franz](https://arxiv.org/pdf/0711.0198), equation (2.1),
uses the residual pivot with the covariance cross-term. Our known-covariance
Gaussian specialization gives the positive-restricted confidence set

$$
\mathcal C_+=\left\{m>0:
(F-ma)^2\le q^2(C_{FF}-2mC_{Fa}+m^2C_{aa})\right\}.
$$

Here q is a normal critical value; the paper's estimated-covariance version
uses a Student-t reference. The single defining inequality handles measured
zeros without dividing by a. At double zero it includes every positive m.
At (F,a)=(-4,4), C=I and q≈1.96, it is empty because
16(1+m)²>q²(1+m²) for all m>0.

An empty or unbounded confidence set is not an invalid statistical output.
Nevertheless, the native point F/a fails zeros and positivity, and the set
does not itself supply the positive conditional distribution requested here.

**Verdict:** demonstrated limitation of the point-plus-set package against
the full contract; successful scalar confidence inference in its own sense.
See the [Fieller evidence note](mass_candidate_fieller_geometry_evidence_2026-09-08.md).

## Independent numerical checks and what they mean

New script: `03_trace_and_evaluator/mass_estimator/coverage_second_pass.py`.
Output: `03_trace_and_evaluator/mass_estimator/results/coverage_second_pass_2026-09-08.json`.
Existing estimator and previous numerical records were not changed.

Seven scalar configurations, including both zeros, each one-zero case,
both signs agreeing, and opposing signs, were repeated at correlations
−0.6, 0, and +0.6: **21 correlated scalar checks**. Their covariance scales
were sigma_F=0.7 and sigma_a=1.3. These are scalar restrictions of the
comparison contract; perpendicular vectors are not a scalar configuration.

The independent oracle integrates the original Cartesian Gaussian posterior
over same-sign quadrants, integrating conditional force analytically and
signed acceleration numerically. It does not reuse the reference estimator's
ratio-angle/radial integral. It truncates acceleration only beyond the observed
absolute reading plus 12 standard deviations; no mass cutoff is imposed.
Results with 300 and 600 quadrature nodes were compared.

| Check | Largest observed discrepancy |
|---|---:|
| Our point versus Cartesian conditional-normal oracle, 21 cases | 2.19e−13 relative |
| Oracle CDF at reported 2.5%, 50%, 97.5% mass quantiles | 2.54e−8 absolute probability |
| Oracle CDF change between 300 and 600 nodes | 1.54e−12 absolute probability |
| Liseo independent-error density versus our scalar density, 7 cases | 4.89e−15 relative |
| Quotient of Katz rules versus our known-direction point, 7 cases | 3.91e−14 relative |
| Correlated zero–zero analytic Cauchy check | 2.96e−9 absolute probability |

All assertions passed. These calculations validate translations and numerical
implementation on the named inputs; the change-of-variables and existence
arguments establish the mathematical identities and domain statements.

A specific correlation witness prevents overextending the quotient formula:

$$
F=1.4,\quad a=1.3,\quad
C=\begin{pmatrix}0.49&0.546\\0.546&1.69\end{pmatrix},
\quad\text{known positive direction}.
$$

Joint inference gives point **0.9125851**. The independent Katz quotient gives
**0.8594843** because it omits the supplied correlation. The joint result was
checked against the separate Cartesian integral. A correlated extension can
work; the independent factorization cannot simply be reused unchanged.

The script also reproduces known-covariance scalar Fieller cases and the
anti-aligned TLS profile witness. These are checks of our explicit
specializations, not implementations of the complete original applications.

## What remains distinctive enough to investigate

The scalar-to-vector step is substantive. In one spatial dimension,
same-sign true means occupy an open region of the Cartesian two-channel
space. Positive conditioning is therefore ordinary conditioning on an event
of positive probability. In d>1, exactly codirectional pairs have dimension
d+1 inside dimension 2d, hence probability zero under a nonsingular Cartesian
Gaussian. Ordinary event-conditioning does not select a unique measure
there; a model or limiting prescription is needed. Our df d(alpha) dOmega
measure is one stated choice.

Consequently, exact scalar precedents do not by themselves establish the
full vector formula. They also do not establish its novelty. A published
latent-direction or latent-vector scale model could settle the broader
coverage claim while using a different measure or point estimate.

The central open comparison still has three independent questions:

1. Does a published positive-scale model consume the same supplied joint
   Gaussian covariance and provide an everywhere defined point and uncertainty?
2. Does a publication use this particular flat-magnitude/common-direction law?
3. Does it select E[f]/E[alpha], or an equivalent point functional?

A method answering the first question succeeds under the agreed contract even
if it answers the other two differently. No extra uniqueness requirement has
been added after inspecting these candidates.

## Diminishing returns: assessment for discussion

**There are diminishing returns in continuing indiscriminately through generic
scalar-ratio and least-squares sources.** This pass adds new evidence, but the
mechanisms are starting to repeat: Katz and Roe–Woodroofe share one restricted
normal posterior; Liseo pins down the scalar ratio law; Butler repeats the
OLS boundary; Bartel's linear specialization repeats the positive TLS boundary.

The four original list entries not newly audited by this two-pass method are
Marsaglia, Kim, Little, and Arratia–Goldstein–Kochman. Their roles in the supplied
list are ratio-distribution mathematics, folded-normal ratios, ratios of
posterior expectations, and size bias. The expected return from those entries
is higher for mathematical ancestry or a specific citation than for finding
a full noisy-vector positive-scale method. This is a prioritization judgment,
not an assertion about their unread full contents.

The higher-value unresolved direction is targeted vector/scale inference with
declared covariance. Mardia's fixed-alignment model from the first pass remains
a useful lead, with its scale-dependent noise model explicitly unresolved
against our fixed-covariance contract. A carefully specified covariance-correct
adaptation would clarify what to search for, but a construction we derive is
not evidence that an earlier publication already supplied it.

This pass stops before undertaking that work. There has been no exhaustive
database search, citation-network audit, or proof of historical priority.

## Provenance and custody

Managed repository main was checked through the installed GitHub connector
again and remains **canvalidk/VD-docs@adfe14bbaafb29856f0a30da088495b7c9d483d2**.
The already-read sources at this commit remain the coherent target snapshot:

- `catalog.yaml`, relevant mass-estimator inventory selection;
- `Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`,
  the construction and policy passages;
- `Newton-analysis/_dscn_mass_estimator/mass_estimator_originality_research_2026-09-05.md`,
  historical research context, not an exhaustive priority finding.

The active local September 7 law-assumed-correct steering remains authoritative
for this investigation over the managed specification's older rejection gate.
That local development and the practical-formula/literature-recheck records
were not listed in the selected managed snapshot, as recorded in the first
review; no managed/local reconciliation was performed here.

Local scope: root `AGENTS.md`, first review, September 7 literature recheck,
the already established steering/practical interpretation, estimator source,
existing literature-equivalence script, and the relevant beginning of
`test_estimator.py`. New evidence consists of the six linked candidate notes,
this review, the comparison script, and its JSON. Source-code hashes are saved
with the numerical output. Existing files, entry production, and VDfirst
were not modified.

Selection followed the six named primary references rather than a new broad
keyword search. Access follow-ups were limited to DOI/publisher/author/repository
routes for those sources. Katz's original full text was unavailable; its
formula is corroborated in the later original research manuscript. Liseo was
read from author-upload full text. Bartel was read from the NIST-hosted
manuscript, not assumed identical in pagination to the journal version.
Relevant source text and equations were inspected; PDF screenshot failures
and image-only untranscribed details are recorded in individual notes.

Managed historical collections, inbox/unclassified material, and unrelated
local notes were not exhaustively searched. No managed files were changed,
no source collection was copied, and no author was contacted. Saving these
files locally does not submit them to VD-docs. No commit or push was made.
