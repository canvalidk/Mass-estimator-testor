# Mass-equation comparison: first review checkpoint

Date: 2026-09-08. Status: local research record awaiting transfer to VD-docs;
not submitted. This is a bounded first pass, paused for the user's review of
the comparison method. It is not a completed literature or priority search.

## What this pass establishes

The comparison now has a repeatable structure: identify the published model,
map its inputs, inspect the information its equations retain, establish a
boundary result, and record a verdict. A candidate's published prescription
and any specialization or adaptation we derive receive separate verdicts.

| Candidate | Result for the agreed full input domain | Decisive evidence | Adaptation status |
|---|---|---|---|
| Our conditional Gaussian/flat-magnitude construction | Full mathematical coverage established below | Positive normalizer and finite positive magnitude means; proper continuous mass law | Numerical implementation has finite precision and quadrature limits |
| Ordinary vector least squares | Demonstrated limitation | Undefined at zero acceleration; zero at perpendicular inputs; negative at anti-alignment | A positivity constraint alone does not supply uncertainty or guarantee an interior optimum |
| Aach–Mester–Dümbgen (2001) | Demonstrated limitation for the published statistic | Its squared dot product gives identical residuals for aligned and anti-aligned pairs | Bayesian or positive-scale adaptations require separate evaluation |
| Leonard (2011) | Demonstrated limitation for the published diffuse posterior | Centered-scatter ratios are undefined at zero–zero, with path-dependent limits | Fixed-origin, positive-scale, known-covariance adaptation unresolved |
| Mardia et al. (2013) | Demonstrated limitation for arbitrary supplied joint covariance | Its observed noise covariance is isotropic and tied to unknown scale | Fixed-alignment conditional density has full finite-reading coverage under its own model; a covariance-correct adaptation remains a serious candidate |

No earlier method satisfying the entire agreed contract has been established
in this pass. That does not establish historical priority. In particular,
Mardia should remain under investigation rather than be dismissed because its
application or prior differs.

## 1. Contract used for this checkpoint

Inputs are two finite measured vectors in a fixed common coordinate frame,
with a supplied finite symmetric positive-definite joint error covariance.
The latent true vectors obey F*=m a*, with m>0. Newton II is assumed correct;
large discrepancies do not suppress this conditional inference.

For definiteness this first pass uses the **joint Gaussian measurement model**
of the current equation. Covariance alone does not specify an error
distribution. The earlier quoted contract omitted this distinction, so it
must be made explicit before claiming equivalent inputs across models.

The required outputs are a finite strictly positive point and an uncertainty
description throughout this domain, including measured zeros, from one
mathematical prescription. Relative direction and sign must be retained.
Interior probability quantiles are allowed; a finite mass standard deviation
is not required. The method must have a stated behavior as well-resolved,
nonzero aligned data approach the uncertainty-free limit.

Additional model choices, including priors, must be declared. A different
prior or point functional is not automatically a failure: competitors need
not reproduce our exact answer. However, replacing a supplied covariance,
discarding the origin, or introducing an unknown interchannel rotation changes
the inference problem. Such changes cannot be hidden in the input mapping.

Here “coverage” means the domain on which the mathematical method provides the
requested result. It does **not** mean a 95% interval has 95% repeated-sampling
coverage. Uncertainty type, frequentist calibration, estimator accuracy,
identification, and historical priority are separate questions.

## 2. Our exact target and its domain argument

Let y=(F,a), u be a unit common direction, and f,alpha>0. Set

$$
x=(fu,\alpha u),\quad Q=(y-x)^T\Sigma^{-1}(y-x),\quad
dP=Z^{-1}e^{-Q/2}\,df\,d\alpha\,d\Omega.
$$

The base measure is flat in positive magnitudes and uniform in direction.
The point and uncertainty are

$$
\widehat m=\frac{E_P[f]}{E_P[\alpha]},\qquad M=f/\alpha,
\qquad [q_{.025}(M),q_{.975}(M)].
$$

This is neither E[M] nor the least-error reconciliation. The September 7
local steering retains this equation while removing the older managed
specification's rejection gate for the present investigation.

### Existence and continuity: our derivation

Work in fixed numerical units for the two channels. Let
lambda=lambda_min(Sigma^(-1))>0. Since ||x||²=f²+alpha² and
||y-x||² >= ||x||²/2 - ||y||²,

$$
e^{-Q/2}\leq
e^{\lambda\|y\|^2/2}e^{-\lambda(f^2+\alpha^2)/4}.
$$

Direction space has finite measure. The bound is integrable on the positive
magnitude quadrant, also after multiplication by f, alpha, or any fixed
nonnegative powers of the magnitudes. The likelihood is strictly positive.
Consequently 0<Z<infinity and both magnitude means are finite and strictly
positive. Their quotient therefore is finite and strictly positive for every
finite y and positive-definite Sigma.

On any compact neighborhood inside this input domain, covariance eigenvalues
remain bounded away from zero and infinity and ||y|| is bounded. A common
Gaussian dominating function then gives continuous integrals and continuous
point estimates by dominated convergence. Approaching a singular covariance
is outside this local argument; no finite positive bound is claimed at that
excluded boundary.

After integrating direction, let q(f,alpha) be the normalized magnitude
density. The mass density is

$$
p_M(m)=\int_0^\infty\alpha q(m\alpha,\alpha)\,d\alpha,
\qquad m>0.
$$

It is positive and proper. Its CDF is continuous and strictly increasing
from zero to one, so every interior quantile is finite and strictly positive.
The same domination proves continuity in the input; strict increase then
gives continuity of interior quantiles. The logarithmic axis singularities
are integrable under the magnitude density, so log-mass moments are finite.
These results supply one prescription through all the named configurations.

### Direction, identification, and tails

For independent isotropic channels, the direction-dependent likelihood
contains exp[u dot v], where

$$
v=fF/\sigma_F^2+\alpha a/\sigma_a^2,
\qquad
\|v\|^2=\frac{f^2\|F\|^2}{\sigma_F^4}
+\frac{\alpha^2\|a\|^2}{\sigma_a^4}
+\frac{2f\alpha F\cdot a}{\sigma_F^2\sigma_a^2}.
$$

Thus integration over unknown common direction retains the **signed** dot
product. General covariance enters the original full Q. Equal point values
in symmetric examples do not demonstrate lost direction: the uncertainty
distribution can change while symmetry fixes the point.

The density extends positively to the magnitude axes. The substitution
f=m alpha gives

$$
p_M(m)=m^{-2}\int_0^\infty f q(f,f/m)\,df
\sim K/m^2,\quad K=\int_0^\infty f q(f,0)\,df>0.
$$

Gaussian domination justifies this limit. The mass mean and second raw
moment diverge at every finite empirical input, including aligned readings.
This does not prevent finite quantiles or the chosen ratio-of-means point.

At measured double zero with independent isotropic channels, s=sigma_F/sigma_a,
the point is s and M has half-Cauchy density 2s/[pi(s²+m²)]. Its 95% interval
is approximately [0.03929s,25.4517s]. This scale comes from the instrument
uncertainties and the declared measure; the null trial does not uniquely
identify the object's mass. For nonzero aligned readings and
Sigma=epsilon² Sigma_0, the existing local concentration argument gives
both the point and interior quantiles tending to ||F||/||a|| as epsilon→0.

**Scope of the proof:** mathematical conditional inference under this
Gaussian/flat-magnitude policy, in finite spatial dimension. It does not
certify a floating-point program on every finite representable input, prove
sampling calibration, or establish the policy uniquely from Newton II.

## 3. First published-method comparisons

### Aach–Mester–Dümbgen

The [2001 paper](https://www.gretsi.fr/data/colloque/pdf/2001_001-0258_13386.pdf),
section 3.1, equation (10), introduces two noisy vectors and a positive
multiplicative factor. With A=||x1||², B=x1 dot x2 and C=||x2||², equation (18) is

$$
D^2=\tfrac12[A+C-\sqrt{(A-C)^2+4B^2}].
$$

The iid Gaussian, common-variance model yields a collinearity statistic and
change decision. Section 3.2's approximate chi-square law concerns that
statistic, not uncertainty for the factor. The relative sign disappears in B².

**Our counterexample:** x1=(1,0), x2=±(1,0), joint covariance I4, gives D²=0
for both signs. In contrast, profiling the positive-scale Gaussian residual
for the negative sign gives

$$
R(m)=\min_s\{\|e-s\|^2+\|-e-ms\|^2\}
=\frac{(1+m)^2}{1+m^2}>1,\quad m>0.
$$

Its infimum 1 is not attained at any finite positive m. Therefore simply
constraining TLS to positive scale also fails to supply an everywhere finite
positive maximum-likelihood estimate. This is our deduction, not the paper's
claim. Verdict: demonstrated limitation; alternative Bayesian constructions
remain unevaluated.

### Leonard

[Leonard (2011)](https://arxiv.org/pdf/1202.0957), equations (13)–(21),
uses centered scatter through r=S12/sqrt(S11 S22) and l=sqrt(S22/S11).
Its published diffuse posterior models an unknown intercept, unrestricted
slope, and unknown channel error variances. It retains signed centered
correlation on its native domain, but zero–zero makes r and l undefined.

**Our stronger counterexample:** with fixed measurement covariance I6, take

$$
a_\epsilon=\epsilon(-1,0,1),\quad
F_\epsilon=(c\epsilon/\sqrt3)(1,-2,1),\quad c>0.
$$

Every path approaches the same zero–zero input, yet r=0 and l=c throughout.
Equation (19) then gives distinct posterior scale families for different c.
Even positive-conditioning leaves their medians at c, using reciprocal
invariance (51). No unique continuous zero–zero extension exists for this
diffuse prescription. Verdict: demonstrated limitation, not impossibility
of a fixed-origin known-covariance adaptation.

The [Leonard evidence note](mass_candidate_leonard_evidence_2026-09-08.md)
records the algebra and a separate centering example that reverses raw
vector-direction sign without changing its posterior.

### Mardia

[Mardia et al. (2013)](https://arxiv.org/pdf/1312.1840), equations (3)–(4),
provides a positive-scale density. Our specialization fixes one match,
rotation I, translation zero, and their noise parameter s=sigma_c>0.
For gamma prior shape k>0 and rate b>0, it becomes

$$
h(c)=c^{r-1}e^{-bc}e^{-\|F-ca\|^2/(4s^2)},\quad r=d/2+k.
$$

**Our proof:** 0<h(c)<=c^(r−1)e^(−bc), also after multiplication by any
nonnegative power of c. Thus the normalized density has finite positive
mean, finite interior quantiles, and continuity through every finite reading.
It retains the signed F dot a. At a=0 it is Gamma(r,b).

However, equation (1) implies independent observed covariances s²I and
s²I/c², tied to unknown c; this is not our fixed supplied joint covariance.
Allowing free rotation also removes a single pair's relative angle.

Verdict: full finite-reading coverage for the stated fixed-alignment
conditional model; demonstrated limitation against the full input contract.
A covariance-correct adaptation remains unresolved. The
[Mardia evidence note](mass_candidate_mardia_evidence_2026-09-08.md) records
the parameterization and numerical boundary caveats.

## 4. Reproducible first numerical matrix

New script: `03_trace_and_evaluator/mass_estimator/coverage_first_pass.py`.
Output: `03_trace_and_evaluator/mass_estimator/results/coverage_first_pass_2026-09-08.json`.
The existing estimator was used without modification or its historical gate.

Six configurations were repeated under four covariance settings: small
isotropic errors, large isotropic errors, anisotropic errors, and a full
positive-definite covariance including directional and interchannel
correlations. The JSON records every matrix and input, source-code hashes,
point, mass quantiles, OLS result, and quadrature refinement differences.

An additional 34 points sample an angle path from alignment through
perpendicularity to anti-alignment and signed paths crossing each zero.
These samples inspect behavior; the existence and continuity argument above
establishes the mathematical domain claim.

Example subset: 2D, independent isotropic errors sigma_F=sigma_a=0.25;
acceleration is (1,0) except where explicitly zero. Units are fixed so the
ordinary aligned quotient is 2.

| Measured F; measured a | OLS | Our point | Conditional 95% mass interval |
|---|---:|---:|---:|
| (0,0); (0,0) | Undefined | 1.000 | [0.0393, 25.45] |
| (0,0); (1,0) | 0 | 0.2070 | [0.00813, 0.7276] |
| (2,0); (0,0) | Undefined | 9.946 | [3.395, 253.2] |
| (2,0); (1,0) | 2 | 2.000 | [1.221, 4.063] |
| (0,2); (1,0) | 0 | 8.565 | [2.887, 218.7] |
| (-2,0); (1,0) | -2 | 34.93 | [9.534, 1315] |

All 24 cases and 34 path points returned finite positive points and interior
quantiles. An independent analytic half-Cauchy null oracle passed. A direct
evaluation of Aach's equation reproduced the sign counterexample.

The initial 128/1024 versus 256/2048 direction/ratio-grid comparison exposed
a 1.21% relative difference in the anti-aligned upper quantile. That case was
refined to 384/8192 versus 512/16384. After targeted refinement, the largest
relative difference over the 24 cases was 3.62e−10 for the point and
2.74e−4 (0.0274%) for any reported quantile. These are observed differences,
not certified error bounds. Path points were sampled at one resolution.

No Mardia numerical result is put alongside this table as if it consumed the
same covariance. Leonard and Aach have algebraic failures; a large grid is
unnecessary to establish those failures. This table demonstrates the case
format and our reference behavior, not an accuracy ranking between competitors.

## 5. Review point and next bounded pass

The main methodological decision for review is to keep **published coverage,
adapted coverage, and exact equation identity** distinct. This avoids making
historical priority depend on adaptations we derived ourselves, while keeping
mathematically promising competitors in the comparison.

The most useful next bounded pass is a covariance-correct positive-scale
Bayesian comparator: specify its latent measure and point functional, show
it consumes the same joint covariance, and determine whether an earlier
publication supplies that prescription. Mardia and its cited methodological
ancestors are a concrete starting point. This is more informative than
adding many loosely related sources before reviewing the criteria.

Butler, Liseo, Katz, and the remaining supplied references have not received
a new equation audit in this pass. The September 7 scalar-equivalence checks
remain existing local evidence, not newly reverified results here. The scope
is deliberately paused before expanding this list.

## 6. Provenance, managed/local relationship, and handoff

Managed repository: `canvalidk/VD-docs`. Main resolved through the installed
GitHub connector to commit **adfe14bbaafb29856f0a30da088495b7c9d483d2**.
The catalog and task sources were fetched at this pinned commit:

- `catalog.yaml`: selected mass-estimator inventory entries; no full inventory audit.
- `Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`:
  sections 7–12 for the construction, and the rejection/remaining-policy
  passages for the relationship to later local steering.
- `Newton-analysis/_dscn_mass_estimator/mass_estimator_originality_research_2026-09-05.md`:
  prior findings, tail argument, and coverage limitations.

The specification is the managed construction basis; the September 5 research
is an earlier research record. The catalog does not list the named September
7 local steering, practical-formula, or literature-recheck records at this
snapshot. They were inspected as active local/unsubmitted developments rather
than assumed stale. The steering explicitly overrides the older rejection
framing for this investigation; this is a scope conflict already recorded
locally, not a managed rewrite.

Local sources read: `AGENTS.md`; the mass-estimator laboratory's `README.md`,
`estimator.py`, and `experiments.py`; and these existing `pass to VD-docs/` records:
`../Mass estimator/mass_estimator_law_assumed_correct_steering_2026-09-07.md`,
`../Mass estimator/mass_estimator_practical_limits_and_uncertainty_2026-09-07.md`,
`mass_estimator_literature_recheck_2026-09-07.md`.
The arrow-data audit was incidentally read in the initial batch but supplies
no finding in this report.

Selection used scoped filenames containing mass/ratio/equation, implementation
symbols concerning covariance/integration/moments, relevant catalog paths,
and direct reading of the three named primary papers. No broad new literature
keyword search, historical-folder sweep, inbox/unclassified audit, or
VDfirst inspection was performed. Source PDF text was checked at the cited
equations; web screenshots failed, and no visual verification is claimed.

New local work consists of this review, the Leonard and Mardia evidence notes,
the comparison script, and its JSON output. Existing source documents and
estimator code were not altered. No submission, managed reorganization,
repository commit, push, or contact with authors was performed.
