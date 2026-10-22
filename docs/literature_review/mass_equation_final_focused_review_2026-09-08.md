# Mass equation: final focused predecessor review

Date: 2026-09-08. Status: local active research record awaiting transfer to VD-docs; **not submitted**. This completes the user's authorized final focused check and contribution draft. Literature expansion stops here pending discussion or specialist review. Mardia was not reopened.

## Result and stopping decision

No inspected publication in this final round was established as giving the complete target equation and algorithm. The round checked six source groups: four older errors-in-variables/model-fitting sources and two sources for the point functional's ancestry. The primary-text gaps for Lindley–El-Sayyad and Dellaportas–Stephens were closed through browser access.

The useful result is a more exact account of what has precedents. Lindley–El-Sayyad's generic hierarchy already contains an important target subcase after an alternative hyperprior is chosen. A conditional two-dimensional specialization of Werman–Keren's integration measure also equals the target measure. Dellaportas–Stephens computes ratio uncertainty by transforming paired posterior draws. Little uses a quotient of posterior expectations, and Arratia–Goldstein–Kochman gives the general weighting identity. Gull already integrates latent Gaussian coordinates analytically.

These findings support a qualified claim for the **particular combined construction, analysis and application**, rather than separate first-ever claims for its prior, algebraic operation, general Bayesian algorithm, or force–acceleration mass estimation. The [contribution statement draft](../Mass%20estimator/mass_estimator_contribution_statement_draft_2026-09-08.md) makes that claim and preserves the user's account of human and AI contributions.

Further broad searches for generic least-squares or Bayesian errors-in-variables methods now have diminishing value: their relevance is established. The next useful review is of the exact measure, readout and resulting benefits by a specialist. This recommendation is a stopping judgment for the present work, not a claim that the literature is exhausted. No specialist was contacted.

## Six source groups and supported verdicts

| Source and primary text | Equation-level finding | Verdict |
|---|---|---|
| [Lindley–El-Sayyad, 1968](https://doi.org/10.1111/j.2517-6161.1968.tb01519.x), pp.190–202 | Eq.(1) already fixes a through-origin relation. Eq.(9) is an exact latent-Gaussian hierarchy with generic hyperprior. A flat latent-variance measure produces our radial prior for dimension greater than two; their explicit choices use log-uniform variance. Their displayed point calculations are large-sample slope approximations. | Exact containment of the target independent homogeneous-error subcase under specified alternative choices is established by our derivation. The authors' instantiated prior and point differ; no complete target publication identified. |
| [Gull, 1989](https://bayes.wustl.edu/sfg/line.pdf), pp.511–518 | Gaussian latent abscissae and range hyperpriors; analytic nuisance integration; an approximate optimization point on p.517. | Established methodological predecessor. Displayed prior and point differ; target combination not identified. |
| [Dellaportas–Stephens, 1995](https://www.jstor.org/stable/2533007), pp.1085–1095 | Joint latent posterior Eq.(4), Gibbs/full conditionals Eq.(7), and a ratio computed from paired posterior draws on p.1090. | General adaptable posterior algorithm and explicit ratio-uncertainty predecessor. Published examples do not instantiate target measure plus quotient-of-expectations readout. |
| [Werman–Keren, 2001](https://www.cs.haifa.ac.il/~dkeren/mypapers/werman-keren-fitting.pdf), pp.528–534 | Integrates latent points using arc length and selects a MAP model. Positive origin-fixed lines with an added uniform-angle convention give our measure in two dimensions. | Conditional measure identity established; literal publication of that specialization is unresolved. MAP point does not equal the target readout in general. |
| [Little, 2012](https://www.scb.se/contentassets/ca21efb41fee47d293bbee5bf7be7fb3/calibrated-bayes-an-alternative-inferential-paradigm-for-official-statistics.pdf), Example 3, pp.318–319 | Explicit ratio of posterior expectations, used as an approximation to a posterior expectation of a ratio in a survey example. | Earlier algebraic usage established. This approximation is not the target interpretation, whose ordinary posterior mass mean diverges. |
| [Arratia–Goldstein–Kochman, 2013 v1](https://arxiv.org/pdf/1308.2729v1), §2.3, Eqs.(14)–(15) | Exact change of joint probability measure by one nonnegative coordinate. | Exact general ancestry of acceleration-weighted interpretation established. No target physical posterior or readout application identified in the inspected sections. |

Detailed evidence, including exact priors, qualification of approximations, access provenance and page anchors:

- [Lindley–El-Sayyad and Dellaportas–Stephens](mass_predecessor_lindley_dellaportas_focused_check_2026-09-08.md).
- [Gull and Werman–Keren](mass_predecessor_gull_werman_focused_check_2026-09-08.md).
- [Point-functional ancestry and loss analysis](mass_point_functional_ancestry_evidence_2026-09-08.md).

“Not identified” refers to these inspections. It is not a conclusion about unread papers, a framework's possible specializations, or another author's unpublished work. These six groups supplement the prior fourteen audited groups; they are not twenty independent complete competitors. Little and Arratia–Goldstein–Kochman were already present as background entries in the user's original list, and received this more focused check.

## Target and exact equivalence criteria

The target assumes finite measured vectors $y=(F,a)$, a known positive-definite joint Gaussian error covariance $\Sigma$, and compatible true vectors $(fu,\alpha u)$, where $f,\alpha>0$ and $u\in S^{d-1}$. It normalizes the Gaussian likelihood against $df\,d\alpha\,d\Omega$, reports $E[f]/E[\alpha]$, and uses the posterior law of $m=f/\alpha$ for uncertainty. No discrepancy rejection gate is applied under the current steering.

The change of variables $v=\alpha u$, $f=m\alpha$ gives

$$
df\,d\alpha\,d\Omega=\|v\|^{2-d}\,dm\,d^dv,
\qquad
\widehat m=\frac{E[m\|v\|]}{E[\|v\|]}.
$$

Consequently the Gaussian observation likelihood, the complete prior after marginalization, the point functional, and the uncertainty law must all be compared. Sharing only a likelihood does not establish an identical estimator. Conversely, a conditional Gaussian latent prior is not by itself evidence of a different marginal prior.

### The variance-mixture connection

For $d>2$ and $v\ne0$, direct substitution $s=\|v\|^2/(2t)$ gives

$$
\int_0^\infty N_d(v;0,tI)\,dt
=\frac{\Gamma(d/2-1)}{2\pi^{d/2}}\|v\|^{2-d}.
$$

In three dimensions the result is $1/(2\pi\|v\|)$. The singularity at zero is locally integrable in Cartesian volume. This is an identity of improper measures, not a proper prior distribution. A flat variance measure $dt$, a log-uniform measure $dt/t$, and a flat standard-deviation measure are different choices. The flat-variance integral diverges in two dimensions; the target measure there is directly flat in $(m,v)$.

Lindley–El-Sayyad Eq.(9), with known homogeneous independent measurement variances and $n=d>2$ coordinate pairs, becomes the target subcase when its generic hyperprior is selected as flat positive slope times flat variance. Its explicit log-uniform variance choices instead yield a different radial power. This exact containment is our calculation from the published family; it is not a claim that the authors made the same choice. Their discussion on p.201 also already considers flat incidental-coordinate priors and alternatives based on radial distance.

### Two-dimensional arc-length identity and different point rules

Represent the coordinates as $d$ points $(a_i,F_i)$ on a positive through-origin line. In fixed nondimensional units, choose uniform line angle $\theta=\arctan m$, and use latent arc length $ds_i=\sqrt{1+m^2}\,dv_i$. Then

$$
d\theta\prod_i ds_i=(1+m^2)^{d/2-1}\,dm\,d^dv.
$$

For $d=2$, this equals the target measure. Werman–Keren supplies the arc-length integration; fixing the origin, restricting the sign and specifying this angle measure are additions in our comparison. Its native affine-line prior measure does not uniquely prescribe this restriction.

Even with that same measure, the point decision remains substantive. At both-zero readings with equal independent isotropic errors, the target mass law is half-Cauchy and its magnitude-mean quotient is 1. The posterior line angle is uniform, so MAP in angle ties at every positive slope. Identical posterior measures do not imply identical reported points; MAP itself depends on the reference parameterization.

## Point-functional interpretation and a correction

Let $P$ be the target posterior, $A=E_P[\alpha]$, and $B=E_P[f]$. Define the acceleration-weighted joint law $dQ=\alpha\,dP/A$. The established coordinate-weighting identity gives exactly

$$
E_Q[m]=B/A.
$$

This is an exact interpretation, not an approximation to the divergent $E_P[m]$. The mass quantiles reported by the method are from $P$, so they must not silently be relabelled as quantiles under $Q$.

A tempting decision-theoretic argument requires correction: the expected weighted squared error $E_P[\alpha(d-m)^2]$ is infinite for every finite action $d$. Near $\alpha=0$, its force term entails an integral proportional to $\int d\alpha/\alpha$. Dropping this divergent constant is not a finite-risk minimization proof. The pointwise loss difference is integrable, however, and its expected regret relative to $m_*=B/A$ is $A(d-m_*)^2$.

There is also a finite nonnegative loss with exactly the same optimum:

$$
L(d;f,\alpha)=f\log\!\frac{f}{d\alpha}-f+d\alpha.
$$

Its logarithm has a dimensionless argument. Gaussian tail bounds and integrability of logarithms at zero give finite expected loss for every $d>0$. The risk derivative is $A-B/d$, so its unique minimizer is $B/A$. This is an interpretation derived in this audit, using a general known divergence form; no historical novelty claim for the loss is made. It does not prove that the physical law uniquely selects this loss or the posterior measure. The [composition theorem](../Mass%20estimator/mass_readout_composition_uniqueness_2026-09-06.md) remains a separate characterization under its explicitly declared assumptions.

## What was verified mathematically and computationally

For the target measure, the Gaussian likelihood is bounded by $K\exp[-c(f^2+\alpha^2)]$ for finite data and SPD covariance. Thus its normalizer and all nonnegative amplitude moments are finite, with a strictly positive normalizer and first moments. This proves finite positive output. Local domination in both data and the open SPD domain gives continuity. The mass density is positive on $(0,\infty)$, so all interior quantiles are finite, positive and continuous.

Substitution $\alpha=f/m$ into the ratio density gives

$$
p_M(m)\sim C/m^2,\qquad
C=Z^{-1}\int_{S^{d-1}}\int_0^\infty fL(fu,0)\,df\,d\Omega>0.
$$

The divergent ordinary mass mean is therefore a property of the stipulated measure throughout this Gaussian domain, not a numerical defect or a zero-only exception.

An independent mathematical check also confirmed the radial algorithm. Set $f=s_F r\sin\theta$, $\alpha=s_a r\cos\theta$, with fixed positive reference scales. If $w=(s_F\sin\theta\,u,s_a\cos\theta\,u)$, $\kappa=w^T\Sigma^{-1}w$, $b=w^T\Sigma^{-1}y$, and

$$
J_k(\kappa,b)=\int_0^\infty r^k e^{-\kappa r^2/2+br}\,dr,
$$

then the angle density is proportional to $\int J_1\,d\Omega$, and

$$
\widehat m=\frac{s_F}{s_a}
\frac{\int_0^{\pi/2}\int_{S^{d-1}}\sin\theta\,J_2\,d\Omega\,d\theta}
{\int_0^{\pi/2}\int_{S^{d-1}}\cos\theta\,J_2\,d\Omega\,d\theta}.
$$

The Jacobian is $s_Fs_a r$, with no extra angular power. These expressions describe the existing analytic radial integration and direction/ratio quadrature. They were checked analytically in this round; the estimator code was not changed and previous numerical sweeps were not repeated. The [third review](mass_equation_comparison_third_review_2026-09-08.md) records the existing fixtures, refinement evidence and independent two-dimensional identity checks. Mathematical coverage does not establish uniform floating-point accuracy or fixed-grid convergence for arbitrarily ill-conditioned inputs. The current reference code supports dimensions 1–3; the measure arguments apply in any fixed finite dimension.

## Search scope, source access and provenance

Selection followed the older Bayesian errors-in-variables citation trail from the already-audited Huard–Mailhot, Kelly and Mantz sources, and the two named quotient/size-bias background sources. The aim was to inspect complete priors, transformations, point decisions and computational steps, not to add generic competitors indefinitely. Searches used the exact source titles and concepts including positive through-origin scale, latent variance/range priors, uniform latent coordinates or arc length, quotient of posterior expectations, and coordinate size bias.

A small supplementary application screen used the queries `"mass estimation" "ratio of" "posterior"`, `"force" "acceleration" "ratio of expectations"`, and `"Bayesian" "inertial mass" "estimation"`. It identified no new exact match. A MEMS parameter-extraction article surfaced but its body was not inspected; it is not counted as an audited or failed competitor. An earlier Werman–Keren conference precursor was noticed bibliographically but not audited separately. No inference of non-equivalence is drawn from either unread source. No Scopus or Web of Science search is claimed.

Lindley–El-Sayyad's original 13-page publisher PDF was obtained through the user's browser and read in full, with decisive equations checked visually. Dellaportas–Stephens's original 11 pages were read through the Imperial-authenticated JSTOR viewer. Gull's eight scanned pages and Werman–Keren's seven pages were read from primary hosted texts. Little's relevant example and AGK's coordinate-weighting sections were inspected, including the original 2013 preprint and corresponding 2018 revision. Source notes distinguish full-article inspection from selected-section inspection. Temporary signed URLs and authentication details are not included in the record.

Managed `canvalidk/VD-docs` main was resolved via the installed GitHub connector to the unchanged coherent snapshot `adfe14bbaafb29856f0a30da088495b7c9d483d2`. The ongoing task's already-read managed sources remain:

- `catalog.yaml`, relevant mass-estimator inventory;
- `Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md`;
- `Newton-analysis/_dscn_mass_estimator/mass_estimator_originality_research_2026-09-05.md`.

Active local sources inspected for this final synthesis include the three comparison reviews, the current `03_trace_and_evaluator/mass_estimator/estimator.py`, the September 7 law-assumed-correct steering and practical/uncertainty context, the September 6 composition characterization, and the three new focused evidence notes linked above. The current user steering governs the no-rejection interpretation over the older managed specification. The managed snapshot is not asserted to contain these newer local developments. This is a scoped review, not an exhaustive audit of managed history, inbox/unclassified material or all local records.

The user's account of the discovery is preserved in the separate contribution draft: the human supplied the existence conviction, push, requirements and part of the algorithmic approach; AI found the equation. Detailed algorithmic attribution requires the development record. No independent human or external specialist verification has been represented as completed.

All new research records are local under `pass to VD-docs/`, awaiting transfer. No managed VD-docs material, VDfirst content, production entry, or estimator code was changed in this final round. Existing unrelated work was left alone. No commit, push, publication or external message was made. The literature search stops with this report and draft.
