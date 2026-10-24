# Joint mass learning and what its uncertainty can promise

**Working contribution 4 — 15 September 2026.** Continues
[weighting and calibration](weighting_and_calibration.md). Newton II is assumed
throughout. These are derivations and research checks, not a selected replacement
for the laboratory estimator.

## Main findings

1. A joint model with one mass and several unknown acceleration vectors has a
   tractable Jeffreys reference measure. Multiplying already fitted single-trial
   mass distributions does not reproduce it. With a common instrument ratio and
   identical observed zero vectors, its mass distribution remains unchanged as
   trials are added.
2. Joint fitting can remove the earlier pooling inconsistency while leaving an
   uncertainty problem. In the repeated-trial example below, the point approaches
   the true mass, but its nominal 95% Bayesian interval approaches **87.87%
   sampling coverage**. Each additional trial brings another unknown vector.
3. The zero-signal problem is stronger. With **noisy measurements of truly zero
   excitation**, the joint posterior can concentrate around a random mass or
   drift toward zero or infinity as trials accumulate. Its central 95% interval
   coverage tends to zero.
   The unchanged distribution at exactly zero observations is not protection
   against this repeated-experiment failure.
4. We can construct mass-confidence sets with an exact 95% coverage guarantee
   under the declared Gaussian model. That guarantee leaves meaningful choices:
   the sets differ in width, shape, and which misleading candidate masses they
   exclude. Coverage does not by itself choose the best uncertainty output.

The distinction is consequential: a proper probability distribution, a good
point estimate, and well-calibrated repeated-experiment uncertainty are separate
achievements.

## 1. State which experiment is being repeated

This contribution studies independent trials with one fixed mass and **a
separate fixed unknown acceleration vector in each trial**:

$$
\mathbf F_i\sim N(m\mathbf v_i,\sigma_{F,i}^2 I_d),\qquad
\mathbf a_i\sim N(\mathbf v_i,\sigma_{a,i}^2 I_d),\qquad i=1,\ldots,N.
$$

The two measurement channels are independent, their positive variances are
known, and \(m>0\). The vectors \(\mathbf v_i\) are nuisance parameters: quantities
needed to describe the experiment, although mass is the quantity of interest.
Their directions and magnitudes may differ. No probability distribution over
the true excitation vectors is assumed in the sampling experiment.

This differs from repeatedly measuring **one unchanged** force/acceleration
pair. With that protocol, Gaussian sample means first combine the measurements
and reduce their variances; there is still one unknown vector. The distinction
must be made before choosing how to pool trials or interpret uncertainty.

A third possibility is a physical distribution of newly generated excitations,
such as \(\mathbf v_i\stackrel{iid}{\sim}N(0,\tau^2I_d)\). That supplies an additional
assumption about the apparatus or protocol. It is not implied by Gaussian
measurement errors. Section 7 returns to this option.

## 2. The joint reference measure and mass distribution

Define

$$
D_i(m)=\sigma_{F,i}^2+m^2\sigma_{a,i}^2,\qquad
c_i(m)=m^2/\sigma_{F,i}^2+1/\sigma_{a,i}^2.
$$

The nuisance blocks of the Fisher information matrix are \(c_iI_d\).
Eliminating them by a Schur complement leaves mass information
\(\sum_i\|\mathbf v_i\|^2/D_i\). Consequently the full joint Jeffreys volume is

$$
\boxed{
d\Pi_J\propto
\sqrt{\sum_i\frac{\|\mathbf v_i\|^2}{D_i(m)}}
\prod_i c_i(m)^{d/2}\;dm\prod_i d\mathbf v_i.
} \tag{1}
$$

For one trial this is the metric-tube measure derived in contribution 3.
For several trials the square root contains a **sum**. Taking a product of
single-trial reference densities gives a different measure. Independence of
the observations does not require the reference measure to factor after a
parameter has been shared.

This remains a reference choice, not a consequence of Gaussianity alone. It
also depends on the complete experimental design. Recomputing a joint Jeffreys
measure when new nuisance parameters are added need not equal Bayesian updating
under a previously fixed prior. A sequential use would need that prior and
nuisance extension declared in advance.

### Exact nuisance integration

Put

$$
Q_i=\frac{\|\mathbf F_i-m\mathbf a_i\|^2}{D_i},\quad
\boldsymbol\eta_i=
\frac{m\mathbf F_i/\sigma_{F,i}^2+\mathbf a_i/\sigma_{a,i}^2}{\sqrt{c_i}},\quad
w_i=\frac{\sigma_{F,i}\sigma_{a,i}}{D_i}.
$$

Completing the Gaussian squares and integrating all nuisance vectors yields
the unnormalized density relative to \(dm\):

$$
q_J(m)=e^{-\frac12\sum_iQ_i}
E\sqrt{\sum_iw_i^2\|\mathbf Z_i+\boldsymbol\eta_i\|^2},
\qquad \mathbf Z_i\stackrel{ind}{\sim}N(0,I_d). \tag{2}
$$

The nuisance Jacobians cancel the product in (1). The remaining norm
expectation is essential; discarding it changes the construction.

For every finite dataset and finite \(N\ge1\), the \(\boldsymbol\eta_i\) are
bounded as mass varies. Equation (2) therefore has positive finite limits

$$
q_J(m)\longrightarrow C_0\quad(m\downarrow0),\qquad
q_J(m)\sim C_\infty m^{-2}\quad(m\to\infty).
$$

The posterior is proper. Exactly the power moments with \(-1<p<1\) are finite;
the ordinary mass mean, inverse-mass mean, and second mass moment still diverge.
All fixed-order absolute log-mass moments are finite. More trials can sharpen
central quantiles without removing those asymptotic tails.

### Common instrument ratio gives a one-dimensional calculation

Suppose \(s=\sigma_{F,i}/\sigma_{a,i}\) is common to all trials. Stack the
standardized observations \(\mathbf F_i/\sigma_{F,i}\) and
\(\mathbf a_i/\sigma_{a,i}\) into \(X,Y\in\mathbb R^K\), with \(K=Nd\).
For \(r=m/s\), define

$$
R_r=\frac{X-rY}{\sqrt{1+r^2}},\qquad
S_r=\frac{rX+Y}{\sqrt{1+r^2}},\qquad
h_K(t)=E\|Z+t e_1\|,\quad Z\sim N(0,I_K).
$$

Then (2) reduces to

$$
\boxed{p_J(r\mid X,Y)\propto
\frac{e^{-\|R_r\|^2/2}}{1+r^2}\,h_K(\|S_r\|).} \tag{3}
$$

In \(\theta=\arctan r\), the Jacobian cancels \((1+r^2)^{-1}\). The script
evaluates the remaining noncentral Gaussian norm using

$$
h_K(t)=\sqrt{\frac2\pi}\int_0^{\pi/2}
\frac{1-\cos^K\beta\,e^{-t^2\sin^2\beta/2}}{\sin^2\beta}\,d\beta.
$$

The apparent singularity at zero is removable. The fraction tends to
\((K+t^2)/2\), which also guides stable numerical evaluation.

### Zero observations expose an important pooling difference

If **all observed vectors are exactly zero**, (3) gives

$$
p_J(r)=\frac{2}{\pi(1+r^2)},\qquad
p_J(\ell)=\frac1{\pi\cosh\ell},\quad \ell=\log r,
$$

for every \(N,d\). The width does not shrink. In comparison, multiplying
single-trial posterior densities in \(\ell\) gives a density proportional to
\(\operatorname{sech}^N\ell\), whose variance tends to \(1/N\). Multiplication
in mass coordinates instead gives \((1+r^2)^{-N}\), concentrating toward zero.
These are distinct pooling/reference prescriptions, not two equivalent ways
to combine the independent likelihoods.

Exactly zero observations and true zero excitation are different events. At
true \(\mathbf v_i=0\), the random observations are noisy and their sampling law
is identical for every mass. No inference construction creates identification
there. The formula above is a diagnostic for a specified observed dataset,
not a claim that every random null dataset produces the same posterior.

## 3. A declared extension of the ratio-of-means readout

For homogeneous channel standard deviations across trials, concatenate the
physical latent acceleration vectors into \(V\). One possible extension is

$$
\widehat m_{\rm stack}=
\frac{E\|F^*_{\rm stack}\|}{E\|a^*_{\rm stack}\|}
=\frac{E[m\|V\|]}{E\|V\|}.
$$

It reduces to the original readout for one trial. Choosing a concatenated norm
is an additional aggregation decision: this is not physical addition of
vectors from different trials, and the original single-pair uniqueness proof
does not by itself select this multi-trial operation.

Under (1), the completed-square Gaussian nuisance law is tilted by a norm.
The additional norm in this readout therefore produces a second moment:
\(E\|Z+S\|^2=K+\|S\|^2\). Hence

$$
\boxed{\frac{\widehat m_{\rm stack}}s=
\frac{\int_0^{\pi/2}e^{-\|R_\theta\|^2/2}
        \sin\theta\,(K+\|S_\theta\|^2)\,d\theta}
     {\int_0^{\pi/2}e^{-\|R_\theta\|^2/2}
        \cos\theta\,(K+\|S_\theta\|^2)\,d\theta}.} \tag{4}
$$

At zero observations this equals \(s\). Both expectations defining it are
finite, despite the divergent unweighted mass mean. If only the instrument
ratios agree while the absolute uncertainties differ, (4) describes a
**standardized-stack** readout. A physical-stack norm needs a weighted norm
expectation instead. The numerical experiment below uses homogeneous scales.

## 4. Three exact mass-confidence constructions

A 95% confidence guarantee means

$$
P_{m,V}\{m\in\mathcal C(\text{observations})\}=0.95
$$

at every fixed positive true mass and every fixed collection of excitation
vectors. It is a repeated-experiment statement, not 95% posterior probability
for the particular observed set. Each construction below has this guarantee
under Section 1's assumptions **with the common instrument ratio used in
Section 2**. The full-residual construction extends more broadly as described
below.

At the true \(r\), \(R_r\sim N(0,I_K)\), independently of
\(S_r\sim N(\sqrt{1+r^2}\,v,I_K)\), where \(v\) is the standardized latent stack.
Write \(\chi^2_{k,p}\) for the \(p\) quantile of a chi-square law with \(k\)
degrees of freedom.

| Construction | Accept a candidate \(r>0\) when | Why coverage is exact |
|---|---|---|
| Full residual | \(T=\|R\|^2\le\chi^2_{K,.95}\) | \(T\sim\chi^2_K\) at the true candidate |
| Directional score | \(Z^2=(R\cdot S)^2/\|S\|^2\le\chi^2_{1,.95}\) | Given nonzero \(S\), \(Z\sim N(0,1)\) |
| Split residual | \(Z^2\le\chi^2_{1,p_Z}\) and \(W=T-Z^2\le\chi^2_{K-1,p_W}\) | The two statistics are independent chi-squares and \(p_Zp_W=.95\) |

The script calls these `omnibus`, `score`, and `split`. For the split it uses
\(p_Z=p_W=\sqrt{.95}\) when \(K>1\). When \(K=1\), the orthogonal component
vanishes identically; the score threshold must use .95 instead.

For each fixed candidate, \(S=0\) has probability zero. The prototype explicitly
includes that candidate in the score and split sets. This convention preserves
pointwise coverage, but can create an isolated accepted mass for exactly
antiparallel data. A continuous-extension convention is another choice and
would give a different set on such data. It should not be hidden in polynomial
algebra or used to select a point estimate.

### What the guarantee leaves undecided

With \(A=\|X\|^2,B=\|Y\|^2,C=X\cdot Y\), full-residual inversion is the quadratic

$$
(B-q)r^2-2Cr+(A-q)\le0,\qquad q=\chi^2_{K,.95}.
$$

Score inversion is the quartic

$$
\{C(1-r^2)+(A-B)r\}^2
\le q(1+r^2)(Ar^2+2Cr+B),\qquad q=\chi^2_{1,.95},
$$

away from the separately treated \(S=0\) case. The split adds the inequality
\(W=(AB-C^2)/\|S\|^2\le\chi^2_{K-1,p_W}\).

These sets can be unbounded, disconnected, or empty. At zero observations
all three are the whole positive domain. They are not always-finite error bars,
and do not evade the weak-identification limitation in contribution 3.

Two nondegenerate examples illustrate the remaining choice. Here \(K=3\)
and the table reports sets of dimensionless mass \(r\):

| Observed Gram summaries | Full residual | Score | Split |
|---|---|---|---|
| \(A=B=100,C=-90\): strongly opposed | Empty | \([0.93265,1.07222]\) | Empty |
| \(A=B=100,C=0\): orthogonal, equal norms | Empty | All positive masses | Empty |

The score is zero at stationary points of \(T\). For \(C<0\), the positive
stationary point is the **maximum** of \(T\), yet the score accepts it. The
orthogonal residual catches information the directional score has omitted.
This is a limitation in excluding misleading masses, not a coverage failure.

The split has a cost. Even with strong informative data it sometimes returns
an empty set; in the regular large-trial regime its empty-set probability
approaches \(1-\sqrt{.95}\), about 2.53%. When the orthogonal test passes, its
directional cutoff is wider than the score-only cutoff.

An empty confidence set is an outcome of the declared mass procedure. None of
these sets gates the ratio-of-means point or decides whether Newton II is
admitted. Such use would be a further policy decision.

The full-residual proof also extends to known anisotropic and correlated
channels: replace each scalar denominator by
\(V_i(m)=\Sigma_{FF,i}+m^2\Sigma_{aa,i}-m(\Sigma_{Fa,i}+\Sigma_{aF,i})\) and use
\(\sum_i(\mathbf F_i-m\mathbf a_i)^TV_i(m)^{-1}(\mathbf F_i-m\mathbf a_i)\).
The score and split formulas above use the narrower standardized model.
Estimating covariance from the data requires a separate calibration argument.

## 5. Repeated-experiment check

The reproducible experiment uses \(d=3\), homogeneous scales, and independent
Gaussian channel errors. Each trial has an unknown excitation vector of the
listed standardized magnitude; directions cycle through the three coordinate
axes. Those directions are not supplied to the estimator. Each row contains
8,192 independent simulated datasets, using seed 20260916.

The Bayesian interval is the central equal-tailed 95% interval from (3).
Its coverage is evaluated by checking whether the posterior CDF at the true
mass lies between .025 and .975. The other columns evaluate the exact
confidence inequalities at the true mass.

| Trials | True \(m/s\) | Acceleration signal/noise per trial | Full-residual coverage | Score coverage | Split coverage | Joint Bayesian interval coverage | Median stacked point / true mass |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 16 | 0.25 | 94.82% | 94.74% | 94.67% | 91.67% | 0.2623 |
| 1 | 2 | 1 | 94.92% | 95.04% | 95.10% | 96.28% | 0.7308 |
| 8 | 2 | 1 | 94.56% | 94.80% | 94.54% | 88.87% | 0.9928 |
| 64 | 2 | 1 | 95.35% | 95.00% | 95.24% | 87.62% | 1.0020 |
| 512 | 2 | 1 | 95.41% | 94.36% | 94.91% | 87.30% | 1.0000 |

At 95% coverage the Monte Carlo standard error is about 0.24 percentage points
per row; these are finite random estimates of the exact coverage results.
The posterior CDF changed by at most \(3.2\times10^{-14}\) when both integration
orders were doubled for the first 32 datasets in each row. This is a scoped
quadrature check, not an error bound over arbitrary data.

### Why the joint Bayesian interval remains too narrow

This result has an analytic explanation. In homogeneous physical units put

$$
V_0=\sigma_F^2+m_0^2\sigma_a^2,\qquad
E_N=\sum_i\|\mathbf v_i\|^2,\qquad K=Nd.
$$

For a regular sequence with positive limiting average excitation
\(E_N/N\to e>0\) and no dominating trial, the locally consistent profile
mass estimate has leading sampling variance

$$
\boxed{\operatorname{Var}(\widehat m)\simeq
\frac{V_0}{E_N}+\frac{K\sigma_F^2\sigma_a^2}{E_N^2}.} \tag{5}
$$

To see the extra term, work in standardized coordinates and let
\(L=\sum_i\|\mathbf v_i/\sigma_a\|^2\). At the true ratio,
\(R=\varepsilon_R\) and
\(S=\sqrt{1+r_0^2}\,v+\varepsilon_S\), with independent standard Gaussian errors.
For the profile estimating equation \(H=R\cdot S=0\),

$$
\operatorname{Var}(H)=(1+r_0^2)L+K,\qquad
E\!\left[\frac{dH}{dr}\right]=-L.
$$

The noise in the estimated direction contributes the \(K\) term. Linearizing
the root and restoring physical units gives (5).

Inverse profile curvature gives only \(V_0/E_N\). In (3), the remaining norm
factor has order-one local log derivatives, while the likelihood curvature
grows as \(N\). The joint Jeffreys posterior consequently inherits the smaller
local variance. Its stacked point has the same leading center in this regular
regime, so point consistency can coexist with interval undercoverage.

For \(m_0=2,\sigma_F=\sigma_a=1,d=3,\|\mathbf v_i\|=1\), the two variances are
\(8/N\) and \(5/N\). The limiting coverage is therefore

$$
2\Phi\!\left(1.95996398454\sqrt{5/8}\right)-1
=0.8787354453.
$$

This explains the simulation without attributing it to integration error.
It also identifies a specific calibration issue to address. Multiplying every
reported interval width by \(\sqrt{8/5}\) would only address this particular
asymptotic example; it is not a general weak-signal correction.

### At true zero excitation, more trials can create false precision

Now set every true \(\mathbf v_i=0\) and keep generating noisy observations.
Their joint sampling law contains **no mass information**: it is exactly the
same for every positive mass, for every number of trials. Nevertheless, the
joint Jeffreys posterior increasingly favors a mass, or a boundary direction,
selected by random fluctuations.

This is different from the artificial dataset in which every observation is
exactly zero. A separate simulation, with 4,096 datasets per row and a fixed
reference true ratio \(r_0=2\), gives:

| Trials | Full-residual coverage | Score coverage | Split coverage | Joint Bayesian interval coverage |
|---:|---:|---:|---:|---:|
| 1 | 94.21% | 95.29% | 94.78% | 96.95% |
| 8 | 94.41% | 95.07% | 94.68% | 75.32% |
| 64 | 95.19% | 94.95% | 95.09% | 39.65% |
| 512 | 95.00% | 95.21% | 95.04% | 19.90% |
| 4,096 | 95.61% | 95.02% | 95.07% | 10.52% |

At 95% coverage the Monte Carlo standard error here is 0.34 percentage points.
The largest CDF change on doubling both quadrature orders for the first 32
datasets per row was below \(10^{-13}\).

The limiting failure can also be derived. At true zero, \(X,Y\) are independent
standard Gaussian stacks. Define

$$
G=(\|X\|^2-\|Y\|^2)/2,\qquad C=X\cdot Y.
$$

As \(K\) increases, \((G,C)/\sqrt K\) tends to a pair of independent standard
normal variables. Up to a term constant in \(\theta\),

$$
\log p_J(\theta\mid X,Y)=
\tfrac12[-G\cos(2\theta)+C\sin(2\theta)]+o_P(1). \tag{6}
$$

The norm factor in (3) has vanishing angular log variation of order
\(K^{-1/2}\). The random angular variation of the profile term instead grows
as \(\sqrt K\). Writing \(\rho=\sqrt{G^2+C^2}\), an interior preferred angle has
local variance \(1/(2\rho)\), so its width shrinks as \(K^{-1/4}\). If a boundary
is preferred, the distance to it typically shrinks as \(K^{-1/2}\).

In distribution, the limiting preferred angle is uniform in \((0,\pi/2)\)
with probability one half, and each boundary receives probability one quarter.
The central 95% posterior interval consequently misses any fixed interior true
angle with probability tending to one. Equation (5), which assumes positive
limiting excitation, must not be extended to this null case.

This disqualifies the joint Jeffreys credible interval as a generally
sampling-calibrated uncertainty policy for the declared separate-vector
protocol. It does not refute Newton II or the ratio-of-means uniqueness theorem.
It exposes how much the chosen nuisance reference can influence mass uncertainty
when the sampling experiment cannot identify mass.

## 6. A first concrete answer about information left unused

The distinction between the confidence procedures is already measurable.
With a new unknown vector per trial and fixed positive average excitation:

- Full residual energy fluctuates on the order of \(\sqrt N\), and the cutoff
  exceeds its null mean \(K\) by that same order. A small mass displacement
  \(\delta\) increases residual energy on the order of \(N\delta^2\). When
  nonempty, its usual local width is of order \(N^{-1/4}\).
- The directional score changes linearly with mass and gives local width of
  order \(N^{-1/2}\), while keeping exact candidate-wise coverage. Its leading
  width includes the extra noise term in (5).

Thus a correct Gaussian assumption and a 95% guarantee do not ensure that a
procedure extracts mass information effectively. The full-residual method
spends sensitivity on fluctuations that are relatively unhelpful for locating
mass. The score concentrates on mass displacement but can overlook large
orthogonal discrepancies, as the examples show.

This is a comparison of these specified uncertainty procedures. It is not
yet a percentage of information used by the original ratio-of-means equation,
nor evidence that its point dominates alternatives. The fixed-dimension,
strong-signal Fisher benchmark in contribution 3 concerns a different regime;
there is no contradiction with the additional term in (5).

## 7. What this suggests investigating next

**First settle the experimental meaning of a batch.** Repeated readings of
one stable pair, separate fixed unknown excitations, and excitations drawn
from a physical distribution justify different joint models. The theorem
about a single readout does not decide among them.

For an independently generated Gaussian excitation distribution with variance
\(\tau^2\), nuisance integration gives the per-trial covariance

$$
\begin{pmatrix}
(\sigma_F^2+m^2\tau^2)I_d & m\tau^2 I_d\\
m\tau^2 I_d & (\sigma_a^2+\tau^2)I_d
\end{pmatrix}.
$$

That is a finite-parameter alternative worth testing if the protocol supports
it. The excitation variance should be learned or independently warranted;
fixing a wrong value can shift the inferred mass even with many trials. At
\(\tau=0\), mass is again unidentified. Repeated readings of one common random
vector would also create cross-trial covariance and need a different model.

**Then compare uncertainty at a common guarantee.** Candidate exclusion power,
width, unbounded-set frequency, and empty-set frequency are practical comparison
targets. The true-zero and near-zero regimes must be included: a method that
sharpens around random masses there has failed an essential calibration check.
The subsequent [conditional uncertainty analysis](conditional_uncertainty_and_design.md)
derives and implements a conditional likelihood-ratio benchmark for this
workflow under known independent isotropic noise. It addresses the score's
stationary-maximum problem while retaining the stated calibration; extending
that guarantee to estimated calibration remains a separate question.

The general strategy of conditioning away nuisance dependence and inverting
tests is developed in [Moreira, *A Conditional Likelihood Ratio Test for
Structural Models* (2003)](https://drphilipshaw.com/Moreira-2003.pdf), especially
Sections 2.1–2.2. That paper also motivates checking power under weak
identification. Our displayed pivots and joint mass marginal are derived here;
the script is not an implementation of Moreira's conditional likelihood-ratio
test. This is a focused methodological reference, not a priority assessment.

## Reproduction, verification, and provenance

[joint_checks.py](joint_checks.py) contains the posterior quadrature, declared
stacked readout, exact confidence inequalities, and a prototype polynomial
inversion retaining unbounded and disconnected sets. From the workspace root:

```powershell
& 'C:/Users/canva/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B '03_trace_and_evaluator/mass_estimator/investigation/joint_checks.py' --replicates 8192 --output '.tools/mass_equation_20260915/joint/results_8192.json'
```

For the separate true-zero study:

```powershell
& 'C:/Users/canva/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B '03_trace_and_evaluator/mass_estimator/investigation/joint_checks.py' --null-study --replicates 4096 --output '.tools/mass_equation_20260915/joint/null_results_4096.json'
```

The completed run passed quantile and noncentral-norm checks, exact zero-data
identities, one-trial agreement with the earlier independent sphere
integration, and confidence-set membership checks against the defining
inequalities on generic sampled Gram matrices and explicit degenerate examples.
Independent review identified and corrected a rank-one roundoff error; direct
membership and set inversion now agree at the checked exact singular candidates.
Monte Carlo coverage uses the inequalities directly, not polynomial roots.
The root routine remains a research prototype: unresolved multiple roots and
extreme scales are outside its numerical validation. The coverage proofs apply
to the exact constructions, not to arbitrary floating-point inversion.

The numerical output records the script SHA-256 and source commit. Supporting
derivation probes and machine output are in
`.tools/mass_equation_20260915/joint/`. No existing point estimator, uncertainty
policy, or managed document was changed by this contribution.

Managed context was checked at the single VD-docs `main` commit
[`8782a293297330f7a354d87ebe62792018fa8866`](https://github.com/canvalidk/VD-docs/commit/8782a293297330f7a354d87ebe62792018fa8866).
The relevant inventory and canonical source reads were:

- [catalog.yaml](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/catalog.yaml), selected by mass-estimator paths and the source map from contributions 1–3;
- [full working specification](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/inertial_mass_estimator_full_working_spec_d1.md);
- [readout audit](https://github.com/canvalidk/VD-docs/blob/8782a293297330f7a354d87ebe62792018fa8866/Newton-analysis/_dscn_mass_estimator/what_forces_the_ratio_of_means_readout.md).

The active local scope was the three preceding investigation contributions,
their checks, and the new script and supporting derivations. These results are
new local working material, not ingested managed context. No exhaustive
literature or wider VD document search is claimed in this pass.
