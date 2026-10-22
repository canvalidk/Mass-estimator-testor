# Mass point functional: quotient, size bias, and loss — 2026-09-08

Status: local active/unsubmitted research note. Saving here does not submit it to VD-docs. This bounded audit inspected **two source groups only**: Little (2012) and Arratia–Goldstein–Kochman (2013 preprint, checked against its 2018 revision). It consumed the task-supplied equation and posterior specification, not a managed VD-docs source; the coordinating review records its managed snapshot separately.

Conclusion: **earlier quotient-of-posterior-expectations usage and exact coordinate-weighting mathematics are established**. Neither audited source is established as publishing the complete VD positive-vector equation. The weighted-mean identity is exact, but the tempting weighted-squared-error justification has infinite expected loss. A finite-regret characterization and a separate, finite nonnegative loss are derived below.

## Two source groups and exact anchors

### 1. Little (2012): an explicit quotient of posterior expectations

[Roderick J. Little, *Calibrated Bayes, an Alternative Inferential Paradigm for Official Statistics*](https://www.scb.se/contentassets/ca21efb41fee47d293bbee5bf7be7fb3/calibrated-bayes-an-alternative-inferential-paradigm-for-official-statistics.pdf), Journal of Official Statistics 28(3), 309–334.

Example 3, equations (4)–(7), printed pp. 318–319 (PDF pp. 10–11), concerns a survey population's precision-weighted mean. Its unnumbered display on printed p. 319 replaces

$$
E\!\left[\frac{\sum_j b_j U_j}{\sum_j U_j}\middle|D\right]
\quad\text{by}\quad
\frac{\sum_j b_j E[U_j\mid D]}{\sum_j E[U_j\mid D]},
$$

where this note abbreviates the paper's fixed sampled-stratum means by $b_j$ and uncertain stratum totals by $U_j$. The accompanying text explicitly identifies the operation as a ratio of posterior expectations and assigns an $O(1/n)$ approximation error.

This verifies an earlier algebraic form and stated use. It does **not** make the VD point value an approximation to $E[m\mid D]$: that expectation diverges for the task's posterior. Nor does this survey calculation specify the VD latent-vector manifold or prior measure.

### 2. Arratia–Goldstein–Kochman: exact joint coordinate weighting

[Arratia, Goldstein and Kochman, *Size Bias for One and All*, arXiv v1](https://arxiv.org/pdf/1308.2729v1), submitted 13 August 2013. Section 2.3, p. 9, equations (14)–(15), defines biasing a joint law by one nonnegative coordinate:

$$
\frac{d\mu^{(i)}}{d\mu}=\frac{x_i}{E_\mu X_i},\qquad
E_{\mu^{(i)}}g(X)=\frac{E_\mu[X_i g(X)]}{E_\mu X_i}.
$$

The assumptions include $0<E_\mu X_i<\infty$. The displayed identity is for bounded measurable $g$; extension to nonnegative functions follows by monotone convergence, as used in the paper's unbounded-function discussion, section 2.2.2, p. 5. The same coordinate formulas appear as (19)–(20), p. 11, in the [2018 v3 revision](https://arxiv.org/pdf/1308.2729v3). Section 2.3 itself credits preceding work; 2013 is therefore a verified availability date, not a claimed first invention of this general mathematics.

This establishes exact mathematical ancestry for the weighting identity. The inspected sections do not instantiate $X_i$ as acceleration, choose the VD posterior, or propose the physical mass summary. Those substitutions below are this audit's derivation.

## Exact application to the task posterior

Write $P=P(\,\cdot\mid y,\Sigma)$ for the posterior with density proportional to

$$
\exp\!\left[-\tfrac12(y-(fu,\alpha u))^T
\Sigma^{-1}(y-(fu,\alpha u))\right]
\,df\,d\alpha\,d\Omega,
\quad f,\alpha>0,\quad u\in S^{d-1},
$$

where $y=(F,a)$ is finite and $\Sigma$ is positive definite. Let $m=f/\alpha$, and let

$$
A=E_P[\alpha],\qquad B=E_P[f],\qquad m_*=B/A.
$$

The following arguments depend on the specified measure, not just the phrase “Bayesian inference.” Since $\|(fu,\alpha u)\|^2=f^2+\alpha^2$, a Gaussian bound gives constants $K,c>0$ such that the normalized amplitude-direction density is at most $K e^{-c(f^2+\alpha^2)}$. Consequently $A,B$ are finite and strictly positive. The joint amplitude density is also continuous and positive at finite points of either axis when extended to the boundary.

Define the acceleration-weighted joint law

$$
dQ=\frac{\alpha}{A}\,dP.
$$

Then, with no independence assumption,

$$
E_Q[m]=\frac{E_P[\alpha m]}{E_P[\alpha]}
=\frac{E_P[f]}{E_P[\alpha]}=m_*.
$$

This is **biasing by acceleration**, not size-biasing mass by itself. The latter would require a finite $E_P[m]$, which is absent here. Directional information remains in $P$; the weighting identity neither creates nor removes the need to specify that vector likelihood.

The interpretation does not require replacing the task's original mass-uncertainty distribution $P\circ(f/\alpha)^{-1}$ by $Q$'s mass distribution. They are distinct distributions. If original-posterior quantiles are reported alongside $m_*$, the point summary should be identified as a weighted mean, not called the ordinary posterior mean of mass.

## Why weighted squared-error Bayes risk fails

The tempting loss is

$$
L_2(d;f,\alpha)=\alpha(d-m)^2
=\frac{(d\alpha-f)^2}{\alpha},\qquad d>0.
$$

For every finite action $d$, $E_P L_2(d)=\infty$. To see this without subtracting infinities, restrict $f$ to a positive compact interval and $\alpha$ to a sufficiently small neighborhood of zero. The posterior density has a positive lower bound there, and $(d\alpha-f)^2$ is bounded below by a positive constant. Integration therefore includes $\int_0^\epsilon d\alpha/\alpha=\infty$. The same region also proves $E_P[f/\alpha]=\infty$.

Equivalently, $E_Q[m^2]=E_P[f^2/\alpha]/A=\infty$. Although $E_Q[m]$ exists, every ordinary squared-error risk under $Q$ is infinite. Completing a formal square and discarding the divergent constant is not a proof that $m_*$ uniquely minimizes finite Bayes risk for $L_2$.

### A valid finite-regret statement

For any fixed finite reference action $d_0$, the **pointwise difference** is integrable:

$$
E_P[L_2(d)-L_2(d_0)]
=(d^2-d_0^2)A-2(d-d_0)B.
$$

In particular,

$$
E_P[L_2(d)-L_2(m_*)]=A(d-m_*)^2\ge0.
$$

This uniquely selects $m_*$ by finite expected regret. It is the expectation of an integrable difference, not a subtraction of two infinite expectations. The criterion and its qualification should be stated explicitly.

## A finite, nonnegative loss with the same exact optimum

The following is a separate **loss interpretation derived in this audit**, not a source claim about Little or Arratia–Goldstein–Kochman:

$$
L_{\mathrm{KL}}(d;f,\alpha)
=f\log\!\left(\frac{f}{d\alpha}\right)-f+d\alpha,
\qquad d,f,\alpha>0.
$$

The logarithm's argument is dimensionless. This is generalized KL divergence between the positive force magnitude $f$ and the force magnitude predicted by $d\alpha$. Nonnegativity follows by writing it as $d\alpha\,[t\log t-t+1]$, with $t=f/(d\alpha)$; it is zero exactly when $d=m$.

For every finite action $d>0$, its posterior expectation is finite. The only apparent new boundary terms are logarithmic: $f|\log f|$ and $f|\log\alpha|$. The Gaussian bound above and integrability of $|\log x|$ near zero establish absolute integrability. No $1/\alpha$ singularity remains.

For any reference action $d_0>0$, its finite posterior risk has difference

$$
R_{\mathrm{KL}}(d)-R_{\mathrm{KL}}(d_0)
=A(d-d_0)-B\log(d/d_0).
$$

Thus $R'_{\mathrm{KL}}(d)=A-B/d$ and $R''_{\mathrm{KL}}(d)=B/d^2>0$. Its unique global minimizer is

$$
\arg\min_{d>0}E_P L_{\mathrm{KL}}(d;f,\alpha)=B/A=m_*.
$$

The exact regret is $B[t-1-\log t]$, where now $t=d/m_*$, again nonnegative with equality only at the optimum. The loss is also $\alpha$ times the Bregman divergence generated by $s\mapsto s\log s-s$, evaluated at $(m,d)$, after consistent choice of mass units. No ancestry or originality claim for that general divergence is made here; this lane did not audit a third source group on losses.

For comparison, ordinary squared force residual $(f-d\alpha)^2$ is finite but generally selects $E_P[f\alpha]/E_P[\alpha^2]$, a different functional. The choice of loss therefore matters; the finite KL interpretation is an available justification, not proof that physical considerations uniquely require it.

## What is established, and where to stop

- **Established predecessor of an algebraic operation:** Little explicitly uses a ratio of posterior expectations, as an approximation in his survey setting.
- **Established predecessor of an exact general identity:** coordinate size bias gives $E_Q[m]=E_P[f]/E_P[\alpha]$ after the stated substitution.
- **Established for the supplied VD posterior by this audit:** infinite weighted-squared-error risk; valid finite regret; and a finite nonnegative loss with the exact point estimate as its unique optimum.
- **Unresolved historical question:** neither audited source has been shown to publish the combined latent-collinearity model, $df\,d\alpha\,d\Omega$ measure, full-covariance treatment, physical point functional, and ratio uncertainty as the same complete method.

This note does not change the user's stated discovery account: the human supplied the existence conviction, direction of the search, and part of the algorithm; the AI found the equation. Bibliographic ancestry of mathematical components is a separate question from those contributions.

Coverage: the cited Little example and surrounding model definitions; AGK sections 2.1–2.3 in the original preprint and corresponding revision sections. No third source group, broad loss-function search, or new numerical sweep was undertaken. Stop after this focused note.
