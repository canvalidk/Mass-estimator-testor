# Roe–Woodroofe: scalar boundary uncertainty audit — 2026-09-08

Status: local active/unsubmitted research note. Saving here does not submit it to VD-docs. This bounded lane consumed the primary paper below and the current task contract; no managed VD-docs source was consumed here. The coordinating review records its managed snapshot separately.

Verdict: **full scalar-boundary uncertainty coverage established** in the Bayesian sense specified below. **Demonstrated limitation as a native answer to the complete two-vector contract**: the inspected Gaussian model does not define the two-vector likelihood, unknown true direction, or mass estimator. A joint-vector adaptation remains **unresolved**. Different prior or coverage semantics are not themselves grounds for failure.

## Source facts and equations

[Roe and Woodroofe, *Setting Confidence Belts*, arXiv:hep-ex/0007048v2, 13 October 2000](https://arxiv.org/pdf/hep-ex/0007048). One-based PDF pages are used below. The original text was accessible; the Gaussian comparison uses section IV.A, not the earlier method discussed in section III.

Equations (1)–(2), pp. 3–4, specify $x=\theta+\Delta$, $\Delta\sim N(0,1)$, $\theta\ge0$, a flat half-line prior, and

$$
p(\theta\mid x)=\frac{\phi(x-\theta)}{\Phi(x)},\qquad\theta\ge0.
$$

Equations (3)–(6), p. 4, construct shortest posterior credible intervals by density threshold. Section V.A, p. 6, distinguishes posterior credibility from unconditional repeated-sampling coverage. Equation (12) bounds the latter by $(1-\epsilon)/(1+\epsilon)$; the reported minimum for nominal 90% intervals is about 86%. An optional upper-limit modification is treated separately. Appendix equation (A2), p. 10, establishes a specified conditional frequentist property using an independent replicate and an error-truncation condition. It is not ordinary unconditional nominal coverage.

## Deductions and explicit reconstruction

The posterior is proper for every finite real $x$, including zero and negative readings, because its normalization $\Phi(x)$ is strictly positive. Its density on $(0,\infty)$ is positive, its tails are Gaussian, and its interior quantiles and mean are therefore finite and positive. The mean, derived by integration, is the Katz expression

$$
E[\theta\mid x]=x+\frac{\phi(x)}{\Phi(x)}>0.
$$

The posterior mode $\max(x,0)$ can be zero. This is not a defect in the supplied posterior: choosing the mean or an interior quantile gives a strictly positive point summary. The paper's Gaussian contribution being assessed is its posterior and interval prescription; this note does not attribute a required mean-based point estimate to its authors.

For credibility $1-\epsilon$, $0<\epsilon<1$, write

$$
x_0=\Phi^{-1}\!\left(\frac1{1+\epsilon}\right),\qquad
d(x)=
\begin{cases}
\Phi^{-1}[1-\epsilon\Phi(x)],&x\le x_0,\\
\Phi^{-1}[(1+(1-\epsilon)\Phi(x))/2],&x>x_0.
\end{cases}
$$

Then the shortest interval is

$$
[\ell(x),u(x)]=[\max(x-d(x),0),\ x+d(x)].
$$

Every endpoint is finite for finite $x$. In the lower branch, $d(x)>-x$, because $1-\epsilon\Phi(x)>1-\Phi(x)=\Phi(-x)$, so the upper endpoint remains strictly positive even for negative readings. The two branches meet at $d(x_0)=x_0$. The displayed case split is the analytic solution of **one density-threshold optimization**, not a special repair for failed zero inputs. An uncertainty interval with lower endpoint zero is allowed by the task contract; that contract requires a positive point estimate, not every interval endpoint to be strictly positive.

At $x=0$, the posterior is half-normal. Its mean is $\sqrt{2/\pi}$, and the shortest interval is

$$
[0,\ \Phi^{-1}(1-\epsilon/2)].
$$

Known error standard deviation $\sigma>0$ is handled by standardizing $x/\sigma$ and multiplying the resulting parameter summaries and endpoints by $\sigma$. This uses the same prescription. It does not introduce a new treatment at zero.

## Coverage interpretation

The posterior interval contains probability $1-\epsilon$ for each observed $x$ under the stated prior and likelihood. That is enough to describe Bayesian uncertainty for the current comparison. It is distinct from a guarantee that repeated intervals cover each fixed true $\theta$ at exactly the nominal rate.

The optional frequentist modification should be treated as a separate published interval variant if tested. It is unnecessary to make the base posterior defined, positive, or continuous at zero. The paper's reported modified 90% coverage rounded to three significant figures is not a proof of an exact 0.9 lower bound.

## Relation to the mass task

The new Gaussian evidence here concerns one observed scalar normal mean. It contains no second noisy denominator, shared unknown spatial direction, or supplied full joint covariance for two measured vectors. Directly calling a force/acceleration quotient a normal measurement of mass would require an additional distributional argument; such a quotient need not be normal.

A known-direction construction can combine two such posteriors to define positive magnitudes, their ratio distribution, and a ratio-of-means point estimate, as derived in the companion Katz evidence note. That is a separate construction. Full correlated vector inference additionally needs a joint likelihood, direction treatment, and prior measure. The paper's scalar success neither supplies nor disproves that generalization.

What this adds: explicit published boundary uncertainty, including exact posterior normalization and shortest intervals, plus a clear distinction between Bayesian credibility and frequentist coverage. Its underlying scalar posterior repeats the Katz mechanism; it is not an independent discovery of a complete two-vector solution. No numerical sweep is needed to prove the scalar existence and continuity statements, and no broad candidate search was undertaken.
