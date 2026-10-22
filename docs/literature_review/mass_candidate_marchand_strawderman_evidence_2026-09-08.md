# Marchand–Strawderman: restricted-parameter and prior-measure audit — 2026-09-08

Status: local active/unsubmitted research note. Saving here does not submit it to VD-docs. This bounded lane used the current comparison contract and the publication identified below. It consumed no managed VD-docs source; the coordinating review records that snapshot separately. This is **one candidate/source group** reached through the Katz citation trail. No additional referenced candidate was independently audited or counted.

Verdict: **demonstrated limitation of direct application of the audited published constructions**. They do not provide the stipulated positive proportionality inference for two vectors. **A complete published predecessor remains unresolved**. The review does establish substantial preceding restricted-Bayes machinery. Its abstract generality should neither be counted as an already-published solution of this exact problem nor dismissed as incapable of producing one.

## Source and access

[Marchand and Strawderman (2004), *Estimation in Restricted Parameter Spaces: A Review*](https://doi.org/10.1214/lnms/1196285377), IMS Lecture Notes–Monograph Series 45, 21–44. The [author-uploaded manuscript](https://www.researchgate.net/publication/38389719_Estimation_in_restricted_parameter_spaces_A_review) was readable in full through its extracted text. It is numbered 1–24 and dated April 26, 2004. **Page anchors below refer to that manuscript**, not inferred final-publication pagination. Publisher/JSTOR retrieval produced no readable body, and the published-PDF mirror was rate-limited.

Theorem 1, manuscript pp. 3–4, treats $X\sim N_p(\theta,I_p)$, with $\theta$ in an ambient convex set $C$ having nonempty interior. Its uniform-prior estimator is

$$
\delta_U(x)=x+\nabla\log g_C(x),\qquad
g_C(x)=\int_C\phi_p(x-\nu)\,d\nu.
$$

It proves squared-error domination over $X$, using equations (1)–(2). Section 4.1, p. 8, discusses posterior means and boundary admissibility. Section 5, pp. 10–13, treats two independent noisy normal vectors with known isotropic variances and a difference restriction $\theta_1-\theta_2\in A$. Definition 1 and equation (4) estimate $\theta_1$. The introduction explicitly separates its point-estimation emphasis from interval estimation.

## Deductions for this task

These arguments are this audit's deductions; they are not results attributed to the review.

### 1. The compatible-pair set is not the theorem's convex region

For spatial dimension $d\ge2$, put

$$
\mathcal M=\{(f u,\alpha u):f>0,\alpha>0,\ u\in S^{d-1}\}
\subset\mathbb R^{2d}.
$$

This parameter set has dimension $d+1<2d$. It has empty ambient interior and zero $2d$-dimensional Lebesgue measure. Substitution of $C=\mathcal M$ into the displayed ambient integral therefore gives $g_C(x)=0$, not a normalized posterior.

It is also nonconvex. With distinct coordinate vectors $e_1,e_2$, both

$$
(F_1^*,a_1^*)=(e_1,e_1),\qquad
(F_2^*,a_2^*)=(2e_2,e_2)
$$

are permitted. Their midpoint has force $(e_1+2e_2)/2$ and acceleration $(e_1+e_2)/2$, which are not collinear. Thus Theorem 1's convexity hypothesis also fails. Adding the common zero pair does not change either conclusion.

Whitening a known positive-definite covariance is an invertible linear map; it preserves nonconvexity and the zero-measure obstruction. Consequently, a covariance change alone does not bring the task within this theorem. A separately specified measure on $\mathcal M$ can define a valid posterior, but that is an additional mathematical specification.

### 2. “Uniform on compatible pairs” leaves a real measure choice

For the embedding $h(f,\alpha,u)=(fu,\alpha u)$, its Euclidean metric is

$$
ds^2=df^2+d\alpha^2+(f^2+\alpha^2)\,ds_{S^{d-1}}^2.
$$

Therefore the induced surface-volume element is

$$
dH_{d+1}=(f^2+\alpha^2)^{(d-1)/2}\,df\,d\alpha\,d\Omega.
$$

For three dimensions the factor is $f^2+\alpha^2$. The VD measure $df\,d\alpha\,d\Omega$ supplied for this audit is therefore different from uniform induced surface measure, as well as from the review's ambient Lebesgue measure. This observation identifies a specification to compare; it does not establish that one measure is universally preferable or new.

### 3. The two-vector difference result does not encode positive proportionality

Section 5's condition $\theta_1-\theta_2\in A$ depends only on a fixed difference. Positive proportionality cannot be characterized that way. For example, the valid pair $(2e_1,e_1)$ and invalid noncollinear pair $(e_1+e_2,e_2)$ have the same difference $e_1$. Any fixed $A$ either accepts both or rejects both.

Thus this is a precise limitation of the published reduction, not a rejection simply because the paper discusses another application. Its estimate of one latent vector also requires an additional choice of scalar functional before it becomes a mass estimate. Posterior averages of compatible vector pairs need not themselves be compatible, as the midpoint example shows.

### 4. Broad Bayesian adaptability is real, but is not publication evidence

To make the distinction concrete, the following is **a newly stated illustrative construction, not an estimator attributed to this review**. Let $v=a^*$, $m>0$, and

$$
y=(F,a),\qquad h(m,v)=(mv,v),\qquad
L_\Sigma(y\mid m,v)=\phi_{2d,\Sigma}(y-h(m,v)).
$$

Choose once, as declared model assumptions, a proper density $\pi(m,v)>0$ on $(0,\infty)\times\mathbb R^d$, with finite prior mean of $m$. For instance, a lognormal mass prior and independent Gaussian $v$ prior with fixed, declared reference scales have these properties. These assumptions are different from the VD measure and are allowed by the present broad contract.

For every finite observation and every positive-definite $\Sigma$, the Gaussian likelihood satisfies $0<L_\Sigma\le B_\Sigma<\infty$. Hence

$$
0<Z(y,\Sigma)=\int L_\Sigma\,d\pi\le B_\Sigma,
\qquad
0<E[m\mid y,\Sigma]
\le\frac{B_\Sigma E_\pi[m]}{Z(y,\Sigma)}<\infty.
$$

The posterior is proper. Its mass marginal assigns positive probability to every positive open interval, so all interior mass quantiles are finite and strictly positive. Local dominated convergence over finite $y$ and positive-definite $\Sigma$ gives continuity of normalization and mean. Continuity of the posterior CDF together with its strict increase gives continuity of the interior quantiles. No special branch for zeros or angles is required.

This likelihood also retains directional signs. For example, after nondimensionalizing and taking $\Sigma=I$, $v\sim N_d(0,\tau^2 I)$, integration over $v$ gives

$$
p(F,a\mid m)\propto
k(m)^{-d/2}
\exp\!\left[-\frac{\|F\|^2+\|a\|^2}{2}
+\frac{\|mF+a\|^2}{2k(m)}\right],\qquad
k(m)=m^2+1+\tau^{-2}.
$$

The squared norm contains $2mF\cdot a$; reversing just one measured vector changes the mass likelihood. The full-covariance prescription remains the joint Gaussian integral above.

This construction supplies an existence argument for a broader class of solutions with different priors and point summaries. **It supplies no date or source establishing that an earlier author published this particular application.** It also makes no minimaxity claim. Its role is to prevent the leap from “these particular published formulas do not directly fit” to “only our equation could meet the contract.”

## Review implication and stopping point

What this adds beyond Katz: a direct check of a broad multivariate restricted-normal theorem, a precise distinction among ambient, induced-surface, and amplitude–direction measures, and an audit of a genuine two-noisy-vector result. None establishes this review as a full published predecessor for positive proportionality inference. The additional proper-prior construction shows why the broad contract and exact VD-measure/equation priority must remain separate questions.

Inspection scope: introduction, sections 2, 4.1 and 5 closely; remaining section headings, scalar developments, and relevant citation entries for scope. Search terms included the exact title, Katz, ratio, additional information, and admissibility. No citation beyond this source group was promoted to a separately reviewed candidate. Stop after this note; no broad expansion was undertaken.
