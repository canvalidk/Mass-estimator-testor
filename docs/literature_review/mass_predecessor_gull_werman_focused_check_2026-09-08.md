# Final focused predecessor check: Gull and Werman–Keren

Date: 2026-09-08. Status: local active evidence, awaiting transfer to VD-docs; NOT submitted. Parent task supplies the managed-context snapshot. This note audits two primary publications only and does not establish historical priority.

## Target and comparison boundary

The target uses measurements $y=(F,a)$, true vectors $(fu,\alpha u)$, $f,\alpha>0$, and $u\in S^{d-1}$, with Gaussian measurement likelihood and latent measure $df\,d\alpha\,d\Omega$. Its point is $E[f\mid y]/E[\alpha\mid y]$; its mass uncertainty is the posterior pushforward of $m=f/\alpha$.

Our change of variables $v=\alpha u,\ f=m\alpha$ gives

$$
 df\,d\alpha\,d\Omega
 =\|v\|^{2-d}\,dm\,dv,
 \qquad
 \widehat m=\frac{E[m\|v\|\mid y]}{E[\|v\|\mid y]}.
$$

Thus the measure is flat in $(m,v)$ for $d=2$, but proportional to $1/\|v\|$ for $d=3$. An earlier general Bayesian framework, or an available choice of priors inside it, is different evidence from an earlier explicit publication of this measure, statistic, or numerical procedure.

## 1. Gull (1989)

Stephen F. Gull, *Bayesian Data Analysis: Straight-line fitting*, pp. 511–518, DOI [10.1007/978-94-015-7860-8_55](https://doi.org/10.1007/978-94-015-7860-8_55). [Author-hosted full text](https://bayes.wustl.edu/sfg/line.pdf). All eight scanned pages were rendered and read visually.

Native construction: independent Gaussian errors in both coordinates and latent line $\hat y_i=a\hat x_i+b$ (p. 511). Pages 512–513 use a prior uniform in $(x_0,y_0,\log R_x,\log R_y)$, with Gaussian latent abscissae of mean $x_0$ and variance $R_x^2$. Gaussian nuisance integration gives a posterior over scale and slope (p. 514). The simplified point recommendation on p. 517 minimizes

$$
 \frac{aV_{xx}-2V_{xy}+a^{-1}V_{yy}}
 {a\sigma_x^2+a^{-1}\sigma_y^2},
$$

with centered data sums $V$. It is an approximation obtained by taking infinite range and dropping a determinant term. The native priors and this optimization statistic differ from the target; no target norm-weighted posterior mean or radial/direction quadrature was found. Analytic nuisance integration is already present.

**Verdict:** demonstrated difference from the audited native construction; exact target publication not established. Its Gaussian latent prior, hyperparameters, intercept, and prior bounds cannot silently be replaced by the target measure. This is an important methodological predecessor, not evidence that Bayesian noisy-line inference or analytic marginalization is new.

## 2. Werman–Keren (2001)

Michael Werman and Daniel Keren, *A Bayesian Method for Fitting Parametric and Nonparametric Models to Noisy Data*, IEEE TPAMI 23(5), 528–534, DOI [10.1109/34.922710](https://doi.org/10.1109/34.922710). [Author-hosted full text](https://www.cs.haifa.ac.il/~dkeren/mypapers/werman-keren-fitting.pdf). All seven pages were inspected.

Sections 1.2–2.1 integrate latent points on a model before selecting its MAP instance:

$$
 p(D\mid M)=\prod_i\int_M p(p_i\mid p)\,p(p\mid M)\,dp.
$$

For lines, independent isotropic Gaussian noise and uniform line/point priors yield perpendicular least squares. The line-prior reference measure is not explicitly fixed. The integral uses arc length. Sections 2–4 use analytic integrals or numerical integration followed by optimization; the target weighted mean and radial/direction procedure were not found.

**Conditional identity — our derivation, not an explicit published claim:** map the vector components to $d$ planar data points $(a_i,F_i)$. In fixed nondimensional units, restrict to origin-fixed positive slopes and choose uniform angle $\theta=\arctan m$. Then

$$
 d\theta\prod_{i=1}^{d}ds_i
 =\frac{dm}{1+m^2}\prod_i\sqrt{1+m^2}\,dv_i
 =(1+m^2)^{d/2-1}\,dm\,dv.
$$

At $d=2$ this equals the target measure; at $d=3$ it differs. The origin restriction and angle prior are choices added here. Uniformity over affine lines does not uniquely specify this conditioning.

**Verdict:** conditional two-dimensional measure identity established; literal publication of that specialization unresolved. The full target estimator is not identified: the published point rule is MAP.

## Independent check of why the point-rule distinction matters

In those units, take $d=2$, $F=a=0$, and equal independent isotropic errors $\Sigma=\sigma^2I_4$. The target density in $(f,\alpha,u)$ is proportional to $\exp[-(f^2+\alpha^2)/(2\sigma^2)]$; symmetry gives $\widehat m=1$. Integrating the equivalent flat $(m,v)$ measure instead displays its mass density:

$$
 p(m\mid y)=\frac{2}{\pi(1+m^2)},\qquad m>0.
$$

The density relative to $d\theta$ is uniform. Thus a MAP line-angle prescription ties at every positive slope, whereas the target norm-weighted posterior statistic selects 1. A prior-measure match alone does not establish an identical point estimator. This calculation is an audit derivation, not a result attributed to either publication.

## Consequence for the contribution draft

These sources reinforce that integrating uncertain latent coordinates is established methodology. The two-dimensional conditional identity further narrows any claim that the latent measure itself is unprecedented. A supportable description is: the work specifies a particular magnitude-and-direction measure, the ratio of posterior magnitude means, its distinct ratio uncertainty, and an implementation for the stated covariance model. Claims of originality for that exact combination or implementation remain provisional; this bounded check did not locate an explicit earlier publication of them.

This lane stops here. No new literature expansion is proposed. The 1999 conference precursor listed on Keren's publication page was noticed bibliographically but not audited. Other references in these two papers were not followed. No claim is made that an unread source lacks the target construction. Recovered primary text sufficed; no author was contacted and no shared signed-in browser was operated.
