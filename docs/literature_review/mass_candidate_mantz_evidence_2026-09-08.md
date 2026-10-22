# Mantz (2016): fresh full-covariance multivariate comparator

Date: 2026-09-08.
Status: local active research record; unsubmitted and awaiting transfer to VD-docs. Saving here does not submit it.
Scope: primary-source model audit plus explicitly labeled mathematical specializations. No skipped candidate was reopened. Managed-context provenance belongs to the parent comparison record.

## Primary source and anchors

[Adam B. Mantz, *A Gibbs Sampler for Multivariate Linear Regression*](https://arxiv.org/pdf/1509.00908), Monthly Notices of the Royal Astronomical Society 457, 1279-1288 (2016), DOI 10.1093/mnras/stv3008. [Author-deposit metadata](https://arxiv.org/abs/1509.00908) dates the first preprint to 3 September 2015 and the reviewed v2 to 2 February 2016. [Publisher record](https://academic.oup.com/mnras/article/457/2/1279/966253) confirms journal publication. The full 11-page author PDF and publisher article text were accessible; no paywall action is required.

Equation anchors, using author-PDF pages: (1)-(2), page 2, provide a full multivariate measurement covariance and Gaussian-mixture latent-covariate distribution; (5), page 4, gives the vector linear relation; (19)-(21), page 5, give the coefficient update. That page states that uniform coefficient priors are the default, with normal priors available. Equation (22) and footnote 14 address intrinsic covariance. The source treats independent objects, but allows correlations between every measured property of a single object.

## Why this does carry the full input covariance

Mantz's measurement equation is

$$
\begin{pmatrix}x_i\\y_i\end{pmatrix}
\sim N_{p+r}\!\left(
\begin{pmatrix}\xi_i\\\eta_i\end{pmatrix},M_i\right),
\qquad
\eta_i\sim N_r(b+B\xi_i,S_{\rm int}).
$$

Use one object, p=r=d, x=a, y=F and M=C, where C is the supplied 2d-by-2d positive-definite measurement covariance in that order. Impose the physical relation directly by setting b=0, B=m I_d and S_int=0. Then

$$
z=\begin{pmatrix}a\\F\end{pmatrix},\qquad
H_m=\begin{pmatrix}I_d\\mI_d\end{pmatrix},\qquad
z\mid m,v\sim N_{2d}(H_m v,C).
$$

This is exactly the required observation model with v=a_* and F_*=m v. All within-vector and cross-vector error correlations are already representable by the published M_i. It would be incorrect to reject this candidate because its objects are astronomical sources or because its notation is multivariate regression rather than two physical vectors.

The scalar-matrix restriction, known zero intercept, deterministic relation and m>0 must be specified as model constraints. The zero-scatter relation is a legitimate deterministic specialization; the observed Gaussian remains nonsingular because C is positive definite. The published Gibbs updates use inverse intrinsic covariance and free regression coefficients, so the unchanged sampler cannot simply be given a zero intrinsic covariance and assumed to implement this specialization. No claim is made that the published package already supports the constrained sampler.

## Our complete alternative with all-data existence proof

For an explicit alternative choose a fixed proper latent distribution v~N_d(mu,T), T positive definite, and a proper positive scale prior pi(m), for example a lognormal with fixed finite parameters and positive log standard deviation. These choices are made once for the model, not repaired for individual readings. Integrating v gives the same Gaussian latent-variable construction as the paper:

$$
p(m\mid z)\propto
\mathbf1_{m>0}\pi(m)\,
\phi_{2d}\bigl(z;H_m\mu,C+H_m T H_m^T\bigr).
$$

Equivalently let L(m)=integral phi_2d(z;H_m v,C) q(v) dv. For every finite z and every positive-definite C,

$$
0<L(m)\le K_C=(2\pi)^{-d}|C|^{-1/2}.
$$

The normalizer Z=integral_0^infinity pi(m)L(m)dm therefore obeys 0<Z<=K_C. This proves propriety everywhere. The posterior density is positive on all m>0, so every interior quantile, including the median, is finite and strictly positive. For lognormal pi, the mean is finite and positive because its numerator is at most K_C times the finite prior mean. Dominated convergence gives continuity over finite data and the open positive-definite covariance domain, including all measured-zero and angle-transition paths. The same density supplies uncertainty throughout.

Directional information is retained. For the illustrative isotropic case mu=0, T=tau^2 I, C=diag(s_a^2 I,s_F^2 I), define D(m)=s_F^2(tau^2+s_a^2)+m^2 tau^2 s_a^2. The log likelihood contains the term

$$
\frac{m\tau^2}{D(m)}\,a\cdot F.
$$

Its nonzero coefficient for m>0 distinguishes aligned, perpendicular and anti-aligned vectors at fixed norms. Full C uses the corresponding complete quadratic form, rather than discarding signs or replacing vectors by lengths. At a=F=0 the displayed posterior remains proper, with uncertainty determined by the declared priors, covariance and determinant factor. It does not claim that zero readings determine a unique mass without modeling assumptions.

This establishes full domain coverage for the explicitly adapted model and posterior summaries. It is not a historical claim that the paper published this particular positive-scale prior, scalar-matrix restriction or summary.

## Exact location of the target inside the same likelihood family

This transformation was derived in the parent comparison and independently checked here. The target uses f,alpha>0, common unit direction u, latent pair (f u,alpha u), and measure df d(alpha) dOmega. Put m=f/alpha and v=alpha u. Then

$$
df\,d\alpha\,d\Omega
=\alpha\,dm\,d\alpha\,d\Omega
=\|v\|^{2-d}\,dm\,d^dv.
$$

Therefore the target itself uses precisely the same full-C errors-in-variables likelihood above, with the improper prior measure

$$
\mathbf1_{m>0}\|v\|^{2-d}\,dm\,d^dv.
$$

For d=1 this is |v|dm dv; for d=2 it is flat dm d^2v; for d=3 it is dm d^3v/||v||. The target point is E[m||v||]/E[||v||], and its uncertainty is the posterior law of m. This is a model-class identity, not proof that Mantz published that prior or point. The preceding proper Gaussian/lognormal example is a different legitimate member of the same likelihood family.

The difference can also be seen in magnitude-direction coordinates: a prior pi(m)q(v)dm dv becomes pi(f/alpha)q(alpha u)alpha^(d-2)df d(alpha)dOmega, generally not the target's flat measure. Neither the differing measure nor a differing valid point summary excludes an alternative from the broad coverage contract.

## Native defaults and verdict

The native regression model estimates a full matrix of signed slopes, intercepts and intrinsic covariance from multiple objects. Its default flat coefficient prior gives no positivity guarantee. With one object and free intercept/slopes the native design matrix is underdetermined, so equation (20)'s inverse is unavailable; optional proper normal priors can address regression regularization but do not themselves impose a common positive scale. These are limitations of applying the default fit unchanged, not of the observation model.

**Full required likelihood coverage is already present in the published multivariate model. Full positive-estimate and quantile coverage is established for our explicit constrained, proper-prior adaptation.** Publication of the target's special measure and point, or of an earlier method already making exactly these positive common-scale choices, remains unresolved from this source.

This candidate materially changes the comparison: full arbitrary covariance and vector observations are not sufficient by themselves to separate the target from prior Bayesian errors-in-variables work. The more specific prior measure, point functional and any proposed advantages now deserve the scrutiny. The lane stops after Kelly and Mantz; no new implementation or outward literature branch was started.
