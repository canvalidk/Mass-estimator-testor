# Mardia positive-scale comparison: first evidence pass

Date: 2026-09-08.
Status: local active research record; awaiting transfer to VD-docs. Saving here does not submit it.
Scope: one candidate, equations and decisive boundary argument; no historical-priority finding. Managed-context snapshot and our equation are handled by the parent comparison record; this note independently audits the external source only.

## Primary source checked

[Mardia, Fallaize, Barber, Jackson and Theobald, *Bayesian alignment of similarity shapes*](https://arxiv.org/pdf/1312.1840), Annals of Applied Statistics 7(2), 989-1009 (2013), DOI 10.1214/12-AOAS615. The PDF identifies a 2013 journal publication and a 6 December 2013 arXiv deposit.

Relevant source locations: model equation (1), PDF page 4; posterior equation (2), page 6; scale conditional equations (3)-(5), page 7; proposal calculation, page 8; posterior medians/intervals, page 10; vector representation without translation, pages 13-14; isotropic-error limitation, page 18. Primary PDF text was inspected; web PDF screenshots failed, so equations were cross-checked between model, expansion and stated parameter definitions.

The published model estimates positive scale, rotation, translation and matching under isotropic Gaussian errors in a common latent space. It uses a gamma scale prior and reports posterior uncertainty. The paper also treats origin-based vectors without translation. Its scale density therefore deserves substantive comparison; a different application or prior does not disqualify it.

## Our specialization and proof

This section derives consequences; it does not attribute an all-cases theorem to the authors.

Set one known match, dimension d >= 1, rotation A = I, translation zero, and fix their noise parameter s = sigma_c > 0. Write measured x = F and y = a. For a gamma prior with shape k > 0 and rate b > 0, their equation (3) specializes to

$$
h(c;x,y)=c^{r-1}e^{-bc}\exp[-\|x-cy\|^2/(4s^2)],\quad
c>0,\qquad r=d/2+k.
$$

The normalized density is h divided by its integral. Expanding the square yields their equation (4) with

$$
\nu=\|y\|^2/(2s^2),\qquad
\delta=x\cdot y/(2s^2)-b,
\qquad p(c\mid x,y)\propto c^{r-1}e^{-\nu c^2/2+\delta c}.
$$

Since 0 < h <= c^(r-1) exp(-bc), its normalizer is positive and finite for every finite x,y. The same bound after multiplication by any nonnegative power of c proves all nonnegative moments finite. The mean is strictly positive; every interior quantile is positive and finite. Dominated convergence proves continuity in finite x,y and positive s, including paths through zeros. The strictly positive density makes interior quantiles continuous as well. This is one density prescription throughout.

At y=0, including x=y=0, the normalized density is exactly Gamma(r,b), independent of x. This also identifies the source of the boundary scale: the gamma prior and c^(d/2) factor. At x=0 but y nonzero, it remains proper. At fixed nonzero norms, changing x dot y changes the density; aligned, perpendicular and anti-aligned data are not collapsed together. The mean is strictly increasing in delta because its derivative is Var(c)>0.

## Demonstrated limits of the published model against the proposed contract

1. **Known covariance is the substantive mismatch.** From equation (1), in the original measured coordinates the independent channel covariances are c sigma^2 I and sigma^2 I/c. Equivalently they are s^2 I and s^2 I/c^2. These are isotropic and tied to the unknown scale. They are not an arbitrary supplied positive-definite joint covariance with unequal directional variances and cross-channel correlation. Simply placing s into a measurement-error input slot does not repair this difference.
2. **Free rotation changes the problem.** For a single pair in d >= 2 with a uniform rotation prior and zero translation, integrating exp[c x dot (Ay)/(2s^2)] over rotations depends only on the two norms. Hence the native freely aligned model cannot retain this pair's measured relative angle. Fixing A=I restores sign sensitivity, as demonstrated above; it is a legitimate specialization worth testing, rather than grounds for rejecting the entire Bayesian approach.
3. **Numerical formula versus density.** The displayed proposal mode on PDF page 8 divides by nu. At y=0, nu=0 and its written expression is 0/0. That does not invalidate the density or existence proof. An implementation would need a stable equivalent calculation or a sampler accepting the gamma boundary. No claim is made here that the authors' supplied implementation already does so.

The operative parameterization in equation (2) specifies the noise prior on sigma_c, independently of c; section 2.2 uses a gamma prior on sigma_c^(-2). Since sigma_c = sqrt(c) sigma, holding sigma_c fixed in equations (3)-(4) is different from holding the latent-space sigma in equation (1) fixed. In (c,sigma) coordinates the same independent prior p(c)p(sigma_c) becomes p(c)p_sigma_c(sqrt(c) sigma)sqrt(c). This derived Jacobian is another reason not to substitute a fixed measured noise scale into the conditional formula casually.

## Verdict for this checkpoint

**Full coverage established analytically for the fixed-alignment conditional scale density under its own noise model and a proper gamma prior. Demonstrated limitation against the proposed arbitrary supplied joint-covariance contract.** A covariance-correct adaptation remains a live competitor; its existence and historical publication status are separate questions. No equivalence to our particular point estimate or distribution is established here.

Next bounded decision: retain this specialization as a positive comparator in the first numerical matrix, clearly label its noise assumptions, and decide whether the review should audit a published covariance generalization before expanding the candidate list.
