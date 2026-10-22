# Final focused predecessor lane: Lindley-El-Sayyad and Dellaportas-Stephens

Date: 2026-09-08.
Status: local active research record; unsubmitted and awaiting transfer to VD-docs. Saving here does not submit it.
Scope: two older source groups in the Bayesian errors-in-variables citation trail, aimed at exact prior/estimator identity. This is not another broad coverage survey. No skipped candidate was reopened. The parent comparison holds managed-source provenance and the contribution draft.

## Target identity being checked

The target observation law is Gaussian with supplied positive-definite joint covariance C and latent pair (f u,alpha u), f,alpha>0, common unit direction u. Its chosen measure is df d(alpha) dOmega, point estimate E[f]/E[alpha], and mass uncertainty law f/alpha. With m=f/alpha and v=alpha u,

$$
df\,d\alpha\,d\Omega
=\|v\|^{2-d}\,dm\,d^dv,
\qquad m>0,
\qquad \widehat m=\frac{E[m\|v\|]}{E[\|v\|]}.
$$

All expectations refer to the normalized target posterior. Matching only its Gaussian observation law is insufficient for an exact-prior or exact-estimator finding.

## Source group 1: Lindley and El-Sayyad (1968)

[D. V. Lindley and G. M. El-Sayyad, *The Bayesian Estimation of a Linear Functional Relationship*](https://academic.oup.com/jrsssb/article/30/1/190/7026880), Journal of the Royal Statistical Society, Series B, 30(1), 190-202 (1968), DOI 10.1111/j.2517-6161.1968.tb01519.x. The publisher's issue date establishes 1968; its later online migration date is not the historical publication date.

Access resolved: the parent retrieved the publisher's original PDF through normal browser download, saved locally as `C:/Users/canva/Downloads/jrsssb_30_1_190.pdf`. All 13 PDF pages, printed pp. 190-202, were read by text extraction; printed pp. 193 and 195-199 were also rendered and the decisive equations and prior passages visually checked. No temporary signed download URL is retained. This supersedes the earlier access-limited entry.

### Native model and exact posterior

Printed p. 190, equation (1), already fixes the intercept at zero: eta=theta xi. It observes independent pairs

$$
x_i\mid\xi_i\sim N(\xi_i,\phi_1),\qquad
y_i\mid\xi_i,\theta\sim N(\theta\xi_i,\phi_2),
$$

with conditional independence within and between pairs. Both error variances are unknown and common across observations. The paper leaves the slope sign unrestricted, explicitly acknowledges that positive sign is known in many applications, and returns to that possibility on p. 198. Thus through-origin geometry and awareness of the positive-slope case are already published features; they should not be listed as absent or as new adaptations.

Printed p. 193, Section 3, assigns xi_i|tau independently N(0,tau), with a generic joint prior pi(theta,phi1,phi2,tau). Equation (8) gives the integrated observation covariance

$$
\begin{pmatrix}
\tau+\phi_1 & \theta\tau\\
\theta\tau & \theta^2\tau+\phi_2
\end{pmatrix}.
$$

Writing s11=sum(x_i^2)/n, s12=sum(x_i y_i)/n and s22=sum(y_i^2)/n, equation (9) gives the exact posterior density, up to a data-only constant,

$$
\pi(\theta,\phi_1,\phi_2,\tau)
D^{-n/2}\exp\left[-\frac{n}{2D}
\{s_{11}(\theta^2\tau+\phi_2)-2s_{12}\theta\tau+s_{22}(\tau+\phi_1)\}\right],
\quad D=\phi_1\phi_2+\tau\phi_2+\tau\phi_1\theta^2.
$$

This retains relative signs through s12 theta. It is an exact generic Bayesian model; the later approximation does not exhaust its mathematical scope. Treating n=d coordinate pairs as the two input vectors gives the target likelihood in the special covariance case C=diag(phi2 I_d,phi1 I_d), with the block order (F,a), after conditioning on known phi1 and phi2. Arbitrary correlated or direction-dependent supplied C is not the instantiated observation law in this paper. Extending the likelihood is straightforward Bayesian modeling, but should be identified as an extension.

### The complete prior matters

The paper's explicit choices are not the target's flat variance mixing measure:

- Printed p. 195, before equation (13): theta, tau, and the pair of error variances are independent; arctan(theta) is uniform on [-pi/2,pi/2), giving density proportional to 1/(1+theta^2); log(tau) is uniform, giving d(tau)/tau.
- Printed p. 195, equation (14): phi=phi1, lambda=phi2/phi1, rho=tau/phi1. A phi^-1 prior independent of the remaining parameters permits the integration producing equation (15).
- Printed p. 198, final paragraph, specifies the Section 4 example more fully: theta has a constant prior, rho has density proportional to rho^-1, the two are independent and independent of lambda, whose prior is left as pi(lambda). At fixed phi this again means d(tau)/tau, not d(tau). Jointly, the phi and rho measures transform to d(phi)d(tau)/(phi tau).

Our independent variance-mixture calculation below shows why this is a substantive distinction. Mixing N_d(v;0,tau I) over d(tau) gives the target radial power ||v||^(2-d) for d>2. Mixing it over d(tau)/tau instead gives Gamma(d/2) pi^(-d/2)||v||^(-d) for nonzero v. Therefore the explicitly selected hyperpriors do not reproduce the target measure.

There is also a useful mathematical caution about treating the printed vague prior as an exact unbounded prescription. With positive fixed measurement variances, the likelihood integrated over xi tends to a strictly positive finite density as tau tends to zero. Its integral against d(tau)/tau consequently diverges at that endpoint. This is an inference from equation (9), not a claim made by the authors. It does not invalidate choosing suitable proper priors: printed p. 201 explicitly permits replacing their improper structural priors with sufficiently smooth proper priors. It does mean that the literal log-uniform expression and the paper's large-sample approximations should not be promoted to a proved exact posterior for every finite input.

**Exact containment result, derived here:** keep the generic hierarchy of equation (9), condition on known phi1 and phi2, set n=d>2, and choose pi(theta,tau) dtheta dtau proportional to 1(theta>0) dtheta dtau. Integrating tau then gives exactly the target joint prior and likelihood for the independent homogeneous-error covariance subcase. Thus the target is an exact specialization of the paper's general family under this different hyperprior. That fact is stronger than mere conceptual similarity, but it is not evidence that the authors selected this hyperprior or the target point functional.

For d=2, the target is directly flat in theta and the latent vector. Printed p. 201 explicitly discusses the earlier use of uniform incidental coordinates, independent of theta, with known error variances, and compares this with uniform radial-distance assumptions. That is direct historical evidence that flat-latent prior choices and their symmetry consequences were already being studied. This passage does not specify the complete target positive-slope prior plus point summary. The flat-variance mixture above must not be used in d=2, where it diverges.

### Published estimates and uncertainty

Printed pp. 196-197, equations (19)-(21), yield a large-sample normal approximation for theta conditional on lambda. Its center theta_lambda is the root, having the same sign as s12, of

$$
\theta^2+t\theta-\lambda=0,
\qquad t=(\lambda s_{11}-s_{22})/s_{12},
$$

and its variance is theta_lambda^2 u^2/n, where u^2=s11 s22/s12^2-1. The text identifies this center with the usual least-squares/maximum-likelihood estimate for known lambda. Printed p. 199, equation (25), gives a numerical mixture of these approximate normal laws over lambda; equation (26) offers an approximate log-slope form.

These are not the target's exact norm-weighted expectation ratio E[theta||xi||]/E[||xi||]. In particular, the displayed approximation is undefined at both-zero observations because s12=0; negative s12 chooses a negative center. The latter is consistent with the paper's native unrestricted-sign model and is not evidence that a positive-restricted exact posterior could not be constructed. The paper expressly confines these calculations to large or moderately large samples. No complete prior and point combination identical to the target was found in the inspected full article.

**Supported verdict:** the explicit priors and published approximate point/uncertainty prescriptions differ, with demonstrated limitations for the broad all-input contract. The general exact hierarchy contains the target d>2 independent homogeneous-error subcase after the specific alternative flat-variance and positive-slope choices described above. Exact identity of the full arbitrary-C target prior, point and mass law is not established as an authored method in this paper. This is a bounded primary-source result, not a historical uniqueness proof.

## Source group 2: Dellaportas and Stephens (1995)

[Petros Dellaportas and David A. Stephens, *Bayesian Analysis of Errors-in-Variables Regression Models*](https://www.jstor.org/stable/2533007), Biometrics 51(3), 1085-1095 (1995), DOI 10.2307/2533007. The [author's publication page](https://petrosdel.github.io/publication/dellaportas-1995-bayesian/) confirms authors, title and journal; Kelly (2007) cites the paper directly in its Bayesian errors-in-variables background.

Access resolved: the parent read all 11 original printed pages, pp. 1085-1095, through the authenticated Imperial JSTOR viewer at the stable article URL above. No Dellaportas PDF was downloaded: the download route required accepting additional terms, whereas the article was readable in the viewer. The following equation and page findings are the parent's direct full-text inspection, communicated to this lane for the shared record. They are not inferred from secondary citations or from the abstract. This supersedes the preliminary access-limited result.

### Framework and displayed algorithms

Printed p. 1086, equation (4), states the joint latent-variable posterior in the form

$$
[\theta,X_a\mid Y,X_c]\propto
[Y\mid X_a,\theta,X_c][\theta,X_a\mid X_c].
$$

Equation (5) states the conditional independence structure, and printed p. 1087, equation (7), gives the Gibbs/full-conditional computation. This is a general framework for inferring model parameters together with uncertain covariates, and is adaptable to other likelihoods and priors. It does not itself select a unique prior on a positive vector proportionality manifold.

Printed p. 1089 discusses an improper flat prior over the five real logistic parameters and uses normal priors. The second example, printed pp. 1091-1093, has ordered latent oxygen values and a nonlinear response: phi1 and phi2 have lognormal priors, other priors include normals and gammas, and phi3 has a flat prior. Thus it would be incorrect to say either that all native priors are proper or that the paper has no treatment of positive parameters. These are fully specified priors for its actual examples, not the target's radial measure.

### Nonlinear quantities and uncertainty

Printed p. 1090 calculates relative log potency rho=theta5/theta4 by dividing paired posterior draws, and then potency=exp(rho). This constructs the posterior law of a ratio while preserving the dependence of numerator and denominator; it is not a quotient of their posterior expectations. Printed p. 1094 discusses uncertainty from the MCMC posterior output.

This directly anticipates a general algorithmic building block of the target's mass-uncertainty calculation: if the posterior in its latent variables can be sampled, transform each paired draw to obtain the ratio law and its quantiles. It therefore weakens any claim that propagating a dependent latent posterior through a ratio is a new algorithm. It does not establish the exact target prior, its norm-weighted point, or its all-input vector prescription.

**Supported verdict:** the exact posterior/Gibbs framework is a generic adaptable predecessor, and the paper explicitly computes a ratio distribution from paired posterior draws. Its displayed models and priors do not supply the target measure ||v||^(2-d) dm d^dv together with E[m||v||]/E[||v||], full known arbitrary-C two-vector likelihood, and one prescription at all readings. The target could be implemented within this general Bayesian sampling framework after its likelihood and prior are supplied, which is generic algorithmic containment rather than a demonstrated identity of the published estimators. No complete target combination was found in the parent's inspection of the full original article. Printed p. 1095 cites Lindley-El-Sayyad (1968), confirming the citation-trail relationship.

## Independently derived variance-mixture identity

This calculation does not depend on obtaining either paper and is not attributed to their authors. For d>2 and v nonzero,

$$
\int_0^\infty \phi_d(v;0,tI_d)\,dt
=\frac{\Gamma(d/2-1)}{2\pi^{d/2}}\|v\|^{2-d}.
$$

To check it, substitute s=||v||^2/(2t) in the integral of (2 pi t)^(-d/2) exp[-||v||^2/(2t)]. The remaining integral is Gamma(d/2-1). Thus the target prior in dimensions d>2 also has the hierarchical representation

$$
v\mid t\sim N_d(0,tI_d),\qquad
\pi(t)\,dt\propto dt,\qquad
\pi(m)\,dm\propto\mathbf1_{m>0}dm.
$$

Both hyperprior measures are improper; this is an identity of measures up to a constant, not a claim that either is a probability distribution. In d=3 the integrated density is exactly 1/(2 pi ||v||). Its point singularity is locally integrable in three-dimensional Cartesian volume. Nonnegative integrands permit the order of these integrations by Tonelli; posterior propriety still needs the target's separate existence argument.

This identity does not extend to a flat variance mixing measure in d=2: the variance integral diverges logarithmically at infinity. The d=2 target measure is directly flat in (m,v), as established by the original Jacobian. Nor is a flat variance measure dt interchangeable with a flat standard-deviation measure or a log-uniform variance measure dt/t; those yield different radial powers. The centering assumption v|t~N(0,tI) also matters: integrating an unknown latent mean is a different model.

The radial function ||v||^(2-d) has zero Laplacian away from zero: for a radial power r^p, the d-dimensional Laplacian is p(p+d-2)r^(p-2). Calling the measure a particular radial or harmonic form is therefore mathematically justified, but this lane has not opened a separate historical source group on harmonic priors.

## Consequence for the contribution draft

The exact comparison should inspect complete hierarchical priors after marginalization, not stop at the displayed conditional Gaussian. A prior difference cannot be established merely by observing that a predecessor uses a latent Gaussian distribution. The target's prior can itself be written that way in d=3.

A defensible contribution statement can describe the particular combination of physical positive proportionality, chosen measure, finite weighted point and uncertainty law, together with its proved domain behavior. It should not call the radial prior itself new, or claim that no earlier equation or algorithm is equivalent. Lindley-El-Sayyad already supplies an exact general hierarchy containing a target covariance subcase under a selectable different hyperprior, and discusses flat incidental priors and coordinate/radial symmetry. Dellaportas-Stephens supplies a general posterior sampling framework and explicitly obtains a ratio law from paired posterior draws. Neither inspected article instantiated the full target combination. Historical priority remains narrower than general Bayesian model containment and is not settled by this bounded check.

This lane is capped at these two groups. Further retrieval is limited to the already requested two primary articles; no outward citation branch or new candidate is authorized within this lane.
