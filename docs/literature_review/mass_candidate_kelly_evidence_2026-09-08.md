# Kelly (2007): fresh Bayesian errors-in-variables comparator

Date: 2026-09-08.
Status: local active research record; unsubmitted and awaiting transfer to VD-docs. Saving here does not submit it.
Scope: one fresh candidate outside the original comparison list, primary equations and a bounded existence/propriety audit. No earlier skipped candidate was reopened. Managed-context snapshot and target-source provenance belong to the parent comparison record.

## Primary source and anchors

[Brandon C. Kelly, *Some Aspects of Measurement Error in Linear Regression of Astronomical Data*](https://arxiv.org/pdf/0705.2774), Astrophysical Journal 665, 1489-1506 (2007), DOI 10.1086/519947. [Author-deposit metadata](https://arxiv.org/abs/0705.2774) identifies the accepted paper and its 18 May 2007 deposit. The full 39-page author PDF was readable. Publisher DOI access was unavailable, but the primary author copy was sufficient; no institutional-login request is needed.

Equation anchors in that PDF: equations (11)-(15), page 5, specify latent covariates, Gaussian regression and known correlated measurement errors; equations (16)-(18), pages 5-6, integrate the Gaussian-mixture latent distribution; section 4.3 generalizes the covariates; section 6.1, page 12, specifies priors. Observations are independent across objects. The default regression prior is uniform over intercept, signed slope and nonnegative intrinsic variance. The paper uses posterior simulation to describe parameter uncertainty.

## Equation-level comparison

Translate one object's measured covariate to a and response to F. Kelly's model is

$$
v\sim\sum_k\pi_k N(\mu_k,\tau_k^2),\qquad
\eta\mid v\sim N(b+mv,\sigma_{\rm int}^2),\qquad
(F,a)\mid(\eta,v)\sim N_2((\eta,v),C).
$$

Here C is known and can contain a nonzero cross-covariance. The signed latent relation retains the distinction between concordant and opposing measurements. No observed-vector length is divided by another in this likelihood.

For one Gaussian mixture component, integrating the latent variables gives the bivariate normal of equations (17)-(18):

$$
E(F,a)=(b+m\mu,\mu),\quad
V(m)=
\begin{pmatrix}
m^2\tau^2+\sigma_{\rm int}^2+C_{FF}&m\tau^2+C_{Fa}\\
m\tau^2+C_{Fa}&\tau^2+C_{aa}
\end{pmatrix}.
$$

This is already a genuine two-noisy-channel model. Its Gaussian mixture prior on the true covariate is a modeling choice, not grounds for dismissing it under a contract that permits alternative legitimate priors.

The physical exact-proportionality specialization fixes b=0 and sigma_int=0, and explicitly parameterizes m>0. The collapsed observed covariance remains positive definite because C is positive definite; a zero intrinsic variance does not make the observed likelihood singular. It does, however, remove a parameter and changes the sampler from its native free-intercept/free-scatter configuration.

## Our decisive propriety calculation

Simply restricting the native flat slope prior to m>0 does not automatically yield an all-data posterior. In the scalar one-object specialization with b=0, sigma_int=0, a proper latent prior v~N(0,tau^2), C=I and measured F=a=0,

$$
|V(m)|=1+\tau^2+\tau^2m^2,
\qquad L(m)=\{2\pi\sqrt{1+\tau^2+\tau^2m^2}\}^{-1}.
$$

Thus L(m) is asymptotic to a positive constant divided by m. A flat prior on m>0 has an infinite normalizer. This is a demonstrated scalar failure of that specific specialization, even with a proper Gaussian covariate prior. It is not a claim that every Kelly model or every fixed-dimensional vector specialization is improper.

The unconstrained native model has an additional one-object issue: with a flat intercept b, integrating over b removes the only response location constraint; the remaining marginal in a is independent of m. Integrating the flat slope then diverges. Native many-object regression should therefore not be described as automatically covering one input pair with every default parameter left free.

## Explicit successful adaptation

Fix a proper latent Gaussian or Gaussian mixture and a proper positive scale prior, for example a lognormal. Fix b=0 and sigma_int=0, and use the posterior median and central posterior quantiles. For any finite scalar readings and any known positive-definite C, the likelihood before latent integration is strictly positive and bounded by the finite maximum of N_2(0,C). Its integral against the proper latent prior is therefore strictly positive and bounded by the same constant. Integrating against the proper scale prior gives a finite positive normalizer.

Consequently the positive-scale posterior is proper for all scalar readings, its median and interior quantiles are finite and strictly positive, and these quantities vary continuously through zero readings and sign changes. A lognormal prior also gives a finite positive posterior mean, by bounding its numerator by the corresponding prior first moment. This is our existence argument for an explicit adaptation, not a theorem or default prior attributed to Kelly.

For physical vectors treated as d paired coordinate observations, Kelly's product likelihood accepts unequal and correlated errors within each coordinate pair but, as published, assumes independence between objects. That route does not represent arbitrary cross-coordinate errors. The companion Mantz audit removes this particular obstacle by treating both entire vectors as one multivariate object. We did not reject the general latent-Gaussian likelihood merely because the application used scalar response variables.

## Verdict and additional value

**Published Bayesian model family established as a relevant two-noisy-channel predecessor. Full requested positive-scale contract is not established for the native default model.** Native positivity is not imposed; the one-object default is underdetermined, and the explicit scalar flat-positive-slope specialization above is improper. **Full scalar coverage is established for the stated proper-positive-prior adaptation.** Arbitrary full vector covariance is addressed separately by the Mantz specialization.

This candidate adds a concrete structural errors-in-variables likelihood and exposes the importance of which regression parameters are fixed and which priors are proper. It does not establish an earlier publication of the target df d(alpha) dOmega measure or E[f]/E[alpha] point. No numerical implementation or onward reference search was performed in this bounded lane.
