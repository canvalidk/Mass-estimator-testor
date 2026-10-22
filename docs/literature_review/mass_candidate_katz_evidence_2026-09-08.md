# Katz: scalar boundary mechanism audit — 2026-09-08

Status: local active/unsubmitted research note. Saving here does not submit it to VD-docs. This bounded audit uses the current task contract and the primary-source material identified below; no managed VD-docs source was consumed in this lane. The coordinating review records its managed snapshot separately. Katz and the Chang corroborating manuscript count as **one candidate**.

Verdict, updated after institutional access: **full scalar-boundary coverage established for Katz's original normal-mean rule**, by the source and posterior argument below. **Demonstrated scope limitation as a stand-alone solution of the two-vector contract**: the paper develops point estimation for a scalar parameter restricted to an interval, including one-parameter exponential families. It does not supply the noisy-vector proportionality model, mass distribution, or quotient of posterior magnitude means. The earlier original-text access gap is now closed.

## Source facts and access

[Katz (1961), *Admissible and Minimax Estimates of Parameters in Truncated Spaces*](https://doi.org/10.1214/aoms/1177705146), Annals of Mathematical Statistics 32, 136–142. Public DOI, Project Euclid, and JSTOR routes initially did not yield readable original text. Later on 2026-09-08 the user's Imperial Library session successfully authenticated JSTOR. All seven original pages were read in the [JSTOR article viewer](https://www.jstor.org/stable/2237613); the normal-mean equation on p.139 was also checked visually against the printed page. No login information is stored in this research record.

### Original-text findings added in the third bounded round

- Pages 136–138 introduce a scalar one-parameter exponential family, squared-error loss, and a parameter bounded below (or, by the stated analogous development, above). Theorem 1 on p.138 gives an admissible estimate of the scalar mean in this family. Its scope is broader than just the normal example.
- Page 137 derives a sequence of Bayes rules with exponential priors on the nonnegative natural parameter and takes the diffuse-prior limit. The resulting rule is the mean under the corresponding flat restricted-parameter measure.
- Page 139 explicitly prints the unit-variance normal rule `delta(x)=x+exp(-x^2/2)/integral_{-infinity}^x exp(-t^2/2)dt`, equivalent to `x+phi(x)/Phi(x)`. Lemma 1 explicitly states strict positivity for every finite reading. The known-sigma expression below is its direct rescaling. The page also treats minimaxity; this audit does not independently certify every step of the original decision-theoretic proof.
- Pages 139–141 give further sufficient admissibility results for scalar interval parameter spaces and an additional restricted-normal shrinkage example. The occurrences of multiple observations are reduced to scalar sufficient statistics; they are not two measured vectors with an unknown common scale and direction.
- Page 142 completes the references. Across pp.136–142, no posterior mass ratio, quotient of two posterior means, arbitrary joint vector covariance model, or the target compatible-pair measure is specified. Thus the original-paper *scope* is now checked directly. The uncertainty and vector-composition deductions below remain our deductions, not newly discovered claims of Katz.

[Chang, Shinozaki and Strawderman, 2017 manuscript](https://m.sc.niigata-u.ac.jp/~hirukawa/seminar/niigata2017_program/Chang.pdf), PDF p. 2, equation (1.4), attributes the known-variance rule to Katz:

$$
\delta(x)=x+\sigma\frac{\phi(x/\sigma)}{\Phi(x/\sigma)}.
$$

The manuscript identifies it as generalized Bayes for a uniform positive-half-line prior and reports admissibility and minimaxity. Its own model, equations (1.1)–(1.3), p. 1, concerns independent nonnegative normal means with common variance and independent variance information. Its shrinkage extension, equations (1.5)–(1.7), p. 2, is separate from the scalar Katz rule; it permits negative estimates, motivating equation (2.6), p. 6. We are not transferring its dominance results to mass ratios.

## Independent reconstruction and deductions

Here $\phi$ and $\Phi$ denote the standard normal density and CDF. For a finite scalar reading $x$, known $\sigma>0$, and a flat prior on $\theta\ge0$, direct normalization gives

$$
p(\theta\mid x)=
\frac{\phi((\theta-x)/\sigma)}{\sigma\Phi(x/\sigma)}
\mathbf1_{\{\theta\ge0\}}.
$$

This posterior is proper for every finite $x$: $0<\Phi(x/\sigma)<1$. Its Gaussian tail makes every positive moment finite, and its positive density on $(0,\infty)$ makes its mean strictly positive. Integration gives the displayed Katz rule. Consequently, a negative reading does not produce a negative estimate, and zero is ordinary:

$$
\delta(0)=\sigma\sqrt{2/\pi}>0.
$$

Every interior posterior quantile is finite and strictly positive:

$$
Q(q\mid x)=x+\sigma\Phi^{-1}
\!\left[1-(1-q)\Phi(x/\sigma)\right],\qquad0<q<1.
$$

The expression is continuous throughout the stated finite-data, positive-variance domain. These conclusions are mathematical; naive floating-point evaluation at extreme negative readings can underflow or cancel and needs stable evaluation of the same function.

The posterior mode is $\max(x,0)$, whereas its mean is $\delta(x)>0$. Rejecting the Katz mean because the constrained mode can be zero would confuse two estimators. The source's decision-theoretic risk result also does not itself establish a confidence-interval coverage claim. The normalized posterior and quantiles above are reconstructed here; the original paper's explicit output is a point-estimation rule and its decision-theoretic analysis.

## Known-direction quotient: a separate construction

Suppose a true common oriented direction is supplied, so the two unknown magnitudes are $f>0$ and $\alpha>0$, with independent signed scalar readings $x_F,x_a$ and known errors $\sigma_F,\sigma_a>0$. A flat prior on the positive quadrant gives the product of the two reconstructed posteriors. It therefore permits

$$
m_*=
\frac{E[f\mid x_F]}{E[\alpha\mid x_a]}
=\frac{\delta(x_F;\sigma_F)}{\delta(x_a;\sigma_a)}>0.
$$

Both numerator and denominator are finite and strictly positive for all finite readings. At zero–zero, $m_*=\sigma_F/\sigma_a$. This is a quotient-of-Katz construction, **not an equation attributed here to Katz**.

For uncertainty in $m=f/\alpha$, its transformed posterior density is

$$
p_m(m)=\int_0^\infty \alpha\,p_f(m\alpha)p_\alpha(\alpha)\,d\alpha,
\qquad m>0.
$$

It is proper, positive on the positive half-line, and has finite positive interior quantiles. Its mean need not exist: since the denominator posterior has positive density at zero,

$$
p_m(m)\sim \frac{p_\alpha(0)E[f]}{m^2}\quad(m\to\infty).
$$

Thus $E[f]/E[\alpha]$ is not $E[f/\alpha]$. The heavy tail does not prevent the point summary or quantile uncertainty from meeting this **restricted** problem.

With correlated scalar errors, the analogous positive-quadrant truncated bivariate Gaussian remains proper and has finite positive component means, but those means generally are not separate univariate Katz functions. The joint truncation depends on the correlation. Such an adaptation needs explicit attribution as a construction, not a verified result of this source.

## Boundary of the vector comparison

The full task supplies measured vectors in common coordinates and their full positive-definite joint Gaussian covariance, with unknown common true direction. It does not supply two directly measured positive magnitudes. Norms of noisy Gaussian vectors are not generally normal; replacing vectors by norms also identifies aligned and anti-aligned pairs. A componentwise positive-orthant prior is different again: a permitted true vector such as $(1,-1,0)$ has a negative component.

These are precise obstacles to direct substitution into the cited scalar rule. They do not exclude a jointly derived latent-direction generalization. Whether an earlier publication already supplies that generalization remains unresolved by this audit.

What this adds: it verifies an earlier, all-finite-reading mechanism behind positive estimation at the scalar boundary, and distinguishes the ratio of posterior means from the generally heavy-tailed posterior ratio. It is not another version of Leonard's zero-scatter failure. The first version relied on the Chang corroboration; the institutional-access update closes that specific evidence gap without changing the scalar comparison verdict. Broader searches are recorded separately in the third review.
