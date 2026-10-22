# Liseo scalar-ratio comparison: second bounded evidence pass

Date: 2026-09-08.
Status: local active research record; awaiting transfer to VD-docs. Saving here does not submit it.
Scope: primary-source equation audit and independent scalar derivation. Managed-context provenance belongs to the parent comparison record. This note makes no exhaustive historical-priority claim.

## Source and access

[Brunero Liseo, *Bayesian and conditional frequentist analyses of the Fieller's problem. A critical review*](https://www.researchgate.net/publication/5182203_Bayesian_and_conditional_frequentist_analyses_of_the_Fieller%27s_problem_A_critical_review), METRON LXI(1), 133-150 (2003). The full article text is readable on this author-upload page, which identifies Brunero Liseo as uploader. The direct PDF download link failed. Equations were read from the article's extracted full text and checked by independent change of variables below, not copied from the previous local literature assessment.

Source anchors: pp. 133-134 define the equal-variance independent Gaussian model and ratio theta=omega_1/omega_2; equation (10), p. 137, gives the Jeffreys prior in ratio/denominator coordinates; equations (12)-(13), p. 138, give its marginal ratio posterior. Other priors produce other posteriors; the identity below concerns equation (12) only. The paper discusses the distinction between Bayesian credible sets and frequentist coverage. It does not state the full supplied joint-covariance vector contract.

## Exact identity checked independently

Use F and a for scalar readings. Let v=sigma^2/n be their common measurement variance. A flat prior on the two true signed means gives posterior

$$
(F_*,a_*)\mid F,a\sim N_2((F,a),vI).
$$

With m=F_*/a_* and t=a_*, the transformation (m,t) -> (mt,t) has absolute Jacobian |t|. Thus the original flat-means prior becomes proportional to |t| in (m,t) coordinates, matching equation (10). The ratio density is

$$
p_R(m)=\int_{-\infty}^{\infty}|t|\,
\phi_2((mt,t);(F,a),vI)\,dt.
$$

Completing the square in t gives

$$
p_R(m)=\frac{e^{-(F^2+a^2)/(2v)}}{\pi(1+m^2)}
\left[1+\frac{s(m)}{\phi(s(m))}
\left(\Phi(s(m))-\tfrac12\right)\right],
\quad s(m)=\frac{a+mF}{\sqrt{v(1+m^2)}}.
$$

This is Liseo's equation (12), with equation (13)'s notation translated. Here phi and Phi denote the standard normal density and distribution function.

Restricting to the physical sign event E={F_*a_*>0} gives

$$
p_+(m)=\frac{\mathbf1_{m>0}p_R(m)}{P(E\mid F,a)}.
$$

To compare directly with the target scalar measure, write F_*=u f, a_*=u alpha with f,alpha>0 and u in {-1,+1}. Flat df d(alpha) and equal weights on u give the same measure as flat dF_* da_* restricted to the same-sign quadrants. Substituting f=m alpha produces Jacobian alpha; summing the two directions reproduces the |t| integral exactly. Therefore the positive-conditioned equation (12) is the target unknown-sign scalar mass distribution, provided these model and base-measure assumptions are indeed the target's choices. Positive conditioning is the comparison specialization, not a claim that Liseo's displayed equation already restricts m>0.

## Unequal and correlated scalar covariance: our derivation

For any known positive-definite 2x2 covariance C, replace vI in the integral by C. A flat original-means prior still gives N_2((F,a),C) and the same Jacobian. This proves a correlated scalar extension of the identity, but this extension was not found explicitly in the reviewed paper.

Every finite reading has 0<P(E|F,a)<1 because the Gaussian density is positive everywhere. The conditioned ratio is supported on (0,infinity), with a strictly positive continuous density, hence a positive finite median and all interior quantiles. The Gaussian pushforward and sign conditioning vary continuously over finite readings and positive-definite C. This is one prescription at zeros and adverse signs.

At F=a=0, write C_FF=v_F, C_aa=v_a and C_Fa=kappa. Before positive conditioning, the ratio is a Cauchy with location kappa/v_a and scale sqrt(det C)/v_a. This follows by expressing F_*=(kappa/v_a)a_*+epsilon with epsilon independent of a_*. After conditioning it is that Cauchy restricted to m>0. It becomes the familiar half-Cauchy of scale sqrt(v_F/v_a) only when kappa=0.

For any finite Gaussian means and positive-definite C, the ratio density has positive asymptotic coefficient times m^(-2): in the Jacobian integral substitute z=mt, so m^2 p_R(m) tends to integral |z| phi_2((z,0);(F,a),C) dz, which is positive and finite. Positive conditioning therefore leaves the posterior mean of m infinite. Finite quantiles remain valid; demanding a finite mean or standard deviation would reject the target scalar law too.

The target point E[f|E]/E[alpha|E] is nevertheless finite and strictly positive, since the Gaussian magnitudes have finite moments and E has positive probability. This point is not E[m|E]. Liseo's displayed ratio posterior alone does not establish that the paper proposed this particular ratio-of-expectations point estimator.

## Verdict and incremental learning

**Exact scalar distributional identity established** for Liseo equation (12) after positive conditioning under equal-variance Gaussian readings and the specified flat-means measure. **Full scalar quantile-based coverage established** for this specialization and our direct correlated-covariance extension. **Full vector coverage unresolved from this source**, and no published identity of the target point estimator is established.

The scalar analogy does not by itself settle higher dimensions: exact vector collinearity is a lower-dimensional subset of the Cartesian pair-of-vectors space, so a choice of measure on that subset is still needed. Simply conditioning a nonsingular Cartesian Gaussian on that probability-zero event does not uniquely specify it.

This source adds a firm equation-level scalar precedent and clarifies the role of prior measure and uncertainty interpretation. Another general scalar normal-ratio paper is unlikely to add much to vector coverage unless it changes one of those ingredients. No references branching from this paper were pursued in this bounded pass.
