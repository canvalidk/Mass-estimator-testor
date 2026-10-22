# Huard and Mailhot (2006): Bayesian input-uncertainty framework

Date: 2026-09-08. Status: local evidence note, awaiting transfer; NOT submitted to VD-docs. Scope: one new published candidate, its general framework, and an explicitly identified mathematical instantiation. The combined review records the managed target snapshot; no managed source is modified here.

Source: David Huard and Alain Mailhot, *A Bayesian perspective on input uncertainty in model calibration: Application to hydrological model "abc"*, Water Resources Research 42, W07416 (2006), DOI 10.1029/2005WR004661. [Institutional full text](https://espace.inrs.ca/9481/1/P1330.pdf).

## What the publication supplies

Section 2.2, Eq. (1), treats true inputs and outputs as vectors. Section 2.3.1, Eq. (3), p. 4, integrates nuisance input/output values to obtain the parameter posterior. Before its final factorization, it has the form

$$
p(\theta\mid\widetilde x,\widetilde y,M)
\propto\iint p(\widetilde x,\widetilde y\mid x,y)
p(y\mid x,\theta,M)p(x,\theta\mid M)\,dx\,dy.
$$

The accompanying paragraph explicitly says that conditional independence of input and output errors is a simplifying assumption, not mandatory. Section 3.2's example is already through the origin, $y=\theta x$, with $\theta\sim U[0,5]$, $x\sim U[0,10]$, Gaussian input error and exact output. Sections 3.3.2-3.3.4 compare input priors and parameter posteriors. Example point estimates emphasize the most probable slope. The publication does not print the lognormal-prior, posterior-mean construction below.

## Our instantiation of that published general framework

This is a concrete competitor construction, not a transcription of the worked example. Fix, once for all inputs, positive constants $\tau,m_0,s$. Use the supplied coordinates without fitting a rotation or translation. Let

$$
x\sim N_d(0,\tau^2I),\qquad
\log(m/m_0)\sim N(0,s^2),\qquad y=mx,
$$

with independent priors, and measured $z=(a^T,F^T)^T$ satisfying

$$
z\mid x,m\sim N_{2d}(B_mx,\Sigma),\qquad
B_m=\begin{bmatrix}I\\mI\end{bmatrix},\quad\Sigma\succ0.
$$

The arbitrary known full joint covariance is admitted by the unfactorized likelihood in the published equation. The proper Gaussian input prior, full-support positive lognormal mass prior and posterior-mean decision are choices made in this audit. Applying the general integral gives the closed-form marginal likelihood

$$
q_m(z)=N_{2d}\!\left(z;0,\Sigma+\tau^2B_mB_m^T\right),
\quad
p(m\mid z)=\frac{q_m(z)\pi(m)}{Z(z)},
\quad Z(z)=\int_0^\infty q_t(z)\pi(t)\,dt.
$$

Take $\widehat m=E[m\mid z]$ and report posterior quantiles. This uses one prescription at zeros, across angles, and for every supplied SPD covariance. Neither a hypothesis-rejection gate nor an angle-dependent repair is introduced.

## Our existence and coverage proof

For fixed finite $z$ and $\Sigma\succ0$, the conditional Gaussian likelihood is strictly positive and bounded above by

$$
C_\Sigma=(2\pi)^{-d}|\Sigma|^{-1/2}.
$$

Integrating against the proper input prior preserves $0<q_m(z)\le C_\Sigma$. Consequently $0<Z(z)\le C_\Sigma<\infty$, including when either or both measured vectors are zero. Since the lognormal prior has finite mean,

$$
0<E[m\mid z]
\le \frac{C_\Sigma E_\pi[m]}{Z(z)}<\infty.
$$

The posterior density is positive on all $(0,\infty)$. Every quantile with probability strictly between zero and one is therefore strictly positive and finite. The same bounds, with dominated convergence locally in $z$ and SPD $\Sigma$, give continuity of the normalizer, posterior mean, and interior quantiles. An existence proof covers the entire finite-data domain; numerical checks can assess implementation rather than establish universality.

Relative signs are used. For the simple covariance $\Sigma=I_{2d}$, write $D=1+\tau^2+m^2\tau^2$. The marginal likelihood is proportional to

$$
D^{-d/2}\exp\left[-\frac{
(1+m^2\tau^2)\|a\|^2+(1+\tau^2)\|F\|^2
-2m\tau^2a^TF}{2D}\right].
$$

Flipping one nonorthogonal measured vector changes this expression. The general covariance expression retains the complete coordinates and cross-covariance. The priors have support at every finite true vector and every positive mass, so the construction imposes no upper mass cutoff. Prior influence when data are weak is part of this specified inference, and is not grounds for excluding a competitor when priors are allowed.

The printed $U[0,5]$ example must be treated separately: its posterior can only support masses within that bound, and its exact-output assumption is outside the SPD test domain. Those restrictions do not refute the general framework. A bounded positive prior with a noisy full-rank likelihood would still give finite positive posterior means for all observations, although it would not recover arbitrarily large, well-measured mass ratios. The full-support instantiation above avoids that upper bound.

## Relationship to the target and verdict

**Full coverage established mathematically for our explicit instantiation of the published general framework.** This result is not a claim that the paper's printed worked example already implements the entire contract. The historical statement supported by the source is that the relevant joint-likelihood, latent-input integration framework was published in 2006. The chosen prior and point prescription remain explicitly attributed to this audit.

The constructed estimator is a posterior mean of mass. It is not generally the target $E[f]/E[\alpha]$, and its prior differs from the target's flat measure $df\,d\alpha\,d\Omega$. With $x=\alpha n$ and $f=m\alpha$, our joint prior in those variables is proportional to

$$
\alpha^{d-2}e^{-\alpha^2/(2\tau^2)}\pi(f/\alpha)
\,df\,d\alpha\,d\Omega,
$$

which is not flat. Conversely, the target's flat measure can formally be represented in the generic Bayesian framework using density $p(x,m)\propto\|x\|^{2-d}$ with respect to $dx\,dm$, $m>0$. This Jacobian identity is our observation, not a prior choice or point statistic found in Huard and Mailhot. It supplies no exact historical identity claim.

What this adds: an earlier general Bayesian method can generate a complete functional competitor once the priors and decision rule are specified. This materially narrows any originality claim based solely on finite positive output, sign use, uncertainty and zero coverage. The specific target measure, pushforward distribution and ratio-of-expectations decision require their own historical comparison.

Coverage: primary text sections 2.2-2.3 and 3.1-3.4 inspected, with exact equation/paragraph locators above. No broad priority search or author contact. Numerical verification of this construction is delegated to the combined review; no numerical results are claimed in this note.
