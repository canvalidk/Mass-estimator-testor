# Repeated-sampling coverage and local information

Working mathematical contribution, 2026-09-15. Scratch only. This continues
the Gaussian compatible-pair investigation and introduces no rejection gate.

## Model and the two probabilities

At a fixed true positive mass m and true acceleration vector a, observe

\[
Y=(\widetilde F,\widetilde a)\sim
N\!\left(\begin{bmatrix}ma\\a\end{bmatrix},\Sigma\right),
\qquad\Sigma>0.
\]

For independent isotropic channels,
`Sigma=diag(sigma_F^2 I_d,sigma_a^2 I_d)`. The covariance is assumed known and
correct. Denote the sampling law by P_(m,a).

A conditional/credible mass interval I(y) satisfies, for example,

\[
P(M\in I(y)\mid y;\text{declared latent weighting})=.95.
\]

Repeated-sampling coverage instead is

\[
C(m,a)=P_{m,a}\{m\in I(Y)\}.
\]

The first probability varies latent states after fixing the data. The second
varies datasets after fixing the true mass and acceleration. One does not
imply the other. Proper Bayesian joint models can give prior-predictive
average coverage identities; those are not pointwise uniform coverage, and
an improper reference weighting does not by itself supply that proper
prior-predictive distribution.

## 1. A coverage obstruction at zero excitation

At a=0,

\[
P_{m,0}=N(0,\Sigma)=P_0\quad\text{for every }m>0.
\]

Suppose a procedure reports an interval with measurable endpoints L(Y),U(Y).
If U is finite P_0-almost surely, then

\[
C(m,0)\le P_0\{U\ge m\}\longrightarrow0
\quad(m\to\infty).
\]

Likewise, if L>0 P_0-almost surely,

\[
C(m,0)\le P_0\{L\le m\}\longrightarrow0
\quad(m\downarrow0).
\]

Therefore any procedure that always supplies a finite bounded interval
strictly inside (0,infinity) has no positive uniform lower coverage over the
unrestricted mass parameter at true zero excitation. This applies to any
method, not specifically the ratio-of-means construction.

The assertion is **not** that every fixed mass has poor coverage, that all
weak-signal intervals are useless, or that a physically bounded mass domain
must have unbounded intervals. It is a uniform statement over unbounded
positive mass when arbitrarily weak excitation is allowed.

### What genuinely uniform coverage would require

If an interval has coverage at least c for every m>0 at null, allowing
endpoints 0 and infinity, then

\[
P_0(U=\infty)\ge c,\qquad P_0(L=0)\ge c.
\]

Hence it must report the entire positive mass domain with probability at
least `max(0,2c-1)`. For 95% coverage that lower bound is 90%. This is a
necessary condition, not a recommendation to install a new output rule.

## 2. Removing exactly zero from the parameter domain does not solve it

For two Gaussian laws with the same covariance, let

\[
D^2=\begin{bmatrix}ma\\a\end{bmatrix}^{\!T}
\Sigma^{-1}\begin{bmatrix}ma\\a\end{bmatrix}.
\]

Their total variation distance, defined as the largest difference in event
probabilities, is exactly

\[
\|P_{m,a}-P_0\|_{TV}=2\Phi(D/2)-1
\le {D\over\sqrt{2\pi}}.
\]

The exact expression follows by dividing the two Gaussian densities: the
set where one exceeds the other is a half-space perpendicular to their
whitened mean displacement. Its probabilities reduce to one-dimensional
normal tails. In particular, for any interval rule,

\[
|C(m,a)-C(m,0)|\le2\Phi(D/2)-1.
\]

For independent isotropic channels,

\[
D=\|a\|\sqrt{m^2/\sigma_F^2+1/\sigma_a^2}.
\]

Choose m_n tending to infinity and a_n=v/m_n^2 for any fixed nonzero vector v.
Then a_n is nonzero, both true channel means tend to zero, and D_n tends to
zero. If U is almost surely finite under P_0,

\[
C(m_n,a_n)\le P_0(U\ge m_n)+\|P_{m_n,a_n}-P_0\|_{TV}
\longrightarrow0.
\]

For the small-mass end, take m_n down to zero and a_n=m_n v. The corresponding
coverage also tends to zero if L is strictly positive almost surely.
Thus for any fixed delta>0, uniform positive coverage fails even on
`m>0, 0<||a||<delta` if interval endpoints are finite and positive.

The null observation law can be used as a limit even if the physical
parameterization excludes a=0. All finite-mean Gaussian laws with the same
positive definite covariance are mutually absolutely continuous, so the
almost-sure finite-endpoint assumption cannot be evaded by removing a single
mean from the parameter domain.

This argument extends to any fixed finite number n of repeated observations
of the same pair: replace Sigma by Sigma/n and D by sqrt(n)D. More data help
at nonzero fixed excitation, but no fixed n cures a requirement uniform over
arbitrarily small excitation and unrestricted mass.

## 3. Local Fisher information with an unknown true acceleration vector

The acceleration vector a is a nuisance parameter: it is part of the true
state needed to explain the measurements, while m is the desired parameter.
For independent isotropic noise the Fisher blocks are

\[
I_{mm}=\|a\|^2/\sigma_F^2,\qquad
I_{ma}=m a^T/\sigma_F^2,
\]

\[
I_{aa}=(m^2/\sigma_F^2+1/\sigma_a^2)I_d.
\]

Removing the information shared with the nuisance score by a Schur complement
gives

\[
\boxed{
I_{m\cdot a}=I_{mm}-I_{ma}I_{aa}^{-1}I_{am}
={\|a\|^2\over\sigma_F^2+m^2\sigma_a^2}.
}
\]

The efficient score is

\[
S_{m\cdot a}={a^T(\epsilon_F-m\epsilon_a)
                 \over\sigma_F^2+m^2\sigma_a^2},
\]

where epsilon_F and epsilon_a are the measurement errors at the true
parameter. Its variance is the displayed efficient information.

At a nonzero regular parameter point, the Cramer-Rao lower bound for a
locally unbiased mass estimator is `1/I_(m.a)`. A leading small-noise/local
standard uncertainty is therefore

\[
\sigma_{m,\mathrm{local}}
={\sqrt{\sigma_F^2+m^2\sigma_a^2}\over\|a\|}.
\]

This is local sampling information, not the exact spread of the induced
latent ratio law. It does not promise finite-sample coverage or a finite
exact variance for every chosen estimator. At a=0 information is zero and
the regular mass-identification calculation is singular.

For log mass, let rho_F=m||a||/sigma_F and rho_a=||a||/sigma_a be the two
channel signal-to-noise ratios. Then

\[
I_{\log m\cdot a}=m^2 I_{m\cdot a}
={\rho_F^2\rho_a^2\over\rho_F^2+\rho_a^2},
\qquad
I_{\log m\cdot a}^{-1}=\rho_F^{-2}+\rho_a^{-2}.
\]

The weaker channel limits relative mass information. Repeating the same
independent experiment n times multiplies information by n. Independent
trials with different unknown true accelerations a_i and one common mass
add their efficient information when each trial has its own nuisance vector.

### Full known covariance

The corresponding exact local expression with correlations/anisotropy is

\[
\boxed{I_{m\cdot a}=a^TV(m)^{-1}a,\qquad
V(m)=\Sigma_{FF}+m^2\Sigma_{aa}
-m(\Sigma_{Fa}+\Sigma_{aF}).}
\]

To see it, the mean Jacobian for a is J_a=[mI;I], and the mass derivative is
J_m=[a;0]. The contrast K=[I,-mI] annihilates J_a and satisfies
`K Sigma K^T=V`, `K J_m=a`. Projection off the nuisance score yields the
displayed expression. It agrees with the existing aligned local uncertainty
implementation after reparameterizing (f,alpha,u) as (m,a).

## 4. What directional information remains in Gaussian summaries

With independent isotropic errors, profiling the unknown true a vector gives

\[
Q_{\mathrm{prof}}(m)
={\|\widetilde F-m\widetilde a\|^2
   \over\sigma_F^2+m^2\sigma_a^2}
={\|\widetilde F\|^2+m^2\|\widetilde a\|^2
       -2m\widetilde F\cdot\widetilde a
  \over\sigma_F^2+m^2\sigma_a^2}.
\]

The two squared norms and their dot product are the maximal invariant under
common orthogonal transformations of the data. With no external preferred
direction, a rotationally/reflection-invariant mass procedure can be written
using these three values. The dot product retains relative directional
agreement; keeping only norms discards it. This is an invariance statement,
not a claim that the three values are sufficient for estimating the full
oriented nuisance vector.

At a true nonzero compatible state a=||a||u, the first-order efficient mass
score uses the **longitudinal** residual
`u dot (epsilon_F-m epsilon_a)`. In independent isotropic noise the perpendicular
components have no additional first-order mass score after fitting direction.
They still participate in finite-data direction estimation and the global
likelihood geometry. Known direction and unknown direction have the same
leading local mass information in this particular isotropic model; that does
not equate their weak-signal conditional laws.

In particular, the measured magnitude ratio has the same leading linear
error `(u dot epsilon_F-m u dot epsilon_a)/||a||` at strongly informative
compatible data. One should not claim a first-order local-information gain
over that simple ratio merely because another procedure uses directions.
Any advantage must be established for specified finite-data behavior, richer
covariance, uncertainty output, or downstream use.

For repeated iid Gaussian measurements of one unchanged latent pair with
known covariance, their full vector sample mean and covariance Sigma/n
retain all mean-parameter likelihood information. Residual sample scatter is
ancillary in that declared model. This gives one precise boundary for the
later information-use comparison; it changes if covariance is estimated,
latent states change, or observations are correlated.

### Do not manufacture information from a parameter-dependent residual

At each candidate m, `r_m=Ftilde-m atilde` has covariance V(m), but r_m itself
is a parameter-dependent transformation of the fixed data. Treating it as an
ordinary fixed observed Gaussian vector and adding a `log det V(m)` term is
not the full model's profile likelihood. The actual profile is proportional
to `exp(-Q_prof/2)` because the original covariance determinant is constant.

Integrating over a using a flat Cartesian reference measure really does
produce an additional factor proportional to
`(sigma_F^2+m^2 sigma_a^2)^(-d/2)` in the isotropic case. That is a
nuisance-volume/reference-measure choice. It must not be
presented as extra sampling information about mass at a true null state,
where the observation law is exactly identical for every mass.

## Verification and limits

The coverage and total-variation results above are analytic. The local
information formula is an exact Fisher/Schur-complement calculation.
`verify_local_information.py` compares the full-covariance contrast formula
against direct nuisance elimination and the existing implementation over
random positive definite covariances in dimensions 1,2,3. It also checks the
isotropic and signal-to-noise expressions. No simulation here claims to
establish uniform coverage, and no rejection rule or estimator revision has
been added.
