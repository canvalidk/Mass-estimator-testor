# E10: Can the pure-noise limit select the starting law?

**Status:** exploration, 9 Nov 2026, Codex; user-authorized follow-up to T42 and
E9. The proposed averaged-null fixed-point test fails for every full-support
prior on positive mass in the class below. This is a conditional obstruction,
not a new ruling, a change of prior, or a claim against the existence goal.
**Checks:** `../tests/test_null_limit.py`; research quadrature in
`../met/null_limit.py`. Production equations are unchanged.

## 1. Standard mathematics and the question

The proposed requirement is **average preservation of the prior**:

\[
\mathcal T\pi=E_0[Q_\pi]=\pi,
\]

where each returned law is normalized before taking its average over independent
noise-only experiments. This tests the limiting distribution of T42's outputs;
it is not a claim about iterating the estimator or sequentially feeding its
outputs back as priors.

The standard Bayesian comparison is **Bayes plausibility**, or the posterior
martingale identity. It averages under the prior predictive distribution of a
specified joint model. Kamenica and Gentzkow (2011), section I.A, printed p.2594,
state that constraint; section I.B gives their signal characterization:
[primary paper](https://web.stanford.edu/~gentzkow/research/BayesianPersuasion.pdf).
For a normalized experiment f_theta(y), the elementary identity is

\[
p_\pi(y)=\int f_\theta(y)\pi(d\theta),\qquad
\int Q_\pi(B\mid y)p_\pi(y)\,dy=\pi(B).
\]

Replacing p_pi by the fixed-null sampling law P0 is an additional condition.
If zero true signal is **declared as known**, f_theta is independent of mass
and every prior is preserved samplewise. Here zero signal is the test world,
not a declaration supplied to the procedure. The working factor need not be
that declared-null likelihood.

A second standard connection is **concave log-mixture optimization**. If
Z_pi(y)=integral L_y(theta) pi(dtheta), then F(pi)=E0 log Z_pi has directional
derivative E0[L_y(theta)/Z_pi(y)]-1 towards a point mass, where integrable.
This is the population counterpart of the mixture-likelihood score conditions
reviewed by Watanabe and Ikeda (2015), sections 3.1--3.2, equations (5)--(7):
[primary paper](https://www.lisl.cs.tut.ac.jp/wkazuho/rep/EntMin-Mixing-MACH-wo-cover.pdf).
Their finite-sample mixture results are not being invoked as a theorem about
this population kernel. The proofs needed here follow below.

Z/T24 concerns exact **observed** zeros. It does not already impose this
averaged condition on random observations with zero true signal (O7).
O6's full-support restriction is also still provisional. Both distinctions
matter to the conclusion.

## 2. The exact pure-null specialization of T42

Use T42's independent scalar Gaussian readings, standardized channel covariance
I2, one common instrument ratio s, flat quadratic factors, and Bartlett power.
No sign/misfit weights or finite-nu ours32 factors are included. Put
theta=arctan(m/s), so positive mass is theta in (0,pi/2).

At pure null, T42 section 3 has kappa=q=0. Its returned random law converges
**in distribution**, not to a deterministic law or necessarily along an
individual infinite sample path, to

\[
Q_\pi(d\theta)\ \propto\ \pi_\theta(d\theta)
\exp[-2R^2\sin^2(\theta-\Phi)],
\qquad R^2\sim\chi^2_2.
\]

Phi is uniform on the projective circle, independently of R. The true mass
has disappeared from this pure-null experiment; the supplied prior has not.
With x=2theta and phi=2Phi modulo 2pi, discard the common factor exp(-R^2):

\[
L_x(T,\phi)=e^{T\cos(x-\phi)},\quad T=R^2,\quad
x\in A=(0,\pi),\quad \phi\sim\operatorname{Unif}(0,2\pi).
\]

This is the axial von Mises structure also used in T70. We will allow any
independent T>=0 with P(T>0)>0 and E exp(aT)<infinity for some a>0.
In particular T=beta R^2 works for **every finite beta>0**. Beta is a
comparison of corrections, not a proposed new input or a fitted choice.

For a probability measure p on A define

\[
Z_p(T,\phi)=\int_A L_y(T,\phi)p(dy),\quad
g_p(x)=E\frac{L_x(T,\phi)}{Z_p(T,\phi)},\quad
(\mathcal Tp)(dx)=g_p(x)p(dx).
\]

Tonelli gives integral g_p dp=1. Thus the test is g_p=1 p-almost everywhere.

## 3. No full-support fixed point on positive mass

**Claim.** Under the stated strength assumptions, no probability measure
with full support on A satisfies mathcal T p=p. Smoothness, a density,
reciprocal symmetry, and a chosen centre are not required.

First justify differentiation. At each x0 in [0,pi], full support supplies
a positive mass b in a delta-neighbourhood of x0 inside A. For
|x-x0|<epsilon, cosine's Lipschitz bound implies

\[
L_x/Z_p\le b^{-1}e^{T(\epsilon+\delta)},\qquad
|\partial_xL_x|/Z_p\le b^{-1}T e^{T(\epsilon+\delta)}.
\]

Choose epsilon+delta<a strictly. The exponential moment dominates both
expressions. Therefore g_p is continuously differentiable locally throughout
A and across its endpoints. If mathcal T p=p, continuity and full support
force g_p=1 throughout A; its extension must have g_p'(0)=0.

But pairing phi and -phi in its endpoint derivative gives

\[
g'_p(0)=E_T\frac1{2\pi}\int_0^\pi
T\sin\phi\,e^{T\cos\phi}
\left[\frac1{Z_p(T,\phi)}-\frac1{Z_p(T,-\phi)}\right]d\phi.
\]

For T>0, phi in (0,pi), and every y in A,

\[
\cos(y-\phi)-\cos(y+\phi)=2\sin y\sin\phi>0.
\]

Consequently Z_p(T,phi)>Z_p(T,-phi), and the derivative is **strictly
negative**. The domination above justifies the interchanges and finiteness.
This contradicts the required zero derivative. Endpoints are used only for
the extension argument; no probability at m=0 or infinity is admitted.

This rules out the whole full-support O6 family for this *additional test*;
it does not invalidate that family under the existing premises.

## 4. Why full-circle uniformity does not rescue the half-Cauchy

On the complete doubled-angle circle, uniform measure U is a fixed point by
rotational symmetry. It is the **unique full-support fixed point** there.
Here is a proof, including the uniqueness qualification.

Use k_x=exp[T(cos(x-phi)-1)], z_p=integral k_x dp and F(p)=E log z_p.
Since -2T<=log z_p<=0, F is finite. For each T,

\[
E_\phi z_p=e^{-T}I_0(T)=z_U.
\]

Jensen yields F(p)<=F(U). A full-support fixed point has g_p=1 everywhere
by the same local regularity argument. Thus
E(z_U/z_p)=integral g_p dU=1. The inequality log r<=r-1 gives
F(U)-F(p)<=0, so equality holds in Jensen. Hence z_p(T,phi) is constant
in phi for almost every positive T. Its Fourier coefficients are proportional
to I_n(T) times the coefficients of p. Since I_n(T)>0 for T>0, all
nonconstant coefficients of p vanish: p=U.

The full-circle law restricted to A is uniform in theta, the half-Cauchy
2s/[pi(s^2+m^2)] in mass. However, restriction and averaging do not commute:

\[
U(\cdot\mid A)
=\frac{E[Q_U(A)Q_U(\cdot\mid A)]}{U(A)},\qquad U(A)=1/2,
\]

whereas our proposed test averages Q_U(.|A) **without** the random weights
Q_U(A). The positive-domain update is exactly that latter law. Reweighting
the noise ensemble to repair the identity would change the question.
This is related to E9's restriction issues, but the obstruction here already
occurs on the simplest single-instrument slice.

Two controls prevent overstatement. Every point mass is fixed. An equal
mixture of any two interior points is also fixed, by the reflection exchanging
them. Reciprocal pairs x and pi-x give swap-symmetric examples. Thus dropping
full support produces many fixed points, not uniqueness. Conversely T=0
almost surely makes **every** prior fixed.

## 5. Corrections, order of averaging, and what is not selected

The obstruction persists for every fixed beta>0 multiplying T42's limiting
tilt, and for bounded positive strengths. It is not a special numerical
failure at the Bartlett coefficient.

The newer **raw flat Gaussian Gram correction has the same null limit**.
Let H_N be the Gram eigenvalue gap. T42 gives

\[
H_N/(2\sqrt N)\Rightarrow R,\quad J_N/N\to1,
\quad\lambda_{\max}(G_N)/N\to1.
\]

The doubled-angle concentration is eta H_N/4. Bartlett gives
H_N^2/(4J_N) -> R^2 in distribution, with its cap inactive asymptotically.
The Gram correction eta_G=H_N/lambda_max gives
H_N^2/(4lambda_max) -> R^2 as well. This does not establish equivalence away
from pure null or for T70/B22 with nontrivial reading weights.

If a different correction instead drives the angular concentration to zero
in probability, its law tends to the supplied prior in total variation:
for fixed strength T, a sufficient bound is
TV(Q_p,p)<=min(1,exp(2T)-1). Boundedness then also gives averaged convergence.
**Every prior is preserved in that case**, so it cannot select a principal
member. No such correction is adopted or claimed forced here.

Likewise, averaging the **bounded unnormalized kernel first** gives
E exp[T(cos(x-phi)-1)]=E[e^{-T}I0(T)], independent of x. Normalizing that
averaged factor times p returns p for every p. This is not the average of
the normalized returned laws. Unbounded exp[T cos] should not be averaged
without checking integrability; its expectation can diverge at T42's strength.

## 6. Numerical witnesses and limits

`met/null_limit.py` uses Gauss--Laguerre quadrature for R^2/2 and a periodic
trapezoid rule for phi; prior integrals are represented by supplied weighted
nodes. All reported calculations are deterministic, with refinement checks,
not a Monte Carlo study. The code is a research reference at ordinary finite
strengths, not an extreme-input-certified production method.

For uniform x on (0,pi), quadrature with 320 prior nodes, 128 radial nodes
and 1,024 phase nodes gives (rounded; refinement is not a rigorous error bound):

| beta | response g(0)=g(pi) | response g(pi/2) | derivative g'(0) |
| --- | --- | --- | --- |
| 0.5 | 1.279552 | 0.892781 | -0.752120 |
| 1 | 1.762437 | 0.794312 | -3.192416 |
| 2 | 2.808427 | 0.691894 | -13.029217 |

At beta=1 the total-variation displacement is approximately 0.0955: the
average law redistributes probability toward the ends of the positive-mass
domain. Endpoints in the table are limiting density **ratios**, not masses
assigned probability. Reflection symmetry preserves the centre of the
averaged law and does not preserve its shape.

An independent small-beta check follows by expanding the likelihood ratio.
Writing v(x)=(cos x,sin x), mu=integral v dp, and T=beta R^2,

\[
g_p(x)=1+4\beta^2(\|\mu\|^2-v(x)\cdot\mu)+O(\beta^3).
\]

For uniform p, mu=(0,2/pi): the endpoint and centre increments divided by
beta^2 tend to 16/pi^2 and 16/pi^2-8/pi, and g'(0)/beta^2 tends to -8/pi.
These independent expansion coefficients check the normalization and factor
of two in the doubled-angle representation.

Validation on 9 Nov: all 14 new checks passed, including an independent
reviewer's run. The focused regression command below passed **101 tests**:

```text
python -m pytest -q tests/test_null_limit.py tests/test_gaussian_information.py
  tests/test_deviance_calibration.py tests/test_starting_law_family.py
  tests/test_reading_qubit.py tests/test_audit.py
```

The command is shown wrapped for readability; pass the paths on one command
line. The full production suite was not rerun for this isolated research
module. Independent mathematical and premise reviews found no substantive
correction. Audit findings remain T6/T21/T27 untested and L24 unused.

## 7. What this says about the goal

| proposed claim | result and scope |
| --- | --- |
| Pure noise approaches a limiting distribution of returned shapes | passes in T42's stated flat Gaussian Bartlett model; it is random and prior-dependent |
| Averaged null neutrality selects an admissible full-support principal prior | fails for every prior in the nontrivial isotropic tilt class proved above |
| The full-circle averaged-update test selects uniform angle | passes among full-support laws on the whole circle; restriction to positive mass fails the same fixed-point test |
| Switching raw Bartlett to raw Gram tempering removes this obstruction | fails on the pure-null scalar slice: the limits coincide |
| Existing Z requires the extra averaged condition | not established; exact observed-zero invariance and random-null calibration differ |
| Another correction can satisfy the project's full commitments | unknown; not tested or selected here |

The failed proposal names a conflict among **average null neutrality,
positive mass, full support, and persistent isotropic random evidence**.
Positivity is bedrock; the averaged requirement is a new candidate, not a
premise to attribute to the user or the corpus. If it becomes a requirement,
the evidence correction is a concrete place to investigate. If the null
evidence vanishes, however, this route alone still does not select pi.

The heterogeneous extension is not needed to refute a universal fixed-point
selector: the full input type includes this one-instrument slice. O6 and O15
remain open. The useful next issue is the intended operational meaning of
null neutrality, not searching harder for a fixed prior inside this class.
