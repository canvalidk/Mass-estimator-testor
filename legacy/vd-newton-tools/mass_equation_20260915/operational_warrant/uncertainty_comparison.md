# Flat area versus noise-tube weighting: an uncertainty comparison

Scratch derivation, 2026-09-15. Single vector pair, three dimensions, independent isotropic Gaussian channel errors. This note compares conditional laws; it does not select an estimator policy or claim that a narrower credible interval has better sampling coverage.

## 1. Model and two measures

Standardize the observations as `U=F_obs/sigma_F`, `V=a_obs/sigma_a`. Write compatible latent states as

\[
X=r\sin\theta\,u,\quad Y=r\cos\theta\,u,
\quad r>0,\quad 0<\theta<\pi/2,\quad u\in S^2.
\]

Let `R=tan(theta)`, `m=sR`, `s=sigma_F/sigma_a`, and `ell=log R`. The two reference weights relative to `df dalpha dOmega` are

* flat area: `rho=1`;
* standardized noise tube: `rho=f^2+alpha^2=r^2`.

Both use the same Gaussian likelihood. Introduce the observed Gram quantities

\[
A=\|U\|^2,\quad B=\|V\|^2,\quad C=U\cdot V,
\qquad
D(\theta)=A\sin^2\theta+B\cos^2\theta+2C\sin\theta\cos\theta.
\]

Thus `h=sqrt(D)=||sin(theta) U+cos(theta) V||`. Integration over the common direction gives a factor `4 pi sinh(rh)/(rh)`. Constants depending only on the observations cancel from both conditional laws.

## 2. Exact radial integration and its interpretation

Define

\[
J_1(h)=\int_0^\infty r e^{-r^2/2}\frac{\sinh(rh)}{rh}\,dr,
\qquad
J_3(h)=\int_0^\infty r^3 e^{-r^2/2}\frac{\sinh(rh)}{rh}\,dr.
\]

The normalized angle densities are proportional to `J1(h(theta))` and `J3(h(theta))`, respectively. Direct integration gives

\[
J_1(h)=\sqrt{\frac\pi2}\frac{e^{h^2/2}\operatorname{erf}(h/\sqrt2)}h,
\quad J_1(0)=1,
\]

\[
\boxed{J_3(h)=(1+h^2)J_1(h)+1.}
\]

For a short verification, put `Q(h)=h J1(h)=integral exp(-r^2/2)sinh(rh) dr`. Integration by parts gives `Q'=hQ+1`, so `Q''=(1+h^2)Q+h`; division by `h` gives the second identity, including its continuous value `J3(0)=2`.

Consequently tube weighting changes the flat mass law by the relative factor

\[
\boxed{K(h)=\frac{J_3(h)}{J_1(h)}
=1+h^2+\frac1{J_1(h)}
=E_{\rm flat}[r^2\mid\theta,U,V].}
\]

This factor is nonconstant whenever the conditional radial second moment varies with the candidate mass. The difference is therefore not merely a normalizing constant away from special data geometries.

### The relative factor is strictly increasing

For `h>0`, differentiating the conditional radial expectation gives

\[
K'(h)=\operatorname{Cov}_h\left(r^2,
r\coth(rh)-\frac1h\right)>0.
\]

The second function is strictly increasing in `r`: its derivative is `coth(z)-z csch(z)^2>0` for `z=rh>0`, since `sinh(z)cosh(z)>z`. The radial distribution is nondegenerate, making the covariance strictly positive. Thus the tube measure gives relatively more mass to angles at which the data support a larger compatible latent radius.

The small-`h` expansions are

\[
J_1(h)=1+\frac{h^2}3+\frac{h^4}{15}+O(h^6),
\qquad
J_3(h)=2+\frac{4h^2}3+\frac{2h^4}5+O(h^6),
\]

\[
K(h)=2+\frac{2h^2}3+\frac{2h^4}{45}+O(h^6).
\]

## 3. Measured zero and an exact ordering beyond zero

At `U=V=0`, the angle is uniform for both measures. Hence both have

\[
p_R(r)=\frac{2}{\pi(1+r^2)},\quad
p_\ell(\ell)=\frac1{\pi\cosh\ell},\quad
E[\ell]=0,\quad \operatorname{Var}(\ell)=\frac{\pi^2}4.
\]

Their **joint** laws differ: the polar radius has a `chi_2` law under flat weighting and a `chi_4` law under tube weighting, independently of the angle. Their radial second moments are `2` and `4`. Their respective magnitude means are `sqrt(2/pi)` and `(3/2)sqrt(2/pi)` in each channel. Their common ratio-of-means point is nevertheless `mhat=s`.

There is also an exact uncertainty comparison whenever the observed channel norms agree, `A=B`. In this case

\[
h^2=A+C\operatorname{sech}\ell.
\]

Both laws are symmetric in `ell`, and both point readouts equal `s`. Because `K` is increasing:

* If `C>0`, the tube-to-flat density ratio decreases with `|ell|`. Therefore `|ell|` under the tube law is stochastically smaller: all its central log-mass intervals are narrower, and all finite increasing moments of `|ell|` are smaller.
* If `C<0`, the ordering reverses. Tube weighting gives wider central log-mass intervals and larger increasing moments of `|ell|`.
* If `C=0`, `h` is constant in angle. The two mass laws are exactly identical half-Cauchys, even if the equal observed magnitudes are large.

The stochastic-order statement follows directly from the monotone density ratio on `|ell|>=0`: after normalization, that ratio crosses one at most once. It applies to equal observed magnitudes with any feasible dot product, not only exactly parallel or antiparallel observations. These are width comparisons, not accuracy comparisons.

## 4. Small-observation expansion through second order

Scale the observations as `(U,V)=delta(U0,V0)` and let `delta -> 0`. Then `A,B,C` are of order `delta^2`. Write `k=0` for flat weighting and `k=2` for tube weighting, and define

\[
c_k=\frac{k+2}{6},\qquad
\overline D=\frac{A+B}{2}+\frac{2C}{\pi}.
\]

The normalized angle density is

\[
p_k(\theta)=\frac2\pi\left[1+c_k\{D(\theta)-\overline D\}\right]
+O(\delta^4).
\]

The remainder is uniform in angle. Gaussian radial domination also permits integration against any fixed power of `|log tan(theta)|`. The tube law's first change from the null angle law is exactly twice the flat law's change.

### Log-mass mean and variance

Using

\[
D-\overline D=\frac{B-A}{2}\cos2\theta
+C\left(\sin2\theta-\frac2\pi\right),
\]

and the null identities

\[
E_0[\ell\cos2\theta]=-1,
\qquad
\operatorname{Cov}_0(\ell^2,\sin2\theta)=-\frac\pi3,
\]

gives

\[
\boxed{E_k[\ell]=\frac{c_k}2(A-B)+O(\delta^4),}
\]

\[
\boxed{\operatorname{Var}_k(\ell)=\frac{\pi^2}4
-\frac{\pi c_k}{3}C+O(\delta^4).}
\]

Explicitly:

| Quantity | Flat area | Noise tube |
|---|---|---|
| `E[log(m/s)]` | `(A-B)/6 + O(delta^4)` | `(A-B)/3 + O(delta^4)` |
| `Var[log(m/s)]` | `pi^2/4 - pi C/9 + O(delta^4)` | `pi^2/4 - 2 pi C/9 + O(delta^4)` |

For the identities, `cos(2theta)=-tanh(ell)` and `sin(2theta)=sech(ell)`. Integration by parts with `p0(ell)=sech(ell)/pi` gives the first. For the second, `E0[ell^2 sech(ell)]=pi/6`, whereas `E0[ell^2]E0[sech(ell)]=pi/2`.

**Operational content:** to second order, unequal observed squared magnitudes shift the center; their dot product changes the log spread. Tube weighting contracts twice as fast under positive observed alignment and expands twice as fast under negative alignment. It is not uniformly a narrower uncertainty law.

### Equal-tail quantile changes

For any fixed `0<p<1`, let `ell_p` be the `p` quantile and `theta_p^0=pi p/2`. Integrating the angle-density expansion and inverting the CDF yields

\[
\boxed{\ell_{p,k}=\log\tan\frac{\pi p}{2}
+\frac{c_k}2(A-B)
-c_k C\frac{1-\cos(\pi p)-2p}{\sin(\pi p)}
+O(\delta^4).}
\]

The center-shift term is common to all quantiles. The dot-product term has opposite signs at complementary quantiles. Positive `C` moves the lower and upper equal-tail endpoints inward; negative `C` moves them outward. The tube changes are twice the flat changes at this order.

### Ratio-of-means readout

The reported point must be calculated separately; it is not generally `exp(E[ell])`. The relevant unnormalized first radial moments have quadratic coefficient `(k+3)/6`, because a `chi_(k+2)` radius satisfies `E[r^3]/E[r]=k+3`. Also

\[
\int_0^{\pi/2}(\sin\theta-\cos\theta)D(\theta)\,d\theta
=\frac{A-B}{3}.
\]

Therefore

\[
\boxed{\frac{\widehat m_k}{s}
=1+\frac{k+3}{18}(A-B)+O(\delta^4).}
\]

The flat coefficient is `1/6`; the tube coefficient is `5/18`. Thus the tube law's quadratic point response is `5/3` times the flat response, while its log-law mean and variance responses are twice the flat responses. These distinct coefficients are a useful safeguard against conflating the readout with the ratio law.

## 5. The intermediate radial weight `r^2+2`

One covariance-inspired residual is `||X cross Y||/sqrt(r^2+2)`. Its fixed small-norm tolerance induces the reference weight `rho=r^2+2`. The angle weight is exactly

\[
J_3(h)+2J_1(h)=(h^2+3)J_1(h)+1.
\]

For the unnormalized first magnitude moments, put `J2(h)=sqrt(pi/2)exp(h^2/2)`. Direct differentiation gives `J4(h)=(h^2+3)J2(h)`, so the first-radial-moment weight is

\[
J_4(h)+2J_2(h)=(h^2+5)\sqrt{\frac\pi2}e^{h^2/2}.
\]

The weak-data formulas consequently become

\[
E[\ell]=\frac{A-B}{4}+O(\delta^4),\qquad
\operatorname{Var}(\ell)=\frac{\pi^2}4-\frac{\pi C}{6}+O(\delta^4),
\]

\[
\frac{\widehat m}{s}=1+\frac7{30}(A-B)+O(\delta^4).
\]

More strongly, let `Zflat,Ztube` be the two joint likelihood integrals with their stated reference weights. The intermediate full conditional law is exactly

\[
P_{r^2+2}=\frac{Z_{\rm tube}P_{\rm tube}+2Z_{\rm flat}P_{\rm flat}}
{Z_{\rm tube}+2Z_{\rm flat}}.
\]

This also holds for their mass marginals. At measured zero it is an equal mixture, and its radius is an equal mixture of `chi_2,chi_4`. For equal observed norms, its `|ell|` law lies between the exact flat/tube stochastic orderings in Section 3. For unequal norms one must not infer that its log variance lies between component variances: their log means can differ, producing an additional between-component variance.

The covariance interpretation needs a precise distinction. For independent unit Gaussian channel errors around candidate means `X,Y`, the cross product's sampling covariance is

\[
\mathcal C(X,Y)=(r^2+2)I-XX^T-YY^T.
\]

The candidate residual `B=X cross Y` is an eigenvector of this candidate covariance field with eigenvalue `r^2+2`. Thus the algebraic identity `B^T C(X,Y)^(-1) B=||B||^2/(r^2+2)` is exact. At fixed compatible true means `xu,yu`, however, the fixed sampling covariance is `(r^2+2)(I-uu^T)+2uu^T`; a random observed cross product need not be perpendicular to the fixed true `u`. Candidate-field scaling is therefore distinct from whitening that observed random vector by its fixed true covariance. The cross-product distribution is non-Gaussian, so neither identity supplies an exact chi-squared pivot or a coverage guarantee.

## 6. Tails and strong-signal scope

At every finite observation, `J1(h(theta))` and `J3(h(theta))` are continuous and strictly positive at both angle endpoints. Since `dtheta/dR=1/(1+R^2)`, both ratio densities tend to a positive constant at zero and decay as a positive constant times `R^(-2)` at infinity. Thus both ordinary ratio means diverge, as do inverse-ratio means; their log moments are finite. A strongly concentrated central interval does not remove these mathematical tails.

At fixed positive true mass with compatible standardized means `t f0 u, t alpha0 u`, `f0,alpha0>0`, and unit observation noise, both weights are locally smooth and differ by a factor whose relative variation over an order-one likelihood neighborhood is `O(1/t)`. They therefore share the first-order local Gaussian law and the leading log-mass variance

\[
\frac1{t^2 f_0^2}+\frac1{t^2\alpha_0^2}+o_p(t^{-2}).
\]

The magnitude coordinates have independent unit-noise longitudinal fluctuations at this order; direction fluctuations are orthogonal nuisance directions. This first-order statement does **not** establish equality of higher-order bias or coverage. No higher-order repeated-sampling ranking is claimed here.

The calculations above compare responses to fixed observed data, especially weak data. They do not prove that either resulting credible interval has nominal frequentist coverage. The measures use the same observations and Gram quantities. Their difference is how the chosen compatible-state measure weights radial uncertainty, not additional measured information supplied by one method.

## 7. Independent numerical check

`check_small_data.py` integrates the original radial kernel directly, without the error-function radial formula or its Taylor expansion. It uses Gauss-Legendre quadrature over radius `[0,14]` and log ratio `[-32,32]` with `160` and `512` nodes. It checks Gram triples `(4,1,1)`, `(1,1,0.7)`, `(1,1,-0.7)`, and `(2,2,0)` at observation scales `0.04`, `0.02`, `0.01`.

For example, at scale `0.01` and Gram triple `(4,1,1)`, the scaled changes `(new-null)/delta^2` are:

| Quantity | Flat numerical / predicted | Tube numerical / predicted |
|---|---|---|
| Log mean | `0.5000041 / 0.5` | `0.9999600 / 1` |
| Log variance | `-0.3490809 / -0.3490659` | `-0.6981654 / -0.6981317` |
| Point readout | `0.5000058 / 0.5` | `0.8333173 / 0.8333333` |

The discrepancies are consistent with the retained fourth-order remainder. The equal-norm orthogonal case is unchanged to approximately `5e-12` after scaling by `delta^(-2)`. Results are saved in `small_data_verification.json`. These checks support the algebra; the proofs above establish the exact identities and the limiting coefficients.
