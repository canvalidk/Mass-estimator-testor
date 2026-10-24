# Composition bridge: what force additivity does and does not fix

Scratch derivation, 2026-09-15. This is an independent mathematical contribution to the investigation, not a proposed change to the estimator or a managed record.

## 1. Setting and exact thin-tolerance result

Work first in fixed dimensionless coordinates. In the Gaussian example these are `X = F / sigma_F`, `Y = a / sigma_a`; a dimensionless magnitude ratio becomes a physical mass after multiplication by `s = sigma_F / sigma_a`. Write the coordinates below as `F,a`, their magnitudes as `f,alpha`, and their common compatible unit direction as `u`. Thus powers or sums of magnitudes in this note do not silently add quantities with different physical units.

Let `L(F,a)` be a continuous, strictly positive ambient Gaussian density with positive-definite covariance. Consider

\[
B_h(F,a)=h(\alpha)F\times a,\qquad h:(0,\infty)\to(0,\infty),
\]

and condition on

\[
E_\varepsilon=\{F\cdot a>0,\ \|B_h(F,a)\|<\varepsilon\}.
\]

For each fixed nonzero `a`, this mismatch is additive in force. It is equivariant under proper rotations; under reflections it transforms as an axial vector. Its positive-dot zero set is precisely positive collinearity for every positive `h`.

**Result.** Subject to the integrability condition below, as `epsilon -> 0`, the conditional law converges to the law on `(fu,alpha u)` with

\[
 dP_h(f,\alpha,u)
 =\frac{L(fu,\alpha u)h(\alpha)^{-2}}
 {Z_h}\,df\,d\alpha\,d\Omega(u),\qquad
 Z_h=\int L(fu,\alpha u)h(\alpha)^{-2}\,df\,d\alpha\,d\Omega.
\]

In particular, the reference weight relative to flat positive magnitudes and solid angle is

\[
\boxed{\rho_h(f,\alpha)=h(\alpha)^{-2}.}
\]

Here and below all magnitude integrals run over `(0,infinity)`.

**Proof with the normalizing constant.** Fix `a=alpha u` and `F=fv`. The ambient volume is `f^2 alpha^2 df dalpha dOmega(u) dOmega(v)`. If `gamma` is the angle between `u,v`, the allowed directions have `gamma < pi/2` and

\[
\sin\gamma<\delta,\qquad \delta=\frac{\varepsilon}{h(\alpha)f\alpha}.
\]

Their cap area is `2 pi (1-sqrt(1-delta^2))` for `delta<1`, and `2 pi` otherwise. Divided by `pi delta^2`, it tends to one as `delta -> 0` and is at most two for every positive `delta`. Continuity of `L` therefore gives, pointwise at positive `f,alpha`, the limiting angular integral `pi epsilon^2 L(fu,alpha u)/(h(alpha)^2 f^2 alpha^2)`. Its radial Jacobian cancels the denominator. The same argument applies after multiplication by any bounded continuous test function.

A sufficient dominated-convergence hypothesis is that, for an envelope `L(F,a) <= C exp[-c(f^2+alpha^2)]`,

\[
\int_0^\infty h(\alpha)^{-2}e^{-c\alpha^2}\,d\alpha<\infty.
\]

The cap bound above supplies the domination. This proves both the conditional limit and

\[
\Pr(E_\varepsilon)=\pi\varepsilon^2 Z_h+o(\varepsilon^2).
\]

This sufficient condition covers all the polynomial examples below. It is not claimed to be a necessary tail condition for every possible `h` and Gaussian covariance.

An equivalent local geometric check is that the positive-collinearity surface element is `(f^2+alpha^2) df dalpha dOmega`, while the two-dimensional normal disk has area `pi epsilon^2/[h(alpha)^2(f^2+alpha^2)]`. These factors cancel in exactly the same way.

## 2. Axes, propriety, and uncertainty tails

The preceding formula concerns the interior `f,alpha>0`. The axes have zero ambient probability at every finite tolerance, but neighborhoods of them can dominate the limit if the putative surface weight is not integrable.

* Near the acceleration axis, propriety requires `integral_0^delta h(alpha)^(-2) dalpha < infinity`. This condition is necessary because a nondegenerate Gaussian is bounded below on any fixed compact set of positive `f` and small `alpha`.
* The force axis adds no singular weight: its magnitude measure remains `df`.
* If `h(alpha) ~ c alpha^beta` as `alpha -> 0`, the axis condition is exactly `beta < 1/2`.
* Polynomial growth or decay of `h` at infinity is harmless under a nondegenerate Gaussian. Faster behavior requires checking the actual integral. If the displayed reference law is improper, it must not be normalized formally; a boundary-dominated tolerance limit needs separate analysis.

Continuity of the residual at `a=0` is a different condition. Extending `B_h(F,0)=0` is continuous for every fixed `F` exactly when `alpha h(alpha) -> 0`. For a power behavior this is `beta>-1`. Thus many continuous, proper alternatives exist; continuity plus propriety allows `-1<beta<1/2` for this local power family.

For tame Gaussian tails, define `H(f,alpha)=integral L(fu,alpha u) dOmega`. The ratio density is

\[
p_R(r)=Z_h^{-1}\int_0^\infty \alpha h(\alpha)^{-2}H(r\alpha,\alpha)\,d\alpha.
\]

If `h(alpha) ~ c alpha^beta` at the acceleration axis, then

\[
p_R(r)\sim C_\infty r^{2\beta-2}\quad(r\to\infty),
\qquad p_R(r)\to C_0>0\quad(r\downarrow0).
\]

Consequently its power moments are finite exactly for `-1<t<1-2 beta`, under the stated power and tail assumptions. In particular, every smooth `h` with `0<h(0)<infinity` retains the flat model's `1/r^2` upper tail, even when its middle quantiles and point readout change substantially.

## 3. Explicit alternatives, including a smooth counterexample

At standardized isotropic measured zero, the likelihood on the compatibility surface is proportional to `exp[-(f^2+alpha^2)/2]`. Let `c0=sqrt(2/pi)` denote the half-normal mean. The following are exact, proper laws.

| Mismatch multiplier | Induced weight | Null magnitude laws | Ratio-of-means `E[f]/E[alpha]` |
|---|---|---|---:|
| `h=1` | `1` | independent `chi_1, chi_1` | `1` |
| `h=1/alpha` | `alpha^2` | independent `chi_1, chi_3` | `1/2` |
| `h=(1+alpha^2)^(-1/2)` | `1+alpha^2` | `f~chi_1`; `alpha` is the equal mixture of `chi_1,chi_3` | `2/3` |

For the second row, `B_h=F cross (a/||a||)` is a force-direction mismatch on nonzero `a`. It is force-additive but has no continuous extension at `a=0`: approach direction matters. Its null ratio density and ratio mean are

\[
p_R(r)=\frac{4}{\pi(1+r^2)^2},\qquad E[R]=\frac2\pi.
\]

The point `1/2` is therefore still not the mean of the ordinary ratio law.

The **third row avoids that continuity objection**. The function `h(||a||)=(1+||a||^2)^(-1/2)` and the whole residual are smooth at every Cartesian `a`, including zero. Since the half-normal second and third moments are `1` and `2 c0`, respectively,

\[
E[f]=c_0,\qquad E[\alpha]=\frac{c_0+2c_0}{2}=\frac32c_0.
\]

Its ratio density is

\[
p_R(r)=\frac1{\pi(1+r^2)}+\frac2{\pi(1+r^2)^2}
=\frac{r^2+3}{\pi(1+r^2)^2}.
\]

It retains an infinite ratio mean, but its ratio-of-means is `2/3`. In physical units the three point readouts are `s`, `s/2`, and `2s/3`. A dimensional form of the smooth multiplier uses a declared acceleration scale, `h(alpha)=[1+(alpha/a0)^2]^(-1/2)`; the table chooses that scale as one in the fixed standardized coordinates.

There is also a short direct proof of this example's tolerance limit including the axes. Write `a=alpha u` and `F=b u+w`, with `w` in the two-dimensional perpendicular plane. Positive dot product is exactly `b>0`; the volume is `alpha^2 db d^2w dalpha dOmega`. The tolerance disk is `||w|| < epsilon sqrt(1+alpha^2)/alpha`. For each positive `alpha`, its integral divided by `pi epsilon^2`, after multiplication by the volume factor `alpha^2`, tends to `(1+alpha^2)L(bu,alpha u)`. A Gaussian envelope bounds the same expression at every `epsilon` by `C(1+alpha^2)exp[-c(b^2+alpha^2)]`, because its integral is at most disk area times the envelope. This bound is integrable on `b,alpha>0`, so dominated convergence proves the full limit even where the disk radius becomes large near `alpha=0`.

More generally, every `h(alpha)=(1+alpha^2)^q`, finite real `q`, is smooth, positive, and yields a proper Gaussian conditional law with weight `(1+alpha^2)^(-2q)`.

**Conclusion:** force additivity, rotation equivariance, smoothness at the axes, propriety, and the same compatibility zero set do not jointly determine the flat reference measure. The smooth third row is an explicit counterexample. It does not establish that any alternative should be preferred.

### Broader one-sided class

The cross-product family is a subclass of the continuous, force-additive, proper-rotation-equivariant vector residuals that vanish exactly when force is parallel to nonzero acceleration. With `u=a/alpha`, the general transverse form is

\[
B(F,a)=p(\alpha)(I-uu^T)F+q(\alpha)F\times u,
\quad p(\alpha)^2+q(\alpha)^2>0.
\]

The two terms are orthogonal and have the same transverse-force norm. Fixed norm tolerance therefore induces

\[
\rho(f,\alpha)=\frac{\alpha^2}{p(\alpha)^2+q(\alpha)^2}.
\]

The restricted family has `p=0, q=alpha h(alpha)`. If the residual is required to transform axially under the full orthogonal group, rather than merely equivariantly under proper rotations, the polar projection term is excluded. This still leaves the radial `h` freedom.

## 4. What exchange symmetry can force

Specify exchange in the fixed standardized coordinates, so it maps `(f,alpha,u)` to `(alpha,f,u)`. Require the **whole conditional law** to be equivariant: applying the same `h` construction to the exchanged likelihood must equal the pushforward of the original conditional law.

The two positive densities then have the common exchanged likelihood, multiplied respectively by `h(alpha)^(-2)` and `h(f)^(-2)`. Equality after normalization implies

\[
h(\alpha)^{-2}=C h(f)^{-2}\quad\text{for all }f,\alpha>0.
\]

Taking `f=alpha` gives `C=1`; hence `h` is constant. A constant only rescales tolerance and leaves the normalized law unchanged. It is enough to impose symmetry of the full posterior at one everywhere-positive exchange-symmetric likelihood, such as isotropic measured zero.

This conclusion depends on the one-sided `h(alpha)` family and a fixed, state-independent tolerance. Exchange symmetry of more general reference weights admits every symmetric weight `rho(f,alpha)=rho(alpha,f)`; it does not force flatness. In the broader `p,q` family above, it forces `p(alpha)^2+q(alpha)^2` proportional to `alpha^2`, which fixes the induced weight without uniquely fixing the residual.

Also distinguish whole-law exchange symmetry from merely declaring the reciprocal readout to be the inverse, or observing a reciprocal-symmetric **mass marginal**. Those weaker statements do not prove the result. For an explicit marginal counterexample, choose any `c>1` and set

\[
w(\alpha)=\tfrac12\left[c^{-1}e^{(1-c^{-2})\alpha^2/2}
 +c e^{(1-c^2)\alpha^2/2}\right],\qquad h=w^{-1/2}.
\]

At standardized measured zero, `f~chi_1` and `alpha` is the equal mixture of `c chi_1` and `c^(-1) chi_1`. The ratio is an equal mixture of half-Cauchys with reciprocal scales, so `p_R(r)=r^(-2)p_R(1/r)`. Nevertheless the joint law is not exchange-symmetric, `h` is nonconstant, and the ratio-of-means is `2/(c+c^(-1))`, below one. This example is proper for the standardized unit-covariance Gaussian model; no claim of propriety for every unrelated covariance is needed.

## 5. Radial transformations and state-dependent tolerance

### A global monotone radial transformation preserves the limit

Let `A=||F cross a||`, and replace the residual's norm by `g(A)`, where `g` is continuous and nondecreasing, `g(0)=0`, and `g(t)>0` for every `t>0`. Then `g(A)<epsilon` is equivalent, up to a zero-probability boundary, to `A<t_epsilon` for a scalar cutoff tending to zero. Thus this is exactly a reparameterization of the sequence of tolerance sets. Whenever the area-conditioned limit exists, it is unchanged.

No positive derivative at zero is needed: `g(t)=t^p` works for every `p>0`. The same statement applies to a global monotone radial transformation of an already weighted `B_h`: it preserves that `h`-weighted limit. Such nonlinear transformations are generally no longer force-additive. A nontrivial zero plateau would change the zero set and is excluded by `g(t)>0` for positive `t`.

### State dependence generally changes the limit

Instead condition on `A < epsilon tau(F,a)`, with a positive continuous tolerance factor whose trace on the compatible surface is `tau_0(f,alpha,u)`. Under corresponding integrability and domination assumptions, the normal disk area gains a factor `tau_0^2`, so

\[
\rho_{\tau}(f,\alpha,u)=\tau_0(f,\alpha,u)^2.
\]

For `B_h` with that tolerance the weight becomes `tau_0^2/h(alpha)^2`. Only the trace of the tolerance factor matters; variations off the compatible surface with the same trace do not change the leading limit. A constant trace preserves the flat limit. A compensating trace `tau_0` proportional to `h(alpha)` restores flatness even for a nonconstant residual multiplier.

Examples for the pure area residual in three dimensions:

| Tolerance factor on compatible pairs | Resulting weight | Interpretation |
|---|---|---|
| `1` | `1` | fixed area tolerance |
| `f alpha` | `f^2 alpha^2` | fixed sine-of-angle tolerance |
| `sqrt(f^2+alpha^2)` | `f^2+alpha^2` | first-order standardized normal-distance tolerance |
| `(1+alpha^2)^(-1/2)` | `(1+alpha^2)^(-1)` | a smooth proper alternative with the same exact compatibility relation |

These functions can vanish at boundary axes; the result uses their positive interior traces and still requires the global integrability check. The normal-distance row states its local leading behavior, not equality with exact distance at finite tolerance.

More generally, a state-dependent radial transform with a common leading power `g_s(t) ~ c(s)t^p`, positive `c(s)` and `p>0`, yields weight `c(s)^(-2/p)` relative to the original area-conditioned measure, provided the approximation and integrability justify passage to the limit. A state-independent coefficient is harmless; a state-dependent coefficient is additional weighting information. If powers themselves vary with state, one cannot assume a common surface law without a separate asymptotic analysis.

## 6. Scope of the bridge

Force additivity restricts the form of a residual; it does not by itself specify how a small residual is measured or how its tolerance is standardized. A fixed norm tolerance applied to a specified residual supplies the missing reference weight. Within the one-sided cross-product family, full exchange equivariance of the resulting conditional law removes the radial multiplier freedom. If more general symmetric state-dependent tolerances are admitted, substantial freedom remains. These are logical implications of stated choices, not grounds for silently adopting one choice as the estimator's policy.
