# One-sided scalar mass confidence benchmark

Independent analytical contribution for the weak-acceleration investigation,
15 September 2026. This scratch review supports contribution 8; it does not
replace the ratio-of-means readout or its probability law.

## Scope and model

A signed scalar reading in a **known physical direction** is observed:

$$
(X,Y)\sim N_2((ma,a),\Sigma),\qquad m>0,\ a>0,
$$

where $\Sigma=\begin{pmatrix}s_F^2&c\\c&s_a^2\end{pmatrix}$ is known and
positive definite. Thus $s_F,s_a>0$, $|c|<s_Fs_a$. The readings can be negative.
This is not an unknown-direction vector model with observed norms substituted
for signed readings.

Write $z=\Phi^{-1}(1-\eta)>0$, with $\eta=.05$ for 95% confidence, and

$$s(m)=\sqrt{s_F^2+m^2s_a^2-2mc},\qquad
T_m=\frac{X-mY}{s(m)}.$$

At the true mass, $T_m\sim N(0,1)$ for every fixed $a$, including arbitrarily
weak acceleration. Therefore the inverted set

$$\boxed{C(x,y)=\{m>0:x-my\le z s(m)\}}$$

has **exact** sampling coverage $P_{m,a}\{m\in C\}=1-\eta$.
This direct proof is a one-sided Fieller-type construction. It does not use a
prior/reference measure, and it does not provide a point estimator.

## Set topology and a lower endpoint

Let $g(m)=x-my-zs(m)$. Then

$$g''(m)=-\frac{z(s_F^2s_a^2-c^2)}{s(m)^3}<0.$$

Hence the rejected set $\{m>0:g(m)>0\}$ is an interval. The accepted set can
be all positive masses, empty, a lower interval, an upper interval, or two
intervals separated by a rejected middle interval. Retain this topology if
claiming exact set coverage. Squaring the defining inequality without retaining
the sign of $x-my$ gives the wrong set.

Define the upward hull by $L=\inf C$ and $[L,\infty)\cap(0,\infty)$, with
$L=+\infty$ and an empty hull when $C=\varnothing$. It contains $C$, and so
its coverage is **at least** $1-\eta$. In generic data:

- If $x<zs_F$, then $L=0$.
- If $x>zs_F$ and $y>-zs_a$, then $L>0$ is the unique positive solution of
  $x-Ly=zs(L)$ and $C=[L,\infty)$.
- If $x>zs_F$ and $y\le-zs_a$, then $C=\varnothing$.

The exact boundary $x=zs_F$ needs separate treatment:

- $y\ge zc/s_F$: $C=(0,\infty)$ and $L=0$;
- $-zs_a<y<zc/s_F$: $C=[L,\infty)$ for a positive root $L$;
- $y\le-zs_a$: $C=\varnothing$.

The boundary has probability zero under the continuous model but matters in
an exact implementation. An empty set is a reported incompatibility with this
one-sided confidence rule; $L=\infty$ is bookkeeping, **not** a claim that the
object has infinite mass. Returning a vacuous lower bound 0 with a separate
empty-set flag is another conservative reporting policy.

For finite positive-definite $\Sigma$, the hull's coverage is strictly greater
than $1-\eta$. Almost surely its noncoverage event is

$$\{L>m\}=\{X>zs_F,\ T_m>z\}.$$

Writing $\lambda_F=ma/s_F$ and
$\rho_T=(s_F^2-mc)/(s_Fs(m))\in(-1,1)$ gives its exact coverage:

$$1-\eta+\int_z^\infty\phi(t)
\Phi\left(\frac{z-\lambda_F-\rho_Tt}{\sqrt{1-\rho_T^2}}\right)dt.$$

The added integral is positive and tends to zero as force signal-to-noise
increases with $\rho_T$ fixed. This is an unconditional sampling guarantee.
Selecting only experiments whose measured force passes a threshold does not
automatically preserve it conditionally on that selection.

## Zero acceleration reading

At $y=0$ with strict $x>zs_F$, the positive lower endpoint is

$$\boxed{L=\frac{c+\sqrt{c^2+s_a^2(x^2/z^2-s_F^2)}}{s_a^2}.}$$

For independent errors this reduces to

$$L=\frac{\sqrt{x^2-z^2s_F^2}}{zs_a}.$$

For strict $x<zs_F$ the lower hull is vacuous ($L=0$). The exact boundary has
$L=2c/s_a^2$ if $x=zs_F,c>0$, and $L=0$ if $c\le0$; the strict formula
above must not silently override this correlated boundary case.

With perfectly known positive force $F$, the scalar residual model reduces to
$Y\sim N(a,s_a^2)$ and

$$C=\{m>0:m(Y+zs_a)\ge F\}.$$

Thus $C=[F/(Y+zs_a),\infty)$ if $Y+zs_a>0$, and $C=\varnothing$ otherwise.
This set and its lower hull have exact $1-\eta$ coverage. At $Y=0$,

$$L=F/(1.644854s_a)\quad\text{for 95% confidence}.$$

## Comparison with the flat-positive-magnitude construction

Under the **same known direction, perfectly known force** information, a flat
reference measure $da$ on $a>0$ gives a lower mass probability bound

$$L_{\rm flat}(y)=F/U_{\rm flat}(y),$$
$$U_{\rm flat}(y)=y+s_a\Phi^{-1}\{1-\eta\Phi(y/s_a)\}.$$

This follows from the positive-truncated Gaussian law for acceleration. At
$y=0$, $U_{\rm flat}=s_a\Phi^{-1}(1-\eta/2)$, and the 95% lower mass bound is
$F/(1.959964s_a)$. It is weaker (smaller) than the exact one-sided confidence
bound at the same observation. This is a comparison of different probability
guarantees, not a mistake in either arithmetic.

In fact this flat-positive lower mass bound has conservative repeated-sampling
coverage throughout this restricted known-force model:

$$U_{\rm flat}(y)>y+zs_a,$$

so $P_a\{L_{\rm flat}(Y)\le F/a\}>1-\eta$. Coverage tends to 1 as
$a/s_a\downarrow0$ and to $1-\eta$ as $a/s_a\to\infty$. This supplies a
useful positive theoretical result: a conditional lower bound from this
construction can have a valid practical confidence interpretation in a clear
restricted domain, while leaving conservatism on the table.

The corresponding uniform guarantee **does not extend automatically** when
force is noisy and both $f,a$ receive flat positive magnitude weighting.
Indeed, for any fixed $a_0>0$ and fixed positive-definite measurement covariance,
let $m\downarrow0$, so $f=ma_0\downarrow0$. Every finite observed pair has a
strictly positive 5th posterior percentile of $f/a$. Under a coupled Gaussian
error realization that percentile converges to a positive limit as $m\to0$,
while the true mass goes to zero. Thus its sampling coverage tends to 0.
Continuity of the posterior percentile follows from the everywhere-positive
smooth Gaussian density restricted to the positive quadrant.

This boundary counterexample is outside a strong-force operating regime. It
does not establish failure within that regime; it establishes why that regime
must be delimited and evaluated, and why a nominal conditional percentage
cannot be advertised as a global sampling guarantee.

## Attribution and limits

- E. C. Fieller (1954), *Some Problems in Interval Estimation*, JRSS B 16(2),
  175–185: https://academic.oup.com/jrsssb/article/16/2/175/7026708 . The publisher
  abstract explicitly describes inversion of standardized linear relations to
  infer ratios. Metadata/abstract inspected; full original PDF not inspected.
- Ulrike von Luxburg and Volker H. Franz (2009), *A Geometric Approach to
  Confidence Sets for Ratios: Fieller's Theorem, Generalizations, and Bootstrap*,
  Statistica Sinica 19(3), 1095–1117: https://arxiv.org/abs/0711.0198 . Abstract
  inspected; it describes reduction to one-dimensional linear-combination
  confidence sets and exact/conservative guarantees.

The precise one-sided topology and formulas above are direct derivations,
not assertions copied from those abstracts. No claim is made that this
confidence construction meets the composition demands imposed on the mass
point readout, or that it supersedes that readout. It benchmarks only how much
one-sided mass information the same scalar observations can support.

The unknown-direction vector model needs a separate benchmark using the same
information as its probability construction. The existing Gaussian residual
chi-square set is one such comparator; substituting observed norms into this
signed scalar formula is not justified.
