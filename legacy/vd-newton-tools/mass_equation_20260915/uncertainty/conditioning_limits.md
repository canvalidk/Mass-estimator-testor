# Two ways to approach exact positive collinearity

Bounded mathematical check, 2026-09-15. Scratch contribution only.

Let X and Y be independent standard Gaussian d-vectors. In physical units the
latent vectors are sigma_F X and sigma_a Y. Write x=||X||, y=||Y||, and
s=sigma_F/sigma_a. The exactly compatible manifold has parameterization

    (X,Y) = (x u,y u), x,y>0, u in S^(d-1).

Two different neighborhoods of this probability-zero manifold produce
different limiting conditional laws. This is a statement about conditioning
and reference measures; neither neighborhood has been selected as physically
required by the mass-readout theorem.

## Angular tolerance

Condition the relative angle of X and Y to be at most epsilon. Under the
independent isotropic Gaussian law, directions and radii are independent.
Consequently the magnitude law is iid chi_d for **every positive angular
tolerance**, not just in its zero limit. With R=M/s=x/y,

\[
p_R(r)={2\over B(d/2,d/2)}{r^{d-1}\over(1+r^2)^d},\qquad r>0.
\]

The limiting compatible measure in magnitude/direction coordinates is
proportional to `(xy)^(d-1) dx dy dOmega`; physical scale constants cancel.
The point is `E[f]/E[alpha]=s`, and

\[
E[R^t]={\Gamma((d+t)/2)\Gamma((d-t)/2)\over\Gamma(d/2)^2},
\quad |t|<d,
\qquad
\operatorname{Var}\log R={1\over2}\psi_1(d/2).
\]

For d=3:

\[
p_R(r)={16\over\pi}{r^2\over(1+r^2)^3},\qquad
F_R(r)={2\arctan r\over\pi}-{\sin(4\arctan r)\over2\pi}.
\]

- Point: s.
- Central 95% mass interval / s: [0.2545000379, 3.9292724999].
- SD(log(M/s)): sqrt(pi^2/4-2) = 0.6836673901.
- E[M]/s = 4/pi = 1.2732395447.
- E[M^2]/s^2 = 3; Var(M)/s^2 = 3-16/pi^2 = 1.3788610617.
- The third moment diverges logarithmically. The mean and second inverse-mass
  moments follow by reciprocity; the moment domain is -3<t<3.
- E[f]/sigma_F = E[alpha]/sigma_a = 2 sqrt(2/pi) = 1.5957691216.

## Standardized Euclidean tube

Instead condition the standardized pair (X,Y) in R^(2d) to lie within Euclidean
distance epsilon of the closed positive-collinearity cone. The metric matters:
this is Euclidean distance **after** dividing channels by their noise scales.

The compatible manifold has metric

\[
ds^2=dx^2+dy^2+(x^2+y^2)\,ds_{S^{d-1}}^2.
\]

Its induced surface volume is therefore
`(x^2+y^2)^((d-1)/2) dx dy dOmega`. Multiplying by the ambient isotropic
Gaussian density gives the thin-tube conditional limit. Boundary pieces have
one higher codimension and disappear from the normalized limit.

Set x=r sin(theta), y=r cos(theta). The radial law is chi_(d+1) and theta is
uniform on (0,pi/2), independently of r. Thus, for every d,

\[
p_R(r)={2\over\pi(1+r^2)},\qquad
q_R(p)=\tan(\pi p/2),\qquad
\operatorname{SD}\log R=\pi/2.
\]

The ratio-of-means point is still s. All positive mass moments of order at
least one diverge, as do reciprocal moments of order at least one.

For d=3:

- Point: s.
- Central 95% mass interval / s: [0.0392901070, 25.4516995794].
- SD(log(M/s)): 1.5707963268.
- E[f]/sigma_F = E[alpha]/sigma_a = 3/sqrt(2*pi) = 1.1968268412.

The mass ratio law equals the existing flat-positive-magnitude null law, but
the full compatible-state law does not: its radial law is chi_4 rather than
chi_2. In d=3 its normalized magnitude density is

\[
g(x,y)={(x^2+y^2)e^{-(x^2+y^2)/2}\over\pi}
=\tfrac12[p_{\chi_3}(x)p_{\chi_1}(y)+p_{\chi_1}(x)p_{\chi_3}(y)].
\]

This is useful for reuse: agreeing on both the point and the complete mass
uncertainty law still does not establish agreement on the full joint latent
force-acceleration law.

## Independent finite-tube verification in d=3

The verification integrates the original iid chi_3 magnitudes and the uniform
relative direction cosine. It does not use the limiting surface measure.

Let a=x^2, b=y^2, and c=X dot Y. Distance squared to the positive cone is

\[
\begin{cases}
\frac{a+b-\sqrt{(a-b)^2+4c^2}}2,&c\ge0,\\
\min(a,b),&c<0.
\end{cases}
\]

Hence at fixed x,y the exact tube acceptance probability is 1 if either norm
is at most epsilon; otherwise it is

\[
\frac12\left[1-\sqrt{(1-\epsilon^2/x^2)(1-\epsilon^2/y^2)}\right].
\]

For fixed positive x,y this is asymptotic to
`epsilon^2 (x^2+y^2)/(4x^2y^2)`. Multiplying by the iid chi_3 radial density
gives exactly the claimed tube-limit surface weight. The total acceptance
probability divided by epsilon^2 tends to 1/2.

Selected direct-quadrature results for the **raw norm ratio**:

| epsilon | Acceptance / epsilon^2 | SD(log ratio) | Lower 2.5% mass quantile / s |
|---:|---:|---:|---:|
| 0.1 | 0.5240213 | 1.5136830 | 0.04593009 |
| 0.02 | 0.5052186 | 1.5830989 | 0.03317209 |
| 0.005 | 0.5013235 | 1.5826004 | 0.03738768 |
| 0.001 | 0.5002657 | 1.5760114 | 0.03889583 |
| 0.0005 | 0.5001329 | 1.5741541 | 0.03909216 |
| Limit | 0.5 | 1.5707963 | 0.03929011 |

The reciprocal upper quantile follows by exchange symmetry. The point stays
1 in standardized units throughout. Spread can overshoot its limiting value;
convergence is not monotone. The finite-tube mean is finite for d=3 but grows
as the tube shrinks, consistent with the divergent limiting mean. The raw
means at epsilon .001 and .0005 are 6.1457592 and 6.5885458.

**Projection caveat:** projecting each accepted pair to its nearest point of
the closed positive cone creates boundary atoms at every positive tube radius.
For negative dot product the best projection sets one magnitude to zero.
Thus projected mass can have atoms at zero and infinity, and a finite-tube
projected log standard deviation is not defined. These boundary probabilities
vanish in the thin-tube limit. The table uses raw norms, whose ratio approaches
the limit without introducing artificial endpoint atoms.

`tube_limits.py` uses a rationalized angular probability and transformed
quadrature to remove the square-root cusp at x=epsilon or y=epsilon. Moments
were checked at 256 and 512 nodes per interval; lower quantiles at epsilon
.005, .001, and .0005 agree between 128 and 256 nodes to the bisection
resolution (about 2.4e-10). Gaussian radial truncation is at 14. These are
numerical checks for the listed cases, not a universal convergence certificate.
Full output is in `tube_verification.json`.
