# O3: what supports the flat weighting and the stacked length

**2026-09-30 (late afternoon).** Record, written by Claude after his "Go
ahead" on O3, the item the lemma record marks as carrying the equation since
ours29. The arguments are Claude's; the premises they run on are his rulings
(L8's conditional form, L10 as restated today, L12, L18). It offers a closure
of O3 and asks for the ruling that would make it one.

**What it records.** Two things needed support after the strong form of axis
splitting fell on 29 Sept: the flat-in-vector weighting of candidates (T6),
which is now ours29's whole per-reading factor $e^{H_i/2}$, and the stacked
length (T9), the readout's size of a series. Both are supported without the
strong form. **T6:** the conditional form of L8 that survived the ruling
("components are independent given the reading's shared state") together with
isotropy (L9) forces a Gaussian candidate weight (Herschel–Maxwell), and L10
as restated today (no signal size is declared) sets its scale to infinity: the
flat weight. Additive state (the series' law depends on the readings only
through the sum of their Gram matrices; T12's exactness) is then a theorem,
not a premise, and the exponential is the only kernel with that property.
**T9:** one joining rule (L12) applied to the components as atoms (L18)
forces the Euclidean join, because a vector's length already joins its
components that way. The morning's candidate ("additive state as a property
of mass") is the weaker route and is not needed.

**Relations.** Answers O3 in `mass_estimator_lemmas_2026-09-29.md` (the
rewritten L8, T6, T9, T11, T12, L10, L12, L18) and the 25 Sept record
`../axis_splitting_and_the_24_september_change_2026-09-25.md`, whose premise
was withdrawn as forcing on 29 Sept and is replaced here. Corrects one
sentence of the lemma record's L8 row ("the conditional form forces
nothing"), offered, not placed. Binds nothing.

---

## 1. The question

With the size law gone (ruling 2 of 30 Sept), each reading's factor in
ours29 is $e^{H_i/2}$, the profile of the Gaussian likelihood over a
candidate weighted **flat as a vector** in noise units. The alternative the
project used until 24 Sept weighted the candidate flat in its two
*magnitudes* (the ours21 kernel $K_d(h)$, different in each dimension). The
two agree in 2D and disagree in 1D and 3D. Under F2 they are rival
definitions of mass at weak signal. The 25 Sept record chose the vector
weighting by the strong form of axis splitting; that form is false by ruling.
So: what forces the vector weighting now?

## 2. T6: the conditional form of L8 does force the kernel

**The premises.**

- **L8, conditional form** (ruled 29 Sept, 21:22, and kept when the strong
  form was struck): a reading's components are independent *given the
  reading's shared state* ($N_i$, its off-line misfit $Q_i$, and the size
  scale $\omega_i$ while one existed).
- **L9** (strong): the candidate weight depends only on the candidate's
  length in noise units within the declared subspace.
- **L10** (restated 30 Sept): no signal size is declared, per reading or per
  series.

**The argument.** Take a reading with $k_i\ge2$ and fix its shared state.
The candidate weight is a measure on $\mathbf v\in\mathbb R^{k_i}$ (the
signal in noise units). By L8's conditional form it is a product over
components; by L9 it is a function of $\lVert\mathbf v\rVert$ alone. A
density on $\mathbb R^k$, $k\ge2$, that is both a product of one-dimensional
factors and rotation-invariant is Gaussian, $\exp(-\lVert\mathbf v\rVert^2/2\omega^2)$
(Herschel–Maxwell). So the candidate weight is Gaussian with some scale
$\omega$, or its $\omega\to\infty$ limit, flat.

A finite $\omega$ is a statement that the signal is typically $\omega$ noise
units. That is not a law on $\omega$, but it is a declared size, and L10's
reason (the definition cannot rely on what the setups are like) excludes it
exactly as it excludes a law: a fixed $\omega$ is the claim that every setup
delivers signals of that size. So $\omega=\infty$: the weight is flat in the
vector, and the per-reading factor is
$$
E_i(m)=\int e^{-\frac12\lVert\text{whitened residual}\rVert^2}\,d\mathbf v
\;\propto\;e^{H_i(m)/2},
$$
the $m$-free constant dropped (T7). **This is T6 with the endpoint selected,
from L8's conditional form, L9 and L10.** The lemma record's L8 row says
the conditional form "forces nothing"; it forces the family, and L10 now
picks the member. That sentence is corrected below.

**1D.** Herschel–Maxwell needs $k\ge2$. A 1D reading is a reading declared
$k=1$, and T25 (a $d$-dimensional reading declared $k=1$ equals the 1D
reading exactly) says the 1D kernel is the $k\ge2$ kernel restricted to one
component: $e^{h^2/2}$. The ours21 kernel in 1D,
$1+\sqrt{\pi/2}\,h\,e^{h^2/2}\operatorname{erf}(h/\sqrt2)$, would break T25.

**The improper members.** The 25 Sept record (§5) already found this
route, "not yet discussed with the user", and listed its limits. One of them
matters here: among *improper* weights, $e^{+c\lVert\mathbf v\rVert^2}$ with
$c>0$ is also isotropic and factorises, and Herschel–Maxwell does not exclude
it. L10 does not either (it prefers ever-larger signals rather than a size).
What excludes it is the coefficient: against the Gaussian likelihood the
weight $e^{c\lVert\mathbf v\rVert^2}$ gives the factor $e^{H/(2(1-2c))}$
(and nothing at all for $c\ge\tfrac12$), and at strong signal the law's
width must be the Gaussian propagation of the declared $\Sigma_i$ (L7,
L13), which is $e^{H/2}$ exactly. So $c=0$. The same fact fixes the
$\tfrac12$ against any other rescaling. (With $\eta$ present a coefficient
above $\tfrac12$ would even be corrected back at strong signal, and one
below it would not; the cap makes $\tfrac12$ the largest coefficient that
leaves strong signal alone. Not needed, since L7 already fixes it.)

**What is new against the 25 Sept record.** Its second route rested on
"independent components", which it called "a modelling claim, not a law",
and left the scale $\omega$ open (it leaned toward a learned one). Since
then: independence given the reading's shared state became his ruling (L8,
29 Sept), the scale question was closed by L10 (30 Sept), and the improper
members are removed by L7 above. Nothing in the route itself is new; its
premises are.

## 3. Additive state is a theorem, and the exponential is the only kernel with it

The morning's candidate was to *assume* that a series' law depends on the
readings only through the sum of their Gram matrices $S=\sum_iM_i$ and derive
the kernel from that. The direction is right and the assumption is not needed:
under §2 the log-law is $\tfrac12\sum_i\mathbf u^\top M_i\mathbf u=\tfrac12\mathbf u^\top S\mathbf u$,
so additive state follows. For the record, the converse holds too, which is
what makes O3 a real fork rather than a matter of taste:

**Claim.** Let the per-reading log factor be $f(h_i^2)$ with
$h_i^2=\mathbf u(\theta)^\top M_i\mathbf u(\theta)$. The series' law depends
on the readings only through $S$ iff $f$ is affine, i.e. the kernel is
$e^{bh^2}$.

**Proof.** Two 1D readings $\lambda(1,0)$ and $\lambda(0,1)$ have the same
Gram sum, $\lambda^2I$, as $\lambda(1,1)/\sqrt2$ and $\lambda(1,-1)/\sqrt2$.
Equal laws for every $\theta$ means
$f(\lambda^2a)+f(\lambda^2(1-a))$ is the same for all $a\in[0,1]$, at every
$\lambda$; writing $x=\lambda^2a$, $f(x)+f(\Lambda-x)$ is constant in $x$ on
$[0,\Lambda]$ for every $\Lambda>0$, so $f'(x)=f'(\Lambda-x)$ for all
$\Lambda>x$, and $f'$ is constant. The converse is immediate. $\square$

Checked: with the two configurations above, the exponential kernel gives
identical laws (difference $0$ to rounding); the ours21 kernel differs by up
to $0.055$ nats at $\lambda=1$ and $0.52$ nats at $\lambda=2.5$.

So "readings combine by adding their Gram matrices" and "the kernel is the
exponential" are the same statement. Which one is the premise is a matter of
which is a property of mass. §2 says neither needs to be: both follow from
L8's conditional form, L9 and L10.

## 4. T9: the stacked length from one joining rule

The readout divides the mean size of the series' candidate force by the mean
size of its candidate acceleration, and the question was how a series'
*size* is formed from its readings' sizes. ours25b answered
$\big(\sum_i\lVert\mathbf a^\ast_i\rVert^2\big)^{1/2}$, the length of the
stacked vector, and supported it by the strong form of axis splitting.

**The premises.** L12 (medium): one joining rule, $\operatorname{size}(A\cup B)=\operatorname{size}(A)\oplus\operatorname{size}(B)$
with one fixed $\oplus$, a series having one size per channel per
candidate. L18 (ruled 28 Sept): the reading structure, which components
belong to which reading, is a declared instrument fact; the atoms of a
series are its components.

**The argument.** Within one reading, the candidate's size is the length of
a vector, and a vector's length is the Euclidean join of its components:
$\lVert\mathbf a^\ast_i\rVert=\big(\sum_c a^{\ast2}_{i,c}\big)^{1/2}$. That is
geometry, not a choice. L12 says there is one $\oplus$ for every join of
sets of components. The within-reading join is Euclidean, so the
across-reading join is Euclidean, and the series' size is
$\big(\sum_i\sum_c a^{\ast2}_{i,c}\big)^{1/2}$: the stacked length. The plain
sum ($\sum_i\lVert\mathbf a^\ast_i\rVert$) would need a second $\oplus$ for
joins across readings, which L12 forbids.

This does not reinstate the strong form of L8. Strong axis splitting said a
3D reading and its three 1D readings give the same *law*; that is false
because the declaration differs (the misfit weight and $\eta$ read it). The
stacked length is a statement about the *size* of a fixed set of
components, and it is the same whether those components are declared as one
reading or three. The rewritten L8 lets the law differ; L12 keeps the size
the same.

## 4b. O8 dissolves under L10 and T7 (added the same afternoon)

**The worry** (his question, 15:56, and Claude's flag before it). L9 says
the candidate weight depends only on the candidate's length **in noise
units**. With an instrument whose uncertainty differs by axis (or is
correlated, a full $\Sigma_i$), noise units and physical space disagree on
which directions are "the same size", and O8 asked which is right.

**Why it no longer arises.** That fork is real for a Gaussian weight of
finite width: a round Gaussian in noise units is squashed in physical space.
But after L10 the weight is flat, and a flat weight is flat under any linear
change of coordinates. The two choices differ only by the Jacobian of the
whitening, restricted to the candidate subspace at $m$:

- flat in **noise units** (the whitened candidate): the Jacobian cancels
  and the factor is $e^{H_i/2}$ exactly;
- flat in **physical** $\mathbf a^\ast$: the factor is
  $e^{H_i/2}\,\det G_i(m)^{-1/2}$, $G_i(m)=B_m^\top\Sigma_i^{-1}B_m$, and
  that determinant depends on $m$. It is a factor on mass attached to every
  reading, a prior counted $N$ times, which T7 forbids.

So T7 picks noise units by itself, whatever $\Sigma_i$ is. The flat
weighting needs no isotropy premise beyond what fixes the Gaussian family
(and at $\omega=\infty$ the family's shape no longer matters).

**Checked by direct quadrature** (2D, a covariance with SDs spanning a factor
20 and correlations across axes and between channels; both integrals done on
grids, not from closed forms): $\log E_{\rm noise}-H/2$ is constant in $m$
to $10^{-8}$ ($-9.19260918$ at $m=0.1$, 0.5, 1, 3, 10, 100), and
$\log E_{\rm phys}-\log E_{\rm noise}$ equals $-\tfrac12\log\det G(m)$ to
six decimals at every $m$, running from $-1.04$ to $-10.88$.

**Two remarks.** The same $m$-dependent determinant appears with isotropic
noise too ($\det G=(m^2/\sigma_F^2+1/\sigma_a^2)^k$): it is the tester's
radius-versus-acceleration nuisance split, so noise units were always forced
by T7 for the flat weight, and O8 was only ever a question about finite-width
weights. With the size law gone, there are none. **O8 closes under L10**
(offered with the rest, §6).

## 5. What the closure rests on, and what it does not do

| piece | rests on | grade |
|---|---|---|
| T6, the Gaussian family | L8 conditional form + L9 (Herschel–Maxwell) | L8 ruled, L9 strong |
| T6, the flat endpoint | L10 as restated (no declared size, fixed or lawful); the improper $c>0$ members removed by L7 | ruled 30 Sept; L7 mild |
| T6 in 1D | T25 | theorem |
| T6 under any $\Sigma_i$ (O8) | L10 + T7: flat is coordinate-free up to a Jacobian, and T7 removes the $m$-dependent one | theorem (§4b) |
| the coefficient $\tfrac12$ | L7, L13 | mild, forced |
| additive state, T11, T12 | corollary of T6 | theorem |
| T9, the stacked length | L12 + L18 (components as atoms) | L12 medium, L18 ruled |

- It does **not** say that L8's conditional form is itself forced. It is his
  ruling of 29 Sept, and the closure inherits its standing. If "components
  independent given the shared state" were ever doubted, O3 reopens at T6.
- It does **not** fix the misfit weight's middle (O12) or $\eta$'s finite-$N$
  form (O14); those are outside the per-reading kernel.
- It does **not** use the additive-state candidate of the morning as a
  premise; §3 records it as a theorem and a characterisation.
- The weakest link is **L12 at "medium"**. If L12 is doubted, T9 falls back
  to being chosen, and the plain-sum readout (ours24) returns as a rival.

## 6. Correction offered to the lemma record (not placed)

- **L8 row:** "the conditional form forces nothing" → "the conditional form,
  with L9, forces the Gaussian candidate family (Herschel–Maxwell); with L10
  it forces the flat endpoint (this record §2)".
- **T6 row:** source L8 (conditional), L9, L10; "support from L8 withdrawn"
  → "supported by L8's conditional form; the strong form was never needed".
- **T9 row:** source L12, L18; "needs its own (O3)" → "forced by one joining
  rule applied to components as atoms (this record §4)".
- **T11, T12:** additive state now a corollary of T6 (§3).
- **O3:** close, on his ruling; note the two dependencies (L8 conditional
  form, L12 at medium).
- **O8:** close (§4b): under L10 the weight is flat, and T7 fixes noise units
  for any covariance; the question survives only for finite-width weights,
  of which ours29 has none.

## 7. Checks

Sandbox, numpy and scipy. §4b: 2D, $\Sigma$ with SDs (2, 0.1, 0.5, 1.5) and
a random correlation, $\mathbf y=(1.3,-0.4,0.7,0.2)$; the noise-unit
integral on a 1601² grid over the whitened candidate plane, the physical one
on an 801² grid placed on the posterior's axes (placement only; the integrand
is the Gaussian likelihood itself), $m\in\{0.1,0.5,1,3,10,100\}$. §3: seven candidate angles $\theta\in[0.05,\pi/2-0.05]$;
laws $\sum_i\log K(\lvert\mathbf y_i\cdot\mathbf u(\theta)\rvert)$ for the
two configurations and two kernels, centred over $\theta$; max absolute
difference reported. Script not kept.

## Sources

- `../mass_estimator_lemmas_2026-09-29.md`: L7, L8, L9, L10, L12, L13, L18,
  T6, T7, T9, T11, T12, T25, O3.
- `../axis_splitting_and_the_24_september_change_2026-09-25.md` (the premise
  replaced).
- `../excitation_weighting_premise_2026-09-24.md` §5 (the Herschel–Maxwell
  second route and its four limits; this record takes that route with the
  conditional L8 and answers the limits one by one).
- `readings_are_evidence_just_less_2026-09-29.md` §10 (the ruling that struck
  the strong form).
- `ours29_misfit_weight_and_bartlett_power_2026-09-30.md` §§3, 5, 11.
- `../ours25b_the_many_reading_readout_2026-09-25.md` (the stacked length).
