# Mass as force paired with the reciprocal of acceleration: units and axes

**2026-10-10.** Record, written by Claude at the user's request ("make a
record on [the idea]. But don't call it conjecture 1. Give an explanation of
how we arrived here analyzing units and axes"). The idea and every premise
marked **his** are the user's, quoted verbatim where possible. The
derivations, checks and the literature in §2 are Claude's. It covers the
**aligned case only**, on his instruction: "For now, we don't need to concern
ourselves too much with the unaligned cases. What we are doing is
understanding the basic form for now, and hopefully the natural extension
will drop out."

**Relations.** It scopes and binds nothing. It bears on F1 (the estimator is
the definition) and on the exact law of `mass_estimator_lemmas_2026-09-29.md`
(L1, T1, T2): it restates $m=|\mathbf F|/|\mathbf a|$ on aligned pairs as an
operation that is defined, rather than as a division by a vector. §8 parks
what was said about unaligned pairs and about O11.

**Revised in place 2026-10-10, afternoon, on his "Yes put it in. But keep
L14 provisional."**: §9 added, correcting the comparison with the wave in
§§1–2 and §5 (the mass corner is the forward pairing, in two versions) and
noting, unruled, the L14/swap trade. Sources moved to §10. Nothing else
touched.

**Revised in place 2026-10-10, later afternoon, on his "Do 1 and 2 please"**:
§9.4 added with the start on §9.3's four items and their checks; §9.1's
second point corrected to match. The results are carried into
`mass_estimator_lemmas_2026-10-10.md`, which supersedes the 29 September
lemma record. Nothing else touched.

---

## 0. The idea

**His**, stated as a conjecture: "every time we have a cancelling of
dimensions, instead of division it must instead be caused by the matrix
multiplication by the reciprocal covector."

Applied to Newton II, **his**: "instead of $m = F/a$ or its equivalent. We
would be interested in $F_i\,a'^j$ where $a'$ would be the reciprocal
vector."

In the aligned case this record arrives at

$$
m=\mathbf F\cdot\mathbf a',\qquad \mathbf a'=\frac{\mathbf a}{|\mathbf a|^2},\qquad |\mathbf a'|=\frac{1}{|\mathbf a|},\qquad (\mathbf a')'=\mathbf a,
$$

so mass is force **paired with the reciprocal of acceleration**, and every
unit that cancels on the way cancels by a pairing. The build-up of §6 makes
this concrete:

$$
m=\mathbf F\cdot\mathbf a'=\Delta\mathbf p\,\underbrace{(\Delta t)'\,\Delta t}_{=1}\cdot(\Delta\mathbf v)'=\Delta\mathbf p\cdot(\Delta\mathbf v)'.
$$

## 1. Where it started: the plane wave

**His**: "for a planar wave, if you take the dot product of a wavenumber and
the displacement vector, it doesn't work. The dimensions don't work. You have
to treat the $\mathrm m^{-1}$ vector as a covector. The method by which you
cancel dimensions is that..."

The phase is a pairing, $\varphi=k_i x^i$. The wavevector carries
$\mathrm m^{-1}$, the displacement carries $\mathrm m$, and the units cancel
because each component of $\mathbf k$ carries the inverse unit of the
matching component of $\mathbf x$. No dot product (no metric) is needed for
the pairing itself.

**His question:** "When the wavevector and displacement are aligned, do we
get a scalar equation?" Yes:

$$
\varphi=|\mathbf k|\,|\mathbf x|\quad\Longrightarrow\quad |\mathbf k|=\frac{\varphi}{|\mathbf x|},
$$

a division of magnitudes (wavenumber is phase per metre; rad/m times m is
rad). It is the same shape as $m=|\mathbf F|/|\mathbf a|$ on aligned pairs.

## 2. Is the idea already established?

Checked by Claude at his request ("is conjecture 1 actually consensus in
physics?").

- **One axis: yes, as mathematics.** The careful formalisations of units
  treat each unit as a one-dimensional space, the inverse unit as its dual,
  and products of units as tensor products. Janyška, Modugno and Vitolo: the
  inverse of a unit $u$ "is just the dual semi-basis of $(u)$". Whitney (1968)
  builds physical quantities the same way. There, "divide by 3 m" literally is
  "apply the covector $(3\,\mathrm m)^{-1}$", and a length cancels against an
  inverse length by contraction. In one axis the reciprocal is unique and
  needs no metric, which is why everyday unit algebra ($Q=\{Q\}[Q]$, dividing
  freely) loses nothing and nobody needs to state the principle.
- **Several axes: the covector reading is standard, the universal claim is
  not.** Wavevector, force and momentum ($p_i=\partial L/\partial v^i$) are
  one-forms in geometric physics (Misner, Thorne and Wheeler draw the phase as
  the count of wave surfaces crossed). No source found states "every
  cancellation is a contraction" as a law. The closest is Hart's
  *Multidimensional Analysis*, on vectors and matrices of mixed units: a
  matrix is multipliable only if its array of units is an outer product
  $\mathbf u\,\tilde{\mathbf v}$, a unit vector times the reciprocal of
  another (his reciprocals are taken per component, not geometrically).
- **Where the idea has content.** It is uncontroversial exactly where it is
  trivial (one axis). With $d$ axes a nonzero vector $\mathbf x$ has a whole
  plane of covectors with $\langle\mathbf k,\mathbf x\rangle=1$, so "the
  reciprocal of a vector" needs something to pick one.

## 3. No tensors: orthogonal axes, and the symmetry of duality

**His**: "We won't need tensors. Why? 1, all the dimensions are orthogonal.
We are talking about perpendicular space dimensions, and associated time
dimensions. 2, If B is the covector of some vector A, vector A is the
covector to some vector B. What I mean is, hopefully m and $\mathrm m^{-1}$
are symmetric. And t and $t^{-1}$ are symmetric. So there shouldn't be any
difference between them."

- **(1) holds for Newtonian space.** On Cartesian axes the metric is the
  identity, $g_{ij}=\delta_{ij}$, so a vector and a covector have the same
  components and the upper/lower index distinction disappears (Cartesian
  tensors). The condition is that the axes are **orthonormal**, not merely
  orthogonal: every axis carries the same length unit. The identity metric is
  then hidden in that shared unit, not absent. Newtonian mechanics supplies
  one length unit for all of space.
- **(2) holds.** In finite dimensions the dual of the dual is the original
  space, so neither m nor $\mathrm m^{-1}$ is preferred, and likewise s and
  $\mathrm s^{-1}$. An earlier remark of Claude's ("force is the covector")
  was a convention and was dropped.
- **The reciprocal of a vector**, on orthonormal axes and in the aligned
  sense: $\mathbf a'=\mathbf a/|\mathbf a|^2$, with
  $\langle\mathbf a',\mathbf a\rangle=1$, $|\mathbf a'|=1/|\mathbf a|$ and
  $(\mathbf a')'=\mathbf a$. Taking the reciprocal twice returns the vector.

## 4. Two magnitudes

**His**: "There are 2 types of magnitudes right. ... the 'sum in quadrature
of components' ... Then we also have the magnitude. And the magnitude of the
covector is defined differently correct? Isn't the magnitude of the dot
product wavevector the same as the covector wavevector because of how we
define the magnitude?"

- **Of a vector:** the root-sum-square of components, $\sqrt{\sum_i x_i^2}$.
  (Not "rms", which would divide by the number of components and bring in a
  $\sqrt d$.)
- **Of a covector:** its strength as a measuring device, the most it returns
  per unit length of input:
  $$|\mathbf k|_\ast=\sup_{\mathbf x\neq0}\frac{|\mathbf k(\mathbf x)|}{|\mathbf x|}.$$
  For the wave this is the most phase per metre, which is what a wavenumber
  means.
- **They agree** because of Cauchy–Schwarz: the root-sum-square norm is its
  own dual. That is why the dot-product wavevector and the covector
  wavevector have the same magnitude, and why $|\mathbf a'|=1/|\mathbf a|$
  holds in either reading.
- **Alignment needs no second metric.** With the covector magnitude defined
  this way, $\mathbf k$ and $\mathbf x$ are aligned exactly when $\mathbf x$
  attains the supremum: $\mathbf k(\mathbf x)=|\mathbf k|_\ast|\mathbf x|$.
  The scalar equation of §1 is that condition.
- **Where they separate.** They agree only because the vector norm is the
  root-sum-square on orthonormal axes. In noise units,
  $|\mathbf x|_\Sigma=\sqrt{\mathbf x^{\mathsf T}\Sigma^{-1}\mathbf x}$, the
  covector magnitude becomes $\sqrt{\mathbf k^{\mathsf T}\Sigma\,\mathbf k}$.
  The estimator measures in noise units (L9), so this is a second metric to
  keep in view (§8).

## 5. The basic form

With §§3–4, for aligned $\mathbf F$ and $\mathbf a$:

$$
m=F_i\,a'^i=|\mathbf F|\,|\mathbf a'|=\frac{|\mathbf F|}{|\mathbf a|},\qquad \mathbf a\cdot\mathbf F'=\frac1m .
$$

Swapping the roles gives the reciprocal, as the symmetry of (2) requires. In
the aligned, one-axis case pairing and division give the same number, so the
basic form is right but not yet distinguishable from what was there before.
Its content has to show when the units are built up, which §§6–7 begin.

## 6. Constant rates, one multiplication at a time

**His**: "we can start by examining the $\mathrm m\cdot\mathrm m^{-1}$ case
which we have an example for. And see how the very basic story of dividing
and multiplying these quantities when they are constant e.g. $dx/dt = x/t$,
makes sense this time." And on the units: "In F we have kg, m... in $a'$ we
have $\mathrm m^{-1}$, s, s."

**The wave ($\mathrm m\cdot\mathrm m^{-1}$).** The wavevector is the gradient
of phase, $\mathbf k=\nabla\varphi$: differentiating by position is what
makes the $\mathrm m^{-1}$ slot. For a plane wave the gradient is constant, so
$\varphi=\mathbf k\cdot\Delta\mathbf x$ exactly, and the derivative coincides
with "phase times the reciprocal of the displacement". The reverse has a
familiar meaning: $2\pi\,\mathbf k'=2\pi\,\mathbf k/|\mathbf k|^2$ is the
wavelength vector, so $\lambda=2\pi/k$ is the reciprocal of the wavevector
(metres per radian against radians per metre).

**Constant velocity ($\mathrm s\cdot\mathrm s^{-1}$).** Two operations, kept
apart:

- **building a rate:** $\mathbf v=\Delta\mathbf x\,(\Delta t)'$, with
  $(\Delta t)'=1/\Delta t$. Metres and inverse seconds are different units,
  so nothing cancels; the reciprocal is attached, and $\mathbf v$ is a
  converter from durations to displacements;
- **using the rate:** $\Delta\mathbf x=\mathbf v\,(\Delta t)$, where
  $\mathrm s^{-1}$ meets s and cancels by pairing.

"$dx/dt=x/t$" is the statement that, for a constant rate, the derivative is
attaching the reciprocal of the interval. It needs **increments**: the
reciprocal is of a duration, not of a clock reading, since clock readings
have an arbitrary zero.

**Building $\mathbf a'$.** One derivative at a time, through velocity, not
straight from position: $x/t^2$ is not $a$ (with constant acceleration
$x=\tfrac12at^2$). So

$$
\mathbf a=\Delta\mathbf v\,(\Delta t)'\quad\Rightarrow\quad \mathbf a'=(\Delta\mathbf v)'\,\Delta t,\qquad (\Delta\mathbf v)'=\frac{\Delta\mathbf v}{|\Delta\mathbf v|^2},
$$

in $(\mathrm s/\mathrm m)\cdot\mathrm s=\mathrm s^2/\mathrm m$, with
$(\mathbf a')'=\mathbf a$. In words: the time taken per unit of velocity
gained.

**Force**, as a constant $d\mathbf p/dt$: $\mathbf F=\Delta\mathbf p\,(\Delta t)'$.

**The pairing, one unit at a time:**

$$
m=\mathbf F\cdot\mathbf a'=\Delta\mathbf p\,(\Delta t)'\,\Delta t\cdot(\Delta\mathbf v)'=\Delta\mathbf p\cdot(\Delta\mathbf v)'.
$$

The time factors cancel first. What is left is momentum change paired with
the reciprocal of velocity change: $p=mv$, with mass as the converter from
velocity changes to momentum changes. kg is the one unit never paired.
Check: $\Delta\mathbf p=m\,\Delta\mathbf v$ gives
$\Delta\mathbf p\cdot\Delta\mathbf v/|\Delta\mathbf v|^2=m$.

**A condition that came in quietly:** $\mathbf F$ and $\mathbf a$ must be
rates over the **same interval** $\Delta t$, or the time factors do not
cancel. Free for constant quantities; named here so it is not lost.

## 7. Units and axes

**His**: "There is only 1 dimension that is giving us trouble. It's the space
dimension. But that sounds confusing. How do we differentiate unit dimensions
and vector... what should we call them?" Terminology adopted ("ok great"):

- **units** (ISO: quantity dimensions): M, L, T, i.e. kg, m, s;
- **axes**: the $d$ perpendicular directions a vector has components along.

Stated with these words: **of the three base units, only L comes with
axes.** M and T are scalar units (one axis each). Hence:

- every pairing of a scalar unit with its reciprocal
  ($\mathrm s\cdot\mathrm s^{-1}$, $\mathrm{kg}\cdot\mathrm{kg}^{-1}$) is
  canonical, since a one-axis reciprocal is unique;
- the only pairing that can require a choice is
  $\mathrm L\cdot\mathrm L^{-1}$, since a $d$-axis reciprocal is unique only
  on aligned pairs.

The direction in velocity, acceleration, momentum and force is all inherited
from their single factor of L: every vector in Newton II is a scalar-unit
quantity times one directional L. (A Newtonian fact. In relativity T acquires
an axis too, as in the phase $\mathbf k\cdot\mathbf x-\omega t$; Newtonian
mechanics is the regime with exactly one directional unit.)

**The factorisation, his question:** "can we break it up through commutative
algebra? Can we say that $F\cdot a' = \mathrm{kg}\,(\mathrm m\,\mathrm m^{-1})(\dots)$?"
Yes.

- **Units:** they form a free commutative group, so
  $$[\mathbf F\cdot\mathbf a']=\mathrm{kg}\,(\mathrm m\,\mathrm m^{-1})\,(\mathrm s^{-2}\,\mathrm s^{2})=\mathrm{kg}.$$
- **Values:** write $\mathbf F=\tilde{\mathbf F}\,\mathrm{kg\,m\,s^{-2}}$
  and $\mathbf a'=\tilde{\mathbf a}'\,\mathrm{m^{-1}\,s^{2}}$, with
  $\tilde{\mathbf F},\tilde{\mathbf a}'$ lists of $d$ plain numbers. Scalar
  units have one axis, so their values commute with everything and come out;
  all $d$ components ride on the m factor, so the dot product sits in that
  bracket:
  $$\mathbf F\cdot\mathbf a'=\mathrm{kg}\,\big(\tilde{\mathbf F}\,\mathrm m\cdot\tilde{\mathbf a}'\,\mathrm m^{-1}\big)\,(\mathrm s^{-2}\,\mathrm s^{2})=(\tilde{\mathbf F}\cdot\tilde{\mathbf a}')\,\mathrm{kg}.$$
- **The caveat lives in the $\mathrm m\cdot\mathrm m^{-1}$ bracket.** As
  units it is 1. As values it is a dot product of axis components:
  commutative on orthonormal axes, but it returns 1 only when a vector meets
  its own reciprocal. Here it meets the reciprocal of acceleration's length
  part, and returns the mass in kg. The $\mathrm s^{-2}\mathrm s^{2}$ bracket
  never has this issue.
- **The split of the number between factors is not fixed.** Only the product
  $\tilde{\mathbf F}\,\mathrm{kg\,m\,s^{-2}}$ is defined, so a factor
  $\lambda$ can move from the kg part to the m part with nothing changing.
  "The length part of $\mathbf F$" means its direction and axis components,
  not a separately meaningful length.

## 8. Parked, for the build-up and after

Raised during the session, deliberately not pursued (his: "is it possible
that you are confusing yourself by going ahead"). Claude's, untested unless
stated.

1. **Unaligned pairs: contracting breaks the symmetry.** Contracting with the
   shortest reciprocal gives $\mathbf F\cdot\mathbf a/|\mathbf a|^2$ one way
   and $|\mathbf F|^2/(\mathbf F\cdot\mathbf a)$ the other; their ratio is
   $\cos^2\theta$. These are the m- and μ-regressions of
   `what_forces_the_ratio_of_means_readout` (each passes only one composition
   law), so a literal contraction brings back the Deming asymmetry.
2. **The outer product keeps it.** $M=\mathbf F\,\mathbf a'^{\mathsf T}$ (his
   $F_i a'^j$ with two indices) satisfies $M\mathbf a=\mathbf F$ at every
   angle, $\lVert M\rVert=|\mathbf F|/|\mathbf a|$, and the reverse
   construction $\mathbf a\,\mathbf F'^{\mathsf T}$ is its pseudo-inverse. Its
   trace is the contraction, and $\operatorname{tr}M/\lVert M\rVert=\cos\theta$.
   Newton II with a scalar mass becomes "$M$ has no rotation part".
3. **One shared direction restores the symmetry.** With
   $\mathbf a'=\hat{\mathbf b}/(\mathbf a\cdot\hat{\mathbf b})$ and
   $\mathbf F'=\hat{\mathbf b}/(\mathbf F\cdot\hat{\mathbf b})$,
   $m=(\mathbf F\cdot\hat{\mathbf b})/(\mathbf a\cdot\hat{\mathbf b})$ and
   $\mathbf a\cdot\mathbf F'=1/m$ for every $\hat{\mathbf b}$. The bisector of
   $\hat{\mathbf F}$ and $\hat{\mathbf a}$ gives $|\mathbf F|/|\mathbf a|$,
   with $|\mathbf a'|=\sec(\theta/2)/|\mathbf a|$ (the shortest reciprocal of
   $\mathbf a$ giving that value). Off alignment, sweeping $\hat{\mathbf b}$
   gives every positive mass and the negatives beyond: dividing $\mathbf F$
   by $\mathbf a$ off alignment means choosing a direction, and that choice
   is the estimation. Possibly the estimator's candidate direction
   $\hat{\mathbf u}$ is this $\hat{\mathbf b}$; the bisector depends on the
   metric and moves in noise units.
4. **O11 / L3′ by units (offered, not ruled).** Before any reading the only
   declared objects carrying force and acceleration units are $\sigma_F$ and
   $\sigma_a$, so by L3 any default centre is $c=k\,\sigma_F/\sigma_a$; the
   swap sends $c\to1/c$ and $s\to1/s$, so $k=1$ and $c=s$. Holds for one
   noise pair only: with several instruments the ratios $s_1/s_2$ are
   dimensionless and admit any homogeneous mean (the T20 no-go), and a full
   $\Sigma_F,\Sigma_a$ admit several rotation-invariant scalars.
5. **Noise units are a pairing too.** $\Sigma^{-1}$ is the only declared
   object that cancels the mixed units of a stacked $(\mathbf F,\mathbf a)$
   reading; read this way, L9's "isotropy in noise units" is the reading's
   only means of cancelling units, not a convenience.
6. **Next, his plan:** "once we have the basic form that we are targeting, we
   can build up the dimensions of F and a one multiplication at a time." §6
   is the first pass; the open question it sets is whether every pairing in
   the build is as canonical as the one-axis case. By §7 only the
   $\mathrm L\cdot\mathrm L^{-1}$ pairing can fail to be.

## 9. Revision, afternoon: the right comparison with the wave

**Revised in place 2026-10-10, afternoon, on his "Yes put it in. But keep
L14 provisional."** Nothing above this section was changed. This section
corrects the comparison with the wave in §§1–2 and §5, and records what
follows from the corrected one.

**His correction:** "Why can't force be the converter. $m = F'(a)$. That
works too no? You just compared them wrong." §§1–2 compare the mass corner
with the wave's *inverse* problem (given $\varphi$ and $\mathbf x$, find
$\mathbf k$), which is where "the reciprocal of a vector is unique only on
aligned pairs" came from. The mass corner is the wave's *forward* problem: a
covector acting on a vector, as $\varphi=\mathbf k(\mathbf x)$. Two
converters are available, and the units fix which gives which ($\mathbf F'$
is in 1/N, so force as the converter gives inverse mass):

$$
m=\mathbf a'(\mathbf F)=\frac{\mathbf F\cdot\mathbf a}{|\mathbf a|^2},\qquad
\frac1m=\mathbf F'(\mathbf a)=\frac{\mathbf F\cdot\mathbf a}{|\mathbf F|^2}.
$$

Each is a forward pairing, defined at every angle. They are the two sides of
the symmetry of §3 (2): the swap $\mathbf F\leftrightarrow\mathbf a$ sends one
to the other with $m\leftrightarrow1/m$. On aligned pairs both equal
$|\mathbf F|/|\mathbf a|$, which is §5.

### 9.1 What changes

- **The mass corner is a total operation, in two versions.** "Vector
  division is undefined" is not the point (`what_forces_the_chain_to_the_estimator`
  already called it rhetorical, not load-bearing). The fact is
  **non-uniqueness**: two forward pairings that agree only on aligned pairs.
  Their ratio is $\cos^2\theta$ (§8.1), and their product
  $\mathbf a'(\mathbf F)\,\mathbf F'(\mathbf a)=\cos^2\theta$ is dimensionless
  and equals 1 exactly when the pair is aligned. Both converters are built
  from a measured vector by taking a reciprocal, and that is where the metric
  enters, once per converter. The wave's $\mathbf k$ arrives ready-made as
  $\nabla\varphi$, which is why the wave never shows the disagreement: the
  phase is only ever computed one way.
- **The exact branch (T1) falls out, holes included.** $\mathbf a'$ covers
  $\mathbf F=\mathbf 0$ (gives $m=0$); $\mathbf F'$ covers $\mathbf a=\mathbf 0$
  (gives $1/m=0$); neither covers $0/0$. The trichotomy's "no mass" becomes
  either "the two totals disagree" (a non-collinear pair) or "they agree on
  a value outside $(0,\infty)$" (anti-aligned, or one of them zero): a
  consistency check between defined values, not a partial function failing
  (corrected the same afternoon, see 9.4; the first wording said
  disagreement only). For an entry,
  which is a behavioural trigger, that is a different procedure. On measured
  inputs it sits with the 29 Sept ruling that mass is never undefined
  (`why_mass_cannot_be_undefined_2026-09-29.md`); on exact inputs T1 is
  untouched.
- **The composition laws are the linearity of the two pairings.** A covector
  is linear in its argument, so
  $\mathbf a'(\mathbf F_1+\mathbf F_2)=\mathbf a'(\mathbf F_1)+\mathbf a'(\mathbf F_2)$
  is lumping (L5), and
  $\mathbf F'(\mathbf a_1+\mathbf a_2)=\mathbf F'(\mathbf a_1)+\mathbf F'(\mathbf a_2)$
  is common force (L14). This explains the 2 Sept finding of
  `what_forces_the_ratio_of_means_readout` that the m-regression passes only
  (L) and the μ-regression only (R): the m-regression *is*
  $\mathbf a'(\mathbf F)$ and the μ-regression *is* $\mathbf F'(\mathbf a)$.
- **What estimation is, in this language** (Claude's reading, not ruled).
  The two pairings coincide exactly on aligned pairs. The candidate cloud is
  the set of pairs on which they coincide, and the ratio of means is the
  average over it that keeps both linearities. So estimation is restoring the
  agreement of the two pairings, and the one non-canonical step of §7, the
  $\mathrm L\cdot\mathrm L^{-1}$ choice of direction, is the variable the law
  averages over.

### 9.2 L14 and the swap: noted, kept provisional on his instruction

A readout that passes lumping and is swap-symmetric ($R(Y,X)=1/R(X,Y)$)
passes common force:

$$
\frac{1}{R(X,Y_1+Y_2)}=R(Y_1+Y_2,X)=R(Y_1,X)+R(Y_2,X)=\frac{1}{R(X,Y_1)}+\frac{1}{R(X,Y_2)}.
$$

In the lemma record the readout swap (T5) is a theorem of the ratio of means
(T4), which is derived from L5 **and** L14, so this is a trade, not a free
derivation: under the premise of §3 (2), that m and $\mathrm m^{-1}$ are
symmetric, the swap becomes a premise and L14 a theorem; otherwise L14 stays
a premise and T5 a theorem. **His ruling:** "keep L14 provisional. The
confusion and the tradeoff with the swap vs the lumping etc. ruling is well
understood and can remain unruled." L14 stands as the lemma record has it.
Nothing is rebound.

### 9.3 Next moves (Claude's, in the order proposed; none begun)

1. Write the exact law in the pairing form and check that it reproduces T1,
   including the reading of "no mass" as disagreement (9.1, second point).
2. Express the misfit $Q$ (L8, L21) in whitened coordinates and see whether
   it is a function of the product $\mathbf a'(\mathbf F)\,\mathbf F'(\mathbf a)$.
   If it is, the misfit weight $g$ gets a units-level reading. Untested.
3. Finish the build-up of §6. The chain lands on
   $m=(\Delta\mathbf v)'(\Delta\mathbf p)$ and $1/m=(\Delta\mathbf p)'(\Delta\mathbf v)$
   with time paired off canonically, which reads a reading as one interval
   over which $\Delta\mathbf p$ and $\Delta\mathbf v$ are taken. May bear on
   O9 (one push measured repeatedly: one reading or several). A guess.
4. O11 by units (§8.4): still parked.

### 9.4 The start on the four, with checks (same afternoon)

On his "ok make a start on the 4 listed", then recorded on his "Do 1 and 2
please". The checks ran on his laptop against the kept
`noise_vs_signal/ours29_reference.py` (`reading_state`, `misfit_weight`).
The check script sat in a scratch folder and is **not yet kept**: item 3 of
the recording plan (beside this record, or in the tester as a test) is
pending. The results are carried into `mass_estimator_lemmas_2026-10-10.md`
(T1 restated, T32, T33, L5, L14, L18, L3′, O9, O11), which supersedes the
29 September lemma record and marks all of it provisional. §9.1's second
point was corrected the same afternoon to match (1) below; its first wording
said disagreement only.

**(1) The exact law in pairing form reproduces T1, with finer structure.**
Two totals and their ratio $P=m_1/m_2=\cos^2\theta$:

| exact pair | $m_1=\mathbf a'(\mathbf F)$ | $m_2=1/\mathbf F'(\mathbf a)$ | $P$ | T1 |
|---|---|---|---|---|
| aligned, $\mathbf F=3\mathbf a$ | 3 | 3 | 1 | one mass |
| anti-aligned | −3 | −3 | 1 | no mass |
| misaligned 30° | 2.598 | 3.464 | 0.75 | no mass |
| perpendicular | 0 | ∞ | 0 | no mass |
| $\mathbf F=\mathbf 0,\ \mathbf a\neq\mathbf 0$ | 0 | undefined | — | no mass |
| $\mathbf F\neq\mathbf 0,\ \mathbf a=\mathbf 0$ | undefined | ∞ | — | no mass |
| $0/0$ | undefined | undefined | — | every mass |

"No mass" is **not** only disagreement: the anti-aligned pair is collinear,
so the two totals agree, on a negative value. $P=1$ means "one projective
value", and T1's "no mass" splits into *disagreement* ($P<1$: the $d-1$
codirectionality conditions failing) and *agreement outside $(0,\infty)$*
(the sign and zero conditions). That is the chain audit's "counting +
positivity" split, computed as two total operations plus a comparison. In 1D
$P\equiv1$ and only the sign condition is left, matching "in $d=1$ it is
zero".

**(2) $Q$ is not a function of the product alone, and should not be.** In
whitened coordinates with $\Sigma=\mathrm{diag}(\sigma_F^2,\sigma_a^2)\otimes I_k$
and $N=I$, the reading's Gram matrix $M$ has $Q=\lambda_-$ (the 30 Sept
eigenvalue dictionary), and the product of the two pairings is its squared
correlation:

$$
P=\frac{M_{12}^2}{M_{11}M_{22}},\qquad Q\,\lambda_+=M_{11}M_{22}\,(1-P),\qquad
Q=\tfrac12\Big[T-\sqrt{T^2-4M_{11}M_{22}(1-P)}\Big],\quad T=\operatorname{tr}M .
$$

Checked against `reading_state` on 300 random readings in $k=1,2,3$:
agreement to $7\times10^{-13}$, and the physical and whitened $P$ coincide
when the channel noise is diagonal. So $Q$ is the dimensionless disagreement
$(1-P)$ converted into noise units of energy by the sizes. The same angle at
three scales: $P=0.770$ throughout, $Q=0.002,\ 0.21,\ 21.1$, so
$g=1,\ 1,\ 0.09$. That is the right behaviour (a small misaligned pair is
consistent with alignment under noise), and it means the misfit weight cannot
be read from the pairings alone: it needs the energy. Two side findings: with
cross-channel correlation $\rho\neq0$, $Q$ sees the *whitened* product, not
the physical one ($P=0.62$ physically; $0.71$ at $\rho=0.5$; $0.97$ at
$\rho=-0.7$), because whitening corrects $\mathbf a$ by $\mathbf F$ through
the conversion $\sigma_a/\sigma_F$; and at equal whitened energies
$Q=s^2(1-\lvert\cos\theta_w\rvert)$ exactly.

**(3) The build-up lands on $p=mv$ and makes the window part of a reading.**
$m=(\Delta\mathbf v)'(\Delta\mathbf p)$ and $1/m=(\Delta\mathbf p)'(\Delta\mathbf v)$,
with time paired off canonically. The quiet condition of §6 is now exact:
for a non-constant push, force and acceleration must be averaged over the
**same window**, or the pair is not a Newton II pair even without noise.
Checked on $F(t)=5(1+\sin2t)$ in a fixed direction, $m=2$: the same window
gives $2.000$; the acceleration window shifted by half gives $1.935$; one
window twice the other gives $2.417$. So a reading is
$(\text{window},\Delta\mathbf p,\Delta\mathbf v)$, and the window belongs to
the reading structure of L18 as a declared instrument fact. **His ruling:**
"I like the window. You can just add that in to the lemma lists."

For **O9** this gives a criterion: two measurements over the same window
have the same signal by construction, so under L19's reason they cannot be
two readings; they are one reading with a combined $\Sigma$. Measurements
over different windows of an "unchanged" push are different signals, hence
separate readings, with L19 as written. The criterion is the window, not the
push. Proposed, not ruled.

**(4) O11 by units: what the swap forces, and where it stops.** Premise: the
default declaration for mass with instruments $(\sigma_F,\sigma_a)$ and the
default for inverse mass with $(\sigma_a,\sigma_F)$ are one declaration read
backwards (the symmetry of §3 (2)). Then $c(\sigma_F,\sigma_a)=1/c(\sigma_a,\sigma_F)$,
and with $c=k\,\sigma_F/\sigma_a$ from L3 this gives $k^2=1$, so $c=s$. The
general form: **$c$ may depend only on declared quantities that invert under
the swap; any swap-invariant dimensionless factor must be 1.** So $c$ cannot
depend on $\rho$, $k$, $d$ or $N$, and for one noise pair, correlated or
not, $c=s$ is forced. Where it stops: several instruments leave every
weighted geometric mean $\prod_g s_g^{w_g}$ (the T20 no-go stands);
anisotropic $\Sigma_F,\Sigma_a$ leave any $h$ of the eigenvalues of
$\Sigma_a^{-1/2}\Sigma_F\Sigma_a^{-1/2}$ with $h(1/\lambda)=1/h(\lambda)$,
among which the determinant ratio $(\det\Sigma_F/\det\Sigma_a)^{1/2d}$ is
the unique power mean. Still parked as far as a ruling goes.

## 10. Sources

- G. W. Hart, *Multidimensional Analysis: Algebras and Systems for Science
  and Engineering* (Springer, 1995); summary at
  https://georgehart.com/research/multanal.html; multipliable matrices as
  described at https://multi-dimpy.readthedocs.io/en/latest/theory/
- J. Janyška, M. Modugno, R. Vitolo, *Semi-vector spaces and units of
  measurement*, arXiv:0710.1313
- H. Whitney, "The Mathematics of Physical Quantities", *Amer. Math. Monthly*
  75 (1968), Parts I–II (cited via arXiv:1503.00564; not read directly)
- C. W. Misner, K. S. Thorne, J. A. Wheeler, *Gravitation* (1973), ch. 2
  (one-forms; cited from Claude's knowledge, not re-checked this session)
- In this discussion: `mass_estimator_lemmas_2026-09-29.md` (F1, L1, L3, L9,
  T1, T2, T20, O11); `what_forces_the_ratio_of_means_readout.md` (the m- and
  μ-regressions); `why_mass_cannot_be_undefined_2026-09-29.md`.
