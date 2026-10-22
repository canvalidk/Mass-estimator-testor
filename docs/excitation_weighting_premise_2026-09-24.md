# Counting force and acceleration as vectors: the weighting the premise forces

**2026-09-24.** Record, filed at the user's request ("write it up as a record
with the premise in; explain how we got to it"). Written by Claude from the
conversation of that evening. It follows two records from the same day in
the same discussion, `old_information_in_the_prior_slot_2026-09-24.md` and
`vonmises_state_2026-09-24.md`. It concerns the 21 September specification
in `mass_estimator_equation/`. It **proposes a change to that equation in 1D
and 3D**. It does not edit that folder and binds nothing.

**In one paragraph.** The estimator averages over every true force and true
acceleration that could lie behind a reading, so it needs a weighting over
those true vectors before any data arrive. The 21 September equation weights
every *length* equally ($df\,d\alpha\,d\Omega$). The user's premise is that
force and acceleration play symmetric roles, and that in Newtonian mechanics
either can take any value. Read for vectors, "no preferred value" means the
weighting does not change when every candidate true acceleration is shifted
by the same vector. The only weighting with that property is the flat one on
the vector itself, the weighting Maxwell–Boltzmann uses for velocities. In
the equation's own variables it is $\alpha^{d-2}\,df\,d\alpha\,d\Omega$.
**In 2D nothing changes. In 3D the current equation favours small true
accelerations, and in 1D large ones.** After the change the equation has
the same simple form in every dimension. The tests so far show no clear
winner, so the change is justified by the argument, not by performance.

---

## 1. The premise

Two statements by the user, 2026-09-24 (quoted with light cleanup):

- **P1 (symmetry).** "There is no difference between $\alpha$ or $f$. You
  are arbitrarily picking one. But there is a symmetry here." Any rule for
  weighting one must apply equally to the other.
- **P2 (no preferred value).** "No, alpha does not depend on anything.
  Neither does net force. In NM, they can have any numbers."

**Reading P2 for vectors** (proposed by Claude, accepted by the user).
Force and acceleration add as vectors. So "no value is preferred" means:
shift every candidate true acceleration $\mathbf a^\ast$ by the same vector
$\mathbf b$, and the weighting before data looks exactly the same. Adding a
push should not change which pushes you thought were likely.

**What the premise is about, and what it is not about.** The user asked:
*"Does that mean, if you have a prior of some amount and you apply it to the
equation, it doesn't change from the prior?"* No. There are two separate
priors:

- **The excitation prior** is the weighting over the true force and true
  acceleration *before* any reading. The premise is about this one.
- **The mass prior** ($m_0$, $\sigma_{m_0}$ in the earlier records) is the
  user's to set. The premise says nothing about it.

Data update the law exactly as before. The shift-invariance is a property of
the starting weighting only, as a flat ruler looks the same wherever you
slide it.

## 2. What it forces

**Step 1: the only shift-invariant weighting is the flat one.** On
$\mathbb R^d$ the only measure unchanged by every shift is ordinary volume,
$d^d\mathbf a^\ast$, up to a constant (Lebesgue measure; the uniqueness of
Haar measure). Any weighting that favours short vectors, or long ones, fails
the test, because after the shift the favoured region has moved.

**Step 2: P1 is respected.** P1 is used in two strengths. The weak one:
choosing $f$ or $\alpha$ changes the excitation weighting only by a factor
in $m$. The strong one: the whole measure, mass factor included, is
unchanged when $f$ and $\alpha$ swap roles. Steps 1–3 need only the weak
one. At fixed mass $\mathbf F^\ast = m\,\mathbf a^\ast$,
so $d^d\mathbf F^\ast = m^d\, d^d\mathbf a^\ast$. Choosing to make the force
flat instead of the acceleration changes the weighting only by a factor that
depends on $m$ alone, and a factor in $m$ alone belongs to the mass prior. So
P1 and P2 together fix the excitation weighting, and which of $f$ or
$\alpha$ is picked does not matter.

**Step 3: in the equation's variables.** With $\mathbf a^\ast=\alpha\mathbf u$
and $f=m\alpha$ (so $df = \alpha\,dm$ at fixed $\alpha$):

$$d^d\mathbf a^\ast\,dm \;=\; \alpha^{d-1}\,d\alpha\,d\Omega\;dm \;=\; \alpha^{d-2}\;df\,d\alpha\,d\Omega .$$

Read the other way, the 21 September measure is

$$df\,d\alpha\,d\Omega \;=\; |\mathbf a^\ast|^{\,2-d}\; d^d\mathbf a^\ast\,dm ,$$

so it weights the true acceleration vector by $|\mathbf a^\ast|^{2-d}$:

| dimension | 21 September ($df\,d\alpha\,d\Omega$) | forced (flat in the vector) | change |
|---|---|---|---|
| 1D | weight $\propto\lvert\mathbf a^\ast\rvert$: disfavours small pushes | flat | multiply by $1/\alpha$ |
| 2D | flat | flat | **none** |
| 3D | weight $\propto 1/\lvert\mathbf a^\ast\rvert$: favours small pushes | flat | multiply by $\alpha$ |

**The picture for 3D.** Hold the mass fixed. True accelerations of length 1
to 1.1 fill a thin shell, and those of length 10 to 10.1 fill a shell with
100 times the volume. The flat weighting gives the outer shell 100 times the
weight. The 21 September measure gives it only 10 times: at fixed mass its
weight per unit length grows like $\alpha$ (from $df=\alpha\,dm$), not
$\alpha^2$. So it packs the inner shell 10 times more densely, a hidden
preference for small pushes. In general, at fixed mass the forced weighting
on $\alpha$ is $\alpha^{d-1}$, the shell factor, while the 21 September one
is $\alpha^{1}$ in every dimension; the two agree only in 2D.

**What the premise leaves open: a factor in $m$ alone.** Steps 1–3 fix the
weighting up to a function of $m$, and that function is the mass prior,
which is open (investigation question 2). For the tester, and for the
formulas below, it is chosen as follows. Write the pair in noise units as one
vector $\mathbf v = r\,\mathbf u\in\mathbb R^d$ with $r^2 = (f/\sigma_F)^2+(\alpha/\sigma_a)^2$,
and $m = s\tan\theta$ with $s=\sigma_F/\sigma_a$. Then

$$d\nu_{\text{new}} \;=\; d^d\mathbf v\,d\theta \;\propto\; r^{\,d-2}\;df\,d\alpha\,d\Omega ,$$

which is $\alpha^{d-2}\,df\,d\alpha\,d\Omega$ times a function of $m$ alone.
This choice:

- is flat in $\theta$ in every dimension;
- matches the 21 September equation exactly in 2D;
- respects P1 in the strong sense, in noise units: swapping $f/\sigma_F$ with $\alpha/\sigma_a$
  sends $\theta\to\pi/2-\theta$ and leaves the measure unchanged.

The plain choice $\alpha^{d-2}\,df\,d\alpha\,d\Omega$ is flat in $m$
instead. It respects P1 in the strong sense only in 2D, and in 1D it leaves
the single-reading law improper (its tail falls only like $1/m$). It is
written out in §3 for comparison.

## 3. The equation after the change, in the measured quantities

Symbols: measured $\mathbf F_i$, $\mathbf a_i$; noise per component
$\sigma_F$, $\sigma_a$; dimension $d$.

**One reading.** Integrating the true acceleration vector out is now a
Gaussian integral in every dimension:

$$p(m)\;\propto\;\frac{1}{\sigma_F^2+m^2\sigma_a^2}\;
\exp\!\left[-\,\frac{\big|\mathbf F-m\,\mathbf a\big|^2}{2\,\big(\sigma_F^2+m^2\sigma_a^2\big)}\right].$$

The exponent measures how well Newton II fits. The front factor is what
"flat in $\theta$" becomes in $m$ ($d\theta/dm$): a half-Cauchy with scale
(and median) $s$. With the flat-in-$m$ completion the front factor would be
$(\sigma_F^2+m^2\sigma_a^2)^{-d/2}$ instead.

**Many readings of the same mass** (the symmetric rule, mass factor counted
once):

$$p(m)\;\propto\;\text{prior}(m)\;
\exp\!\left[-\,\frac{\sum_i\big|\mathbf F_i-m\,\mathbf a_i\big|^2}{2\,\big(\sigma_F^2+m^2\sigma_a^2\big)}\right].$$

Readings combine by adding residuals, in every dimension. A $d$-dimensional
reading is exactly $d$ one-dimensional readings, one per axis. The
doubled-angle vector state of `vonmises_state_2026-09-24.md`, which was
exact only in 2D, becomes exact in 1D and 3D as well.

**The readout is unchanged in form**: ratio of means,
$\widehat m = E[f]/E[\alpha]$. The one new ingredient is the expected true
acceleration size at mass $m$:

$$\hat{\mathbf a}_i(m)=\frac{m\,\sigma_a^2\,\mathbf F_i+\sigma_F^2\,\mathbf a_i}{\sigma_F^2+m^2\sigma_a^2},\qquad
\sigma_e(m)=\frac{\sigma_F\,\sigma_a}{\sqrt{\sigma_F^2+m^2\sigma_a^2}},\qquad
\bar\alpha(m)=\sum_i \sigma_e(m)\,\bar\chi_d\!\left(\frac{|\hat{\mathbf a}_i(m)|}{\sigma_e(m)}\right),$$

$$\widehat m=\frac{\int_0^\infty m\,\bar\alpha(m)\,p(m)\,dm}{\int_0^\infty \bar\alpha(m)\,p(m)\,dm},$$

where $\bar\chi_d$ is the mean length of a unit Gaussian vector centred at
distance $h$ (the noncentral chi mean):

$$\bar\chi_1(h)=\sqrt{\tfrac{2}{\pi}}\,e^{-h^2/2}+h\,\operatorname{erf}\tfrac{h}{\sqrt2},\qquad
\bar\chi_3(h)=\sqrt{\tfrac{2}{\pi}}\,e^{-h^2/2}+\Big(h+\tfrac1h\Big)\operatorname{erf}\tfrac{h}{\sqrt2},$$

and $\bar\chi_2$ is the Rice mean.

For comparison, the 21 September law in 3D has a per-reading kernel
$\propto e^{h^2/2}\operatorname{erf}(h/\sqrt2)/h$ in the doubled-angle
variable. The extra factor, about $1/h$ for a strong reading, is the
$|\mathbf a^\ast|^{-1}$ weighting of §2 showing through.

## 4. The same move as Maxwell–Boltzmann

The user's observation: this is what is done for the Maxwell–Boltzmann
distribution.

- **The base weighting.** Maxwell–Boltzmann starts from equal weight per box
  of *velocity* space, $d^3\mathbf v$, and multiplies by
  $e^{-mv^2/2kT}$. The $4\pi v^2$ in the speed distribution comes only from
  converting "per box of velocity" into "per unit of speed". The shell
  factor is $v^{d-1}$: Maxwell in 3D, Rayleigh in 2D, half-normal in 1D.
- **The mistake the old weighting makes, in gas language.** At fixed mass
  the forced weighting on $\alpha$ is $\alpha^{d-1}$, exactly Maxwell's
  shell factor. The 21 September weighting is $\alpha^1$ in every dimension:
  in 3D that is replacing Maxwell's $v^2$ with $v$, which puts too many
  molecules near rest. The power reads $\alpha^{d-2}$ in the $(f,\alpha)$
  variables only because $df=\alpha\,dm$.
- **The hierarchical rule is Maxwell–Boltzmann for accelerations.** The
  tester's hierarchical rule (commit `0c3235f`) gives each reading's true
  pair a Gaussian $\mathbf v_i\sim N(0,\omega^2 I)$ with one shared,
  learned $\omega$. That is a Maxwell–Boltzmann law with $\omega^2$ playing the
  role of $kT/m$. The flat weighting is its infinite-temperature limit
  ($\beta=\omega^2/(1+\omega^2)\to1$). With many readings, that limit is
  what makes the law about $1/\beta$ too concentrated, and so it
  undercovers.

## 5. Where the argument is thin, and a second route found while writing

**The weak spot.** In statistical mechanics the flat $d^3\mathbf v$ is not
only "no preferred value". Two physical facts back it:

- Liouville's theorem: the dynamics preserves phase-space volume.
- Galilean boosts: shifting every velocity by the same vector is a real
  symmetry of the physics, and it is exactly the shift test.

For accelerations, shifting every value by the same vector is what moving
to a uniformly accelerating frame does to *coordinate* accelerations (those
from tracking position). That is not a symmetry of this measurement: in the
new frame $\mathbf F = m\mathbf a$ needs a fictitious force $-m\mathbf A$
that no force sensor registers. An accelerometer, which measures proper
acceleration, does not shift at all. Worse, zero acceleration *is*
physically distinguished in Newtonian mechanics: by the first law,
force-free bodies do not accelerate in an inertial frame. So for
us the translation step rests on P2 read as a statement of ignorance, without
the dynamical backing Maxwell–Boltzmann has. A careful examiner would press
here.

**A second route that avoids the weak spot** (Claude, found while writing
this record; not yet discussed with the user). Maxwell's own first
derivation (1860, after Herschel 1850) did not use translation. It used two
assumptions:

- **Isotropy.** No direction is preferred. This is already in the 21
  September equation, as the uniform $d\Omega$.
- **Independent components.** The push along $x$ tells you nothing about the
  push along $y$.

The Herschel–Maxwell theorem says that for $d\ge2$ these two force a centred
Gaussian, $N(0,\omega^2 I)$, with one unknown scale. That is the hierarchical
prior of §4, centred at zero, which is consistent with zero being special
(the first law). The flat weighting is its $\omega\to\infty$ limit.

So the two routes agree on the base weighting, and both reject the 21
September weighting in 3D: $|\mathbf a^\ast|^{-1}$ is neither flat nor
Gaussian. The second route also says more. It suggests the *repeated-reading*
estimator should carry a finite, learned $\omega$ rather than the flat
limit, and that is what the coverage tests of §6 favour.

The second route has limits of its own:

- **1D.** In 1D isotropy is only $\pm$ symmetry and there is nothing for
  independence to act on, so the theorem gives nothing directly. It does
  apply if the 1D reading is one axis of a push that really lives in 2D or
  3D: that component inherits the Gaussian.
- **Improper weightings.** The theorem forces a Gaussian among *proper*
  laws. Among improper weightings, $e^{c|\mathbf a|^2}$ with $c\ge0$ also
  passes both tests: $c=0$ is the flat one, and $c>0$ favours ever-larger
  pushes. The 21 September $|\mathbf a^\ast|^{-1}$ still fails, since it
  does not split into independent components.
- **Independence.** Independence of components was the criticised step of
  Maxwell 1860, which he replaced with a collision argument in 1867. For an
  experimenter's pushes it is a modelling claim, not a law.
- **Which vector.** The tester places the Gaussian on the noise-unit pair
  $\mathbf v$, which keeps P1. A Gaussian placed on $\mathbf a^\ast$ alone
  would not.

## 6. What the tests say

All runs are in the tester repo `canvalidk/Mass-estimator-testor`. "Within
×k" is the percentage of series whose point estimate lies within a factor k
of the true mass. Masses are in units of $s=\sigma_F/\sigma_a$.

**One reading** (`presets/reference_single_reading.yaml`, commit `cd53203`,
1D and 3D, 30 cells × 2000 series, fixed push size, fixed direction). The
two weightings trade places at weak signal and agree once the push is strong:

| cell (SNR 0.5) | 21 Sept, within ×2 | forced, within ×2 |
|---|---|---|
| 1D, $m=s/4$ | 10.9% | 7.1% |
| 1D, $m=s$ | 81.6% | 87.6% |
| 3D, $m=s/4$ | 10.1% | 12.7% |
| 3D, $m=s$ | 82.9% | 78.0% |

- From acceleration SNR 4 upward, every cell agrees within about 1
  percentage point.
- Worst-case regret of the ratio of means: 4.8 vs 10.6 points at ×1.25, and
  8.3 vs 8.6 at ×2 (forced vs 21 Sept).
- 95% interval coverage: 94.6–99.7% for the forced weighting, 94.7–99.5% for
  21 Sept.
- These worlds use a fixed push size, which is neither weighting's world, so
  a trade with no clear winner is the expected outcome. The smaller worst
  case at ×1.25 is the one edge, and it belongs to the forced weighting.

**Fifty readings sharing only the mass** (`presets/reference_repeated.yaml`,
3D, 6 cells × 1000 series).

- Point accuracy is the same within about 1 point in every cell.
- Both weightings undercover at weak signal. 80% intervals cover 57.2–80.1%
  (21 Sept) and 52.0–79.8% (forced).
- The change of base weighting does not fix this, and on this evidence
  should not be expected to.

**The same fifty readings with the Maxwell–Boltzmann (hierarchical) prior**
(`presets/excitation_hierarchy.yaml`, commit `0c3235f`, 12 cells × 500
series).

- In the weak cells ($m=s/4$ and $s$, acceleration SNR 1 or Gaussian RMS 1),
  80% intervals cover 81.4–87.6%, against 56.4–60.4% (21 Sept, symmetric
  rule) and 52.8–56.8% (forced, flat limit).
- Worst-case regret at ×1.25: 2.4 points, against 7.4 and 9.0.

**Summary.** The base weighting is decided by the argument. The tests
neither clearly favour nor oppose it. What the tests do favour, for repeated
readings, is a finite temperature on top of it (§5, second route).

**Checks made for this record.** The single-reading law of §3 matches the
tester's `cartesian` reference to $3\times10^{-15}$ in 1D, 2D and 3D. The
flat-in-$m$ completion matches direct quadrature over the true acceleration
components to $1\times10^{-14}$. The tester's `flat` reference matches direct
quadrature of the 21 September measure $df\,d\alpha\,d\Omega$ in 3D to
$10^{-8}$.

## 7. How we got here

1. **Morning.** The old × new split of the exact combination put the old
   readings' law in the prior slot (first record of the day).
2. **Afternoon.** The user asked whether the uncertainties of $f$ and $a$
   could be written in $\lambda$'s natural form. Result: each reading is a
   von Mises law in the doubled angle, and readings add as vectors (second
   record). In 1D and 3D a leftover factor (about $h^{2-d}$ for strong readings)
   spoiled the exactness.
3. **The user asked whether it is wrong to think of force and acceleration
   in 3D,** since what is measured is a projection. Chasing that leftover
   factor showed that under a weighting flat in the vector's *components*
   it disappears. A $d$-dimensional reading then becomes exactly $d$
   projection readings. This was built into the tester as the `cartesian`
   reference (`cd53203`) and tested: no clear winner. It was rendered in his symbols,
   with the prior as $(m_0,\sigma_{m_0})$ and $\lambda=1/\sigma_{\ln m}^2$
   at $m_0=s$.
4. **The user: every freedom is a prior.** The freedoms split into a mass
   prior (a factor in $m$ alone) and an excitation prior (on the true
   pushes, with data-dependent effect). He asked to look at the excitation
   pair.
5. **The excitation pair.** The Gaussian hierarchical prior explained the
   undercoverage with many readings (the flat weighting is its $\beta\to1$
   limit) and repaired it (`0c3235f`). Declared oracles (`2128823`) showed
   that the weak-signal bottom is limited by information and by the mass
   prior, not by method.
6. **"You sound like you let go of the wall."** The user judged, correctly, that this had
   drifted into choosing weightings by test performance, and that the
   project's method is to find what is *forced*. He supplied P1 and P2
   (§1).
7. **The forcing argument** (§2) followed from P1 and P2 in a few lines.
   It picks the weighting of step 3, which was reached there for its
   elegance. So step 3 now has a reason, not only a test result.
8. **The user asked whether the premise meant the prior is not changed by
   data.** The answer is §1: it is about the excitation weighting before
   data.
9. **The user connected it to Maxwell–Boltzmann** (§4). Claude named the
   weak spot. Writing this record turned up the Herschel–Maxwell route
   (§5).

## 8. What changes, what does not, what is owed

- **Changes (proposed):** the per-reading law and the conditional mean
  $\bar\alpha(m)$ in 1D and 3D (§3).
- **Unchanged:**
  - everything in 2D;
  - Newton II as the premise;
  - the uniform direction;
  - the ratio-of-means readout;
  - how readings are combined.
- **The 21 September text** says flat "does not mean flat in mass, log mass
  or Cartesian vector volume". Under this record the last clause would be
  reversed for the excitation. Whether and when to amend that specification
  is the user's decision; it lives on VD-Newton's surface.
- **The tester's default** is still the `flat` reference. Presets choose the
  reference explicitly, so nothing else changes until a preset asks.
- **Owed:**
  - the mass-only factor, i.e. the mass prior (investigation question 2);
  - whether repeated-reading estimation should adopt the finite-temperature
    (hierarchical) prior as its default, and the hyperprior $b$ on $\beta$;
  - an answer to the accelerating-frame objection of §5;
  - the two-cart composition test under the new weighting;
  - a second seed;
  - real data.
- **A pointer that now makes sense.** The unsanctioned literature search in
  `literature_review/` described the target as a Gaussian errors-in-variables
  likelihood with an improper prior $\|\mathbf v\|^{2-d}$. That is the power
  in the table of §2, and the change here sets it to zero.

## 9. Reading list

Cited from general knowledge; not verified against the texts. Nothing here
is in `papers/`.

**Ignorance priors from symmetry.**
- Jaynes (1968), *Prior probabilities*, IEEE Trans. Systems Science and
  Cybernetics SSC-4, 227–241. The method of transformation groups: a
  location parameter with no preferred value gets the flat prior. §2 applies
  it to a vector.
- Jaynes (2003), *Probability Theory: The Logic of Science*, ch. 12.
- Jeffreys (1961), *Theory of Probability*, 3rd ed.

**Maxwell–Boltzmann and its derivations.**
- Herschel (1850), review of Quetelet, *Edinburgh Review* 92; Maxwell
  (1860), *Illustrations of the dynamical theory of gases*, Phil. Mag. The
  isotropy + independence derivation.
- Maxwell (1867), *On the dynamical theory of gases*, Phil. Trans. The
  collision argument that replaced it.
- Bryc (1995), *The Normal Distribution: Characterizations with
  Applications*, Lecture Notes in Statistics 100. The Herschel–Maxwell
  theorem stated carefully.

**Why the flat (infinite-temperature) limit fails with many readings.**
- Neyman & Scott (1948), *Consistent estimates based on partially
  consistent observations*, Econometrica 16, 1–32. The incidental-parameter
  problem: one latent push per reading.
- Lindley (1947), *Regression lines and the linear functional relationship*,
  JRSS Suppl. 9; Solari (1969), JRSS B 31, 372–375. The functional
  (flat-latent) errors-in-variables model and its pathologies.
- Kendall & Stuart, *The Advanced Theory of Statistics*, vol. 2, chapter on
  functional and structural relationships. Structural means Gaussian
  latents: our hierarchical rule.
- Berger, Liseo & Wolpert (1999), *Integrated likelihood methods for
  eliminating nuisance parameters*, Statistical Science 14, 1–28.

## 10. Reproduction

```sh
# tester repo canvalidk/Mass-estimator-testor
python -m pytest -q tests/test_cartesian.py tests/test_hierarchical.py
python -m met run presets/reference_single_reading.yaml   # §6, one reading
python -m met run presets/reference_repeated.yaml         # §6, fifty readings
python -m met run presets/excitation_hierarchy.yaml       # §6, hierarchical
```

Code: `met/law.py`, `ReadingSet.reading_terms(..., reference="cartesian")`
and `cartesian_mean_radius`; `met/analyst.py`, `_hierarchical_log_density`.
