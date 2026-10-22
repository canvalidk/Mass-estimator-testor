# Readings are evidence, just less: tempering by the off-line misfit

**2026-09-29 (night).** Record, written by Claude at the user's request ("Yes
go ahead", to the offer of a record after "Yes investigate various ways. Try
to make sense of it."). The property (§1.5) and the ruling (§1.7) are his,
quoted in §8. The derivation, the tests and the reading of them are Claude's.

**Terminology.** Ruled by the user tonight: the true, errorless vector is the
**signal**, not a push ("No more using the word push for this. Its signal.").
This record uses *signal* throughout. Earlier records keep their wording.

**What it records.** ours28 has a defect separate from the $\sqrt N$ tilt: a
reading whose force and acceleration are misaligned far beyond their declared
noise (a reading the premises call impossible) is treated as strong signal,
and the law collapses onto the least-bad mass with unbounded confidence as
$\sigma\to0$. The user ruled that such readings are evidence, just less. Under
that ruling a reading's factor must be *flattened* by its off-line misfit,
$E_i(m)^{g(Q_i)}$. Z fixes $g(0)=1$. Bounded confidence plus "evidence, not
none" fixes the tail as $g\sim\kappa/Q$. The middle is free. Tested: every
form with a $1/Q$ tail fixes the impossible reading (one such reading drops
ours28's coverage from 97% to 8.5%; tempered, 97–98%), and forms that stay at
$g\approx1$ inside the noise band cost 2–6% width on normal data.

**Relations.** Scopes `mass_estimator_lemmas_2026-09-29.md` (L8, L13/T23, L20,
R1, T24, T25, O6, O7), `alignment_what_one_reading_can_say_about_noise_2026-09-29.md`
§3 (whose "the residual is weakly pro-signal" is overridden here by ruling,
outside the noise band), `null_series_tilt_and_property_N_2026-09-29.md` (not
touched: see §6), `ours26_learned_omega_is_a_hint_2026-09-27.md` (its
tempering identity is the mechanism, now driven by each reading's own misfit
rather than a learned $\omega$), and the 7 Sept ruling (not a gate: §6). It
binds nothing.

---

## 1. How it was reached

The evening's thread, in order. 1D and 2D, unit noise, $c=1$, ours27 and ours28
as filed, sandbox.

### 1.1 Signal, error, and the triangle

Measured $=$ signal $+$ error, as vectors, one triangle per channel. The
signal-to-noise ratio is $|\text{signal}|/\sigma$ (in noise units with a full
$\Sigma$). In ours28 $\omega$ is the signal's typical size in noise units, and
the typical sizes obey Pythagoras: legs $\omega$ and $1$, hypotenuse
$\sqrt{1+\omega^2}$. With $\tan\varphi=\omega$,

$$
\beta=\frac{\omega^2}{1+\omega^2}=\sin^2\varphi,\qquad 1-\beta=\cos^2\varphi .
$$

$\beta$ is the signal's share of the measured variance, the shrink factor
(best-guess signal $=\beta\times$ measured, along the candidate line), and the
leftover uncertainty of the signal. The arcsine law on $\beta$ is uniform in
$\varphi$. One realized triangle is not right-angled; the Pythagoras holds for
typical sizes.

### 1.2 Does the back-calculated error add up?

For a candidate $(m,\mathbf a^\ast)$ the error vectors are exact:
$\mathbf e_F=\tilde{\mathbf F}-m\mathbf a^\ast$,
$\mathbf e_a=\tilde{\mathbf a}-\mathbf a^\ast$, size
$Q(m,\mathbf a^\ast)=\mathbf e^\top\Sigma^{-1}\mathbf e$. Averaged over the
equation's own law, the error energy per reading should equal the declared 2
(1D, two channels). 100 readings, 40 series:

| data | ours27 | ours28 |
|---|---|---|
| pure noise | 1.93 | **1.72** |
| signal of 1 noise unit | 1.96 | 1.99 |
| signal of 3 noise units | 2.02 | **2.65** |

ours28 books noise as signal on null data and signal as noise at 3 units.
The best single guess always leaves less than declared (ours27 about 1,
ours28 1.2–1.5): any fit absorbs noise. Pure noise and fixed-size series are
not ours28's own world, so this is a property it is not strictly obliged to
meet there. Under R1 it is a candidate property. **Not ruled.**

### 1.3 The averaged true pair is not a Newton II pair

Every candidate has true force $=m\times$ true acceleration, codirectional.
The averages over the law of mass are not: the force average is weighted by
$m$ and leans toward the measured force.

| measured $\tilde F$, $\tilde a$ | estimated true force | estimated true acceleration |
|---|---|---|
| 1D: $+1$, $-1$ (ours27) | $+0.206$ | $-0.206$ |
| 1D: $+1$, $-1$ (ours28) | $+0.069$ | $-0.069$ |
| 2D: $(2,0)$, $(0,2)$, both | | $25^\circ$ apart |
| 2D: $(3,0)$, $(1,2)$ | | $8^\circ$ (ours27), $10^\circ$ (ours28) |

They agree only when the measured vectors are aligned or the law of mass has
collapsed. Recorded as a fact about the averages; nothing yet uses it.

### 1.4 The plug-in signal share does not win

Requiring the along-line error to carry its declared energy $k$ gives
$\hat\beta(m)=\max\!\big(0,\,1-k/H(m)\big)$: the single-reading maximum
likelihood, and the positive-part James–Stein factor. Slotted in per reading,
the null-series tilt (1D, median of 60 series; 20 at $N=1000$):

| $N$ | ours28 | ours28 with $\hat\beta$ |
|---|---|---|
| 10 | 1.05 | 1.07 |
| 100 | 4.40 | 5.28 |
| 1000 | 15.6 | 16.3 |

As the null-series record's §5 predicts for any fixed per-reading rule. Pooled
over readings, $1-kN/\sum H$, it is ours28b's shared scale (L19).

### 1.5 His property

"If the vectors are aligned and the uncertainty decreases, signal increases.
If the directions are unaligned and the uncertainty decreases, signal
decreases." ours28's signal share (mean $\beta$ over the law), 2D, $k=2$:

| $\sigma$ | aligned $(2,0),(1,0)$ | perpendicular $(2,0),(0,1)$ | off-line misfit $Q$, perpendicular |
|---|---|---|---|
| 2 | 0.28 | 0.27 | 0.3 |
| 1 | 0.41 | 0.35 | 1 |
| 0.5 | 0.82 | 0.75 | 4 |
| 0.2 | 0.98 | 0.97 | 25 |
| 0.1 | 0.99 | 0.99 | 100 |

The first half holds. The second fails: ours28's signal share reads only the
energy along the best line, $H$, and never the off-line misfit $Q$.

### 1.6 A defect, but not the measured overconfidence

It is not the source of the coverage failures or the $\sqrt N$ tilt: 1D has no
misalignment and carries the tilt in full, and the own-world coverage tests
contain no readings misaligned beyond the noise. It is a defect of its own:
readings with misalignment far beyond $\sigma$ have probability zero under
L1 and R1, the definition must still answer (never undefined, ruled 29 Sept),
and what ours27 and ours28 answer (a law collapsing on the least-bad mass
with unbounded confidence) was chosen by nobody. Three ways to speak there,
each costing a premise, were put to him: misfit as evidence the $\sigma$'s are
small (R1), as evidence Newton II failed (7 Sept), or a trust flag only.

### 1.7 His ruling

"The readings are evidence. Just less."

## 2. What the ruling forces

Write each whitened reading as $(T_i,C_i)$ (the alignment record's §1). Along
the best line it has energy $\max_\theta H_i=\tfrac12(T_i+|C_i|)$; off it,

$$
Q_i=\tfrac12\big(T_i-|C_i|\big),
$$

which is $m$-free. Under the model $Q_i$ is at most $\chi^2_{\nu}$ with
$\nu=d-1$ (smaller under noise, as the fit absorbs part of it).

**Less means flatter, not smaller.** Multiplying $E_i(m)$ by any function of
$Q_i$ changes nothing in the law, since $Q_i$ does not depend on $m$. The only
way for a reading to count for less is a power below 1:

$$
p(m)\propto\pi(m)\prod_iE_i(m)^{\,g(Q_i)},\qquad g(0)=1,\ g\ \text{decreasing}.
$$

This is the tempering identity of the learned-$\omega$ records, driven by the
one number the law currently discards, per reading, from that reading alone.

**The two ends.**

- **$g(0)=1$.** A reading with no off-line misfit keeps full weight. An exact
  zero reading has $T_i=|C_i|=0$, so $Q_i=0$, $g=1$ and $E_i$ constant: Z
  holds (T24 unchanged). Exact zero components add nothing to
  $|\tilde{\mathbf F}|^2$, $|\tilde{\mathbf a}|^2$ or $\tilde{\mathbf F}\cdot\tilde{\mathbf a}$,
  so $Q_i$ is unchanged and T25 survives (analytic).
- **The tail.** Fix a misaligned pair and let $\sigma\to0$. $T_i$, $|C_i|$ and
  $Q_i$ all grow like $1/\sigma^2$, and at large energy the reading's log-tilt
  is $|C_i|/2$. The tempered tilt $g(Q_i)\,|C_i|/2$ stays bounded iff
  $g=O(1/Q)$. A faster tail drives the reading's evidence to zero, which is
  "none", not "less". So the ruling plus bounded confidence fix the tail's
  power: $g\sim\kappa/Q$. With $\kappa=\nu$ the limit is

$$
\lim_{\sigma\to0}\ \text{tilt}_i=\nu\,\frac{a_i}{1-a_i},\qquad a_i=\frac{|C_i|}{T_i},
$$

  set by geometry alone: zero for a reading no mass can explain at all,
  unbounded as the reading aligns. The coefficient $\kappa$ is not forced;
  $\nu$ has the reading "the share of the residual the declared noise accounts
  for". **Sharpened in §10: $\nu=k_i-1$, the number of codirectionality
  conditions inside the declared subspace, and $Q_i$ is computed there.**

**What is excluded.**

- Scale-free discounts, such as a function of $a_i$ alone: they do not change
  as $\sigma$ shrinks, so they cannot stop the collapse.
- A heavy-tailed error model (Student-$t$) would give "less evidence with
  misfit" by itself, but it lets an exact zero reading tilt the law (the 1.6%
  of L13), breaking Z. So T23's forcing of Gaussian noise holds for the
  along-line shape and says nothing about how the residual is treated.

**Where it does nothing.** In 1D, $|C_i|=T_i$, so $Q_i\equiv0$ and $g\equiv1$.

## 3. Forms tested

| name | $g(Q)$ |
|---|---|
| ours28 | $1$ |
| smooth | $1/(1+Q/\nu)$ |
| kink | $\min(1,\ \nu/Q)$ |
| p-value | $P(\chi^2_\nu>Q)$ (tail faster than $1/Q$: evidence goes to zero) |
| Cauchy $3\nu$ | $1/\sqrt{1+(Q/3\nu)^2}$ |
| Cauchy $9\nu$ | $1/\sqrt{1+(Q/9\nu)^2}$ |

## 4. Results (2D, $k=2$, ours28 base)

**4.1 The impossible pair as $\sigma\to0$.** $\tilde F=(2,0)$,
$\tilde a=(0,1)$, fixed; tilt of the log law in nats.

| $\sigma$ | $Q$ | ours28 | smooth | kink | p-value | Cauchy $3\nu$ | Cauchy $9\nu$ |
|---|---|---|---|---|---|---|---|
| 1 | 1 | 0.51 | 0.26 | 0.51 | 0.16 | 0.49 | 0.51 |
| 0.5 | 4 | 3.79 | 0.76 | 0.95 | 0.17 | 2.27 | 3.46 |
| 0.1 | 100 | 148 | 1.46 | 1.48 | 0.00 | 4.44 | 13.3 |
| 0.005 | 40 000 | 60 000 | 1.50 | 1.50 | 0.00 | 4.50 | 13.5 |

Here $a=0.6$, so $\nu a/(1-a)=1.5$: the $\kappa=\nu$ forms reach it exactly;
the Cauchy forms reach $\kappa a/(1-a)$ with $\kappa=3\nu,9\nu$.

**4.2 The price on own-world data** (signal of $r$ noise units along random
directions, true mass drawn over the central range, $\sigma$'s correct).
Coverage of the 95% interval / median log-width:

| world | ours28 | kink | Cauchy $3\nu$ | Cauchy $9\nu$ | smooth | p-value |
|---|---|---|---|---|---|---|
| $N=20$, $r=1$ | 94.5 / 2.53 | 95.5 / 2.68 | 95.5 / 2.60 | 94.5 / 2.54 | 98.0 / 3.15 | 98.5 / 3.58 |
| $N=20$, $r=3$ | 97.5 / 0.44 | 97.0 / 0.48 | 98.0 / 0.46 | 97.5 / 0.44 | 99.0 / 0.55 | 99.5 / 0.64 |
| $N=100$, $r=1$ | 94.2 / 0.82 | 96.7 / 0.85 | 95.8 / 0.84 | 94.2 / 0.83 | 98.3 / 1.01 | 98.3 / 1.19 |

Forms that keep $g\approx1$ through the noise band cost 0–6% width. Forms
that already discount at $Q\approx\nu$ cost 25–45% and over-cover. (This 2D
baseline is not the no-borrowing record's 84% case, which was 3D with 200
readings; the numbers differ, not the picture.)

**4.3 Rescue.** 20 readings at $r=3$ plus one reading at the impossible pair
with $\sigma=0.1$ ($Q=100$):

| | coverage | median log error |
|---|---|---|
| ours28 | **8.5%** | 0.83 |
| smooth | 99.0% | 0.083 |
| kink | 98.0% | 0.087 |
| p-value | 100% | 0.084 |
| Cauchy $3\nu$ | 97.0% | 0.089 |
| Cauchy $9\nu$ | 91.0% | 0.085 |

With a mildly misaligned extra reading ($\sigma=0.5$, $Q=4$) ours28 already
covers 97% and tempering changes little: the model absorbs what it can.

## 5. Reading

- **Both ends are forced; the middle is free.** $g(0)=1$ by Z. $g\sim\kappa/Q$
  at large $Q$ by the ruling plus bounded confidence. Between them (misfits
  of a few $\nu$) the shape is an F2 item, but a narrow one: it acts only on
  readings between normal and impossible.
- **The own-world test prefers "no discount inside the noise band".** That is
  also what the alignment record found the model to license there. The
  override of its §3 is then confined to where the model licenses nothing.
- **A reading's evidence is scaled by how much of its residual its instrument
  can vouch for.** A perfectly measured misaligned pair is worth
  $\kappa\,a/(1-a)$ nats, a finite amount fixed by its geometry.

## 6. What it does not do

- **The $\sqrt N$ tilt.** Under noise $Q_i$ sits below $\nu$ (mean about 0.43
  in 2D), where $g\approx1$. The null-series problem is untouched.
- **The readout (T10).** Only the law was tested.
- **3D.** Not tested; $\nu=2$ there.
- **A gate.** $g>0$ always: every reading keeps some evidence and nothing
  halts. Compatible with the 7 Sept ruling as far as its letter goes; whether
  its intent covers a soft discount is his call.

## 7. Parked (delivered, not dropped)

1. **Axis splitting (L8) collides.** A 2D reading has a residual; its 1D axis
   readings have none ($Q=0$, $g=1$). **Resolved by ruling the same night: see
   §10.** The strong form of L8 is false; the difference is intended.
2. **The alignment record §3** should be marked as overridden outside the
   noise band by this ruling.
3. **$\kappa$ and the middle shape.** $\kappa=\nu$ is natural, not forced. The
   kink in $\min(1,\nu/Q)$ vs a smooth form with the same ends.
4. **The back-calculated-error property (§1.2)** as a named property: ours27
   meets it, ours28 misses both ways. Not ruled.
5. **The averaged true pair (§1.3)** is not codirectional. Whether any output
   should report a true pair, and which, is open.
6. **Tester presets** for §4 (these numbers are sandbox), and the 3D case.

## 8. His words, 29 Sept (verbatim)

> No more using the word push for this. Its signal.

> So we should have that, if the vectors are aligned and the uncertainty
> decreases, signal increase. If the directions are unaligned and the
> uncertainty decreases, signal decreases.

> So. Could the fact that the measure weights such highly misaligned
> readings equally be the source of the overconfidence?

> So this is a defect yes? We can start by fixing this?

> The readings are evidence. Just less.

> Yes investigate various ways. Try to make sense of it.

## 9. Checks

Claude's sandbox, numpy and scipy; scripts not kept. The whitened reading
state $(T_i,C_i)$ as in the alignment record; ours28 per-reading factor
$\log{}_1F_1(\tfrac12;1+\tfrac k2;\tfrac{H}{2})$ with its large-argument
asymptotic beyond $H/2=500$ (matched to $10^{-8}$ at the switch); law on 4001
points uniform in $\theta$ (half-Cauchy, $c=1$).

- §1.2: 1D, 40 series of 100 readings per cell; posterior-expected error
  $Q+(1-2E\beta+E\beta^2)H+E\beta$ per reading (ours27: $\beta\equiv1$),
  averaged over the law.
- §1.3: 1D and 2D, exact quadrature over the law; ours28 with $E[\beta\mid H]$
  per mass node.
- §1.4: 1D, 60 series (20 at $N=1000$), seed 7.
- §1.5: 2D, $k=2$, 8001-point grid, $\beta$ posterior on a 4000-point grid.
- §4.1: exact, no sampling. §4.2: 200 series per cell (120 at $N=100$), seed
  21, signal of fixed length along uniform random directions, true
  $\theta_0$ uniform on $[0.2,\ \pi/2-0.2]$. §4.3: same seed and design.

## 10. Sharpened later the same night (21:13–21:25)

**10.1 A side path, closed.** The principled-looking alternative to $g(Q)$ is
a per-reading shared nuisance for reliability, like $\omega_i$: a spherical
error inflation $s_i^2=1+\tau_i^2$ in whitened space, averaged over a prior.
Checked: it keeps Z exactly (zero reading, spread $0.0$), and caps the
impossible pair at a geometry-only limit ($0.693$ nats at $a=0.6$, the
coefficient forced by the model). But it tempers by the misfit **at each
mass**, $Q_i(m)$, where the ruling is about the misfit **at the best mass**,
$Q_{\min,i}$. Misfit at a mass is evidence about the mass; misfit at the best
mass is evidence about the reading; "just less" is the second. Unable to tell
them apart, the inflation model caps every reading's evidence against every
far mass, and a perfectly aligned, perfectly measured pair grows only like
$\tfrac12\log(\text{energy})$ (3.5 nats at $\sigma=0.05$ where ours28 gives
about 1000). The $g(Q_{\min})$ form leaves an aligned reading untouched by
construction. Closed.

**10.2 The collision, stated sharply.** Under ours27 the explained energy adds
over axes, so a $d$-dimensional reading is literally $d$ 1D readings sharing
$m$, and $Q_{\min}$ is exactly the disagreement between its axis readings
about the mass. The impossible pair is two 1D readings, one certain of
$m=\infty$ and one certain of $m=0$. So the tempering acts on a contradiction
packaged inside one reading and cannot act on the same contradiction spread
over two readings, where the law still collapses onto the stronger.

**10.3 His ruling: one reading is one completeness claim.** In his words:

> yes so the clear conclusion is. 1 3D readings is not the same as 3 1D
> readings? maybe the better thing to focus on is, when we make a 3D reading,
> we assume that the true vector couldn't have been in the 4th or 5th
> dimensions. And when we make a 1D reading, that same thing can be allowed.

A $k_i$-dimensional reading declares that the signal has no component outside
its declared subspace $N_i$. In Newtonian mechanics a 3D reading is complete.
A 1D reading of a 3D world is a projection: silent about the other
components, it makes no claim that can be contradicted. What the claim buys
is exactly the codirectionality conditions the chain audit counted: a
complete $k$-dimensional reading imposes $k-1$ of them, and $Q_{\min}$ is
their violation. A 1D reading imposes none and is aligned with some mass
automatically.

So "one 3D reading is not three 1D readings" means: same data, different
declaration. Three 1D readings not declared one signal are three projections,
and their disagreement is disagreement about $m$, the law's business. Three
1D readings declared projections of one signal *are* the 3D reading. Under
ours27 the definition could not read the declaration, which is why axis
splitting looked too elementary to be wrong: there it was an identity. It
stops being one as soon as anything reading-level exists ($\omega_i$ in
ours28, $Q_i$ here). The reading structure (L18, an instrument fact) is what
does the work.

**10.4 Corrections to §2 and §4.** $\nu=k_i-1$, not $d-1$, and $Q_i$ is the
residual inside the declared subspace at the best mass. A 3D reading declared
$k=1$ imposes no condition and is not tempered, which agrees with T25 (it
equals the 1D reading exactly). In the 2D, $k=2$ tests the two coincide, so
§4's numbers stand. $\nu$ now has a meaning: the number of conditions the
reading's completeness claim imposes, each costing one noise unit of residual
on average; $(k-1)/Q$ is the share of the reading's self-contradiction that
the declared noise accounts for. That is forced by meaning, not by a
property. A projection, a reading with axes declared exactly zero, and a
plain 1D reading give the same law: the law sees only the $k-1$ conditions
inside the measured subspace.

**10.5 Consequences carried to the lemma record (edited in place the same
night, on his "Go ahead").** L8 rewritten: the strong form is false, the
conditional form (components independent given the reading's shared state
$N_i,\omega_i,Q_i$) forces nothing. New premise L21 (this record's ruling)
and theorem T30 (the tempered law). O3 reopened: the flat-in-vector weighting
(T6) and the stacked length (T9) had axis splitting as their sole support and
need their own; the 25 Sept axis-splitting record's premise is withdrawn as a
forcing argument (offered, not placed there). O12: $\kappa$ and the middle
shape. O13: contradiction between readings, which the law cannot see and the
resolution reports.

## Sources

- `mass_estimator_lemmas_2026-09-29.md`: L1, L8, L13, L20, R1, T23, T24, T25,
  O6, O7.
- `alignment_what_one_reading_can_say_about_noise_2026-09-29.md` §§1, 3.
- `null_series_tilt_and_property_N_2026-09-29.md` §5.
- `ours26_learned_omega_is_a_hint_2026-09-27.md` §1 (tempering identity).
- `why_mass_cannot_be_undefined_2026-09-29.md` (never undefined).
- `mass_estimator_law_assumed_correct_steering_2026-09-07.md` (no gate).
- `ours28b/ours28b_no_borrowing_between_readings_2026-09-28.md` (L19; the
  84% case).
- His messages of 29 Sept, 20:04–20:50 (§8) and 21:12–21:25 (§10).
