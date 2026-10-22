# The noise-measure view: one missing property behind the open items

**2026-09-29 (evening).** Record, written by Claude at the user's request
("Could you write a markdown on this way of looking at the problem"). The
view is his, reached in conversation the same afternoon; his words are in
§7. The reconstruction, the numbers and the map are Claude's.

**What it says.** The equation has no notion of how much of a reading is
push and how much is noise. Every open item in the lemma record, from the
push scale to the readout to the centre of $\pi$, is a place where that
missing notion shows through. So the open items may not be several
freedoms to close one by one; they may be one property, not yet stated,
whose absence leaves all of them loose.

**Relations.** It scopes `mass_estimator_lemmas_2026-09-29.md` (its §4 and
§5), `why_mass_cannot_be_undefined_2026-09-29.md` (a correction to it is
owed, §6), the ours28b naming record, and the coverage failures in
`ours26/README.md`. It binds nothing. It is a way of looking, like
`mass_estimator_default_weighting_direction_2026-09-23.md`, not a result.

---

## 1. How it was reached

The thread started at `==^`, the compatibility relation physicists use to
say whether two measurements disagree. His stance: `==^` predates us; we
named it, we cannot redefine it, and we have to trust that the practice
does not lead to contradiction.

**1.1 The contradiction.** Four things could not all stand:

- (A) at $0/0$ the output is a point $c$ and a proper law $\pi$ with a
  declared centre (his ruling of 29 Sept);
- (B) F1 read literally: the output law *is* the object's mass;
- (C) the law's width is the result's uncertainty under `==^`;
- (D) an independent route (density, $\rho V$) can give a mass far from $c$,
  and Newtonian mechanics has no contradictions.

With $c=1\,$kg and a density result of 1000 kg, (C) makes `==^` fail. So one
of the four gives. His resolution, in his words: the combination rule is
already fixed, it is inverse-variance weighting, and for the combined answer
to be about 1000 kg the $0/0$ result must enter with weight zero. **The
uncertainty quoted for the $0/0$ result is "unresolved", not the width of
the law around $c$.** (C) is what gives.

Underneath: a *result* is $\pi\times$data; a *contribution* is the data
factor alone. Inverse-variance weighting is the Gaussian case of the product
rule, and the product counts $\pi$ once (T7). The $0/0$ contribution is
exactly 1. Combining two *results* would count $\pi$ twice. Even done
correctly $\pi$ is still there once: half-Cauchy at $c=1\,$kg times
$N(1000,100)$ peaks at 980 kg, a 2% tilt from $\pi$ alone, inside the
100 kg width. That is the price of a declared centre.

**1.2 His question: does $0/0$ really carry no information?** A real
instrument with $\sigma=0.1$ almost never prints $0/0$; it prints
$0.1,\,-0.1,\,0.2,\ldots$ Those should be the readings that carry little.
Three cases:

| case | what the law does | status |
|---|---|---|
| exact $0/0$ | data factor $\equiv1$: nothing | theorem (T24) |
| one noise-sized reading, e.g. $(0.1,-0.1)$ at $\sigma=1$ | tilts by 0.5%; $O(\varepsilon^2)$ | continuous with the exact case |
| **a series of $N$ noise-sized readings** | **tilts by $O(\sqrt N)$ in log, in a direction the noise picked** | **the current equation does not converge to "we know nothing"** |

Simulated (ours27, 1D, unit noise, true pushes zero, 400 series each):

| $N$ | tilt of the law, $\max_m-\min_m$ of $\log$ data factor, median | 90th percentile |
|---|---|---|
| 1 | 0.5 | 1.8 |
| 10 | 2.8 | 5.6 |
| 100 | 9.2 | 17 |
| 1000 | 29 | 56 |

The mechanism is the cross term $2m\sum_i\tilde F_i\tilde a_i/(m^2+1)$ in
$\sum_iH_i(m)$. Force noise and acceleration noise are independent, but
their sample correlation over $N$ draws is of size $\sqrt N$, never zero,
and the equation reads it as evidence that the pushes lie along some slope.
ours28 tempers this by about a factor of 3 in the log (its per-reading
evidence is ${}_1F_1$, not exponential) and does not remove the growth.

**1.3 His proposal, property N.** The result on a null series must remain
a true statement: if the point wanders, the quoted uncertainty must grow
with it. Checked on the same series with the data-only compatibility set
$\{m:\sum_iQ_i(m)\le\chi^2_{N,0.95}\}$:

| $N$ | series that excluded any mass at 95% | median fraction of the mass range excluded |
|---|---|---|
| 10 | 11% | 0 |
| 100 | 15% | 0 |
| 1000 | 14% | 0 |

The misfit $\sum Q(m)$ sits at $N\pm O(\sqrt N)$ for every $m$ and the
threshold's margin is also $O(\sqrt N)$, so the set stays essentially the
whole half-line at every $N$. **The resolution grows in step with the
wandering point.** On a null series the law *prefers* (a chance
correlation), the data *exclude* nothing, and "point $\widehat m$,
resolution $(0,\infty)$" is true at every $N$. The $0/0$ case is the limit
where even the preference vanishes.

**1.4 His diagnosis.** What is missing is a noise measure: a function that
says, of the aligned pair the equation can always find, how likely it is to
be representative of a push rather than of noise, so that its evidence is
weighted accordingly. Its absence is why the equation is overconfident and
why ours28b, which estimates a signal-to-noise ratio, is not.

## 2. The noise measure, in the equation's own terms

Each reading contributes an explained energy $H_i(m)$: the part of the
whitened reading lying along the candidate slope $m$. Under pure noise
$\sum_iH_i(m)\approx N$ for every $m$ (one unit of noise energy per reading,
whatever the line), with chance fluctuations of size $\sqrt{2N}$. Under a
real signal of total push $R^{\ast2}$ in noise units,
$\sum_iH_i(m_{\rm true})\approx N+R^{\ast2}$. So

$$
\text{signal excess}=\sum_iH_i(\widehat m)-N,
\qquad
\text{chance level}=\sqrt{2N},
$$

and the noise measure is the excess against the chance level. An
alignment whose excess is below the chance level is what noise produces on
its own. The $\sqrt N$ tilt of §1.2 is this quantity, read as a defect to
tolerate instead of a number to measure.

**It is cross-reading by construction.** In 1D a single reading
$(\tilde F,\tilde a)$ is always exactly aligned with some mass,
$m=\tilde F/\tilde a$. One reading cannot tell a push from noise that
happens to point somewhere. The only evidence that a series is noise is
that its readings disagree about $m$. So any noise measure shares
something across readings.

**The three equations in this light.**

| | push scale | what the equation can learn about noise |
|---|---|---|
| ours27 | flat ($\omega\to\infty$) | nothing: every push is infinitely large in noise units, the overconfidence in one line |
| ours28 | per reading, arcsine | nothing across readings: each reading's scale is its own nuisance; tempers by ~3, growth stays $\sqrt N$ |
| ours28b | one scale per noise pair, shared | the series' signal-to-noise; withdrawn under L19 |

## 3. The map: what the noise measure touches

| freedom / open item | how it touches it |
|---|---|
| push scale (L10′, $G$, O3) | it *is* a push-scale statement |
| no borrowing (L19) | forbids the only route a noise measure can take |
| $\pi$'s centre $c$ (L3′, O11) | $c$ is the answer at zero signal-to-noise; the knowledge-free default $c=\sigma_F/\sigma_a$ is the noise measure evaluated at zero signal. One object, two faces. |
| $\pi$'s shape and $G$'s shape (O6) | "uniform in the angle the noise defines" is already noise-referenced; a noise measure adds a third shape (the law on shared energy) or, if forced, pins the others |
| the readout (T29, O1, O2) | the readout fails Z because a candidate cloud has noise-sized size with no data; size $=$ explained energy minus its noise expectation passes Z. An O2 candidate. |
| declared subspace $N_i$, $k$ (O7) | the null level of $\sum H$ depends on $k$; T25's "$k$ matters" is the noise measure showing through |
| axis splitting restated (L8) | a reading's "shared state" would include its noise budget |
| coverage (76.5% at 21 weak readings, `ours26/README.md`) | the measured symptom |

## 4. Two readings of "touches everything"

**It is the missing foundation.** The freedoms were left open because
nothing forced them, and nothing forced them because the equation had no
notion of how much of a reading is push. State that as a property and
several may close together: the push scale, the readout repair, what $c$
means. The map is some evidence: the same quantity, $\sum H-N$, appears in
most rows.

**It is one more free shape.** A shared energy scalar needs a law of its
own, what the series' signal energy was expected to be before data, and
that law is declared. Then several small freedoms have been traded for one
larger one.

What decides between them is whether a *property* forces the noise
measure's form, as axis splitting forced the flat-in-vector weighting.
Candidates:

- **Null convergence.** A series of instrument noise converges to the
  starting law. His sentence; the current equation violates it (§1.2).
- **Z extended.** A reading indistinguishable from noise tends to add
  nothing. Z as ruled covers exact nulls only.

If one of these forces a unique noise measure, it is the foundation. If it
leaves a shape free, it is a trade.

**The weakening of L19 that keeps its point.** One scalar may be shared,
the series' total signal energy; nothing else. No direction, no per-reading
size, no spread. A reading is never told where its push points, only
whether the series as a whole has pushed. That is narrower than ours28b's
per-noise-pair scale.

**Where it lives is a separate choice.** In the law: the shared scale
enters the combination, the law's width widens on quiet series, the
coverage failure should close (the ours28b route). In the trust, not the
law: L19 stays intact for the law, the noise measure sits in the output
beside the resolution and governs the result's weight in combination; the
law's own intervals stay overconfident but nothing downstream is misled.

## 5. What holds regardless of the choice

- **The output has a third part.** Law, point, and *resolution*: the
  data-only compatibility set. Case II §9 of the 5 Sept two-cases record had
  this five-line output; property N makes it a requirement.
- **The resolution is computed from the data alone.** $\pi$ and $G$ stay
  out of it. That is what makes the result a true statement whatever is
  declared.
- **`==^` and combination act on contributions and resolution, never on
  the law's width.** A null series carries weight zero in any combination.
- **Property N is independent of the 7 Sept ruling.** The resolution is
  descriptive. An empty set means "no mass fits at level $\alpha$", read
  under 7 Sept as improbable error, never rejection. Nothing halts.

## 6. Parked (from the conversation; delivered, not dropped)

1. Whether a single shared energy scalar breaks the "F = 1, then 50000000,
   then 0.00001" case that withdrew ours28b, or whether the per-pair scale
   was what broke. Needs the no-borrowing record reread and a preset.
2. Under 7 Sept the noise measure must not become a gate: "excess below
   chance level" means weight lightly, never reject.
3. Z and the noise measure: an exact null reading adds 0 to $\sum H$ and 1
   to $N$, so it lowers the estimated excess. Either right (a null reading is
   evidence the series is quiet) or a breach of Z's intent. Ruling needed.
4. Whether the $\sqrt N$ *preference* on a null series is acceptable once
   the result is true (property N), or whether a further property should
   suppress it.
5. **Correction owed to `why_mass_cannot_be_undefined_2026-09-29.md`.** Its
   A5 (the cost of making "undefined" pass `==^`) is overstated: the real
   `==^` already has a slot, "unconstrained", which passes against
   everything at no cost and propagates as interval arithmetic. Its §2.3
   called the identity reading "a renaming that buys nothing"; against
   practice that is unfair. What still rules out "undefined" is F1, Z and
   his ruling, not `==^`. The record's §3 (never quote $\pi$'s width as the
   result's uncertainty) is confirmed, and is now forced by practice rather
   than chosen. Offered; not made.

## 7. His words, 29 Sept (verbatim)

On combination:

> Clearly, if we combine the two pieces of information: 0/0 from the mass
> measurement and 1000kg from rho*volume, we should have the answer be
> around about 1000. But the method we use to combine them is picked
> already. It's how sigma's are combined using inverse weighting. Then for
> those to hold, the uncertainty quoted by for the 1kg result needs to be
> that inverse weighting making everything consistent and not the width of
> the law of 1kg?

On whether $0/0$ carries no information:

> If we had, a stationary object and really have uncertainty be say 0.1,
> what are the odds that we actually get 0/0 form that? I mean isn't it
> more likely that the information makeup is more like 0.1, -0.1, 0.2...
> Then THAT should be that thing that conveys little info. The more and
> more we are sure that in all these cases all the data comes out of the
> noise, then we become sure that we don't know anything about the actual
> mass of the thing.

On property N:

> If the uncertainty also grows with the point estimate, we might have a
> result that is a true statement.

On the noise measure:

> We need to come up with an understanding of a noise measure. In other
> words some function that allows us to say, yes although we can find some
> aligned pair for the data, it is unlikely for the found pair to be
> representative of the signal and thus we need to recognize that although
> this measure is possible, when considering how trustworthy it is among
> other data, we need to recognize it is likely unreliable evidence. If
> that doesn't exist in the equation, that would explain why we are
> overconfident and why the ours28b equation actually isn't. Because it
> solves the issue by estimating the noise to signal ratio.

> So we really do have a problem that touches all the freedoms we identify
> in the equation.

## 8. Checks

All in Claude's sandbox, scipy and numpy, 1D, unit noise, $c=1$, ours27
unless stated. Scripts not kept.

- **Tilt on null series (§1.2):** for each $N\in\{1,10,100,1000\}$, 400
  series of $N$ readings with $\tilde F_i,\tilde a_i\sim N(0,1)$
  independent; $\log$ data factor $=\tfrac12\sum_i(m\tilde F_i+\tilde a_i)^2/(m^2+1)$
  on a 401-point log-spaced mass grid over $[e^{-4},e^{4}]$; reported
  $\max_m-\min_m$, median and 90th percentile over series.
- **Single reading:** $(0.1,-0.1)$ gives a range of 0.0048.
- **Resolution on null series (§1.3):** same series;
  $\sum_iQ_i(m)=\sum_i(\tilde F_i^2+\tilde a_i^2)-\sum_iH_i(m)$ (exact in
  1D); threshold $\chi^2_{N,0.95}$; 601-point grid over $[e^{-6},e^{6}]$.
- **$\pi$'s residue (§1.1):** $\arg\max_m\big[-(m-1000)^2/(2\cdot100^2)-\ln(1+m^2)\big]=979.6$.
- **ours28 tempering:** ${}_1F_1(\tfrac12;\tfrac32;x)\approx1+x/3$ against
  $e^{x}\approx1+x$ at small $x$, so the log-slope is one third.

## Sources

- `mass_estimator_lemmas_2026-09-29.md`: F1, L3′, L8, L10′, L19, L20, T7,
  T12, T24, T25, T27, T29, O1–O3, O6, O7, O11.
- `why_mass_cannot_be_undefined_2026-09-29.md` §§2.3, 3, 4.
- `mass_estimation_zero_zero_two_cases.md`, Case II §9.
- `mass_estimator_law_assumed_correct_steering_2026-09-07.md` (the 7 Sept
  ruling).
- `ours28b/` naming record and `ours28b_no_borrowing_between_readings_2026-09-28.md`
  (L19; ours28b withdrawn).
- `ours26/README.md` (the coverage failures).
- His messages of 29 Sept, 14:44–17:06 (§7).
