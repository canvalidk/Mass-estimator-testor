# Scattered observations, 30 September morning

**2026-09-30 (morning).** Record, written by Claude at the user's request
("should we make a record the other ideas in this convo? Maybe a scattered
one that just list our observations?"). A list, not an argument. Each item
says whose it is. Numbers are from Claude's sandbox (1D, unit noise, $c=1$,
ours27 and ours28 as filed); scripts not kept.

**Relations.** Scopes `the_noise_measure_view_2026-09-29.md` (its §4 fork
and parked item 3), `null_series_tilt_and_property_N_2026-09-29.md` §5 (the
dichotomy), `alignment_what_one_reading_can_say_about_noise_2026-09-29.md`,
`ours28b/ours28b_no_borrowing_between_readings_2026-09-28.md` §§3–4, the
lemma record (L8, L16, L19, L20, T6, T9, T25, O3, O6, O13) and the 27 Sept
refusal of the sandwich (`ours26/README.md`). Binds nothing.

---

## 1. The argument as read under F1 (the session's starter, his paste)

The estimator is mass because: Newton II is three operations over the same
quantities and an entry must run; force and acceleration must be obtainable
without each other (shared-quantity argument, 1 Aug); so the mass corner
receives non-codirectional pairs and the exact rearrangement is inert on the
domain traces run over; nobody has ever held an exact mass, so there is no
meta-mass behind the readout; and because the estimator is mass it must
compose as mass composes, which forces the ratio of means. Three joints
flagged as open:

- The 1 Aug authority ruling ("outside the VD's authority, like
  differentiation") is consistent with F1 only if the estimator is as unique
  as differentiation. So uniqueness (F2) is what makes F1 and the authority
  ruling consistent, and under F1 each residual shape choice is a different
  definition of mass at weak signal.
- Persistence (O10) is unargued: every many-reading formula assumes one
  shared mass across readings; under "mass = estimator" that assumption is
  part of the definition and needs its own argument (guess: from the
  trackability/lumping chain).
- The T14 finding becomes definitional: if the starting law is part of the
  estimator it is part of the definition of mass, so a starting law carrying
  knowledge of the object would make mass depend on what was already known.
  The missing premise, "the default starting law comes from the instruments
  alone", keeps mass from being observer-relative.

## 2. $\eta$: the Bartlett power (Claude)

- Under F1 "overconfident" means: the law misdescribes how far its own point
  moves under repeats. That is property N in his words and the Bartlett
  identity (score variance $=$ curvature) in statistical ones. It needs no
  meta-mass.
- ours28's per-reading size law is not "no expectation". Independent
  $\omega_i$ from the arcsine (median 1) says about half the readings carry
  signal above one noise unit; a quiet series has prior weight of order
  $\varepsilon^N$. L19's reason (readings from any regime in any combination:
  no distribution behind the sequence of signal sizes) says $\prod_iE_i$ is a
  *working* likelihood, not a true one. The $\sqrt N$ tilt is that declared
  population reading chance correlation as signal.
- A working likelihood is honest only when raised to the power that matches
  its curvature to its score variance. Mass is a scalar, so one equation, one
  number:
  $$\eta=\min\!\Big(1,\ \frac{\mathcal H}{J}\Big),\qquad
  \mathcal H=-\sum_i\partial^2_\theta\ell_i(\hat\theta),\quad
  J=\sum_i\big(\partial_\theta\ell_i(\hat\theta)\big)^2,$$
  both at the best mass; law $\pi\,(\prod_iE_i^{g(Q_i)})^{\eta}$.
- **Z holds exactly.** A zero reading has zero score and zero curvature;
  50 appended zeros left $\eta$, law and point unchanged to $10^{-8}$. This
  separates $\eta$ from every push-scale device (ours28b and learned $\omega$
  moved 1.80 → 1.65 on one zero reading): $\eta$ responds only to readings
  disagreeing about $m$, never to how large their signals are.
- Evaluated at the maximum, $\eta$ is reparametrisation-invariant, so
  swap-invariant. The cap at 1 makes one reading give ours28 exactly
  ($J=0$) and means $\eta$ only ever widens. Strong signals have $\eta\approx1$.
- **Closed form from the carried state** (ours27, 1D): with
  $C_i=(\tilde a_i^2-\tilde F_i^2)+2i\tilde F_i\tilde a_i$, $C=\sum C_i$,
  $\hat\theta=\arg C/2$,
  $$\eta^{-1}=\frac{\sum_i|C_i|^2-\operatorname{Re}\!\big(e^{4i\hat\theta}\,\overline{\textstyle\sum_iC_i^2}\big)}{8\,|C|},$$
  matched to numerics at $10^{-8}$. The state grows from $(C,T)$ to
  $(C,T,\sum|C_i|^2,\sum C_i^2)$; combination stays exact.
- **Numbers** (coverage of the law's 95% interval, ours28 raw → tempered;
  300 series per cell; widths in $\ln m$):

  | world | ours28 | with $\eta$ | median $\eta$ |
  |---|---|---|---|
  | $N=200$, $r=1$ | 89.3 / 0.97 | 93.0 / 1.15 | 0.85 |
  | $N=100$, $r=0.5$ | 84.3 / 3.66 | 91.3 / 4.46 | 0.61 |
  | one $r=6$ + 200 nulls | 85.3 / 2.12 | 97.0 / 3.25 | 0.58 |
  | own world (arcsine), $N=20$ / $200$ | 95.0 / 96.0 | 95.7 / 96.7 | 1.00 |
  | mixed regimes $r\sim e^{U(-4,4)}$, $N=50$ | 95.0 | 95.7 | 1.00 |
  | null, $N=10$ / $100$ / $1000$ | 92 / 71 / 37 | 92 / 79 / 79 | 1.0 / .59 / .17 |

  Null tilt of the law, ours28: 1.2 / 4.5 / 12.2 nats at $N=10$ / $100$ /
  $1000$; tempered 1.2 / 2.6 / 2.0. The $\sqrt N$ is gone; a floor of about
  2 nats and 80% coverage remains.
- **It is ours25's sandwich by another route.** The 27 Sept refusal had
  three grounds: bolted on from outside (now derived from property N, L19's
  reason and the per-reading impossibility of the no-borrowing record);
  breaks axis splitting (whose strong form he ruled false 29 Sept; the
  per-reading scores are what "one reading is one completeness claim"
  licenses); the carried state cannot hold it (answered above for ours27,
  1D). Whether the derivation counts as forcing is his call.
- Not tested: 2D/3D with $g(Q)$, the readout, tester presets.

## 3. Rotation kept, axis splitting dropped (his clarification, Claude's check)

Rotation invariance (L16, T22) is untouched and $\eta$ preserves it, because
each reading's score and curvature in $m$ are frame-free once the covariance
is full. This is also why $\eta$ must use one score per reading, not one per
axis: per-axis scores would bring back ours26's frame dependence. Axis
splitting in its strong form is what changes the answer: splitting one
reading into three changes its entry in $J$ from $(\psi_1+\psi_2+\psi_3)^2$
to $\psi_1^2+\psi_2^2+\psi_3^2$. Under the rewritten L8 that difference is
intended (same data, different declaration). The weak form stays true and
forces nothing. The cost is unchanged: O3.

## 4. What is forced in the $\eta$ chain (Claude, at his question)

1. **Z rules out counting readings.** Any shared quantity must be unchanged
   by an exact zero reading, so it can only be built from per-reading terms
   that vanish on a zero reading, and $N$ cannot appear. This kills the view
   record's "$\sum H-N$" measure and settles its parked item 3: that measure
   breaks Z. Forced.
2. **Mass is a scalar, so one shared number is the ceiling**, not a
   weakening of L19. Forced as a count.
3. **L19 rules out reading signal size across readings.** Z still permits
   size statistics ($\sum T_i^2/\sum T_i$ is zero-invariant and measures how
   quiet a series is), but they infer signal size collectively, which L19
   bans. What survives is how readings bear on the shared $m$: scores and
   curvatures. **Forced only after a ruling on what "borrowing" means.**
4. **Property N picks $\mathcal H/J$** among score-and-curvature quantities,
   for long series. Forced in the limit.
5. Finishing: power form natural not forced (other finite-$N$ forms agree
   asymptotically); evaluation at the maximum forced (invariance); cap at 1
   natural (one reading gives ours28).

## 5. Null evidence does not accumulate under Z and L19 (Claude; checked pattern)

The permitted collective statistic is coherence, not size. Under pure noise
$|\sum_iC_i|^2/\sum_i|C_i|^2$ has mean $1.00$, $1.00$, $0.96$ at $N=10$,
$100$, $1000$ (2000 series each). A long null series looks no more null than
a short one. If general, exact convergence of the law to $\pi$ on null
series is unreachable inside the law under Z and L19; what is reachable is
bounded tilt (§2). "The more sure we are it is all noise, the more we know we
know nothing" would then live in the resolution, not the law. Not a proof.

## 6. What the instrument declares (his, sharpened)

His statement: force and acceleration have the same dimensions; the number
of dimensions of each vector is known; the dot product depends on the
dimension implicitly but not on how many components are zero, since zeros
are accounted for by rotation; so the missing information is vector
dimensionality only.

- Instrument facts before data: **same vector space** for $\mathbf F$ and
  $\mathbf a$ (a type fact; makes codirectionality well-typed); **$k_i$**
  (the completeness claim); **$\Sigma_i$**. Nothing about signal size or
  direction beyond $N_i$.
- **Zeros carry nothing because a zero component is not frame-invariant**:
  $(4,0,0)$ rotated has none. A law reading zeros would break L16. T25 is
  this confirmed.
- **So $k_i$ cannot be inferred from the data** (zero-counting is a frame
  artefact) and must be declared. The completeness claim is an instrument
  fact by necessity. Small forced move.
- The dot product sharpened: its value needs no $k$; whether it is large
  *for noise* does. Alignment fraction of noise: $1.00$, $0.87$, $0.71$ in
  $k=1,2,3$. Everything one reading can say about noise is a function of
  $(T_i,|C_i|,k_i)$; $\nu=k_i-1$ is a special case.
- The law never sees the world's $d$; a 3D reading declared $k=1$ equals a
  1D reading. The definition is dimension-agnostic; readings have dimensions.

## 7. His $d$-intuition, split (his statement, Claude's split)

"Two aligned vectors: in 1D the direction is certain; in 2D both could be
errored just right to align again; so confidence in direction falls with
$d$ for a given aligned pair."

- Chance alignment gets **rarer** with $d$, not commoner (automatic in 1D). A
  perfectly aligned pair is more surprising under noise as $d$ grows.
- For a **fixed energy along the line**, noise's share grows with $d$
  ($H_{\max}\approx3.0$ under noise in 2D against 1 in 1D), so the same $H$
  is weaker evidence about $m$ in higher $k$. This is the effect the law
  sees, since the law never reads $Q$ (a perfect and a mildly misaligned pair
  with equal $H$ get the same factor). ours28 has it exactly: the
  $(1-\beta)^{k/2}$ factor, the $1+k/2$ in ${}_1F_1$.
- Direction is pinned less tightly in higher $d$ for fixed size (angular
  error $\sim\sqrt{k-1}\,\sigma/\text{size}$), but direction is integrated
  out; only its effect on $m$ matters.

## 8. The shape being sought (Claude, at his question)

Given data and a candidate $(m,\mathbf a^\ast)$ the likelihood is fixed
(L13). What is sought is the measure on candidates, which by L9 is one
radial function per $k$, $\rho_k(r)$, $r$ the candidate size in noise
units; $E_i(m)=\int N(\text{data}\mid m,\mathbf a^\ast)\rho_k(|\mathbf a^\ast|)\,d\mathbf a^\ast$.
The $k$-dependence is in the equation (§7). The size part is where L19's
reason bites: the only law on $r$ that expects nothing is scale-free, which
is improper at the noise end with an $m$-free divergence, so every reading
contributes nothing. This is the impasse *inside the law*; see §10.

## 9. A candidate to half-close O3 (Claude; unchecked)

Herschel–Maxwell conditional on the reading's scale: given $\omega_i$,
components independent and isotropic $\Rightarrow$ Gaussian with scale
$\omega_i$. Averaging over $G$ gives ours28's family: $\rho_k$ is a scale
mixture of Gaussians, $G$ free. So the conditional L8 forces the *class*
(T6 as a family) and leaves the mixing law, which is declared anyway. T9
still owed. Not checked against the lemma record's "forces nothing".

## 10. The reframe (his correction, Claude's identification)

His words: "if the truth isn't where you looked, you looked in the wrong
place as a definitional axiom. You're trying to force this ideal, it
doesn't work, but the correct answer still exists. So the only remaining
conclusion is, you were wrong that the solution has this idyllic property."

The property identified: **the law carries its own uncertainty.** Every
"impossible" of the morning (no per-reading rule is calibrated; the
$\sqrt N$ cannot be removed inside L19; null convergence is out of reach
under Z and L19) is a true theorem about the law's width, and a wall only if
that width is the result's uncertainty. The records already gave that up:
(C) broke in the `==^` thread; the output is law, point and resolution; the
view record §4 named "in the trust, not the law" and it was not taken.

Consequences, offered:

- The $\sqrt N$ tilt on a null series is a preference; the resolution
  $(0,\infty)$ beside it says what the preference is worth (property N).
  Under F1 that is what mass is on that data.
- The noise measure ("how trustworthy among other data") is a combination
  question; combination acts on resolution and contributions, never on the
  law's width. So it lives beside the law. It is collective, a function of
  the series, and that is allowed: L19 governs how readings join *in the
  law*. The null-series dichotomy was false: sharing in the trust layer is
  not borrowing.
- The coverage tests of the past week measured the law's interval against
  the truth. Under this reading they measured the wrong object; the standard
  becomes the resolution's coverage and the point's behaviour. **A change to
  what the tester counts as a failure; his call.**
- Candidates for the noise measure beside the law: the resolution set;
  $|C|/(2\sqrt N)$ against its Rayleigh null; $\eta=\mathcal H/J$, now read
  as trust in the law's width rather than a repair to the law.
- Next step proposed: the 1 kg vs 1000 kg combination case, each candidate
  as the quoted uncertainty, to see which makes inverse-variance weighting
  give about 1000 kg.

## 11. His words, 30 Sept (verbatim)

> So to clarify, we have still have ratatiol invariance, what we are
> rejecting is axis splitting.

> Can you see natural and forced moves beyond what is written in the notes
> so far? is etta that?

> We were gonna start by establishing what we do know from the instrument.
> And I think the answer is simple. First, the force and acceleration
> equations have the same dimensions. We already know the number of
> dimensions of each vector too. In the same way that we don't really
> mention the dot product depends on the number of dimensions on the vector
> right? But then it doesn't depend on the number of dimensions which are
> zero. Since those already are accounted for by rotational invariance?
> Which means the missing info is just vector dimensionality only?

> It feels like you keep hitting your head on the same thing. You keep
> looking for some 'perfect' solution. And getting demotivated when you
> can't find it. But you aren't allowed to believe it can't be done. If it
> can't be done, just allow yourself to fail today.

## 12. Checks

Sandbox, numpy and scipy. Per-reading log factor $\log{}_1F_1(\tfrac12;\tfrac32;H/2)$
via Dawson's function (ours28, $k=1$), $H/2$ for ours27; law on 2000 points
uniform in $\theta\in(0,\pi/2)$ (half-Cauchy, $c=1$); $\hat\theta$ by
Brent on the full circle then bracketed; scores and curvatures by central
differences, $h=10^{-4}$.

- §2 null table: 200 series per $N$ (60 at $N=1000$), seed 7. Coverage
  table: 300 series per world, seed 21, true $\theta_0\sim U(0.2,\pi/2-0.2)$,
  signal of the stated length along random signs; own-world cell rechecked
  on a 20001-point local grid (96.0 / 96.7). Null coverage: $\theta_0$ drawn
  from $\pi$, 300 series (120 at $N=1000$), seed 5. Z check: 30 readings
  with a 0.3 offset plus 50 exact zeros, seed 3. Closed form: two series,
  seed 1, agreement $2\times10^{-8}$ and $4\times10^{-9}$.
- §5: 2000 series per $N$, seed 0.
- §7 alignment fractions and $H_{\max}$: from the alignment record §§2, 4.

## Sources

- The session's starter (his paste, 30 Sept, 10:09).
- `the_noise_measure_view_2026-09-29.md` §§1.1, 2, 4, 5, 6.
- `null_series_tilt_and_property_N_2026-09-29.md` §§1, 5, 6.
- `alignment_what_one_reading_can_say_about_noise_2026-09-29.md` §§2, 4, 5.
- `readings_are_evidence_just_less_2026-09-29.md` §10.
- `ours28b/ours28b_no_borrowing_between_readings_2026-09-28.md` §§1, 3, 4.
- `ours28b/ours27_one_scale_per_noise_pair_2026-09-28.md` §1 (Bartlett).
- `ours26_learned_omega_is_a_hint_2026-09-27.md` §3.
- `../ours26/README.md` (the 27 Sept refusal of the sandwich).
- `../the_combination_equation_is_exact_2026-09-26.md` §1 (the state).
- `../mass_estimator_lemmas_2026-09-29.md`: F1, F2, L8, L9, L13, L16, L19,
  L20, T6, T9, T22, T25, T27, O3, O6, O10, O13.
- His messages of 30 Sept, 10:09–11:04.
