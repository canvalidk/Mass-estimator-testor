# The null-series tilt, and property N

**2026-09-29 (evening).** Record, written by Claude at the user's request
("Okay do number 1"), developing §1.2–1.3 of
`the_noise_measure_view_2026-09-29.md`. The question and property N are
his; the derivation, the numbers and the structural conclusion in §5 are
Claude's.

**What it records.** A series of $N$ readings that are pure instrument
noise (true pushes zero) does not leave the law of mass alone. Under
ours27 the law tilts by about $\sqrt N$ in log; under ours28 by about half
that; the point moves by a factor of 6 from $c$ at $N=100$, and the law's
95% width shrinks below the starting law's. The tilt has an exact form: it
is the carried state $C$ of the combination record, and its size under
noise is Rayleigh with scale $2\sqrt N$. The data-only resolution stays the
whole half-line at every $N$, so the result is a true statement (property
N). A fix at the level of the law has to be collective; no per-reading
tempering can remove the $\sqrt N$.

**Relations.** Scopes the view record, `the_combination_equation_is_exact_2026-09-26.md`
(whose state $C$ is the tilt), `ours26/ours26_learned_omega_is_a_hint_2026-09-27.md`
(whose tempering identity is the mechanism a fix would use), the lemma
record (L13, L19, L20, T7, T24, T27, O2) and `ours26/README.md` (the
coverage failures). Binds nothing.

---

## 1. The tilt has an exact form

1D, unit noise, ours27. The log of the data factor at mass $m=\tan\theta$,
$\theta\in(0,\pi/2)$, is

$$
L(\theta)=\tfrac12\sum_i\frac{(m\tilde F_i+\tilde a_i)^2}{m^2+1}
=\tfrac14 T+\tfrac14\,|C|\cos(2\theta-\varphi_0),
$$

with $T=\sum_i(\tilde F_i^2+\tilde a_i^2)$ and
$C=\big(\sum_i\tilde a_i^2-\sum_i\tilde F_i^2\big)+2i\sum_i\tilde F_i\tilde a_i$,
$\varphi_0=\arg C$. These are the $(C,T)$ of the carried state in the
combination record: **the law's whole dependence on the data is one
complex number, and the tilt is its modulus.** The constant $T/4$ is
$m$-free; subtracting the expected noise energy changes nothing in the
law. What tilts it is $|C|$.

Two things follow at once.

- **Positivity halves the tilt.** $m>0$ means $2\theta$ runs over half a
  circle, so the range of the cosine is $\tfrac14|C|(1+|\cos\varphi_0|)$,
  not $\tfrac12|C|$. Mean factor $1+2/\pi$ over $2$, about $0.82$.
- **Under pure noise $|C|$ is Rayleigh.** Each of $\sum\tilde a^2-\sum\tilde F^2$
  and $2\sum\tilde F\tilde a$ has variance $4N$, so $|C|/(2\sqrt N)$ is
  Rayleigh with unit scale, median $\sqrt{\ln4}=1.18$. Chance correlation
  between force noise and acceleration noise is the whole effect.

Checked: grid and formula agree to $10^{-2}$ at $N=1000$ (grid
resolution); the median of $|C|/(2\sqrt N)$ is $1.21$–$1.23$ at
$N\ge100$; $T/(2N)$ is $1.00$.

**The noise measure, in this form.** A signal of total push $R^{\ast2}$
adds $R^{\ast2}$ to $|C|$ (and to $T$). So $|C|/(2\sqrt N)$ against its
Rayleigh null is exactly "is the alignment we found bigger than chance".
It is a function of the carried state, so it combines exactly.

## 2. Measured, ours27 against ours28

1D, unit noise, $c=1$, readings $\tilde F_i,\tilde a_i\sim N(0,1)$
independent. ours28 is the reference implementation (arcsine size law,
$k=1$).

**The law's tilt** ($\max_m-\min_m$ of log law minus log $\pi$, median over
200 series):

| $N$ | ours27 | ours28 | ratio |
|---|---|---|---|
| 1 | 0.49 | 0.18 | 0.37 |
| 10 | 2.77 | 1.25 | 0.45 |
| 100 | 9.03 | 4.49 | 0.50 |

ours28 tempers by a factor of 2–3 and the factor weakens as $N$ grows
(readings with larger $H_i$ are tempered less). Growth from $N=10$ to
$100$: ours27 $\times3.3$, ours28 $\times3.6$. Both are $\sqrt N$ or a
little faster.

**The point and the width** (median over 30 series at $N=10$, 12 at
$N=100$; the readout is T10 as implemented):

| $N$ | equation | median $\lvert\ln\widehat m\rvert$ | 90th pct | median 95% width of the law (log) |
|---|---|---|---|---|
| 0 | either | 0 | 0 | 6.47 ($\pi$ alone) |
| 10 | ours27 | 0.75 (factor 2.1) | 1.47 (4.3) | 5.11 |
| 10 | ours28 | 0.25 (1.3) | 0.87 (2.4) | 6.08 |
| 100 | ours27 | 1.88 (6.5) | 2.65 (14) | 3.72 |
| 100 | ours28 | 0.97 (2.6) | 1.55 (4.7) | 4.35 |

A hundred readings of instrument noise move ours28's point by a typical
factor of 2.6 from $c$ and narrow its 95% interval from 6.5 to 4.4 in log,
in a direction the noise chose. The definition has learned something from
nothing.

## 3. Dimension does not help the law

In $d$ dimensions with identity noise the form of §1 is unchanged, with
$\sum\tilde F_i\cdot\tilde a_i$ for the cross term and each sum
$\chi^2_{dN}$. So the tilt grows like $\sqrt{dN}$: measured $1.34\sqrt N$
in 2D against $0.95\sqrt N$ in 1D. Each reading's perpendicular component
(the $\chi^2_{d-1}$ residual of the old gate) sits in $Q_i$, not $H_i$,
and the law never sees $Q$. The information that would let a single
reading call itself noise is present in $d\ge2$ and discarded by the law.
That is a per-reading fact and is developed in the companion record on
alignment; it does not change §5.

## 4. Property N holds: the result is true

His property: on a null series the result must stay a true statement,
the quoted uncertainty growing with the wandering point. The data-only
compatibility set $\{m:\sum_iQ_i(m)\le\chi^2_{N,0.95}\}$ on the same
series (view record §1.3): 11–15% of series exclude any mass, and the
median fraction of the mass range excluded is 0, at $N=10$, $100$ and
$1000$. $\sum Q(m)=T-\sum H(m)$ sits at $N\pm O(\sqrt N)$ for every $m$
and the threshold's margin is also $O(\sqrt N)$, so the set does not
tighten. **The law prefers; the data exclude nothing; "point
$\widehat m$, resolution $(0,\infty)$" is true at every $N$.**

So the tilt is not a false statement. It is a preference with no data
behind it, and the output type that carries the resolution beside the law
(view record §5) is what keeps it honest.

## 5. What a fix must be

- **Per-reading tempering changes the coefficient, never the power.**
  ours28 tempers each reading's factor on its own. But each reading of a
  null series looks the same to it (noise-sized, $|y_i|^2\approx2$), so it
  tempers each by about the same amount, and the sum of $N$ tempered
  chance tilts is still $\sqrt N$. That is the table in §2: the ratio sits
  near $\tfrac12$ and the growth is unchanged. No law of the form
  $\prod_i g(H_i)$, with $g$ fixed, can remove the $\sqrt N$, because
  $\sum_i\log g(H_i(m))$ has the same chance-correlation structure as
  $\sum_iH_i(m)$.
- **A fix that works is collective.** The tilt is $|C|$; the chance level
  is $2\sqrt N$; a law that is tempered by their ratio (the tempering
  identity of the learned-ω record: raise the likelihood to a power
  $\beta$, which leaves the peak and widens the law by $1/\beta$) can make
  a null series converge to $\pi$. Any such $\beta$ is a function of the
  series, which is what L19 forbids. **So the dichotomy is exact: keep L19
  and keep the $\sqrt N$, or share one scalar and lose it.** The narrowest
  scalar is $|C|/(2\sqrt N)$, or equivalently the signal excess against
  $T$; it carries no direction and no per-reading size.
- **Which world this is in.** For ours27 the flat push says every push is
  infinitely large in noise units, so a null series is outside its declared
  world: a tier-3 test. For ours28 the arcsine size law gives small pushes
  positive density, so a null series is inside its world, in a
  low-probability region. T27's calibration is an average over the world
  and can hold while this region is overconfident. Whether that counts as a
  tier-2 failure depends on whether calibration is demanded conditionally
  on the series. Not ruled.

## 6. The zero-signal limit is $c$

As the readings shrink along any fixed direction, $|C|\to0$, the law tends
to $\pi$ and the point to $c$ (the undefined record §4 has the
direction-dependent limits of the *best-guess* readout; T10's cloud
readout goes to $c$ when $c=s$). Whatever noise measure is built must
reduce to $c$ continuously as $|C|/(2\sqrt N)\to0$. That is the sense in
which $c$ is the noise measure evaluated at zero signal, and the sense in
which L3′ ($c$ from the instruments) is a statement about the same object.

## 7. Open

- **O-a.** Rerun §2 in the tester as presets (this record's numbers are
  from Claude's sandbox and the reference implementation). Add a 2D case.
- **O-b.** The conditional-calibration question of §5: is a definition
  allowed to be overconfident on a region of its own world?
- **O-c.** Which scalar: $|C|/(2\sqrt N)$, or $(\sum H(\widehat m)-N)/\sqrt{2N}$,
  or the learned $\beta$ of the ω records. They are the same object at
  leading order; the differences are what a property would have to decide.
- **O-d.** Under 7 Sept the scalar tempers, never gates. A form that
  vanishes smoothly at zero excess and saturates at 1 for strong signal is
  the shape to look for.

## 8. Checks

Sandbox, numpy and scipy; ours28 via `ours27/ours28_reference.py`
unmodified (`flat=True` for ours27), `nodes` 241–361, `lm_half` 5–6,
`tnodes=121`. Scripts not kept.

- §1: 300 series per $N$; grid $\max-\min$ of $L$ on 1801 log-spaced masses
  over $[e^{-9},e^{9}]$ against the closed form, max discrepancy
  $9.5\times10^{-4}$ ($N=1$) to $1.5\times10^{-2}$ ($N=1000$).
- §2 tilt: 200 series per $N$, seed 7. §2 point and width: seed 11.
- §3: 300 series, 2D, unit noise on both axes.
- §4: from the view record (seed 2, 400 series).
- §6: undefined record §4 and §8.

## Sources

- `the_noise_measure_view_2026-09-29.md` §§1.2–1.3, 2.
- `the_combination_equation_is_exact_2026-09-26.md` §1 (the state $C,T,M$).
- `ours26/ours26_learned_omega_is_a_hint_2026-09-27.md` §1 (tempering identity).
- `ours27/ours28_reference.py`; `ours27/ours28_reference_implementation_and_readout_2026-09-28.md`.
- `mass_estimator_lemmas_2026-09-29.md`: L13, L19, L20, T7, T10, T24, T27,
  O2; §6 (tiers).
- `why_mass_cannot_be_undefined_2026-09-29.md` §4.
