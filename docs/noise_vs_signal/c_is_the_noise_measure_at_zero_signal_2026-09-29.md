# $c$ is the noise measure at zero signal, and the noise measure does not fix $c$

**2026-09-29 (night).** Record, written by Claude at the user's request
("Ok now do number 3"), developing §6 of
`null_series_tilt_and_property_N_2026-09-29.md`. Short by design.

**What it records.** The centre $c$ of the starting law is what the
definition returns when the noise measure reads zero: it is the noise
measure's own zero. But the noise measure cannot supply $c$. Under pure
noise the law's preferred mass is spread over the whole half-line with no
concentration at any instrument scale, and Z forbids $c$ from depending on
which instruments a series happens to contain. So $c$ is declared before
data, and L3′ ("from the instruments' SDs") is a convention for what to
declare, never a derivation. This settles the direction of O11 without
settling its value.

**Relations.** Scopes the lemma record (L3′, L20, T14, T20, O11) and the
two records in this folder. Binds nothing.

---

## 1. The requirement

Whatever noise measure gets built, as it goes to zero the point must go
to $c$ continuously. With one noise pair the cloud readout does this
(T14: zero readings give $s$, when $c=s$), and the best-guess readout does
not (undefined record §4: direction-dependent limits). Continuity at zero
signal is a test any candidate readout has to pass, and it is the same
test as Z taken to its limit.

## 2. Noise picks no scale

**One noise pair.** Units (L3) and the swap (L4, T5) force the zero-signal
point to be $c$, and the only scale in play is $s=\sigma_F/\sigma_a$. So
$c=s$ is forced *if $c$ is to come from the instruments at all*. That is
T14.

**Several noise pairs.** Two instruments with $s_1=1$ kg and $s_2=100$ kg,
fifty pure-noise readings each, 2000 series, ours27: the law's preferred
mass falls near $s_1$ in 19% of series, near $s_2$ in 22%, between them in
20%, and at $m\to0$ or $m\to\infty$ in 39%. The direction $\varphi_0$ of
the carried state under noise is uniform, so the law prefers *something*
in every series and it is never the same thing. The noise measure reads
zero on average and carries no scale.

Units and swap alone allow any $c$ with $c(1/\mathbf s)=1/c(\mathbf s)$:
every weighted geometric mean of the $s_g$ qualifies. Nothing in the noise
selects among them.

## 3. Z closes it

L20: a reading that adds no information does not change the answer. Add
one null reading from a new instrument with a new $s$. If $c$ were any
function of the series' instruments, the answer at zero signal would move.
So $c$ cannot depend on the instrument set of the series. That is T20 (the
mixed-instrument no-go) reached from the other side: **$c$ is declared
before data, and it is not a function of the data's instruments.**

## 4. What this does to L3′ and O11

- L3′ as written, "$\pi$'s scale comes from the instruments' SDs alone", is
  not a derivable property. It is a *convention for the declaration*: with
  one instrument known in advance, declare $c=s$. With several, declare a
  $c$; the geometric mean is a convention that respects the swap, and no
  more than that.
- O11 splits. *Where $c$ comes from* is settled: declared, before data,
  outside the noise measure. *What to declare* stays open, and is the same
  question as O6 (the shapes): a parameter of the declaration, to be
  argued as forced or recorded as an alternative under F2.
- The 2% residue at 1000 kg (view record §1.1) is the cost of the
  declaration, and it is the only footprint $c$ leaves once data arrive.

## 5. Checks

Sandbox, numpy, ours27, 1D, seed 5. Per instrument $g$: $\sigma_F=s_g$,
$\sigma_a=1$; whitened explained energy summed over both; peak of the log
data factor on 1601 log-spaced masses over $[e^{-8},e^{8}]$; edge counted
at $|\ln m|>7.9$. Histogram over 2000 series: 416 / 378 / 405 / 438 / 363
for [low edge, near $s_1$, between, near $s_2$, high edge], "near" meaning
within $e^{\pm2}$. Script not kept.

## Sources

- `null_series_tilt_and_property_N_2026-09-29.md` §§1, 6.
- `the_noise_measure_view_2026-09-29.md` §§1.1, 3.
- `why_mass_cannot_be_undefined_2026-09-29.md` §4.
- `mass_estimator_lemmas_2026-09-29.md`: L3, L3′, L4, L20, T5, T14, T20,
  O6, O11.
