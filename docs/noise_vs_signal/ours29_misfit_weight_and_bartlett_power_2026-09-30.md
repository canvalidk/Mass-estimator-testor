# ours29: the flat law, each reading's own misfit, and the Bartlett power

**2026-09-30 (afternoon).** Record, written by Claude after his three rulings
of 11:57 ("1. Yes good. 2. Good. 3. Yes. but we can report both for now. Go
ahead."). The rulings, the reason for L19 and the correction that led here
are his, quoted in §12. The equation's assembly, the code, the derivations
and the numbers are Claude's. It continues
`scattered_observations_2026-09-30.md`, which covers the morning up to 11:04.

**What it records.** The working equation moves from ours28 to **ours29**:
each reading contributes the flat (ours27) factor, weighted by its own
off-line misfit (L21, kink form), and the product is raised to the
**Bartlett power** $\eta=\min(1,\mathcal H/J)$, computed from the readings'
own residuals at the series' best line. ours28's per-reading size law $G$ is
dropped. Z forces two details the earlier sandwich (ours25) got wrong: $J$
is the raw sum (no Bessel factor, which counts readings), and the best line
is sought on the whole projective line of slopes (not on $m>0$ only). The
record also carries the afternoon's other results: the combination case, the
$\lambda$ move applied to the size law, the eigenvalue dictionary, and the
reason for L19 as he stated it.

**Relations.** Continues `scattered_observations_2026-09-30.md` (§§2, 10).
Scopes `mass_estimator_lemmas_2026-09-29.md` (proposed changes in §11,
offered, not placed), `readings_are_evidence_just_less_2026-09-29.md` (its
kink form is chosen here; O12 stays open), `null_series_tilt_and_property_N_2026-09-29.md`
§5 (the dichotomy, resolved by ruling 1), `ours28b/ours28b_no_borrowing_between_readings_2026-09-28.md`
(L19's scope, §2 here), `ours26_learned_omega_is_a_hint_2026-09-27.md`
(whose "the width problem lives in L10" this answers, §5), and the ours28
records in `../ours27/` (ours28's size law is dropped). The tester's ours25
(`calibrate: sandwich`) is kept there as it was. Binds nothing.

**Files.** `ours29_reference.py` (reference implementation: full $\Sigma_i$,
declared subspaces, readout), `ours29_reference_checks.py` (§7),
`ours29_reference_worlds.py` (§8), beside this record. In the tester
(`canvalidk/Mass-estimator-testor`, commit `645dfc7`): `misfit` and `calibrate: bartlett` on the
likelihood product, `tests/test_ours29.py`, `presets/ours29_noise_vs_signal.yaml`
(§9).

---

## 1. The rulings (30 Sept, 11:57)

Claude's three questions, his answers:

1. **L19's scope.** "Nothing about any reading's signal may be inferred from
   the others (unchanged, ironclad). The width of the mass law may be
   computed from each reading's own residual, the Birge move. Yes or no?"
   **"Yes good."**
2. **The move.** "Flat per-reading factors plus $\eta$, dropping ours28's
   per-reading arcsine as a declared population of setups. $g(Q)$ tempering
   stays either way." **"Good."**
3. **The standard.** "The tester counts a failure when the quoted width (the
   law widened by $\eta$) fails to cover, not when the raw law does. Or both
   get reported." **"Yes. but we can report both for now."**

## 2. The reason for L19, and its scope

His words (11:37): "L19 is no borrowing correct? That 1 is ironclad for a
reason. So if that reason is upheld we could wiggle it a bit. [...] We are
looking for the equation that can act as the definition of mass, so it can't
rely on things like a repeated experimental setup."

**The reason.** Readings can come from any experiment, in any order, with
nothing in common but the object. So the definition may not treat reading
$j$ in a way that presupposes reading $i$ resembles it. Anything that says
the signals in a series share a size, a direction or a spread is a claim that
the setups resemble each other.

**What the reason rules out.**

- ours28b's shared scale: every signal has the same typical size. Out, as
  ruled 28 Sept.
- **ours28's per-reading arcsine $G$.** Independent $\omega_i$ from one
  declared law says the setups are exchangeable draws from one population of
  experiments. With $\omega_i$ half-Cauchy (median 1), a long series is a
  priori sure to carry signal on about half its readings; a quiet series has
  prior weight of order $\varepsilon^N$ (scattered record §2). That is a
  claim about the ensemble of setups, weaker than a shared size and of the
  same kind. Dropped (ruling 2).
- The flat factor makes no cross-reading claim. Its defect is per reading
  (it overstates each reading's precision when the signal is weak) and
  setup-independent.

**What upholds it: $\eta$.** $\eta$ uses each reading's own score
$\ell_i'$ and curvature $\ell_i''$ at the series' best line, and adds them.
Nothing about reading $i$'s signal is inferred from reading $j$. The only
thing treated as shared is $m$, which L1 already shares. What $\eta$ reads is
whether the readings **support or contradict one another about $m$** by
more than their declared noise allows (his reading, 11:46). The scope,
ruled: **nothing about any reading's signal is inferred from the others; the
width of the mass law may be computed from each reading's own residual.**

**Precedent.** This is the practice the Particle Data Group uses when
combining measurements of one quantity from unrelated experiments: if they
disagree by more than their error bars allow, the combined error is inflated
by $S=\sqrt{\chi^2/(N-1)}$ (the Birge ratio), and only when $S>1$. $\eta$ is
the same move made on the law of mass rather than on Gaussian points, and
the cap at 1 is the "only when $S>1$". Like `==^`, it predates us.

**The caveat, stated before someone else does.** For the flat base, raising
the likelihood to $\eta$ gives exactly the law one would get by assuming every
signal shares one Gaussian size, $\omega^2=\eta/(1-\eta)$ (the tempering
identity). Same function, different premise: the sample mean is also the
Gaussian MLE, and taking a mean does not assume Gaussian errors. What
separates them is Z (§5): the shared-size model reads a zero reading as
evidence of small signals and moves (ours28b, 1.80 → 1.65 on one zero
reading); $\eta$ does not move.

## 3. The equation

A reading $i$ declares $\mathbf y_i=(\tilde{\mathbf F}_i;\tilde{\mathbf a}_i)$,
$\Sigma_i$ and $N_i$ ($d\times k_i$), as in the lemma record's type. The
starting law $\pi$ is declared (half-Cauchy at $c$). Nothing else is.

**Per reading**, for every candidate slope on the whole projective line,
$m=c\tan\phi$, $\phi\in(-\pi/2,\pi/2]$ (negative and infinite $m$ included):

$$
H_i(\phi)=\mathbf b^\top G^{-1}\mathbf b,\quad
C(\phi)=\begin{pmatrix}c\sin\phi\,N_i\\ \cos\phi\,N_i\end{pmatrix},\
G=C^\top\Sigma_i^{-1}C,\ \mathbf b=C^\top\Sigma_i^{-1}\mathbf y_i ,
$$

the explained energy (ours27's). $E_i$ is the whitened energy of the data
inside the span of all candidates, $\big(\begin{smallmatrix}N_i&0\\0&N_i\end{smallmatrix}\big)$,
and

$$
Q_i=E_i-\max_\phi H_i,\qquad \nu_i=k_i-1,\qquad
g_i=\min\!\Big(1,\frac{\nu_i}{Q_i}\Big),\qquad
\ell_i(\phi)=\tfrac12\,g_i\,H_i(\phi).
$$

**Series.** At $\hat\phi=\arg\max_\phi\sum_i\ell_i$ on the whole line,

$$
\mathcal H=-\sum_i\ell_i''(\hat\phi),\qquad J=\sum_i\ell_i'(\hat\phi)^2,\qquad
\eta=\min\!\Big(1,\frac{\mathcal H}{J}\Big)\quad(\eta=1\ \text{if}\ J=0).
$$

**Law.** $p(m)\propto\pi(m)\exp\!\big(\eta\sum_i\ell_i(m)\big)$ on $m>0$.

**Candidates and readout.** Each reading keeps its own flat cloud,
$\alpha_i\mid m\sim N(G^{-1}\mathbf b,G^{-1})$; the readout is T10 on the
tempered law, $\widehat m=\int mA\,p\,dm/\!\int A\,p\,dm$,
$A(m)=E[\lVert\mathbf a^\ast_{\rm stack}\rVert\mid m]$. Neither $g_i$ nor
$\eta$ touches $A$.

**Output** (ruling 3): the law's median and 95% interval (the quoted width),
the raw law's interval ($\eta=1$), the point, and $\eta$ itself. In the
tester, "report both" is two estimators in one preset, with and without
`calibrate: bartlett`, seeing the same readings.

**Declared:** $\pi$ (shape and centre, as before) and the middle of $g$'s
shape (O12). $G$ is gone.

## 4. The eigenvalue dictionary

His question (11:18): "Is it not the case that where there is energy, there
are eigenvalues?" It is. For noise that is the same in every spatial
direction, $\Sigma_i=S_i\otimes I_{k}$ (any correlation between the force and
acceleration channels), whiten each reading to $(\mathbf x_i,\mathbf y_i)$
and form its $2\times2$ Gram matrix

$$
M_i=\begin{pmatrix}\mathbf x_i\cdot\mathbf x_i&\mathbf x_i\cdot\mathbf y_i\\
\mathbf x_i\cdot\mathbf y_i&\mathbf y_i\cdot\mathbf y_i\end{pmatrix}.
$$

| the records' word | in eigen terms |
|---|---|
| a reading | $M_i$, plus $k_i$ |
| energy $T_i$ | $\operatorname{tr}M_i$ |
| alignment $\lvert C_i\rvert$ | eigenvalue gap $\lambda_+-\lambda_-$ |
| alignment fraction $a_i$ | eccentricity $(\lambda_+-\lambda_-)/(\lambda_++\lambda_-)$ |
| energy along the best line, $\max H_i$ | $\lambda_+$ |
| off-line misfit $Q_i$ (L21's input) | $\lambda_-$ |
| the reading's own best mass | the principal eigenvector's angle |
| explained energy $H_i(\theta)$ | the Rayleigh quotient $\mathbf u(\theta)^\top M_i\mathbf u(\theta)$, $\mathbf u=(\sin\theta,\cos\theta)$ |
| the series' state $(C,T)$ | $S=\sum_ig_iM_i$: its spectrum and frame |
| the law's tilt | $\tfrac12$ gap of $S$ |
| a null reading / null series | Wishart$_2(k,I)$ / Wishart$_2(kN,I)$ |
| "is there signal?" | sphericity of $S$ |
| $\mathcal H$ | gap of $S$ |
| $J$ | $\sum_i\big(g_i[M_i]_{+-}\big)^2$, the readings' off-diagonal entries in $S$'s frame |
| a zero reading | the zero matrix |
| $k$ | the Wishart degrees of freedom |

Checks: $\operatorname{tr}M_i=T_i$ and gap $=\lvert C_i\rvert$ to $10^{-12}$;
the median eccentricity of Wishart$_2(k,I)$ is 1.00, 0.87, 0.71, 0.38 for
$k=1,2,3,10$, reproducing the alignment record's table in closed-form terms.

**Correction to the conversation.** Claude said "$\lambda_-=2\,Q_{\min}$" at
11:18. It is $\lambda_-=Q_{\min}=\tfrac12(T_i-\lvert C_i\rvert)$, as the
evidence record defines $Q$. The code uses the correct one.

**Scope.** The $2\times2$ form needs $\Sigma_i=S_i\otimes I_k$, and the
series rows need one channel covariance for the whole series, so that every
reading's candidate angle is the same $\theta$ for the same $m$ (with several
instruments each reading has its own angle, and $S$ is not one matrix). For a general
$\Sigma_i$ the same objects exist (the energy in the span, the best line, the
off-line residual, the series' best line and its curvature) but they are not
the spectrum of one $2\times2$ matrix. The reference implementation works
with the general form.

**What becomes automatic.** Rotation invariance (Gram matrices do not see the
frame), Z for the law (the zero matrix adds nothing), and $k$'s role (degrees
of freedom). **What the rename must not do.** The principal axis of $S$ is the
state, not the readout: it is the Deming root the spec demoted, and the
readout stays the ratio of means forced by composition. And on $m>0$ the
principal axis falls outside the allowed range on about half of null series,
which is the fibre structure the chain audit found, seen again.

## 5. What is forced and what is chosen

| piece | status | why |
|---|---|---|
| flat per-reading factor | follows from ruling 2 | L19's reason excludes a declared population of setups; with no size law the factor is the flat one, read as a working likelihood |
| $g(0)=1$, the $\kappa/Q$ tail | forced (T30) | Z; bounded evidence that is "less, not none" |
| $\kappa=\nu=k_i-1$ | forced by meaning (T30) | the reading's codirectionality conditions |
| the kink middle, $\min(1,\nu/Q)$ | **chosen** | flat inside the noise band, which the own-world test preferred; O12 stays open |
| $\eta$'s form $\mathcal H/J$ | forced in the long-series limit | property N: the width must track how far the point moves (Bartlett identity) |
| one $\eta$ per series | forced as a count | $m$ is one number, so one calibration equation |
| **raw $J$** (no $N/(N-1)$) | **forced by Z** | an exact zero reading has $\ell_i'=\ell_i''=0$; a Bessel factor counts it and moves. ours25's sandwich does |
| **the whole projective line** | **forced** | on $m>0$ the best line can sit at a boundary where the curvature measures nothing; ours25 then falls back to the prior (3 of 40 pure-noise 3D series of 12 readings in one check) |
| evaluation at the best line | forced | the ratio is parametrisation-free only where the score vanishes; gives the swap for free |
| cap at 1 | natural | never sharpen from the series; one reading gives its own law |
| candidates unchanged | forced by L19 | shrinking each cloud by $\eta$ would infer signal size from the series |
| the power form at finite $N$ | **chosen** | other finite-$N$ forms (a Gaussian with the sandwich variance, a bootstrap) agree only asymptotically |

**The hint answered.** `ours26_learned_omega_is_a_hint_2026-09-27.md` ruled
learned $\omega$ a hint and located the width problem in L10: "a solution
should come from what L10 should say, not from fitting $\omega$." ours29's
answer: L10 says no size law is declared; the flat factor is a working
likelihood, not the claim that every signal is infinitely large; and how far
that working likelihood overstates the precision is measured on the series by
$\eta$, from the readings' disagreement about $m$, never from their sizes.
That is why learned $\omega$ was the right shape (§2's caveat) and the wrong
object (it reads sizes, and fails Z).

## 6. The afternoon's other results

**6.1 The combination case** (11:12, before the rulings). A null series at
the noise scale $s=1$ kg, combined with an independent density result $B$.
The exact product rule (contributions multiplied, $\pi$ once) is the reference.

- $B=1000$ kg: every rule gives about 1000 kg. The pull is $-0.010$ in
  $\ln m$ at every $N$, the residue of $\pi$'s declared centre, and the null
  series adds nothing. The case does not discriminate.
- $B=1$ kg $\pm20\%$, at the noise scale: rerun on ours29's base (1D, $g=1$),
  pull of the combined median in $\ln m$:

  | $N$ | $\eta$ | raw law | tempered by $\eta$ |
  |---|---|---|---|
  | 10 | 1.00 | +0.036 | +0.036 |
  | 100 | 0.27 | +0.196 | +0.053 |
  | 1000 | 0.06 | **+0.440** | +0.031 |

  Pure noise moves an independent 1 kg measurement by 55% at $N=1000$ under
  the raw law; tempered, by 3%.
- The data-only resolution (property N) gives zero pull but throws away
  precision on good series (computed at 11:12 on ours28's base, 1D: true 2
  kg, $N=20$, signal 3: law width 1.00 in
  $\ln m$, resolution width 10.2; $N=200$: 0.27 against 2.45). So the
  resolution is the `==^` check, not the weight: in a contradiction (a 2 kg
  series against 1000 kg) the product returns 990 kg with a narrow width and
  hides it; the resolution and the (point, width) summary both show 25σ.
- For a near-Gaussian contribution, the power $\eta$ equals quoting
  $\sigma/\sqrt\eta$, so "in the law" and "beside the law" give the same
  combination. The ruling that matters is which combination rule the
  practice uses across results; ruling 3 moves the tester's standard to the
  quoted width.

**6.2 The $\lambda$ move on the size law** (11:18, his question: "can we do
something like how we ruled out lambda?"). With $\lambda$, the exact
combination showed a fitted prior factor was carried state. Same move here:
for the flat base, the power $\eta$ is exactly a candidate measure, one
shared Gaussian size with $\omega^2=\eta/(1-\eta)$. So the combination factor,
read back as a shape, is a shared scale set by the calibration requirement,
not by a declared $G$: $G$ disappears as $\lambda$ did, and $\eta$ is state
(it starts at 1 and every reading updates it; for ours27 in 1D it is carried
exactly by $(C,T,\sum|C_i|^2,\sum C_i^2)$). Inside the law it would collide
with L19 as then worded; ruling 1 is what lets it stand, because it is
computed from residuals, not read as a size.

**6.3 A Z-safe noise statistic.** $\text{gap}(S)/(2\sqrt J)$ is
Rayleigh(1) under pure noise at every $N$ (medians 1.49, 1.20, 1.18 at
$N=10$, 100, 1000; Rayleigh 1.18), and a zero reading leaves it alone. It is
the null-series record's O-c statistic with $\sqrt N$ replaced by $\sqrt J$.
Not used by ours29; recorded as the natural "is there signal" test in these
terms.

## 7. Exact checks (reference implementation)

`python ours29_reference_checks.py`, full covariance with cross-channel
correlation where stated:

- **Both steps off is ours27:** readout and median equal the ours28
  reference's flat branch to $10^{-10}$ on its two test series.
- **1D $\eta$** against the carried-state closed form: $7\times10^{-9}$.
- **A weak 3D series** (40 readings of pure noise, full covariance):
  $\eta=0.403$; 4 of 40 readings weighted below 1 (min 0.62). Quoted 95%
  interval $[1.30,57.8]$ against the raw law's $[1.87,34.0]$ and ours27's
  $[1.63,15.0]$ ($c=2$).
  - each reading rotated with its covariance: $\eta$, readout, median
    unchanged to $3\times10^{-8}$;
  - swap: readout and median products 1.000000, $\eta$ unchanged;
  - units (force $\times3$): ratios 3.000000;
  - ten exact zero readings appended: $\eta$ and the law unchanged exactly;
    **the readout moves, 3.5983 → 3.6017** (T29, below).
- **Declared subspace:** a 3D reading declared $k=1$ gives the 1D reading's
  law and readout exactly; a misaligned 3D reading ($Q=2.34$, $g=0.853$)
  declared $k=1$ has $Q=0$, $g=1$.
- **One reading:** $\eta=1$ for aligned, weak and zero readings.
- **The impossible pair** ($\tilde{\mathbf F}=(2,0)$, $\tilde{\mathbf a}=(0,1)$,
  2D): $Q=1/\sigma^2$, and the tilt is exactly **1.500 nats at every
  $\sigma\le1$**, the limit $\nu a/(1-a)$. On the flat base the limit is
  reached as soon as $Q\ge\nu$ (on ours28's base it was approached, 1.48 at
  $\sigma=0.1$).
- **Totality:** perpendicular, anti-parallel, force-only and
  acceleration-only readings all finite; the all-zero experiment gives
  median $c$ and $\eta=1$.
- **Axis splitting:** 30 3D readings against their 90 axes, $\eta$ 0.458
  against 0.445, intervals slightly different, as the rewritten L8 intends.

In the tester, `tests/test_ours29.py` (13 tests, all pass with the existing
861): the closed form against the numerical path and a 200 001-point
brute-force grid in 1D, 2D and 3D; $g$ against the Gram-eigenvalue rule; Z;
the sandwich's Bessel factor moving by exactly $(12/11)/(24/23)$ when 12 zero
readings are appended while $\eta$ stays; rotation, swap, units; the
impossible pair ($1.5$ against $1.5/\sigma^2$ untempered); $\eta$ against the
sandwich up to $N/(N-1)$ where the peak is inside $m>0$; the readout weight
unchanged; the refusals. The tester and the reference agree to $10^{-8}$ on
null series in 1D, 2D and 3D.

## 8. Worlds the tester cannot generate (reference, law only)

`python ours29_reference_worlds.py`, seed 30. Coverage of the 95% interval:

| world | ours27 | ours29 raw law | ours29 quoted | median $\eta$ | series |
|---|---|---|---|---|---|
| (a) rescue: 20 readings of signal 3 plus one impossible reading ($\sigma=0.1$, $Q=100$), 2D | 3.3 | 95.3 | 95.3 | 1.00 | 150 |
| (b) one reading of signal 6 among 100 pure-noise readings, 3D | 65.0 | 62.0 | 100.0 | 0.19 | 100 |
| (c) every reading its own instrument ($\sigma_F,\sigma_a$ each $10^{U(-2,2)}$) and its own signal (0, 0, 0.5, 1, 3 or 10), 3D, $N=40$ | 92.0 | 93.0 | 94.0 | 1.00 | 100 |
| (d) pure noise, planted mass drawn from $\pi$ itself, 2D, $N=100$ | 53.0 | 50.0 | 83.0 | 0.18 | 100 |
| (d) the same, 3D | 42.0 | 40.0 | 75.0 | 0.17 | 100 |

- **(a)** is the misfit weight's work: $\eta$ stays at 1 and the raw law
  already covers. One impossible reading otherwise takes ours27 to 3%.
- **(b)** is $\eta$'s: one informative reading drowned in pure noise. The
  raw law is too narrow because the null readings add score noise and no
  curvature (T26); $\eta$ widens it, here to over-coverage.
- **(c)** is the world L19's reason is about. Nothing is shared, so there
  is nothing to correct and $\eta$ stays at 1: the wiggle costs nothing
  where the readings do not resemble each other.
- **(d)** is the null floor (§13, item 4). Against a mass drawn from $\pi$
  the quoted law covers 75–83%, not 95%: the tilt is bounded, not zero. It
  falls with dimension.

## 9. The tester study

`presets/ours29_noise_vs_signal.yaml`: 1D, 2D, 3D; planted mass 0.25, 1, 4
($s=1$); signal 0, 0.5, 1, 3 noise units per reading, random directions;
$N=20$, 200; 200 series per cell (MCSE about 1.5 points at 95%). Four
estimators on the same readings: ours27, ours29's raw law (misfit weight
only, $\eta=1$), ours29 (quoted), and ours25 (the sandwich).

**Coverage of the 95% interval** (percent; mean over the three planted
masses, worst mass in brackets). "Raw law" is ours29 with $\eta=1$ (misfit
weight only); "quoted" is ours29. Last column: median log width, raw law →
quoted.

| $d$ | signal | $N$ | ours27 | ours29 raw law | **ours29 quoted** | ours25 | width |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 20 | 73.5 [64] | 73.5 [64] | **81.0** [76] | 82.5 [75] | 4.98 → 5.70 |
| 1 | 0 | 200 | 39.3 [37] | 39.3 [37] | **84.5** [82] | 72.2 [71] | 4.51 → 5.83 |
| 1 | 0.5 | 20 | 82.8 [79] | 82.8 [79] | **87.7** [86] | 86.8 [84] | 3.73 → 4.31 |
| 1 | 0.5 | 200 | 77.0 [66] | 77.0 [66] | **94.0** [94] | 91.8 [87] | 0.85 → 1.58 |
| 1 | 1 | 20 | 89.0 [87] | 89.0 [87] | **91.3** [91] | 91.8 [92] | 1.46 → 1.76 |
| 1 | 1 | 200 | 91.3 [88] | 91.3 [88] | **97.0** [96] | 97.0 [96] | 0.39 → 0.48 |
| 1 | 3 | 20 | 95.7 [94] | 95.7 [94] | **97.0** [96] | 97.7 [96] | 0.42 → 0.44 |
| 1 | 3 | 200 | 94.0 [93] | 94.0 [93] | **95.3** [95] | 95.5 [96] | 0.13 → 0.13 |
| 2 | 0 | 20 | 64.7 [62] | 65.0 [62] | **83.3** [82] | 80.2 [78] | 4.81 → 5.91 |
| 2 | 0 | 200 | 32.7 [30] | 32.8 [28] | **83.5** [79] | 69.8 [67] | 4.42 → 5.83 |
| 2 | 0.5 | 20 | 79.2 [72] | 80.8 [73] | **88.8** [84] | 87.5 [84] | 3.72 → 4.81 |
| 2 | 0.5 | 200 | 70.8 [58] | 71.5 [57] | **94.7** [93] | 91.5 [84] | 0.94 → 2.46 |
| 2 | 1 | 20 | 89.2 [85] | 91.0 [88] | **95.3** [94] | 94.3 [92] | 1.52 → 2.06 |
| 2 | 1 | 200 | 83.5 [74] | 84.7 [74] | **93.7** [94] | 93.2 [92] | 0.42 → 0.56 |
| 2 | 3 | 20 | 94.2 [94] | 95.5 [94] | **96.2** [95] | 95.8 [96] | 0.45 → 0.46 |
| 2 | 3 | 200 | 93.0 [91] | 93.8 [92] | **94.7** [93] | 94.0 [92] | 0.14 → 0.14 |
| 3 | 0 | 20 | 59.0 [58] | 60.5 [60] | **81.8** [78] | 76.2 [74] | 4.69 → 6.04 |
| 3 | 0 | 200 | 25.5 [24] | 28.5 [26] | **84.3** [82] | 70.5 [70] | 4.47 → 5.93 |
| 3 | 0.5 | 20 | 75.2 [68] | 75.5 [68] | **88.0** [84] | 89.0 [86] | 3.26 → 4.83 |
| 3 | 0.5 | 200 | 60.8 [48] | 62.2 [48] | **94.0** [92] | 89.8 [80] | 0.93 → 3.05 |
| 3 | 1 | 20 | 83.8 [76] | 84.3 [76] | **93.3** [92] | 93.2 [91] | 1.67 → 2.62 |
| 3 | 1 | 200 | 78.3 [65] | 80.2 [68] | **94.2** [94] | 94.2 [92] | 0.41 → 0.63 |
| 3 | 3 | 20 | 94.5 [90] | 95.2 [92] | **96.0** [94] | 97.3 [96] | 0.45 → 0.47 |
| 3 | 3 | 200 | 92.0 [88] | 93.0 [90] | **94.0** [93] | 95.2 [94] | 0.14 → 0.14 |

**Reading.**

- **Wherever there is signal and $N=200$, the quoted width covers 93.7–97%**
  in every dimension, where the raw law covers 62–94% (62–91% below signal
  3). The worst cell before,
  3D at signal 0.5, goes from 60.8% (worst mass 48%) to 94.0%.
- **At strong signal $\eta\approx1$ and nothing changes** (widths equal to two
  decimals).
- **Pure noise:** the raw law's coverage collapses with $N$ (28–39% at
  $N=200$; ours27 26–39%); the quoted width holds at 81–85% at both $N$. This is the null
  floor: bounded, not removed (§13, item 4).
- **Weak signal with few readings** (0.5, $N=20$): 88–89%. With 20 readings
  $J$ is a noisy estimate, biased low, and the cap cannot help; this is the
  finite-$N$ item (§13, item 2).
- **ours25 is close where the series has signal** and worse on pure noise
  (70–72% at $N=200$): on $m>0$ its search finds no interior peak on a
  share of null series and it falls back to the prior's width about a
  misplaced point. The whole-line search is what separates them there.
- **The misfit weight alone** changes little in these worlds (raw law against
  ours27: at most +3 points). Nothing here is misaligned beyond its noise;
  its work is §8 (a).

## 10. The readout

The readout is T10 on the tempered law, with each reading's candidate cloud
unchanged. $\eta$ never moves the law's peak, but it moves the readout,
because the readout averages over the law.

**On pure noise $\eta$ stops the point wandering.** Planted mass 1 (the
prior's centre), quartiles of the readout:

| $d$ | $N$ | ours27 | ours29 |
|---|---|---|---|
| 1 | 20 | 0.40 / 0.98 / 2.30 | 0.49 / 0.99 / 1.83 |
| 1 | 200 | 0.19 / 0.92 / 3.98 | 0.62 / 0.98 / 1.48 |
| 2 | 20 | 0.28 / 0.81 / 2.15 | 0.47 / 0.96 / 1.41 |
| 2 | 200 | 0.14 / 0.83 / 3.64 | 0.55 / 0.98 / 1.53 |
| 3 | 20 | 0.30 / 1.02 / 3.18 | 0.63 / 1.02 / 1.69 |
| 3 | 200 | 0.16 / 0.95 / 6.58 | 0.61 / 0.99 / 1.50 |

Under ours27 the interquartile range widens with $N$ (in 3D from a factor 10
to a factor 40); under ours29 it stays between a factor 2.4 and 3.7 and does
not grow with $N$.
This is the requirement of `c_is_the_noise_measure_at_zero_signal_2026-09-29.md`
§1, met: as the evidence of signal goes to zero, the point goes to the
zero-signal value continuously. It was not designed in; it follows from the
readout being an average over a law that $\eta$ has flattened.

**Accuracy** (within a factor 2 of the planted mass, ratio of means): on
pure noise at $N=200$, 16–21% → 31–33%; at signal 0.5, $N=200$, up by 1–6
points; at signal 3, within 0.3 points.

**At strong signal in 2D and 3D the readout moves a little.** Series by
series against ours27 at signal 3: median $|\log|$ difference 1.5–1.7% at
$N=20$, 0.6% at $N=200$ (1D: zero). This is the misfit weight, not $\eta$
($\eta=1$ there): with the kink at $Q=\nu$, the ordinary residual of a
normal reading exceeds $\nu$ often ($P(\chi^2_1>1)=32\%$ in 2D,
$P(\chi^2_2>2)=37\%$ in 3D), so about a third of ordinary readings are
weighted slightly below 1. It costs nothing measurable here, but it is
evidence on O12: a kink placed at the mean of $\chi^2_\nu$ is not "flat
inside the noise band", only flat inside half of it.

**Z (O1) is untouched.** Ten exact zero readings appended to the weak 3D
series of §7 leave $\eta$ and the law unchanged and move the readout
3.5983 → 3.6017: the stacked cloud counts the zero readings (T29). ours29
neither causes nor cures it.

## 11. What this changes in the lemma record (offered, not placed)

- **T8 (the law):** ours29's law (§3) replaces ours28's.
- **L10′ (per-reading scale under $G$):** struck by ruling 2. **L10**
  restated: no size law is declared; the flat factor is a working likelihood
  (§5, the hint answered).
- **L4′ (the size law's symmetry)** and the $G$ row of §4 (declared shapes):
  moot with $G$ gone. The "one rule, two uses" remark keeps only its $\pi$
  half.
- **L19:** scope added, as ruled (§2).
- **New T31 (ours29's $\eta$):** raw $J$ forced by Z; the whole projective
  line; evaluation at the best line; cap at 1; candidates unchanged by L19.
- **T18** (the flat law's width error): measured per series by $\eta$.
  **T26** (limits of per-reading rules) stands, and is why $\eta$ is
  series-level. **T27** (own-world calibration) has no declared world left to
  apply to; the standard is ruling 3's.
- **O3** becomes more load-bearing, not less: ours29 uses the flat weighting
  (T6's endpoint) and the stacked length (T9), whose support the rewritten L8
  withdrew. **O6** shrinks to $\pi$'s shape. **O12** stays (the kink is
  chosen). **O1** is untouched (§10).

## 12. His words, 30 Sept (verbatim)

> So can we do something like how we ruled out lambda? Can we use the
> combination to then infer which shape inside the estimator would amount to
> be the combination factor? And think about this only if it matters: Is it
> not the case that where there is energy, there are eigenvalues?

> So what now? Does this imply a reframing pass where we rename the objects
> and reinterpret the problem we are facing?

> Ok looks like we are getting one connection after another. But L19 is no
> borrowing correct? That 1 is ironclad for a reason. So if that reason is
> upheld we could wiggle it a bit. It's important that you understand what
> that reason is because its integral to what we are going for here. We are
> looking for the equation that can act as the definition of mass, so it
> can't rely on things like a repeated experimental setup.

> Okay I don't understand this well. What I think you might be saying is, if
> we do have multiple readings and we use the weightings appropriately, the
> way we combine the results might very much resembles assuming some of their
> set up is shared? But that comes out of the idea that regardless of how we
> combine the information, 2 measurements will either support or contradict
> eachother?

> Alright that makes sense. Then what do we do first? Make the reframing
> record? Start with the next equation? what makes sense?

> 1. Yes good. 2. Good. 3. Yes. but we can report both for now. Go ahead.

## 13. Open

1. **O12**, the middle of $g$: the kink is chosen, not forced, and §10
   shows its threshold at $\nu$ touches about a third of ordinary readings in
   2D and 3D. A threshold at a high quantile of $\chi^2_\nu$, or a smooth
   form with the same two ends, would leave ordinary readings alone.
2. **The finite-$N$ form of $\eta$.** The power is asymptotically forced; at
   small $N$ the raw $J$ is biased low (by about $(N-1)/N$ when the model is
   right), which the cap absorbs. A Z-safe small-sample correction would
   count effective readings, $(\sum_i c_i)^2/\sum_i c_i^2$ with $c_i$ the
   curvature contributions (zero readings contribute nothing), not $N$. Not
   built.
3. **Several instruments.** The ours28b calibration theorem makes sharing
   exact only within classes of equal sensitivity (per noise pair). One
   $\eta_g$ per noise pair may be what property N asks for. Not checked.
4. **The null floor.** $\eta$ stops the law's tilt growing but leaves it
   bounded (about 2 nats; coverage 75–83% against a mass drawn from $\pi$ in
   2D and 3D, §8 (d), and near 80% in the morning's 1D check). Under Z and L19, null evidence does not
   accumulate (scattered record §5), so this may be the most the law can do;
   "we know nothing" then lives in the resolution. Not proven.
5. **O1**, the readout under Z (§10).
6. **O3**, now load-bearing (§11).

## 14. Checks and reproduction

- Reference: `ours29_reference.py` (numpy only), checks and worlds scripts
  beside it; the checks need `../ours27/ours28_reference.py` for check 1.
  Projective line on 2001 (reading) and 4001 (series) points, refined by
  golden section; curvature and score by central differences,
  $h=10^{-4}$ in $\phi$; law on 721 points over $\ln c\pm9$.
- Tester (commit `645dfc7`): `python -m met run presets/ours29_noise_vs_signal.yaml`;
  `python -m pytest -q tests/test_ours29.py`. The study in §9 was run in
  Claude's sandbox in two parts split by dimension (cells do not depend on
  where they sit in the grid); per-cell data are identical to a single run.
- §6.1: sandbox scripts (numpy, scipy), 1D, seed 11; not kept.
- §6.3 and §4's Wishart medians: sandbox, 2000 and 20 000 draws; not kept.
