# Axis splitting: what the 24 September change is, and what forces it

**2026-09-25.** Record, filed at the user's request ("write up the record"),
written by Claude from the conversation of the night of 24–25 September. It
follows `excitation_weighting_premise_2026-09-24.md` in the same discussion
and concerns the 21 September specification in `mass_estimator_equation/`.
It edits neither and binds nothing.

**In one paragraph.** The 24 September equation is the 21 September one
with every candidate true pair $(f,\alpha)$ reweighted by $r^{d-2}$, where
$r$ is the pair's length in noise units. That one weight is the whole
change. It makes the law of mass the same in every dimension: the 21
September law's 2D form. On accuracy the two equations tie, because where
they differ the answer is dominated by the pull toward $s=\sigma_F/\sigma_a$.
What the new one gains is consistency. The same measurement gives the same
mass law however its axes are grouped, and the old equation fails this.
That requirement, **axis splitting**, is elementary: Newton II is three
component equations sharing one mass. Taken as a premise, it forces the
flat-in-the-vector weighting (without the translation step the previous
record flagged as weak) and it forces the mass factor that the previous
record left open. With many readings both equations converge to the true
mass, and the new one's too-narrow intervals now have an exact formula.

---

## 1. The premise

**Axis splitting (A).** A $d$-dimensional reading, and its $d$ axes treated
as $d$ one-dimensional readings of the same mass, must give the same law
of mass.

- **Why it is elementary.** $\mathbf F^\ast=m\,\mathbf a^\ast$ is
  $F^\ast_k=m\,a^\ast_k$ for each axis $k$: three equations sharing one
  mass. Whether the channels are recorded as one vector or as three
  numbers is bookkeeping. The mass should not depend on it.
- **The user's judgement** (2026-09-25): "axis splitting seems too
  elementary to be wrong."
- **Where it came from.** The previous record already stated it as a
  *consequence* of the flat-vector weighting ("a $d$-dimensional reading is
  exactly $d$ one-dimensional readings"). This record reads it the other
  way, as a requirement, and asks what it forces.
- **The same kind as (L) and (R).** The ratio-of-means readout was forced
  by two composition laws of Newtonian mass: lumping at common acceleration,
  and common force (`what_forces_the_ratio_of_means_readout.md`). Axis
  splitting is a third composition law: the estimator must split as
  Newton II splits.

"Combined" means the symmetric rule: the readings' likelihoods multiplied,
each with its latent pair integrated under the radius split, and the mass
factor counted once (contribution 12; the tester's `likelihood_product`,
nuisance `radius`). §6 shows why this is the combination rule to use: it is
consistent.

## 2. The whole change is one weight

Write the candidate pair in noise units, $f=\sigma_F\,r\sin\theta$,
$\alpha=\sigma_a\,r\cos\theta$, so $m=s\tan\theta$. The 24 September measure
is $d^d\mathbf v\,d\theta$ with $\mathbf v=r\mathbf u$, and

$$d^d\mathbf v\,d\theta \;=\; r^{d-1}\,dr\,d\Omega\,d\theta \;=\; \frac{r^{\,d-2}}{\sigma_F\sigma_a}\;df\,d\alpha\,d\Omega .$$

So the new readout is the old one with each candidate pair reweighted:

$$\widehat m_{24}\;=\;\frac{E_{21}\!\left[f\,r^{d-2}\right]}{E_{21}\!\left[\alpha\,r^{d-2}\right]},\qquad r^2=\Big(\frac{f}{\sigma_F}\Big)^2+\Big(\frac{\alpha}{\sigma_a}\Big)^2 .$$

| dimension | weight on each candidate pair |
|---|---|
| 1D | $1/r$ |
| 2D | $1$ (no change) |
| 3D | $r$ |

Checked by brute-force quadrature in the 21 September coordinates to
$10^{-9}$ (tester `tests/test_cartesian.py`). In 2D, all 4000 series in
every cell of the head-to-head run agree exactly.

## 3. What the weight does

**The law of mass becomes dimension-free.** Integrating the true pair out
leaves a kernel in $h=\lvert\mathbf x\sin\theta+\mathbf y\cos\theta\rvert$
($\mathbf x=\mathbf F/\sigma_F$, $\mathbf y=\mathbf a/\sigma_a$):

| | 1D | 2D | 3D |
|---|---|---|---|
| 21 Sept | $1+\sqrt{\pi/2}\,h\,e^{h^2/2}\operatorname{erf}(h/\sqrt2)$ | $e^{h^2/2}$ | $\sqrt{\pi/2}\,e^{h^2/2}\operatorname{erf}(h/\sqrt2)/h$ |
| 24 Sept | $e^{h^2/2}$ | $e^{h^2/2}$ | $e^{h^2/2}$ |

So a reading's 24 September law of mass *is* the 21 September law of the 2D
reading with the same $\lvert\mathbf x\rvert$, $\lvert\mathbf y\rvert$ and
$\mathbf x\cdot\mathbf y$. Relative to that 2D form, the old 1D law was
sharper and the old 3D law flatter. The readout keeps one dimension
dependence: $E[\alpha\mid m]$ is the mean length of a $d$-dimensional
Gaussian vector (the noncentral chi mean).

**The only behavioural change is the pull toward $s$.** At weak signal
every estimate is pulled toward $s=\sigma_F/\sigma_a$, the median of the
zero-reading law. The new law pulls **less in 3D and more in 1D**. It
changes nothing at zero signal, at $m=s$ on a noise-free reading (the
self-reciprocal point), or from acceleration SNR about 4 upward. That one
mechanism explains every row of the previous record's single-reading
table.

## 4. What axis splitting forces

**(a) The 21 September equation fails it.** One reading, same numbers,
median estimate ($\sigma_F=0.5$, $\sigma_a=0.1$, so $s=5$):

| reading | 21 Sept | 24 Sept |
|---|---|---|
| $F=0.3$, $a=0.15$ as 1D | 3.12 | 3.46 |
| the same, written as a 3D vector $(x,0,0)$ | 3.73 | 3.46 |
| $\mathbf F=(0.3,0.2,-0.4)$, $\mathbf a=(0.15,-0.05,0.02)$ as one 3D vector | 3.97 | 3.72 |
| the same three axes as three 1D readings | 3.46 | 3.72 |

The x-axis alone gives 3.12 (21 Sept) and 3.46 (24 Sept). That differs from
the last two rows for both equations, as it should: dropping the y and z
readings is less data, not another grouping of the same data.

**(b) The excitation weighting must be flat in the vector.** Splitting
gives each axis its own latent scalar $a^\ast_k$. For the vector reading to
agree, the weighting on $\mathbf a^\ast$ must be a product over axes. The
21 September weighting $\lvert\mathbf a^\ast\rvert^{2-d}$ is not.

- Product over axes plus isotropy (already in the 21 September equation):
  by Herschel–Maxwell, a centred Gaussian $N(0,\omega^2 I)$, or among
  improper weightings $e^{c\lvert\mathbf a\rvert^2}$.
- No preferred excitation scale then leaves the flat one, $\omega\to\infty$.

This is the previous record's second route with its criticised step
(independent components) replaced. It is no longer a modelling claim about
the pushes. It is a consistency demand on the inference: the vector must
agree with its components. **It does not use the translation step**, so the
accelerating-frame objection of the previous record's §5 does not reach it.

**(c) The mass factor is forced too.** The premise fixed the weighting only
up to a factor in $m$ alone. Under the symmetric rule each reading's
likelihood is $e^{h_k^2/2}$ and the mass factor is counted once. So a
$d$-dimensional reading and its $d$ axes agree only if the single-reading
mass factor is the same in every dimension.

- **The 2D anchor.** In 2D the 21 September equation already agrees with the
  premise, and there "flat in $\mathbf a^\ast$ with flat in $m$" is the same
  measure as "flat in $\mathbf F^\ast$ with flat in $1/m$" and as "flat in
  $\mathbf v$ with flat in $\theta$". Keeping it there fixes the factor:
  uniform in $\theta$, the record's completion.
- **Consequence.** Given axis splitting, the 24 September equation adds **no
  freedom** the 21 September one did not already have. The one remaining
  free function is the 2D mass factor, which is investigation question 2 as
  before.

**Why 2D never had this problem.**

$$d^d\mathbf F^\ast\,d(1/m)\;=\;m^{\,d-2}\;d^d\mathbf a^\ast\,dm .$$

The premise read with force and read with acceleration give the same
measure only when $d=2$.

**(d) The rival completion fails.** The 23 September direction asked that
the default weighting reference no instrument. With the premise, the only
completion satisfying that is

$$(f\alpha)^{(d-2)/2}\,df\,d\alpha\,d\Omega=\sqrt{d^d\mathbf a^\ast dm\cdot d^d\mathbf F^\ast d(1/m)} ,$$

the geometric mean of the two readings of the premise. It has no $\sigma$
in it. It fails axis splitting (its mass factor changes with $d$) and it
behaves badly (§5). **In 1D and 3D, axis splitting and the no-instrument
requirement cannot both hold.** This record takes axis splitting.

**Not split-invariant: the readout.** The ratio of means uses vector
lengths, and the sum of axis lengths is not the vector length. The laws
agree to $10^{-15}$, but the readouts of one 3D reading and its three axes
differ by 0.3–0.4% (median) and 1.1–1.4% (90th percentile) over 200 random
readings. If axis splitting is a law, the readout owes an account of this.

## 5. Forced and chosen, measured

Taken apart, the forced factor $\alpha^{d-2}$ and the completion are each
large, and they cancel. In 3D the forced factor alone with a flat-in-$m$
completion (the previous record's §3 comparison law, "plain") drags
readings at $m=s$ to a median of 0.65–0.79 of the truth. The completion
brings them back to 1.00. So the new equation stays close to the old one
*because of* the completion. That is why it matters that axis splitting
forces the completion rather than leaving it to choice.

3D, one reading, 4000 series per cell, ratio of means
(`presets/sept24_completions_3d.yaml`):

| cell | within ×2 (%): 21 Sept / plain / 24 Sept / σ-free | median $\widehat m/m$: 21 Sept / plain / 24 Sept / σ-free |
|---|---|---|
| $m=s/4$, SNR 2 | 34.3 / 53.2 / 37.9 / 28.2 | 2.56 / 1.89 / 2.44 / 2.71 |
| $m=s$, SNR 1 | 79.5 / 72.3 / 74.7 / 85.0 | 1.00 / 0.68 / 1.00 / 1.00 |
| $m=4s$, SNR 1 | 70.7 / 50.7 / 71.4 / 65.8 | 0.70 / 0.51 / 0.72 / 0.61 |

In 1D (`presets/sept24_completions_1d.yaml`; the plain law is improper
there), the σ-free completion overshoots: median $\widehat m/m$ is 1.18 at
$m=4s$, SNR 2. At SNR 8 it still trails, with 76.2% of series within ×2
against 89.5% (21 Sept) at $m=s/4$. The 21 September equation sits between
the record's completion and the σ-free one in both dimensions.

## 6. Many readings: consistency and the width formula

**Both equations are consistent** (true noise supplied, independent
isotropic channels). For a candidate angle $\theta$, each reading gives
$\mathbf w=\mathbf x\sin\theta+\mathbf y\cos\theta\sim N(\lambda(\theta)\,\hat{\mathbf e},\,I_d)$
with

$$\lambda(\theta)=r^\ast\cos(\theta-\theta_0),\qquad r^{\ast2}=\frac{\lvert\mathbf F^\ast\rvert^2}{\sigma_F^2}+\frac{\lvert\mathbf a^\ast\rvert^2}{\sigma_a^2},\qquad \tan\theta_0=\frac{m_0}{s},$$

where $r^\ast$ is the true pair's noise-unit length and $\theta_0$ the true
angle. The expected log-likelihood is $E[\log K(\lvert\mathbf w\rvert)]$,
increasing in $\lambda$ for any kernel $K$ increasing in $h$. So it peaks at
$\theta_0$ for every reading, whatever its excitation. This answers the 23
September direction's "prove or refute convergence" for the symmetric rule.
It also explains why rules that carry a mass factor per reading (eq. 18,
eq. 20, the $(f\alpha)^k$ family) settle away from the truth: each reading
adds the same tilt.

**The spread and the interval width.** At $\theta_0$,
$\mathbf w_\perp=\mathbf x\cos\theta-\mathbf y\sin\theta\sim N(0,I_d)$ is
independent of $\mathbf w$. With the per-reading log-likelihood $G(h)$:

$$J=E\big[G'(h)^2\big],\qquad H=E\Big[-G''(h)-G'(h)\,\frac{d-1-h^2}{h}\Big],\qquad h\sim\chi_d(r^\ast).$$

The estimate's variance is $J/(NH^2)$ and the posterior's is $1/(NH)$. For
the 24 September law $G=h^2/2$, so $J=d+r^{\ast2}$ and $H=r^{\ast2}$ exactly:

$$\frac{\text{true variance}}{\text{posterior variance}}\;=\;1+\frac{d}{r^{\ast2}} .$$

**This is the undercoverage** the previous record found and the flat
(infinite-temperature) limit was blamed for, now as a formula. It is the
functional errors-in-variables sandwich in this problem's variables (cf.
Fuller 1987, cited from general knowledge, unverified). The 21 September
law has no closed form. Its factor is smaller in 3D and larger in 1D, and
the two spreads differ by at most 4%.

**Tested** (`presets/sept21_vs_sept24_repeated.yaml`, new excitation per
reading, random direction, 500 series per cell). At $N=200$, 80% interval
coverage, observed (predicted):

| cell | 21 Sept | 24 Sept |
|---|---|---|
| 1D, $m=s$, SNR 1 | 64.4 (66.3) | 68.4 (70.5) |
| 3D, $m=s$, SNR 1 | 60.4 (61.4) | 57.6 (58.2) |
| 3D, $m=s/4$, SNR 1 | 53.2 (52.2) | 51.2 (48.8) |
| 3D, $m=4s$, SNR 3 | 78.4 (79.7) | 78.0 (79.6) |

Medians are within about 3% of the truth in every cell at $N=200$. The
weakest 1D cell ($m=s/4$, SNR 1) is still approaching its prediction from
above (73.0 → 66.6 → 63.0 against 59.1 at $N$ = 10, 50, 200).

## 7. Is it better?

**Not on accuracy.** One reading, 4000 series per cell, ratio of means
within ×2 (%) (`presets/sept21_vs_sept24_single.yaml`):

| cell | 21 Sept | 24 Sept |
|---|---|---|
| 3D, $m=s/4$, SNR 2 | 34.3 | 37.9 |
| 3D, $m=s$, SNR 1 | 79.5 | 74.7 |
| 1D, $m=s/4$, SNR 1 | 19.1 | 14.8 |
| 1D, $m=s$, SNR 1 | 77.8 | 82.4 |

- Which one wins at weak signal depends on whether the true mass sits near
  $s$. That is a question about the mass prior, not the method.
- From SNR 4 upward they agree within about one point.
- Series by series they differ by about 5% (median) and at most 20%.
- The one edge is the worst case: the largest regret across all cells is
  5.6 / 6.9 / 8.5 points (24 Sept) against 11.7 / 8.9 / 8.8 (21 Sept), at
  ×1.25 / ×1.5 / ×2.

**On consistency, yes.** It never gives two answers to one measurement.
**On form, yes.** It is a Gaussian in every dimension, readings add
residuals, and the interval width has an exact formula. The change is
adopted because it is forced, not because it scores higher.

## 8. How we got here

1. The user asked for the newest equation to be implemented and tested
   against the 21 September one, to isolate the difference. The tester
   commits the previous record cites (`cd53203`, `0c3235f`, `2128823`,
   `9f27343`) are not on GitHub, so it was reimplemented independently
   (§10).
2. The identity of §2 isolated the difference as one weight. Its effect,
   the pull toward $s$, explained the previous record's table.
3. Splitting the weight into the forced $\alpha^{d-2}$ and the completion
   showed two large pieces cancelling. The σ-free completion suggested by
   the 23 September direction was found and tested; it behaves badly.
4. Noticing that the new law of mass does not depend on $d$ led to axis
   splitting as a requirement. It forces the completion (§4).
5. The user asked whether the new equation is better (§7), how the equation
   is used in 1D (a cart on a track: signed readings, direction $\pm1$),
   and whether the property survives nonzero y and z components (§4a). He
   judged axis splitting too elementary to be wrong.

## 9. What changes, what does not, what is owed

- **Proposed:** amend the 21 September specification in 1D and 3D to the
  24 September law, with axis splitting stated as the reason for both the
  weighting and its completion. That text lives on VD-Newton's surface; the
  call is the user's.
- **Unchanged:** 2D; Newton II as the premise; the uniform direction; the
  ratio-of-means readout; the symmetric rule.
- **Owed:**
  - the 2D mass factor itself (investigation question 2), which now decides
    weak-signal accuracy in every dimension;
  - the interval width for many readings: either the finite-temperature
    (hierarchical) prior, or widening by the estimated $\sqrt{J/H}$;
  - an account of the readout's split-dependence (§4);
  - whether axis splitting extends to correlated or unequal per-axis noise
    (the tester's 24 September law covers independent isotropic noise only);
  - a second seed; real data.

## 10. Reproduction

Tester `canvalidk/Mass-estimator-testor`, three commits on top of `cfd6714`
(`3760f6f`, `48980d3`, `fb9b78b`). They were handed over as a bundle and are
not yet pushed. The runs are at `48980d3`. `reference: cartesian` is the 24
September law; `angle_power` expresses the other completions.

```sh
python -m pytest -q tests/test_cartesian.py
python -m met run presets/sept21_vs_sept24_single.yaml
python -m met run presets/sept24_completions_3d.yaml
python -m met run presets/sept24_completions_1d.yaml
python -m met run presets/sept21_vs_sept24_repeated.yaml
```

The large-N predictions of §6 come from a short quadrature script
(noncentral chi expectations of $G'$ and $G''$), not yet in the tester.
The tester's `flat` reference agrees with the standalone
`mass_estimator_equation/mass_estimator.py` to $2\times10^{-10}$ in 1D, 2D
and 3D.
