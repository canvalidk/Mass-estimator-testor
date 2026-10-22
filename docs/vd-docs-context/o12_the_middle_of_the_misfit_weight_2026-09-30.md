# O12: the middle of the misfit weight

**2026-09-30 (late afternoon).** Record, written by Claude after his "Ok do
that too" on O12. The arguments and numbers are Claude's; the premises are
his rulings (L21, T30's reading of $\kappa$, R1). It offers a closure of O12
and **retracts** the reading Claude gave of the ours29 study earlier the same
day ("the data say the threshold is set too low").

**What it records.** The misfit weight $g(Q)$ had two forced ends ($g(0)=1$;
$g\sim\nu/Q$ in the tail) and a free middle. The middle is not free: T30
already fixed $\kappa=\nu$ **by meaning**, reading $\nu/Q$ as the share of a
reading's self-contradiction that the declared noise accounts for. A
quantity cannot mean that only asymptotically. Read at every $Q$, and capped
at 1 because noise cannot account for more than all of it, it gives
$g=\min(1,\nu/Q)$: the kink. It is the same rule as $\eta$ one level down
(expected disagreement over observed, capped at 1), and the same convention
as the Particle Data Group's scale factor. Its cost on ordinary data is real
and measured: about 7% extra width at strong signal in 2D and 3D, coverage
94–96%, accuracy unchanged. That is the cost of a mean-based share, not an
error; the "threshold too low" reading measured ours29 against ours27, which
is not the standard.

**Relations.** Answers O12 in `../mass_estimator_lemmas_2026-09-29.md`
(T30, L21, R1, T31, O14). Scopes `readings_are_evidence_just_less_2026-09-29.md`
§§2, 5 ("the middle is free") and `ours29_misfit_weight_and_bartlett_power_2026-09-30.md`
§§10, 13 (the "evidence on O12" paragraph, retracted here). Binds nothing.

---

## 1. What was open

T30: $g(0)=1$ (a reading with no misfit keeps full weight); bounded
confidence on a misaligned pair as $\sigma\to0$, with evidence "less, not
none", forces a tail $g\sim\kappa/Q$; and $\kappa=\nu=k_i-1$ was fixed "by
meaning, not by a property": $(k-1)/Q$ is the share of the reading's
self-contradiction the declared noise accounts for (each of the $k-1$
codirectionality conditions costs one noise unit of residual on average). The
shape between $Q\approx0$ and $Q\gg\nu$ was left free (O12), and ours29
chose the kink.

## 2. The meaning forces the middle

The argument that fixed $\kappa$ says what $\nu/Q$ **is**: the declared
noise's expected misfit over the observed misfit. That is a statement about
the ratio at a given $Q$, not about its behaviour as $Q\to\infty$. If $g$ is
"how much of the reading's contradiction its instrument vouches for", then
at every $Q$

$$
g(Q)=\min\!\Big(1,\ \frac{\nu}{Q}\Big),
$$

the cap because a share cannot exceed the whole. Any other middle (the
smooth $1/(1+Q/\nu)$, a threshold at a quantile $q_\alpha$ of $\chi^2_\nu$, the
shifted kink $\nu/(Q-t+\nu)$) keeps the meaning in the tail and gives it up
in the middle, so under it $\nu/Q$ would mean two things. The quantile and
shifted forms also add a parameter ($\alpha$ or $t$) that nothing supplies,
which under F2 is a family of rival definitions.

**It is the same rule as $\eta$.** $\eta=\min(1,\mathcal H/J)$ is the
series-level version: $\mathcal H$ is what the readings' disagreement about
$m$ should be if the declared noise is right (Bartlett), $J$ is what it is.
$g$ is the reading-level version: $\nu$ is what the reading's internal
disagreement should be, $Q$ is what it is. Both cap at 1. And both are the
PDG convention: the scale factor $S=\sqrt{\chi^2/\text{dof}}$, applied
whenever $S>1$, i.e. with the threshold at the mean.

## 3. The cost, measured, and why it is not an error

Under the model (R1: declared $\sigma$ correct), at strong signal $Q$ is
$\chi^2_\nu$ and independent of the reading's score along the best line. A
weight $g(Q)$ then costs efficiency $E[g]^2/E[g^2]$:

| $\nu$ | form | $E[g]$ | efficiency | share of ordinary readings weighted below 1 |
|---|---|---|---|---|
| 1 (2D) | **kink at $\nu$** | 0.849 | 0.915 | 32% |
| 1 | smooth $1/(1+Q/\nu)$ | 0.656 | 0.860 | all |
| 1 | kink at the 95% quantile | 0.987 | 0.995 | 5% |
| 2 (3D) | **kink at $\nu$** | 0.852 | 0.929 | 37% |
| 2 | smooth | 0.596 | 0.881 | all |
| 2 | kink at the 95% quantile | 0.989 | 0.997 | 5% |

In the tester study (signal 3), ours29's raw law is 5–9% wider than ours27's
in 2D and 3D (1D: identical), covers 94–96%, and its points are as accurate
(within 0.3 points). So the kink is **slightly conservative on ordinary
data**, never overconfident.

**Why the earlier reading was wrong.** Claude wrote that "a kink at the mean
of $\chi^2_\nu$ is not flat inside the noise band, only flat inside half of
it", and called it evidence the threshold was too low. The comparison was
with ours27, whose readouts it moves by about 1.5%. ours27 is not the
standard: under ruling 3 the standard is coverage, which the kink meets.
What the table shows is the price of a mean-based share with few degrees of
freedom (one or two per reading, against about $N$ for $\eta$), the same
price the PDG pays when it applies $S$ to a handful of measurements.
Retracted.

**The cost could be recovered, and that is O14's business.** With $g$ in
place the per-reading Bartlett ratio at strong signal is
$E[g]/E[g^2]\approx1.09$ (3D): the weighted law is underconfident by about
that factor, which is exactly what an uncapped $\eta$ would read and undo
(width $\times1/\sqrt{1.09}\approx0.96$). The cap at 1 (natural, not forced;
T31) is what keeps the cost. Arithmetic, not tested; recorded under O14.

## 4. Offered closure (not placed)

- **T30:** the middle is forced by the meaning that fixed $\kappa$:
  $g=\min(1,\nu/Q)$ at every $Q$. Status "forced by meaning" for the whole
  function, not only $\kappa$.
- **O12:** close. Note the cost (≈7% width at strong signal in 2D/3D, 94–96%
  coverage) and its one possible remedy, the cap on $\eta$ (O14).
- **ours29 record §10:** the "evidence on O12" paragraph is retracted by
  this record (§3); not edited there.
- **§4 of the lemma record (declared shapes):** $g$'s middle leaves the
  table; the only declared shape left is $\pi$'s.

## 5. Checks

Sandbox, scipy: $E[g]$ and $E[g^2]$ by quadrature against the $\chi^2_\nu$
density; share weighted below 1 is $P(\chi^2_\nu>\nu)$ for the kink. Widths
and coverage from the tester run of `presets/ours29_noise_vs_signal.yaml`
(ours29 record §9), signal-3 cells: ours29-law over ours27 median log width
1.052–1.085 in 2D and 3D, 1.000 in 1D. Script not kept.

## Sources

- `../mass_estimator_lemmas_2026-09-29.md`: L21, T30, T31, R1, O12, O14, §4.
- `readings_are_evidence_just_less_2026-09-29.md` §§2, 4, 5, 10.4 (the
  meaning of $\kappa$).
- `ours29_misfit_weight_and_bartlett_power_2026-09-30.md` §§2 (the PDG
  precedent), 9, 10, 13.
