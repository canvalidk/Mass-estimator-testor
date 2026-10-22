# What ours24 and ours25b change in the earlier results

**2026-09-26.** Record, filed at the user's request ("write the record please
explaining all the pieces that changes"), written by Claude from the
conversation of the night of 25–26 September. The user had asked: "check if
the new equation breaks our earlier results. Because I was under the
impression that we had some exact algebraic proofs and results." It follows
`ours25b_the_many_reading_readout_2026-09-25.md` in the same discussion. It
edits nothing and binds nothing.

**Method, and its limit.** Three readers went through the 19 flat records of
this discussion, the 21 September equation (`mass_estimator_equation/README.md`)
and the equation note's d1, and listed every exact, algebraic or closed-form
claim, about 80 in all, with what each depends on. Claude then rechecked
every claim flagged as changed, and the new results below, by computation
(§10). Claude did **not** re-derive every line of every record. Items marked
*not re-derived* rest on the readers' reasoning only.

**In one paragraph.** Two changes are in play, and they damage different
things. The 24 September law change (ours21 → ours24) altered the law in 1D
and 3D, so it moves formulas that were worked out with the old per-dimension
law: a literature identity, a warrant, some closed forms and one worked
example. None of these is a proof a chain step rests on. The 25 September
readout change (ours24 → ours25b) touched only the many-reading readout, so it
breaks exactly one proof: that the readout's acceleration weight adds reading
by reading. That proof has an exact replacement that is simpler than the
original: the whole estimator is carried by four numbers. Everything about the
law, including the reading of λ, $m_0$ and $\sigma_{\ln m}$ as carried state, is
untouched. Two new exact results come with ours25b.

---

## 1. Two changes, two kinds of damage

| change | what it altered | what it can break |
|---|---|---|
| ours21 → ours24 (24 Sept) | the law, in 1D and 3D only (2D identical) | formulas computed with the old 1D or 3D law; nothing in 2D; nothing about the readout's form |
| ours24 → ours25b (25 Sept) | the ratio-of-means readout, for more than one reading only | anything that adds up per-reading acceleration sizes; nothing about the law, the median, the intervals, or one reading |

So a result survives both if it is about the law (or the median and
intervals) in 2D, or about one reading in any dimension with the readout's
form alone.

## 2. The one proof that broke

**What `old_information_in_the_prior_slot_2026-09-24.md` §2.2 proved.** The
readout's acceleration weight splits over readings,

$$
A_{1:N+1}(m)=A_{1:N}(m)+E[\alpha_{N+1}\mid m],
$$

so the exact state carried from old readings to a new one is two curves,
$\big(p_{1:N},A_{1:N}\big)$.

**Why it breaks.** Under ours25b the weight is the expected length of the
stacked vector, $A=\sigma_a\cos\theta\,\mu_{M}\big(H(\theta)\big)$ with
$H^2=\sum_i h_i^2$ and $M=Nd$. A length of a stacked vector is not a sum of
the parts' lengths, so the split fails.

**The exact replacement: four numbers.** Write each reading's von Mises vector
(from `vonmises_state_2026-09-24.md`) and its energy in noise units,

$$
c_i=\frac{(Q_i-P_i)+2i\,D_i}{4},\qquad T_i=P_i+Q_i,\qquad
P_i=|\mathbf x_i|^2,\ Q_i=|\mathbf y_i|^2,\ D_i=\mathbf x_i\cdot\mathbf y_i ,
$$

and add them up: $C=\sum_i c_i$, $T=\sum_i T_i$, $M=Nd$. Then, exactly,

$$
\log p(\theta)=\operatorname{Re}\!\big[C\,e^{-2i\theta}\big]+\text{const},\qquad
H(\theta)^2=\frac{T}{2}+2\operatorname{Re}\!\big[C\,e^{-2i\theta}\big],
$$

so the law and the readout are both functions of $(C,T,M)$ and $s$ alone,
in every dimension. A new reading adds its own $c$, $T$ and $d$. The carried
state went from two curves to four real numbers.

Checked by rebuilding $\widehat m$ from $(C,T,M,s)$ alone and comparing with the
tester's `ratio_of_means_rss` ($\sigma_F=0.7$, $\sigma_a=0.3$, random readings):

| readings | from $(C,T,M,s)$ | tester |
|---|---|---|
| 4 in 3D | 1.7926966274 | 1.7926966274 |
| 6 in 1D | 1.5865390071 | 1.5865390071 |
| 3 in 2D | 2.0698883108 | 2.0698883108 |

**What leaned on the broken half, and what happens to it.**

- `vonmises_state` §2, "what the vector does not carry": it says $C$ cannot
  carry the readout, citing the additive weight. **Reversed:** with $T$ and
  $M$, it can, exactly.
- Contribution 12 (`mass_estimator_regimes_priors_and_updating_2026-09-23.md`),
  quoting contribution 11 §7: "carry a second function for the readout,
  $E[\alpha\mid m]$". **Shrinks:** something extra must still be carried, but
  it is two numbers ($T$, $M$), not a function.
- The tester's sequential rule: its `exact` carry, and the test that it equals
  the symmetric rule on all readings, are built on the additive weight. They
  are correct for ours24. For ours25b the sequential rule is currently refused.
  **Owed:** a `(C, T, M)` carry.
- `old_information` §3 table, ratio-of-means column: numbers (§3).

**What leaned on the other half, and stands.** The same record's §2.1, "the
exact law of $N+1$ readings is the old readings' law times the new reading's
likelihood, so the old law sits exactly in the prior slot", is about the law
and is untouched. So are the things built on it:

- $m_0$ and $\sigma_{\ln m}$ as the geometric mean and log-spread of the carried
  law, starting at $(s,\pi/2)$ and updated by each reading;
- λ as a reparametrisation of $\sigma_{\ln m}$ ($\lambda=0\leftrightarrow\pi/2$,
  $\lambda^\ast=0.35\leftrightarrow1.457$), and the tilt family's
  $\sigma_{\ln m}$ values;
- λ as one reading never taken (the 25 September chain);
- the von Mises combination of the law, exact in every dimension under ours24.

Two cautions that predate ours25b and are unchanged: $(m_0,\sigma_{\ln m})$
were a two-number fit to the old law, never an exact state (about 1% with
strong old readings, 1.6–16% with weak ones; the exact two-number state of the
law is $C$); and neither they nor $p$ alone ever carried the readout.

## 3. Other things ours25b changes

None of these is a proof. Each is a formula or a number stated for the plain
sum over readings (contribution 11, eq. 28).

- **`excitation_weighting_premise_2026-09-24.md` §3**, the readout formula
  $\bar\alpha(m)=\sum_i\sigma_e\,\mu_d(h_i)$: for ours25b it is
  $\sigma_e\,\mu_{Nd}(H)$. Its §8 line "unchanged: the ratio-of-means readout,
  how readings combine" is out of date for more than one reading.
- **`axis_splitting_and_the_24_september_change_2026-09-25.md` §2**, the
  identity $\widehat m_{24}=E_{21}[f\,r^{d-2}]/E_{21}[\alpha\,r^{d-2}]$: written for
  one pair, and still right there. With many readings under ours25b the sizes
  in it are the stacked ones.
- **Many-reading ratio-of-means numbers** in: the 24 September record's
  `reference_repeated` and `excitation_hierarchy` point columns; contribution
  12's N = 50 repetition study (including eq. 18's readout of about $0.08\,s$);
  `old_information` §3. Medians, intervals and coverage in the same tables
  are law-only and unchanged. The size of the change at N = 50 is under 0.06%
  at the median (25 September studies).
- **Contribution 11's combined value 2.008902161134** (and the scaled
  1.997267503687): still true of ours21 with eq. 18. The example uses two
  instruments with different noise ratios, which is outside ours25b's type.
- **`what_forces_the_ratio_of_means_readout.md` §6**, "pool first moments
  $\sum\bar f_i/\sum\bar\alpha_i$ from separately fitted readings": already
  superseded by ours24's joint law; ours25b is a further step away.
- **`mass_estimator_additional_derivations_2026-09-05.md` §3**, the
  associative bookkeeping $(p,pA,pB)$: still true for mixtures of cases of one
  law, which is what it is about. It no longer gives reading-by-reading
  bookkeeping of the readout.

**Not affected:** the composition theorem
(`mass_readout_composition_uniqueness_2026-09-06.md`). It is about one
positive pair, and its §6.4 says pooling across readings is outside it.
ours25b supplies that piece; it does not contradict the theorem.

## 4. What the 24 September law change moved (1D and 3D)

**4.1 The 1D literature identity moves; it does not break.** Note d1 §6 says
the 1D law is Liseo's (2003) eq. (12), the posterior under the Jeffreys prior
for Fieller's problem, restricted to positive ratios. That is true of ours21.
ours24's 1D law, in noise units ($s=1$),

$$
p(m)\;\propto\;\frac{1}{1+m^2}\,\exp\!\left[-\frac{(x-m\,y)^2}{2\,(1+m^2)}\right],
$$

is Liseo's eq. (14): the posterior under the Berger–Bernardo reference prior
$\pi(\lambda,\theta)\propto(1+\theta^2)^{-1/2}$, flat in the nuisance mean
(Liseo 2003, pp. 137–138, in `papers/`), again restricted to positive ratios.
So the 1D law is still a known default, just a different one.

**4.2 The warrant for the flat weighting.** Note d1 §7 cites contribution 6
(the flat measure as the limit of Gaussian conditioning on a small
$|\mathbf F\times\mathbf a|$). That warrants ours21's measure, which agrees with
ours24's only in 2D. ours24's warrant is axis splitting (25 September record).
d1 §5's measure $df\,d\alpha\,d\Omega$ is likewise ours21's.

**4.3 The F/0 and 0/a closed forms** (spec d1 §17/§18;
`mass_estimation_zero_zero_two_cases.md` §2):
$\widehat m\approx\sigma_F\sqrt{2/\pi}/|a|$ and
$|F|/(\sigma_a\sqrt{2/\pi})$. They reproduce ours21 in 1D and hold only as limits
under ours24. Ratio of means, $\sigma_F=\sigma_a=1$:

| d | 0/5: ours21 / ours24 / formula | 5/0: ours21 / ours24 / formula | 0/50 (all) |
|---|---|---|---|
| 1 | 0.1596 / 0.1629 / 0.1596 | 6.2666 / 6.1370 / 6.2666 | 0.0160 |
| 2 | 0.1631 / 0.1631 / 0.1596 | 6.1304 / 6.1304 / 6.2666 | 0.0160 |
| 3 | 0.1670 / 0.1633 / 0.1596 | 5.9897 / 6.1243 / 6.2666 | 0.0160 |

**4.4 The 21 September worked example** (3D, $\mathbf F=(10,0,0)$,
$\mathbf a=(2,0.1,0)$, $\sigma_F=0.5$, $\sigma_a=0.1$). The point is the same to
12 digits under both laws, 4.993753906248, and equals the median. The
interval and log SD move in the 4th digit:

| | 95% interval | log SD |
|---|---|---|
| ours21 (converged) | [4.344240, 5.740238] | 0.071023 |
| ours24 | [4.344624, 5.739731] | 0.070979 |

Why the point does not move: for one reading, the law and the readout's weight
both depend on $\theta$ only through $h(\theta)$, which is symmetric about the
law's mode, the total-least-squares angle. When the law sits far from both
ends of $(0,\pi/2)$, every version's median and ratio of means both land on
$s\tan\theta_{\rm mode}$, up to what the ends cut off. Here that is
$5\,(\sqrt{1+800^2}-1)/800=4.9937539$. The tester's test still pins ours21's
numbers, which remain correct for ours21.

**4.5 Zero–zero.** The intermediate step "at zero readings $f$ and $\alpha$ are
independent half-normals" (`zero_zero` §§3–4) holds under ours24 only in 2D.
The conclusion, $\widehat m=s$, survives in every dimension. So does
$\kappa_0=\Gamma(3/4)^2/\sqrt\pi=0.847213$ (`additional_derivations` §4): it
depends only on the uniform angle, and was checked to 12 digits.

**4.6 Identities about rejected rules.**

- "eq. 18 = product with acceleration as nuisance, flat-mass prior once"
  (`old_information` §2.3; contribution 12 §7): under ours24's `cartesian2`
  split the two differ by $\cos^{N(d-2)}\theta$, so it holds only in 2D (and
  for `cartesian1` in every dimension). The per-reading factors that reject
  eq. 18 and eq. 20 survive.
- Contribution 12's $m^{2N-2}$ count between the two nuisance references
  becomes dimension-dependent under ours24. *Not re-derived.*

**4.7 Formulas specific to ours21's channels.** Katz's known-direction
truncated-normal means (d1 §4) use independent channels, which ours24 couples
in 1D and 3D. *Not re-derived.* The readout audit's measure comparison (flat
against cone volume, under 5% apart, §5) never included ours24's measure.
*Untested.*

**4.8 What ours24 improved.** `vonmises_state`'s leftover factor in 1D and 3D,
and its failures there, are artefacts of ours21; under ours24 the von Mises
form is exact in every dimension (the 24 September premise record already
says so).

## 5. What stands

- **The chain's steps:** Newton II and the three-case check for exact inputs;
  the disagreement rule; σ_F/σ_a at zero readings for any estimator of the
  type; averaging over candidates; the ratio-of-means readout for one pair
  (any two of lumping, common force and symmetry; the 6 September theorem);
  axis splitting and the flat push weighting; the combination rule and its
  convergence argument; λ as a fake reading.
- **Law-level results:** the half-Cauchy null law, its intervals and log SD
  $\pi/2$, the log-moment formulas and infinite mean of $m$; the profile
  $Q(m)=|\mathbf F-m\mathbf a|^2/(\sigma_F^2+m^2\sigma_a^2)$, now the exponent of the
  law; the prior-slot split; the von Mises form and its addition; the width
  factor $1+d/r^{\ast2}$ (per reading; $1+M/R^{\ast2}$ summed).
- **Readout properties:** the force/acceleration swap, unit scaling,
  $\widehat m=E[\alpha m]/E[\alpha]$, zero readings giving $s$ for any N
  under both poolings, the reciprocal-invariance survey (every candidate passes
  it, so it selects nothing).
- **Law-free results:** Gleser–Hwang, Fieller's and Koschat's results, the
  loss characterisation of $B/A$ (`additional_derivations` §1).

## 6. New exact results that come with ours25b

**6.1 Equal energies give exactly $s$.** If the readings carry equal total
energy in the two channels, $\sum_i|\mathbf x_i|^2=\sum_i|\mathbf y_i|^2$ in noise
units, then the ratio of means and the median are both exactly
$\sigma_F/\sigma_a$. This holds for any number of readings, any split of the
energy among them, and any alignment. Reason: with $P=Q$ the swap
$\mathbf x\leftrightarrow\mathbf y$ leaves $(P,Q,D)$ unchanged, so by the
force/acceleration swap $\widehat m/s=s/\widehat m$. It generalises the
zero-reading result, which is the case $P=Q=0$. Under ours24 it holds only
when each reading has equal energies on its own. Checked ($s=7/3$ except the last row, $s=4$):

| readings | $D$ | ours25b | median | ours24 |
|---|---|---|---|---|
| 5 in 1D, anti-aligned | −0.93 | 2.3333333333 | 2.3333333333 | 2.333483 |
| 5 in 3D, anti-aligned | −1.83 | 2.3333333333 | 2.3333333333 | 2.319562 |
| 5 in 2D | −2.53 | 2.3333333333 | 2.3333333333 | 2.345670 |
| 3 in 3D, unequal per reading | +8.46 | 4.0000000000 | 4.0000000000 | 4.011060 |

**6.2 The four-number state** of §2.

## 7. What the equation note's next draft must not copy

- d1 §5: the measure $df\,d\alpha\,d\Omega$ (ours21's; ours24's is flat in the
  push vector).
- d1 §6: "the 1D law is Liseo's eq. (12)". For ours24 it is eq. (14) (§4.1).
- d1 §7: contribution 6 as the warrant (§4.2).
- d1 §10: "several readings of one mass need a declared joint model" is no
  longer open: ours24's law and ours25b's readout supply it.
- `vonmises_state` §2's "does not carry the readout", contribution 12's
  "carry a second function", `old_information` §2.2's additive weight, and
  the excitation premise record's §3 readout formula (§§2–3 above).

## 8. Owed

- The ours25b sequential carry $(C,T,M)$ in the tester, with a test that it
  equals the symmetric rule on all readings.
- The items marked *not re-derived* or *untested* (§4.6–4.7).
- A full line-by-line re-derivation of every record; this audit rechecked
  only what was flagged.

## 9. How we got here

1. After ours25b was filed, the user asked whether the change breaks the
   earlier exact results.
2. Three readers extracted the claims; Claude rechecked the flagged ones and
   found the Liseo move (reading the paper in `papers/`) and the equal-energy
   result (while explaining why the worked example's point did not move).
3. He asked whether it had changed an earlier result (yes, one: §2), whether
   other results depended on it (yes; the law half stands, the readout half
   shrinks), and whether $m_0$ and $\sigma_{\ln m}$ are gone (no: they
   describe the law). Then he asked for this record.

## 10. Reproduction

Tester `canvalidk/Mass-estimator-testor`, commit `e288110` on top of
`19d3b41` (ours25b, committed together with the `stacked_state` and ours26 work). Every number above comes from `met.analyst.estimate` with
`law: {reference: flat | cartesian2, nuisance: radius}` and, for more than one
reading, the symmetric rule (`likelihood_product` with the uniform-angle
prior). The four-number rebuild used Gauss–Legendre quadrature in $\theta$
(6000 nodes) with `met.law.stacked_mean_radius`. $\kappa_0$ was computed by
quadrature over the uniform angle. The scripts were run in Claude's
workspace and are not kept.
