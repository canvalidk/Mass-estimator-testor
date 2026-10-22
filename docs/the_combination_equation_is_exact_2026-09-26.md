# The combination equation is exact: four numbers carry everything

**2026-09-26, written 2026-09-27.** Record, filed at the user's request ("Make a
record of 1. as well please"), written by Claude from the conversation of 26
September. It delivers what
`what_ours24_and_ours25b_change_in_earlier_results_2026-09-26.md` §8 listed as owed,
"the ours25b sequential carry $(C,T,M)$ in the tester, with a test that it equals
the symmetric rule on all readings", and reports the rerun of the old-information
study that it made possible. It edits nothing and binds nothing.

**In one paragraph.** When new readings arrive, what must be kept from the old ones
to get exactly the answer the rule would give on all readings together? Under
ours24's readout the answer was two curves: the old law and the old acceleration
weight. ours25b broke the second, because a stacked length is not a sum of lengths.
The audit replaced it on paper with four real numbers per series. Those numbers
are now built and tested: the old readings' $(C,T,M)$ plus the new reading's gives
the law, the median, every interval and ours25b's readout, all equal to the
all-readings answer to $10^{-14}$ over 21,600 series. Combining is adding, in any
order and grouping. The rerun also measures what the earlier compact carries
lose:
- keeping the old law but not the old weight moves the readout by up to about 6%
  (typical per cell);
- the two-number fits of the law move the median by up to about 19%, and interval
  endpoints by far more.

---

## 1. The equation

**Setting.** ours25b: one $\sigma_F$ and one $\sigma_a$ for the whole series,
$s=\sigma_F/\sigma_a$, either cartesian reference, the radius nuisance, readings
multiplied, the uniform-angle prior counted once. Noise units
$\mathbf x_i=\tilde{\mathbf F}_i/\sigma_F$, $\mathbf y_i=\tilde{\mathbf a}_i/\sigma_a$,
angle $m=s\tan\theta$.

**Each reading's contribution.** For a reading in $d$ dimensions, with
$P_i=|\mathbf x_i|^2$, $Q_i=|\mathbf y_i|^2$ and $D_i=\mathbf x_i\cdot\mathbf y_i$:

$$
c_i=\frac{(Q_i-P_i)+2i\,D_i}{4}\in\mathbb C,\qquad T_i=P_i+Q_i,\qquad M_i=d .
$$

**The state of a set of readings** is the sum of its readings' contributions,
$(C,T,M)=\sum_i(c_i,T_i,M_i)$: one complex and two real numbers, four real numbers
in all.

**The combination equation.**

$$
\boxed{\;(C,T,M)_{A\cup B}\;=\;(C,T,M)_A+(C,T,M)_B\;}
$$

for any two disjoint sets of readings of the same series. The empty set has state
$(0,0,0)$.

**Everything is read off the state.** With

$$
H(\theta)^2=\frac{T}{2}+2\,\mathrm{Re}\!\big[C\,e^{-2i\theta}\big]
\qquad\big(=P\sin^2\theta+Q\cos^2\theta+2D\sin\theta\cos\theta\big),
$$

the law and the readout are

$$
p(\theta)\;\propto\;\exp\!\big(\mathrm{Re}\big[C\,e^{-2i\theta}\big]\big),\qquad
\widehat m\;=\;s\,\frac{\displaystyle\int_0^{\pi/2}\sin\theta\;\mu_M\big(H(\theta)\big)\,p(\theta)\,d\theta}
{\displaystyle\int_0^{\pi/2}\cos\theta\;\mu_M\big(H(\theta)\big)\,p(\theta)\,d\theta},
$$

with $\mu_M$ the mean length of an $M$-dimensional Gaussian vector (the ours25b
record, §1). The law needs only $C$. The readout also needs $T$ and $M$. The prior
is not in the state: it is applied once, when the answer is read off.

**Algebraic facts.**

- The state is a sum, so combining is commutative and associative, and the order
  in which readings arrive, or how they are batched, cannot matter.
- The empty state gives the zero-reading answer: a uniform law in $\theta$ and
  $\widehat m=s$.
- Two states can be added only if they share $s$ and $\sigma_a$, one instrument.
  With per-axis noise (ours26) the state becomes one $(C,T,M)$ per noise pair
  (`ours26_per_axis_noise_2026-09-26.md` §5).

## 2. What it replaces

| carry | what it keeps of the old readings | exact for |
|---|---|---|
| eq. 28's split (`old_information` §2.2) | the old law and the old weight $A_{\rm old}(m)$ as two curves | ours24; **broken under ours25b** |
| `vonmises_state` | $C$ | the law (under either cartesian reference, every $d$) |
| `curve_only` | the old law as a curve | the law |
| `lognormal_fit`, `tilt_fit` | two numbers: the old law's centre and log-spread ($m_0$, $\sigma_{\ln m}$) | nothing: fits |
| **`stacked_state`** | **$(C,T,M)$** | **the law and ours25b's readout, every $d$** |

Four real numbers replace two curves, and carry more. The von Mises record
(`vonmises_state_2026-09-24.md` §2) said $C$ cannot carry the readout; with $T$
and $M$ beside it, it can.

**Credit.** $C$ is a rearrangement of the classical errors-in-variables sums
$S_{xx}$, $S_{yy}$, $S_{xy}$, so carrying the law this way is not new. $T$ and $M$
are there for the readout, which is the project's
(`credit_and_originality_of_the_mass_estimator_2026-09-26.md`).

## 3. Checks in the tester

`combine: {rule: sequential, carry: stacked_state}`. The last reading is "new"; the
others are "old", carried as their state. `tests/test_stacked_state.py` has 26
checks:

- **equal to the rule on all readings:** law, median, 95% interval and
  `ratio_of_means_rss`, relative $10^{-12}$, for 1D, 2D and 3D, $N=2,6,21$, under
  `cartesian1` and `cartesian2`;
- **order and grouping:** adding one reading at a time from the left, from the
  right, or as two batches gives the same $C$, $T$, $M$ ($10^{-13}$);
- **against the per-reading law:** $H^2/2$ equals the summed per-reading
  likelihood up to a constant ($10^{-10}$), and the state's weight equals the
  stacked weight computed from the readings ($10^{-12}$);
- **instruments:** states with different $\sigma_a$ refuse to add;
- **refusals:** a non-cartesian reference, an old rule other than the likelihood
  product, and eq. 28's plain `ratio_of_means`, which the state does not carry.

## 4. The rerun: `presets/old_info_stacked_state.yaml`

**Design.** The old-information studies of 24 September (`old_info_as_prior`,
`old_info_vonmises`) were run under the 21 September law. This reruns them under
ours25b:
- 1D, 2D and 3D;
- masses $s/4$, $s$ and $4s$;
- acceleration SNR 1 and 3;
- $N=2$, 6 and 21, a new push each reading;
- 54 cells × 400 series.

`full` applies the rule to all readings. Each carry is compared with it series by
series. The distance is $|\log(\text{carry}/\text{full})|$, and the table gives
the median over a cell's 400 series, as a range across cells. For intervals it is
the worse endpoint. For every carry except `stacked_state`, the readout is compared
with `full`'s plain `ratio_of_means`, so that gap is exactly what dropping
$A_{\rm old}$ costs.

| carry | readout | median | 95% interval | worst single series |
|---|---|---|---|---|
| **stacked_state** | $\le 9\times10^{-16}$ | $\le 1.2\times10^{-15}$ | $\le 1.7\times10^{-15}$ | $3.9\times10^{-14}$ |
| vonmises_state | 0–0.055 | $\le 1.3\times10^{-15}$ | $\le 1.6\times10^{-15}$ | 0.79 (readout) |
| curve_only | 0–0.055 | $\le 4\times10^{-16}$ | $\le 9\times10^{-16}$ | 0.79 (readout) |
| tilt_fit | 0–0.073 | 0–0.151 | 0–0.51 | 5.06 (interval) |
| lognormal_fit | 0–0.110 | 0–0.172 | 0–0.59 | 4.47 (median) |

**Reading it.**

- **The combination equation is exact** to rounding in every cell, for every
  output. Scores are identical: within-×1.25 differs from `full` by 0.0 points in
  every cell.
- **Keeping the law is not enough for the readout.**
  - `vonmises_state` and `curve_only` carry the law exactly (so their medians and
    intervals match), yet their readout is off by up to 0.055 in $\log m$ (about
    6%) at the cell median, and by up to 0.79 in a single series.
  - The gap is largest with few, weak readings: 0.002–0.055 at SNR 1, 0–0.015 at
    SNR 3, and it falls with $N$.
  - Within-×1.25 scores move by up to 1.5 points.
  - This is the price of dropping $A_{\rm old}$. It is the same for both carries,
    since both keep the law exactly.
- **The two-number fits lose the law itself.**
  - Their medians are off by up to 0.17 in $\log m$ (about 19%). Their interval
    endpoints are off by up to 0.59, a factor of 1.8.
  - Scores move by up to 7.7 points.
  - They are worst at small mass ($m=s/4$), even at SNR 3 (0.165 for the
    lognormal median). They are best at $m=s$, SNR 3, $N=21$ ($4\times10^{-4}$).
  - The tilted family beats the lognormal in most cells (46 of 54 on the median,
    46 on the readout, 49 on the interval), but neither comes close to exact.

**An incidental finding, not about carries.** In the same run, `full`'s 95%
intervals cover as little as 76.5%: 3D, $N=21$, SNR 1, at $m=s/4$ and at $m=s$
(400 series; standard error about 1.1 points). Coverage falls with $N$ and with
weaker signal. At SNR 1 the cell average goes from 95.8% ($N=2$) to 91.1% ($N=6$)
to 86.3% ($N=21$). At SNR 3 it goes from 95.3% to 93.9%. This is the many-reading undercoverage that ours25's sandwich was
built to correct. The sandwich in turn conflicts with axis splitting (the ours25b
record, §8), so calibrated many-reading intervals remain open.

## 5. What changes in earlier records

- **`old_information_in_the_prior_slot_2026-09-24.md`.**
  - §2.2's additive weight is exact for ours24 only, as the audit found. The exact
    ours25b state is now built and tested.
  - Its §3 table was computed under the 21 September law. The rerun above is the
    ours25b version. The carries' ranking is unchanged (exact state, then law-only,
    then tilt, then lognormal); the sizes differ.
- **`vonmises_state_2026-09-24.md` §2** ("what the vector does not carry"): the
  readout is carried once $T$ and $M$ are added, as the audit said. Now tested.
- **Contribution 12** ("carry a second function for the readout"): it is four
  numbers, not a function.
- **The audit's §8 owed item:** delivered.

## 6. Was it a surprise?

When the rerun came back exact, the user asked "wait that worked? did it really?"
It did, but the exactness was not in doubt once the algebra of §1 held. The law and
the weight are exact functions of $(C,T,M)$, so a correct implementation must
match to rounding. What the tester adds is that the code is right: the carry
matches an independent computation from every reading.

The informative part is §4: how much the earlier compact carries were losing. The
answer is a few percent through the readout, and up to about 19% through the fits.

## 7. Owed

- **Per-noise-pair states for ours26:** the equation extends (one $(C,T,M)$ per
  pair), but the tester's carry is not built for it.
- **The sequential rule with ours25's sandwich** is refused. The sandwich
  estimates spread from per-reading scores, which the state as it stands does not
  keep. That is one more place where ours25 does not fit the rest.
- **ours21 (the flat reference) has no such state.** Its per-reading kernel is not
  $e^{h^2/2}$ in 1D and 3D, so the series does not collapse. In 2D it coincides
  with the cartesian references. This form of the combination equation belongs to
  the cartesian references.
- **The many-reading undercoverage** (§4), which belongs in the open-problems list.

## 8. How we got here

1. The audit (26 September) found the one broken proof, the additive weight. It
   replaced it on paper with $(C,T,M)$, checked by rebuilding three answers from
   the four numbers, and listed the tester carry as owed.
2. The user asked: "before we thought the combination equation gave us a
   non-exact expression. Now with ours25b, can we rerun the combination equation
   on the exact version of it?" Then: "Yes exactly. Give that one a go."
3. The carry was built (`StackedState` in `met/law.py`), tested, and run against
   the old carries.
4. On the exact result he asked "wait that worked? did it really?", then "were
   quite far away from getting there right...".
5. On 27 September, asked what in the session was unrecorded, Claude named this
   result first. He asked for this record.

## 9. Reproduction

Tester `canvalidk/Mass-estimator-testor`:
- commit `e288110`, which holds the `stacked_state` work (`ad26e21` in Claude's
  sandbox) together with ours25b and ours26;
- the rerun itself was made in Claude's sandbox on `2af8033` with the
  `stacked_state` changes uncommitted;
- 54 cells × 400 series, 24 minutes.

- `StackedState` and `ReadingSet.stacked_cond_alpha` in `met/law.py`; the
  `stacked_state` carry and `combined_stacked_state` in `met/analyst.py`.

```sh
python -m pytest -q tests/test_stacked_state.py
python -m met run presets/old_info_stacked_state.yaml
```
