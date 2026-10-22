# ours26: the estimator with a different noise on each axis

**2026-09-26.** Record, filed at the user's request ("Ok will you make a record
about ours26"), written by Claude from the conversation of the evening of 26
September. It names ours26 and records how it was built in the tester. It follows
`ours25b_the_many_reading_readout_2026-09-25.md` in the same discussion. It edits
nothing and binds nothing.

**Names used here.**

| name | what it is |
|---|---|
| ours24 | the 24 September law (tester `reference: cartesian2`) with the plain-sum readout |
| ours25b | ours24's law with the stacked-length readout (tester `ratio_of_means_rss`) |
| **ours26** | ours25b applied to readings with a different noise SD on each axis |

**In one paragraph.** Every equation up to ours25b assumed one force SD and one
acceleration SD for every axis. The user had not noticed that assumption and did
not want it ("we don't care about elegance… sticking to the wall"). So the chain
was walked again with a separate $\sigma_{F,k}$, $\sigma_{a,k}$ on each axis $k$,
keeping every move that stays forced. Axis splitting still fixes the law: each axis
is its own one-dimensional reading, so the law is a sum of Beale terms, one per axis.
The readout keeps its definition, the ratio of the stacked vectors' mean lengths.
Its weight is now a *weighted* noncentral chi mean, with no closed form, but its
Laplace transform is closed, so one numerical integral gives it to about $10^{-10}$.
One thing comes off the wall: the prior's centre. With several axis ratios there is
no single instrument ratio. The symmetries that fixed it before now leave a whole
family, so the centre is declared, like the prior's shape. ours26 gives ours25b
exactly when the axes agree. In the tester, a 3D series and its axes agree to
$5\times10^{-11}$, and the open centre moves the answer by up to about 50% on weak
data and under 1% on strong data.

---

## 1. The equation (ours26)

**Types.**

- $N\ge1$ readings of one object in $d\in\{1,2,3\}$ dimensions: measured force
  $\tilde F_{ik}$ and acceleration $\tilde a_{ik}$, reading $i$, axis $k$, any
  values including zero.
- Noise SDs $\sigma_{F,ik}>0$ and $\sigma_{a,ik}>0$, one pair per component,
  independent errors. The tester's worlds use one pair per axis, the same in every
  reading, but the equation allows any.
- A declared prior centre $c>0$ (§4).
- Output: a point $\widehat m\in(0,\infty)$ and a law over $(0,\infty)$.

**Noise pairs.** Group the $M=Nd$ components by their noise pair
$(\sigma_{F,g},\sigma_{a,g})$. In the tester's worlds a group is one axis across
all readings. For each group $g$:

$$
s_g=\frac{\sigma_{F,g}}{\sigma_{a,g}},\qquad
\tan\theta_g=\frac{m}{s_g},\qquad
x_c=\frac{\tilde F_c}{\sigma_{F,g}},\quad y_c=\frac{\tilde a_c}{\sigma_{a,g}}\quad(c\in g),
$$

$$
P_g=\sum_{c\in g}x_c^2,\qquad Q_g=\sum_{c\in g}y_c^2,\qquad D_g=\sum_{c\in g}x_c\,y_c,\qquad M_g=\#g .
$$

**The law of mass.**

$$
p(m)\;\propto\;\frac{c}{c^2+m^2}\;
\exp\!\left[-\sum_{i=1}^{N}\sum_{k=1}^{d}
\frac{\big(\tilde F_{ik}-m\,\tilde a_{ik}\big)^2}{2\left(\sigma_{F,ik}^2+m^2\sigma_{a,ik}^2\right)}\right]
$$

or, by noise pair, with one angle per pair,

$$
p(m)\;\propto\;\frac{c}{c^2+m^2}\;\exp\!\Big[\tfrac12\sum_g H_g^2\Big],\qquad
H_g^2=P_g\sin^2\theta_g+Q_g\cos^2\theta_g+2D_g\sin\theta_g\cos\theta_g .
$$

When every axis has the same pair, there is one group, and this is ours24's law
term by term.

**The readout.**

$$
\boxed{\;
\widehat m\;=\;\frac{E\big[\,|\mathbf F^\ast|\,\big]}{E\big[\,|\mathbf a^\ast|\,\big]}
\;=\;\frac{\displaystyle\int_0^\infty m\,A(m)\,p(m)\,dm}{\displaystyle\int_0^\infty A(m)\,p(m)\,dm}\;}
$$

with $\mathbf F^\ast,\mathbf a^\ast\in\mathbb R^{M}$ the stacked true vectors, as in
ours25b, and the weight

$$
A(m)=E\big[\,|\mathbf a^\ast|\;\big|\;m\big]
=E\sqrt{\sum_g \sigma_{a,g}^2\cos^2\theta_g\;Y_g},\qquad
Y_g\sim\chi'^2_{M_g}\!\big(H_g^2\big)\ \text{independent}.
$$

At fixed $m$ each component's true push, in its own noise units, is
$v_c\sim N(w_c,1)$ with $w_c=x_c\sin\theta_g+y_c\cos\theta_g$. The true
acceleration is $\alpha_c=\sigma_{a,g}\cos\theta_g\,v_c$, and a group's squared
pushes add to a noncentral chi-square with $M_g$ degrees of freedom and
noncentrality $H_g^2$. The force side needs no separate weight: within one
candidate $\mathbf F^\ast=m\,\mathbf a^\ast$, since
$\sigma_{F,g}\sin\theta_g=m\,\sigma_{a,g}\cos\theta_g$.

**How $A$ is computed.** With one group,
$A=\sigma_a\cos\theta\,\mu_M(H)$, ours25b's weight. With several groups there is no
closed form, but the Laplace transform of the weighted sum is closed:

$$
A=\frac{1}{2\sqrt\pi}\int_0^\infty\Big[1-\prod_g\big(1+2t\,c_g\big)^{-M_g/2}
\exp\!\Big(-\frac{t\,c_g\,H_g^2}{1+2t\,c_g}\Big)\Big]\,t^{-3/2}\,dt,
\qquad c_g=\sigma_{a,g}^2\cos^2\theta_g ,
$$

from $\sqrt X=\frac{1}{2\sqrt\pi}\int_0^\infty(1-e^{-tX})\,t^{-3/2}\,dt$. The
integral is one trapezoid sum in $\log t$ (§6).

**Where each piece comes from.**

| piece | status under per-axis noise |
|---|---|
| divide average force size by average acceleration size | forced, unchanged (composition) |
| flat weighting over each push, in that axis's own noise units | forced: axis splitting makes each axis its own 1D reading, and the combination rule rules out a per-reading factor in $m$ (§2) |
| readings (and axes) multiply, prior counted once | forced, unchanged (combination rule) |
| size of the series = length of the stacked vector | forced, unchanged (axis splitting on the answer) |
| the weight $A$ as a weighted noncentral chi mean | follows from the above; computed numerically |
| equal belief for equal angles | shape still open; **centre now open too** (§4) |

## 2. The law: why each axis is a Beale term

Axis splitting says a vector reading and the same reading taken axis by axis give
the same mass. With per-axis noise, axis $k$ taken alone is a one-dimensional
reading with its own ratio $s_k$. Under ours24's one-dimensional law its latent
push $v_k$, in its own noise units
($f_k/\sigma_{F,k}=v_k\sin\theta_k$, $\alpha_k/\sigma_{a,k}=v_k\cos\theta_k$), is
weighted flat.

Integrating $v_k$ out leaves

$$
\int \exp\!\Big[-\tfrac12\big(x_k-v\sin\theta_k\big)^2-\tfrac12\big(y_k-v\cos\theta_k\big)^2\Big]dv
\;\propto\;\exp\!\Big[-\tfrac12\big(x_k\cos\theta_k-y_k\sin\theta_k\big)^2\Big]
=\exp\!\Big[-\frac{(\tilde F_k-m\tilde a_k)^2}{2(\sigma_{F,k}^2+m^2\sigma_{a,k}^2)}\Big],
$$

with no factor that depends on $m$. The other ways to weight the push, for example
flat in the physical acceleration $\alpha_k$ at fixed mass, leave a factor
$\cos\theta_k$ per axis, which depends on $m$. The combination rule (nothing per
reading depends on mass) excludes them, as it did for ours24. So the per-axis
law is not a new choice. It is ours24's law, applied to one-dimensional readings
whose instruments differ.

In the tester this is literal: `ReadingSet.reading_terms` on per-axis readings
splits them into their axes and sums the one-dimensional terms, and a test checks
the result against the Beale formula above ($4\times10^{-14}$).

## 3. The readout: nothing new to decide

ours25b's argument for the stacked length used axis splitting on the answer and one
joining rule. It never used the noise. So the readout's definition carries over
unchanged: $E|\mathbf F^\ast|/E|\mathbf a^\ast|$, pooled on each candidate before
averaging. Only its value changes, because the components now weigh differently:
$\sigma_{a,g}\cos\theta_g$ differs between groups, both through $\sigma_{a,g}$ and
through the angle.

The same weighted mean also covers a case ours25b refused: SDs that differ between
the readings of a series (different instruments), even when every axis agrees.

## 4. The prior's centre comes off the wall

**What fixed it before.** With one ratio $s$, the zero-reading argument pinned
the centre at $s$. Any zero-data answer must be a mass built from the SDs, so it
scales with the force unit and against the acceleration unit, which leaves $k\,s$.
Swapping force and acceleration must invert it, which leaves $k=1$.

**What those conditions leave now.** A zero-data centre $c(s_1,\dots,s_d)$ must:

1. scale as a mass: $c(\lambda s)=\lambda\,c(s)$;
2. invert under the swap: $c(1/s)=1/c(s)$;
3. not care about the order of the axes.

These no longer pick one answer. For any $p$,

$$
c_p(s)=\sqrt{\,\mathcal M_p(s)\;\mathcal M_{-p}(s)\,},\qquad
\mathcal M_p(s)=\Big(\tfrac1d\sum_k s_k^{\,p}\Big)^{1/p},
$$

meets all three, because $\mathcal M_p(1/s)=1/\mathcal M_{-p}(s)$. $p=0$ is the
geometric mean. $p\to\infty$ gives $\sqrt{\max_k s_k\,\min_k s_k}$. So the centre
is not forced, and it joins the prior's shape as an open, declared choice.

**Which angle.** "Equal belief for equal angles" also needs one angle. With
several ratios there are several angles $\theta_k$, and a prior uniform in one of
them is not uniform in another. Choosing the centre $c$ is choosing the angle
$\tan\theta=m/c$.

**In the tester.**

- `first_reading` (the old centre) is refused when the first reading's axes have
  different ratios.
- The new `first_reading_geometric` declares $p=0$, the geometric mean of the
  first reading's axis ratios. It equals `first_reading` when the axes agree.
- A number is always allowed.
- After `split_axes`, "the first reading" still means the original one, so the
  centre does not change when a series is split.

**The exception.** If $\sigma_{F,k}/\sigma_{a,k}$ is the same on every axis but the
SDs are not, every axis has the same angle. The centre is forced again, at that
common ratio, and `first_reading` is allowed. The readout is still the weighted
one, because $\sigma_{a,k}$ differs.

## 5. What holds, what is lost

**Kept, and tested.**

- **ours25b exactly on its own ground:** with one SD per channel on every axis,
  ours26 equals ours25b to $10^{-12}$ (1D, 2D, 3D; one and five readings; both
  cartesian implementations). Two old presets were rerun on the old and new code
  with the same data. All 684 saved arrays are bit-identical.
- **Axis splitting, law and readout:** a 3D series with per-axis noise and the same
  series as $3N$ one-dimensional readings give the same answer, with the geometric
  centre or a numeric one.
- **Units:** scaling the force readings and force SDs by $\lambda$ scales
  $\widehat m$ by $\lambda$.
- **The force/acceleration swap** gives $1/\widehat m$ (to $10^{-9}$, the
  quadrature).
- **Only sums matter, per noise pair:** the answer depends on the data only
  through $P_g,Q_g,D_g,M_g$ for each noise pair. Mixing the readings by a rotation,
  separately on each axis, leaves it unchanged to $10^{-10}$.

**Consequence of the last point.** ours25b's exact combination state $(C,T,M)$
survives, one copy per noise pair:

$$
H^2(\theta)\;\to\;\sum_g\Big(\tfrac{T_g}{2}+2\,\mathrm{Re}\big[C_g\,e^{-2i\theta_g}\big]\Big),
$$

and combining sets of readings is still adding states, pair by pair. The tester's
`stacked_state` carry is not built for it yet.

**Lost.**

- **Rotation of the lab frame.** Rotating diagonal per-axis noise makes it
  correlated, which is outside the model. The instrument's axes are physical now,
  so no rotation symmetry is expected.
- **One reading in $Nd$ dimensions.** A series is no longer one big reading. It is
  one block per noise pair, each with its own angle.
- **Zero readings give the centre.** The median still does, since the law of an
  empty series is the prior. The readout does not, because the axes weigh
  differently. With axis ratios 5, 2 and 0.8 and the geometric centre 2, zero
  readings give $\widehat m=1.759$ for one reading, 1.730 for five and 1.720 for
  fifty. Under ours25b the zero-reading readout was exactly $s$ for any $N$. The
  mechanism is the padding drift of the ours25b record (§7): zero readings still
  carry noise-sized unknown pushes, and here the axes weigh them unequally.
- **Not built for per-axis readings:** ours24's plain sum of sizes, laws other
  than the cartesian radius-nuisance product (split the axes first to use them),
  and the `stacked_state` carry.

## 6. Numerics

- **The weight.** The trapezoid rule in $\log t$, step 0.4, 181 nodes centred on
  $t=1/E[X]$, with both tails added in closed form. The integrand is analytic in
  $\log t$ and decays like $e^{-|\log t|/2}$, so the rule converges fast.
- **Accuracy of the weight:**
  - against ours25b's closed form (one group, $M$ from 1 to 150, $h$ from 0 to 300,
    weights from $10^{-3}$ to 50): $\le 10^{-11}$;
  - two groups, against polar-coordinate integration: $6\times10^{-11}$;
  - two groups with weights $10^4$ apart, against nested adaptive quadrature:
    $7\times10^{-11}$.
- **On the posterior grid.** $A$ is computed at nodes no more than 0.02 apart in
  $\log m$ and carried to the 2001-point grid by a cubic spline, about
  $5\times10^{-9}$.
- **The whole readout** for one 2D reading with two ratios, against an independent
  Gauss–Legendre quadrature over $\log m$ with polar integration at each node:
  $3\times10^{-11}$.
- **Cost.** A full run of the preset below (18 cells × 5 estimators × 400 series)
  took 18 minutes.

## 7. First study: `presets/ours26_per_axis_noise.yaml`

**Design.**

- **World:** 3D, a new push per reading, random directions. Force SDs 1 on every
  axis; acceleration SDs 0.2, 0.5 and 1.25, so the axis ratios are 5, 2 and 0.8.
  Acceleration SNR is measured against the RMS of the axes' SDs.
- **Grid:** mass 0.5, 2 and 8; SNR 1, 3 and 10; N = 1 and 5; 400 series per cell.
- **Estimators:**
  - `ours26`: the geometric centre, 2;
  - `centre_low` and `centre_high`: ours26 with the centre moved to the extreme
    axis ratios, 0.8 and 5;
  - `ours26_axes`: ours26 on the split series;
  - the known-excitation oracle as a ceiling.

**Axis splitting.** Across all 18 cells, `ours26_axes` differs from `ours26` by at
most $4.9\times10^{-11}$ in the readout and $1.4\times10^{-14}$ in the median.

**How much the open centre matters.** The table gives the within-×1.25 success
rate of the readout for each estimator, and the typical (median, series by series)
$|\log|$ shift of `centre_low` and `centre_high` from `ours26`.

| $m$ | SNR | N | ours26 | centre 0.8 | centre 5 | oracle | shift, 0.8 | shift, 5 | 95% coverage |
|---|---|---|---|---|---|---|---|---|---|
| 0.5 | 1 | 1 | 0.11 | 0.27 | 0.06 | 0.05 | 0.42 | 0.37 | 0.96 |
| 0.5 | 1 | 5 | 0.23 | 0.27 | 0.21 | 0.30 | 0.20 | 0.09 | 0.91 |
| 0.5 | 3 | 1 | 0.29 | 0.34 | 0.27 | 0.35 | 0.17 | 0.06 | 0.97 |
| 0.5 | 3 | 5 | 0.42 | 0.41 | 0.42 | 0.47 | 0.06 | 0.02 | 0.94 |
| 0.5 | 10 | 1 | 0.60 | 0.60 | 0.60 | 0.61 | 0.03 | 0.01 | 0.96 |
| 0.5 | 10 | 5 | 0.91 | 0.91 | 0.91 | 0.95 | 0.006 | 0.001 | 0.94 |
| 2 | 1 | 1 | 0.27 | 0.16 | 0.27 | 0.33 | 0.36 | 0.31 | 0.98 |
| 2 | 1 | 5 | 0.39 | 0.33 | 0.41 | 0.55 | 0.13 | 0.10 | 0.93 |
| 2 | 3 | 1 | 0.54 | 0.47 | 0.53 | 0.72 | 0.08 | 0.06 | 0.97 |
| 2 | 3 | 5 | 0.88 | 0.89 | 0.88 | 0.99 | 0.01 | 0.01 | 0.96 |
| 2 | 10 | 1 | 0.95 | 0.96 | 0.96 | 1.00 | 0.006 | 0.006 | 0.95 |
| 2 | 10 | 5 | 1.00 | 1.00 | 1.00 | 1.00 | 0.001 | 0.001 | 0.94 |
| 8 | 1 | 1 | 0.36 | 0.35 | 0.42 | 0.84 | 0.04 | 0.10 | 0.96 |
| 8 | 1 | 5 | 0.70 | 0.69 | 0.70 | 0.99 | 0.005 | 0.02 | 0.94 |
| 8 | 3 | 1 | 0.79 | 0.79 | 0.79 | 1.00 | 0.003 | 0.01 | 0.96 |
| 8 | 3 | 5 | 0.99 | 0.99 | 0.99 | 1.00 | 0.000 | 0.002 | 0.97 |
| 8 | 10 | 1 | 0.99 | 0.99 | 0.99 | 1.00 | 0.000 | 0.001 | 0.95 |
| 8 | 10 | 5 | 1.00 | 1.00 | 1.00 | 1.00 | 0.000 | 0.000 | 0.96 |

**Reading it.**

- **The centre matters where the data are weak, as the prior's shape does.**
  - At SNR 1 with one reading, moving the centre from 2 to either extreme shifts the
    typical answer by 0.31–0.42 in $\log m$, a factor of about 1.4–1.5.
  - At SNR 3 with five readings it is under 6%.
  - At SNR 10 it is at most 3%, and under 1% with five readings.
- **Which centre scores best follows the true mass.** At $m=0.5$ the low centre
  wins; at $m=8$ the high one does. That is the prior pulling toward itself, not
  evidence for a centre, because the true mass decides the winner.
- **The intervals** cover 91–98%. The lowest (91%, about three standard errors
  below 95) is at $m=0.5$, SNR 1, N = 5. It has not been looked into.
- **The gap to the oracle** is what it was before per-axis noise, largest at high
  mass and low SNR (0.36 against 0.84).

## 8. What is owed

- **The centre.** It is open, like the shape. Nothing here chooses $p$. A
  principled choice would need a new premise, for example about what zero readings
  should return. As §5 shows, no centre makes the zero-reading readout equal the
  centre for every $N$.
- **The per-pair combination state.** It exists (§5) but the tester does not carry
  it yet.
- **Correlated errors.** Between axes within a channel, a correlation that is the
  same in every reading, with the two channels sharing principal axes, reduces to
  ours26 by rotating to those axes. That is a sketch, not checked. Correlation
  between the force and acceleration channels is not addressed.
- **Mixed instruments** within a series: the equation and the tester's law and
  readout support them, but no tester world generates them yet.
- **The low coverage cell,** and ours25's sandwich width under per-axis noise (not
  tried).
- **A record** of the per-pair state, if the user wants the combination equation
  written for ours26.

## 9. How we got here

1. After the combination-equation rerun, the user asked whether correlated
   errors and different noise ratios are the same thing. Claude answered that they
   are not, and in doing so said that every equation so far assumed one noise pair
   for all axes.
2. He had not known that ("I thought we already allowed different errors on
   different axis"), and asked if it was an easy change.
3. Claude's first answer weighed what the change cost in elegance, since the
   four-number state and the closed-form weight would go. He objected: "we don't
   care about eleganceeee!!!! All we care about sticking to the waaaallllll!!!!"
4. Claude re-sorted the moves by what stays forced: the law and the readout stay
   forced, and the prior's centre comes off the wall.
5. He named the result: "Yes we can call that ours26."
6. It was built in the tester the same evening. Writing this record turned up two
   things that corrected the first build's README:
   - the answer still depends only on per-pair sums, so the state is not lost;
   - the zero-reading readout is not the centre.

## 10. Reproduction

Tester `canvalidk/Mass-estimator-testor`: commit `e288110`. It holds ours26 (two
commits in Claude's sandbox, on top of `ad26e21`) together with the ours25b and
`stacked_state` work.

- **Where it lives.**
  - The law: `ReadingSet.reading_terms` (per-axis SDs, shape (B, N, d)).
  - The weight: `weighted_stacked_mean` and `ReadingSet.stacked_cond_alpha`,
    in `met/law.py`.
  - The centre options, refusals and oracle: `met/analyst.py`.
  - The world setting `noise: {acceleration_sd: {per_axis: [...]}}`:
    `met/world.py`.
- **Tests.** `tests/test_ours26.py` checks against closed forms, polar
  integration, nested quadrature, an independent quadrature of the whole readout,
  and the Beale formula. It also covers the reduction to ours25b, axis splitting,
  per-pair sums, the zero-reading result, units, the swap, per-axis pooling, the
  oracle, the world's draws, and the refusals. All 842 tests pass.

```sh
python -m pytest -q tests/test_ours26.py
python -m met run presets/ours26_per_axis_noise.yaml
```
