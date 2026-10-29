# The mass estimator equation: what it is, and why each part is there

**Note:** `_note_mass_estimator_equation`, draft d3, 2026-10-29. To be filed in VD-docs as
`Newton-analysis/_note_mass_estimator_equation/the_mass_estimator_equation_d3.md`, beside d2.
Written by Claude at the user's request when he ruled the switch to ours30: "That involves all the
stuff, writing the python, making an explainer document for the algebraic equation, setting up
the file space." The equation, the premises and the rulings are taken from the tester's living
layer (`canvalidk/Mass-estimator-testor`: `definition/equation.md`, `definition/premises.md`,
`rulings.md`, `problems/`), which is where the work is done; this note explains, it decides
nothing.

**What this is.** The equation the tester now treats as the definition of inertial mass,
**ours30**, stated once in full, then written out for the common case where it can be done by
hand, then justified part by part, with a status on every part:

- *ruled*: the user's decision, with its ruling ID in the tester's `rulings.md`;
- *forced*: follows from a named premise or ruling;
- *provisional*: the best sort so far, open to his check against the goal;
- *input*: not the function's to decide; the caller declares it, as a branch is declared to
  the complex logarithm (R-2026-10-24-3);
- *open*: not settled; the item is named.

**What this is not.**

- Not the argument that an equation is needed at all (`_note_why_mass_estimation`).
- Not the lemma list. The premises and theorems live in the tester's `definition/`; this note
  cites them.
- Not the story of how ours30 was found. That is in the tester, `problems/O13.md`.
- Not the uniqueness proof. That is the open half of the deliverable (§7).

**Relation to d2.** d2 (26 September) stated ours25b. Five equation states have come since:

| state | date | what changed |
| --- | --- | --- |
| ours26 | 26 Sept | per-axis noise allowed |
| ours27 | 27 Sept | any error covariance per reading |
| ours28 | 28 Sept | a per-reading size law (withdrawn 30 Sept, R-2026-09-30-6) |
| ours29 | 30 Sept | each reading weighted by its own misfit; the series' law raised to the Bartlett power η |
| **ours30** | **29 Oct** | **the misfit is taken at the reading's best positive mass, with numerator k − ½** |

Along the way three other things changed that d2 stated differently: the point is now the
law's **median**, not the ratio of means (28 Oct, provisional); composition (lumping, common
force) is a commitment about the **law**, not the point (R-2026-10-28-1, provisional); and the
strong form of axis splitting, which d2's many-reading step used, was ruled false
(R-2026-09-29-9). d2's §8 (no always-finite interval can guarantee coverage) stands.

**Direction.** Not yet ruled (owed for d2 as well).

---

## 0. Where this note starts

One vector cannot be divided by another, and independently measured force and acceleration
never line up exactly. So mass has to be estimated, and the estimator **is** the definition:
there is no mass behind it (F1, R-2026-09-27-1). The task in the tester is to find the function
that is **unique**: everything the readings cannot decide is moved into the input, the way the
branch is an input to the complex logarithm, and what is left inside the function must be forced
(R-2026-10-24-3; TASK.md). ours30 is the current answer to "what is left inside".

## 1. Types

**Inputs, per reading i** (any number of readings of one object, sharing one mass):

- measured force $\tilde{\mathbf F}_i$ and acceleration $\tilde{\mathbf a}_i$ in $\mathbb R^d$,
  $d\in\{1,2,3\}$, any values, zero included; write $\mathbf y_i=(\tilde{\mathbf F}_i;\tilde{\mathbf a}_i)$;
- an error covariance $\Sigma_i$ ($2d\times2d$, positive definite: every declared $\sigma>0$,
  any correlation, also between the channels);
- a candidate subspace $N_i$ ($d\times k_i$, orthonormal columns): the directions the true
  signal could have taken; $N_i=I_d$ when nothing is known, so $k_i=d$;
- a window: the time interval both channels are averaged over (carried, not used in the formula).

**Declared once:** the starting law $\pi$, half-Cauchy with centre $c$,
$\pi(m)=c/(c^2+m^2)$ on $m>0$. By default $c=s=\sigma_F/\sigma_a$ for one noise pair.

**Output:** a law $p(m)$ on $m>0$; the point $\widehat m$ (its median); the 95% interval
(quoted width, η applied) and the raw interval (η = 1) beside it; η; each reading's weight
$g_i$ and misfit $Q_{+,i}$. Never "undefined" on measured inputs. The law has no mean (its tail
is the half-Cauchy's), so "mass ± SD" is not an output.

**Not this function:** exact inputs ($\sigma=0$) go to a three-way check (§2.6).

## 2. The equation, ours30

### 2.1 The explained energy of each reading

Put every candidate slope on the projective line, $m=c\tan\phi$, $\phi\in(-\pi/2,\pi/2]$
(negative and infinite $m$ included). For reading $i$:

$$
C(\phi)=\begin{pmatrix}c\sin\phi\,N_i\\ \cos\phi\,N_i\end{pmatrix},\qquad
G(\phi)=C^\top\Sigma_i^{-1}C,\qquad \mathbf b(\phi)=C^\top\Sigma_i^{-1}\mathbf y_i,\qquad
H_i(\phi)=\mathbf b^\top G^{-1}\mathbf b .
$$

$C(\phi)$ maps a true acceleration $\alpha\in\mathbb R^{k_i}$ in the subspace to the true pair
$(m\,N_i\alpha;\,N_i\alpha)$, scaled so that every slope is treated alike. $H_i(\phi)$ is the part of
the reading's whitened energy that a true pair on slope $\phi$ can explain, with the true
acceleration fitted freely (no size law: L10). $E_i$ is the whitened energy of $\mathbf y_i$ in the
span of all candidates; $E_i-H_i(\phi)$ is the error the reading needs at slope $\phi$.

### 2.2 The misfit and the reading's weight

$$
Q_{+,i}=E_i-\max_{0\le\phi\le\pi/2}H_i(\phi),\qquad
\kappa_i=k_i-\tfrac12,\qquad
g_i=\min\!\Big(1,\ \frac{\kappa_i}{Q_{+,i}}\Big)\quad(g_i=1\text{ if }Q_{+,i}=0),
$$

$$
\ell_i(\phi)=\tfrac12\,g_i\,H_i(\phi).
$$

$Q_{+,i}$ is the smallest error that reconciles the reading with Newton II **with a positive
mass**: the best slope is searched over $m\in[0,\infty]$ only, the two ends as limits. $\kappa_i$
is the most misfit the declared noise produces in expectation. $g_i$ is the share of the
reading's self-contradiction the declared noise can account for, capped at the whole. **This is
the one change from ours29**, which had $Q_i=E_i-\max_{\text{all }\phi}H_i$ and $\nu_i=k_i-1$.

### 2.3 The Bartlett power

At the series' best line on the **whole** projective line,
$\hat\phi=\arg\max_\phi\sum_i\ell_i(\phi)$:

$$
\mathcal H=-\sum_i\ell_i''(\hat\phi),\qquad J=\sum_i\ell_i'(\hat\phi)^2,\qquad
\eta=\min\!\Big(1,\ \frac{\mathcal H}{J}\Big)\quad(\eta=1\text{ if }J=0\text{ or }\mathcal H\le0).
$$

$\mathcal H$ is how sharply the readings, taken together, pin the slope; $J$ is how much they
disagree about it. If the declared noise were the whole story, the two would match on average
(Bartlett's identity); when the readings disagree more than the noise allows, $\eta<1$ and the law
widens. One η per series, because $m$ is one number.

### 2.4 The law

$$
p(m)\ \propto\ \pi(m)\,\exp\!\Big(\eta\sum_i\ell_i(m)\Big)\qquad\text{on }m>0,
$$

with $\ell_i(m)=\ell_i(\phi)$ at $\phi=\arctan(m/c)$.

### 2.5 The point and the interval

$$
\widehat m=\text{the median of }p,\qquad
[\,m_{2.5\%},\ m_{97.5\%}\,]\ \text{of }p\ \text{(quoted)},\qquad
\text{the same of }p\text{ with }\eta=1\ \text{(raw)}.
$$

### 2.6 Exact inputs

When $\sigma=0$, compute the two pairings $m_1=\mathbf F\cdot\mathbf a/\lvert\mathbf a\rvert^2$ and
$m_2=\lvert\mathbf F\rvert^2/(\mathbf F\cdot\mathbf a)$. **One mass** if $m_1=m_2\in(0,\infty)$;
**no mass** if they disagree (the vectors are not parallel) or agree outside $(0,\infty)$ (anti-parallel,
or one vector zero); **every mass** at exact $0/0$, the only place "undefined" survives.

## 3. The common case written out: isotropic noise, full subspace

Take one noise pair per reading, $\Sigma_i=\mathrm{diag}(\sigma_F^2 I,\ \sigma_a^2 I)$, no restriction on
direction ($k_i=d$), and $c=s=\sigma_F/\sigma_a$. Work in noise units,
$\mathbf x_i=\tilde{\mathbf F}_i/\sigma_F$, $\mathbf y_i=\tilde{\mathbf a}_i/\sigma_a$, and write the slope as
$t=m/s=\tan\phi$. Three numbers per reading carry everything:

$$
A_i=\lvert\mathbf x_i\rvert^2,\qquad B_i=\lvert\mathbf y_i\rvert^2,\qquad C_i=\mathbf x_i\cdot\mathbf y_i .
$$

**The error a reading needs at slope $t$.** A true pair on slope $t$ is $(t\mathbf u,\mathbf u)$ in noise
units. Fitting $\mathbf u$ gives

$$
K_i(t)=\frac{A_i+t^2B_i-2tC_i}{1+t^2},\qquad H_i(t)=A_i+B_i-K_i(t)=\frac{t^2A_i+2tC_i+B_i}{1+t^2}.
$$

**The misfits.** Over every slope, the smallest error is the smaller eigenvalue of the Gram matrix
$\begin{pmatrix}A_i&C_i\\C_i&B_i\end{pmatrix}$:

$$
Q_i=\tfrac12\big(A_i+B_i-\Delta_i\big),\qquad \Delta_i=\sqrt{(A_i-B_i)^2+4C_i^2}.
$$

Over positive slopes only (T32):

$$
Q_{+,i}=\begin{cases}Q_i,& C_i>0\ \text{(the best slope is already positive)},\\[2pt]
\min(A_i,B_i),& C_i\le0\ \text{(the best positive slope is } t\to0\text{ or }t\to\infty\text{)}.\end{cases}
$$

So the two misfits differ exactly when force and acceleration point more than 90° apart. Then
the cheapest way to keep $m>0$ is to declare the smaller channel's whole reading an error.

**The weight.** $g_i=\min\big(1,(d-\tfrac12)/Q_{+,i}\big)$.

**The Bartlett power in closed form.** With $t=\tan\theta$,
$H_i=\tfrac12(A_i+B_i)+\tfrac12\,\mathrm{Re}\big(Z_ie^{-2\mathrm i\theta}\big)$, where
$Z_i=(B_i-A_i)+2\mathrm iC_i$. Then

$$
\hat\theta=\tfrac12\arg\Big(\sum_ig_iZ_i\Big),\qquad
\mathcal H=\Big\lvert\sum_ig_iZ_i\Big\rvert,\qquad
\ell_i'(\hat\theta)=-\tfrac12\,g_i\,\mathrm{Im}\big(\bar Z_ie^{2\mathrm i\hat\theta}\big),
$$

and $\eta=\min(1,\mathcal H/J)$ with $J=\sum_i\ell_i'(\hat\theta)^2$. (In the tester's code the same
quantities carry the names P, Q, D for $A$, $B$, $C$.)

**The law.**

$$
p(m)\ \propto\ \frac{s}{s^2+m^2}\ \exp\!\Big(\frac\eta2\sum_ig_i\,\frac{t^2A_i+2tC_i+B_i}{1+t^2}\Big),
\qquad t=\frac ms,\ m>0 .
$$

Everything in ours30 for this case is in that line, plus the median and the two quantiles.

## 4. Why each part is there

| part | what puts it there | status |
| --- | --- | --- |
| Mass is the output of an estimator; nothing behind it | F1 | ruled, R-2026-09-27-1 |
| Newton II assumed true; mismatch is measurement error; no rejection gate | "see what must be true if the law is correct" | ruled, R-2026-09-07-1 |
| σ > 0 in the type; never undefined on measured input | L2; the all-zero experiment returns a law and a point | ruled, R-2026-09-25-1, R-2026-09-29-2 |
| Σᵢ, the reading structure, the window | what the instrument declares; the readings cannot decide them | input; R-2026-09-28-1, R-2026-09-28-7, R-2026-10-10-4 |
| The subspace Nᵢ | the same | input, proposed (L18's second clause: ruling owed) |
| Gaussian noise model | forced along the line by Z (T23); the data need not be Gaussian | forced; R-2026-09-27-12; L13's wording owed |
| Each reading's factor $e^{H_i/2}$, the true acceleration fitted freely | L8's conditional form with L9 forces the Gaussian candidate family; no size law is declared (L10), which gives its flat endpoint (T6) | forced; L10 ruled R-2026-09-30-6 |
| Readings joined only through m | no borrowing between readings (L19) | ruled, R-2026-09-28-6 |
| The misfit is $Q_+$, at the best **positive** mass | the smallest error under the law as assumed, $m>0$ included (R-2026-09-07-1); "misfit at the best mass is evidence about the reading" (L21) | ruled, R-2026-10-29-1 |
| The weight's numerator $\kappa=k-\tfrac12$ | the declared noise's expected misfit (O12's meaning) read without a size (L10): the most it produces over every allowed size, $k-1$ for L8's conditions and ½ for L1's sign | ruled, R-2026-10-29-1; the value is numerical, its proof is owed (tester brief B02); under a general covariance, open |
| The form $\min(1,\kappa/Q_+)$ | a share, at every Q, capped at the whole | forced by the meaning; ruled, R-2026-09-30-10 |
| η = min(1, ℋ/J) | the width must track how far the point moves under repeats (property N); one η because m is one number | forced in the long-series limit; property N ruled R-2026-09-29-4 |
| J without a Bessel factor | an exact zero reading has $\ell'=\ell''=0$ and must leave η alone (Z) | forced from Z, ruled R-2026-09-27-5 |
| η at the whole-line best fit | the ratio ℋ/J is parametrisation-free only where the score is zero; on $m>0$ that point can sit on a boundary | forced (T31); checked against R-2026-09-07-1 on 29 Oct |
| η's cap at 1; the power form at finite N | natural, not forced | open (O14) |
| The law on $m>0$ | Newton II with a positive mass (L1) | ruled (bedrock) |
| π half-Cauchy | symmetric under $m/c\leftrightarrow c/m$, unit-covariant, declared before data; the shape among those is free | input; its shape open (O6) |
| π's centre $c=s$ | the principal branch for one noise pair | input; L3′ not yet ruled (O11) |
| The point is the median | composition belongs to the law, not the point (R-2026-10-28-1); a point that does not depend on the coordinate the law is read in (L23) is a quantile, and the swap picks the median (T35) | provisional; L23 is proposed, not ruled |
| Quoted interval with η, raw beside it | the standard is the quoted width; report both for now | ruled, R-2026-09-30-7 |
| All-zero readings return $c$ | Z, units and the swap | derived, not ruled |
| The exact branch | Newton II on exact vectors, in pairing form | T1 (10 Oct) |

**On the change, in one paragraph.** ours29 measured each reading's misfit over every slope,
negative mass included. So a reading whose force points against its acceleration needed no error
at all (a negative mass fits it) and kept full weight, and its confidence grew without bound as
its noise shrank: the law ran to $m\to0$ or $m\to\infty$. That contradicted the 7 September
ruling, which already said what the misfit is: assume the law, find the smallest error. The law
includes $m>0$. The 7 September record had even computed the right number on its own example:
$\mathbf F=(-\sqrt6,0)$, $\mathbf a=(\sqrt6,0)$, unit noise, "the infimum discrepancy is 6, approached by
keeping one channel and moving the other toward zero". That is $Q_+=6$; ours29's $Q$ is 0 there.
The numerator then changes for a reason of its own: once the sign is counted, the noise can
produce more misfit than before. When one channel carries no signal, noise decides the sign half
the time, and a wrong sign costs a whole extra unit, so the most the noise produces is
$k-\tfrac12$, not $k-1$ (O13, 29 Oct).

## 5. Worked examples

All numbers are from the tester's `met/full.py` (`ours30`, `ours29`) as of 29 October 2026;
c = 1 unless stated.

**(a) One 2D reading, anti-parallel:** $\mathbf F=(1,0)$, $\mathbf a=(-1,0)$, equal SDs.

| SD | | Q | Q₊ | g | 95% interval | evidence span |
| --- | --- | --- | --- | --- | --- | --- |
| 0.1 | ours29 | 0 | 100 | 1 | [0.0005, 1950] | 50 nats |
| 0.1 | ours30 | 0 | 100 | 0.015 | [0.026, 39.2] | 0.75 nats |
| 0.01 | ours29 | 0 | 10⁴ | 1 | [5·10⁻⁶, 2·10⁵] | 5000 nats |
| 0.01 | ours30 | 0 | 10⁴ | 0.00015 | [0.026, 39.2] | 0.75 nats |

Under ours29 the law escapes to both ends as the noise shrinks. Under ours30 the reading stays
worth a fixed 0.75 nats: evidence, just less (L21).

**(b) One 1D reading, anti-parallel:** F = 2, a = −1, SD 0.1. ours29 left every 1D reading at
g = 1, so it returned a point of **291** with interval [55, 7960]: a confident, enormous mass from
a reading whose force opposes its acceleration. ours30: Q₊ = 100, g = 0.005, point **1.59**,
interval [0.045, 47.9], close to the starting law.

**(c) One 2D reading at 30°:** F = (2, 0), a = (cos 30°, sin 30°), SD 0.1. Acute, so Q₊ = Q = 20.9.
ours29: g = 0.048, point 2.143, [0.83, 14.9]. ours30: g = 0.072, point 2.177, [1.00, 9.97]. An
acute reading is discounted less, by the factor (k − ½)/(k − 1) = 1.5 here.

**(d) Four 3D readings that disagree about m** (the tester's STATE.md example; σ_F = 0.2,
σ_a = 0.1, c = s = 2): F = (4, 0.3, 0), (6.1, −0.2, 0.1), (1, 0.8, 0), (3, 0.1, −0.2);
a = (2, 0, 0.1), (2, 0.1, 0), (0.6, −0.4, 0), (1, 0, 0). ours30: point **2.6451**,
[2.0911, 3.4004], η = 0.12, the third reading weighted 0.079. ours29: 2.6454, [2.0914, 3.4009],
the third reading 0.064. All four readings are acute, so only the discount changes.

**(e) The all-zero experiment:** three 3D readings of exactly zero, σ_F = 0.2, σ_a = 0.1. The
point is exactly c: 0.700000 at c = 0.7, 2.000000 at c = 2. Every g is 1 (Q₊ = 0).

**(f) Real data** (Skrede's robot dataset, Zenodo 10.5281/zenodo.11096791; tester
`tools/robot_run.py`): 1-baseline, 0.2 s windows. ours30: median 0.9213, [0.9145, 0.9282];
ours29: 0.9213, [0.9137, 0.9290]. The static calibration mass is 0.9323, outside both. That miss is
about the declared noise budget of that dataset, not the equation (tester `problems/O16.md`,
28 Oct): its noise per window was declared as if samples were independent, and they are not.

## 6. What the switch moves, and what it does not

- **Moves:** readings whose force and acceleration are more than 90° apart now get a bounded,
  never-zero weight, in 1D too. Acute readings are discounted less (1.5× at k = 2, 1.25× at k = 3);
  under the declared model at strong signal that raises efficiency from 0.915 to 0.954 (k = 2) and
  from 0.929 to 0.954 (k = 3), and the share of ordinary readings weighted below 1 falls from
  0.32 to 0.22 and from 0.37 to 0.29.
- **Does not move:** η's rule, the law's form, the starting law, the point's rule, the exact
  branch. Where every reading is acute with misfit at most k − 1 (so g = 1 under both), ours30 and
  ours29 are the same number.
- **Kept:** ours29 stays callable in the tester (`ours29()`, and the batch pipeline without
  `positivity`), so every number recorded under it still reproduces.

## 7. Open

- **The value k − ½** is shown numerically for isotropic channels (a seeded grid over every pair
  of channel strengths: maxima 0.4997, 1.4987, 2.5015 for k = 1, 2, 3); the proof, and the
  companion fact that the sign-blind maximum is exactly k − 1, are queued (brief B02). Under a
  general covariance (correlated channels) the value is not known yet.
- **ours30's coverage** has not been measured; only its efficiency at strong signal.
- **The point** (O2): the median rests on L23, proposed, not ruled.
- **π's shape** (O6) and, with several noise pairs, its centre (O11).
- **η at small N** (O14; brief B04) and the null floor (O15).
- **Errors shared across readings** (O16) have no slot in the type yet.
- **Isotropy of inertia** (L22, O17): a scalar mass is assumed; it is read as a clause of L1.
- **The uniqueness theorem:** that every freedom left is an input, never a choice inside the
  function. The tester's `definition/ladder.md` tracks which commitment removes which freedom.

## Sources

In the tester (`canvalidk/Mass-estimator-testor`):

- `definition/equation.md` (ours30, stated); `definition/premises.md` (L1–L23);
  `definition/theorems/` (T1, T6, T8, T23, T30, T31, T32, T35); `definition/ladder.md`.
- `rulings.md`: R-2026-09-07-1, R-2026-09-29-8, R-2026-09-29-9, R-2026-09-30-10, R-2026-10-24-3,
  R-2026-10-28-1, R-2026-10-29-1 and the others cited in §4.
- `problems/O13.md`, "29 Oct: the numerator sorted"; `problems/O2.md`; `problems/O16.md`.
- `met/full.py` (`ours30`, `ours29`, `reading_state`); `met/analyst.py` (`positivity: true`);
  `tests/test_ours30.py`, `tests/test_sign_numerator.py`; `tools/robot_run.py`,
  `tools/eta_sign_check.py`.

In VD-docs: d2 of this note; `mass_estimator_law_assumed_correct_steering_2026-09-07.md` (the
anti-aligned example); `o12_the_middle_of_the_misfit_weight_2026-09-30.md` (the meaning of the
numerator); `noise_vs_signal/ours29_misfit_weight_and_bartlett_power_2026-09-30.md` (ours29).
