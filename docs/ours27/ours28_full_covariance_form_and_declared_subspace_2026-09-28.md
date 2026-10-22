# ours28 in full-covariance form, and the declared candidate subspace

**2026-09-28 (late night).** Record, filed at the user's request ("continue
to work on the solution"), written by Claude after reading the five records
filed in this folder during the evening: the betweenness/fork record and the
four Codex records. It edits nothing and binds nothing. Where it disagrees
with an earlier record of mine it says so; it does not review the others'.

**In one paragraph.** ours28 has a two-line full-covariance form: per
reading, $H_i(m)=T_i-Q_i(m)$ with $Q_i$ ours27's own misfit and $T_i$ the
reading's total whitened energy, and the evidence is
${}_1F_1(\tfrac12;\tfrac k2+1;\tfrac{H_i}{2})$ in place of $e^{H_i/2}$, where
$k$ is the dimension of the reading's **declared candidate subspace**.
Rotation covariance, the zero-reading property and the reduction to the
per-axis form are then automatic (§1). That subspace is a declared input of
the reading, like $\Sigma_i$: $k=d$ says nothing about the push's direction
is known; $k=1$ along $\mathbf n$ says it was along $\mathbf n$. With it,
Codex's heterogeneous counterexample resolves: undeclared, ours27 and ours28
both under-cover (width ratios 1.66 and 1.48, Codex's 1.61 reproduced);
declared, both calibrate (1.01 and 0.95, coverage 94–95%) (§3). The same
declaration is what my known-empty record §3 was missing: the widening from
precise transverse nulls is a solid-angle effect, present in any units,
independent of the transverse instruments' precision, and correct for an
undeclared direction; that record's "physically backwards" verdict and its
"one $\omega$ per noise pair" candidate are withdrawn (§2). On the fork
(grouping bookkeeping or physics, Z, properness) I recommend: grouping is
physics, because a reading is one push; Z and properness are kept; C3 in its
narrow form is folded into L12 as a premise (§4). Calibration is judged in
the definition's own world, which any proper measure passes by theorem;
planted-mass worlds are robustness checks whose failures name an undeclared
input (§5). What ours28 is now, and what is still owed, in §6.

---

## 1. ours28 in full-covariance form

**Per reading.** Data $\mathbf y_i=(\tilde{\mathbf F}_i;\tilde{\mathbf a}_i)\in\mathbb R^{2d}$,
declared error covariance $\Sigma_i$, declared candidate subspace with
orthonormal basis $N_i\in\mathbb R^{d\times k}$ (the push $\mathbf a^\ast_i=N_i\boldsymbol\alpha_i$,
$\boldsymbol\alpha_i\in\mathbb R^k$). Put $C_m=B_mN_i$ with $B_m=(mI_d;I_d)$.

$$
T_i=\mathbf y_i^\top\Sigma_i^{-1}\mathbf y_i,\qquad
H_i(m)=\mathbf y_i^\top\Sigma_i^{-1}C_m\big(C_m^\top\Sigma_i^{-1}C_m\big)^{-1}C_m^\top\Sigma_i^{-1}\mathbf y_i .
$$

$H_i$ is the reading's explained energy: the squared length of the whitened
data projected onto the whitened candidate subspace. For $k=d$,
$T_i-H_i(m)=Q_i(m)$, ours27's Mahalanobis misfit (checked, $3\times10^{-15}$).

**The evidence.** With the push isotropic in the whitened subspace with
scale $\omega_i$, $\beta_i=\omega_i^2/(1+\omega_i^2)$, the reading's
marginal is $N(0,I+\omega_i^2P_m)$ with $P_m$ the projector, so
$\log L_i=-\tfrac k2\log(1+\omega_i^2)-\tfrac12T_i+\tfrac12\beta_iH_i(m)$
(the inverse is $I-\beta_iP_m$, the determinant $(1+\omega_i^2)^k$, both
$m$-free). Averaging $\beta_i$ under the size-angle rule:

$$
\boxed{\;E_i(m)={}_1F_1\!\Big(\tfrac12;\;\tfrac k2+1;\;\tfrac{H_i(m)}{2}\Big),\qquad
p(m)\propto\pi(m)\prod_iE_i(m)\;}
$$

ours27 is $E_i=e^{H_i/2}$ (the $\beta_i\to1$ endpoint). Any other symmetric
size law replaces the $\tfrac12$ by its $a$ (Jeffreys-check record §2).

**Checks (numerical, one 2D reading, random correlated $\Sigma$).**

- Diagonal $\Sigma$: $H=\sum_ch_c^2(\theta_c)$, the per-axis form of the
  ours28 record, to $2\times10^{-16}$.
- Rotation of data and $\Sigma$ together: $H$ unchanged to $2\times10^{-15}$.
- Zero reading: $H\equiv0$, so $E_i$ is constant (Z, exact, any $\Sigma$).
- $T-H$ against ours27's $\mathbf r^\top V^{-1}\mathbf r$: $3\times10^{-15}$.

**What the form makes visible.** Every reading contributes through one
number, $H_i(m)$, and the size measure acts only as the function
$H\mapsto\log E(H)$ that turns explained energy into evidence. ours27's
function is linear; ours28's grows like $H/2-\tfrac{k+1}{2}\log H$. The
readout, the stacked length, T7 and composition are untouched.

## 2. The declared candidate subspace, and a correction to my known-empty record

**The type.** A reading declares $(\mathbf y_i,\Sigma_i,N_i)$. $N_i=I_d$
($k=d$) declares nothing about the push's direction. $k<d$ declares that the
push lies in the span of $N_i$: a cart on a track is $k=1$ along the track.
"Known-empty" in the candidate-size record §5 is $k<d$ with $N_i$ a subset
of the axes; this is its general form, with any direction.

**The solid-angle fact.** Take one reading with the push along axis 1 and
two transverse axes that read zero within their noise, and compare $k=3$
(undeclared) with $k=1$ (declared). In their own noise units the transverse
readings are $N(0,1)$ whatever their $\sigma$, so $H$ has the same law and
$E(H)$ differs only by $k$: ${}_1F_1(\tfrac12;\tfrac52;\cdot)$ against
${}_1F_1(\tfrac12;\tfrac32;\cdot)$. The widening (known-empty record §2: at
$r=4$, 3.65 against 3.04 in 95% log-width) is therefore **independent of the
transverse instruments' precision**. In physical units it is the same fact:
under any isotropic push measure, "the push pointed exactly along one axis"
is a small-solid-angle event whose likelihood ratio against "the push was
small" is $\propto1/\text{scale}^2$, with the transverse $\sigma$ entering
only as a data constant. Codex's zero-reading record §4 states the mechanics
("$\delta=\epsilon$, $\epsilon\to0$: the factor remains $t^{-1/2}$").

**So the verdict of `ours28_known_empty_vs_measured_empty_2026-09-28.md` §3
was wrong**, and I withdraw it: the widening is not "physically backwards"
and it is not a units fork. It is the correct cost of not declaring a
direction. If the push was aligned with an axis on purpose, that is a
declaration, and the type carries it ($k=1$). The candidate repair in that
record, one $\omega$ per noise pair, is withdrawn with it: per-pair scales are
frame-dependent when the noise is isotropic (three independent scales in an
arbitrary frame are not rotation-invariant), so they would jump at isotropy.
The record's §2 (the declaration buys width under ours28 and bought nothing
under ours27) stands.

**What is left of the units fork.** T7 forces the size measure to be
isotropic in whitened (noise) units; the push-units bend is forced out three
ways. The scenario that seemed to argue for physical units no longer does.
One scenario is still open: pushes that are physically isotropic under an
anisotropic instrument, where "isotropic in whitened units" is a different
direction law from the physical one. Its test is the own-world test of §5
run in a physically isotropic world; not done.

## 3. Codex's heterogeneous counterexample, resolved by the declaration

`ours27_limits_of_scalar_and_tensor_size_repairs_2026-09-28.md` §1: $d=2$,
$\sigma_F=(1,10)$, $\sigma_a=(1,1)$, $m_0=1$, one fixed push with
standardized lengths $r^2=(0.1,9.9)$, repeated $N$ times; the signal sits
mostly on the axis least sensitive to mass, so one scalar shrinkage cannot
match the sampling variance (predicted SD ratio 1.615).

Run here with $\pi$ the half-Cauchy centred at $\sqrt{s_1s_2}$, 200 series
per row, the push's direction $\mathbf n=\mathbf a^\ast/\lVert\mathbf a^\ast\rVert$
either undeclared ($k=2$) or declared ($k=1$, $N_i=\mathbf n$). Coverage of
the 95% interval / ratio of the sampling SD of the median to the law's SD:

| $N$ | ours27, undeclared | ours28, undeclared | ours27, declared | ours28, declared |
|---|---|---|---|---|
| 84 | 74.0 / 1.65 | 85.5 / 1.34 | 93.5 / 0.97 | 96.0 / 0.88 |
| 336 | 73.0 / 2.18 | 78.5 / 1.83 | 94.0 / 1.02 | 94.5 / 0.93 |
| 1344 | 78.5 / 1.66 | 84.0 / 1.48 | 93.5 / 1.01 | 95.0 / 0.95 |

- Codex's number is reproduced (1.66 against 1.615) for ours27 undeclared.
- **ours28 undeclared does not fix it** (1.48). A per-reading scalar
  tempers the whole reading's misfit uniformly, and with rank-one data per
  reading a per-reading tensor is the same scalar; the direction structure
  cannot be learned within a reading, and learning it across readings
  breaks Z and meets the improper-measure problem Codex found.
- **Declared, both calibrate.** ours27 declared sits at 1.01, which is T18
  with $M=Nk=N$ and $R^{\ast2}=10N$ ($\sqrt{1.1}=1.05$); ours28 declared at
  0.95, slightly wide, coverage 95%.

**Reading.** The defect is not in the size measure; it is a fixed direction
left undeclared, in a world where the direction never varies. The estimator
was told the push could point anywhere and the world never let it. That is
the same relation as a wrong $\sigma$, and the morning ruling ("the
$\sigma$'s are correct") extends to it: **the declared inputs are correct**,
and the push's direction structure is one of them. When the direction is
unknown, the isotropic hedge is what "unknown" means, and §5 says what
calibration then means.

The ladder of declarations, all inputs to E10, none learned: $k=d$ (nothing
known about direction) → $N_i$ with $k<d$ (direction known) → a declared push
covariance $R_i$ (direction and size known approximately; Codex's tensor, as
an input rather than a learned state) → old readings in the prior slot
(`old_information_in_the_prior_slot_2026-09-24.md`).

## 4. The fork

`ours27_c3_from_betweenness_and_the_scale_fork_2026-09-28.md` §3.2: grouping
as bookkeeping, the zero-reading property Z, and a proper size measure cannot
all hold, and the ruling owed is whether **which components were measured
together is bookkeeping or physics**.

**Recommendation: physics, because a reading is one push.** A reading's
components are the components of one physical event, one applied force at
one moment; two readings are two pushes. A "reading" made of $x$ of reading 1
and $y$ of reading 2 corresponds to no push. Codex's axis-splitting record
says the same from the other side: representation invariance holds when the
shared state travels with the recoding, and independent elimination is an
extra rule. Under this ruling:

- **L8** reads: a reading and its components give the same law when the
  components carry the reading's shared state (its $N_i$ and its
  $\beta_i$). That is what ours28 does, and what §1's form makes automatic.
- **Z** is kept, per reading, exactly.
- **Properness** is kept, so the own-world calibration of §5 holds.
- **C3 shrinks to its narrow form and stays a premise:** the readout sees a
  candidate through its two joined sizes, not through the profile of
  per-reading sizes (the $\sum\lVert\mathbf v_i\rVert^4/(\sum\lVert\mathbf v_i\rVert^2)^2$
  weight of the betweenness record §1.2). The cross-reading regrouping that
  would kill it is exactly what "a reading is one push" denies. I would fold
  the narrow form into L12, which already says a series has one size per
  channel per candidate: *and the readout is a function of the law of those
  two numbers*. One premise, stated once, whose content is that the series
  is measured as one thing. What the betweenness derivation still gives:
  betweenness, own frames and signs, and same-law-same-mass reduce C3 to
  that one line.

What is given up: the Kingman route to "one Gaussian scale per series" (it
needs regrouping across readings), and the shared-scale form's four-number
state. What is gained: Z, properness, the known-empty behaviour of §2, and
the declaration ladder of §3, which needs the reading to be a unit.

## 5. What calibration means, and where the shapes live

The betweenness record §2 proposes the non-realist yardstick: feed the
definition data generated by its own measure, and its intervals must be
calibrated. For a proper measure that is a theorem (a correctly specified
model's credible intervals cover at their nominal rate on average under the
prior), so every proper ours28 passes it and the improper ours27 fails it.
Checked here in the anisotropic instrument of §3, $k=2$: 95.0% at $N=6$,
94.0% at $N=21$ (300 series; the $N=84$ row read 91% because my mass grid was
coarser than the posterior at large $\omega$; not a finding).

Two consequences.

- **Planted-mass worlds are robustness checks**, and every failure in one
  so far has named an undeclared input: an infinite push (L10, the ours28
  record §1), a fixed direction (§3 here), an aligned push (§2 here). That
  is the pattern to expect: the definition is a procedure on declared
  inputs, and a world that fixes what the declaration left open is a world
  the declaration does not describe.
- **The two shapes are declared.** $g$ for the mass and $a$ for the size
  both pass the own-world test, both pass every property (they are the
  prior-shape bends the matrix found free), and no yardstick internal to the
  project separates them. T20 already says $\pi$ is declared before data as
  part of the specification. The size shape goes in the same slot. The
  current declaration is one rule for both: uniform in the noise-defined
  angle, the half-Cauchy on $m/s$ and on $\omega$. Whether a declared shape
  is "a rival definition of mass" is then the same question as whether a
  second instrument is, and I would say it is not: it is the experiment's
  description, and the definition is what is done with it.

## 6. Where ours28 stands

**The equation.** §1, with the readout of ours25b/26/27 unchanged
(ratio of mean stacked sizes over the candidates, which now include
$\beta_i$).

**Its declared inputs.** Readings $\mathbf y_i$; error covariances
$\Sigma_i$; candidate subspaces $N_i$; the starting law's shape (half-Cauchy)
and centre (T20); the size law's shape (uniform angle). All assumed correct
(the morning ruling, extended).

**Premises it rests on, beyond the lemma record's L1–L9, L11–L14.** L4′
(push–noise symmetry; Jeffreys-check record §4); grouping is physics (§4
here); C3 narrow, folded into L12 (§4 here). L10 is gone. L13 stays forced
given the flat push and is untested for the tempered evidence.

**What passes.** Composition, swap, units, rotation (any $\Sigma_i$), axis
splitting with the state carried, Z exactly, consistency (within the
evening's Monte Carlo), own-world calibration (theorem), planted-mass
calibration in every world whose fixed features are declared.

**Still owed, in order.**

1. The matrix rerun on ours28 in the form of §1, with two new columns: the
   declared subspace ($k$ against embedding) and the size-angle rule bent.
2. L13 under ours28: the Gaussian-noise forcing used the flat push; redo
   with the tempered evidence.
3. The physically-isotropic-push world under an anisotropic instrument
   (§2, last paragraph), as an own-world test.
4. The readout point and its interval under ours28. Nothing tonight touched
   the point.
5. The tester build: the reading type with $N_i$; `ratio_of_means_rss` on
   the §1 law; the five cells of the ours28 record and Codex's example as
   presets.
6. The next lemma record, carrying the sweep's corrections, the matrix's,
   the ours28 record's §8, the Jeffreys-check record's L4′ and this record's
   §4 and §6.

## 7. Reproduction

Scripts in Claude's workspace on 28 September, not kept; numpy and scipy.
`o28full.py`: $H_i(m)$ by the projector formula, the four checks of §1.
`codexD.py`: Codex's example with $N_i=I$ and $N_i=\mathbf n$, 1601-point
grid in $\log m$ on $[-2.5,2.5]$, ${}_1F_1$ from `scipy.special.hyp1f1`,
half-Cauchy at $\sqrt{10}$. `ownworld.py`: masses drawn half-Cauchy, per
reading $\omega_i=\tan\varphi_i$ with $\varphi_i$ uniform, pushes
$N(0,\omega_i^2I_k)$ in the whitened subspace, 3601-point grid on $[-9,9]$,
the evidence by quadrature in $\varphi$ (240 points) in log space.

## 8. How we got here

1. The user: "new documents just dropped... continue to work on the
   solution... worry about solving it yourself not correcting others."
2. Read the five records. The betweenness record's fork and Codex's
   heterogeneous counterexample were taken as the two things that bear on
   ours28; Codex's zero-reading record §4 was taken as correcting my
   known-empty §3.
3. ours28 was rewritten in the projector form to have a single object to
   test; its checks passed first time.
4. Codex's example was run undeclared and declared. A first run used a flat
   log-mass law and small $N$, where the readings' weak sensitivity
   ($J=0.12$ per reading) made everything biased and wide; rerun with the
   half-Cauchy at larger $N$.
5. The own-world check overflowed ${}_1F_1$ and truncated the half-Cauchy
   tails; redone in log space on a wider grid.
