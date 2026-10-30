# The proof: from the type to ours30, and what is left

**What this is.** The second half of the deliverable. `equation.md` states ours30; this file is
the argument that ours30 is *the* function, written as one chain from the type and the premises
to the equation, with every step marked. Started 2026-10-29 by Claude on Can's "Ok go ahead".
It is a draft with holes, on purpose: the holes (§5) are the queue. Edit in place; git holds
the history.

**Statuses.** Each step carries one:

- **ruled**: Can's decision; the ruling ID is in `../rulings.md`.
- **forced**: follows from named premises by an argument that is *written* in this repo, with
  the file named. "Carried from the records" does not count as written here.
- **forced (argument in archive)**: the argument exists in `../docs/` but has not been written
  into the theorem file. A writing task, not a research task.
- **input**: the readings cannot decide it; it is declared with them (`../TASK.md` §2).
- **provisional**: Claude's sort, from the goal, open to Can's check (R-2026-10-27-1).
- **gap**: no argument, or an argument that does not reach the claim. Typed in §5.

**The rule for this file.** A step may cite a test as its check, but a test is not a proof.
Where the proof is missing the step says *gap*, however well the test passes.

**What is to be proved.** Let $\mathcal F$ be the functions from the type (`type.md`: readings
with $\mathbf y_i,\Sigma_i,N_i$, structure, window; a declared $\pi$) to an output (a law on
$m>0$ and a point). Claim: the premises in `premises.md`, with the rulings, select one member
of $\mathcal F$, and it is ours30. Two halves:

- **(A) ours30 satisfies the premises.** Mostly done, by tests (§3).
- **(B) nothing else does.** The uniqueness half. §2 is the chain that would prove it; §4 says
  how much of it is proved today.

Agreement with experiment (`../TASK.md` §5) is a check on the result, not part of this proof;
it has not been made (O16, brief B01).

---

## 1. Before the chain: the two branches of the type

**Step 0. Exact inputs ($\sigma=0$).** From L1 alone: an exact pair $(\mathbf F,\mathbf a)$ has one
mass iff the two pairings $m_1=\mathbf F\cdot\mathbf a/\lvert\mathbf a\rvert^2$ and
$m_2=\lvert\mathbf F\rvert^2/(\mathbf F\cdot\mathbf a)$ agree on a value in $(0,\infty)$; no mass if they
disagree ($\mathbf F\nparallel\mathbf a$) or agree outside $(0,\infty)$; every mass at exact $0/0$ (T1).
*forced* from L1; the algebra is in T1's statement; `tests/test_full.py::test_exact_branch`. The
one place "undefined" survives.

**Step 0′. Measured inputs return a law and a point, never "undefined".** The type: every
$\sigma>0$ (L2, *ruled* R-2026-09-25-1); on measured inputs mass is never undefined, the
all-zero experiment included (*ruled* R-2026-09-29-2). Why a *law*: a measured pair never fixes
$m$, so the honest output is what the readings say about every $m$; the point is a summary of it
(F1, *ruled* R-2026-09-27-1: nothing behind the readout). *ruled.*

Everything below is about the measured branch.

---

## 2. The chain

Each step: the claim; what it rests on; where the argument is written; status.

### 2.1 What a reading can mean

**Step 1. The candidates.** Under L1 with a scalar mass (L22), the true pair behind reading $i$
is $(\mathbf F^\ast_i,\mathbf a^\ast_i)=(m\,N_i\alpha_i,\ N_i\alpha_i)$ for some $\alpha_i\in\mathbb R^{k_i}$:
one direction shared by both channels, inside the declared subspace, with $m>0$ (T2).
Rests on L1 (bedrock), L22 (that $m$ is a scalar), L18 (that $N_i$ is a correct declaration).
*forced* from L1 and L22; trivial once L22 is granted. **L22 is provisional** (Claude, 27 Oct:
a clause of L1, "Newton's mass is a scalar"); **$N_i$ as an input is proposed, not ruled** (L18's
second clause). → gaps G10, G11.

**Step 2. The noise model is Gaussian with the declared $\Sigma_i$.** The model attaches to the
data a density $\phi_{\Sigma_i}(\mathbf y_i-\text{true pair})$. L13 says the model uses only the declared
second moments and is not a claim about the world's errors (R-2026-09-27-12). That the family is
*Gaussian* and not another second-moment family: T23 gives two routes. The Z route: only the
Gaussian leaves an exact-zero reading's factor constant in $m$ (L20); written as a test,
`tests/test_theorems.py::test_t23_only_gaussian_noise_leaves_a_zero_reading_alone`, under the
flat weight. The consistency route (L11) under the tempered evidence is O4, unwritten.
*forced (test only)* from L20 and L11. The derivation behind the Z route is not written out;
the consistency route is open. → gap G1.

**The frame for steps 1–5 (T36, 29 Oct).** The reading space $W_i$, with the noise metric
$\Sigma_i^{-1}$ forced by units; the hypothesis $m$ as the subspace $L_i(m)$ (step 1's candidates);
$H_i$ as projection onto it and $K_i=E_i-H_i$ as the residual; $k_i$ residual dimensions per reading,
$\nu_i=k_i-1$ when the reading picks its own line. The raw law reads one residual,
$\sum_iK_i(m)$. T36 is geometry and a dimension count; it removes no hole, and it states G13's and
G3's questions (what the power on $E_i$ and η are for) in the form in which they are now asked (O14).

**Step 3. The candidate weight, and the per-reading factor.** Integrating the candidates out
needs a weight $w(\alpha_i)$ on the true acceleration. Four premises fix it:

- L8 (conditional form, ruled R-2026-09-29-9): components independent given the reading's
  shared state, so $w$ is a product over components;
- L9 (strong): $w$ depends on $\alpha_i$ only through its whitened length;
- a density on $\mathbb R^{k}$, $k\ge2$, that is both a product and rotation-invariant is
  Gaussian, $e^{-\lVert\alpha\rVert^2/2\omega^2}$ (Herschel–Maxwell);
- L10 (ruled R-2026-09-30-6): no size is declared, so no $\omega$: the weight is the
  $\omega\to\infty$ endpoint, flat in the vector;
- L7 with L13 (the width at strong signal is the Gaussian propagation of $\Sigma_i$) removes the
  improper members $e^{+c\lVert\alpha\rVert^2}$ and fixes the coefficient $\tfrac12$;
- T7 (from L11): nothing attached to a reading may depend on $m$, which puts the weight in noise
  units and drops its $m$-free constant.

Then $E_i(m)=\int\phi_{\Sigma_i}\big(\mathbf y_i-C_m\alpha\big)\,d\alpha\ \propto\ e^{H_i(m)/2}$, with
$H_i$ the explained energy of `equation.md` (T6, T8's per-reading factor). *forced (argument in
archive)*: the argument is written in full in
`../docs/vd-docs-context/o3_what_supports_the_flat_weighting_and_the_stacked_length_2026-09-30.md`
§2, not in T06.md. **Two things it does not cover.** (i) Herschel–Maxwell needs $k\ge2$. The 1D
factor is carried by "a reading declared $k=1$ is the 1D reading" (T25). T25 is a theorem about
the kernel in use; as a *forcing* step it needs the premise that declaring $k=1$ of a
$d$-dimensional reading is the same reading as a 1D one, which L8's "same data, different
declaration" says but no row states as a premise. (ii) It is stated for one reading's $\Sigma_i$
with the conditional form of L8 "given $N_i$, $Q_i$": whether the weight may depend on $Q_i$ is what
T30 settles later, so the two steps lean on each other. → gaps G2, G17.

**Step 4. Noise units, and nothing else.** The whitening in step 3 is by $\Sigma_i$ and the
candidate map $C_m$ is scaled so every slope is treated alike (`equation.md` §Per reading). That
the law is covariant under a change of units (L3) and under a rotation of the reading's frame
with its covariance (L16, T22) are checks: `tests/test_full.py::test_rotation_with_covariance_changes_nothing`,
`::test_swap_gives_the_reciprocal`. *forced* from L3, L16; the covariance of $H_i$ under both is
two lines of algebra, not yet written in T22. → part of G15.

### 2.2 Joining readings

**Step 5. The product.** The joint factor of a series is $\prod_iE_i(m)$ (T12). Rests on L13
(errors independent across readings: the joint density factorises) and L19 (no borrowing: no
latent quantity shared across readings besides $m$, so no hierarchical or pooled factor). *forced*
from L13 and L19. The argument is one sentence and is not written anywhere but here. **A premise
is doing silent work**: L13's "independent across readings" is a declaration about the instrument
(errors shared across readings have no slot, O16). As a forcing step this is fine: the function
is defined on what is declared. → noted, not a gap.

**Step 6. The starting law $\pi$, declared once.** $\pi$ enters once, not per reading (T7, from
L11: a per-reading prior would make the law depend on $N$ at fixed total signal). It is declared
before data (T20). Its family is forced: symmetric under $m/c\leftrightarrow c/m$ (L4), covariant under
units (L3), proper (T27's reason: an improper $\pi$ gives no law on null input), so
$\pi(m)=m^{-1}g(\ln(m/c))$ with $g$ even and integrable (O6). **Its shape within the family is an
input** (R-2026-10-24-3: the readings cannot decide it; the half-Cauchy, uniform in the angle the
noise defines, is the principal branch, T28). **Its centre is an input**; $c=s$ for one noise pair
is the principal branch (L3′, O11, the units argument of 10 Oct), *not ruled*. *input*, with the
family *forced (argument in archive)*: the family derivation is spread over O6, T28 and the
Jeffreys-check record, not written in one place. → gaps G8, G9.

### 2.3 How much each reading counts

**Step 7. The weight is a power on the reading's factor.** L21 (*ruled* R-2026-09-29-8): the law
is $\pi\prod_iE_i^{g(Q_i)}$, $g(0)=1$, $g$ decreasing, $g>0$. The *form* (a power, one per reading, a
function of the reading's own misfit) is the ruling's statement. **Why a power and not another
tempering** is not argued anywhere: L21 was ruled in that form. → gap G13.
**30 Oct (T36, T37):** a power constant in $m$ is a declared $\Sigma_i/g$ in disguise (homogeneity).
With the σ's believed, the cut's exact per-reading factor is Student-type (T37 (4)), and the ruled
power form is its first order in the along-line misfit: equal near the reading's own line, not far
from it. Which is the body is a ruling. → G19.

**Step 8. What the misfit is.** The smallest error that reconciles the reading with the law as
assumed, $m>0$ included (R-2026-09-07-1), so $Q_{+,i}$ (T32: $\lambda_-$ of the whitened Gram
matrix when $\mathbf x\cdot\mathbf y>0$, else $\min(A,B)$; for a general $\Sigma_i$, by search over
$\phi\in[0,\pi/2]$). *ruled* R-2026-10-29-1; the closed form is T32's, written there with its
proof (Codex, 26 Oct) and tested (`tests/test_sign_numerator.py`). *forced.* Scope note: T32's
closed form is for isotropic uncorrelated channels; the general case is defined by the search and
tested against a dense grid, with no closed form claimed.

**Step 9. The two ends and the middle of $g$.** $g(0)=1$ from Z (an exact-zero reading has
$Q_+=0$ and must not change the law). The tail $g\sim\kappa/Q_+$ from the requirement of bounded,
nonzero confidence on a fixed pair no positive mass fits as its noise shrinks (T30). The middle,
$g=\min(1,\kappa/Q_+)$, from the *meaning* of $\kappa/Q_+$ as the share of the reading's
self-contradiction the declared noise accounts for, read at every $Q_+$ and capped at the whole
(*ruled* R-2026-09-30-10, O12). The ends: *forced (argument in archive)* for $g(0)=1$ (Z, one line);
**the bounded-confidence requirement itself is not derived from a named premise.** T30 states it
as a requirement; L21 gives "not none" ($g>0$) and "less" (decreasing), not "bounded". The
candidate source is R-2026-09-07-1 ("it must simply be an unlikely thing to happen": a reading no
mass fits is an improbable event under the law, and an improbable event's evidence is finite),
but that is a reading, not an argument. → gap G13.
**30 Oct: the bounded requirement derived, and the middle in conflict.** With each reading's noise
size its shared state (L8) and Σ its covariance (T37, Claude's formalisation of R-2026-10-30-1), the
bounded, non-vanishing requirement forces a power-law tail on the noise size, tail weight
$(\nu+2\alpha)/Q$ (T37 (1)): this is the argument G13 lacked. But the same model puts the clean end
above 1 (T37 (3)), against the ruled cap. And the cap conflicts with R1 without T37's model: by
homogeneity $g<1$ is the law reading the misfit as evidence the σ's were too small, which R1 says
the misfit is not. See §4, "The rulings at steps 9–10 conflict". → G18.

**Step 10. The numerator $\kappa$.** O12's meaning, "the declared noise's expected misfit", read
without a declared size (L10): the most the noise produces over every allowed signal size.
That reading is *ruled* (R-2026-10-29-1). Its value, $k-\tfrac12$ for isotropic uncorrelated
channels, is **numerical** (a seeded grid: `tests/test_sign_numerator.py::test_the_numerator_is_the_supremum_over_sizes`;
Codex's bounds $k-\tfrac12\le S_k\le k$ in O13). Under a general covariance the value is **unknown**.
→ gaps G4 (the proof, brief B02), G5 (general $\Sigma$).
**30 Oct:** with Σ the covariance, the tail coefficient is $\nu+2\alpha>k+1$; $k-\tfrac12$ needs a noise size
with no mean (T37 (2)). So the meaning ruled for $\kappa$, the most misfit the declared noise produces,
is the floor's meaning (§4). → G18; G4 and G5 stand only if G18 keeps the ruled meaning.

### 2.4 How much the series counts

**Step 11. One power on the whole series, $\eta$.** The flat factor overstates the series'
precision when the signals are weak (T18: the width error is a function of the true signals, which
L10 says are undeclared). The correction must come from the readings' disagreement about $m$ and
nothing else (L19's scope, *ruled* R-2026-09-30-5: the Birge move is allowed) and must leave an
exact-zero reading alone (Z). Property N (*ruled* R-2026-09-29-4): on a null series the quoted
width must be a true statement, which fixes the ratio $\mathcal H/J$ in the long-series limit
(Bartlett's identity: the curvature of the summed log-evidence equals the variance of its score
when the model is right). One $\eta$ per series because $m$ is one number (L1). The raw $J$
(no Bessel factor) because a factor that counts readings moves on appended zeros (Z). The
evaluation point: the series' best line on the whole projective line, because the ratio is
parametrisation-free only where the score vanishes, and on $m>0$ the best line can sit at a
boundary (T31; checked against R-2026-09-07-1 on 29 Oct). Candidates unchanged by $\eta$ (L19).
*forced in the limit* from property N, L19, L20, L1; T31's statement carries the argument in
compressed form; nothing is written out as a proof. **Not forced:** the cap at 1, and the power
form at finite $N$ (T31 says "natural"). **Under F2 a choice inside the function is an alternative
mechanics**, so this is the one component of the law with an unsorted freedom in it. → gap G3
(a sort; brief B04 supplies the finite-$N$ facts it needs).
**30 Oct:** η is a power on every reading at once, so under L24 (candidate) it is admissible only as
the cut posterior of a mean-1 noise size shared by the series, which is a correlation between
readings the type does not declare. Brief B12 item 4 checks it. → G21.

**Step 12. The law.** $p(m)\propto\pi(m)\exp(\eta\sum_ig_iH_i(m)/2)$ on $m>0$ (T8). Steps 3–11
assemble it. The closed form for one noise pair (`equation.md`, note d3 §3) is checked
against the numerical path (`tests/test_ours30.py::test_the_closed_forms_in_equation_note_d3`).

### 2.5 The point

**Step 13. The point is a functional of the law.** L19: readings are joined only through $m$;
any point that joins them through their signal sizes borrows (the cloud readout did, through the
stacked length: `tests/test_readout_weight.py`). So the point is a function of the law of $m$ and the
declared inputs. *forced* from L19; the argument is written in O2 (28 Oct, §1). **Not from Z:**
Z is ruled for the law only (R-2026-10-29-2); that the median does not move on null readings is a
consequence of this step, not a requirement it was chosen to meet.

**Step 14. No point meets composition at finite noise.** L5 and L14 as constraints on the point
(masses add; inverse masses add) cannot both hold for any functional of spread laws (Jensen; the
law has no mean, T19) and on data no point meets either (`tests/test_point_composition.py`,
`tests/test_composition_law.py`). So they are commitments about the *law*, which exact inputs
(T1) and consistency (T17) deliver for any consistent point (*ruled, provisional*
R-2026-10-28-1). The law-level statement is checked by coverage. *provisional*; the no-go is
written in O2 (28 Oct, §§2–3).

**Step 15. The swap.** Swapping force and acceleration (with their covariances) and reading
$m$ as $1/m$ must give the reciprocal answer: force and acceleration are equal information
(*ruled* R-2026-09-24-1, R-2026-09-07-1). The *law* has this symmetry: $H_i$ is invariant under
the swap with $\phi\mapsto\pi/2-\phi$, $\pi$ by L4, $Q_+$ and $\kappa$ unchanged, $\eta$ unchanged
(`tests/test_sign_numerator.py::test_swap_and_units`, `tests/test_ours29.py::test_rotation_swap_and_units`).
So the point must commute with $m\mapsto c^2/m$. **This replaces T5's old derivation** (T5 was a
theorem of T4, the ratio of means, which is no longer the point). *forced* from R-2026-09-24-1 and
the law's symmetry; the two-line argument is not written in T05.md. → part of G15 (G7).

**Step 16. The median.** Among functionals of the law that commute with the swap, respect units
(L3) and return $c$ on the all-zero experiment (step 17), a family remains: every location
functional of $\ln m$ that commutes with negation (median, geometric mean, …). L23 (*proposed*,
28 Oct: the point commutes with every increasing reparametrisation of mass, as T31 already
requires of $\eta$) makes the point a quantile, and the swap makes it the median (T35, written
with its proof and test). *provisional on L23.* Without L23 the choice in that family is
"unknown yet" (R-2026-10-27-1), and the median is the default by T35's argument applied to the
sub-family that is parametrisation-free. → gap G6 (ruling).

**Step 17. The all-zero experiment returns $c$.** For the median this is immediate: the all-zero
law is exactly $\pi$ (T24: Z for the law) and $\pi$'s median is $c$ (L4). The older argument through
Z at the point (the answer depends neither on $N$ nor on $\Sigma_i$, so by L3 it is $kc$ and by the
swap $k=1$; undefined-mass record §4) assumed point-level Z, which is not ruled (R-2026-10-29-2),
and is not needed. *forced* from L20 (law), L4, step 13; `tests/test_point_composition.py::test_point_returns_the_centre_on_the_all_zero_experiment`.

**Step 18. The interval.** The quoted interval is the law's central 95% with $\eta$ applied; the
raw law's beside it (*ruled* R-2026-09-30-7). The level 95% is a reporting convention, not part
of the definition: the law is the output and the interval summarises it. *ruled.*

---

## 3. Half (A): ours30 satisfies the premises

| premise | holds for ours30 | checked by |
| --- | --- | --- |
| L1 (one $m>0$, shared) | the law lives on $m>0$; one $m$ per series | by construction |
| L2 (never undefined) | T13 | `tests/test_full.py::test_every_measured_input_gives_a_mass`; the all-zero case in `tests/test_point_composition.py` |
| L3 (units) | law, point, interval scale by $\lambda/\mu$ | `tests/test_sign_numerator.py::test_swap_and_units`, `tests/test_ours29.py::test_rotation_swap_and_units` |
| L4 (swap, law) | $\pi$ symmetric; the whole law symmetric (step 15) | same |
| L5, L14 (composition, law level) | computed and measured laws cover alike | `tests/test_composition_law.py` |
| L8, L9, L10 | the per-reading factor is the flat Gaussian endpoint | `tests/test_theorems.py::test_t6_t7_*`; `tests/test_cartesian.py` |
| L11 (consistency) | T17 | `tests/test_theorems.py::test_t17_*` |
| L13 (Gaussian model) | by construction | T23's test |
| L16 (frame) | T22 | `tests/test_full.py::test_rotation_with_covariance_changes_nothing` |
| L17 (same candidate law, same mass) | the point is a law functional (step 13) | by construction since 28 Oct |
| L18 (declared inputs) | used as declared | by construction |
| L19 (no borrowing) | nothing crosses readings but $m$; $\eta$ reads disagreement only | `tests/test_readout_weight.py` (the old readout's failure), step 13 |
| L20 (Z, for the law) | T24. (The median also does not move on null readings, as a consequence of step 13; point-level Z is not ruled, R-2026-10-29-2, O7) | `tests/test_full.py::test_a_zero_reading_leaves_the_law_alone_and_moves_the_readout`; `tests/test_point_composition.py::test_point_is_the_median_and_meets_z_exactly` |
| L21 (evidence, just less) | $0<g\le1$, decreasing in $Q_+$ | `tests/test_sign_numerator.py::test_no_gate_*` |
| L22 (scalar mass) | the candidates share one direction | by construction (T2) |
| L23 (parametrisation-free point) | the median | `tests/test_point_composition.py::test_t35_*` |
| property N | $\eta$ removes the $\sqrt N$ null tilt | `tests/test_ours29.py`; O15 open on whether the floor is the best possible |
| R-2026-09-07-1 (smallest error under the law) | $Q_+$ | `tests/test_sign_numerator.py::test_q_plus_*` |
| T30's bounded confidence | anti-parallel pairs bounded in 1D–3D | `tests/test_sign_numerator.py::test_t30_bounded_tilt_*` |
| **R1 (the σ's correct; misfit is not evidence about them)** | **fails** (30 Oct): by homogeneity $g<1$ declares $\Sigma_i/g$, i.e. reads the misfit as evidence the σ's were too small | `tests/test_scale_reading.py::test_tempering_a_reading_is_declaring_its_noise_larger` |
| **L24 (candidate)** | **fails** for $g$ (T37 (2), (3)); η pending | `tests/test_noise_size.py`; brief B12 |

Not checked: ours30's **coverage** on its own world (only ours29's was measured; ours30 changes
the weights on acute pairs). A study, not a proof step; queued with B04's cells.

---

## 4. Half (B): what the chain proves today

Reading §2 with the statuses:

- **Forced, with the argument written in this repo:** steps 0, 8, 13, 15 (two lines owed), 17,
  T35's conditional (step 16).
- **Forced, argument in the archive, not yet in the theorem file:** steps 3 (T6), 6 (the family
  of $\pi$), 9's ends, 11 (T31, compressed), and most of the consistency theorems (T12, T17,
  T18, T19, T22, T24, T25: 31 of 35 theorem files say "carried from the records").
- **Ruled:** steps 0′, 7 (the form of the weight), 8, 9's middle, 10's meaning, 14 (provisional),
  18, the inputs of step 6.
- **Provisional:** L22 (step 1), step 14, step 16 (L23).
- **Gaps that are research, not writing:** step 2's consistency route; step 9's bounded
  requirement; step 10's value and its general-covariance form; **step 11's finite-$N$ form**.

So the honest statement of what is proved: *on isotropic uncorrelated channels, any function
that satisfies the premises as ruled, with L22, L23 and the finite-$N$ form of $\eta$ taken as
given, is ours30 up to the declared inputs $\Sigma_i,N_i$, the structure, the window, and $\pi$'s
shape and centre.* The three things "taken as given" are the uniqueness half's open work, and
two of them are rulings rather than mathematics.

### The rulings at steps 9–10 conflict (30 Oct; Claude's sort, provisional)

Two groups of rulings meet at the misfit weight:

- **The σ's are believed.** R1 (28 Sept): "the declared σ's are assumed correct… not whether misfit
  is evidence about the σ's". R-2026-09-29-8 ("evidence, just less"), chosen over the option "misfit as
  evidence the σ's are small" (readings-are-evidence record §1.6). R-2026-10-30-1: a misaligned
  reading is possible under the same σ's, only unlikely.
- **The weight's middle and numerator.** R-2026-09-30-10 (O12): $g=\min(1,\kappa/Q)$, the share of
  the reading's self-contradiction the declared noise accounts for, capped at the whole.
  R-2026-10-29-1: $\kappa$ the most misfit the declared noise produces over every size.

They cannot all hold. By T36's homogeneity any $g<1$ is a declared $\Sigma_i/g$: the second group
makes the misfit evidence that the σ's were too small, which the first group rules out. With the
σ's formalised as the covariance (T37), the second group's cap and $\kappa$ are also unreachable.
By `TASK.md` §4 one premise is mis-stated. **The sort: the mis-stated one is O12's meaning.** "A share
of the self-contradiction the declared noise accounts for, capped at the whole" presupposes that the
declared noise is the most there can be: the floor, the option declined on 29 Sept. The first group
is three rulings over two months saying the same thing; the second group's meaning is one reading of
it that turned out to encode its opposite. Can rules whether the sort holds. → G18.

**What is proved, restated.** Until G18 is ruled: *on isotropic uncorrelated channels, any function
that satisfies the premises as ruled, with R1 read as not constraining the misfit weight, and with
L22, L23 and η's finite-$N$ form taken as given, is ours30 up to the declared inputs.* If R1
constrains the weight (the sort above, and L24), **steps 9–11 are open and ours30 is not the
function**; T37 gives the form the reading level would take (a cut posterior of a mean-1 noise size
with a power-law tail), with its tail index an input whose principal branch is unknown (G20).

**What is no longer load-bearing.** Switching the point to the median (28 Oct) took these off
the path from premises to ours30: T3, T4 (ratio of means), T9 (stacked size), T10 (the cloud
readout), T11 (sufficiency through the stacked sums), T14, T15 (the readout at zero readings and
equal energies), T21 (C3) and its owed step O5, T29 (the readout fails Z), and L12 (one joining
rule) with the series-size question. They remain true statements about the cloud readout, kept
for comparison, and they come back only if the ratio of means does. The audit should stop
counting L12 and O5 as owed by the definition. → G12.

**Two rungs are not fixed in VD itself** (`ladder.md`): R3 (NML2.2: persistence, lumping) is
unwritten and R5 (Newton III) is provisional. L5, L14 and the use of Newton III (F4) inherit
that. The proof here is relative to the version as it stands; it cannot be stronger than the
rungs. → G11.

---

## 5. The holes, typed

Each is one of: **write** (the argument exists, put it in the file: Claude), **derive** (a bounded
derivation with a named deliverable: a brief for Codex), **sort** (a fork to sort against the
goal: Claude), **rule** (Can), **VD-Newton** (a rung to write there: Can). Order is the chain's.

| id | step | hole | type | where |
| --- | --- | --- | --- | --- |
| G10 | 1 | L22, scalar mass, as a clause of L1 | rule | O17 |
| G11 | 1, 4 | $N_i$ as an input (L18's second clause); R3 unwritten, R5 provisional | rule; VD-Newton | L18; `ladder.md` |
| G1 | 2 | T23: the Z route's derivation written out; the consistency route under tempered evidence | write; derive | O4 → brief B05 |
| G2 | 3 | T6's argument moved from the O3 record into T06.md | write | T06 → brief B07 |
| G17 | 3 | the 1D factor: "a reading declared $k=1$ is the 1D reading" needs stating as a premise or deriving from L8 and L18 | sort | T25, L8 |
| G15 | 4, 5, 15, and the consistency theorems | 31 theorem files with the proof "carried from the records": write them out in the chain's order, shortest first (T22, T24, T25, T12, T19, T17, T16, T18, T13) | write; the pure algebra ones can be briefs | theorems/ → briefs B06 (T16–T19, T24, T25), B10 (T5, T12, T13, T22) |
| G8 | 6 | the family $m^{-1}g(\ln m/c)$ derived in one place (L3, L4, properness) | write | O6, T28 → brief B08 |
| G9 | 6 | the centre: L3′'s first half, then $c=s$ | rule | O11 |
| G13 | 7, 9 | why the weight is a power on $E_i$, and which premise makes bounded confidence a requirement | sort | T30, L21. **30 Oct, provisional:** the power is a noise scale read once (L21's cut, O14); bounded confidence derived from a noise size with a power-law tail (T37 (1)), which supersedes the morning's "from Σ as a floor" |
| G18 | 9, 10 | R1 and its kin (the σ's believed) against O12's cap and κ's meaning (R-2026-09-30-10, R-2026-10-29-1) | sort (done, provisional: O12's meaning is the mis-stated one) → rule | §4; O14; T37 |
| G19 | 7 | L21's ruled power form, or the cut's exact factor (T37 (4)), as the body | rule | T37 |
| G20 | 9 | T37's tail index α: an input, but TASK §5's guard asks that practice supply it. **Sorted 30 Oct, provisional (T37, G20 section):** α = ν_eff/2, the effective degrees of freedom of the stated uncertainty (GUM Annex G), on calibration certificates; the inverse-gamma shape forced by L13 applied to how σ was obtained. Practice normalises the precision ($E[1/v]=1$), which restates T37 (2) and (3): ours30 is a σ declared with ν_eff = ½; the cap holds exactly in 1D | sort (done) | T37 |
| G21 | 9, 11 | L24's test: which temperings are cut posteriors of a mean-1 noise size; η's series version | derive | L24 → brief B12 |
| G4 | 10 | $S_k=k-\tfrac12$, and $\sup E[Q]=k-1$, proved | derive | brief B02 |
| G5 | 10 | $\kappa$ under a general $\Sigma_i$ (correlated channels) | derive, after B02 | O13 |
| G3 | 11 | $\eta$'s cap at 1 and its power form at finite $N$: forced, input, or a mis-stated premise (property N stated in the limit only?) | sort | O14; brief B04 supplies the finite-$N$ facts. **30 Oct, partly (O14):** the cap is the series floor and the power form is a scale; the finite-$N$ form is still open. **30 Oct (O15):** property N is stated in the limit, and ours30 meets it there (0.93–0.95 at $N=80$–200, signal 1–2, 3D); only the finite-$N$ form is left |
| G7 | 15 | T5 for the median: the law's swap symmetry, two lines | write | T05 → brief B10 |
| G6 | 16 | L23 | rule | O2 |
| G12 | 4 | the audit and the problem files still count L12, O5, T9–T10's freedoms as the definition's | write | `tools/audit.py`, O2, `ladder.md` → brief B09 |
| G14 | — | agreement with experiment (`../TASK.md` §5) | derive (a study) | brief B01 |

Not holes, but noted: step 5's reliance on L13's "independent across readings" (a declaration;
O16 is where it may gain a slot); the 95% level (a convention, step 18).

---

## 6. How this file is used

- When a hole closes, its row moves to the step it belongs to as a status change, and the
  row is struck here (kept, so references resolve).
- A new freedom found anywhere in the tester is a new row here before it is anything else.
- Briefs for the *derive* rows, and for the *write* rows whose argument already exists, are
  written from this table (`../briefs/`), so Codex works on what the proof needs. The *sort*
  rows are Claude's: a sort is what a brief must have before it is written, so it cannot be
  one. The *rule* rows are Can's, listed in `../rulings.md`'s owed table.
- `docs/STATE.md`'s "Read first" points here after `ladder.md`.

## History

- 2026-10-29: created (Claude). Eighteen steps, sixteen holes.
- 2026-10-29 (later): point-level Z removed as a premise (R-2026-10-29-2); steps 13, 16, 17 and §3 reworded.
- 2026-10-29 (night): the frame for steps 1–5 named (T36), on R-2026-10-29-3.
- 2026-10-30: steps 7, 9, 10, 11 and §§3–4 brought up to T37 and L24: the rulings at steps 9–10 conflict (R1 against O12), sorted provisionally; what is proved restated; holes G18–G21.
