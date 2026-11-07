# The proof: from the type to ours30, and what is left

**What this is.** The second half of the deliverable. `equation.md` states ours30; this file is
the argument that ours30 is *the* function, written as one chain from the type and the premises
to the equation, with every step marked. Started 2026-10-29 by Claude on Can's "Ok go ahead".
It is a draft with holes, on purpose: the holes (§5) are the queue. Edit in place; git holds
the history.

**3 Nov, information obstruction (Codex).** T41 proves that no measurable estimator
satisfies literal L11 over the unrestricted reading type: its universal promise
on divergent-energy histories would imply consistent discrimination of two
proof mixtures that remain insufficiently separated. The current full premise
list therefore has no model under those quantifiers. L11 is strong; its full-type information-growth condition
is unsorted (G25, O19). T42 characterizes the flat Gaussian Bartlett procedure's
attainable scope; T53 gives the matching estimator-independent threshold
for a specified uniform Gaussian energy class. T52 rules out even divergent
oracle-information sums as a universal unknown-allocation guarantee.
The conditional chain below cannot establish existence
or uniqueness until this conflict is sorted. No premise or ruling is replaced.

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

- **(A) ours30 satisfies the premises.** Fails as currently stated (§3, G22).
  T41 additionally proves literal unrestricted L11 unattainable by any estimator.
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
*Gaussian* and not another second-moment family was attributed to two routes in T23.
**B05 (30 Oct, Codex): neither establishes the forcing claim.** Z admits elliptical
non-Gaussian laws and even an independent-channel non-Gaussian counterexample; the
zero-line integral does not imply full spherical symmetry. Under current ours30,
positive-misfit weights also give a Gaussian-world counterexample to literal L11
with fixed nonzero weak signal and increasing N. T23 writes out both arguments and
names their tests. The narrower fixed-N amplitude consistency result includes
non-Gaussian bounded-scale mixtures. The Gaussian remains the model stipulated by
L13; its uniqueness from L20 and L11 is a **gap**, not "forced (test only)". O4 remains
open; no premise row or equation changes here. → gap G1.

**The frame for steps 1–5 (T36, 29 Oct).** The reading space $W_i$, with the noise metric
$\Sigma_i^{-1}$ forced by units; the hypothesis $m$ as the subspace $L_i(m)$ (step 1's candidates);
$H_i$ as projection onto it and $K_i=E_i-H_i$ as the residual; $k_i$ residual dimensions per reading,
$\nu_i=k_i-1$ when the reading picks its own line. The raw law reads one residual,
$\sum_iK_i(m)$. T36 is geometry and a dimension count; it removes no hole, and it states G13's and
G3's questions (what the power on $E_i$ and η are for) in the form in which they are now asked (O14).

**Step 3. The candidate weight, and the per-reading factor.** Integrating the candidates out
needs a weight $w(\alpha_i)$ on the true acceleration. The archived argument uses:

- L8 (conditional form, ruled R-2026-09-29-9): components independent given the reading's
  shared state, so $w$ is a product over components;
- L9 (strong): $w$ depends on $\alpha_i$ only through its whitened length;
- a proper density on $\mathbb R^{k}$, $k\ge2$, that is both a product and rotation-invariant is
  Gaussian, $e^{-\lVert v\rVert^2/2\omega^2}$ in noise coordinates (Herschel–Maxwell);
- L10 (ruled R-2026-09-30-6): no size is declared, so no $\omega$: the weight is the
  $\omega\to\infty$ endpoint, flat in the vector;
- the archive's strong-signal width condition, with L13, removes the improper members
  $e^{+\tau\lVert v\rVert^2}$ and fixes the coefficient $\tfrac12$. **B07 finding (Codex, 30 Oct):**
  literal L7 specifies point calibration and MAE continuity, not this width condition;
- T7 (from L11): nothing attached to a reading may depend on $m$, which puts the weight in noise
  units and drops its $m$-free constant.

With T7's noise-volume convention, integrating against $dv$ gives a factor proportional to
$e^{H_i(m)/2}$; physical $d\alpha$ would add $\det(C_m^\top\Sigma_i^{-1}C_m)^{-1/2}$.
$H_i$ is the explained energy of `equation.md` (T6, T8's per-reading factor).
**B07, 30 Oct:** T06 now contains the characterization proof, exact Gaussian integral and
coefficient test. The flat coefficient is given explicitly by L10 as written, but its
independent derivation from L7/L13 remains *gap*: the width condition attributed to L7 is
absent from its row. The tests verify the conditional width calculation, not that missing
premise. T36's measure/π question is unchanged. → G2.
**Two further things it does not cover.** (i) Herschel–Maxwell needs $k\ge2$. The 1D
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
`::test_swap_gives_the_reciprocal`. *forced* from L3, L16. **B10, 31 Oct (Codex):**
T22 now writes out the covariance-congruence proof, including the declared subspace,
and an exact diagonal-discarding counterexample; T05 handles the two swap coordinates.
This part of G15's writing is done.

### 2.2 Joining readings

**Step 5. The product.** The joint factor of a series is $\prod_iE_i(m)$ (T12). Rests on L13
(errors independent across readings: the joint density factorises) and L19 (no borrowing: no
latent quantity shared across readings besides $m$, so no hierarchical or pooled factor). *forced*
from L13 and L19. **B10, 31 Oct (Codex):** T12 now writes the product proof with the prior
counted once, compact Gram-state scope, common-candidate derivative state, and a two-reading
example with $\eta=32/225$. The union recomputes $\eta$ and generally is not the product of
its separately tempered laws. **A premise
is doing silent work**: L13's "independent across readings" is a declaration about the instrument
(errors shared across readings have no slot, O16). As a forcing step this is fine: the function
is defined on what is declared. → noted, not a gap.

**Step 6. The starting law $\pi$, declared once.** $\pi$ enters once, not per reading (T7, from
L11: a per-reading prior would make the law depend on $N$ at fixed total signal). It is declared
before data (T20). Its family is forced: symmetric under $m/c\leftrightarrow c/m$ (L4), covariant under
units (L3), proper (T27's reason: an improper $\pi$ gives no law on null input), so
$\pi_c(m)=m^{-1}g(\ln(m/c))$ for nonnegative even $g$ of integral one, within the
absolutely continuous class (O6). For a general starting law, which the type also permits,
the exact form is $M=ce^X$ with $X$ any symmetric probability law; absolute continuity is
not forced by those three conditions. **Its shape within the family is an
input** (R-2026-10-24-3: the readings cannot decide it; the half-Cauchy, uniform in the angle the
noise defines, is the principal branch, T28). **Its centre is an input**; $c=s$ for one noise pair
is the principal branch (L3′, O11, the units argument of 10 Oct), *not ruled*. *input*, with the
family *forced* by the three named conditions. **B08, 31 Oct (Codex):** O6 now contains
the unit and reciprocal Jacobians, necessity and converse, the general-law qualification,
and what counted once does and does not add. T28 derives its uniform-angle half-Cauchy
member under the historical model without selecting it from the family. The three-member
null tests are algebraic because `met/full.py` has no general-prior API. G8's writing is
done; the centre question G9 and O6's principal-member question remain open.

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
channels, is **proved** (T38, 31 Oct: B02, `tests/test_positivity_supremum.py`): a supremum not attained
at any finite signal size, splitting exactly as $\tfrac12(k-1)+\tfrac12k$. Under a general covariance the
value is **unknown**. → gap G5 (general $\Sigma$); G4 closed.
**30 Oct:** with Σ the covariance, the tail coefficient is $\nu+2\alpha>k+1$; $k-\tfrac12$ needs a noise size
with no mean (T37 (2)). So the meaning ruled for $\kappa$, the most misfit the declared noise produces,
is the floor's meaning (§4). → G18; T38 is true either way, but it is load-bearing, and G5 worth
deriving, only if G18 keeps the ruled meaning.

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
**B10, 31 Oct (Codex):** T13 proves $0\le H_i\le E_i<\infty$, so the multiplier is bounded
above and below by positive constants. A proper starting law gives a proper output law
and a finite positive median, including the all-zero case. This proves totality for the
current median convention without using the old cloud readout.

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
For a physical channel exchange the centre also becomes $1/c$, giving $m\mapsto1/m$.
The fixed-centre form instead uses the calibrated channel swap
$(F,a)\mapsto(ca,F/c)$, with its covariance transformed too, giving $m\mapsto c^2/m$.
**B10, 31 Oct (Codex):** T05 now proves the image-law and median identities for both forms,
without T4. For $\eta$ it states the needed correspondence of selected maximizers
(a unique maximizing projective line suffices); it does not resolve an arbitrary tie
rule. The current half-Cauchy gives a unique median; the proof also states the
equivariant median-interval convention. This checks the declared median's symmetry,
not the choice of median itself, which stays conditional on L23 at step 16. → G7's
writing done; the remaining parts of G15 are unchanged.

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
| L11 (consistency) | **fails** (B05, 31 Oct): at fixed per-reading signal as N grows the point converges to a wrong mass (truth 2, signal 1: 1.73 in 1D, 1.81 in 3D), because $g(Q_+)$ is correlated with the score. Holds for fixed N as the signal grows (T23). T17 is the flat weight's → G22 | `tests/test_theorems.py::test_t23_gaussian_ours30_fails_literal_total_signal_consistency`; `tests/test_weight_consistency.py` |
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
  of $\pi$), 9's ends, 11 (T31, compressed), and some of the consistency theorems (19 of 38 theorem files still say "carried from the records", ten
  of them on the path; B10 wrote T5, T12, T13, T22). T16–T19, T24, T25 written 31 Oct (B06), with
  scope conditions; T17 holds for fixed $N$ only (G22).
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

**And ours30 fails L11 (31 Oct, B05; G22).** The sign of F·a entered the weight through $Q_+$
(step 8, R-2026-10-29-1, to meet T30 on anti-parallel pairs). A weight that sees the sign is
correlated with the score, so at fixed per-reading signal the point is inconsistent; a weight that is a
function of the Gram eigenvalues alone (ours29's) is not. So T30 with R-2026-09-07-1 and L11, as
stated, cannot all hold for a power-form weight: steps 8–10 are open on this ground too, independently
of G18. The statement above holds only with L11 read as fixed-$N$ amplitude consistency (T23).

**3 Nov (Codex): ours32 repairs this population bias, with a different L11 limit.**
T39 shows that its $Q_+$ is a mass-independent normalizing constant: its criterion
is $-\sum_i(\nu_i+k_i-1)\log(\nu_i+K_i(m))/2$. Under correctly whitened Gaussian
scalar-scale-mixture noise its population optimum is the true mass. With iid
fixed nonzero signal and finite-mean noise scales the continuous ours32 law and
median concentrate there as $N\to\infty$. However, one allowed 1D reading at
$\nu=1/2$ has limiting median $1.8219669\ldots$ when the true mass is 2 and its
amplitude tends to infinity ($c=s=1$). This is a separate tail-integrability
failure of literal L11, not G22's directional weighting. No input restriction
or new equation is adopted. → G23; O18; T39.

The stronger T39 example uses two 1D readings at ordinary ν = 10 with fixed
opposite perpendicular errors of size 3. Their true line remains the global
mode, but $\eta_R\to1/90$ makes the combined tail index $2/9<1$ and the median
tends to $1.3107868\ldots$ instead of 2. Without η the same law concentrates.
The failure persists on a positive-probability Gaussian error set. Thus the
finite-N global power with the current median conflicts with L11 on the cut factors; a lower
bound on the declared ν alone cannot fix it. The repeated-signal asymptotic
Bartlett argument remains valid in its named regime.

**3 Nov, continuation (Codex): conditional repairs exist, with unforced global structure.**
T40 studies corrections $\exp[-F_{\eta,a}(D(m))]$ of the total deviance
$D=\widehat L-L$. Two different semigroups have the same local slope η and
raw far-deviance rate, and conditionally preserve the concentrating raw law
on T39's ν = 10 family. A fixed nonlinear map does not thereby have the
correct Gaussian local limit: local deviance is O(1), not infinitesimal.
Letting an additional scale a grow with repeated readings, while it stays
bounded on fixed-N amplification, conditionally supplies both limits and
preserves a bounded-null-tilt reference under stated null scalings.

These are compatibility witnesses, not the next adopted equation. Their
family, scale and observational state are not selected. Continuous Z-safe
scales can be computed from individual factor oscillations without a latent
size law; L19's permission for this aggregate is still unsorted. Full
Gaussian calibration is a benchmark studied here, not a new ruling about
property N; η's cap can make its reference limit conservative. → G24.

**3 Nov, continuation (Codex): the unrestricted consistency target is itself impossible.**
T41 uses correct independent Gaussian errors and allowed deterministic signals
of length $b/\sqrt i$. Their quadratic energy diverges, but their quartic
energy is finite. Two sign-mixture sampling laws remain mutually contiguous;
universal consistency for every fixed sign stream would imply consistency
under both mixtures, a contradiction. L16's independent sign-convention
invariance strengthens the example to the all-positive amplitude schedule.
The proof applies to every estimator, independently of the law and readout.

This identifies the false implication in literal L11: unbounded total push
need not mean unbounded distinguishing information about mass. T17/T39/T40's
named regimes remain valid, as do G22/G23's separate failures there. G24
cannot repair this universal obstruction. A sufficient Gaussian scatter
condition illustrates attainable growth without selecting a new premise,
gate or input. Because L11 is strong, this is evidence against consistency
of the current premise list, under GOAL's guard. → G25, O19.

**3 Nov, further continuation (Codex): a sharp attainable Gaussian regime.**
T42 proves that the flat Gaussian factors with their actual Bartlett cap and
fallback have law concentration and median consistency precisely when
$\Delta_N^2/(N+\Delta_N)\to\infty$ as N grows, with a fixed regular prior,
common noise ratio, independent Gaussian errors and arbitrary deterministic
signal amplitudes/signs. No iid signal or leverage condition is needed for
this consistency result. The converse gives a finite random-precision law
and a wrong-angle event of positive probability; it proves median failure
as well as failure of law concentration.

The observable capped angular precision is $\eta\mathcal H$, not
$\mathcal H^2/J$ alone. It can diverge even while η vanishes. A scaled-deviance
comparison extends T40 to $a_N\eta_N\to\infty$; the stated oscillation scales
make both research families TV-equivalent to the Gaussian Bartlett law for
arbitrary deterministic arrays as N grows. Their global choice and L19
permission remain G24. The theorem characterizes this procedure rather than
selecting the full-type replacement for L11 (G25); T39's finite-ν amplitude
failures remain. No gate, physical input, production equation or ruling changed.

**3 Nov, calibration and cut continuation (Codex).** T43 separates a universal
Gaussian fitted-point limit from fitted-J variance estimation. A sequence
with one squared signal amplitude $4N$ and N − 1 null signals meets T42's
growth condition, yet its quoted 95% interval covers only 0.9272811852 in
the limit. A Z-safe observed-spectrum control reaches nominal Gaussian
asymptotic coverage without a count input. Its finite-N body remains
unforced. G26 records the leverage condition missing from a universal
fitted-J argument, rather than treating the problem as finite-N only (G3).

T44 separates factor averaging from a frozen-noise modular cut. Averaging
factors and then normalizing weights the noise law by its component mass
normalizer; preserving the specified noise marginal instead requires a
mixture of normalized whole-series conditionals. That control fixes the
single-reading ν < 1 amplitude failure, with a different temperature
placement, but reintroduces the positive-misfit sign bias at repeated weak
Gaussian signal. G19 therefore includes normalization order as well as
power form. Neither the cut interpretation nor a new calibration is adopted.

**3 Nov, conditional-score and positive-weight continuation (Codex).**
T45 extends the archived directional score test to arbitrary declared
Gaussian covariance and subspace. Conditioning on each reading's own
candidate projection gives an exact common-mass score variance, without
inferring a population of signals. Inversion is Z-safe and has exact
pointwise coverage in every nonempty Gaussian experiment. Locally it gives
T43's spectral width, but globally it can be disconnected or unbounded:
this is a confidence set, not the required proper-law central interval.
O15's model-specific corollary is also proved directly: at every fixed N,
proper-law central intervals have infimum coverage zero over arbitrarily
weak signals. That corner does not invalidate identifying-regime calibration.

T46 shows that T44's weak-signal bias is generic for the stated fixed
positive-misfit quadratic-weight class. Every nonconstant monotone w(Q+)
has a positive Gaussian-null weighted crossmoment and an explicit nonzero
identifying signal with the wrong population line. The corresponding
normalized-cut mean-precision class has the same limit under its stated
joint LLN. Changing only that monotone weight cannot repair the class;
nonlinear residual factors such as ours32 remain outside the theorem.
For the uniform angular prior, even a nonconcentrating nonnegative global
temperature cannot repair the weighted-quadratic median in the explicit
witness: it stays between the wrong fitted angle and the prior's centre.
G19/G22 and the global-law selection in G26/G3 stay open. No equation,
output, premise or ruling is replaced.

**3 Nov, CDF/residual continuation (Codex).** T47 gives a probability
recalibration alternative to a global likelihood power. A law-only operator
natural under every increasing reparameterization, with a full Gaussian
variance-inflation benchmark, is uniquely the normal-score CDF warp within
that class. It preserves the raw median exactly, retains bounded-positive-eta
concentration, and gives regular iid Gaussian calibration. Vanishing eta
needs moving-tail control; ordinary BvM and diverging corrected curvature
alone do not imply it. The operator-class reduction is not supplied by L23
or L17 and remains unforced.

T48 derives the Student residual exponent (nu+k)/2 from fully normalized
Gaussian integration under a variance-independent declared-metric latent
measure. It repairs the raw nu<1 single-reading amplification failure and
the stated fixed-N regime. Bounded relative evidence also survives for a
fixed incompatible reading with a positive residual floor. The same outer
Bartlett power still cancels the prefactor change exactly on T39's pair.
Using T47 instead provides a conditional proper-law construction compatible
with these named regimes and regular iid Gaussian uncertainty. The measure,
ambient conditioning and calibration class are still unsorted G19/G24;
G26's full-type variance target and composition arguments remain open. No
definition, input, output, premise or ruling is replaced.

**4 Nov, operator/subspace continuation (Codex).** T49 derives canonical
graph metric and volume with arbitrary fixed covariance and distinguishes
the intrinsic scalar curve from its ambient matrix chart. The former gives
a dimension-independent half-Cauchy under common uncorrelated noise; the
correlated scalar construction changes shape while retaining median s.
Different shrinking tubes around the Haar-null scalar curve give different
limits. Physical frame covariance does not supply the full active symmetry
that would select ambient Haar, and neither graph volume nor noisy score
rank selects a latent-vector measure or an input prior.

T50 connects the common-isotropic Gaussian profile to standard TLS/Bingham
geometry and extends T45's scalar conditional score to a matrix tangent
family. Its exact vector statistic uses observed covariance rank and
respects zero padding. A scalar tangent recovers the scalar benchmark;
using every operator tangent tests a different family. The comparison
assumes full signal vectors in a common known frame and Gaussian world
errors. It supplies an uncertainty benchmark, not a proper matrix law or
a change of physical inertia. O6/G19/G26 remain open with these distinctions
made explicit; no premise or production body is changed.

**4 Nov, composition and consistency scope (Codex).** T51 proves that a
universal law-only CDF warp commuting with independent addition or parallel
sum must be the identity. This is a conditional result about an extra axiom,
not a consequence attributed to L5/L14. Their provisional recorded scope is
pushforward of a supplied joint law and coverage consistency of a new
measurement. Correct univariate calibration leaves cross-body dependence
unsettled. In the stronger centered-Gaussian benchmark, exact nominal
coverage of both parts and their sum forces the cross-covariance; conservative
coverage alone does not. O2/G27 records the already owed joint/protocol
construction without importing universal convolution commutation.

T52 strengthens the information obstruction using unknown fixed positive
sparse supports, without sign invariance. Almost every support under one
proof measure has infinitely many fixed-amplitude signals and divergent
oracle information, yet the two mass mixtures retain a uniformly positive
testing error. T53 supplies a positive counterpart: uniform point consistency
over a fixed compact interior angle interval and arbitrary Gaussian
exact-energy nuisance vectors is possible precisely when
Delta_N^2/(N+Delta_N) diverges. This is now an estimator-independent sharp
threshold for that class, rather than only T42's method-specific threshold.
It does not characterize every individual fixed stream. O19 proposes
separating fixed identifying noiseless limits from explicitly quantified
long-series classes; the full-type statement remains unsorted. No premise,
grade, type input, production equation or ruling is replaced.

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
goal: Claude), **goal** (the sort turns on what the goal means and the page does not say: a
question to Can, asked as one; `bounds.md`), **VD-Newton** (a rung to write there: Can). Order is
the chain's. (Until 7 Nov the fourth type was **rule (Can)**; R-2026-11-07-1: the deductions are
Claude's to make, Can clarifies the goal.)

Rows G23–G27, and the clauses of G8, G19 and G22 that cite T39–T53, were harvested 7 Nov from the exploration
branch `codex/ours32-consistency` (Codex, 3–4 Nov): input, not checked by Claude. The branch's later rows and clauses
(the joint, score-law and anchor line after T53) were left on it.

| id | step | hole | type | where |
| --- | --- | --- | --- | --- |
| G10 | 1 | L22, scalar mass, as a clause of L1 | sort (Claude); goal if the version's L1 does not settle whether inertia is a scalar | O17 |
| G11 | 1, 4 | $N_i$ as an input (L18's second clause); R3 unwritten, R5 provisional | sort (Claude: TASK §5's inputs-cost test); VD-Newton | L18; `ladder.md` |
| G1 | 2 | T23's Gaussian-forcing claim: B05 delivered counterexamples to both routes as stated; further justification is open. **Merged 31 Oct;** L13 is now stipulated, conditionally forced by Herschel–Maxwell given full sphericity and independent channels | write; derive; rule | O4 → brief B05, merged |
| G2 | 3 | T6's argument: **written 31 Oct (B07, merged)**. Its remainder: the coefficient ½ (excluding the improper $e^{+\tau\lVert v\rVert^2}$, $0<\tau<\tfrac12$) does not follow from L7, which calibrates a point (B07 found the brief mis-sorted); L10 states $e^{H/2}$ as a clause. **Claude, 31 Oct, provisional:** by T36's homogeneity the $\tau$-member is the flat factor under a declared $(1-2\tau)\Sigma$, so the coefficient is forced exactly when the factor may not alter the declared noise law (R1 with L24, candidate); it rides on G18, not a new fork | write done; the remainder attaches to G18 | T06 §3; L24 |
| G17 | 3 | the 1D factor: "a reading declared $k=1$ is the 1D reading" needs stating as a premise or deriving from L8 and L18 | sort | T25, L8 |
| G15 | 4, 5, 15, and the consistency theorems | 31 theorem files (19 by 4 Nov: B06 wrote T16–T19, T24, T25, with conditions T16, T17, T25's statements lacked; B07 T06; B08 T28; B10 T5, T12, T13, T22. Of the 19, nine are comparison-readout only (B09); **ten are on the path: T1, T2, T7, T8, T20, T26, T27, T30, T31, T33**, of which T30 and T31 are sorts, G13 and G3) with the proof "carried from the records": write them out in the chain's order, shortest first (T22, T24, T25, T12, T19, T17, T16, T18, T13) | write; the pure algebra ones can be briefs | theorems/ → briefs B06 (T16–T19, T24, T25), B10 (T5, T12, T13, T22) |
| G8 | 6 | **Written 31 Oct (B08, merged):** L3, L4 and properness give the even log-mass family; no member selected. **Claude, provisional:** L11 adds a full-support clause. T49 gives the half-Cauchy as canonical scalar-graph curve volume in common uncorrelated noise and a correlated extension, but choosing that reference remains additional. Ambient Haar cannot be restricted to the scalar curve by density substitution; existing frame covariance does not select the larger active group | family derivation done; geometric reference conditional, no principal member selected | O6, T28, T49; `tests/test_starting_law_family.py`; `tests/test_graph_geometry.py` |
| G9 | 6 | the centre: L3′'s first half, then $c=s$ | sort (Claude); goal: whether "no knowledge of the object" is part of what the definition promises | O11 |
| G13 | 7, 9 | why the weight is a power on $E_i$, and which premise makes bounded confidence a requirement | sort | T30, L21. **30 Oct, provisional:** the power is a noise scale read once (L21's cut, O14); bounded confidence derived from a noise size with a power-law tail (T37 (1)), which supersedes the morning's "from Σ as a floor" |
| G18 | 9, 10 | R1 and its kin (the σ's believed) against O12's cap and κ's meaning (R-2026-09-30-10, R-2026-10-29-1) | sort (done, provisional: O12's meaning is the mis-stated one) → **open: stands until forced** (R-2026-11-01-1; not a ruling) | §4; O14; T37 |
| G19 | 7 | L21's power form, T37's factor average, normalized whole-series cut (T44), or normalized Gaussian residual integration (T48), as the body. T44's frozen-marginal correction reintroduces T46's positive-misfit bias. T48's Student exponent depends on a specified latent volume and retained ambient residuals. T49 shows graph-hypothesis volume is a different measure and cannot select the latent noise-scale normalization. B17 refutes the proposed T36 forcing: a scalar V^b weight need not change the metric, and b=-1 also has unit residual total. Fixed-V normalization with the declared scale marginal would select b=0 but needs an additional rule/sort. Relative normalization at fixed V, mixing V and temperature placement are not interchangeable | **sort/rule for the latent measure; open: stands until forced** (R-2026-11-01-1); cut, latent measure and retained data/model scope remain unforced; no adoption. **4 Nov, a candidate for what forces it (Claude, not a sort):** L11 with T30. The power form with a sign-blind $g(Q)$ fails T30 (ours29); with a sign-seeing $g(Q_+)$ it fails L11 (G22). The cut factor sees the sign through $H^+$ and $Q_+$ and still keeps the population optimum at the truth (2.00 ± 0.01 at signal 1–2, 1D/3D; `tests/test_weight_consistency.py`) | T37, T44, T46, T48, T49; O18; `tests/test_cut_normalization.py`; `tests/test_positive_weight_obstruction.py`; `tests/test_residual_mixture.py`; `tests/test_graph_geometry.py`; `tests/test_latent_measure_exponent.py` |
| G20 | 9 | T37's tail index α: an input, but TASK §5's guard asks that practice supply it. **Sorted 30 Oct, provisional (T37, G20 section):** α = ν_eff/2, the effective degrees of freedom of the stated uncertainty (GUM Annex G), on calibration certificates; the inverse-gamma shape forced by L13 applied to how σ was obtained. Practice normalises the precision ($E[1/v]=1$), which restates T37 (2) and (3): ours30 is a σ declared with ν_eff = ½; the cap holds exactly in 1D | sort (done) | T37 |
| G21 | 9, 11 | L24's test: which temperings are cut posteriors of a mean-1 noise size; η's series version | derive | L24 → brief B12 |
| G22 | 8–10, 11 | **The sign in the weight against consistency** (31 Oct, from B05). T30 with R-2026-09-07-1 put the sign of F·a into the misfit ($Q_+$); L11 forbids a weight correlated with the score at the truth. Eigenvalue-only weights keep the population optimum at the truth (T38's step 1 symmetry); $g(Q_+)$ moves it (2D/3D: −10%/−9% at signal 1, −3%/−2% at signal 2). **Not inherited by ours32 (4 Nov):** T37's cut factor keeps the optimum at the truth (G19). **Also at zero signal (4 Nov, B13 merge):** the weighted cross-moment is positive, so ours30's pure-noise law concentrates on $m=c$ as $N$ grows (O15). T46 proves the scalar obstruction for every fixed nonconstant bounded monotone w(Q+) quadratic weight, and a separate analytic argument forces constant w for universal Gaussian stationarity even without monotonicity. Nonlinear residual factors remain outside this class. Forced, input, or which premise is mis-stated | sort (Claude) → rule; the fixed-quadratic class obstruction is proved, not a selected repair | T30, T32, T46, L11; `tests/test_weight_consistency.py`; `tests/test_positive_weight_obstruction.py`; T23 |
| G23 | L11, finite-ν probes; G3's finite-N power | **Tempered tail index against amplitude consistency** (3 Nov, Codex). T39 proves raw nu<1 and eta-caused nu=10 failures for current factors. T48's specified Student residual measure repairs raw fixed-N amplification, but its prefactor change is canceled exactly by the outer Bartlett power on the opposed pair. T47's median-preserving CDF warp instead retains concentrating raw laws when eta has a positive lower bound. Together they give a conditional compatibility witness, not the selected replacement or unrestricted growth theorem | failures and conditional named-regime repair proved; measure/class selection unsorted; no premise changed | T39, T47, T48; O18; `tests/test_cut_consistency.py`; `tests/test_cdf_calibration.py`; `tests/test_residual_mixture.py` |
| G24 | G3/G23, global calibration shape | **Local η does not select a global correction**. T40's scalar-deviance semigroups need unforced experiment state; T42 extends their stated Gaussian TV comparison to vanishing eta. T47 derives a unique normal-score warp under a law-only/full-reparameterization class and full Gaussian benchmark, without an extra scale. That class reduction is unforced, and vanishing eta needs raw-tail control beyond ordinary BvM. T51 rules out any nonidentity universal CDF warp under an additional independent-composition interchange axiom, which is not recorded VD scope. Research state and class assumptions are not newly licensed physical inputs under TASK §5 | conditional class uniqueness, calibration and limitations proved; stronger composition route closed; global selection and permitted state unsorted; no adoption | T40, T42, T47, T51; O14, O19, O2/G27; `tests/test_deviance_calibration.py`; `tests/test_gaussian_information.py`; `tests/test_cdf_calibration.py`; `tests/test_composition_calibration.py` |
| G25 | L11, identifying growth | **L11 restated 7 Nov (Claude, provisional): $\widehat m_N\to m$ when $\Delta_N/\sqrt N\to\infty$**, the literature's rate (B19), = T42's and T53's threshold; the old wording struck in its slot. Sorted as a wrong premise (the deduction misread F1, not the goal). Still open: the full type ($k>1$, general $\Sigma_i$) and the phase-two bin. The branch's text follows. **Latent energy or oracle information alone does not justify universal consistency**. T41 disproves literal L11; T52 extends the obstruction to fixed positive sparse supports with every positive-power energy sum and known-magnitude information divergent, without L16. This challenges a strong premise and the consistency of the present list. T53 proves T42's energy threshold sharp for all estimators under specified uniform Gaussian energy-class quantifiers. It is not a criterion for each individual stream or the selected full-type replacement; T39's finite-ν failures remain | contradiction and sharp uniform Gaussian threshold proved; proposed experiment-scoped obligations recorded in O19; full-type statement unsorted; L11 statement and grade retained | T41, T42, T52, T53; O19; `tests/test_signal_identifiability.py`; `tests/test_gaussian_information.py`; `tests/test_sparse_identification.py`; `tests/test_energy_minimax.py` |
| G26 | 11, long-series uncertainty | **Raw fitted J need not estimate sampling variance under concentrated leverage** (T43). A Z-safe spectral control attains Gaussian asymptotic coverage; T45 gives an exact full scalar-type Gaussian score benchmark. T50 extends this to vector scores in a specified common-frame/full-subspace operator comparison, using observed covariance rank. These exact sets are not proper-law intervals, and noisy conditional rank does not establish identifying growth. T49 separates the signal-independent graph metric from Fisher/score covariance. The selected law variance and finite-N body remain unsorted; marginal calibration alone does not settle the joint obligation G27 | counterexample, Gaussian attainability and exact scalar/operator tests proved; proper-law selection and full-type asymptotic law calibration open; no adoption | T43, T45, T49, T50, T51; T18, T31; O14, O15, O17, O2; `tests/test_leverage_calibration.py`; `tests/test_conditional_score.py`; `tests/test_operator_score.py` |
| G27 | L5/L14, joint composition protocol | **A supplied joint law has a unique pushforward; calibrated marginal laws do not supply that joint**. This makes O2's existing bridge obligation explicit: multiple-body identities, shared/separate observations and cross-covariances, propagated joint prior, and the coverage protocol for a new composite measurement. T51 separates this from the unrecorded universal calibration-interchange axiom and derives a cross-covariance constraint under an explicitly stronger centered-Gaussian benchmark. R3/R5 remain provisional/unwritten as in G11 | joint/protocol construction open; stronger interchange route sorted; no new premise or type slot | T51; O2, G6/G11/G26; R-2026-10-28-1; `tests/test_composition_bridge.py`; `tests/test_composition_calibration.py` |
| ~~G4~~ | 10 | ~~$S_k=k-\tfrac12$, and $\sup E[Q]=k-1$, proved~~ **closed 31 Oct: T38** (B02); load-bearing only if G18 keeps O12's meaning | derive | brief B02, merged |
| G5 | 10 | $\kappa$ under a general $\Sigma_i$ (correlated channels) | derive, after G18 (only if it keeps O12's meaning) | O13 |
| G3 | 11 | $\eta$'s cap at 1 and its power form at finite $N$: forced, input, or a mis-stated premise (property N stated in the limit only?) | sort | O14; brief B04 supplies the finite-$N$ facts. **30 Oct, partly (O14):** the cap is the series floor and the power form is a scale; the finite-$N$ form is still open. **30 Oct (O15):** property N is stated in the limit, and ours30 meets it there (0.93–0.95 at $N=80$–200, signal 1–2, 3D); only the finite-$N$ form is left. **31 Oct correction:** that is finite-$N$ only; ours30's point is inconsistent at fixed signal (G22), so its coverage must fall as $N$ grows. **31 Oct (B04 merged):** the leverage identity $E[J]=C(1-1/n_{\rm eff})$ holds for fixed nonnegative curvatures only; ours30's observed curvatures are signed, so $n_{\rm eff}$ is not a count and the target has no finite mean. Standard remedy (HC2 leverages from Fisher, not observed, information) is a pointer, unsorted. **4 Nov (O14):** η's target identified: on long series η → Δ/(Δ+M), which is exactly the empirical-Bayes shrinkage $\hat\beta=1-M/H^2$ of ours28b's declared shared size law (Lindley–El-Sayyad); so η is that shrinkage recovered from the identity's failure rather than from a declaration. $\hat\beta$ itself fails Z (a zero reading adds $k$ to $M$). Constraint, provisional: the finite-$N$ form must estimate Δ/(Δ+M) continuously in the data; the shared/per-reading fork (the leverage case of Codex's theorem 43, branch `codex/ours32-consistency`) is where the form is open |
| G7 | 15 | ~~T5 for the median: write the law's swap symmetry~~ Written in T05, with physical/calibrated maps and maximizer-selection condition explicit | write (done, 31 Oct, Codex; merged 4 Nov) | T05; `tests/test_short_law_proofs.py` → B10 |
| G6 | 16 | L23 | sort (Claude; a draft sort of 6 Nov, L23 against mass's ratio-scale type, is on branch `claude/b21-brief`) | O2 |
| ~~G12~~ | 4 | ~~the audit and the problem files still count L12, O5, T9–T10's freedoms as the definition's~~ **closed 4 Nov (B09):** the nine theorems, L12, O5 and the two ladder rows marked comparison-readout only, not removed; `tools/audit.py` tags T21 as not load-bearing | write | brief B09, merged |
| G14 | — | agreement with experiment (`../TASK.md` §5): a check on the survivor, not on a probe; no venue yet. The robot is a phase-one testbed, not this check (R-2026-10-31-1); B01 merged 31 Oct as a robot input check | derive (a study), later | — |

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

- 2026-11-04 (composition and consistency continuation): Codex, T51 separates joint-law propagation from stronger calibration interchange and gives a conditional Gaussian covariance constraint; G27 exposes the existing O2 joint/protocol obligation. T52 rules out unknown-allocation oracle-information promises; T53 gives the sharp uniform Gaussian energy-class threshold. O19 records a proposed experiment-scoped contract without changing L11, grades, type or rulings.

- 2026-10-29: created (Claude). Eighteen steps, sixteen holes.
- 2026-10-29 (later): point-level Z removed as a premise (R-2026-10-29-2); steps 13, 16, 17 and §3 reworded.
- 2026-10-29 (night): the frame for steps 1–5 named (T36), on R-2026-10-29-3.
- 2026-10-30: steps 7, 9, 10, 11 and §§3–4 brought up to T37 and L24: the rulings at steps 9–10 conflict (R1 against O12), sorted provisionally; what is proved restated; holes G18–G21.
- 2026-10-31: B01, B02 (T38, G4 closed), B04 and B05 merged; ours30 fails L11 (G22, new); L13 stipulated.
- 2026-11-03: Codex, T39: ours32's residual identity and repeated fixed-signal consistency; an allowed single-reading heavy-tail counterexample opens G23. Numerical projection and whole-line-search fixes preserve the probe's continuous equation.
- 2026-11-03 (continuation): Codex, T40: a state-free local/tail obstruction and two conditional dual-regime repairs; semigroup and common tail constant still do not force a rule. G24 records the unsorted experiment state and global shape.
- 2026-11-03 (further continuation): Codex, T41: literal total-energy consistency is unattainable for every estimator over the unrestricted type; G25/O19 expose the missing identification condition. No premise or equation replaced.
- 2026-11-03 (information-growth continuation): Codex, T42: necessary and sufficient flat-Gaussian/Bartlett consistency criterion, critical random-law converse and a vanishing-power comparison for T40. G25 gains an attainable submodel; G24's choices and the full-type premise remain unsorted.
- 2026-11-03 (calibration and cut continuation): Codex, T43/T44: an identifying Gaussian undercoverage counterexample and a Z-safe spectral comparison (G26); factor averaging versus normalized conditional cuts, amplitude repair and a separate weak-signal bias (G19/G22). No production or ruling change.
- 2026-11-03 (conditional-score and positive-weight continuation): Codex, T45/T46: exact full Gaussian score tests and the proper-law corner boundary; generic fixed positive-misfit quadratic-weight bias. G26 gains a test benchmark, O15's model-specific corollary is proved, and G19/G22's class obstruction is strengthened. No production, output, premise or ruling change.
- 2026-11-03 (CDF/residual continuation): Codex, T47/T48: law-level normal-score calibration with conditional uniqueness and a vanishing-eta tail boundary; normalized Student residual integration, raw amplification repair and persistence of the outer-power failure. A proper-law compatibility witness, with measure/class and full-type forcing still open. No production, input, output, premise or ruling change.
- 2026-11-04 (operator/subspace continuation): Codex, T49/T50: canonical graph measures and scalar/tube restriction, covariance and Fisher/fiber distinctions; TLS/Bingham operator profile and an exact conditional vector-score benchmark. Prior, latent measure and law-calibration choices remain open; scalar physical type unchanged.
