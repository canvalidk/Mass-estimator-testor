# The mass estimator's properties and lemmas, 10 October

**2026-10-10.** Record, written by Claude at the user's request ("The mass
estimator lemmas are in a supersession folder right? Just copy the old lemma
list and edit the new version. That way you can write all of it as accepted
and just put a provisional above all of it. Do 1 and 2 please."). It is the
29 September record copied and edited: every row keeps its 29–30 September
state, and the 10 October changes are marked **10 Oct** inside the row.

**Supersedes `mass_estimator_lemmas_2026-09-29.md`**, which superseded the
27 September record. It binds nothing else. The records whose results it
carries are not corrected by it.

**PROVISIONAL.** Everything marked 10 Oct is written as accepted, on his
instruction, and is provisional as a whole. The rulings behind it are the
ones quoted in `mass_as_force_paired_with_reciprocal_acceleration_2026-10-10.md`
("I like the window. You can just add that in to the lemma lists"; "keep L14
provisional. The confusion and the tradeoff with the swap vs the lumping etc.
ruling is well understood and can remain unruled"). In particular L14's
status is unchanged and the L14/swap trade is unruled; O9's window criterion
and O11's units argument are written in as proposed, not as closed.

**What changed on 10 Oct** (source: the units-and-axes record, §§6–9):

- §1: a reading declares its **window**; its force and acceleration are
  averages over that one window (T33).
- L5, L14: read as the linearity of the two pairings $\mathbf a'(\cdot)$ and
  $\mathbf F'(\cdot)$. L18: the window joins the reading structure. L3′: the
  units argument fixes its constant.
- T1: restated in pairing form, with the product $P$ classifying the exact
  branch. T32 (the misfit is the pairings' disagreement in noise units) and
  T33 (same window) added.
- O9: the window criterion proposed. O11: the units argument written in.
- Sources: the units-and-axes record.

**The 30 September revisions of the 29 September record** (L8 rewritten,
L21 and T30 added, 29 Sept night; T8 to ours29, L10 restated, L10′, L4′ and
$G$ withdrawn, T31 added, 30 Sept afternoon; O3 and O8 closed, late
afternoon; O10 and O12 closed, 16:26) are folded in as that record's final
state; its header carries their wording and sources.

**What it is for.** The project's goal is the estimator function that can act
as the definition of mass. So the definition *is* this list: the properties
mass must have, the declared inputs, and what they force. The equation is
what the list produces. Every item below says whether it is **ruled** (the
user's words, dated), a **premise** (graded, not ruled), a **theorem** (with
its sources), **declared** (an input), **proposed** (not ruled), **struck**,
or **open**.

**What it folds in.** The 27 Sept lemma record (L1–L14, T1–T20) and every
correction proposed since:

- the ours26 sweep §4 (errors E1–E6, gaps G1–G4);
- the perturbation matrix §6;
- the ours28 record §8;
- the Jeffreys-check record §4 (L4′);
- the betweenness record §4;
- the full-covariance / declared-subspace record §§2–5;
- the no-borrowing record §5;
- the C3 narrow-form record §9;
- the reference implementation record;
- the rulings of 27–28 Sept;
- the two findings of 29 Sept (zero components contribute nothing; the
  readout fails Z).

**Numbering.** L1–L14 and T1–T20 keep their numbers so old references still
work. Struck items keep their slot. New items start at L15 and T21.
10 Oct adds T32 and T33; no L is added (the window enters §1 and L18).

**Grades** (premises only). *bedrock*: denying it changes the problem.
*strong*: a real premise, well motivated. *medium*: plausible, narrower.
*mild*: technical.

---

## 0. The frame (inherited, not re-argued here)

| | statement | status |
|---|---|---|
| F1 | **The estimator is the definition.** What mass is equals how it is read from data; there is no further mass behind the readout ("there is no super meta mass that we can't achieve"). The goal is to answer "what are the estimators estimating?", not to build a good estimator. | ruled, 27 Sept |
| F2 | **Uniqueness is the stake.** Any choice the properties leave open counts as an alternative definition of mass, i.e. an alternative Newtonian mechanics. "The choices only differ where nobody measures" is not an answer. | ruled, Sept |
| F3 | **Inherited from `_note_why_mass_estimation`:** provenance independence, laws as predictions, the authority ruling. Named here so they are not left out silently (sweep G3); not restated. | inherited |
| F4 | **Commitment order is suspended** for this work: Newton III "and more" may be used. | ruled, 27 Sept |

## 1. The type

**Input.** A series of readings $i=1,\dots,N$, each declaring:

- data $\mathbf y_i=(\tilde{\mathbf F}_i;\tilde{\mathbf a}_i)\in\mathbb R^{2d}$, any values including zero, in the reading's own frame;
- an error covariance $\Sigma_i\in\mathbb R^{2d\times2d}$, positive definite, any correlation, including between the channels;
- a candidate subspace $N_i\in\mathbb R^{d\times k_i}$ with orthonormal columns (default $I_d$: nothing known about the push's direction);
- membership: which components belong to which reading (the reading structure).
- **10 Oct:** the reading's **window** $\Delta t_i$, over which both its force and its acceleration are averaged: one window per reading (T33). Part of the reading structure (L18).

Also declared for the series:

- the starting law $\pi$: its shape and its centre $c$;
- ~~the size law $G$ on each reading's push scale~~ (withdrawn 30 Sept: no
  size law is declared, L10); instead the middle of the misfit weight $g$
  (T30, O12).

**Output.** A law of mass on $(0,\infty)$ and a point $\widehat m$; since 30
Sept also $\eta$ and the raw law's interval beside the quoted one (ruling 3).
The misfit is kept with the result (sweep G2, d2 P3: an output clause, not a
premise the mathematics uses).

## 2. The properties of mass (premises)

| | name | statement | status | change since 27 Sept |
|---|---|---|---|---|
| L1 | **Newton II, one object** | $\mathbf F^\ast_i=m\,\mathbf a^\ast_i$ with one $m>0$ shared by every reading: one object, one inertial frame, over the series. Assumed, not tested: no rejection gate. | bedrock; "not tested" ruled 7 Sept; **persistence ruled a premise 30 Sept, 16:26** ("a perfectly fine lemma. Doesn't need proving") | Addressing and persistence written in (sweep G1); persistence no longer owed an argument (O10 closed). |
| L2 | **Noisy type** | every declared $\sigma>0$ ($\Sigma_i$ positive definite). As d2 reads it, the type promises a mass for every measured input, including measured zero–zero. | bedrock; $\sigma>0$ ruled 25 Sept; "a mass for every measured input" is d2's reading, **now in tension with O1** | "So $\widehat m\in(0,\infty)$" moved to T13 (sweep E5). |
| L3 | **Units covariance** | force units $\times\lambda$, acceleration units $\times\mu$ $\Rightarrow$ law and point scale by $\lambda/\mu$. The only mass scale allowed in is the declared $\pi$, which transforms the same way. | bedrock | none |
| L3′ | **Knowledge-free default** | with no knowledge of the object, $\pi$'s scale comes from the instruments' SDs alone; knowledge of the object enters only as old information in the prior slot. | **proposed** (sweep E1), not ruled | new; T14 rests on it. **10 Oct:** the units argument (O11) fixes the constant: given that the default declaration for mass with $(\sigma_F,\sigma_a)$ and for inverse mass with $(\sigma_a,\sigma_F)$ are one declaration read backwards, $c$ may depend only on declared quantities that invert under the swap, so one noise pair, any $\rho$, gives $c=s$. What remains to rule is L3′'s first half, that the default references the instruments only. |
| L4 | **Swap, law half** | $\pi$ is symmetric under $m/c\leftrightarrow c/m$. | strong | none |
| L4′ | ~~Push–noise symmetry~~ | **withdrawn 30 Sept** with the size law it constrained: there is no declared size law in ours29 (L10). The finding it recorded (every symmetric size law calibrated in ours28's own world, every asymmetric one failed) stands as a fact about ours28. | — | withdrawn with $G$ (ruling 2 of 30 Sept) |
| L5 | **Lumping** | $R(X_1+X_2,Y)=R(X_1,Y)+R(X_2,Y)$: glue two bodies at common acceleration and their masses add. | strong | **10 Oct:** read as the linearity of the pairing $\mathbf a'(\cdot)$ in its argument, $\mathbf a'(\mathbf F_1+\mathbf F_2)=\mathbf a'(\mathbf F_1)+\mathbf a'(\mathbf F_2)$ (units-and-axes record §9.1). The m-regression is $\mathbf a'(\mathbf F)$ and passes only this law. |
| L6 | ~~C3~~ | **struck as a premise.** Now theorem T21. | — | derived; ruled 28 Sept, 23:41 |
| L7 | **Calibration** | $R=f/\alpha$ on sure inputs; continuous at constants **in mean absolute error**; defined on every law with finite means. | mild | topology named (sweep E4) |
| L8 | **One reading is one completeness claim** (was: axis splitting, restated) | a $k_i$-dimensional reading declares that the signal has no component outside its declared subspace $N_i$; that claim imposes $k_i-1$ codirectionality conditions, whose violation $Q_i$ (the off-line residual at the best mass, inside $N_i$) is a property of the reading. Components are independent **given the reading's shared state** ($N_i$, $\omega_i$, $Q_i$). **The strong form is false:** a reading and its axes read as separate 1D readings do not give the same law, because separate 1D readings are projections, impose no condition, and cannot be contradicted. Same data, different declaration; the reading structure (L18) does the work. | ruled 29 Sept, 21:22 ("1 3D readings is not the same as 3 1D readings"; "when we make a 3D reading, we assume that the true vector couldn't have been in the 4th or 5th dimensions. And when we make a 1D reading, that same thing can be allowed") | **Rewritten 29 Sept, night.** The old strong form failed under ours28 (4.0807 against 3.7328) and is now false by ruling. ~~The conditional form forces nothing~~ (corrected 30 Sept): the conditional form, with L9, forces the Gaussian candidate family (Herschel–Maxwell), and with L10 the flat endpoint (T6); the stacked length rests on L12 instead (T9). O3 closed. |
| L9 | **Isotropy, in noise units** | a reading's push weight depends only on the push's length in whitened (noise) units, within its declared subspace. | strong; units fixed by T7 (the push-units bend dies three ways) | Units fork closed except one scenario (O8). |
| L10 | **No size law is declared** (restated 30 Sept) | the definition declares nothing about how large a signal is in noise units, per reading or per series: L19's reason (readings from any setup, in any combination) excludes a declared population of setups as much as a shared size. Each reading's factor is therefore the flat one, $e^{H_i/2}$, read as a **working likelihood**, not as the claim that every signal is infinitely large; how far it overstates the series' precision is measured on the series by $\eta$ (T31), from the readings' disagreement about $m$, never from their sizes. | strong; ruling 2 of 30 Sept | **restated.** The 28 Sept striking stands against the old reading ($R^{\ast2}=\infty$ as a claim); the learned-$\omega$ hint record's "the width problem lives in L10" is answered here |
| L10′ | ~~Push scale is a per-reading nuisance~~ | **struck 30 Sept** (ruling 2): independent $\omega_i$ from one declared law $G$ say the setups are exchangeable draws from one population of experiments, which L19's reason forbids (a quiet series had prior weight of order $\varepsilon^N$). | — | struck; see the ours29 record §2 |
| L11 | **Consistency** | $\widehat m_N\to m$ whenever the total true push $R^{\ast2}\to\infty$. | strong | none |
| L12 | **One joining rule** | $\operatorname{size}(A\cup B)=\operatorname{size}(A)\oplus\operatorname{size}(B)$, one fixed $\oplus$; a series has one size per channel per candidate. | medium | The narrow-C3 clause once proposed for here is now a theorem (T21). |
| L13 | **Gaussian noise model** | the model is zero-mean Gaussian with the declared covariance, independent across readings. This is a statement about the **model**, not a claim that the world's errors are Gaussian: the model uses only declared second moments. | **was declared, now forced**: by consistency under the flat push (matrix §4) and by Z under ours28 (Student-$t_5$ lets a zero reading tilt the law by up to 1.6%) | Moved to theorem-supported (T23). **Check with your 27 Sept line "any insistence that incoming data be Gaussian is of course wrong":** I read it as compatible, since the model is forced and the data need not be Gaussian. Confirm. |
| L14 | **Common force** | $1/R(X,Y_1+Y_2)=1/R(X,Y_1)+1/R(X,Y_2)$: reduced mass (uses Newton III, allowed by F4). | strong | **10 Oct:** read as the linearity of $\mathbf F'(\cdot)$, $\mathbf F'(\mathbf a_1+\mathbf a_2)=\mathbf F'(\mathbf a_1)+\mathbf F'(\mathbf a_2)$; the μ-regression is $\mathbf F'(\mathbf a)$ and passes only this law. **Status unchanged, by his ruling.** Noted, unruled: with the readout swap (T5) taken as a premise, L14 follows from L5, $1/R(X,Y_1+Y_2)=R(Y_1+Y_2,X)=R(Y_1,X)+R(Y_2,X)$; since T5 is now a theorem of T4, that is a trade, not a derivation (units-and-axes record §9.2). |
| L15 | **Betweenness** | mixing two candidate laws cannot give a readout outside both of their readouts. | ruled **provisionally**, 28 Sept ("i think we can accept betweenness provisionally") | new |
| L16 | **Frame independence** | each reading may be written in its own frame and sign convention; a fortiori the whole series may be rotated, or boosted in a Galilean way. | ruled 28 Sept ("Yes of course. Bookkeeping can't matter"); the rotation property was his, 27 Sept | New as a named premise; it forced ours27's full covariance (T22). |
| L17 | **Same candidate law, same mass** | whichever instruments produced a candidate law, it gives one mass: the readout's score cannot read the declared instrument facts. | ruled 28 Sept ("Sure") | new |
| L18 | **Declared inputs are correct** | the $\sigma$'s are correct (ruled 28 Sept, morning). The reading structure is a declared instrument fact like $\Sigma_i$ (ruled 28 Sept, 23:41: "yes they can be a proprty of the instrument"). **Proposed extension:** $N_i$ too, i.e. the push's declared direction structure is an input, not learned (full-covariance record §3). **10 Oct:** the reading's **window** is part of the reading structure: one window per reading, over which both channels are averaged (T33). | ruled for $\sigma$ and for the reading structure; **$N_i$ proposed**; window ruled 10 Oct ("I like the window. You can just add that in") | new; window added 10 Oct |
| L19 | **No borrowing between readings** | nothing inferred from one reading about signals (size, direction, spread) is used on another. Readings are joined only through the mass they share. Readings may come from any regime in any combination ("F = 1, then F = 50000000, then F = 0.00001"). **Reason** (his, 30 Sept): the equation must act as the definition of mass, "so it can't rely on things like a repeated experimental setup". **Scope** (ruled 30 Sept, 11:57): the *width* of the mass law may be computed from each reading's own residual (the Birge move): whether the readings support or contradict one another about $m$ is about $m$, which L1 shares, not about their signals. | ruled 28 Sept, night; reason and scope 30 Sept | scope added 30 Sept; the reason also struck L10′ |
| L20 | **Z: a reading that adds no information does not change the answer** | his wording, 27 Sept: "a reading that adds no information from other dimensions must not change the answer". Applies to the law **and to the point**. | ruled as a property, 27 Sept. **Scope ruling owed (O7).** | New as a named premise; T20 and the matrix used it unnamed. **Holds for the law (T24); fails for the current readout (T29).** |
| L21 | **Readings are evidence, just less** | a reading's evidence decreases with its off-line misfit $Q_i$: the law is $\pi(m)\prod_iE_i(m)^{g(Q_i)}$ with $g(0)=1$ and $g$ decreasing. Misfit at the best mass is evidence about the reading; misfit at a mass is evidence about the mass, and stays in $E_i$. Not a gate: $g>0$. | ruled 29 Sept, 20:46 ("The readings are evidence. Just less.") | new; see the evidence record §§2, 10 |

**Rulings that constrain the tests, not the equation:**

- **R1.** The declared $\sigma$'s are assumed correct. The question is whether
  the equation responds correctly, not whether misfit is evidence about the
  $\sigma$'s. (28 Sept)
- **R2.** Candidate pairs are hypotheticals the average runs over, not
  quantities recovered. Knowing $m$ leaves both slots unknown. (27 Sept)

## 3. What follows (theorems)

Sources in brackets. **[flat]** marks results that hold only for ours27's flat
push; **[one pair]** marks results that need one noise pair.

| | name | statement | from | status after 29 Sept |
|---|---|---|---|---|
| T1 | Exact branch | $\sigma=0$: one mass, no mass, or every mass; "undefined" only at exact $0/0$. **Pairing form (10 Oct):** compute two totals, $m_1=\mathbf a'(\mathbf F)=\mathbf F\cdot\mathbf a/\lvert\mathbf a\rvert^2$ and $m_2=1/\mathbf F'(\mathbf a)=\lvert\mathbf F\rvert^2/(\mathbf F\cdot\mathbf a)$, and their ratio $P=m_1/m_2=\cos^2\theta$. One mass iff $P=1$ and the common value is in $(0,\infty)$. No mass in two ways: **disagreement**, $P<1$ (the $d-1$ codirectionality conditions fail; perpendicular gives $m_1=0$, $m_2=\infty$), or **agreement outside $(0,\infty)$** (anti-aligned: $P=1$ and the value is negative; $\mathbf F=\mathbf 0$: $m_1=0$, $m_2$ undefined; $\mathbf a=\mathbf 0$: $m_2=\infty$, $m_1$ undefined). Every mass at $0/0$: both undefined. In 1D $P\equiv1$ and only the sign and zero conditions remain. | L1 | **restated 10 Oct** in pairing form: the chain audit's counting + positivity split, computed as two totals and a comparison. Checked on seven exact pairs (units-and-axes record §9.4). |
| T2 | Candidates | $\mathbf F^\ast=f\mathbf u$, $\mathbf a^\ast=\alpha\mathbf u$, $f,\alpha>0$ | L1 | unchanged |
| T3 | Separation | $R(X,Y)=R(X,1)\,R(1,Y)$ | L5, L14, L7 | unchanged |
| T4 | Ratio of means | $R=E[X]/E[Y]=E[\alpha m]/E[\alpha]$ | T3, L5, L14, T21, L7 | unchanged; C3 now enters as T21 |
| T5 | Swap, readout half | $R(Y,X)=1/R(X,Y)$ | T4 | unchanged |
| T6 | Push weight | L8 and L9 force the Gaussian family $e^{-\lVert\mathbf v\rVert^2/2\omega^2}$ in noise units (Herschel–Maxwell). The flat weight is its $\omega\to\infty$ endpoint, and **since 30 Sept it is the one in use** (L10: no size law, so no $\omega$ to average). | L8 (conditional form), L9, L10, L7, **T7** (sweep E2) | Back to the flat weight, now as a working likelihood tempered by $\eta$. **Supported again, 30 Sept:** L8's conditional form and L9 force the Gaussian family; L10 (no declared size, fixed or lawful) selects the flat endpoint; L7/L13 remove the improper $e^{+c\lVert\mathbf v\rVert^2}$ members; T25 carries it to 1D; T7 fixes noise units for any $\Sigma_i$ (O8). The strong form of L8 was never needed. See the O3 record §§2, 4b. |
| T7 | Prior once | nothing attached to a reading may depend on $m$. The push weight's determinant is $m$-free for every $\omega$ (ours28 §3), which forces noise units. | L11 | strengthened |
| T8 | The law (**ours29**) | $p(m)\propto\pi(m)\exp\!\big(\eta\sum_i g_i H_i(m)/2\big)$: each reading's flat factor $e^{H_i/2}$ ($H_i$ the explained energy, whitened data projected on the whitened candidate subspace), weighted by its own misfit $g_i=g(Q_i)$ (L21, T30), the product raised to the Bartlett power $\eta$ (T31). ours27 is $g_i=\eta=1$. ours28's ${}_1F_1$ factors (the arcsine size law) are withdrawn with L10′. | L1, L13, T6, T7, L10, L21, T31 | **Changed** from ours28, 30 Sept (ruling 2) |
| T9 | Stacked size | $\lVert\mathbf a^\ast\rVert=\big(\sum_i\lVert\mathbf a^\ast_i\rVert^2\big)^{1/2}$ per candidate, before averaging | L12, L18 | Unchanged. **Supported again, 30 Sept:** within a reading a vector's length is the Euclidean join of its components; L12 (one joining rule) makes the join across readings the same one, with components as the atoms (L18). The plain sum would need a second rule. Rests on L12, graded medium: if L12 is doubted, T9 is chosen and ours24's plain sum returns as a rival. See the O3 record §4. |
| T10 | The readout | $\widehat m=\int mA\,p\,dm\big/\!\int A\,p\,dm$ with $A(m)=E[\lVert\mathbf a^\ast\rVert\mid m]$ over the candidate cloud | T4, T8, T9 | **Fails L20 (T29). Open (O1).** |
| T11 | Sufficiency | data enter through $(P_g,Q_g,D_g,M_g)$ per noise pair | T8, T9, T6 | **Returns under ours29** for the law's shape, weighted by $g_i$; $\eta$ needs in addition the readings' second moments ($\sum_i\lvert C_i\rvert^2$, $\sum_iC_i^2$ in 1D). Under ours28 the state was each reading's $H_i(\cdot)$. |
| T12 | Combination | exact, as a product of per-reading factors | T8 | Under ours29 the law's state adds (the Gram matrices, weighted by $g_i$) and $\eta$'s does too; the tempered law of a union is not the product of the tempered laws of its parts, since $\eta$ is recomputed on the union. **Additive state is a corollary of T6, and characterises it:** the law depends on the readings only through $\sum_iM_i$ iff the kernel is exponential (O3 record §3; the ours21 kernel differs by up to 0.52 nats between two series with the same Gram sum). |
| T13 | Totality | proper $\pi$ $\Rightarrow$ $\widehat m\in(0,\infty)$; the totality checks pass under ours28 and ours29 (perpendicular, anti-parallel, force only, acceleration only; the all-zero experiment gives median $c$) | T8, T10, L2 | Absorbs L2's old clause (E5). **Breaks at the all-zero experiment under the O1 repair.** |
| T14 | Zero readings give $s$ | readout at zero readings is $s$ | L2, L3′, L4, T5 | Scoped to this readout, [one pair], and $c=s$; rests on L3′, which is unruled (sweep E1: at $c=3s$ the readout is 2.095). **Contested by the O1 repair.** |
| T15 | Equal energies | $P=Q\Rightarrow\widehat m=\text{median}=s$ | T11, T5, L4 | [one pair] [flat]; unchecked under ours28 |
| T16 | A tilt is an unseen reading | $\lambda\sin2\theta_c$ at any centre $c$ equals one reading with $P=Q$, $D=2\lambda$, at the level of the law | T8 | Broadened (sweep E3); **[flat]**: under ours28 a reading's evidence is ${}_1F_1$, not exponential; unchecked |
| T17 | Convergence | $\widehat m\to m$ as $R^{\ast2}\to\infty$ | T8, T7 | Unchanged; ours28 checked by Monte Carlo (ours28 record §5) |
| T18 | Width | $\mathrm{Var}_{\rm true}/\mathrm{Var}_{\rm law}=1+\sum_c\theta_c'^2\big/\sum_c\theta_c'^2r^{\ast2}_c$, $\theta_c'=s_c/(s_c^2+m^2)$ | T8 [flat], L13 | The flat factor's width error, as a function of the true signals. Under ours29 it is **measured per series** by $1/\eta$ (T31), from the readings alone; $\eta$ equals its inverse in the long-series limit. Per-component form (E6). |
| T19 | No mean | $E_p[m]=\infty$ under the half-Cauchy | T8 | unchanged |
| T20 | Mixed-instrument no-go | $\pi$ is declared before data, as part of the specification | T8, T14 | Unchanged; the size law that shared this slot is withdrawn (30 Sept). |
| T21 | **C3** | the readout depends only on the law of the two joined sizes $(X,Y)$: no direction weight, no weight on the per-reading size profile | L15, L16, L17, L18, composition | **new, derived 28 Sept.** The betweenness representation step (Dekel, borrowed) is owed in writing (O5). |
| T22 | Full covariance | the error model must be a full $\Sigma_i$ per reading; ours26's diagonal model gives a different answer in a rotated frame | L16 | new (ours27) |
| T23 | Gaussian noise forced | by consistency given the flat push; by Z under ours28 | L11, L20 | new; the consistency route under tempered evidence is owed (O4) |
| T24 | Z holds for the law | an exact zero reading has $H_i\equiv0$, so $E_i$ is constant for any $\Sigma_i$, any $k_i$ | T8 | new; exact |
| T25 | Zero components carry nothing; $k$ does | measured exact zeros on some axes of a signal reading leave $H_i$ exactly as without those axes, at every $m$, whatever their $\sigma$. The law changes only through $k_i$. A 3D reading declared $k=1$ equals the 1D reading exactly. | T8 | **new, 29 Sept** |
| T26 | Limits of per-reading rules | no per-reading size rule has the right width at every push size; none is never overconfident (a pure-noise reading adds score noise and no curvature) | L19, T8 | new (no-borrowing record). **Why $\eta$ is series-level** (T31): the width cannot be fixed reading by reading. |
| T27 | Own-world calibration | a proper law's intervals are calibrated on average in its own world; ours27 (improper) is not. The planted-mass yardstick condemns every proper starting law at extreme masses. | T8 (ours28) | new (betweenness record §2, full-covariance record §5). **Under ours29 there is no declared world of signals to apply it to** (L10); the standard is ruling 3 of 30 Sept: the quoted width ($\eta$ applied) must cover, and the raw law is reported beside it. |
| T28 | Mass and push scale are orthogonal | Fisher-orthogonal under ours28's model; the reference prior on mass is uniform in $\theta$ (the half-Cauchy) for any size law | T8 (ours28) | new (Jeffreys-check record §1); with no size law in ours29 it says only that $\pi$'s reference form does not depend on one |
| T29 | **The readout fails Z** | any readout that averages a push size over each reading's candidate cloud changes when a null reading is appended: the cloud has non-zero size, and that size depends on $m$ through the noise units. ours28, 3D, +50 null readings: 2.2524 → 2.2406; law unchanged to $10^{-14}$. Using the size of each reading's **best-guess (conditional-mean) candidate** passes Z exactly, passes swap and the exact limit, and leaves the all-zero experiment with no readout. | T10, L20 | **new, 29 Sept.** Uniqueness of the repair not shown. |
| T30 | **The tempered law's two ends** | $g(0)=1$ by L20 (an exact zero reading has $Q_i=0$; T24, T25 survive). Bounded confidence on a fixed misaligned pair as $\sigma\to0$, with the evidence not going to zero, forces $g\sim\kappa/Q$; the limit tilt is $\kappa\,a_i/(1-a_i)$, $a_i=|C_i|/T_i$. $\kappa=\nu=k_i-1$ reads as the share of the reading's self-contradiction the declared noise accounts for (forced by meaning, not by a property). **The same meaning fixes the middle (30 Sept, O12 closed):** read at every $Q$ and capped at 1 because a share cannot exceed the whole, $g=\min(1,\nu/Q)$; the smooth, quantile-threshold and shifted forms keep the meaning only in the tail, and the parametric ones add an unsupplied parameter. It is $\eta$'s rule (expected over observed disagreement, capped at 1) at reading level, and the PDG scale-factor convention. Cost on ordinary data, measured: about 7% width at strong signal in 2D/3D (efficiency 0.92–0.93), coverage 94–96%, accuracy unchanged; recoverable only by lifting $\eta$'s cap (O14). Scale-free discounts excluded; a heavy-tailed error model excluded by L20 (T23 forces the along-line shape only). Inert in 1D and for $k_i=1$. Tested 2D: rescue 8.5% $\to$ 91–100% coverage; own-world cost 0–6% width for forms flat inside the noise band. | L21, L20 | **new, 29 Sept.** Does not touch the $\sqrt N$ tilt. |
| T31 | **The Bartlett power** | $\eta=\min(1,\mathcal H/J)$ with $\mathcal H=-\sum_i\ell_i''$, $J=\sum_i\ell_i'^2$ at the series' best line $\hat\phi$ on the **whole projective line** of slopes (negative and infinite $m$ included), $\ell_i=g_iH_i/2$; $\eta=1$ when $J=0$. Forced: $\mathcal H/J$ by property N in the long-series limit (the width must track how far the point moves under repeats); one $\eta$ per series because $m$ is one number; **the raw $J$** (no $N/(N-1)$) by Z, since an exact zero reading has $\ell_i'=\ell_i''=0$ and a Bessel factor counts it (ours25's sandwich moved on appended zeros); the whole line because on $m>0$ the best line can sit at a boundary where the curvature measures nothing (ours25 then fell back to the prior); evaluation at $\hat\phi$ because the ratio is parametrisation-free only where the score vanishes (the swap follows); candidates unchanged by L19 (shrinking each cloud by $\eta$ would infer signal size from the series). Natural, not forced: the cap at 1; the power form at finite $N$ (O14). Removes the $\sqrt N$ growth of the null tilt and leaves a bounded floor (O15). Tested 1D–3D: quoted width 93.7–97% wherever there is signal at $N=200$; the point stops wandering on pure noise. | L19 (scope), L20, L1, property N | **new, 30 Sept.** The PDG's Birge-ratio inflation is the same move on Gaussian points. For the flat base the tempered law coincides with a shared-size model, $\omega^2=\eta/(1-\eta)$: same function, different premise; Z separates them. |
| T32 | **The misfit is the pairings' disagreement in noise units** | for $\Sigma_i=\mathrm{diag}(\sigma_F^2,\sigma_a^2)\otimes I_k$ and $N_i=I$: with $M_i$ the reading's whitened Gram matrix (30 Sept dictionary), $Q_i=\lambda_-(M_i)$, and the product of the two pairings is $M_i$'s squared correlation, $P=M_{12}^2/(M_{11}M_{22})$. So $Q_i\,\lambda_+=M_{11}M_{22}(1-P)$ and $Q_i=\tfrac12\big[T_i-\sqrt{T_i^2-4M_{11}M_{22}(1-P)}\big]$: the dimensionless disagreement $1-P$ converted into noise units of energy by the sizes. The physical and whitened $P$ coincide when the channel noise is diagonal; with cross-channel correlation $\rho$ the law sees the whitened one ($P=0.62$ physically, $0.71$ at $\rho=0.5$, $0.97$ at $\rho=-0.7$). $g$ is not a function of $P$ alone and must not be: the same angle at scales $0.1,1,10$ gives $Q=0.002,0.21,21.1$ and $g=1,1,0.09$ (a small misaligned pair is consistent with alignment under noise). At equal whitened energies $Q=s^2(1-\lvert\cos\theta_w\rvert)$ exactly. | T8, L21, T30, the 30 Sept eigenvalue dictionary | **new, 10 Oct.** Checked against `noise_vs_signal/ours29_reference.py` (`reading_state`) on 300 random readings, $k=1,2,3$: agreement to $7\times10^{-13}$. |
| T33 | **Same window** | a reading's force and acceleration are averages over one window: the time pairing $(\Delta t)'\Delta t=1$ in $m=\mathbf F\cdot\mathbf a'=\Delta\mathbf p\cdot(\Delta\mathbf v)'$ needs the same $\Delta t$ in both. For a non-constant push, mismatched windows break Newton II with no noise at all: $F(t)=5(1+\sin2t)$ in a fixed direction, $m=2$, gives $\mathbf a'(\mathbf F)=2.000$ on one window, $1.935$ with the acceleration window shifted by half, $2.417$ with it doubled. Free for constant signals. Two measurements over the same window share one signal by construction, so under L19's reason they cannot be two readings: one reading, combined $\Sigma$. Measurements over different windows of an unchanged push are separate readings (O9). | the units build-up (units-and-axes record §6), L18, L19 | **new, 10 Oct**; the window is declared (§1, L18). |

## 4. What is declared (the free shapes)

| input | current declaration | forced parts | free part |
|---|---|---|---|
| $\pi$ | half-Cauchy, $\pi(m)=c/(c^2+m^2)$ | symmetric (L4), counted once (T7), declared before data (T20), unit-covariant (L3) | shape: any $m^{-1}g(\ln(m/c))$ with $g$ even; centre with several noise pairs (L3′ unruled) |
| $N_i$ | $I_d$ unless the push's direction is known | — | the declaration ladder: $k=d$ → $k<d$ → a declared push covariance → old readings in the prior slot |

**One rule, one use** (was two). "Uniform in the angle the noise defines" gives
the half-Cauchy on $m/s$. Its second use, the arcsine on $\beta$, went with the
size law on 30 Sept, and $g$'s middle left this table the same day (T30, O12
closed). **$\pi$ is the one declared shape left in the definition.**

**F2 asks a question the records have not answered.** The full-covariance
record §5 argued that a declared shape is "the experiment's description", not
a rival definition of mass. Under F2, that is the "nobody cares" move unless
it is argued. It needs your ruling (O6).

## 5. Open, in the order the definition needs them

| | item | what would close it |
|---|---|---|
| **O1** | **The readout under Z** (T29). Three things cannot all stand as now written: Z at the point (L20); a readout at the all-zero experiment (L2 as d2 reads it); a readout averaging sizes over the candidate cloud (T10). The best-guess repair keeps Z and drops the second. **Untouched by ours29** (ten zero readings still move the readout, 3.5983 → 3.6017, while $\eta$ and the law do not). | your ruling on the all-zero experiment; then O2 |
| O2 | **Uniqueness of the repaired readout.** Does composition (L5, L14) still single out a ratio of means when sizes are best guesses? Does any other readout pass Z? Can any Z-passing readout give a value at the all-zero experiment? | proof, or a counterexample family |
| ~~O3~~ | **Closed 30 Sept** (his "Yes I like it. I agree."). The flat-in-vector weighting (T6), ours29's whole per-reading factor, follows from L8's conditional form and L9 by Herschel–Maxwell, with L10 setting the scale to infinity and L7 removing the improper members; additive state (T12) is then a theorem, and the exponential is the only kernel with it. The stacked length (T9) follows from L12 with components as atoms (L18). The 25 Sept axis-splitting premise is replaced, not needed. **Dependencies kept in view:** L8's conditional form (ruled, inherited, not itself forced) and L12 (medium). | closed; `noise_vs_signal/o3_what_supports_the_flat_weighting_and_the_stacked_length_2026-09-30.md` |
| O4 | L13 under the tempered evidence (the consistency route). | quadrature, as in matrix §4 |
| O5 | The betweenness representation step, written out (T21's borrowed step). | proof on candidate space |
| O6 | **$\pi$'s shape** (was: the two shapes; $G$ is gone). Forced by some property of mass, or a parameter of an alternative Newtonian mechanics (F2)? | a forcing argument, or a ruling that records it as an alternative |
| O7 | **Z's scope.** Your proposal of 29 Sept: two kinds of zero, a whole null reading and zeros alongside a signal. T25 says zero components contribute nothing and what matters is the declared $k$. Is "embedding in more dimensions" information (then Z is about whole readings), or should Z cover it (then ours28's $k$-dependence is a defect)? | ruling |
| ~~O8~~ | **Closed 30 Sept** with O3. Physically isotropic pushes under an anisotropic instrument: the fork between noise-unit and physical isotropy is real only for finite-width weights. Under L10 the weight is flat, a flat weight changes under a linear map only by its Jacobian, and flat-in-physical carries $\det G_i(m)^{-1/2}$, a factor on mass per reading that T7 forbids; so noise units, for any $\Sigma_i$. Checked by quadrature under a strongly anisotropic, correlated $\Sigma$. | closed; `noise_vs_signal/o3_what_supports_the_flat_weighting_and_the_stacked_length_2026-09-30.md` §4b |
| O9 | Repeated measurements of one unchanged push: are they one reading (one scale) or several? L19 treats readings as separate pushes. **10 Oct, proposed:** the criterion is the window (T33). Same window: one reading, one signal, combined $\Sigma$. Different windows: separate readings, L19 as written. That the push is "unchanged" does not enter: it would be a statement about signals across readings. | ruling on the window criterion |
| ~~O10~~ | **Closed 30 Sept, 16:26, by ruling:** persistence is a premise of L1 ("a perfectly fine lemma. Doesn't need proving"). | closed |
| O11 | L3′, and with it T14. **10 Oct, the units argument, offered.** Premise: the default declaration for mass with instruments $(\sigma_F,\sigma_a)$ and for inverse mass with $(\sigma_a,\sigma_F)$ are one declaration read backwards (the symmetry of m and $\mathrm m^{-1}$, units-and-axes record §3). Then $c(\sigma_F,\sigma_a)=1/c(\sigma_a,\sigma_F)$, and with $c=k\,\sigma_F/\sigma_a$ from L3, $k^2=1$. In general: **$c$ may depend only on declared quantities that invert under the swap, and any swap-invariant dimensionless factor is 1.** So $c$ cannot depend on $\rho$, $k$, $d$ or $N$, and one noise pair, correlated or not, forces $c=s$ (T14's value). Where it stops: several instruments leave every weighted geometric mean $\prod_g s_g^{w_g}$ (T20 stands); anisotropic $\Sigma_F,\Sigma_a$ leave any $h$ of the eigenvalues of $\Sigma_a^{-1/2}\Sigma_F\Sigma_a^{-1/2}$ with $h(1/\lambda)=1/h(\lambda)$, of which the determinant ratio $(\det\Sigma_F/\det\Sigma_a)^{1/2d}$ is the unique power mean. | ruling on the swap-of-declarations premise and on L3′'s first half; then the several-instrument and anisotropic cases |
| ~~O12~~ | **Closed 30 Sept, 16:26** ("O12 is good also"). The middle of $g$ is forced by the meaning that fixed $\kappa$ (T30). The "evidence" line of the afternoon (threshold too low) is retracted in the O12 record: it measured against ours27, not against coverage. | closed; `noise_vs_signal/o12_the_middle_of_the_misfit_weight_2026-09-30.md` |
| O13 | **Contradiction between readings.** Two readings declared separate that contradict each other about $m$ collapse the law onto the stronger with unbounded confidence; the law cannot see the contradiction, the data-only resolution (property N) can. Same family as the $\sqrt N$ tilt. **Partly answered by T31:** contradiction about $m$ is exactly what $\eta$ reads, and it widens the law; the resolution still reports it. | ruling on whether anything more is owed |
| O14 | **$\eta$ at finite $N$** (T31). The power form is asymptotically forced; with few readings the raw $J$ is biased low (weak signal at $N=20$ covers 88–89%). A Z-safe small-sample correction would count effective readings, $(\sum_ic_i)^2/\sum_ic_i^2$ with $c_i$ the curvature contributions, never $N$. With several instruments, one $\eta_g$ per noise pair may be what property N asks for (the ours28b calibration theorem: sharing is exact only within classes of equal sensitivity). | derivation, then the tester |
| O15 | **The null floor.** $\eta$ stops the null tilt growing but leaves it bounded (about 2 nats; the quoted width covers 75–85% against a mass drawn from $\pi$). Under Z and L19 null evidence does not accumulate (a long null series looks no more null than a short one), so this may be the most any law can do, with "we know nothing" living in the resolution. | a proof, or a counterexample law |

## 6. Standards for tests (from the rulings, for the tester)

- **Tier 1, properties.** Every L and T above that is stated as exact gets a
  pass/fail check: swap, units, rotation with cross-correlation, Z (law *and*
  point), axis splitting with state carried, $k=1$ equals 1D, totality,
  lumping, common force, betweenness.
- **Tier 2, own world.** Data drawn from the definition's own declared
  measures; intervals must be calibrated (T27). A failure here is a bug in
  the code or in a derivation.
- **Tier 3, other worlds.** Planted masses, fixed-length pushes, undeclared
  direction structure, shared push size. These say what the definition *does*
  when a declaration is wrong, not whether the definition is right. Every
  tier-3 failure so far has named an undeclared input.

## 7. Sources

- `mass_estimator_lemmas_2026-09-29.md` (superseded by this record; it superseded the 27 Sept record)
- `mass_as_force_paired_with_reciprocal_acceleration_2026-10-10.md` §§6–9 (the window and T33; L5, L14 as pairings; T1 in pairing form; T32; O9, O11; 10 Oct). Its §9.4 holds the checks, run on his laptop against `noise_vs_signal/ours29_reference.py`; the check script is not yet kept (item 3 of the 10 Oct recording plan, pending).
- `ours26/ours26_sweep_of_the_flat_records_2026-09-27.md` §§4–5
- `ours27/ours27_full_covariance_2026-09-27.md`
- `ours27/ours27_perturbation_matrix_2026-09-27.md` §§3, 6
- `ours27/ours28_candidate_push_scale_as_nuisance_2026-09-28.md` §§1–4, 8
- `ours27/ours28_size_measure_jeffreys_check_2026-09-28.md` §§1, 4
- `ours27/ours27_c3_from_betweenness_and_the_scale_fork_2026-09-28.md` §§1–4
- `ours27/ours28_full_covariance_form_and_declared_subspace_2026-09-28.md` §§1–6
- `ours28b/ours28b_no_borrowing_between_readings_2026-09-28.md` §§1, 5
- `ours27/ours28_reference_implementation_and_readout_2026-09-28.md`
- `ours27/ours28_c3_narrow_form_and_the_reading_structure_2026-09-28.md` §§4, 9
- `../_note_mass_estimator_equation/the_mass_estimator_equation_d2.md` §3.1 (P1–P8, step 1)
- `noise_vs_signal/readings_are_evidence_just_less_2026-09-29.md` §§2, 10 (L8, L21, T30, O3, O12, O13; 29 Sept, night).
- `noise_vs_signal/o3_what_supports_the_flat_weighting_and_the_stacked_length_2026-09-30.md` §§2–4b, 6 (L8, T6, T9, T11, T12, O3, O8; 30 Sept, late afternoon).
- `noise_vs_signal/o12_the_middle_of_the_misfit_weight_2026-09-30.md` §§2–4 (T30, O12, §4; 30 Sept, 16:26).
- `noise_vs_signal/ours29_misfit_weight_and_bartlett_power_2026-09-30.md` §§1–3, 5, 7–11 (T8, L10, L10′, L4′, L19, T31, T18, T26, T27, §4, O1, O3, O6, O12, O13, O14, O15; 30 Sept, afternoon), and `noise_vs_signal/scattered_observations_2026-09-30.md` §§2, 4, 10 behind it.
- 29 Sept, this session: T25 and T29, computed with
  `ours27/ours28_reference.py` (numpy, scipy). T25: $H_i$ on a 7-point mass
  grid for the 3D reading $F=(4,0,0)$, $a=(2,0,0)$ against the 1D reading
  $(4,2)$, unit noise (difference 0.0); then the law and readout with the
  transverse SDs at 10, 1, 0.1, 0.01 and 0.001 (law identical in every case).
  T29: 1, 5 and 50 exact zero readings appended to one 1D and one 3D reading
  ($c=2$, unit noise), ours27 and ours28; and 20 appended to six weak 3D
  readings with the best-guess size computed by sampling the per-reading
  scale weights (3000 draws per mass node). Scripts not kept.
