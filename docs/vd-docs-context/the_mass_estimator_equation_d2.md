# The mass estimator equation: what it is, and why each part is there

**Note:** `_note_mass_estimator_equation`, draft d2, 2026-09-26. Written by
Claude at the user's request, from the conversations of 25 and 26 September
and the records they produced.

**What this is.** The equation we feature, **ours25b**, stated once with its
types, then justified as a chain of moves from Newton II. Each move is tagged:

- *assumed*: a starting point the user has accepted;
- *proved*: follows from what comes before it;
- *ruled*: the user's decision;
- *open*: not yet settled.

The method is the user's: look for the moves that are forced, "put your hand on
the left wall and walk". A choice made because it scored well in a test is not
a move in this chain. Where the chain has not reached the wall yet, it says so.

**What this is not.**

- Not the argument that an equation is needed at all. That is the other note,
  `_note_why_mass_estimation`.
- Not the literature or credit account. That is
  `../_dscn_mass_estimator/credit_and_originality_of_the_mass_estimator_2026-09-26.md`.
- Not the story in date order. That is d1.

**Relation to d1.** d1 (23 September) told the story up to the 21 September
equation, here called ours21. Since then the law changed (24 September) and so
did the many-reading readout (25 September). The following d1 claims are **not**
carried into this draft (from
`what_ours24_and_ours25b_change_in_earlier_results_2026-09-26.md` §7):

- d1 §5's weighting $df\,d\alpha\,d\Omega$. It is ours21's; the current one is
  flat in the push vector (step 7).
- d1 §6's "the 1D law is Liseo's eq. (12)". For the current law it is Liseo's
  eq. (14), the reference-prior posterior.
- d1 §7's warrant for the weighting (lab contribution 6). It warrants ours21's
  weighting, not the current one.
- d1 §7's verdict "the weighting is not forced", its §9 row "declared choice
  with a partial warrant", and its §10 open item on the three requirements
  behind the flat weighting. The push part is now forced, given axis splitting
  (step 7). Only the shape of the starting belief stays open (step 8).
- d1 §10's "several readings of one mass need a declared joint model", listed
  as open. The current law and readout supply it.

d1 §8 (no always-finite interval can guarantee coverage) stands unchanged.

**Direction.** Not yet ruled by the user.

---

## 0. Where this note starts

`_note_why_mass_estimation` argues three things, taken as given here:

1. $m=\mathbf F/\mathbf a$ is not an operation. One vector cannot be divided by
   another.
2. Force and acceleration must each be obtainable without the other, or
   $\mathbf F=m\mathbf a$ and $\mathbf a=\mathbf F/m$ could not be used to
   predict. So independently measured pairs do not line up.
3. Lining them up is therefore part of finding the mass, not a step before it.

This note builds the operation that does the lining up and returns the mass.

## 1. Types

**Inputs.**

- A dimension $d\in\{1,2,3\}$. A cart on a track is $d=1$.
- $N\ge1$ readings of one object, all sharing one mass. Reading $i$ is a
  measured net force $\tilde{\mathbf F}_i\in\mathbb R^d$ (in N) and a measured
  acceleration $\tilde{\mathbf a}_i\in\mathbb R^d$ (in m/s²), in the same axes.
  Any values are allowed, including zero.
- Uncertainties $\sigma_F>0$ and $\sigma_a>0$: one per channel for the whole
  series, the same on every axis and in every reading. Each is a per-component
  standard uncertainty.
- **The declared noise model:** the errors are zero-mean Gaussian, independent
  across components, channels and readings, with those standard uncertainties.
  This is a declared input, as it was for ours21. The law, the readout and the
  arguments of steps 9 and 10 all use it.
- **Nothing else goes in.** The function has no built-in scale of its own, so
  changing the units of force or acceleration changes the answer only as a mass
  changes.

**The "$>0$" is part of the type.** We cannot promise that the acceleration is
nonzero; that is the problem. We can promise that each instrument has noise.
That promise is what keeps the function from ever returning "undefined"
(steps 1, 2 and 4).

**Output.** A point $\widehat m\in(0,\infty)$ and a law over $(0,\infty)$,
from which the median, the intervals and the spread of $\log m$ come. The
readings and uncertainties travel with the result. The output is never
"undefined".

**Not this function.** Exact inputs ($\sigma=0$) go to a separate check with
three outcomes: *one mass* (nonzero vectors pointing the same way), *no mass
fits* (exactly one zero, or different directions), *every mass fits* (both
zero). "Undefined" survives only there, for exact zero–zero.

**Working quantities.**

- Noise units: $\mathbf x_i=\tilde{\mathbf F}_i/\sigma_F$ and
  $\mathbf y_i=\tilde{\mathbf a}_i/\sigma_a$. Each channel's noise is then one
  unit per component.
- The instrument ratio $s=\sigma_F/\sigma_a$. It is a mass.
- A **candidate** is a possible true pair behind a reading. In noise units it is
  a push vector $\mathbf v_i\in\mathbb R^d$, shared by both channels, and an
  angle $\theta\in(0,\pi/2)$, shared by every reading:

$$
\frac{\mathbf F^\ast_i}{\sigma_F}=\mathbf v_i\sin\theta,\qquad
\frac{\mathbf a^\ast_i}{\sigma_a}=\mathbf v_i\cos\theta,\qquad
m=s\tan\theta .
$$

The length $r_i=|\mathbf v_i|$ is how hard reading $i$ was pushed, measured in
both channels together. The angle carries the mass.

## 2. The equation (ours25b)

**Names.**

| name | what it is |
|---|---|
| ours21 | the 21 September equation (`../_dscn_mass_estimator/mass_estimator_equation/`) |
| ours24 | the 24 September law, with the old many-reading readout |
| **ours25b** | ours24's law with the stacked-length readout (25 September) |

For one reading, ours24 and ours25b are the same equation.

**The three sums.**

$$
P=\sum_{i=1}^{N}|\mathbf x_i|^2,\qquad
Q=\sum_{i=1}^{N}|\mathbf y_i|^2,\qquad
D=\sum_{i=1}^{N}\mathbf x_i\cdot\mathbf y_i,\qquad
M=Nd .
$$

**The law of mass.**

$$
p(m)\;\propto\;\frac{s}{s^2+m^2}\;
\exp\!\left[-\,\frac{\sum_i\big|\tilde{\mathbf F}_i-m\,\tilde{\mathbf a}_i\big|^2}{2\left(\sigma_F^2+m^2\sigma_a^2\right)}\right].
$$

- The first factor is the starting belief, counted once (step 8).
- The exponent adds one misfit per reading: how far the reading is from
  satisfying $\mathbf F=m\mathbf a$, in units of the noise that such a misfit
  would carry.

In the angle, the same law is

$$
p(\theta)\;\propto\;\exp\!\left[-\tfrac12\left(P\cos^2\theta+Q\sin^2\theta-2D\sin\theta\cos\theta\right)\right],\qquad \theta\in\left(0,\tfrac{\pi}{2}\right).
$$

**The readout.** Stack the $N$ readings' true force vectors into one vector
$\mathbf F^\ast\in\mathbb R^{M}$, and the true accelerations into
$\mathbf a^\ast\in\mathbb R^{M}$. Then

$$
\boxed{\;
\widehat m\;=\;\frac{E\big[\,|\mathbf F^\ast|\,\big]}{E\big[\,|\mathbf a^\ast|\,\big]}
\;=\;s\;\frac{\displaystyle\int_0^{\pi/2}\sin\theta\;\mu_M\big(H(\theta)\big)\,p(\theta)\,d\theta}
{\displaystyle\int_0^{\pi/2}\cos\theta\;\mu_M\big(H(\theta)\big)\,p(\theta)\,d\theta}
\;}
$$

with

$$
H(\theta)^2=P\sin^2\theta+Q\cos^2\theta+2D\sin\theta\cos\theta,
\qquad
\mu_M(h)=\sqrt2\,\frac{\Gamma\!\left(\frac{M+1}{2}\right)}{\Gamma\!\left(\frac{M}{2}\right)}\;
{}_1F_1\!\left(-\tfrac12;\tfrac{M}{2};-\tfrac{h^2}{2}\right).
$$

- $\mu_M(h)$ is the mean length of an $M$-dimensional Gaussian vector with unit
  noise whose centre is at distance $h$ from the origin. At a fixed angle the
  stacked true push is such a vector, centred at
  $\mathbf X\sin\theta+\mathbf Y\cos\theta$. Here $\mathbf X$ and $\mathbf Y$
  are all the readings' $\mathbf x_i$ and $\mathbf y_i$ stacked into two vectors
  in $\mathbb R^{M}$, so $|\mathbf X\sin\theta+\mathbf Y\cos\theta|=H(\theta)$.
- For $M=1,2,3$, $\mu_M$ is the folded-normal, Rice and 3D closed forms.

**Plain reading.** Weight every candidate true pair by how well it explains all
the readings. Average the size of the true force, and the size of the true
acceleration, over those candidates. Divide.

**What the output reports.**

- The point $\widehat m$ above.
- The median and central intervals of $m$ under $p$.
- The spread of $\log m$.

The ordinary mean of $m$ under $p$ is infinite, so "mass ± SD" is not
available.

## 3. The chain

### 3.1 Premises carried in

| # | premise | status and source |
|---|---|---|
| P1 | **Newton II.** $\mathbf F^\ast=m\,\mathbf a^\ast$ for one object, one interval, one inertial frame, with $m>0$. | *assumed* |
| P2 | **The input type** of §1: $\sigma_F,\sigma_a>0$, the declared Gaussian noise model, and no built-in scale. | *assumed*; $\sigma>0$ as part of the type is the user's, 25 Sept; the Gaussian model is carried from ours21 |
| P3 | **Newton II is assumed true here, not tested.** Any mismatch between the readings is measurement error. The old rejection gate is not part of the equation; the misfit is kept with the result. | *ruled*, 7 Sept (`mass_estimator_law_assumed_correct_steering_2026-09-07.md`) |
| P4 | **Force and acceleration count as equal information.** Swapping their roles, with $m$ read as $1/m$, gives the reciprocal answer. | *assumed*; the user, 7 and 24 Sept |
| P5 | **Lumping.** Glue two bodies together so they share one acceleration, and their masses add. | *assumed*; a commitment of the dictionary's lumping material |
| P6 | **Axis splitting.** A reading, and the same reading taken as one-dimensional readings (one per axis), must give the same law of mass and the same answer. Newton II is $d$ component equations sharing one mass; whether a reading is recorded as one vector or as $d$ numbers is bookkeeping. Step 11 uses it in its **strong form**: one candidate has one size, however its axes are grouped. | *assumed*; proposed 25 Sept, judged by the user "too elementary to be wrong"; the strong form is how the ours25b record reads it |
| P7 | **No preferred direction in space, and no preferred push size before data.** | *assumed*; the first is already in ours21 (uniform direction), the second from the 24–25 Sept records |
| P8 | **One joining rule.** When two sets of readings of the same object are joined, the size of the whole comes from the sizes of the parts by one fixed rule. The old plain sum assumed this too. | *assumed*; 25 Sept |

### 3.2 The steps

**Step 1. Measured and exact inputs are separate cases.** A measured mismatch
still gives a mass: force read with no acceleration, or the two vectors pointing
different ways. Only exact zero–zero may end at "undefined", and that is on the
exact branch, outside this function.

- **The disagreement rule's distinction.**
  - *No mass fits* is a disagreement with Newton II. On measured data it is
    repaired as measurement error (P3), not ended at "undefined".
  - *Every mass fits* is agreement without identification. For exact inputs it
    ends at "undefined", and no mass may be manufactured there.
- **Measured zero–zero is not settled by that rule.** The rule itself (its
  §10) leaves a mass at measured zero–zero as a separate, stronger proposal.
  What settles it here is P2: with $\sigma>0$ a zero reading is not an exact
  zero, and the function's type promises a mass for every measured input.
  Steps 2 and 4 show that the promise can be kept.
- **Status.**
  - The exact/measured split and the repair of mismatches: *proved*
    (`disagreement_rule_mass_estimator_trace.md`, 29 Aug; the 2 Sept chain
    audit §6).
  - A mass at measured zero–zero: *assumed* (P2, the user's, 25 Sept; it is
    also spec d1 §21's position).

**Step 2. At zero readings, any estimator of this type returns
$s=\sigma_F/\sigma_a$.**

- With every reading zero, the only mass-sized quantity the function receives
  is $s$, since it has no built-in scale (P2). So a function that always returns
  a mass must return $c\,s$, with $c$ a pure number (it may depend on $N$ and
  $d$).
- Swap the channels (P4). The same function then returns $c/s$ for the inverse
  mass, and P4 demands that this equal $1/(c\,s)$. So $c^2=1$, and $c=1$.
- This holds for every construction, not only this one. The instruments appear
  in the answer at zero readings because of dimensions, not by a choice we made.
  The answer describes the apparatus, not the object: an object that is neither
  pushed nor accelerating tells us nothing about its inertia.
- *proved*; 25 Sept conversation.

**Step 3. Candidates.** Every true pair that points the same way, with positive
sizes: $\mathbf F^\ast=f\,\mathbf u$ and $\mathbf a^\ast=\alpha\,\mathbf u$,
with $f,\alpha>0$ and a common unit direction $\mathbf u$. In noise units this
is the push vector and angle of §1.

- *proved*; from P1.

**Step 4. Average over the candidates; do not pick the best one.**

- **Picking the best fit is not total, even with $\sigma>0$.** With force read
  as nonzero and acceleration as zero, the misfit
  $|\tilde{\mathbf F}|^2/(\sigma_F^2+m^2\sigma_a^2)$ keeps falling as $m$ grows,
  so the "best" mass is infinitely heavy, which is not a mass. At zero–zero
  readings every mass fits perfectly, so there is no single best. With
  anti-aligned readings the best positive mass runs to $0$ or to infinitely
  heavy.
- **Averaging is total.** With $\sigma>0$ the averages of $f$ and $\alpha$ over
  the weighted candidates are always finite and positive, so their ratio is
  always a mass.
- *proved*; spec d1 §7, the 27 Aug median record §2, the 2 Sept chain audit §5,
  with P2.

**Step 5. The readout divides the average force size by the average
acceleration size.**

- It is forced by any two of three requirements:
  - **lumping** (P5);
  - **common force:** one force on two bodies, so the accelerations add, so
    the inverse masses add;
  - **the swap symmetry** (P4).
- The third follows from the other two:
  - **Lumping and the swap give common force**, in one line:
    $1/R(X,Y_1+Y_2)=R(Y_1+Y_2,X)=R(Y_1,X)+R(Y_2,X)$. This uses the swap at the
    level of the readout: $R(Y,X)=1/R(X,Y)$ for every law in the theorem's
    domain, which reaches as far as the composition laws themselves do.
  - **Lumping and common force** force the readout (the 6 Sept theorem, which
    uses exactly these two), and that readout is symmetric.
- This note leads with lumping and the swap. Lumping is the most physical, and
  the swap is already a premise.
- Neither composition law alone is enough; each lets a whole family of readouts
  through.
- **The theorem's other conditions:**
  - it is exact on sure inputs;
  - it is continuous at constants;
  - its domain is closed, covering every joint law of positive sizes with
    finite means;
  - the readout depends only on the joint law of the two sizes (C3).
  
  The first two are mild. The last two are real premises. C3 is the one a
  counterexample in the 6 Sept record turns on (a readout that also uses the
  direction). The domain is wider than the laws this equation produces.
  Both are listed in §6.
- **Why it composes:** on measured data the estimator is what "mass" means, so
  the number it reports must add as masses add. That is a coherence requirement
  between dictionary entries, not an optimality claim.
- **An identity:** $\widehat m=E[\alpha\,m]/E[\alpha]$. It is an average of the
  mass, weighted by how hard each candidate accelerates. The weight cancels the
  $1/\alpha$ that makes the plain average of $m$ infinite.
- *proved, under the theorem's premises*;
  `what_forces_the_ratio_of_means_readout.md` (2 Sept),
  `mass_readout_composition_uniqueness_2026-09-06.md` (6 Sept), and the
  any-two-of-three observation (25 Sept).

**Step 6. The weighting over candidates has two parts.** Before any data, a
candidate is a push for each reading ($\mathbf v_i$) and one angle ($\theta$).
The weighting is written as:

- a **push part** over each $\mathbf v_i$;
- a **mass part** over $\theta$, the starting belief.

This is a way of writing, not a claim. Step 9 says why the push part may not
depend on $\theta$.

**Step 7. The push part is flat in the push vector (axis splitting applied to
the law).**

- **Axis splitting (P6)** gives each axis its own one-dimensional push. For the
  vector reading to agree with its axes, the push weighting must be a product
  over axes.
- **No preferred direction (P7)** makes it depend only on the push's length.
- **A product over axes that depends only on length** is, for $d\ge2$, a
  Gaussian $e^{-|\mathbf v|^2/2\omega^2}$ or an improper $e^{c|\mathbf v|^2}$
  (Herschel–Maxwell).
- **No preferred push size (P7)** leaves only the flat one.
- **1D is fixed the same way.** A 1D reading's weighting is the per-axis factor
  of that product, so a 1D reading must match one axis of a 2D or 3D reading.

What this changes:

- ours21 weighted each push by $r^{2-d}$. That favoured small pushes in 3D,
  favoured large ones in 1D, and agreed with the flat weighting only in 2D.
  Making it flat is the whole of the ours21 → ours24 change.
- Axis splitting also forces the mass part to be the same in every dimension:
  a reading and its axes share one starting belief, counted once.

A further candidate route, proposed 25 September and in no record yet, would
reach the same result without axis splitting. No push size is preferred, so the
weighting is a power of the push's length. A zero push is not special, so it
has finite, nonzero weight at zero, which leaves only the flat one.

This route has two weak spots, which is why it is not relied on (§6):

- "Not special" picks the flat weighting only when it is read in the push
  vector. The 2 Sept record used the same phrase for the magnitudes, where it
  gave ours21.
- The 24 Sept record notes that zero acceleration is physically special, by
  the first law.

- *proved, given P6 and P7*; `excitation_weighting_premise_2026-09-24.md`,
  `axis_splitting_and_the_24_september_change_2026-09-25.md` §4.

**Step 8. The starting belief.**

- **Forced:**
  - it is centred on $s$ (step 2);
  - it is symmetric between $m/s$ and $s/m$ (P4);
  - it is the same in every dimension (step 7);
  - it is counted once, so it washes out as readings accumulate (step 9).
- **Open:** its shape.
- **The current choice** is equal belief for equal angles $\theta$, which is
  the $s/(s^2+m^2)$ factor in the law.
  - It is half of a Cauchy curve with scale $s$: median $s$, about 95% of it
    between $0.039\,s$ and $25.45\,s$, spread of $\log m$ equal to $\pi/2$.
  - It is also the law at zero readings.
  - Its reasons so far:
    - it is uniform in the angle that the noise defines;
    - measured from it, it adds no fake data (step 10);
    - in 1D it makes the law a known default (see §5).
  
  None of these is yet a proof. The first depends on choosing that angle as the
  coordinate, and the second is measured from the choice itself.
- **What rides on the shape:** only the answers at weak signal. Where many
  readings converge does not depend on it.
- *open* for the shape. The forced parts listed above are *proved*.

**Step 9. The combination rule: readings multiply their evidence, and nothing
attached per reading may depend on mass.**

- **The argument.**
  - The log of the joint law is the sum of the readings' log-evidence, plus the
    log of the starting belief, plus the log of any factor $c(\theta)$ attached
    to each reading.
  - As long as the readings keep being pushed, so that the total true push
    grows, the evidence and the per-reading factors both grow in proportion to
    $N$; the starting belief does not.
  - So as $N$ grows, the law peaks where the average log-evidence plus
    $\log c(\theta)$ is largest.
  - The average log-evidence peaks at the true angle (the 25 Sept convergence
    proof, which holds for any per-reading evidence that grows with $h$).
  - Adding $\log c$ moves the peak wherever $\log c$ is not flat at the true
    angle. A mass-dependent per-reading factor therefore drags the answer off
    the true mass, except at the one or two masses where it happens to be flat.
    Eq. (20), for instance, is right only when the true mass is exactly $s$.
  - With no push at all (true zero excitation) the law never concentrates,
    under any rule.
- **What it rejects**, all at once:
  - the old eq. (18) and eq. (20). They multiply whole single-reading laws, so
    they bring an extra mass factor with every reading: $\cos^2\theta$, which
    pulls toward $m\to0$, and $\sin2\theta$, which pulls toward $s$
    (contribution 12 §7);
  - the $(f\alpha)^k$ family;
  - a tilt $\lambda$ used as part of the equation, which would come in once per
    reading.
- **What it fixes:** the push must be measured in both channels together, by
  the length of $\mathbf v_i$. Measuring it by the acceleration alone
  ($\mathbf a^\ast_i=\sigma_a\cos\theta\,\mathbf v_i$) or by the force alone
  attaches a factor of $\theta$ to every reading.
- *proved*; `axis_splitting_and_the_24_september_change_2026-09-25.md` §6,
  `mass_estimator_regimes_priors_and_updating_2026-09-23.md` §7 (the measured
  drift of eq. 18 and eq. 20).

**Step 10. What a tilt is: a reading never taken.**

- **The algebra.** Each reading's log-evidence, as a function of the angle, is
  exactly
  $\tfrac12\big(P_i\sin^2\theta+Q_i\cos^2\theta+2D_i\sin\theta\cos\theta\big)$
  plus a constant, in every dimension. The tilt
  $\lambda\sin2\theta=2\lambda\sin\theta\cos\theta$ is therefore exactly the
  evidence of any reading with equal force and acceleration energies in noise
  units ($P=Q$) and $D=2\lambda$.
- **What that means for $\lambda$:**
  - The fitted $\lambda^\ast=0.35$ is one such reading of an object of mass
    $s$. The weakest one is an aligned reading at signal-to-noise
    $\sqrt{0.7}\approx0.84$ in each channel. In 2D and 3D, stronger readings
    with the same $D$ but less alignment give the same tilt.
  - A negative $\lambda$ is an anti-aligned reading, which pushes estimates away
    from $s$.
- **As a starting belief, counted once,** a tilt is a claim to have seen data.
  That claim belongs in the prior slot as old information, if it is real.
- **Attached to every reading,** it adds one fake reading per real one, all
  pointing at $s$, and step 9 rejects it.
- **The other rejected factors are not exactly readings.** Eq. (20)'s
  per-reading factor $\sin2\theta$ is a product, not a tilt of this form. But it
  also peaks at $s$ and pulls the same way, which is why eq. (20) settles
  between the truth and $s$.
- This does not settle step 8: "no fake data" is measured from the
  equal-angle belief itself.
- *proved*; exact for the ours24/ours25b law in every dimension, checked by
  direct integration (25 Sept).

**Step 11. With many readings, the size of the series is the length of the
stacked vector (axis splitting applied to the answer).**

- **One reading:** its force size is the length of its true force vector.
- **A 2D reading and its two axes:** by the strong form of P6, and P8, the
  joining rule $\oplus$ must satisfy $|F_1|\oplus|F_2|=\sqrt{F_1^2+F_2^2}$ for
  every candidate. So $p\oplus q=\sqrt{p^2+q^2}$.
- **Any series:** by P8, its size is $\sqrt{\sum_i|\mathbf F^\ast_i|^2}$, the
  length of the stacked vector. The same holds for the acceleration.
- **The old plain sum** $\sum_i|\mathbf F^\ast_i|$ fails at the first step.
- **When to pool:** the pooling is done on each candidate, before averaging.
  Pooling the averages instead breaks lumping.
- **Composition still holds.** Within one candidate every reading has the same
  mass, so both composition laws still pass. They cannot tell the plain sum from
  the stacked length; axis splitting can.
- **Steps 7 and 11 are one move:** axis splitting applied to the law, then to
  the answer.
- **Only the weak form is unchecked.** The proof uses the strong form of P6
  (one candidate has one size, however its axes are grouped). A weaker form, in
  which only the final answers agree, is sketched and not checked (§6).
- *proved, given the strong form of P6, and P8*;
  `ours25b_the_many_reading_readout_2026-09-25.md` §3.

**Result.** Steps 3–11 give ours25b, with step 8's shape set to the current
choice.

### 3.3 The chain at a glance

| step | what | status |
|---|---|---|
| 1 | measured and exact inputs are separate; measured inputs always give a mass | split and repair proved; a mass at measured zero–zero assumed (P2) |
| 2 | zero readings give $s$, for any estimator of this type | proved |
| 3 | candidates: true pairs pointing the same way, positive sizes | proved |
| 4 | average over candidates, do not pick the best | proved |
| 5 | readout: average force size ÷ average acceleration size | proved under the theorem's premises (any two of lumping, common force, swap; C3 and the domain) |
| 6 | weighting = push part × mass part | a way of writing |
| 7 | push part flat in the push vector | proved, given P6 and P7 |
| 8 | starting belief: centre, symmetry, dimension-free, counted once | those proved; **shape open** |
| 9 | readings multiply; nothing per reading depends on mass | proved, while the pushes keep coming |
| 10 | a tilt is a reading never taken | proved |
| 11 | series size = stacked length, per candidate | proved, given the strong form of P6, and P8 |

## 4. What follows

- **One reading.** ours25b is ours24. The median and every interval come from
  the law alone, so they are ours24's for any $N$.
- **Zero readings give $s$**, for any $N$ and $d$ (step 2).
- **Equal energies give exactly $s$.** If $P=Q$, meaning the readings carry
  equal total energy in the two channels in noise units, then both the readout
  and the median are exactly $s$. This holds for any $N$, any alignment and any
  split of the energy among the readings. It is step 2 generalised: the swap
  leaves $(P,Q,D)$ unchanged, so $\widehat m/s=s/\widehat m$.
  (`what_ours24_and_ours25b_change_in_earlier_results_2026-09-26.md` §6.1.)
- **Symmetry and units.** Swapping force and acceleration gives
  $1/\widehat m$. A change of force or acceleration units rescales
  $\widehat m$ as a mass.
- **Grouping never matters.** A series and its axis readings give the same
  answer. So does the series in any recording frame, in any order, however the
  readings are grouped. $N$ readings in $d$ dimensions are exactly one reading
  in $Nd$ dimensions, and the answer depends on the data only through $P$, $Q$,
  $D$, $M$ and $s$. (ours25b record §5.)
- **Adding a reading means carrying four numbers.** With $\mathrm i$ the
  imaginary unit and $k$ running over the readings,

  $$
  C=\sum_{k}\frac{(Q_k-P_k)+2\mathrm iD_k}{4},\qquad T=\sum_{k}(P_k+Q_k),\qquad M=Nd,
  $$

  with $\log p(\theta)=\operatorname{Re}\!\big[C\,e^{-2\mathrm i\theta}\big]+\text{const}$
  and $H(\theta)^2=T/2+2\operatorname{Re}\!\big[C\,e^{-2\mathrm i\theta}\big]$.
  Here $P_k$, $Q_k$ and $D_k$ are reading $k$'s share of the three sums. A new
  reading adds its own share of $C$, its $P_k+Q_k$ to $T$, and $d$ to $M$.
  (`what_ours24_and_ours25b_change_in_earlier_results_2026-09-26.md` §2.) The
  old readings' law sits exactly where a starting belief sits; that is the
  prior slot.
- **The peak of the law in $\theta$** is the total-least-squares (Deming) slope
  in noise units, when it lies inside $(0,\pi/2)$. Anti-aligned data put it at
  an end. The readout is a different summary of the same law.
  (`vonmises_state_2026-09-24.md` §2.)
- **With many readings it converges to the true mass**, as long as the total
  true push grows (step 9).
- **Intervals are conditional.**
  - A 95% interval holds 95% of the law. It is not a repeated-experiment
    guarantee, and no interval that is always finite can be one (Gleser and
    Hwang 1987; d1 §8).
  - With many readings the law is too narrow. The true variance of the estimate
    is $1+M/R^{\ast2}$ times the law's variance, so intervals are too narrow by
    the square root of that. Here $R^{\ast2}=\sum_i r_i^{\ast2}$ is the total
    squared true push in noise units. (Axis-splitting record §6; ours25b record
    §5.5.)
- **An honest caveat: padding.**
  - Appending readings that read exactly zero, such as a coasting interval,
    leaves the law and the median unchanged but moves the readout slightly,
    toward $s\,E[\sin\theta]/E[\cos\theta]$. Under the premises this is forced:
    a zero reading still has a noise-sized unknown push, and the stacked length
    counts it.
  - The same effect makes a 1D reading differ from $(x,0,0)$ in 3D.
  - The dimension $d$ is physics, the number of directions the push can take.
    Axis splitting is only about bookkeeping. (ours25b record §7.)

## 5. Where the parts come from

In short, the models are borrowed; the reasons for choosing them, and the exact
readout, are the project's.

- The law is Lindley and Beale's (1966), with a starting belief pinned to the
  instruments' noise ratio.
- In 1D it is Liseo's reference-prior posterior, restricted to $m>0$.
- Uniform-angle priors are known in data units. Measuring the angle in noise
  units, and so centring it on $s$, is the project's.
- Lindley already saw one of the drifts step 9 rejects: a per-point factor
  pulling toward $m\to0$.
- The selection principles are not found in the literature so far: composition
  forcing the readout, and axis splitting forcing the weighting and the
  pooling.

For the part-by-part account, what was read and how, and what can honestly be
claimed, see
`../_dscn_mass_estimator/credit_and_originality_of_the_mass_estimator_2026-09-26.md`.

## 6. Open

- **The shape of the starting belief** (step 8).
- **The weak form of axis splitting** (step 11). It is sketched only.
- **The candidate route to step 7 without axis splitting** ("no preferred push
  size, and a zero push is not special"). It is unchecked, and has the two weak
  spots named in step 7.
- **The readout theorem's own premises** (step 5).
  - C3: whether a readout that composes could use more than the law of the two
    sizes.
  - The domain: the theorem is proved over all positive laws with finite means,
    wider than the laws this equation produces.
- **Interval width with many readings.**
  - The fix could be a finite-temperature (hierarchical) push weighting, or
    widening by the sandwich factor.
  - ours25's sandwich correction currently depends on how readings are grouped.
- **Instruments with different noise ratios**, per axis or per reading. These
  are one corridor, since splitting a reading into axes turns the first case
  into the second. Correlated force and acceleration errors are also untested
  under the current law.
- **Tester work.**
  - The $(C,T,M)$ carry for sequential updating.
  - A world with coasting readings.
  - Rerunning the earlier many-reading ratio-of-means tables with ours25b.
  - A second seed.
- **Real data.**

---

## Sources

Records in `../_dscn_mass_estimator/`, in the order the chain uses them:

| record | used in |
|---|---|
| `mass_estimator_law_assumed_correct_steering_2026-09-07.md` | P3 |
| `disagreement_rule_mass_estimator_trace.md` (29 Aug) | step 1 |
| `inertial_mass_estimator_full_working_spec_d1.md` (30 Aug) §7; `median_based_estimator_inertial_mass.md` (27 Aug) §2; `what_forces_the_chain_to_the_estimator.md` (2 Sept) §5 | step 4 |
| `what_forces_the_ratio_of_means_readout.md` (2 Sept); `mass_readout_composition_uniqueness_2026-09-06.md` | step 5 |
| `excitation_weighting_premise_2026-09-24.md`; `axis_splitting_and_the_24_september_change_2026-09-25.md` | steps 7, 9; §4 |
| `mass_estimator_regimes_priors_and_updating_2026-09-23.md` §7 | step 9 |
| `old_information_in_the_prior_slot_2026-09-24.md`; `vonmises_state_2026-09-24.md` | step 10; §4 |
| `ours25b_the_many_reading_readout_2026-09-25.md` | §2; step 11; §4 |
| `what_ours24_and_ours25b_change_in_earlier_results_2026-09-26.md` | relation to d1; §4 |
| `credit_and_originality_of_the_mass_estimator_2026-09-26.md` | §5 |

Steps 2 and 10, and the any-two-of-three observation in step 5, come from the
conversation of 25 September. The records of 25–26 September list them (ours25b
record §1 and §8; the changes record §§2 and 5; the credit record §§2 and 4),
but this draft is the first place their arguments are written out.

The ours25b record's §1 credits any-two-of-three to the 6 September theorem. That
theorem uses lumping and common force; the swap routes are the 25 September
observation.

**Reproduction.** Tester `canvalidk/Mass-estimator-testor`, commit `e288110`
on top of `19d3b41` (ours25b, committed together with the `stacked_state` and ours26 work): the readout is `ratio_of_means_rss`, and the tests are
`tests/test_stacked.py`. The following were rechecked for this draft by an
independent implementation of the formulas in §2 (26 Sept):

- zero readings give $s$ in 1D, 2D and 3D;
- equal energies give $s$ to about $10^{-15}$;
- the swap gives $\widehat m\cdot\widehat\mu=1$;
- the 21 September worked example gives $4.9937539062$;
- the 25 September worked reading gives $3.8855$ as a 3D vector and as three
  axes;
- a Monte Carlo draw of candidates agrees with the readout formula within
  sampling error.
