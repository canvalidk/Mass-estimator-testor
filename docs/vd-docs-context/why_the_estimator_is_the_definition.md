# Why the Estimator Is the Definition, and the Law Is Its Limit

Standalone note in the `vector_division_mass_estimation` family. Authored
in-repo 2026-08-01 from a working session with Can, after the E10 blind entry
writing pass.

Not a draft of the note. `why_the_mass_entry_must_be_long_joint_measurement_d3.md`
remains `current: true`; this is a companion, filed here on the same grounds as
`why_reconciliation_cannot_be_upstream_of_E10.md`.

**What it does.** The session produced three claims about null trials by
reasoning informally. This note checks them against the actual estimator
equation in `vector_division_inertial_mass_measurement_d1.md` §4, finds
agreement, and then turns the agreement into an argument about which of two
candidate definitions of inertial mass is the real one.

---

## 1. Two candidates for the definition

**A — the law, rearranged.** For one material object,

$$
m=\frac{|\mathbf F|}{|\mathbf a|}
$$

after checking that the two are codirectional.

**B — the estimator.** d1 §4's errors-in-variables fit: the positive scalar for
which the smallest uncertainty-weighted corrections to both measured vectors
make them satisfy $\mathbf F^{*}=m\mathbf a^{*}$.

These are not two notations. d1 §3 shows that any rotation-invariant estimator
built from the two measured vectors has the form

$$
\widehat m=\frac{|\widetilde{\mathbf F}|}{|\widetilde{\mathbf a}|}\,g(\cos\theta),
$$

and that the exact law fixes only $g(1)=1$. A is the choice $g\equiv1$. B is a
different choice, determined by the uncertainties. **They agree exactly on
exactly aligned exact inputs and nowhere else.** So at most one of them is the
definition, and the cases where they disagree are the cases that decide which.

---

## 2. The check: the estimator has a closed form here

For a single trial with isotropic covariances $\Sigma_F=\sigma_F^2 I$ and
$\Sigma_a=\sigma_a^2 I$, minimising d1 §4's objective over the latent
acceleration and then over $m$ gives stationary points satisfying

$$
\sigma_a^{2}\,s\,m^{2}
-\left(\sigma_a^{2}p^{2}-\sigma_F^{2}q^{2}\right)m
-\sigma_F^{2}\,s
=0,
$$

where $p=|\widetilde{\mathbf F}|$, $q=|\widetilde{\mathbf a}|$ and
$s=\widetilde{\mathbf F}\cdot\widetilde{\mathbf a}$. This is the Deming /
orthogonal-distance slope in the form the mass problem takes.

The product of its roots is

$$
\frac{-\sigma_F^{2}s}{\sigma_a^{2}s}
=-\frac{\sigma_F^{2}}{\sigma_a^{2}}<0 .
$$

**Result 1 — the estimator always returns exactly one positive mass.** The
product of the roots is negative for every measured pair, independent of $s$, of
the angle, and of the magnitudes. One root is positive, one negative; the
positivity constraint $m>0$ selects the first, and it always exists.

This is the session claim — *we estimate the most likely codirectional pair, not
the most likely pair, so every sample yields a mass greater than zero* — and it
is not an approximation or a statistical tendency. It is a property of the
quadratic. Numerically: 200,000 null realisations at each of five
$(\sigma_F,\sigma_a)$ settings gave a positive root in 100.00% of cases.

It also disposes of a case d1 §6 leaves to judgement. Antiparallel measured
vectors do not break the estimator; they return a positive mass with a large
residual, which is what §6 says the correct response is.

---

## 3. What a null trial returns

Take a trial with no signal: the object contributes nothing, so
$\widetilde{\mathbf F}$ and $\widetilde{\mathbf a}$ are pure instrument noise,
with typical magnitudes $p\sim\sigma_F$ and $q\sim\sigma_a$.

Substituting $p=\sigma_F$, $q=\sigma_a$ kills the middle coefficient exactly:

$$
\sigma_a^{2}p^{2}-\sigma_F^{2}q^{2}
=\sigma_a^{2}\sigma_F^{2}-\sigma_F^{2}\sigma_a^{2}=0 ,
$$

leaving $\sigma_a^{2}s\,m^{2}=\sigma_F^{2}s$, so

$$
\boxed{\;\widehat m=\frac{\sigma_F}{\sigma_a}\;}
$$

independent of the angle between the two noise vectors.

**Result 2 — a null trial returns the ratio of the two noise floors.** Simulated
medians over 200,000 draws:

| $\sigma_F$ | $\sigma_a$ | $\sigma_F/\sigma_a$ | median $\widehat m$ | IQR |
|---|---|---|---|---|
| 1.0 | 0.5 | 2.000 | 2.005 | 0.83 – 4.84 |
| 1.0 | 1.0 | 1.000 | 0.998 | 0.41 – 2.42 |
| 2.0 | 0.25 | 8.000 | 8.002 | 3.33 – 19.29 |
| 0.1 | 0.05 | 2.000 | 2.008 | 0.83 – 4.87 |
| 3.0 | 7.0 | 0.429 | 0.429 | 0.18 – 1.04 |

This is a detection limit, and it is the physically correct answer: *I could not
resolve a force below $\sigma_F$ or an acceleration below $\sigma_a$, so masses
of order $\sigma_F/\sigma_a$ are what this trial cannot distinguish from
nothing.* The estimate is broad and weak, exactly as d1 §6 predicts for the case
where the uncertainties are comparable to the measured magnitudes — and a null
trial is the extreme of that case, since the magnitudes *are* the uncertainties.

**Which statistic.** The median is the one that equals $\sigma_F/\sigma_a$, and
it is what the table reports. **The mean does not exist.** Occasionally the two
noise vectors come out nearly perpendicular; then
$s=\widetilde{\mathbf F}\cdot\widetilde{\mathbf a}\to0$, the quadratic's leading
coefficient $\sigma_a^2 s$ collapses, and $\widehat m$ blows up like $1/s$. The
resulting tail falls as $1/M$ — measured at a factor of 9.7 to 9.9 per decade
across six decades — so $\int M\cdot M^{-2}\,\mathrm dM$ diverges logarithmically
and the sample mean never settles. At $N=10^4,10^5,10^6,10^7$ it reads 39.2,
13.7, 21.8, 23.6 while the median holds at 2.00 throughout. Report the median and
the interval. The mean of a null trial is not a number.

---

## 4. The disagreement, and what adjudicates it

On a null trial the two candidates return different things, and the difference
is not a matter of convention:

```text
A (law)       0/0        no value; error
B (estimator) sigma_F / sigma_a   a positive mass, broad, with a residual
```

There is an independent standard here that neither candidate supplies: ordinary
detection-limit practice, which d1 already cites in the form of JCGM 100:2008. A
null measurement is a legitimate measurement outcome. It yields a bound set by
the sensitivity of the apparatus. It is not an error, and no experimentalist
reports it as one.

**B agrees with that practice. A contradicts it.**

**Result 3 — therefore the estimator is the definition, and the law is its
limit.** If one accepts that a sensitivity bound is a valid way of determining a
mass, then the operation that produces it is the operation that determines mass.
A is what B reduces to when the uncertainties vanish and the pair is aligned; it
is a special case that happens to be the one classroom problems are built from.
It cannot be the definition, because it disagrees with a legitimate
determination on the case where a legitimate determination exists and A returns
nothing.

Note the shape of the argument. It is not that B is more convenient or more
rigorous. It is that A and B *disagree in truth value* on a reachable input, so
at most one is the definition, and an outside standard settles which.

### The obvious objection, and the reply

*A null trial tells you about your apparatus, not about the object. A is right to
return nothing.*

The premise is correct — §3 shows the answer is $\sigma_F/\sigma_a$, which is a
fact about the instruments. But that is what a detection limit is, and detection
limits are results. The reply is continuity. B returns object-facts when there is
signal, apparatus-facts when there is none, and graded mixtures in between, with
the residual and the interval saying which regime you are in. A has two states:
it works, or it errors, with nothing between them and no way to tell a marginal
trial from a null one. The information content of a trial varies continuously;
only B varies continuously with it.

The same reply covers the non-aligned case, where A also errors and B returns a
broad estimate. The null trial is simply the cleanest instance, because there is
no temptation to blame the data.

---

## 5. Where B fails, and why that is the right place to fail

B is not singularity-free. The quadratic's three coefficients are
$\sigma_a^2 s$, $-(\sigma_a^2p^2-\sigma_F^2q^2)$ and $-\sigma_F^2 s$. All three
vanish when $p=q=s=0$ — that is, when the *readings* are exactly
$\mathbf 0$ and $\mathbf 0$. There the equation degenerates to $0=0$ and the
argmin is the whole ray $m\in(0,\infty)$.

The condition for reaching that point is what matters.

- With $\sigma_F,\sigma_a>0$, an exactly zero reading is a probability-zero
  event. The degenerate point exists and is never visited. Every realisable null
  trial lands in §3.
- Only with $\sigma_F=\sigma_a=0$ does the exactly zero reading become certain
  rather than accidental. That is the case where the failure is forced.

**Result 4 — the singularity has been moved, not removed.** Both candidates
carry a $0/0$. A's sits at $(\mathbf F,\mathbf a)=(\mathbf 0,\mathbf 0)$, which
is the commonest physical state there is — every object at rest, every object in
uniform motion. B's sits at $(\sigma_F,\sigma_a)=(0,0)$, a perfectly precise
null measurement.

And it is not removable by taking a limit. Holding the ratio fixed while
improving both instruments leaves the estimate exactly where it was:

| precision | $\sigma_F/\sigma_a$ | median $\widehat m$ |
|---|---|---|
| $\times1$ | 2.0 | 1.993 |
| $\times10$ | 2.0 | 1.995 |
| $\times100$ | 2.0 | 1.997 |
| $\times1000$ | 2.0 | 2.003 |

So $\lim_{(\sigma_F,\sigma_a)\to(0,0)}\widehat m$ is path-dependent: approach
along a fixed ratio $r$ and the answer is $r$. It is a genuine $0/0$, of the same
algebraic type as A's, relocated from the data to the instrument.

**Why the relocation is a cure.** A perfectly precise measurement is not a
physical state. We should expect a floor on the joint precision of the two
channels: an acceleration is a second derivative of position, so resolving it
arbitrarily well requires position information that the uncertainty principle
does not permit alongside the momentum information the force channel needs. The
exact form of that floor is not settled here and is not needed — thermal and
instrumental floors sit far above the quantum one in any case. What is needed is
only that $\sigma_F=\sigma_a=0$ lies outside the domain where classical mechanics
is the applicable theory.

A definition whose only singularity lies outside its own theory's domain of
validity is in good shape. It is not a defect that a classical definition of
mass breaks down exactly where the classical description breaks down; it is what
one would want. Candidate A, by contrast, breaks down on a book resting on a
table.

*Flagged, not claimed:* the quantum argument is stated as an expectation. It
would be worth someone deriving an actual bound on $\sigma_F\sigma_a$ or on the
attainable $\sigma_F/\sigma_a$, because if such a bound exists it says there is a
smallest mass this operation can fail to detect, which would be a real result and
not merely a boundary disclaimer.

---

## 6. Consequence: what a free body is, operationally

Carried from the working session and now grounded in the equation.

A free body in an inertial frame has zero net force and zero acceleration, so it
is the null trial of §3. Its mass estimate is $\sigma_F/\sigma_a$ — a fact about
the apparatus and not about the object. This gives an operational test that
distinguishes it from a badly measured accelerating body, which is also
uninformative but for a contingent reason.

**Improve one instrument and see whether the answer moves.**

| | $\sigma_F$ | $\sigma_a$ | median $\widehat m$ |
|---|---|---|---|
| free body, force sensor improved | 2.0 → 0.25 | 0.5 fixed | 4.00 → 0.50 |
| free body, both improved, ratio fixed | ↓ ×1000 | ↓ ×1000 | 2.00 → 2.00 |
| real mass $m_0=3$, both improved ×64 | ↓ | ↓ | 3.0000, IQR 1.26 → 0.019 |

For a real mass, tightening either channel narrows onto the same value. For a
free body the estimate tracks $\sigma_F/\sigma_a$ and moves in *opposite
directions* depending on which instrument you improved. That is observable, and
it is a sharper claim than "a free particle carries no information about its
mass", which is not directly checkable.

Two guards on this, unchanged from the session.

1. **It runs one way.** Free implies apparatus-tracking; apparatus-tracking does
   not imply free, since a badly instrumented accelerating body also tracks the
   apparatus — but only until the instruments are improved, which is what the
   test exploits. Inside an inertial frame the biconditional does hold, and the
   frame condition is Newton I's own precondition, which is why this result sits
   in Newton I's periphery rather than reaching into it.
2. **Consequence, never definition.** If `free-body` were defined by Newton II's
   silence, Newton I would reference a headword that does not exist yet at that
   point in the log.

This is not a Newton I / Newton II collapse. It is the opposite: it says that
exactly where Newton I applies, Newton II is silent. The two laws partition
rather than nest, which is independent support for the pedagogy note's refusal to
teach Newton I as the zero-force case of Newton II — support that note does not
currently have. It also means a free photon and a free cart are indistinguishable
to Newton II, so the free case never has to adjudicate masslessness.

---

## 7. What this note does not establish

- The closed form in §2 is for a **single trial with isotropic covariances**.
  The product-of-roots argument and Result 1 depend on that form. The general
  anisotropic, correlated, multi-trial case should be checked before Result 1 is
  quoted without qualification, though nothing suggests it fails there.
- Everything here is **inside a declared policy**. §3's $\sigma_F/\sigma_a$ is
  what the Gaussian errors-in-variables policy returns; a different policy
  returns something else. That is d1 §3's non-uniqueness, and the null trial is
  where it bites hardest, because the data constrain nothing and the answer is
  therefore *all policy*. This is an argument that d5's requirement to record the
  policy matters most exactly where a naive reading says there is nothing to
  record.
- It does **not** establish the retention requirements. Those rest on the
  self-confirmation argument in d3 §§3–4, a separate chain.
- It does **not** license reporting a bare point value from a null trial. The
  number is only interpretable next to its interval and residual; without them
  $\sigma_F/\sigma_a$ reads as a mass and is not one.

---

## 8. Consequences for E10 and the packet

1. **`ZeroZeroUnderdetermined` is an exact-branch result only.** It is correct
   for stipulated exact zeros, where there is no uncertainty and no estimation.
   It is wrong for measured null inputs, which return a positive broad estimate.
   The packet currently states it without that restriction.
2. **No zero clause is needed in the entry.** A measured null reading is never
   literally $(\mathbf 0,\mathbf 0)$; it is a small vector pair that is not
   codirectional, so it takes the estimation branch as an ordinary small-signal
   case. The exact zeros reach $0/0$ through the first clause. The two are
   already separated correctly with no special casing.
3. **The `if` in E10 is a capability branch, not an operation branch.** §1 and §4
   together say the estimator is total and the law is its aligned zero-noise
   limit. The exact clause is a fast path through the general operation, written
   first for the reader's sake. A later editor who keeps it and deletes the
   general clause lands on horn one of the reconciliation dilemma.
4. **The free-body result is wall or Newton-analysis material**, not definiens
   material: free-particle vocabulary is on E10's forbidden-mentions list and the
   law writing guide §5 bars a triplet entry from reaching outward.

---

## Ratification status

**Checked, not merely asserted.** Results 1–4 were derived from d1 §4's objective
and confirmed numerically (200,000 null realisations per setting; the control run
with a real mass converges to $m_0=3.0000$ with the IQR shrinking as
$1/\text{precision}$). The arithmetic is reproducible from §2's quadratic.

**Argued, not ratified.** The claim that B is therefore the definition and A its
limit. It depends on accepting detection-limit practice as a legitimate
determination of mass, which is stated here as a premise, not defended.

**Expectation, flagged.** The quantum-boundary reading of $\sigma_F=\sigma_a=0$.
The structural point — that B's singularity is unreachable while A's is ordinary
— does not depend on the quantum argument being right in detail, only on perfect
precision being unattainable for some reason.

---

## Intended catalog disposition

```text
why_the_estimator_is_the_definition.md
  note_name: vector_division_mass_estimation
  draft: null          # not a draft — standalone companion
  current: null
  role: primary
  informs: -> E10 entry writing and packet (link only)
           -> Newton I periphery / inertial-mass wall (free-body result, sec 6)
```

`why_the_mass_entry_must_be_long_joint_measurement_d3.md` remains `current: true`
for this note family.
