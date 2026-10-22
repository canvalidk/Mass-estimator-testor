# Why Reconciliation Cannot Be Upstream of E10

Standalone correction note against `ENTRY_PACKET_E10_INERTIAL_MASS_MATERIAL_OBJECT_d4.md`.
Authored 2026-08-01, in-repo, from a working session with Can.

Not a draft of this note. `why_the_mass_entry_must_be_long_joint_measurement_d3.md`
remains the current draft; this is a companion that supplies one argument d3 does
not make and settles part of an open question d3 and the d4 packet both leave
standing. Filed here because it is about the same object — the mass corner of
Newton II — and belongs with the material it corrects.

---

## What this note corrects

The d4 packet establishes that empirical mass determination runs a joint
Newton-II-constrained fit rather than a magnitude quotient. It does not establish
where that fit *lives*. Decision 11 leaves estimator ownership capability-named,
and open question 3 asks whether E12 is the joint estimation operation, the
mass-address wall, both, or neither.

Underneath that open question sits an assumption neither d3 nor d4 challenges:
that reconciliation *could* in principle be a prior operation — performed by
metrology, by an evaluator pre-pass, or by a separate numbered entry — which then
hands E10 a conforming pair. On that picture E10's definition can stay exact and
short, and the messy part is somebody else's problem. Call this the **deferral
picture**. It is the picture most textbook treatments use, and it is what makes
the mass entry look as though it ought to be repairable by better wording.

The deferral picture is not merely incomplete. It is unavailable. Not because
deferring is bad practice, but because there is no consistent way for any
upstream operation to deliver what it would have to deliver.

---

## The argument

### 1. The triplet's three operations, and what each requires of its inputs

Newton II is addressed by three entries over the same three quantities:

```text
E8   inertial-acceleration-material-object
E9   net-force-material-object
E10  inertial-mass-material-object
```

Two of the three run as total operations on their input types:

$$
\mathbf a = \frac{\mathbf F_{\mathrm{net}}}{m},
\qquad
\mathbf F_{\mathrm{net}} = m\mathbf a .
$$

Read them as operations rather than as rearrangements of one equation, because
that is how the entries are used.

**E9 takes acceleration as an input** and returns net-force. For that to be
runnable, acceleration must be identifiable *without* reference to net-force. If
what counted as the object's acceleration depended on the direction of the net
force, then supplying E9's input would require E9's output.

**E8 takes net-force as an input** and returns acceleration. Symmetrically,
net-force must be identifiable without reference to acceleration.

Neither quantity may be defined by reference to the other. This is a requirement
of the triplet's internal operability, not a methodological preference.

### 2. The channels are in fact independent

The requirement is also satisfied in practice, which matters because it removes
the escape route of saying that real measurement never separates them anyway.

A net force is determined from the force side alone: summed component forces from
load cells, spring extensions, known field strengths and calibrated interactions.
No trajectory data enters.

An acceleration is determined from the trajectory side alone: time-indexed
position samples differentiated, or accelerometer output. No force data enters.

So the two quantities arrive from genuinely separate channels, each carrying its
own instrument model, calibration and uncertainty account.

### 3. Independently determined vectors are generically non-conforming

Two vectors produced by unrelated channels, each carrying its own error, will not
in general stand in the exact relation

$$
\mathbf F_{\mathrm{net}} = m\mathbf a,
\qquad m>0 ,
$$

for any positive scalar. Exact codirectionality of two independently estimated
vectors is a measure-zero coincidence. Displayed alignment, when it occurs, is an
artefact of finite precision rather than evidence that the pair is exact.

This is d3's result and is not re-argued here.

### 4. The dilemma

Now ask what an upstream reconciliation step would have to hand E10.

**Horn one — it delivers the two quantities as its channels determined them.**
Then the pair is generically non-conforming, and an exact definition of mass has
nothing to apply to. The definition is not wrong; it is *inert over the entire
domain on which traces actually run*. Deferral in this form does not divide the
labour. It declines to do any.

**Horn two — it delivers a reconciled pair.** Then fixing the acceleration
required the direction of the net force, and fixing the net force required the
direction of the acceleration. Both quantities are now defined partly by
reference to the other. E8 and E9 immediately become inoperable: each would need
its own output in order to obtain its input. The reconciliation buys E10 a clean
pair by destroying the other two corners of the triplet.

There is no third delivery. Whatever an upstream step passes on is either
unreconciled — and then useless to an exact definition — or reconciled, and then
the reconciliation has already corrupted the quantities the other entries depend
on.

### 5. Conclusion

Reconciliation cannot be upstream of E10, because anything upstream acts on the
net force and the acceleration *before* E8 and E9 receive them, and any such
action changes what those quantities are.

Therefore the estimation of a conforming pair is **internal to mass
determination**. It is a constitutive part of what recovering a mass from
empirical inputs consists in, not a preprocessing step that happens to precede
it.

Stated as the result:

> Inertial mass is the positive scalar relating an *estimated* codirectional
> pair, and the estimation is part of the determination rather than prior to it.

---

## Relation to the argument in d3 §3

d3 §3 already rules out defining a measurement as one that has been reconciled
with Newton II. Its objection is **epistemic**: reconciled inputs make the force
and acceleration channels mutually dependent, so the law is built into the data
used to apply it and can appear to confirm itself.

The objection above is **structural**: reconciled inputs make E8 and E9
inoperable, because each entry's input would presuppose its output.

The structural objection is the stronger of the two for three reasons.

1. It does not depend on any position about confirmation, evidence or the status
   of the law. It goes through even for someone who is untroubled by
   self-confirmation.
2. Its premises are internal to the triplet. E10's constraints are derived from
   what E8 and E9 require, which is where a triplet's constraints ought to come
   from.
3. It closes the deferral picture, which the epistemic objection leaves open. A
   defender of deferral can accept d3 §3 and reply that metrology will simply be
   careful. Horn one blocks that reply: careful independent metrology produces
   exactly the non-conforming pair that an exact definition cannot consume.

d3 §3 is not superseded. The two objections are independent and both hold.

---

## Consequences

### Open question 3 is settled by dissolution

The joint fit cannot be owned by an entry that runs *before* E10 — that is the
ruled-out case above, whatever the entry is numbered. But the remaining
alternatives were both malformed, because both assumed the estimator must have a
home *inside the dictionary*.

It does not. **The estimation capability exists, but it lies outside the scope of
the VD's authority — as differentiation does.**

E8 returns an acceleration. Obtaining one in practice means differentiating
position data. The VD does not define differentiation, does not number it, and
does not own it; it is external mathematical competence the user brings. Nothing
about that makes E8 defective or incomplete. The same holds here: estimating a
latent codirectional pair is an external competence, invoked within the mass
determination, not a node in the trace.

So:

- E12 as an upstream estimation entry: **ruled out** (the dilemma).
- A separate measurement entry feeding E10 a conforming pair: **ruled out** (the
  dilemma).
- A dictionary capability E10 owns and invokes: **ruled out** (nothing in the
  dictionary to own).
- E10 teaching the estimator itself: **ruled out** (outside the authority).

What E10 must do is **name the operation, not supply it**. Decision 11's
capability-naming was the right instinct, but not as a deferral pending the E12
ratification — a capability-name is simply what a reference to an external
competence looks like, and it stays that way permanently.

This also dissolves the pedagogical tension. A first-year uses acceleration
without being able to derive a numerical differentiation scheme. A first-year can
likewise hold "the true pair must be estimated" without meeting orthogonal
distance regression. That is not a gap in the entry; it is the ordinary
relationship between a dictionary and the mathematics its users bring.

**One rider.** Unlike the textbook case of differentiation, estimation here is
non-unique in a way d1 records: magnitude-ratio, projection and componentwise
estimators disagree, and the $g(\cos\theta)$ result shows the disagreement is not
a detail. So the trace must record *which policy was applied*, where a trace
invoking differentiation is not usually asked to. The packet already provides for
this — policy and provenance are both in decision 6's joint-result contents — so
nothing new is required, but the requirement's source is worth naming: it comes
from the non-uniqueness of the external competence, not from the mass corner.

### E10's payload becomes derived rather than asserted

The argument yields what the entry must contain, rather than asserting a length:

1. the constitutive relation;
2. that its two inputs are independently identified and therefore generically
   non-conforming;
3. that the reconciliation is located inside the mass determination.

Whether that is best delivered as one long definiens, a short definiens plus
IRIL, or some other split is a wording question and is not decided here.

### Terminological care when stating the result

The result should not be stated as a ratio *between two vectors* — that is the
operation the whole note denies exists. Once codirectionality has been secured by
estimation, what is read off is the ratio of the two **magnitudes**, or
equivalently the positive scalar relating the estimated pair.

---

## What this note does not establish

It proves that reconciliation is internal to mass determination. It does **not**
prove that the residual, the original measurements, the provenance or the
independence status must be retained in the result.

That is a separate theorem, and it rests on the separate self-confirmation
argument in d3 §§3–4: the latent outputs satisfy Newton II by construction, so
the discrepancy is the only thing still carrying empirical content. The two
chains should be kept apart. Conflating them is what lets requirements that
belong to the evaluator migrate into the definitional core.

Nor does it establish anything about cross-trial invariance, or about what
licenses treating the recovered scalar as a property of the material object
rather than of the trial. That gap is noted, not closed.

---

## Ratification status

Recorded so that later passes can tell what has been settled from what is merely
written down.

**Ratified 2026-08-01.** The E8/E9 input-independence requirement; that any
reconciliation making acceleration depend on net-force direction is incompatible
with the three operations sharing the same quantities; and the conclusion that
estimation is an innate factor of determining mass rather than a prior operation.
This is the argument. It is settled and should not be reopened.

Also ratified the same day: that the estimation capability exists but falls
**outside the scope of the VD's authority**, on the differentiation precedent —
which settles open question 3 by dissolving it rather than by assigning an owner.
E10 names the operation; it does not supply it, and no dictionary entry owns it.

**Written up but not ratified.** The two-horn presentation; horn one (inert over
the trace domain) as a horn distinct from the circularity horn; the
structural-versus-epistemic comparison with d3 §3; the partial settlement of open
question 3; and the three-item payload. These are expositions and consequences of
the ratified argument, not additional commitments. Treat them as revisable
without treating the argument above as revisable.

**Refuted in session, recorded so it is not raised again.** That deferring
measurement to metrology is a coherent division of labour which VD merely cannot
afford. Horn two refutes it. The division of labour was never coherent —
conventional treatments do not notice because they do not trace, so the three
directions of $\mathbf F_{\mathrm{net}} = m\mathbf a$ are never required to share
the same quantities explicitly.

---

## Intended catalog disposition

```text
why_reconciliation_cannot_be_upstream_of_E10.md
  note_name: vector_division_mass_estimation
  draft: null          # not a draft — standalone companion
  current: null
  role: primary
  informs: -> E10 entry writing and packet (link only)
```

`why_the_mass_entry_must_be_long_joint_measurement_d3.md` remains `current: true`
for this note. This file is deliberately outside the draft line, on the same
grounds as `why_the_mass_entry_cannot_be_short.md`: it is not a successor pass
over d3 but a differently aimed document, and ORGANIZATION.md directs that a
`dN` be assigned only when the order is confident.
