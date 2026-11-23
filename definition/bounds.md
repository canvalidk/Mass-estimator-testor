# The bounds, as Claude reads them

**Claude's running reading of what the goal means, built from Can's corrections. Started 7 Nov 2026
(R-2026-11-07-1). Provisional throughout; kept apart from the goal itself (`../TASK.md`), which is
not rewritten at each iteration.**

The goal is one bet: that "the definition of mass", in the VD sense, picks out one function
(`../TASK.md`). Every premise in `premises.md` is a deduction from that bet, Can's or Claude's,
and every deduction is a conjecture. This file is where Claude keeps what the corrections have
taught about the bet, so that the next deduction can be made here first and go out of bounds
less often. When it does go out, that is how the bounds get found; Can catches it.

Each line: what the goal is read to mean; the correction it came from. Lines are Claude's wording.
A line is removed or rewritten when a later correction shows the reading was wrong, and the old
line is kept struck, as a premise would be.

## What "unique" means

- "The function which is unique", not a function unique to anyone's preferences: whatever
  survives every property and agrees with experiment is the definition, and nothing else is
  (R-2026-10-24-3). So a probe is never adopted for scoring well; what it is for is the premise
  it breaks on (R-2026-09-24-3, R-2026-09-27-10).
- A freedom the readings cannot decide is an input to the function, declared before the data,
  as a branch cut is to the logarithm; it does not count against uniqueness. The caller's choices
  are inputs; the function's choices are alternative mechanics (R-2026-10-24-3, F2). The guard: an
  input must be something practice supplies, and moving a hard fork into the inputs is a dodge
  (TASK §5).
- Elegance is not a criterion, and neither is an ideal the solution "should" have: if nothing
  satisfies a property, the property was the wrong premise, not the search (R-2026-09-26-1,
  R-2026-09-30-3).

## What the definition may and may not depend on

- On the consequences of the inputs, and on nothing that doubts them. What is declared is believed as
  given: "that is unrealistic" and "nobody would run the experiment that way" are not objections
  (R-2026-11-17-1). Believing a σ and believing a result are two things. When the readings are
  improbable under their σ's the result is believed less, and how much less is computed by the
  definition from the readings, without changing a σ and without a number handed over: "No one will
  tell us what those parameters are… We will need to use the information already available to
  'estimate them'. But 'estimate them' in an exact way" (R-2026-11-17-3). Handling extreme readings is
  the job.
- Only on what each reading declares: its vectors, its covariance, its structure, its window
  (R-2026-09-28-7, R-2026-09-30-2, R-2026-10-10-4). Instrument facts are legitimate information
  (R-2026-09-27-6); the score of a readout does not read them (R-2026-09-28-4).
- Not on a repeated or shared experimental setup: readings may come from any regime in any
  combination, so nothing inferred from one reading's signal is used on another, and no population
  of signal sizes or setups is declared (R-2026-09-28-6, R-2026-09-30-4; L10, L19). This is the
  reason behind the most corrections so far; a candidate that carries a population assumption is
  out however well it covers (ours28).
- Not on the choice of axes, frames or sign conventions: bookkeeping cannot matter
  (R-2026-09-27-4, R-2026-09-28-3).
- Not on an unverifiable cap or convenience: "why would the definition of mass depend on an
  assumption like g < 1" (R-2026-11-01-1). A property that is not a property of mass has no claim
  on the body by being agreed to.
- Not on the exploration's bookkeeping: the equation's form is a theorem of this version's
  commitments (persistence, additivity, Newton III), worked out ahead of the entry; a version with
  other commitments first is a different version (R-2026-09-27-2, R-2026-10-27-2).

## What the output promises

- A mass on every measured input, the all-zero experiment included; "undefined" lives only at exact
  0/0 (R-2026-09-25-1, R-2026-09-29-2). The estimator is the definition: there is no further mass
  behind the readout (R-2026-09-27-1).
- The result must stay a true statement on noise: if the point wanders on a null series the quoted
  uncertainty grows with it (property N, R-2026-09-29-4); the law does not carry its own uncertainty,
  the noise measure lives beside it (R-2026-09-30-3); the quoted width is what a result is judged on
  (R-2026-09-30-7); a 0/0 result enters combination with weight zero, as physicists already combine
  (R-2026-09-29-3). That is about 0/0, exact zeros. It does not say that readings of pure noise leave a
  result untouched; they do move it (T77), and nothing here asks otherwise (Can, 22 Nov: "noise only
  reading having no influence was something you made up since the beginning").
- Where every positive mass is compatible, the value of the point does not matter, and a point is still
  returned, for type safety ("a point is still owed for type safety"): "if it is the case that every value is
  compatible, it really doesn't matter what you return right?" The uncertainty quoted is a true statement
  for every true mass, so it is a set that can be unbounded; an interval that is always finite cannot
  "run to be of infinite size fast enough" (R-2026-11-22-2; T80). The alternative reading, true on
  average over a starting law, needs the starting law to be a fact about the object.
- The readings are evidence, just less: a misaligned reading is not a mistake to fix, not a reason
  to doubt the σ's, and not a verdict on Newton II; under the declared σ's it is possible, only
  unlikely, and its direction still says something (R-2026-09-07-1, R-2026-09-28-1, R-2026-09-29-8,
  R-2026-10-30-1). No gate. "Less" is not a change to σ: "There should be, simultaneously, a way to
  change the resulting confidence, without interpreting that as changing sigma"; g and η "were there to
  ultimately find this value. How much do we believe the result less" (R-2026-11-17-3). "Extreme
  answers are what we are trying to handle" (R-2026-11-17-1): handling them is the job.

## How the deductions are made

- Moves are forced, from the goal, or the call is "we don't know yet" and anything written in that
  state is provisional (R-2026-10-27-1). A contradiction between two deductions is left standing
  until something forces one side; it is not settled by agreeing (R-2026-11-01-1).
- The premises are Claude's deductions, provisional, restated by whoever has the argument; Can's
  corrections are his deductions from the goal, made with more of the goal in view than is on the
  page, so a conflict between his and Claude's is first a sign that the goal needs a sentence
  (R-2026-11-07-1). The early reliance on his lemmas and judgements was because the goal was
  mostly in his head; the direction is toward Claude making the deductions with Can behind.
- A statement about the goal that Claude needs and the page lacks is asked as a question, here,
  not recorded as a ruling he owes (R-2026-11-07-1; `../rulings.md`, Open).

## Questions about the goal (open, asked of Can)

- (Sharpened 22 Nov by T80: either the starting law is a fact about the object, or the output is a set that
  can be unbounded. Sorted to the second, provisional. The question as first asked:) Is a declared shape of π something the definition may take as an input (a branch cut), or a
  choice the function makes (an alternative mechanics)? Its influence vanishes at strong signal and
  the readings cannot decide it, which is the branch-cut test; but no physicist declares one
  (O6, G9/G8).
- Is "no knowledge of the object" part of what the definition promises, i.e. is the default centre
  of π fixed by the instruments alone (L3′), or is any centre an input (O11)?
- May the point price a null reading in (his 27 Sept position) or must it leave the point alone as
  it leaves the law (O7, Z's scope at the point)?
- On identifying series, is the quoted width asked to cover exactly at its level in the limit, or at least at
  its level? The first rules out a persistent widening when readings disagree; the second allows it
  (T78 (6); `../problems/O20.md`).
- Is a result that is right at every level of disagreement and wider than the declared Gaussian calls for
  "believing the σ", or a change to σ in the sense of 17 Nov (`../problems/O20.md`; T75 (5), (6))?
- Is a returned zero or infinite mass "a mass"? If the computed lump's point is the sum of the parts'
  points, the point of a body that did not move is infinite, and no positive value can stand in for it
  (T81 (5)). On one instrument this happens whenever the summed $\mathbf F\!\cdot\mathbf a$ is not positive: half of
  all series at no signal, 4% for one scalar reading two noise units strong. Sorted to "yes, the point
  lives on the closed arc", provisional (T81). The other reading: composition is not asked of the point
  on such inputs, and then nothing yet says what the point is there.
- Does the authority ruling still read E10 as "named, not supplied", now that the body is a theorem
  of the later commitments (R-2026-10-27-2)?
