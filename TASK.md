# The task

**Find the function that is unique.** Not a function that is unique given someone's
preferences, and not a good estimator. The function which is unique: what survives when
everything that is not forced has been taken out of it. That function, if it agrees with
experiment, is the definition of inertial mass. Everything else depends on the situation,
and so belongs to the input.

This file explains the principle, because it is the task. Read it before working on any
open item. (Ruled by Can, 2026-10-24: `rulings.md`, R-2026-10-24-3.)

## 1. The logarithm

Is the logarithm a unique function on the complex plane? No. exp(z) = exp(z + 2πik) for
every integer k, so exp forgets the winding number, and "the number whose exp is w" has
infinitely many answers.

The mistake is to conclude that there is no logarithm, or that there are infinitely many
equally good ones and the choice is a matter of taste. The right conclusion is that the
logarithm needs one piece of information that its argument does not carry: which branch.
Make that piece an input, and the function is unique:

```
log : (w, branch) -> z     unique
log : w -> z               not a function
```

Nothing was chosen inside the function. The freedom was moved to where it belongs, the
caller, who knows which branch the problem needs. A principal branch is a convention for
the common case. It does not make other branches wrong, and nobody says complex analysis
has several logarithms.

Two features of the example carry over exactly:

- **The function cannot manufacture the missing piece.** No computation on w alone recovers
  k, because exp destroyed it. An attempt to pick the branch from w (the nearest to the
  real axis, the smallest modulus, …) is a bolted-on clause, and it shows up as a
  discontinuity somewhere.
- **Give it nothing and it returns the input.** The only thing log can say about the
  imaginary part beyond what w fixes is which branch you asked for.

## 2. The same move for mass

Mass from measured force and acceleration is not a function of the readings alone. Two
things are missing from any series of readings:

- what the measured numbers' errors are (the instrument), and
- what was believed about the mass before any reading (the prior).

So the definition takes them as inputs. The type in `definition/type.md` is that move
already made:

| Fork the search hit | Can the readings decide it? | Where it went |
| --- | --- | --- |
| The noise of each reading | No: one reading cannot tell signal from noise | input Σᵢ (L18: the declared σ's are correct) |
| Which components were measured together | No: same numbers, different declaration | input, the reading structure (ruled 28 Sept, 23:41) |
| Which directions the signal could have taken | No | input Nᵢ (proposed, L18) |
| The averaging window | No | input Δtᵢ (10 Oct) |
| The starting law π: its centre and its shape | No: by L19 nothing about signals may be inferred from the series, and by Z a null series must return what it was given | input π. The half-Cauchy at c = s is the principal branch |
| How strongly the series contradicts itself | Yes: computed from the readings' own residuals | inside the function: η (T31) |
| How misaligned a single reading is | Yes: its own off-line misfit | inside the function: g(Q) (T30) |
| The readout | Must be yes | inside the function: O2 |

The all-zero experiment is the null input. The function must return the centre of π
there, as log returns the branch, because there is nothing else it could honestly say.
That is why "mass is never undefined" (29 Sept) and "the value is c" are the same
statement once π is an input.

## 3. The procedure at a fork

Every time the work reaches a choice that the premises do not settle, sort it before doing
anything else. There are three outcomes and no fourth:

1. **Forced.** Some premise, applied carefully, decides it. Find the argument; write it into
   the theorem it belongs to. This is "hand on the left wall".
2. **Input.** The readings cannot decide it, because the information is not in them. Move it
   to the type, declared before data, and make it travel with the result. It is now
   situation, not definition.
3. **Wrong premise.** Neither: the readings should decide it, and no argument does. Then one
   of the premises is mis-stated. Find which. This is where breakthroughs come from (30
   Sept: the property that the law carries its own uncertainty was the wrong premise, and
   dropping it produced ours29).

What is *not* an outcome: picking the option that scores best in a study, or asking Can to
choose. A study can check a sort; it cannot make one. Can rules on whether a sort is right,
not on which side of a fork to take.

## 4. Believing it exists

The search assumes the unique function exists. That assumption is licensed, and it is the
one assumption that is never dropped. It is what lets you choose when two large networks of
argument contradict each other ("A says B is wrong, B says A is wrong"): if the function
exists, then at least one network contains a mis-stated premise, and the contradiction tells
you to look for it rather than to give up or to split the difference.

The assumption directs the search. It never discounts a proof. If someone proves that three
requirements cannot hold together, the conclusion is outcome 3, that one of them is
mis-stated; it is never a clause bolted on to make them hold.

## 5. The guards

Moving things into the input is powerful enough to abuse. Two checks stop that:

- **Inputs cost.** An input must be something real measurements need and real practice
  supplies: an instrument's specification, a declared prior, a window. It must be declared
  before the data and reported with the result. If a fork is being moved to the input only
  because it is hard, that is a dodge; sort it as outcome 1 or 3 instead.
- **It must agree with experiment.** The definition is whatever is unique *and* reproduces
  what practice does when the inputs are known: detection limits on null measurements,
  combination with other routes (the density route, a spring, weighing), and real data.
  Until the function has been run on real data, this check has not been made. The robot
  data in `legacy/vd-newton-lab/raw/robot/` is the first place to make it.

## 6. What this does to the open items

| Item | Under this principle |
| --- | --- |
| F2 ("any free choice is an alternative mechanics") | True of choices the function makes. A choice the caller makes is an input, and two labs with different π no more have different mechanics than two people using different branches of log have different analysis. |
| O6, π's shape | An input. The open question becomes: what is the admissible family (already known: m⁻¹g(ln(m/c)), g even) and which member is the principal branch. |
| O11, π's centre | An input, with c = s as the principal branch for one noise pair. |
| O7, the scope of Z | Follows from the type: Z is about what the function may do with information, so what counts as "no information" is fixed by what the inputs declare. |
| O2, the readout | Must be inside the function, so it must be forced. The spec is now sharp: compose, ignore what carries no information, and return the input's centre on null input. |
| The uniqueness theorem | Show that every freedom the properties leave is prior information (so an input by L19), never a choice internal to the function. "Every freedom is a prior" (Can, 24 Sept), proved. |

Where the analogy strains: a branch cut changes no observable, and π does change answers at
weak signal. So the result is not a bare number. A mass is (m̂, law, π, Σ), as a complex log
is (value, branch). That is the type being honest, not a defect.

## 7. The order of commitments (27 Oct)

VD appends commitments in order, and an entry is complete when appended: there is no
deferred body in the log. The entry for mass is `mass := mass_estimator(F, a)`, written out
in full at Newton II. The "tbd" is ours, now: we are determining what that body must say
before it is written. What determines it is the commitments that *will* follow in this
version (`mass is additive`, Newton III, …): the body has to be the one function that those
later commitments leave, so the estimator's form is a theorem of the whole version, worked
out ahead of the entry and written into it. That is F4's content and its scope: commitment
order is suspended for the *work* of determining the body, using commitments later in this
version's order; it is not suspended in the log. Energy is not among them; a version with
energy before mass is a different version (R-2026-10-27-2).

The uniqueness theorem is then a ladder: each later commitment narrows the family of
admissible bodies, and the definition is reached when one member is left. At the Newton II
stage alone the family is larger; that more general function is not the one sought. Two
consequences. A freedom the ladder does not remove is, by §3, an input or a mis-stated rung,
never a choice. And whether a rung is in reach is a fact about the version's order, to be
looked up, not a ruling to request.

## Sources

- Can, 2026-10-24, in conversation: the logarithm and branch cut; "the function which is
  unique"; believing it exists (`rulings.md`, R-2026-10-24-3).
- Can, 2026-09-24: "every freedom is a prior" (R-2026-09-24-3).
- Can, 2026-10-27: the order of commitments; the estimator as a theorem of what follows the
  headword (R-2026-10-27-2); moves are forced or "unknown yet" (R-2026-10-27-1).
- Can, 2026-09-30: dropping the property that the law carries its own uncertainty
  (R-2026-09-30-3), the worked example of outcome 3.
