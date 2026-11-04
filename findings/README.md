# Findings: how an agent reports what it found

Every run by another agent (Codex and the like), briefed or exploring, ends with **one findings file**:
`findings/<branch name>.md` (slashes as dashes: `findings/codex-b15.md`). It is the only thing Claude
reads first. Theorem files, attempts, code and tests are the evidence behind it; Claude opens them only
when a card points there. A finding that is not on a card is treated as not found.

## The file

```
# Findings: <branch>

**Run:** <brief id, or "exploration: <target>">  **Base:** <commit>  **Head:** <commit>
**Suite:** <passed / failed / xfailed, or which subset and why>  **Audit:** unchanged | <what changed>

## Headline
Up to five card ids, most consequential for the definition first, one line each.

## Stop
Exploration only: the last break named (a) a premise of the definition or (b) part of this run's own
construction. If (b), the run stopped there. One line.

## Cards
<one card per finding, in the headline's order, then the rest>
```

## A card

At most twelve lines. No narration: the argument lives in the theorem or problem file the card names.

```
### F<n>: <the claim in one sentence, with its numbers>
- **Kind:** break | no-go | identity | construction | measurement | correction
- **Bears on:** <probe (ours33 …), premise ids (L11 …), hole ids (G23 …), rulings (R-…)>
- **Status:** proved (test: `tests/…::…`) | numerical (preset or tool + seed) | conjecture
- **Check:** one number anyone can reproduce in one command, and the command.
- **Standard name:** what the mathematics is called, with one citation.
- **Does not show:** the nearest stronger claim that is *not* established.
- **Where:** `definition/theorems/T…md` §…, `problems/O….md` (date heading)
- **Absorb as (proposal):** which living-layer line would change and to what, in one sentence.
  Claude decides; the agent does not edit that line.
- **Retracts:** earlier card or theorem this one undercuts, or "none".
```

**Kinds, and what each means for the definition.** A *break* is a probe failing a premise (name both).
A *no-go* is an impossibility that holds for every function in a stated class. An *identity* is an exact
equality between two things the project treats as different. A *construction* is a function or law
built to meet named requirements; it is a benchmark unless a card says which fork it sorts and why that
sort is forced. A *measurement* is a number from a study, with its seed. A *correction* says an earlier
statement in the repo is wrong and gives both versions.

## What an agent may write

| Slot | Agent | Claude |
| --- | --- | --- |
| `findings/<branch>.md` | **writes** | reads first |
| new `definition/theorems/T*.md`, new sections in `problems/O*.md` (dated heading) | writes | absorbs |
| `met/`, `tests/`, `tools/`, `presets/` | writes | merges with tests |
| one line in `docs/STATE.md`'s log | writes | keeps |
| `definition/proof.md`, `equation.md`, `premises.md`, `ladder.md`, `type.md`, theorem index status | **does not** | writes from cards |
| `GOAL.md`, `TASK.md`, `CLAUDE.md`, `AGENTS.md`, STATE's main page, `rulings.md` | **does not** | writes |

`python tools/gate.py <branch>` lists any file a branch changes outside the agent's slots and checks
that `findings/<branch>.md` exists with every card's fields filled. A branch that fails the gate is
returned unread.

## Example (from `codex/ours32-consistency`, 3 Nov, rewritten into this format by Claude)

```
### F1: Literal L11 has no model: no estimator is consistent on every allowed history with divergent push
- **Kind:** no-go
- **Bears on:** L11, G25, O19; every probe
- **Status:** proved (test: `tests/test_signal_identifiability.py`)
- **Check:** r_i = b/sqrt(i), b = 1/2: TV(mixtures) <= 0.1624, so test error >= 0.4188 at every N.
  `python -m pytest tests/test_signal_identifiability.py -q`
- **Standard name:** Le Cam two-point method with mixtures (Polyanskiy & Wu); functional EIV consistency
  (Gleser 1981; Gallo 1982)
- **Does not show:** what the right growth condition is (T42/T53 give it for a Gaussian class only)
- **Where:** `definition/theorems/T41.md`; `problems/O19.md` (3 Nov)
- **Absorb as (proposal):** L11 restated as consistency under identifying information growth, as
  composition was restated on 28 Oct; grade unchanged.
- **Retracts:** none (T17's fixed-N scope already noted)
```
