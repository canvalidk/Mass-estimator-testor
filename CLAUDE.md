# Mass estimator tester

This repo holds the VD mass estimator, the function the project claims is the
definition of inertial mass, together with its proof and every test of it.

- `definition/`: the deliverable (type, premises, one file per theorem, the equation)
- `problems/`: one file per open item
- `rulings.md`: every ruling of Can's, dated, verbatim, append-only
- `met/`, `tests/`, `presets/`, `reference/`: the code, its tests, the studies
- `docs/`: archive, the records the living layer was built from
- `briefs/`: work queued for other agents (Codex), one sorted fork each; the protocol is its README
- `legacy/`: earlier programs, history, never run

**Read `TASK.md` first.** The task is to find the function which is unique: every fork
the readings cannot decide becomes an input (as a branch cut is an input to the complex
logarithm), and what survives, if it agrees with experiment, is the definition. At a
fork, sort it as forced, input, or wrong premise; do not pick by test score and do not
ask Can to choose a side.

**Bring your mathematics.** Before working a fork, say what it is called in standard mathematics or
statistics and what that theory already knows; the project's vocabulary is not a reason to
re-derive known results, and Can should not have to supply them (GOAL.md, standing rule 6).

**The equation the work is at is ours33** (`definition/equation.md`, `met/ours33.py`; R-2026-11-17-2). Read every
result, brief and premise against it. ours30 and the probes before it are history.

**The direction pursued is E13** (R-2026-11-22-1; `problems/O20.md`): harvest from the exploration branches into main's own
slots, each piece checked first and its branch named, until E13 has landed; nothing on the law-and-point line is struck in advance.

**The consequences of the inputs.** The definition computes what follows from what is declared, and a doubt about a
declared input is outside the problem (R-2026-11-17-1). Every σ is believed as given. How much less an improbable set of
readings is believed is computed by the definition from the readings, and handling such readings is the job
(R-2026-11-17-3).

**The fluke** (R-2026-11-24-4, R-2026-11-24-5). The test is calculated with the law held true, always. So some mass is
true, a candidate mass is compared with the other masses and with nothing else, and the χ² of that comparison reaches zero at
the best mass every time. The cut is fixed by how sure the statement is to be, before the readings are seen, and the kept set
is never empty. The curve's height says how large a fluke each mass requires. Its least (the leftover) is the fluke that
happened whichever mass is true: it tests nothing, because the law failing is not an alternative, and it is a measurement of
how misaligned the readings are. The definition computes that amount of fluke from the readings, finds the masses that need
the least fluke beyond it, and uses the amount to weight how reliable the result is. How the amount enters is what is being
looked for (`problems/O14.md`, 24 Nov).

**Then `docs/STATE.md`**, then run `python tools/audit.py`. Do not reconstruct
the state from the records; the living layer is the state.

Write only into the living layer: a theorem file, a problem file's attempts list, a
ruling appended to `rulings.md`, a test. Do not write a new record at the end of a
conversation.

## Rules that are easy to break

The full rules are in README.md, "How work gets into the tester". The ones that
matter most:

- Code lands as commits with tests. Every study is a preset with its seed.
- Never edit a record to agree with the code, or the code to agree with a record.
  A disagreement is a finding: add a `xfail(strict=True)` test that names both
  numbers, and leave both as they are.
- Rulings are Can's. Quote them verbatim with the date. Do not write your own
  suggestion as if it were decided; mark what is yours. A ruling is his deduction from the
  goal, provisional like every premise (R-2026-11-07-1): the premises are Claude's to deduce
  and restate, with the argument; what Can holds is the goal, and a sort that turns on what
  the goal means is a question to him.
- Only decisions go to VD-docs, through `outbox/VD-docs/`. Work stays in `docs/`.
- When you change what `docs/STATE.md` describes, update it in the same commit,
  and add a line to its log at the end of a session.
- **Runs by other agents (Codex and the like) get one sorted fork and a named deliverable:**
  a test, a counterexample, or a derivation, written into a named theorem or problem file.
  Not an open avenue to explore. Findings go into that file's attempts and one line of
  STATE.md's log; study narration does not go on STATE.md's main page. Such runs are queued
  as files in `briefs/`, branch from `main` as `codex/<brief id>`, and are merged by Claude as
  delivered (`briefs/README.md`, "How a brief moves"). Exploration branches (`explore/README.md`)
  are input: never merged; parts are harvested by need.
- A theorem is "proved" only when its file names the test that checks it. Mark
  what is yours as yours; "forced" means forced from a named premise.

## Commands

```sh
python -m pip install -r requirements.txt   # or requirements-lock.txt: the exact versions last run
python tools/audit.py                    # what is out of order in the living layer
python -m pytest -q                      # about 2.5 minutes; MET_SLOW=1 adds the slow ones
python -m met list                       # presets and the question each answers
python -m met run presets/ours29_noise_vs_signal.yaml --replicates 50   # --workers N|max: cells side by side
```

`docs/STATE.md` has a short example that calls the frontier probe, ours33, from Python.

The golden-output tests in `tests/test_reference_outputs.py` compare printed numbers,
so a different numpy or CPU can move a value that is only rounding noise. Printed
differences that are zero in exact arithmetic are checked against a bound (ZERO_BOUND), not
their digits; anything else that moves, check which value before treating it as real. Never
edit `reference/expected/` to make the test pass.
