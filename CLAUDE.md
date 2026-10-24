# Mass estimator tester

This repo holds the VD mass estimator, the function the project claims is the
definition of inertial mass, together with its proof and every test of it.

- `definition/`: the deliverable (type, premises, one file per theorem, the equation)
- `problems/`: one file per open item
- `rulings.md`: every ruling of Can's, dated, verbatim, append-only
- `met/`, `tests/`, `presets/`, `reference/`: the code, its tests, the studies
- `docs/`: archive, the records the living layer was built from
- `legacy/`: earlier programs, history, never run

**Read `TASK.md` first.** The task is to find the function which is unique: every fork
the readings cannot decide becomes an input (as a branch cut is an input to the complex
logarithm), and what survives, if it agrees with experiment, is the definition. At a
fork, sort it as forced, input, or wrong premise; do not pick by test score and do not
ask Can to choose a side.

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
  suggestion as if it were decided; mark what is yours.
- Only decisions go to VD-docs, through `outbox/VD-docs/`. Work stays in `docs/`.
- When you change what `docs/STATE.md` describes, update it in the same commit,
  and add a line to its log at the end of a session.
- A theorem is "proved" only when its file names the test that checks it. Mark
  what is yours as yours; "forced" means forced from a named premise.

## Commands

```sh
python -m pip install -r requirements.txt
python tools/audit.py                    # what is out of order in the living layer
python -m pytest -q                      # about 3.5 minutes; MET_SLOW=1 adds the slow ones
python -m met list                       # presets and the question each answers
python -m met run presets/ours29_noise_vs_signal.yaml --replicates 50
```

`docs/STATE.md` has a five-line example that calls ours29 from Python.

The golden-output tests in `tests/test_reference_outputs.py` compare printed numbers,
so a different numpy can move a value that is only rounding noise. Check
which value moved before treating a failure as real.
