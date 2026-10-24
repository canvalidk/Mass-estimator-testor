# Mass estimator tester

This repo runs every test and study of the VD mass estimator, the function the
project claims is the definition of inertial mass. The code is `met/`; the records
are `docs/`; earlier programs are in `legacy/` (history, never run).

**Start with `docs/STATE.md`.** It says what is current, what is open, and which
five files to read first. Do not reconstruct the state from the records in date
order; that is what the page is for.

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

## Commands

```sh
python -m pip install -r requirements.txt
python -m pytest -q                      # about 3.5 minutes; MET_SLOW=1 adds the slow ones
python -m met list                       # presets and the question each answers
python -m met run presets/ours29_noise_vs_signal.yaml --replicates 50
```

`docs/STATE.md` has a five-line example that calls ours29 from Python.

The golden-output tests in `tests/test_reference_outputs.py` compare printed numbers,
so a different numpy can move a value that is only rounding noise. Check
which value moved before treating a failure as real.
