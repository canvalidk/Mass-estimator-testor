# reference/: implementations written beside the records

These are byte-for-byte copies of the reference implementations written beside
the mass-estimator records (the private VD-docs repository,
`Newton-analysis/_dscn_mass_estimator/`), carried into the tester on 30 September
so that the studies they back can be rerun from here. They are independent of
`met/`: nothing in `met/` imports them.

**They are never edited here.** The folder keeps the records' layout, so each
script's own imports work unchanged (`noise_vs_signal/ours29_reference_checks.py`
imports `../ours27/ours28_reference.py`). `SHA256SUMS` holds the hash of every
script as carried over, and `tests/test_reference_outputs.py` fails if any file
differs from it. `.gitattributes` stores this folder with no line-ending
conversion. A change to the math belongs in the records first, then a new copy
here.

## What is here

| file | record that defines it | run from its folder | time |
|---|---|---|---|
| `ours27/ours28_reference.py` | `ours27/ours28_reference_implementation_and_readout_2026-09-28.md` | (a module) | |
| `ours27/ours28_reference_checks.py` | the same, §2 and §4 (seed 7) | `python ours28_reference_checks.py` | 45 s |
| `ours27/ours28_reference_sim.py` | the same, §3 | `python ours28_reference_sim.py N omega series seed` | minutes to tens of minutes per row |
| `noise_vs_signal/ours29_reference.py` | `noise_vs_signal/ours29_misfit_weight_and_bartlett_power_2026-09-30.md` | (a module) | |
| `noise_vs_signal/ours29_reference_checks.py` | the same, §7 (seed 7) | `python ours29_reference_checks.py` | 30 s |
| `noise_vs_signal/ours29_reference_worlds.py` | the same, §8 (seed 30) | `python ours29_reference_worlds.py` | about 20 min |
| `noise_vs_signal/ours28b/ours28b.py` | `noise_vs_signal/ours28b/ours28b_one_scale_per_noise_pair_2026-09-28.md` (withdrawn the same night by the no-borrowing ruling; kept as record) | (a module) | |
| `noise_vs_signal/ours28b/ours28b_checks.py` | the same, §2 (seed 7) | `python ours28b_checks.py` | about 23 min |
| `mass_estimator_equation/mass_estimator.py` | `mass_estimator_equation/README.md` (the 21 September equation) | (a module) | |
| `mass_estimator_equation/example.py` | the same, its worked example | `python example.py` | under 1 s |

## How they are tested

- `tests/test_reference_outputs.py`: each checks script, and the ours29 worlds
  study, is run again and its printed output compared with `expected/`, the
  output it gave when carried over (slow ones with `MET_SLOW=1`). Text must match exactly; numbers to 1e-6
  relative (last-digit floating-point differences between machines, nothing
  more). The same test checks `SHA256SUMS`. This is repeatability, not
  correctness.
- `tests/test_reference_records_*.py`: that the code reproduces the numbers its
  records state, each assertion carrying the record's words verbatim, to the
  precision the record prints. Where the code does not reproduce a stated
  number, the test is kept and marked `xfail(strict=True)` with both numbers:
  the code and the record are both left as they are. Where a record's number
  cannot be regenerated from what was kept (a seed not recorded, a script not
  kept), the test file's docstring lists it, with the record's words.
- Slow tests (the ours29 worlds study, ours28b's checks, every ours28b call with
  the scale law) run only with `MET_SLOW=1`.

## Found when they were carried over (30 September)

- `ours28b_checks.py` is slow: 23 minutes on an idle machine (on 30 Sept, with
  other jobs sharing the machine, it had not finished after 50). Its output was
  identical in the tester's layout and the records' as far as both runs went.
- The ours29 record's §8 table is reproduced exactly by
  `ours29_reference_worlds.py` (all five worlds, every column).
- The ours28 record's §3 study table cannot be regenerated exactly:
  `ours28_reference_sim.py` takes the seed as an argument, and the table does
  not record the seed of any row.
- Numbers the code does not reproduce, and properties that do not hold as
  stated, are listed in the docstrings of the `tests/test_reference_records_*.py`
  files and marked `xfail(strict=True)` there.
