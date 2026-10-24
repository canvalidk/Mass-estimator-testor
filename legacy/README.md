# legacy/: earlier mass-estimator programs, kept as history

Nothing in here is run, maintained or imported by `met/`. The tests ignore this
folder (`pytest.ini` sets `testpaths = tests`). Code that is still wanted gets
rebuilt in `met/` with tests, under the rules in the top-level README; it is not
revived here.

## vd-newton-lab/

The mass-estimator laboratory from VD-Newton (`canvalidk/VD-Newton`),
`03_trace_and_evaluator/mass_estimator/`, as it stood at VD-Newton commit
`ecd1a1c` (2026-09-23), plus one script that was never committed
(`investigation/default_weighting_consistency_check.py`). 131 files, moved here
2026-10-24 and removed from VD-Newton, on his ruling that the estimator's work
lives in the tester and VD-Newton goes back to entry production. The tester
replaced this lab in September ("we can start over essentially with the tests").

- Paths inside it are the old paths with the prefix
  `03_trace_and_evaluator/mass_estimator/` dropped. Records and the
  reproducibility inventory (`docs/mass_estimator_reproducibility_inventory_2026-10-01.md`)
  cite the old prefix; the rest of the path still resolves here.
- `raw/robot/` and `raw/cart/` hold the raw robot and cart trial data. The
  robot data is the Zenodo dataset 10.5281/zenodo.11096791 (CC-BY 4.0,
  attribution in `raw/robot/record.json`). The cart data came from a course web
  page with no stated licence, so `raw/cart/` is gitignored and stays local, as
  does the paper `Gleser-Nonexistence1001-1987.pdf`; this repo is public.
  `SHA256SUMS` still lists them.
- Its `results/` folder is tracked here, despite the repo-wide `results/` rule:
  it is history, not output of this tester.
- `SHA256SUMS` holds the hash of every file as it left VD-Newton
  (`sha256sum -c SHA256SUMS` from inside the folder).
- Its own `README.md` is the lab's README, unchanged. Its links to VD-docs and
  to other VD-Newton paths are as they were.

## vd-newton-tools/

The scripts and notes from VD-Newton's gitignored `.tools/` folder that were not
already in the tester: 184 files from 28 folders, moved or copied 2026-10-24. The
study outputs they produced stayed in VD-Newton's `.tools/`. Its own README says
exactly what came, what stayed and why.
