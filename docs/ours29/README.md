# ours29/ — the ours29 equation

**2026-09-30 (afternoon).** Written by Claude on the user's instruction: "make
it outside it. And it can stay mostly empty. We can keep the ours29 folder
exclusive to it for now and keep working in noise."

## What this folder is

The folder for ours29, the working equation since 30 September: ours27's
flat factor per reading, weighted by the reading's own misfit, with the
product raised to the Bartlett power $\eta$. Like `ours26/` and `ours27/`, it
is a grouping inside `_dscn_mass_estimator`, not a nested discussion.

It is kept **exclusive to the equation itself** and mostly empty for now.
The investigation continues in `../noise_vs_signal/`.

## Where ours29 is written down

- `../noise_vs_signal/ours29_misfit_weight_and_bartlett_power_2026-09-30.md`:
  the record that defines ours29, with the rulings of 30 Sept (11:57), the
  checks and the tester study.
- `../noise_vs_signal/ours29_reference.py` (with `_checks.py` and
  `_worlds.py`): the reference implementation.
- In the tester (`canvalidk/Mass-estimator-testor`, commit `645dfc7`): `misfit` and
  `calibrate: bartlett`, `presets/ours29_noise_vs_signal.yaml`.

Nothing has been moved here. Those files stay in `noise_vs_signal/`.
