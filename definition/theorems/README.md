# Theorems

One file per theorem. Each names the premises and theorems it rests on,
the test that checks it (or `none yet`), and its status. `[flat]` marks results that hold
only for the flat push; `[one pair]` marks results that need one noise pair.

| id | name | from | test |
| --- | --- | --- | --- |
| [T1](theorems/T01.md) | Exact branch | L1 | none yet |
| [T2](theorems/T02.md) | Candidates | L1 | none yet |
| [T3](theorems/T03.md) | Separation | L5, L14, L7 | none yet |
| [T4](theorems/T04.md) | Ratio of means | T3, L5, L14, T21, L7 | partial |
| [T5](theorems/T05.md) | Swap, readout half | T4 | yes |
| [T6](theorems/T06.md) | Push weight | L8 (conditional form), L9, L10, L7, **T7** (sweep E2) | partial |
| [T7](theorems/T07.md) | Prior once | L11 | none yet |
| [T8](theorems/T08.md) | The law (ours29) | L1, L13, T6, T7, L10, L21, T31 | yes |
| [T9](theorems/T09.md) | Stacked size | L12, L18 | yes |
| [T10](theorems/T10.md) | The readout | T4, T8, T9 | yes |
| [T11](theorems/T11.md) | Sufficiency | T8, T9, T6 | yes |
| [T12](theorems/T12.md) | Combination | T8 | yes |
| [T13](theorems/T13.md) | Totality | T8, T10, L2 | yes |
| [T14](theorems/T14.md) | Zero readings give $s$ | L2, L3′, L4, T5 | yes |
| [T15](theorems/T15.md) | Equal energies | T11, T5, L4 | none yet |
| [T16](theorems/T16.md) | A tilt is an unseen reading | T8 | none yet |
| [T17](theorems/T17.md) | Convergence | T8, T7 | none yet |
| [T18](theorems/T18.md) | Width | T8 [flat], L13 | partial |
| [T19](theorems/T19.md) | No mean | T8 | none yet |
| [T20](theorems/T20.md) | Mixed-instrument no-go | T8, T14 | none yet |
| [T21](theorems/T21.md) | C3 | L15, L16, L17, L18, composition | none yet |
| [T22](theorems/T22.md) | Full covariance | L16 | none in met/ (only reference/noise_vs_signal/ours29_reference_checks.py rotates a covariance) |
| [T23](theorems/T23.md) | Gaussian noise forced | L11, L20 | none yet |
| [T24](theorems/T24.md) | Z holds for the law | T8 | yes |
| [T25](theorems/T25.md) | Zero components carry nothing; $k$ does | T8 | yes |
| [T26](theorems/T26.md) | Limits of per-reading rules | L19, T8 | none yet |
| [T27](theorems/T27.md) | Own-world calibration | T8 (ours28) | none yet |
| [T28](theorems/T28.md) | Mass and push scale are orthogonal | T8 (ours28) | none yet |
| [T29](theorems/T29.md) | The readout fails Z | T10, L20 | yes |
| [T30](theorems/T30.md) | The tempered law's two ends | C_i | yes |
| [T31](theorems/T31.md) | The Bartlett power | L19 (scope), L20, L1, property N | yes |
| [T32](theorems/T32.md) | The misfit is the pairings' disagreement in noise units | T8, L21, T30, the 30 Sept eigenvalue dictionary | partial |
| [T33](theorems/T33.md) | Same window | the units build-up (units-and-axes record §6), L18, L19 | none yet |
