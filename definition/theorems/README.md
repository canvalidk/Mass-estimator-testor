# Theorems

One file per theorem. Each names the premises and theorems it rests on,
the test that checks it (or `none yet`), and its status. `[flat]` marks results that hold
only for the flat push; `[one pair]` marks results that need one noise pair.

| id | name | from | test | scope |
| --- | --- | --- | --- | --- |
| [T1](theorems/T01.md) | Exact branch | L1 | none yet | |
| [T2](theorems/T02.md) | Candidates | L1 | none yet | |
| [T3](theorems/T03.md) | Separation | L5, L14, L7 | none yet | comparison readout only, since 28 Oct |
| [T4](theorems/T04.md) | Ratio of means | T3, L5, L14, T21, L7 | partial | comparison readout only, since 28 Oct |
| [T5](theorems/T05.md) | Swap, readout half | T4 | yes | |
| [T6](theorems/T06.md) | Push weight | L8 (conditional form), L9, L10, L7, **T7** (sweep E2) | partial | |
| [T7](theorems/T07.md) | Prior once | L11 | none yet | |
| [T8](theorems/T08.md) | The law (ours29; ours30 since 29 Oct) | L1, L13, T6, T7, L10, L21, T31 | yes | |
| [T9](theorems/T09.md) | Stacked size | L12, L18 | yes | comparison readout only, since 28 Oct |
| [T10](theorems/T10.md) | The readout | T4, T8, T9 | yes | comparison readout only, since 28 Oct |
| [T11](theorems/T11.md) | Sufficiency | T8, T9, T6 | yes | comparison readout only, since 28 Oct |
| [T12](theorems/T12.md) | Combination | T8 | yes | |
| [T13](theorems/T13.md) | Totality | T8, T10, L2 | yes | |
| [T14](theorems/T14.md) | Zero readings give $s$ | L2, L3′, L4, T5 | yes | comparison readout only, since 28 Oct |
| [T15](theorems/T15.md) | Equal energies | T11, T5, L4 | none yet | comparison readout only, since 28 Oct |
| [T16](theorems/T16.md) | A tilt is an unseen reading | T8 | none yet | |
| [T17](theorems/T17.md) | Convergence | T8, T7 | none yet | |
| [T18](theorems/T18.md) | Width | T8 [flat], L13 | partial | |
| [T19](theorems/T19.md) | No mean | T8 | none yet | |
| [T20](theorems/T20.md) | Mixed-instrument no-go | T8, T14 | none yet | |
| [T21](theorems/T21.md) | C3 | L15, L16, L17, L18, composition | none yet | comparison readout only, since 28 Oct |
| [T22](theorems/T22.md) | Full covariance | L16 | none in met/ (only reference/noise_vs_signal/ours29_reference_checks.py rotates a covariance) | |
| [T23](theorems/T23.md) | Gaussian noise forced | L11, L20 | none yet | |
| [T24](theorems/T24.md) | Z holds for the law | T8 | yes | |
| [T25](theorems/T25.md) | Zero components carry nothing; $k$ does | T8 | yes | |
| [T26](theorems/T26.md) | Limits of per-reading rules | L19, T8 | none yet | |
| [T27](theorems/T27.md) | Own-world calibration | T8 (ours28) | none yet | |
| [T28](theorems/T28.md) | Mass and push scale are orthogonal | T8 (ours28) | none yet | |
| [T29](theorems/T29.md) | The readout fails Z | T10, L20 | yes | comparison readout only, since 28 Oct |
| [T30](theorems/T30.md) | The tempered law's two ends | C_i | yes | |
| [T31](theorems/T31.md) | The Bartlett power | L19 (scope), L20, L1, property N | yes | |
| [T32](theorems/T32.md) | The misfit is the pairings' disagreement in noise units | T8, L21, T30, the 30 Sept eigenvalue dictionary | partial | |
| [T33](theorems/T33.md) | Same window | the units build-up (units-and-axes record §6), L18, L19 | none yet | |
| [T34](theorems/T34.md) | Gain is the geometric mean of the pairings | L1, T1 | yes | |
| [T35](theorems/T35.md) | The parametrisation-free point is the median | L23 (proposed), T5 | yes | |
| [T36](theorems/T36.md) | The reading space, its metric, and one residual | the type, L3, L1, L22, L13 | yes | |
| [T37](theorems/T37.md) | The per-reading noise size: what believing Σ and T30 force | L18, L13, L8, L21, T30 | yes | |
| [T38](theorems/T38.md) | The supremum of the expected positive-mass misfit is k − ½ | T32, L13, L9, L10 | yes | |
| [T39](T39.md) | Ours32: repeated-signal consistency and the global-power tail obstruction | L1, L11, L13 (conditional noise), L18, T31, T36, T37 | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T40](T40.md) | Local curvature, tail preservation and experiment state | T31, T39, L11, L20, L3, L16; stated Gaussian-local benchmark | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T41](T41.md) | Divergent total signal need not identify mass | L1, L2, L13, L18, unrestricted type/L19; L11 under test; L16 for corollary | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T42](T42.md) | Sharp information growth for the flat Gaussian Bartlett law | T18, T31, T40, T41; stated Gaussian submodel and fixed-prior conditions | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T43](T43.md) | Leverage-sensitive interval calibration and a Z-safe spectral control | T18, T31, T41, T42; stated Gaussian/prior conditions | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T44](T44.md) | Factor averaging versus normalized whole-series cuts | T37, T39, provisional L21 cut; L19, L20; specified frozen noise module | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T45](T45.md) | Exact Gaussian score sets and the proper-law coverage boundary | L1, L13, L18, full Gaussian type; T43; proper positive output law | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T46](T46.md) | Positive-misfit quadratic weighting creates weak-signal bias | L13, L18; specified fixed weight class, T30, T42, T44 | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T47](T47.md) | Normal-score CDF recalibration and its tail boundary | T39–T40, T43, T45; explicit law-only/naturality class and Gaussian benchmark | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T48](T48.md) | Normalized Gaussian residual mixing and the latent-volume fork | T36–T37, T39; specified variance-independent declared-metric latent measure | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T49](T49.md) | Graph geometry, scalar curve volume and measure distinctions | T36, T45, T48; canonical Grassmann metric, fixed declared covariance | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T50](T50.md) | Gaussian operator graphs, TLS/Bingham profile and a vector score test | T45, T49; specified matrix experiment, full subspaces and selected tangents | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T51](T51.md) | CDF recalibration, independent composition and joint covariance | T47; explicit additional universal interchange axiom; specified centered-Gaussian benchmark | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T52](T52.md) | Infinite active information does not reveal an unknown signal allocation | T41; unrestricted fixed positive supports; proof-only sparse Gaussian mixtures | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T53](T53.md) | Sharp uniform energy-class consistency in rank-one Gaussian estimation | T42, T52; specified exact-energy Gaussian nuisance classes | yes | harvested 7 Nov from `codex/ours32-consistency`; not checked by Claude |
| [T69](T69.md) | The maximal-invariant likelihood and its limits | L13, T36; one-instrument Gaussian experiment (Moreira's reduction to a noncentral Wishart) | yes | B21, merged 17 Nov as delivered; not checked by Claude |
| [T70](T70.md) | The reading as a qubit; misalignment is imaginarity | T36, L9, L13, T30, L20; exact identities on 1D–3D readings, a no-go in the limit | yes | Claude, 7 Nov; row added 17 Nov (missed then) |
| [T71](T71.md) | What makes inertia a real scalar | L16, L22 (provisional), T36, T70; a linear inertia commuting with the stated symmetries | yes | B23, merged 17 Nov as delivered; not checked by Claude |
| [T72](T72.md) | Fisher consistency does not force sign-dependent reading statistics to cancel | L13, L18, L3, L16, T36, T39, T46; independent Gaussian scalar slice; a counterexample, not an adopted factor | yes | B25, merged 17 Nov as delivered; checked by Claude 17 Nov (O18) |
| [T73](T73.md) | The evidence curve with no starting law | T36, L13, L18, L2, L20, T24, T12 | yes | harvested 22 Nov from `explore/e13-unnormalized-evidence` (Can); checked by Claude |
| [T74](T74.md) | Evidence curves of two bodies compose with no starting law | T73, L5, L14, T5, T12 | yes | harvested 22 Nov from `explore/e13-unnormalized-evidence` (Can); checked by Claude; (6), the composite with a kept promise, from its insight 7 |
| [T75](T75.md) | Exact conditional tests of a candidate mass on the full type; the scale-free null; the two joined | T36, T45, T73, L13, L18, L20, L19, L10 | yes | harvested 22 Nov from `explore/e13-unnormalized-evidence` (Can); reimplemented and checked by Claude; a research reference, not in `equation.md` |
| [T76](T76.md) | The weighted conditional scores: one exact family, two closed forms, no best member | T45, T75, L19, L10, L20, L11 | yes | harvested 22 Nov from `explore/e13-unnormalized-evidence` (Can; insight 9); checked by Claude |
| [T77](T77.md) | Noise-only readings beside a fixed result: what cannot be promised, and what must be carried | L11, L13, L18, T45, T76 | yes | harvested 22 Nov from `explore/e2-null-combination` (Codex); checked by Claude, the fixed-history form read only |
| [T78](T78.md) | On one scalar instrument four requirements have a model, and a family of them | L20, L11, T30, T45, T42, T43, L3, L4 | yes | harvested 22 Nov from `explore/e12-scalar-mass-proof` (Codex; its scalar proof); rewritten and checked by Claude, the limit theorem read and simulated; a witness, not an equation |
| [T79](T79.md) | T45's signed ranks are not one CDF: how much of a law the readings alone supply | T45, T73, L13, L18 | yes | harvested 22 Nov from `explore/e12-scalar-mass-proof` (Codex; the conditional-ranks note, first two sections); checked by Claude |
| [T80](T80.md) | A true statement about mass must be able to keep every large mass: what forces the kind of output | L2, L10, property N (R-2026-09-29-4), T45, T73, T75 | yes | Claude, 22 Nov, on R-2026-11-22-2; (1) is T45 §4 with its proof checked; the sort is provisional; sorts G28 |
| [T81](T81.md) | The point read off the curve: the least composes, composition at the point forces it, and it can sit at zero or infinite mass | T73, T74, L5, L14, L3, L11 | yes | harvested 23 Nov from `explore/e14-ours34-walls` (Antigravity with Can; its E14-1 and E14-2); proofs read and numbers recomputed by Claude; not adopted in `equation.md` |
| [T82](T82.md) | The discount under one unknown common scale: forced there and only there; what it costs where the declared σ is right | T45, T73, T76, T77, T81; R-2026-11-17-3 | yes | harvested 24 Nov from `explore/e14-ours34-walls` (its E14-5, E14-6, E14-7), checked by Claude on main's code; (5)–(8) and the sorts are Claude's, provisional |
