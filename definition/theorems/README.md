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
| [T8](theorems/T08.md) | The law (ours29; ours30 since 29 Oct) | L1, L13, T6, T7, L10, L21, T31 | yes |
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
| [T34](theorems/T34.md) | Gain is the geometric mean of the pairings | L1, T1 | yes |
| [T35](theorems/T35.md) | The parametrisation-free point is the median | L23 (proposed), T5 | yes |
| [T36](theorems/T36.md) | The reading space, its metric, and one residual | the type, L3, L1, L22, L13 | yes |
| [T37](theorems/T37.md) | The per-reading noise size: what believing Σ and T30 force | L18, L13, L8, L21, T30 | yes |
| [T38](theorems/T38.md) | The supremum of the expected positive-mass misfit is k − ½ | T32, L13, L9, L10 | yes |
| [T39](T39.md) | Ours32: repeated-signal consistency and the global-power tail obstruction | L1, L11, L13 (conditional noise), L18, T31, T36, T37 | yes |
| [T40](T40.md) | Local curvature, tail preservation and experiment state | T31, T39, L11, L20, L3, L16; stated Gaussian-local benchmark | yes |
| [T41](T41.md) | Divergent total signal need not identify mass | L1, L2, L13, L18, unrestricted type/L19; L11 under test; L16 for corollary | yes |
| [T42](T42.md) | Sharp information growth for the flat Gaussian Bartlett law | T18, T31, T40, T41; stated Gaussian submodel and fixed-prior conditions | yes |
| [T43](T43.md) | Leverage-sensitive interval calibration and a Z-safe spectral control | T18, T31, T41, T42; stated Gaussian/prior conditions | yes |
| [T44](T44.md) | Factor averaging versus normalized whole-series cuts | T37, T39, provisional L21 cut; L19, L20; specified frozen noise module | yes |
| [T45](T45.md) | Exact Gaussian score sets and the proper-law coverage boundary | L1, L13, L18, full Gaussian type; T43; proper positive output law | yes |
| [T46](T46.md) | Positive-misfit quadratic weighting creates weak-signal bias | L13, L18; specified fixed weight class, T30, T42, T44 | yes |
| [T47](T47.md) | Normal-score CDF recalibration and its tail boundary | T39–T40, T43, T45; explicit law-only/naturality class and Gaussian benchmark | yes |
| [T48](T48.md) | Normalized Gaussian residual mixing and the latent-volume fork | T36–T37, T39; specified variance-independent declared-metric latent measure | yes |
| [T49](T49.md) | Graph geometry, scalar curve volume and measure distinctions | T36, T45, T48; canonical Grassmann metric, fixed declared covariance | yes |
| [T50](T50.md) | Gaussian operator graphs, TLS/Bingham profile and a vector score test | T45, T49; specified matrix experiment, full subspaces and selected tangents | yes |
| [T51](T51.md) | CDF recalibration, independent composition and joint covariance | T47; explicit additional universal interchange axiom; specified centered-Gaussian benchmark | yes |
| [T52](T52.md) | Infinite active information does not reveal an unknown signal allocation | T41; unrestricted fixed positive supports; proof-only sparse Gaussian mixtures | yes |
| [T53](T53.md) | Sharp uniform energy-class consistency in rank-one Gaussian estimation | T42, T52; specified exact-energy Gaussian nuisance classes | yes |
| [T54](T54.md) | Physical joint Gaussian law and calibrated small-noise composition | O2/G27, T45, T50; specified acquisition, supplied joint prior and Gaussian profile | yes |
| [T55](T55.md) | Joint profile sensitivity, sandwich covariance and scalar-inflation limits | T49, T51, T54; specified repeated-block regime and centered-Gaussian benchmark | yes |
| [T56](T56.md) | A conditional score law can have a globally inconsistent median | T45, T50, T54, T55; explicit exponential-score proposal and correlated Gaussian witness | yes |
| [T57](T57.md) | A vanishing profile anchor identifies a score law without changing its local width | T54–T56; additional oscillation-anchor family and stated repeated-block regime | yes |
| [T58](T58.md) | Fixed power anchors can lose coverage on identifying weak histories | T42, T43, T53, T57; fixed Gaussian weakening history and normalized-law proof | yes |
| [T59](T59.md) | Spectral half-power anchoring on all identifying scalar Gaussian histories | T43, T56, T58; common covariance, stated spectral power class and fixed regular prior | yes |
| [T60](T60.md) | Frame-safe joint trace anchoring on arbitrary identifying Gaussian histories | T54, T55, T59; common physical covariance, full shared-acceleration protocol and regular joint prior | yes |
| [T61](T61.md) | Joint Gaussian calibration with changing covariance shapes and arbitrary precisions | T54, T55, T60; uniformly comparable covariance shapes and stated information growth | yes |
| [T62](T62.md) | Posterior-only sequential closure versus additive carried state | T12, T56, T60; explicit additional sequential-closure condition | yes |
| [T63](T63.md) | Degenerating instrument shapes defeat a total-energy anchor | T56, T60, T61; fixed identifying Gaussian covariance history and a fixed proper prior | yes |
