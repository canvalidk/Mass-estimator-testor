# E12 · A partial proof for one scalar instrument

[Read the proof](explore/e12-scalar-proof.md). One scalar Gaussian instrument already admits an uncountable family of laws with exact measured-zero neutrality, sharp-threshold consistency, persistent spread on every fixed incompatible series and an asymptotically calibrated interval. Two members widen the project's 5-and-6 example by different amounts. The full definition is still open.

| Premise or question | Finding |
| --- | --- |
| O14, per-reading discount | An exactly fitted reading among three conflicting ones drives the interval to zero. The rule is point-consistent on every identifying array, yet its raw interval has limiting coverage zero on pure noise. [Proof](explore/e12-per-reading.md) |
| O14, L19, information about conflict | Two readings agreeing on 1 and two reading 4 and ¼ give identical pooled relative Gaussian evidence at every error scale. Individual relative mass curves recover the missing information. [Proof and sampling version](explore/e12-pooled-information.md) |
| L11, T30, scalar existence | A pooled construction attains the sharp identifying threshold for arbitrary signed pushes. An explicit extension treats axis and balanced-scatter conflicts. The evidence retained at finite noise is unselected. [Construction](explore/e12-pooled-witness.md) · [Extension](explore/e12-scalar-proof.md) |
| G28 and G6, report and point | A companion interval covers at least its stated level for every finite Gaussian experiment, and tends to that level on identifying series. Several reports and point readouts meet the tested requirements. [Proof](explore/e12-report.md) |
| T12 and G27, composition | Bounded final log contrast on every fixed conflict and multiplication of final factors prevent singleton law concentration. Same-object pooling from richer state and physical composition of distinct bodies remain separate obligations. [Proof](explore/e12-product-obstruction.md) |

Branch exploration from main `d2a3352`, 18 November 2026. Run `python -m pytest -q -k e12 tests` for the branch checks.
