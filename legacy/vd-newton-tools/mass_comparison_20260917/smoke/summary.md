# Refined mass-estimator comparison — simulation results

512 trials per signal cell; 46 cells; 11 point estimators. Identical observations are used for every method. All losses are smaller-is-better. These are results for the declared Gaussian 3D model.

Negative paired differences favor our ratio-of-means equation. ± values below are one empirical Monte Carlo standard error of a paired difference, not uncertainty in one object's mass. No global winner is calculated across the chosen cells.

## Original force-SNR-8 path

Mass changes along this path: true mass = 8 / acceleration SNR.

| Accel. SNR | True mass | Our squared log error | Norm ratio | Paired difference ± MCSE | Our factor-two success | Norm success |
|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 32 | 1.8715 | 3.2621 | -1.3906 ± 0.0718 | 8.59% | 4.30% |
| 0.5 | 16 | 0.7124 | 1.3994 | -0.6871 ± 0.0365 | 51.56% | 19.92% |
| 1 | 8 | 0.2714 | 0.4820 | -0.2106 ± 0.0151 | 80.27% | 59.38% |
| 2 | 4 | 0.1868 | 0.2040 | -0.0173 ± 0.0117 | 90.43% | 88.87% |
| 3 | 2.67 | 0.1358 | 0.1145 | 0.0213 ± 0.0058 | 93.75% | 96.29% |
| 5 | 1.6 | 0.0621 | 0.0536 | 0.0085 ± 0.0013 | 99.61% | 99.80% |
| 8 | 1 | 0.0326 | 0.0306 | 0.0020 ± 0.0002 | 100.00% | 100.00% |

## All equations around the earlier transition

Unbounded squared log error and bounded factor-ten-capped squared log error are displayed separately. The capped comparison includes every failed output.

| Estimator | Log², a=1 | Log², a=2 | Log², a=3 | Capped log², a=1 | Failed outputs, a=1 |
|---|---:|---:|---:|---:|---:|
| Our ratio of means | 0.2714 | 0.1868 | 0.1358 | 0.2714 | 0.00% |
| Flat-law geometric | 0.3332 | 0.3217 | 0.2054 | 0.3332 | 0.00% |
| Flat-law median | 0.3032 | 0.2248 | 0.1475 | 0.3032 | 0.00% |
| Flat-law root-loss point | 0.3793 | 0.3729 | 0.2276 | 0.3793 | 0.00% |
| Tube ratio of means | 0.2721 | 0.1903 | 0.1372 | 0.2721 | 0.00% |
| Vector-length ratio | 0.4820 | 0.2040 | 0.1145 | 0.4820 | 0.00% |
| Floored vector-length ratio | 0.4505 | 0.1857 | 0.1145 | 0.4505 | 0.00% |
| Positive TLS / profile | ∞¹ | ∞¹ | ∞¹ | 1.5282 | 16.80% |
| Forward OLS | ∞¹ | ∞¹ | ∞¹ | 2.3025 | 16.80% |
| Reverse OLS | ∞¹ | ∞¹ | ∞¹ | 1.5278 | 16.80% |
| Noise-corrected moment | ∞¹ | ∞¹ | ∞¹ | 3.3769 | 53.12% |

¹ Under the declared infinite penalty for invalid positive masses, population risk is infinite for these boundary/signed/unavailable-output rules. This remains true if a finite high-signal sample contains no invalid output. JSON retains the finite sample diagnostics and the separate population-risk annotation.

## Entire probability laws

CRPS evaluates log mass; 95% coverage is the fraction of intervals containing true mass. Coverage should be considered together with interval scores and widths in the full output.

| Accel. SNR | Flat log-CRPS | Tube log-CRPS | Flat 95% coverage | Tube 95% coverage | Flat mean log width | Tube mean log width |
|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 0.7553 | 0.7489 | 94.34% | 94.34% | 4.0525 | 4.0437 |
| 0.5 | 0.4669 | 0.4651 | 97.46% | 97.46% | 3.9000 | 3.8885 |
| 1 | 0.3414 | 0.3432 | 97.85% | 97.85% | 3.5449 | 3.5277 |
| 2 | 0.2768 | 0.2784 | 98.24% | 98.05% | 2.7021 | 2.6751 |
| 3 | 0.2153 | 0.2158 | 95.51% | 95.12% | 1.9006 | 1.8729 |
| 5 | 0.1395 | 0.1394 | 94.92% | 94.73% | 1.0570 | 1.0425 |
| 8 | 0.1015 | 0.1014 | 94.92% | 94.73% | 0.7339 | 0.7275 |

## Holding true mass fixed

Entries are our squared log error minus the norm ratio's. Negative favors our equation. Columns increase total signal strength √(f²+a²), with mass held fixed. For masses 16 and32, this selected range still leaves acceleration weak; it does not map their eventual high-signal limit.

| True mass | Signal 0.5 | Signal 1 | Signal 2 | Signal 4 | Signal 8 | Signal 16 |
|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | -0.2087 ± 0.0519 | -0.1546 ± 0.0447 | -0.0568 ± 0.0286 | -0.1027 ± 0.0140 | -0.0199 ± 0.0105 | 0.0082 ± 0.0026 |
| 1 | -0.1471 ± 0.0270 | -0.1484 ± 0.0312 | -0.0572 ± 0.0218 | 0.0450 ± 0.0101 | 0.0090 ± 0.0009 | 0.0005 ± 0.0000 |
| 4 | -0.0930 ± 0.0475 | -0.0786 ± 0.0445 | -0.0193 ± 0.0276 | -0.1070 ± 0.0132 | -0.0208 ± 0.0118 | 0.0139 ± 0.0034 |
| 16 | -0.0484 ± 0.0872 | -0.0435 ± 0.0856 | -0.1215 ± 0.0705 | -0.8777 ± 0.0627 | -0.6886 ± 0.0366 | -0.2278 ± 0.0152 |
| 32 | -0.0278 ± 0.1082 | -0.0330 ± 0.1081 | -0.2022 ± 0.0970 | -1.4299 ± 0.0985 | -1.3911 ± 0.0718 | -0.7283 ± 0.0368 |

## Reproducibility and scope

All 46 per-cell refinement checks passed. Among the first 32 readings per cell, doubling integration order changed points by at most 1.8e-08 relatively and log interval endpoints by 5.96e-06 absolutely. These are numerical checks on that subset, not a uniform error bound on every possible reading.

The full JSON includes all point losses, quantiles, tolerance curves, paired differences, proper distribution scores, interval calibration and the two stated held-out prediction tasks. Empirical squared-mass-loss MCSEs for the norm ratio are descriptive only: their population variance is infinite. Consult comparison_guide.md for equations, assumptions and primary sources.

Training seed 2026091601; independent validation seed 2026091702. Python 3.12.14, NumPy 2.3.5. Source hashes are in results.json. Managed context: `8782a293297330f7a354d87ebe62792018fa8866`. No real-data validation is claimed.
