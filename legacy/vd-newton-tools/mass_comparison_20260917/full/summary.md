# Refined mass-estimator comparison — simulation results

32,768 trials per signal cell; 46 cells; 11 point estimators. Identical observations are used for every method. All losses are smaller-is-better. These are results for the declared Gaussian 3D model.

Negative paired differences favor our ratio-of-means equation. ± values below are one empirical Monte Carlo standard error of a paired difference, not uncertainty in one object's mass. No global winner is calculated across the chosen cells.

## Original force-SNR-8 path

Mass changes along this path: true mass = 8 / acceleration SNR.

| Accel. SNR | True mass | Our squared log error | Norm ratio | Paired difference ± MCSE | Our factor-two success | Norm success |
|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 32 | 1.9404 | 3.3176 | -1.3772 ± 0.0089 | 8.94% | 3.31% |
| 0.5 | 16 | 0.7527 | 1.4331 | -0.6804 ± 0.0045 | 47.99% | 18.89% |
| 1 | 8 | 0.2882 | 0.4961 | -0.2079 ± 0.0020 | 79.22% | 59.07% |
| 2 | 4 | 0.1904 | 0.2052 | -0.0148 ± 0.0017 | 89.00% | 89.13% |
| 3 | 2.67 | 0.1377 | 0.1195 | 0.0181 ± 0.0008 | 93.64% | 95.89% |
| 5 | 1.6 | 0.0631 | 0.0552 | 0.0079 ± 0.0002 | 99.01% | 99.43% |
| 8 | 1 | 0.0333 | 0.0312 | 0.0022 ± 0.0000 | 99.96% | 99.97% |

## All equations around the earlier transition

Unbounded squared log error and bounded factor-ten-capped squared log error are displayed separately. The capped comparison includes every failed output.

| Estimator | Log², a=1 | Log², a=2 | Log², a=3 | Capped log², a=1 | Failed outputs, a=1 |
|---|---:|---:|---:|---:|---:|
| Our ratio of means | 0.2882 | 0.1904 | 0.1377 | 0.2882 | 0.00% |
| Flat-law geometric | 0.3386 | 0.3184 | 0.2052 | 0.3386 | 0.00% |
| Flat-law median | 0.3177 | 0.2281 | 0.1497 | 0.3177 | 0.00% |
| Flat-law root-loss point | 0.3808 | 0.3674 | 0.2271 | 0.3808 | 0.00% |
| Tube ratio of means | 0.2886 | 0.1937 | 0.1390 | 0.2886 | 0.00% |
| Vector-length ratio | 0.4961 | 0.2052 | 0.1195 | 0.4959 | 0.00% |
| Floored vector-length ratio | 0.4663 | 0.1873 | 0.1178 | 0.4663 | 0.00% |
| Positive TLS / profile | ∞¹ | ∞¹ | ∞¹ | 1.4628 | 16.58% |
| Forward OLS | ∞¹ | ∞¹ | ∞¹ | 2.2298 | 16.58% |
| Reverse OLS | ∞¹ | ∞¹ | ∞¹ | 1.4587 | 16.58% |
| Noise-corrected moment | ∞¹ | ∞¹ | ∞¹ | 3.3376 | 52.47% |

¹ Under the declared infinite penalty for invalid positive masses, population risk is infinite for these boundary/signed/unavailable-output rules. This remains true if a finite high-signal sample contains no invalid output. JSON retains the finite sample diagnostics and the separate population-risk annotation.

## Entire probability laws

CRPS evaluates log mass; 95% coverage is the fraction of intervals containing true mass. Coverage should be considered together with interval scores and widths in the full output.

| Accel. SNR | Flat log-CRPS | Tube log-CRPS | Flat 95% coverage | Tube 95% coverage | Flat mean log width | Tube mean log width |
|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 0.7772 | 0.7707 | 93.98% | 93.75% | 4.0259 | 4.0166 |
| 0.5 | 0.4817 | 0.4797 | 96.50% | 96.33% | 3.8709 | 3.8589 |
| 1 | 0.3477 | 0.3493 | 97.35% | 97.21% | 3.5113 | 3.4934 |
| 2 | 0.2776 | 0.2791 | 96.67% | 96.32% | 2.6666 | 2.6390 |
| 3 | 0.2158 | 0.2163 | 95.39% | 95.08% | 1.8801 | 1.8523 |
| 5 | 0.1403 | 0.1403 | 95.04% | 94.80% | 1.0537 | 1.0391 |
| 8 | 0.1026 | 0.1026 | 95.05% | 94.88% | 0.7348 | 0.7284 |

## Holding true mass fixed

Entries are our squared log error minus the norm ratio's. Negative favors our equation. Columns increase total signal strength √(f²+a²), with mass held fixed. For masses 16 and 32, this selected range still leaves acceleration weak; it does not map their eventual high-signal limit.

| True mass | Signal 0.5 | Signal 1 | Signal 2 | Signal 4 | Signal 8 | Signal 16 |
|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | -0.1844 ± 0.0067 | -0.1228 ± 0.0059 | -0.0520 ± 0.0039 | -0.1112 ± 0.0019 | -0.0163 ± 0.0016 | 0.0129 ± 0.0004 |
| 1 | -0.1949 ± 0.0040 | -0.1725 ± 0.0039 | -0.0848 ± 0.0033 | 0.0415 ± 0.0013 | 0.0094 ± 0.0001 | 0.0005 ± 0.0000 |
| 4 | -0.1700 ± 0.0066 | -0.1155 ± 0.0059 | -0.0417 ± 0.0037 | -0.1072 ± 0.0019 | -0.0190 ± 0.0017 | 0.0127 ± 0.0005 |
| 16 | -0.1434 ± 0.0116 | -0.0649 ± 0.0108 | -0.1224 ± 0.0089 | -0.8843 ± 0.0078 | -0.6819 ± 0.0045 | -0.2229 ± 0.0019 |
| 32 | -0.1317 ± 0.0143 | -0.0456 ± 0.0135 | -0.1847 ± 0.0121 | -1.4278 ± 0.0122 | -1.3777 ± 0.0090 | -0.7186 ± 0.0045 |

## Reproducibility and scope

All 46 per-cell refinement checks passed. Among the first 32 readings per cell, doubling integration order changed points by at most 1.8e-08 relatively and log interval endpoints by 5.96e-06 absolutely. These are numerical checks on that subset, not a uniform error bound on every possible reading.

The full JSON includes all point losses, quantiles, tolerance curves, paired differences, proper distribution scores, interval calibration and the two stated held-out prediction tasks. Empirical squared mass, relative and held-out residual loss MCSEs for the norm ratio (including their paired differences) are descriptive only: their population variance is infinite. Consult comparison_guide.md for equations, assumptions and primary sources.

Training seed 2026091601; independent validation seed 2026091702. Python 3.12.14, NumPy 2.3.5. Source hashes are in results.json. Managed context: `8782a293297330f7a354d87ebe62792018fa8866`. No real-data validation is claimed.
