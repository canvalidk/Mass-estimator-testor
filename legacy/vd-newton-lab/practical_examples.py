"""Regenerate the compact convergence and uncertainty examples."""
from dataclasses import asdict
from pathlib import Path
import json
import math
from estimator import estimate, isotropic_covariance
from practical_formulas import summarize, aligned_local_uncertainty


def main():
    rows = []
    for sigma in (1, .5, .25, .125):
        e = estimate([2, 0], [1, 0], isotropic_covariance(2, sigma, sigma),
                     direction_order=256, ratio_order=4096)
        rows.append(dict(sigma=sigma, **asdict(summarize(e))))
    cov = isotropic_covariance(2, .6, .2)
    strong = summarize(estimate([20, 0], [10, 0], cov,
                               direction_order=768, ratio_order=8192))
    local = aligned_local_uncertainty([20, 0], [10, 0], cov)
    null = summarize(estimate([0, 0], [0, 0], isotropic_covariance(2, 1, 1)))
    output = dict(
        model="Newton II assumed; Gaussian errors; flat positive magnitudes; uniform unknown 2D direction",
        uncertainty="Conditional equal-tail 95% interval and SD(log M); no finite exact SD(M)",
        aligned_force_2_acceleration_1=rows,
        resolved_example=dict(force=20, acceleration=10, sigma_force=.6,
                              sigma_acceleration=.2, **asdict(strong),
                              local_mass_spread=local, local_relative_spread=local/2,
                              local_log_interval=[2*math.exp(-1.959963984540054*local/2),
                                                  2*math.exp(1.959963984540054*local/2)]),
        zero_zero=asdict(null))
    target = Path(__file__).resolve().parent / 'results' / 'practical_examples.json'
    target.write_text(json.dumps(output, indent=2), encoding='utf-8')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
