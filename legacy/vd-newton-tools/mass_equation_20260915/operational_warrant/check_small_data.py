"""Independent radial/log-ratio quadrature for the small-data comparison."""
import json
import math
from pathlib import Path

import numpy as np


def gauss_interval(count, lower, upper):
    points, weights = np.polynomial.legendre.leggauss(count)
    return lower + (points + 1) * (upper - lower) / 2, weights * (upper - lower) / 2


def integrate(a, b, c, radial_count=160, log_count=512):
    radius, radius_weights = gauss_interval(radial_count, 0.0, 14.0)
    ell, ell_weights = gauss_interval(log_count, -32.0, 32.0)
    sin_theta = 1 / np.sqrt(1 + np.exp(-2 * ell))
    cos_theta = 1 / np.sqrt(1 + np.exp(2 * ell))
    sech = 1 / np.cosh(ell)
    h2 = a * sin_theta**2 + b * cos_theta**2 + c * sech
    hr = radius[None, :] * np.sqrt(np.maximum(h2[:, None], 0.0))
    kernel = np.ones_like(hr)
    positive = hr > 0
    kernel[positive] = np.sinh(hr[positive]) / hr[positive]
    angular_weights = ell_weights * sech / 2
    result = {}
    for name, power in (("flat", 0), ("tube", 2)):
        radial_weights = radius_weights * radius ** (power + 1) * np.exp(-radius**2 / 2)
        angle_density = kernel @ radial_weights
        first_radial_moment = kernel @ (radial_weights * radius)
        mass = angular_weights * angle_density
        normalizer = mass.sum()
        mean_log = float(mass @ ell / normalizer)
        variance_log = float(mass @ (ell - mean_log)**2 / normalizer)
        numerator = angular_weights @ (first_radial_moment * sin_theta)
        denominator = angular_weights @ (first_radial_moment * cos_theta)
        result[name] = {"mean_log": mean_log, "variance_log": variance_log,
                        "readout": float(numerator / denominator)}
    return result


def main():
    base = integrate(0, 0, 0)
    rows = []
    # Gram triples include positive and negative observed alignment.
    for a0, b0, c0 in ((4, 1, 1), (1, 1, 0.7), (1, 1, -0.7), (2, 2, 0)):
        for epsilon in (0.04, 0.02, 0.01):
            values = integrate(*(epsilon**2 * x for x in (a0, b0, c0)))
            for name, power in (("flat", 0), ("tube", 2)):
                coefficient = (power + 2) / 6
                expected = {"mean_log": coefficient * (a0 - b0) / 2,
                            "variance_log": -coefficient * math.pi * c0 / 3,
                            "readout": (power + 3) * (a0 - b0) / 18}
                observed = {key: (values[name][key] - base[name][key]) / epsilon**2
                            for key in expected}
                rows.append({"gram": [a0, b0, c0], "epsilon": epsilon, "measure": name,
                             "expected_quadratic_coefficients": expected,
                             "observed_scaled_changes": observed})
    output = {"method": "Independent Gauss-Legendre integration over radius and log ratio",
              "radial_domain": [0, 14], "log_ratio_domain": [-32, 32],
              "nodes": [160, 512], "null": base, "rows": rows}
    destination = Path(__file__).with_name("small_data_verification.json")
    destination.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    for row in rows:
        if row["epsilon"] == 0.01:
            print(row)


if __name__ == "__main__":
    main()
