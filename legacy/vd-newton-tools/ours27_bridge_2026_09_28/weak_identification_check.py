"""Weak-identification check, not a calibration target or candidate selection.

Nonzero R^2=.01 avoids assigning an identifiable mass to the exactly null
experiment. All cells share noise draws; instruments and laws are unchanged.
This tests finite central probability intervals, not Fieller confidence sets.
"""
from pathlib import Path
import json
import numpy as np
from local_information_probe import summary, cdf_at


def main():
    rng = np.random.default_rng(28092871)
    n, M, R2 = 10000, 63, .01
    ex = rng.normal(size=(n, M))
    ey = rng.normal(size=(n, M))
    v = np.full(M, np.sqrt(R2/M))
    rows = []
    for ratio in [.001, .01, .1, .25, 1., 4., 10., 100., 1000.]:
        theta = np.arctan(ratio)
        s = summary(np.sin(theta)*v+ex, np.cos(theta)*v+ey)
        row = {'mass_over_instrument_scale': ratio}
        for adjusted, label in [(False, 'ours27'), (True, 'local_information_diagnostic')]:
            cdf = cdf_at(s, theta, adjusted)
            refined = cdf_at(tuple(a[:80] for a in s), theta, adjusted, nodes=160)
            err = float(np.max(np.abs(refined-cdf[:80])))
            assert err < 2e-9, err
            row[label] = {'central_95_interval_coverage': float(np.mean((cdf >= .025)&(cdf <= .975))),
                          'cdf_refinement_max_error': err}
        rows.append(row)
    out = {'purpose': 'Demonstrate weak-signal coverage failure of bounded central intervals; not a candidate ranking.',
           'source_commit': 'b4bfc2040031ef010416334d0c5f4aa2ad05cbd8',
           'trials_per_cell': n, 'components': M, 'total_standardized_signal_squared': R2,
           'regular_information_expression': R2**2/(R2+M),
           'angle_prior': 'uniform on (0,pi/2), equivalently half-Cauchy at instrument scale',
           'rows': rows}
    Path(__file__).with_name('weak_identification_results.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
