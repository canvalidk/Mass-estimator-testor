r"""The 21 September equation: what its records state, checked against the kept code.

Records (read-only, in the records folder Newton-analysis/_dscn_mass_estimator/):
  mass_estimator_equation/README.md  (the 21 September specification and its example)
  axis_splitting_and_the_24_september_change_2026-09-25.md  (section 10, the claim that
      the tester's flat reference agrees with the standalone script; and the 21 Sept
      column of the section 4(a) median table)
Code:
  reference/mass_estimator_equation/mass_estimator.py  (byte-for-byte copy; example.py
      is not imported because it prints, its inputs are rebuilt here line for line)
  met.analyst.estimate  (the tester), law {reference: flat, nuisance: radius}, rule single

(a) mass_estimator.py against the README's worked example (point, interval, log SD,
    last refinement change) and the other numbers the README states for the function:
    its default settings, the inverse-mass outputs, and the measured-zeros special case.
    Not a duplicate of an existing test: tests/test_readouts.py::test_equation_readme_example
    checks the tester's own implementation (met), not mass_estimator.py, and at rel 2e-4
    for the interval and log SD; tests/test_reference_outputs.py compares example.py's
    printout with the output captured from the script itself (repeatability, not the
    README's numbers).

(a') Properties the README states without printing a number, tested on mass_estimator.py:
    invalid inputs raise ValueError and unresolved quadrature raises ArithmeticError;
    the direction order is unused in the analytic 3D, known-direction and 1D
    calculations; the equal 1/2 weights of the two 1D directions; the half-Cauchy law
    at measured zeros (its quantiles at several interval contents); unit changes;
    rotations; channel reciprocity; the known-direction truncated-normal means.
    The README names the last four among its verification checks but not their content
    or inputs. Their content is read here from the equation, and each such assertion
    carries the note (consequence: ...). Inputs: the README example where it applies, the
    section 4(a) readings of the 25 September record, and otherwise inputs chosen here
    (named as such in each test).
    Tolerances where the record prints no number: exact equality where the two sides
    run the same arithmetic; 1e-12 where only rounding separates them; and the README's
    own numerical tolerance, rtol 2e-4, where the code's quadrature is compared with a
    closed form.
    Observed here, not asserted (the README calls refinement "stability check, not a
    certified error bound"): the README example written in 2D, F=(10, 0), a=(2, 0.1),
    same sigmas, gives mass 4.999529133519121 at default settings (numerical_change
    1.43e-4) against 4.993753906248473 with direction_order=200 (the tester gives
    4.993753906248474), a relative difference of 1.2e-3.

(a'') The 21 Sept column of the 25 September record's section 4(a) median table, first
    three rows (inputs stated there): the tester's flat single-reading median, and the
    standalone script through its central interval at vanishing content.

(b) The cross-implementation claim, section 10 of the 25 September record:
    the tester's flat single-reading ratio-of-means point against mass_estimator.py's
    `mass`, in 1D, 2D and 3D. The record does not fix the inputs, the settings, the
    tester commit, or whether "to $2\times10^{-10}$" is relative or absolute. So:
      - the inputs are chosen here: seeded random readings, default_rng(100 + d), five
        per dimension, drawn in the order F = 3 N(0, 1)^d, a = N(0, 1)^d,
        (sigma_F, sigma_a) = U(0.3, 2)^2, with covariance diag(sigma_F^2 I_d, sigma_a^2 I_d);
        plus the README's own 3D example. This is the first input set tried; it was not
        changed after seeing the results;
      - both implementations run at their default settings;
      - the two readings of "to $2\times10^{-10}$" are tested separately (relative and
        absolute; the absolute reading depends on the units, which the record does not fix).
    Where the code does not give the stated agreement the case is xfail(strict=True)
    with the numbers it gives. Measured here, not asserted: with mass_estimator.py
    started at ratio_order=8192 instead of its default 512, it agrees with the tester to
    at most 1.1e-14 relative (8.6e-14 absolute) on all sixteen cases, and refining the
    tester's grid (grid_points 8001) moves the tester by at most 9e-16 relative.

Not reproducible from what was kept:
  - the inputs, settings and tester commit behind the 2e-10 statement. The record
    (axis_splitting_and_the_24_september_change_2026-09-25.md, "## 10. Reproduction")
    gives only "The tester's `flat` reference agrees with the standalone" /
    "`mass_estimator_equation/mass_estimator.py` to $2\times10^{-10}$ in 1D, 2D" /
    "and 3D." (three source lines). (b) tests the statement on inputs chosen here instead.
  - the README's nine verification checks as they were run: mass_estimator_equation/README.md,
    "## Version check and provenance", says "Verification scripts and rendering intermediates stay"
    / "in workspace scratch space, outside this presentation folder." Those scripts
    were not kept with the records. Some of the named properties are tested in (a')
    on inputs chosen here. Not tested here, from the lines
      "correlated-null quantiles and log spread, known-direction truncated-normal"
      "means, agreement with the latest isotropic backend and the original full"
      "covariance implementation, an independent integral in the original"
      "magnitude/direction coordinates, rotations, unit changes, channel reciprocity,"
      "anisotropic 3D inputs, invalid inputs and explicit numerical nonconvergence."
    : the correlated-null quantiles and log spread (no values stated), the agreement
    with the isotropic backend and the full-covariance implementation (neither is kept),
    the independent integral, and the anisotropic 3D inputs (no inputs or values stated).
  - section 4(a) of the 25 September record, the row
    "| the same three axes as three 1D readings | 3.46 | 3.72 |": the combined rule is
    described in "## 1. The premise" as
      "each with its latent pair integrated under the radius split, and the mass"
      "factor counted once (contribution 12; the tester's `likelihood_product`,"
    and the tester's likelihood_product needs an explicit prior factor. The record does
    not give that factor or its center for the 21 September law, so the row is not
    tested. The 24 Sept column of that table is the 24 September law, outside this
    family, and is not checked here.
"""

import dataclasses
import inspect
import math
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import norm

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "reference" / "mass_estimator_equation"))

from mass_estimator import estimate_mass  # noqa: E402

from met.analyst import estimate  # noqa: E402

FLAT = {"reference": "flat", "nuisance": "radius"}
SINGLE = {"rule": "single", "prior": []}


def printed(decimals):
    """Tolerance for a number printed with `decimals` decimals."""
    return 0.5 * 10.0 ** -decimals + 1e-12


def _example():
    """example.py's inputs, as it builds them."""
    force = [10.0, 0.0, 0.0]         # N: measured net force vector
    acceleration = [2.0, 0.1, 0.0]   # m/s^2: measured acceleration vector
    sigma_force = 0.5               # N: per-coordinate standard uncertainty
    sigma_acceleration = 0.1        # m/s^2: per-coordinate standard uncertainty
    covariance = np.diag([sigma_force**2] * 3 + [sigma_acceleration**2] * 3)
    return force, acceleration, covariance


# ---------------------------------------------------------------- (a) the README

def test_worked_example():
    force, acceleration, covariance = _example()
    result = estimate_mass(force, acceleration, covariance)

    # README.md, "## How to use the Python function":
    #   "For this example, the rounded outputs are **4.993754 kg**, a **95% conditional"
    assert result.mass == pytest.approx(4.993754, abs=printed(6))
    # README.md, "## How to use the Python function":
    #   "For this example, the rounded outputs are **4.993754 kg**, a **95% conditional"
    assert result.probability == 0.95
    # README.md, "## How to use the Python function":
    #   "interval [4.344146, 5.740363] kg**, and **log SD 0.071035**. At the default"
    assert result.interval[0] == pytest.approx(4.344146, abs=printed(6))
    # README.md, "## How to use the Python function":
    #   "interval [4.344146, 5.740363] kg**, and **log SD 0.071035**. At the default"
    assert result.interval[1] == pytest.approx(5.740363, abs=printed(6))
    # README.md, "## How to use the Python function":
    #   "interval [4.344146, 5.740363] kg**, and **log SD 0.071035**. At the default"
    assert result.log_standard_deviation == pytest.approx(0.071035, abs=printed(6))
    # README.md, "## How to use the Python function":
    #   "settings, the final scaled numerical refinement change is about `6.61e-5`."
    #   (6.61e-5 = 0.0000661, printed to 7 decimals)
    assert result.numerical_change == pytest.approx(6.61e-5, abs=printed(7))


def test_worked_example_refinement_bookkeeping():
    """Read from the README's text (it prints no final resolution): the final ratio order
    is the start doubled at least once and at most max_refinements times, and the last
    check passed at rtol."""
    result = estimate_mass(*_example())
    # README.md, "## How to use the Python function":
    #   "Every successful call compares at least two grid resolutions. Starting"
    #   "`direction_order=24` and `ratio_order=512` are doubled as needed, up to"
    #   "`max_refinements=3`; `rtol=2e-4` applies to all returned numerical summaries."
    #   "`numerical_change`, `direction_order` and `ratio_order` report the last check"
    #   "and final resolution. The direction order is unused in the analytic 3D,"
    #   (read from the quoted text, not a printed number: start 512, doubled, at least two
    #   resolutions compared, so at least one doubling; at most 3 doublings)
    assert result.ratio_order in (512 * 2, 512 * 4, 512 * 8)
    # README.md, "## How to use the Python function":
    #   "`max_refinements=3`; `rtol=2e-4` applies to all returned numerical summaries."
    #   "`numerical_change`, `direction_order` and `ratio_order` report the last check"
    #   (read from the quoted text, not a printed number: the reported last check is
    #   within rtol)
    assert result.numerical_change <= 2e-4


def test_default_settings():
    defaults = {k: p.default for k, p in inspect.signature(estimate_mass).parameters.items()
                if p.default is not inspect.Parameter.empty}
    # README.md, "## Types and inputs":
    #   "| `probability` | Probability content of the central equal-tail interval; default `0.95`. |"
    assert defaults["probability"] == 0.95
    # README.md, "## How to use the Python function":
    #   "`direction_order=24` and `ratio_order=512` are doubled as needed, up to"
    assert defaults["direction_order"] == 24
    # README.md, "## How to use the Python function":
    #   "`direction_order=24` and `ratio_order=512` are doubled as needed, up to"
    assert defaults["ratio_order"] == 512
    # README.md, "## How to use the Python function":
    #   "`max_refinements=3`; `rtol=2e-4` applies to all returned numerical summaries."
    assert defaults["max_refinements"] == 3
    # README.md, "## How to use the Python function":
    #   "`max_refinements=3`; `rtol=2e-4` applies to all returned numerical summaries."
    assert defaults["rtol"] == 2e-4


def test_inverse_mass_outputs():
    result = estimate_mass(*_example())
    # README.md, "## Uncertainty: exactly what is returned":
    #   "$\widehat\mu=1/\widehat m$ and $I_\mu=[1/I_c^{\rm upper},1/I_c^{\rm lower}]$."
    #   "These are `result.inverse_mass` and `result.inverse_mass_interval`."
    assert result.inverse_mass == 1.0 / result.mass
    # README.md, "## Uncertainty: exactly what is returned":
    #   "$\widehat\mu=1/\widehat m$ and $I_\mu=[1/I_c^{\rm upper},1/I_c^{\rm lower}]$."
    #   "These are `result.inverse_mass` and `result.inverse_mass_interval`."
    assert result.inverse_mass_interval == (1.0 / result.interval[1], 1.0 / result.interval[0])


@pytest.mark.parametrize("d", [1, 2, 3])
@pytest.mark.parametrize("sf,sa", [(1.0, 1.0), (2.0, 0.5), (0.3, 3.0)])
def test_measured_zeros(d, sf, sa):
    """Two measured zero vectors, independent isotropic errors (sigmas chosen here: the
    README states the result for any). Exact values (sigma_F/sigma_a, pi/2) are asserted
    to 1e-12; the printed interval to its 6 decimals."""
    s = sf / sa
    result = estimate_mass(np.zeros(d), np.zeros(d), np.diag([sf**2] * d + [sa**2] * d))
    # README.md, "## Scope and special cases":
    #   "  measured zero vectors and independent isotropic errors, $\widehat m="
    #   "  \sigma_F/\sigma_a$; $M/(\sigma_F/\sigma_a)$ is half-Cauchy. Its central"
    assert result.mass == pytest.approx(s, abs=1e-12)
    # README.md, "## Scope and special cases":
    #   "  95% interval is approximately `[0.039290, 25.451700]` times that scale and"
    assert result.interval[0] / s == pytest.approx(0.039290, abs=printed(6))
    # README.md, "## Scope and special cases":
    #   "  95% interval is approximately `[0.039290, 25.451700]` times that scale and"
    assert result.interval[1] / s == pytest.approx(25.451700, abs=printed(6))
    # README.md, "## Scope and special cases":
    #   "  its log SD is $\pi/2$. This apparatus-dependent output does not establish"
    assert result.log_standard_deviation == pytest.approx(math.pi / 2, abs=1e-12)


# ---------------------------------------------------------------- (a') README properties without printed numbers

INVALID = {
    # README.md, "## Types and inputs":
    #   "| $\widetilde{\mathbf F}$ / `force` | Finite measured **net force vector**, shape `(d,)`; SI unit N. Components may have either sign or be zero. |"
    "force not finite": {"force": [math.nan, 0.0, 0.0]},
    # README.md, "## Types and inputs":
    #   "| $\widetilde{\mathbf a}$ / `acceleration` | Finite measured acceleration vector, shape `(d,)`, in the same axes; SI unit m/s². |"
    "acceleration not the force's shape": {"acceleration": [2.0, 0.1]},
    # README.md, "## Types and inputs":
    #   "\longrightarrow\mathrm{MassEstimate},\qquad d\in\{1,2,3\}."
    "d = 4": {"force": [1.0, 2.0, 3.0, 4.0], "acceleration": [1.0, 2.0, 3.0, 4.0],
              "covariance": np.eye(8)},
    # README.md, "## Types and inputs":
    #   "| $\Sigma$ / `covariance` | Symmetric positive-definite matrix, shape `(2*d, 2*d)`, for the measurement errors in `[F1,...,Fd,a1,...,ad]` order. |"
    "covariance wrong shape": {"covariance": np.diag([0.25] * 3 + [0.01] * 2)},
    # README.md, "## Types and inputs": (the same row as above)
    #   "| $\Sigma$ / `covariance` | Symmetric positive-definite matrix, shape `(2*d, 2*d)`, for the measurement errors in `[F1,...,Fd,a1,...,ad]` order. |"
    "covariance not symmetric": {"covariance": np.diag([0.25] * 3 + [0.01] * 3)
                                 + np.triu(np.full((6, 6), 1e-3), 1)},
    # README.md, "## Scope and special cases":
    #   "  positive-mass Newton II. Zero covariance is not an empirical input."
    "covariance zero": {"covariance": np.zeros((6, 6))},
    # README.md, "## Types and inputs":
    #   "| $\Sigma$ / `covariance` | Symmetric positive-definite matrix, shape `(2*d, 2*d)`, for the measurement errors in `[F1,...,Fd,a1,...,ad]` order. |"
    "covariance not positive definite": {"covariance": np.diag([0.25] * 3 + [-0.01] * 3)},
    # README.md, "## Scope and special cases":
    #   "- **Known direction:** if externally established, pass an oriented unit vector"
    "known_direction not a unit vector": {"known_direction": [2.0, 0.0, 0.0]},
}


@pytest.mark.parametrize("change", list(INVALID), ids=list(INVALID))
def test_invalid_inputs_raise_value_error(change):
    """The README example with one input made invalid (the invalid values are chosen here;
    what counts as invalid is the README's input table, quoted at each entry of INVALID)."""
    force, acceleration, covariance = _example()
    args = dict(force=force, acceleration=acceleration, covariance=covariance)
    args.update(INVALID[change])
    # README.md, "## How to use the Python function":
    #   "Invalid inputs raise `ValueError`; unresolved quadrature or numerical range"
    with pytest.raises(ValueError):
        estimate_mass(**args)


def test_unresolved_quadrature_raises_arithmetic_error():
    """The README example with rtol=1e-15 and max_refinements=1 (chosen here): the one
    allowed refinement cannot meet that rtol (its change is the example's 6.61e-5)."""
    # README.md, "## How to use the Python function":
    #   "Invalid inputs raise `ValueError`; unresolved quadrature or numerical range"
    #   "failures raise `ArithmeticError` rather than returning an unchecked number."
    with pytest.raises(ArithmeticError):
        estimate_mass(*_example(), rtol=1e-15, max_refinements=1)


def _direction_order_cases():
    force, acceleration, covariance = _example()
    return {
        # the README example: independent isotropic 3D errors, the analytic 3D path
        "analytic 3D": (force, acceleration, covariance, None),
        # section 4(a) of the 25 September record, row 1 (F=0.3, a=0.15, sigma_F=0.5, sigma_a=0.1)
        "1D": ([0.3], [0.15], np.diag([0.5**2, 0.1**2]), None),
        # the README example with a known direction chosen here
        "known direction": (force, acceleration, covariance, [1.0, 0.0, 0.0]),
    }


@pytest.mark.parametrize("case", ["analytic 3D", "1D", "known direction"])
def test_direction_order_is_unused(case):
    force, acceleration, covariance, known = _direction_order_cases()[case]
    a = estimate_mass(force, acceleration, covariance, known_direction=known, direction_order=24)
    b = estimate_mass(force, acceleration, covariance, known_direction=known, direction_order=100)
    # README.md, "## How to use the Python function":
    #   "and final resolution. The direction order is unused in the analytic 3D,"
    #   "known-direction and 1D direction calculations. Refinement is an observed"
    #   (every returned field except the echoed direction_order is identical; 100 is chosen here)
    assert dataclasses.replace(a, direction_order=0) == dataclasses.replace(b, direction_order=0)


@pytest.mark.parametrize("force,acceleration,covariance", [
    # section 4(a) of the 25 September record, row 1
    ([0.3], [0.15], [[0.5**2, 0.0], [0.0, 0.1**2]]),
    # the README example's x components
    ([10.0], [2.0], [[0.5**2, 0.0], [0.0, 0.1**2]]),
    # chosen here, with correlated channels
    ([-1.2], [0.7], [[0.25, 0.01], [0.01, 0.04]]),
])
def test_one_dimension_two_directions_with_equal_weight(force, acceleration, covariance):
    a = estimate_mass(force, acceleration, covariance)
    b = estimate_mass(-np.asarray(force), -np.asarray(acceleration), covariance)
    # README.md, "## The equation":
    #   "uniform probability measure. In one dimension the two directions, $+1$ and"
    #   "$-1$, each have weight $1/2$. Flat means flat in the two **physical positive"
    #   (consequence: reversing both measured vectors swaps the two directions' terms,
    #   which leaves the result unchanged exactly when their weights are equal)
    assert a == b


@pytest.mark.parametrize("d", [1, 2, 3])
@pytest.mark.parametrize("sf,sa", [(1.0, 1.0), (2.0, 0.5), (0.3, 3.0)])
@pytest.mark.parametrize("c", [0.5, 0.8, 0.95, 0.99])
def test_measured_zeros_half_cauchy_quantiles(d, sf, sa, c):
    """Half-Cauchy: P(X <= t) = (2/pi) arctan t, so the central interval of content c is
    [tan(pi (1-c)/4), tan(pi (1+c)/4)]. Sigmas and contents chosen here. The code's
    density in the ratio angle is exactly uniform here, so only rounding separates the
    two sides: asserted to 1e-12 relative."""
    s = sf / sa
    result = estimate_mass(np.zeros(d), np.zeros(d), np.diag([sf**2] * d + [sa**2] * d),
                           probability=c)
    # README.md, "## Scope and special cases":
    #   "  measured zero vectors and independent isotropic errors, $\widehat m="
    #   "  \sigma_F/\sigma_a$; $M/(\sigma_F/\sigma_a)$ is half-Cauchy. Its central"
    # README.md, "## Uncertainty: exactly what is returned":
    #   "I_c=[q_{(1-c)/2},q_{(1+c)/2}],\qquad"
    assert result.interval[0] / s == pytest.approx(math.tan(math.pi * (1 - c) / 4), rel=1e-12)
    # (the same two quotes)
    # README.md, "## Scope and special cases":
    #   "  \sigma_F/\sigma_a$; $M/(\sigma_F/\sigma_a)$ is half-Cauchy. Its central"
    # README.md, "## Uncertainty: exactly what is returned":
    #   "I_c=[q_{(1-c)/2},q_{(1+c)/2}],\qquad"
    assert result.interval[1] / s == pytest.approx(math.tan(math.pi * (1 + c) / 4), rel=1e-12)


@pytest.mark.parametrize("force_unit,acceleration_unit", [(1000.0, 1.0), (1.0, 100.0), (1000.0, 100.0)])
def test_unit_changes(force_unit, acceleration_unit):
    """The README example in other units (chosen here): force in mN (x1000) and/or
    acceleration in cm/s^2 (x100); each sigma scales with its channel."""
    base = estimate_mass(*_example())
    force, acceleration, _ = _example()
    covariance = np.diag([(0.5 * force_unit) ** 2] * 3 + [(0.1 * acceleration_unit) ** 2] * 3)
    other = estimate_mass(np.asarray(force) * force_unit, np.asarray(acceleration) * acceleration_unit,
                          covariance)
    k = force_unit / acceleration_unit
    # README.md, "## Types and inputs":
    #   "| $\widehat m$ / `result.mass` | A strictly positive scalar in force-unit / acceleration-unit; N divided by m/s² gives kg. |"
    # README.md, "## Version check and provenance":
    #   "magnitude/direction coordinates, rotations, unit changes, channel reciprocity,"
    assert other.mass == pytest.approx(base.mass * k, rel=1e-12)
    # (the same two quotes)
    # README.md, "## Types and inputs":
    #   "| $\widehat m$ / `result.mass` | A strictly positive scalar in force-unit / acceleration-unit; N divided by m/s² gives kg. |"
    # README.md, "## Version check and provenance":
    #   "magnitude/direction coordinates, rotations, unit changes, channel reciprocity,"
    assert other.interval == pytest.approx((base.interval[0] * k, base.interval[1] * k), rel=1e-12)
    # README.md, "## Uncertainty: exactly what is returned":
    #   "The log spread is dimensionless and independent of the chosen reference mass"
    assert other.log_standard_deviation == pytest.approx(base.log_standard_deviation, abs=1e-12)


def test_rotation():
    """The README example, both vectors turned by one orthogonal matrix (chosen here:
    QR of default_rng(0).normal(size=(3, 3))). The covariance is isotropic in each channel,
    so it is unchanged by the turn."""
    base = estimate_mass(*_example())
    force, acceleration, covariance = _example()
    q, _ = np.linalg.qr(np.random.default_rng(0).normal(size=(3, 3)))
    turned = estimate_mass(q @ np.asarray(force), q @ np.asarray(acceleration), covariance)
    # README.md, "## The equation":
    #   "Here $S^{d-1}$ is the unit sphere in $d$ dimensions and $d\Omega$ is its"
    #   "uniform probability measure. In one dimension the two directions, $+1$ and"
    # README.md, "## Version check and provenance":
    #   "magnitude/direction coordinates, rotations, unit changes, channel reciprocity,"
    #   (consequence: a uniform direction and an isotropic covariance are unchanged by the turn)
    assert turned.mass == pytest.approx(base.mass, rel=1e-12)
    # (the same quotes and consequence as just above)
    # README.md, "## The equation":
    #   "Here $S^{d-1}$ is the unit sphere in $d$ dimensions and $d\Omega$ is its"
    # README.md, "## Version check and provenance":
    #   "magnitude/direction coordinates, rotations, unit changes, channel reciprocity,"
    assert turned.interval == pytest.approx(base.interval, rel=1e-12)
    # (the same quotes and consequence as just above)
    # README.md, "## The equation":
    #   "Here $S^{d-1}$ is the unit sphere in $d$ dimensions and $d\Omega$ is its"
    # README.md, "## Version check and provenance":
    #   "magnitude/direction coordinates, rotations, unit changes, channel reciprocity,"
    assert turned.log_standard_deviation == pytest.approx(base.log_standard_deviation, abs=1e-12)


def test_channel_reciprocity():
    """The README example with the force and acceleration channels exchanged (vectors and
    covariance blocks). In the boxed equation this exchanges f and alpha, and the flat
    measure df d(alpha) is unchanged by that, so the point becomes E_P[alpha]/E_P[f]."""
    base = estimate_mass(*_example())
    force, acceleration, _ = _example()
    swapped = estimate_mass(acceleration, force, np.diag([0.1**2] * 3 + [0.5**2] * 3))
    # README.md, "## The equation":
    #   "=\frac{E_P[f]}{E_P[\alpha]}"
    # README.md, "## Version check and provenance":
    #   "magnitude/direction coordinates, rotations, unit changes, channel reciprocity,"
    #   (consequence: the exchanged point is the reciprocal of the original)
    assert swapped.mass * base.mass == pytest.approx(1.0, rel=1e-12)
    # README.md, "## Uncertainty: exactly what is returned":
    #   "Use the **same** law $P$ to define the random scalar $M=f/\alpha$. Let"
    # README.md, "## Version check and provenance":
    #   "magnitude/direction coordinates, rotations, unit changes, channel reciprocity,"
    #   (consequence: M becomes 1/M, so the interval becomes the reversed reciprocal)
    assert swapped.interval == pytest.approx((1.0 / base.interval[1], 1.0 / base.interval[0]), rel=1e-12)
    # README.md, "## Uncertainty: exactly what is returned":
    #   "Use the **same** law $P$ to define the random scalar $M=f/\alpha$. Let"
    # README.md, "## Version check and provenance":
    #   "magnitude/direction coordinates, rotations, unit changes, channel reciprocity,"
    #   (consequence: log M changes sign, so its SD is unchanged)
    assert swapped.log_standard_deviation == pytest.approx(base.log_standard_deviation, abs=1e-12)


def _truncated_normal_mean(mu, sigma):
    """Mean of N(mu, sigma^2) restricted to (0, inf)."""
    return mu + sigma * norm.pdf(mu / sigma) / norm.cdf(mu / sigma)


@pytest.mark.parametrize("u", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.6, 0.8, 0.0], [-1.0, 0.0, 0.0]])
def test_known_direction_truncated_normal_means(u):
    """The README example with a known direction u (chosen here). With u fixed and a
    diagonal covariance, Q splits into (f - u.F)^2/sigma_F^2 + (alpha - u.a)^2/sigma_a^2 plus
    terms free of f and alpha, so f and alpha are independent normals truncated to (0, inf).
    The code integrates numerically; asserted to the README's rtol 2e-4."""
    force, acceleration, covariance = _example()
    result = estimate_mass(force, acceleration, covariance, known_direction=u)
    mean_f = _truncated_normal_mean(float(np.dot(u, force)), 0.5)
    mean_a = _truncated_normal_mean(float(np.dot(u, acceleration)), 0.1)
    # README.md, "## Scope and special cases":
    #   "- **Known direction:** if externally established, pass an oriented unit vector"
    #   "  as `known_direction`. This replaces $d\Omega$ by a point mass at that"
    #   "  direction. The function integrates only the two magnitudes. Do not infer"
    # README.md, "## Version check and provenance":
    #   "correlated-null quantiles and log spread, known-direction truncated-normal"
    #   "means, agreement with the latest isotropic backend and the original full"
    # README.md, "## How to use the Python function":
    #   "`max_refinements=3`; `rtol=2e-4` applies to all returned numerical summaries."
    assert result.mean_force_magnitude == pytest.approx(mean_f, rel=2e-4)
    # (the same quotes as just above)
    # README.md, "## Version check and provenance":
    #   "correlated-null quantiles and log spread, known-direction truncated-normal"
    #   "means, agreement with the latest isotropic backend and the original full"
    assert result.mean_acceleration_magnitude == pytest.approx(mean_a, rel=2e-4)
    # README.md, "## The equation":
    #   "=\frac{E_P[f]}{E_P[\alpha]}"
    # README.md, "## How to use the Python function":
    #   "`max_refinements=3`; `rtol=2e-4` applies to all returned numerical summaries."
    assert result.mass == pytest.approx(mean_f / mean_a, rel=2e-4)


# ---------------------------------------------------------------- (a'') section 4(a) of the 25 September record

AXIS_TABLE = [
    # axis_splitting_and_the_24_september_change_2026-09-25.md, "## 4. What axis splitting forces":
    #   "median estimate ($\sigma_F=0.5$, $\sigma_a=0.1$, so $s=5$):"
    #   "| reading | 21 Sept | 24 Sept |"
    #   "| $F=0.3$, $a=0.15$ as 1D | 3.12 | 3.46 |"
    pytest.param([0.3], [0.15], 3.12, id="1D"),
    #   "| the same, written as a 3D vector $(x,0,0)$ | 3.73 | 3.46 |"
    pytest.param([0.3, 0.0, 0.0], [0.15, 0.0, 0.0], 3.73, id="3D (x,0,0)"),
    #   "| $\mathbf F=(0.3,0.2,-0.4)$, $\mathbf a=(0.15,-0.05,0.02)$ as one 3D vector | 3.97 | 3.72 |"
    pytest.param([0.3, 0.2, -0.4], [0.15, -0.05, 0.02], 3.97, id="3D vector",
                 marks=pytest.mark.xfail(strict=True, reason=(
                     "record states 3.97 (21 Sept column, 2 decimals); code gives 3.9646, which "
                     "rounds to 3.96: tester flat single median 3.9646497011704036 (3.964649757675297 "
                     "at grid_points 8001); mass_estimator.py central interval at content 1e-6 "
                     "[3.964644834858348, 3.9646563101202186] at default settings"))),
]


@pytest.mark.parametrize("force,acceleration,stated", AXIS_TABLE)
def test_axis_splitting_table_21_sept_median_tester(force, acceleration, stated):
    median = float(estimate(force, acceleration, 0.5, 0.1, ["median"], law=FLAT, combine=SINGLE)["median"])
    # axis_splitting_and_the_24_september_change_2026-09-25.md, "## 4. What axis splitting forces":
    #   "**(a) The 21 September equation fails it.** One reading, same numbers,"
    #   "median estimate ($\sigma_F=0.5$, $\sigma_a=0.1$, so $s=5$):"
    #   "| reading | 21 Sept | 24 Sept |"
    #   "| $F=0.3$, $a=0.15$ as 1D | 3.12 | 3.46 |"
    #   "| the same, written as a 3D vector $(x,0,0)$ | 3.73 | 3.46 |"
    #   "| $\mathbf F=(0.3,0.2,-0.4)$, $\mathbf a=(0.15,-0.05,0.02)$ as one 3D vector | 3.97 | 3.72 |"
    #   (the 21 Sept column of this case's row, printed to 2 decimals)
    assert median == pytest.approx(stated, abs=printed(2))


@pytest.mark.parametrize("force,acceleration,stated", AXIS_TABLE)
def test_axis_splitting_table_21_sept_median_standalone(force, acceleration, stated):
    """The standalone script has no median output. Its central interval of content c,
    [q_(1-c)/2, q_(1+c)/2], contains the median q_1/2 for every c; at c = 1e-6 (chosen
    here) both ends must then lie within the printed precision of the stated median."""
    d = len(force)
    result = estimate_mass(force, acceleration, np.diag([0.5**2] * d + [0.1**2] * d), probability=1e-6)
    # axis_splitting_and_the_24_september_change_2026-09-25.md, "## 4. What axis splitting forces":
    #   "**(a) The 21 September equation fails it.** One reading, same numbers,"
    #   "median estimate ($\sigma_F=0.5$, $\sigma_a=0.1$, so $s=5$):"
    #   "| reading | 21 Sept | 24 Sept |"
    #   "| $F=0.3$, $a=0.15$ as 1D | 3.12 | 3.46 |"
    #   "| the same, written as a 3D vector $(x,0,0)$ | 3.73 | 3.46 |"
    #   "| $\mathbf F=(0.3,0.2,-0.4)$, $\mathbf a=(0.15,-0.05,0.02)$ as one 3D vector | 3.97 | 3.72 |"
    #   (the 21 Sept column of this case's row, printed to 2 decimals)
    assert result.interval[0] == pytest.approx(stated, abs=printed(2))
    # axis_splitting_and_the_24_september_change_2026-09-25.md, "## 4. What axis splitting forces":
    #   "**(a) The 21 September equation fails it.** One reading, same numbers,"
    #   "median estimate ($\sigma_F=0.5$, $\sigma_a=0.1$, so $s=5$):"
    #   "| reading | 21 Sept | 24 Sept |"
    #   "| $F=0.3$, $a=0.15$ as 1D | 3.12 | 3.46 |"
    #   "| the same, written as a 3D vector $(x,0,0)$ | 3.73 | 3.46 |"
    #   "| $\mathbf F=(0.3,0.2,-0.4)$, $\mathbf a=(0.15,-0.05,0.02)$ as one 3D vector | 3.97 | 3.72 |"
    #   (the 21 Sept column of this case's row, printed to 2 decimals)
    assert result.interval[1] == pytest.approx(stated, abs=printed(2))


# ---------------------------------------------------------------- (b) tester against the script

STATED = 2e-10


def _inputs(d, k):
    """Case k of dimension d: k = 0..4 seeded random, k = "readme" the README example (3D)."""
    if k == "readme":
        force, acceleration, _ = _example()
        return np.array(force), np.array(acceleration), 0.5, 0.1
    rng = np.random.default_rng(100 + d)
    for _ in range(k + 1):
        f = rng.normal(size=d) * 3
        a = rng.normal(size=d)
        sf, sa = rng.uniform(0.3, 2, size=2)
    return f, a, sf, sa


@lru_cache(maxsize=None)
def _both(d, k):
    f, a, sf, sa = _inputs(d, k)
    script = estimate_mass(f, a, np.diag([sf**2] * d + [sa**2] * d)).mass
    tester = float(estimate(f, a, sf, sa, ["ratio_of_means"], law=FLAT, combine=SINGLE)["ratio_of_means"])
    return script, tester


CASES = [(d, k) for d in (1, 2, 3) for k in range(5)] + [(3, "readme")]

# Measured on this checkout (both implementations at their default settings).
CONTEXT = ("; inputs chosen here, not by the record; both codes at default settings. Measured "
           "here: with mass_estimator.py started at ratio_order=8192 it agrees with the tester "
           "to 1e-14 relative on every case")
RELATIVE_FAILS = {
    (3, 0): "record states agreement to 2e-10 (read as relative); code gives relative difference 7.4e-10 "
            "(mass_estimator.py 8.193959350182142, tester 8.193959344154637)" + CONTEXT,
    (3, 3): "record states agreement to 2e-10 (read as relative); code gives relative difference 5.9e-10 "
            "(mass_estimator.py 10.941674097778959, tester 10.941674091277036)" + CONTEXT,
}
ABSOLUTE_FAILS = {
    (2, 0): "record states agreement to 2e-10 (read as absolute); code gives absolute difference 7.2e-10 "
            "(mass_estimator.py 9.927649806640245, tester 9.927649807358035)" + CONTEXT,
    (3, 0): "record states agreement to 2e-10 (read as absolute); code gives absolute difference 6.0e-9 "
            "(mass_estimator.py 8.193959350182142, tester 8.193959344154637)" + CONTEXT,
    (3, 1): "record states agreement to 2e-10 (read as absolute); code gives absolute difference 1.6e-9 "
            "(mass_estimator.py 13.117150805044396, tester 13.117150803463602)" + CONTEXT,
    (3, 3): "record states agreement to 2e-10 (read as absolute); code gives absolute difference 6.5e-9 "
            "(mass_estimator.py 10.941674097778959, tester 10.941674091277036)" + CONTEXT,
}


def _params(fails):
    return [pytest.param(d, k, id=f"{d}d-{k}",
                         marks=[pytest.mark.xfail(strict=True, reason=fails[(d, k)])] if (d, k) in fails else [])
            for d, k in CASES]


@pytest.mark.parametrize("d,k", _params(RELATIVE_FAILS))
def test_flat_reference_agrees_with_the_script_relative(d, k):
    script, tester = _both(d, k)
    # axis_splitting_and_the_24_september_change_2026-09-25.md, "## 10. Reproduction":
    #   "The tester's `flat` reference agrees with the standalone"
    #   "`mass_estimator_equation/mass_estimator.py` to $2\times10^{-10}$ in 1D, 2D"
    #   "and 3D."
    #   (read as a relative difference)
    assert abs(tester - script) <= STATED * abs(script)


@pytest.mark.parametrize("d,k", _params(ABSOLUTE_FAILS))
def test_flat_reference_agrees_with_the_script_absolute(d, k):
    script, tester = _both(d, k)
    # axis_splitting_and_the_24_september_change_2026-09-25.md, "## 10. Reproduction":
    #   "The tester's `flat` reference agrees with the standalone"
    #   "`mass_estimator_equation/mass_estimator.py` to $2\times10^{-10}$ in 1D, 2D"
    #   "and 3D."
    #   (read as an absolute difference in the mass, in the inputs' force / acceleration units)
    assert abs(tester - script) <= STATED
