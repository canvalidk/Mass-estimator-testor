"""Checks with explicit authority, not scores silently promoted into requirements.

A Check is a small callable plus its reason. Add or replace entries by passing a
different check list to inspect_case; no registry or subclass is needed. A check
can state unresolved or inapplicable instead of inventing an expected answer.
"""
from dataclasses import dataclass, field
from fractions import Fraction
import hashlib
import json
import math

import numpy as np

from .transforms import append_null, channel_units, orthogonal_frames, reversed_readings


@dataclass
class Context:
    inputs: object
    masses: np.ndarray
    evaluation: dict
    diagnostics: dict
    adapter: object
    tolerance: float = 1e-7
    _cache: dict = field(default_factory=dict)

    def evaluate(self, inputs, masses=None):
        masses = self.masses if masses is None else np.asarray(masses)
        key = hashlib.sha256((json.dumps(inputs.to_dict(), sort_keys=True, default=str)
                              + repr(masses.tolist())).encode()).hexdigest()
        if key not in self._cache:
            self._cache[key] = self.adapter(inputs, masses, variants=("quoted",), include_size=False)
        return self._cache[key]


@dataclass(frozen=True)
class Check:
    id: str
    title: str
    kind: str
    source: str
    expected: str
    run: object
    assumptions: tuple = ()

    def __call__(self, context):
        result = dict(id=self.id, title=self.title, kind=self.kind, source=self.source,
                      expected=self.expected, assumptions=list(self.assumptions))
        result.update(self.run(context))
        return result


def _finding(context, error, observed, detail=""):
    return dict(status="pass" if math.isfinite(error) and error <= context.tolerance else "fail",
                observed=observed, tolerance={"dimensionless_error": context.tolerance}, detail=detail)


def _unavailable(detail):
    return dict(status="inapplicable", detail=detail, observed={}, tolerance={})


def _law_errors(base, other, mass_scale=1.0):
    a, b = np.asarray(base["log_kernel"]), np.asarray(other["log_kernel"])
    difference = a - b
    # A constant factor does not change a normalized probability law.
    kernel_error = float(np.max(np.abs(difference - difference[len(difference) // 2])))
    cdfs = []
    for law, shift in ((base, 0), (other, math.log(mass_scale))):
        x, density = np.asarray(law["log_mass"]) - shift, np.asarray(law["density"])
        if (not np.all(np.isfinite(density) & (density >= 0)) or not np.all(np.diff(x) > 0)):
            raise ValueError("native law grids must increase and carry finite nonnegative density")
        area = (density[:-1] + density[1:]) * np.diff(x) / 2
        cdf = np.r_[0, np.cumsum(area)]
        if not cdf[-1] > 0:
            raise ValueError("native law has no positive area")
        cdfs.append((x, cdf / cdf[-1]))
    grid = np.union1d(cdfs[0][0], cdfs[1][0])
    native_error = float(np.max(np.abs(np.interp(grid, *cdfs[0]) - np.interp(grid, *cdfs[1]))))
    return dict(sampled_log_kernel_error_modulo_constant=kernel_error, native_cdf_maximum_error=native_error)


def _point_error(base, other, scale=1):
    return abs(math.log(other["point"] / (base["point"] * scale)))


def projection(context):
    errors = []
    for reading, p in zip(context.inputs.readings, context.diagnostics["readings"]):
        norm = 1 + np.max(np.abs(p["candidate_force"]))
        errors.append(np.max(np.abs(p["candidate_force"] - context.masses[:, None] * p["candidate_acceleration"])) / norm)
        errors.append(np.max(np.abs(p["declared_cost"] - p["outside_subspace_cost"] - p["within_subspace_cost"]))
                      / (1 + p["data_energy"]))
        errors.append(max(0, float(np.max(p["residual_cost"] - p["declared_cost"]))) / (1 + p["data_energy"]))
        errors.append(np.max(np.abs(p["candidate_acceleration"] - p["conditional_mean"] @ reading.subspace.T))
                      / (1 + np.max(np.abs(p["candidate_acceleration"]))))
    error = float(max(errors, default=0))
    return _finding(context, error, {"maximum_scaled_constraint_or_energy_error": error})


def explained_energy(context):
    internal = context.evaluation.get("per_reading", [])
    if len(internal) != len(context.inputs.readings):
        return _unavailable("Adapter does not expose per-reading explained-energy curves.")
    errors = [float(np.max(np.abs(a["explained_energy"] - b["explained_energy"])) / (1 + b["data_energy"]))
              for a, b in zip(internal, context.diagnostics["readings"])]
    return _finding(context, max(errors), {"maximum_scaled_energy_error": max(errors)},
                    "Independent QR projection checks the computation of H, not the physical adequacy of the law.")


def frame_invariance(context):
    base = context.evaluation["laws"]["quoted"]
    changed = context.evaluate(orthogonal_frames(context.inputs))["laws"]["quoted"]
    law, point = _law_errors(base, changed), _point_error(base, changed)
    return _finding(context, max(*law.values(), point), dict(law, absolute_log_point_ratio=point))


def unit_covariance(context):
    if (context.inputs.prior.get("family") != "half_cauchy"
            or set(context.inputs.prior) != {"family", "center"}):
        return _unavailable("This unit transform is currently justified for the half-Cauchy center; other prior declarations need their own transform.")
    converted, scale = channel_units(context.inputs)
    changed = context.evaluate(converted, context.masses * scale)["laws"]["quoted"]
    base = context.evaluation["laws"]["quoted"]
    law, point = _law_errors(base, changed, scale), _point_error(base, changed, scale)
    return _finding(context, max(*law.values(), point), dict(law, mass_unit_scale=scale,
        absolute_log_scaled_point_ratio=point))


def order_invariance(context):
    if len(context.inputs.readings) < 2:
        return _unavailable("A single reading has no nontrivial ordering.")
    base = context.evaluation["laws"]["quoted"]
    changed = context.evaluate(reversed_readings(context.inputs))["laws"]["quoted"]
    law, point = _law_errors(base, changed), _point_error(base, changed)
    return _finding(context, max(*law.values(), point), dict(law, absolute_log_point_ratio=point))


def null_law(context):
    base = context.evaluation["laws"]["quoted"]
    changed = context.evaluate(append_null(context.inputs))["laws"]["quoted"]
    errors = _law_errors(base, changed)
    return _finding(context, max(errors.values()), dict(errors, null_readings_added=5))


def null_point(context):
    base = context.evaluation["laws"]["quoted"]
    changed = context.evaluate(append_null(context.inputs))["laws"]["quoted"]
    error = _point_error(base, changed)
    return _finding(context, error, {"null_readings_added": 5, "before": base["point"], "after": changed["point"],
                    "absolute_log_point_ratio": error}, "A failure is the existing O2/T29 issue; it is not suppressed by the inspector.")


def reciprocal_median(context):
    if (context.inputs.prior.get("family") != "half_cauchy"
            or set(context.inputs.prior) != {"family", "center"}):
        return _unavailable("This symmetry check currently establishes the prior's reciprocal symmetry only for the declared half-Cauchy family.")
    if len(context.inputs.readings) != 1:
        return _unavailable("This exact reciprocal-symmetry check is currently for a single full-subspace reading.")
    r = context.inputs.readings[0]
    d = r.force.size
    sf, sa = r.covariance[0, 0], r.covariance[d, d]
    expected_covariance = np.diag([sf] * d + [sa] * d)
    c = context.inputs.prior.get("center")
    if r.subspace.shape[1] != d or not np.array_equal(r.covariance, expected_covariance) or c is None:
        return _unavailable("Requires exact equal whitened energies, isotropic independent errors and c=s.")
    # Near-symmetry is not enough in a law with a nearly flat central CDF: a tiny
    # real asymmetry can move the median greatly. Establish this exact identity
    # on the declared floating-point numbers as rationals, rather than assuming
    # its premise from a loose relative-energy tolerance.
    force_energy = sum(Fraction(float(v)) ** 2 for v in r.force) / Fraction(float(sf))
    acceleration_energy = sum(Fraction(float(v)) ** 2 for v in r.acceleration) / Fraction(float(sa))
    if (force_energy != acceleration_energy
            or Fraction(float(c)) ** 2 * Fraction(float(sa)) != Fraction(float(sf))):
        return _unavailable("Exact reciprocal symmetry is not established by these declarations; near-equal energies are insufficient.")
    median = context.evaluation["laws"]["quoted"]["median"]
    error = abs(math.log(median / c))
    return _finding(context, error, {"mathematical_median": c, "implemented_median": median,
        "absolute_log_median_ratio": error}, "A near-flat CDF can satisfy a probability tolerance while missing its unique median (O13).")


def compatibility_role(context):
    return dict(status="unresolved", observed={"sampled_minimum_total_declared_cost":
                float(np.min(context.diagnostics["total_declared_cost"]))}, tolerance={},
                detail="This is a sampled minimum, not a certified global infimum. No accepted requirement yet selects how absolute incompatibility must change the complete result; no rejection criterion is imposed.")


DEFAULT_CHECKS = (
    Check("projection", "Compatible-fit identities", "algebra", "Covariance propagation and orthogonal projection",
          "Candidates satisfy Newton II and the declared subspace; nested projection costs agree.", projection),
    Check("explained_energy", "Independent explained-energy calculation", "implementation", "T8; independent whitened QR projection",
          "Estimator H agrees with an independent least-squares projection.", explained_energy),
    Check("frames", "Independent orthogonal frames", "physical", "L16; definition/premises.md",
          "Transforming each reading, covariance and subspace preserves sampled log kernels, native-grid CDFs and the point.", frame_invariance),
    Check("units", "Units covariance", "physical", "L3; definition/premises.md",
          "Sampled log kernels and native-grid CDFs transform with mass units; the point scales by force_scale/acceleration_scale.", unit_covariance),
    Check("order", "Reading order", "algebra", "T12; product combination",
          "Reordering declared independent readings preserves sampled log kernels, native-grid CDFs and point.", order_invariance,
          ("No cross-reading errors are represented by this input type.",)),
    Check("null_law", "Null padding: law", "physical", "L20 and T24; definition/premises.md",
          "Appending exact measured-zero readings preserves sampled log kernels and native-grid CDFs.", null_law),
    Check("null_point", "Null padding: point", "physical", "L20, T29 and problems/O2.md",
          "Appending exact measured-zero readings leaves the point unchanged.", null_point),
    Check("reciprocal_median", "Reciprocal mathematical median", "numerical", "Reciprocal symmetry; problems/O13.md",
          "With equal whitened energies and c=s, the mathematical median is c.", reciprocal_median),
    Check("compatibility_role", "Required response to absolute incompatibility", "research", "problems/O13.md; TASK.md",
          "The required relation between absolute correction cost and the complete output remains to be derived.", compatibility_role),
)


def run_checks(context, checks=DEFAULT_CHECKS):
    results = []
    for check in checks:
        try:
            results.append(check(context))
        except Exception as error:
            # An incomplete adapter/check is an operational error, never a pass.
            results.append(dict(id=getattr(check, "id", getattr(check, "__name__", "custom")),
                title=getattr(check, "title", "Custom check"), status="error", kind="execution",
                source=getattr(check, "source", "custom"), expected=getattr(check, "expected", ""),
                observed={}, tolerance={}, assumptions=[], detail=f"{type(error).__name__}: {error}"))
    return results
