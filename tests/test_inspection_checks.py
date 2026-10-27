"""The inspector reports justified failures and permits local experimental checks.

These are framework tests, not new pass criteria for the mass equation. In
particular, reporting the already known L20 point failure is a successful check
of the inspector, rather than an expected failure of this test module.
"""
from copy import deepcopy

import numpy as np
import pytest

from met.inspection.adapters import ours29_adapter
from met.inspection.case import Case, Inputs
from met.inspection.case import load_cases
from met.inspection.checks import Check, DEFAULT_CHECKS
from met.inspection.run import inspect_case
from met.inspection.transforms import orthogonal_frames


GRID = np.exp(np.linspace(-4.0, 4.0, 41))


def selected_checks(*identifiers):
    by_id = {check.id: check for check in DEFAULT_CHECKS}
    return tuple(by_id[identifier] for identifier in identifiers)


def quoted_only_adapter(inputs, masses, *, variants=("quoted",), include_size=False):
    """The diagnostics under test need one law, so omit unused variants and A."""
    return ours29_adapter(inputs, masses, variants=("quoted",), include_size=False)


def findings(report):
    return {finding["id"]: finding for finding in report["checks"]}


@pytest.fixture(scope="module")
def null_case():
    # The same seeded full-covariance world as the established T24/T29 check in
    # tests/test_full.py; the synthetic truth stays outside the case inputs.
    rng = np.random.default_rng(43)
    readings = []
    for index in range(4):
        mixing = rng.normal(size=(6, 6))
        covariance = mixing @ mixing.T / 6 + 0.3 * np.eye(6)
        push = rng.normal(size=3) * 0.7
        y = np.r_[2 * push, push] + np.linalg.cholesky(covariance) @ rng.normal(size=6)
        readings.append(dict(id=f"r{index}", force=y[:3].tolist(), acceleration=y[3:].tolist(),
                             covariance=covariance.tolist()))
    return Case.from_dict(dict(name="known_null_padding", question="Does L20 hold at both outputs?",
                               readings=readings, prior={"family": "half_cauchy", "center": 2.0}))


@pytest.fixture(scope="module")
def null_report(null_case):
    return inspect_case(null_case, masses=GRID, adapter=quoted_only_adapter,
                        checks=selected_checks("null_law", "null_point"))


def test_independent_null_checks_distinguish_law_from_point(null_report):
    checks = findings(null_report)
    assert set(checks) == {"null_law", "null_point"}
    assert checks["null_law"]["status"] == "pass"
    assert checks["null_point"]["status"] == "fail"
    assert "L20" in checks["null_point"]["source"]
    assert "O2" in checks["null_point"]["detail"]
    observed = checks["null_point"]["observed"]
    assert observed["before"] != observed["after"]
    assert observed["absolute_log_point_ratio"] > checks["null_point"]["tolerance"]["dimensionless_error"]


def test_known_physical_failure_is_returned_as_data_without_aborting(null_report):
    """The open O2 requirement must remain visible in an otherwise usable result."""
    assert findings(null_report)["null_point"]["kind"] == "physical"
    assert np.isfinite(null_report["evaluation"]["laws"]["quoted"]["point"])
    assert len(null_report["diagnostics"]["readings"]) == 4
    assert all(check["status"] != "error" for check in null_report["checks"])


@pytest.fixture(scope="module")
def reciprocal_report():
    case = Case.from_dict(dict(name="opposite_equal_energy", question="Numerical median versus reciprocal symmetry",
        prior={"family": "half_cauchy", "center": 1.0}, readings=[dict(id="opposite", force=[1.0, 0.0],
        acceleration=[-1.0, 0.0], covariance=(0.1 ** 2 * np.eye(4)).tolist())]))
    return inspect_case(case, masses=GRID, adapter=quoted_only_adapter,
                        checks=selected_checks("reciprocal_median"))


def test_reciprocal_symmetry_check_passes_since_the_plateau_fix(reciprocal_report):
    """O13 proves the unique mathematical median is one. Until 27 Oct the grid gave
    2.3416 here (the flat-CDF defect); met.full._plateau_centre fixed it."""
    check = findings(reciprocal_report)["reciprocal_median"]
    assert check["status"] == "pass"
    assert check["kind"] == "numerical"
    assert "O13" in check["source"]
    assert check["observed"]["mathematical_median"] == 1.0
    assert check["observed"]["absolute_log_median_ratio"] < 1e-6
    assert check["observed"]["implemented_median"] == reciprocal_report["evaluation"]["laws"]["quoted"]["median"]


def test_reciprocal_median_check_still_fails_on_a_wrong_median():
    """The detector must still fire: an adapter that reports the old defective median."""
    def defective(inputs, masses, **kw):
        out = quoted_only_adapter(inputs, masses, **kw)
        out["laws"]["quoted"]["median"] = 2.3416363971870715
        return out
    case = Case.from_dict(dict(name="opposite_equal_energy", prior={"family": "half_cauchy", "center": 1.0},
        readings=[dict(id="opposite", force=[1.0, 0.0], acceleration=[-1.0, 0.0],
                       covariance=(0.1 ** 2 * np.eye(4)).tolist())]))
    report = inspect_case(case, masses=GRID, adapter=defective, checks=selected_checks("reciprocal_median"))
    check = findings(report)["reciprocal_median"]
    assert check["status"] == "fail"
    assert check["observed"]["absolute_log_median_ratio"] > 1e-2


@pytest.fixture(scope="module")
def correlated_report():
    covariance = np.array([[0.25, 0.03, 0.04, -0.01], [0.03, 0.36, 0.02, 0.05],
                           [0.04, 0.02, 0.16, 0.01], [-0.01, 0.05, 0.01, 0.25]])
    case = Case.from_dict(dict(name="correlated_restricted_series", prior={"family": "half_cauchy", "center": 1.7},
        readings=[dict(id="line", force=[3.1, 0.4], acceleration=[1.3, 0.15], covariance=covariance.tolist(),
                       subspace=[[1.0], [0.0]]),
                  dict(id="plane", force=[1.2, 2.6], acceleration=[0.6, 1.0], covariance=(1.3 * covariance).tolist()),
                  dict(id="tilted_line", force=[0.8, -0.6], acceleration=[0.25, -0.2],
                       covariance=(0.8 * covariance).tolist(), subspace=[[0.6], [0.8]])]))
    return inspect_case(case, masses=GRID, adapter=quoted_only_adapter,
                        checks=selected_checks("frames", "units", "order"))


def test_frames_units_and_order_preserve_correlated_restricted_series(correlated_report):
    checks = findings(correlated_report)
    assert set(checks) == {"frames", "units", "order"}
    for identifier, check in checks.items():
        assert check["status"] == "pass", (identifier, check)
        assert check["source"] and check["expected"]
        assert check["observed"] and check["tolerance"]


def test_investigator_truth_cannot_reach_or_change_estimator(null_case, null_report):
    first_truth = {"mass": 2.0, "latent_acceleration": [[1.0, 0.0, 0.0]], "injected_errors": [0.1]}
    second_truth = {"mass": 900.0, "latent_acceleration": [[-5.0, 4.0, 3.0]], "injected_errors": [-200.0]}
    received = []

    def spy(inputs, masses, *, variants=("quoted",), include_size=False):
        assert isinstance(inputs, Inputs)
        assert not hasattr(inputs, "investigator_truth")
        assert not hasattr(inputs, "truth")
        assert "truth" not in inputs.to_dict()
        received.append(inputs)
        return quoted_only_adapter(inputs, masses, variants=variants, include_size=include_size)

    # One real evaluation also proves truth cannot change the native estimator's
    # result; the second truth is exercised through the same boundary cheaply.
    marked_case = deepcopy(null_case)
    marked_case.investigator_truth = first_truth
    marked = inspect_case(marked_case, masses=GRID, adapter=spy, checks=())
    base, with_truth = null_report["evaluation"], marked["evaluation"]
    assert with_truth["laws"]["quoted"]["point"] == base["laws"]["quoted"]["point"]
    assert with_truth["laws"]["quoted"]["median"] == base["laws"]["quoted"]["median"]
    np.testing.assert_array_equal(with_truth["laws"]["quoted"]["log_kernel"], base["laws"]["quoted"]["log_kernel"])
    assert marked["investigator_truth"] == first_truth

    def boundary_spy(inputs, masses, *, variants=("quoted",), include_size=False):
        assert isinstance(inputs, Inputs)
        assert not hasattr(inputs, "investigator_truth")
        assert not hasattr(inputs, "truth")
        received.append(inputs)
        return {"adapter": "boundary_spy", "laws": {"quoted": {"point": inputs.prior["center"]}}}

    marked_case.investigator_truth = second_truth
    second = inspect_case(marked_case, masses=GRID, adapter=boundary_spy, checks=())
    assert second["investigator_truth"] == second_truth
    assert len(received) == 2
    assert received[0].to_dict() == received[1].to_dict() == null_case.inputs.to_dict()


def test_checks_can_be_selected_replaced_and_supplied_as_plain_callables(null_case):
    seen = []

    def local_projection(context):
        assert isinstance(context.inputs, Inputs)
        assert not hasattr(context, "investigator_truth")
        seen.append("replacement")
        return dict(status="pass", observed={"readings": len(context.inputs.readings)}, tolerance={})

    replacement = Check("projection", "A locally replaced diagnostic", "algebra", "Local test premise",
                        "Report the number of declared readings.", local_projection)

    def tentative_check(context):
        seen.append("plain_callable")
        return dict(id="local_open_question", status="unresolved", detail="A local research question can stay open.")

    def adapter(inputs, masses, *, variants=("quoted",), include_size=False):
        return {"adapter": "minimal", "laws": {"quoted": {"point": 2.0}}}

    result = inspect_case(null_case, masses=GRID, adapter=adapter, checks=(replacement, tentative_check))
    checks = findings(result)
    assert set(checks) == {"projection", "local_open_question"}
    assert checks["projection"]["title"] == "A locally replaced diagnostic"
    assert checks["projection"]["observed"] == {"readings": 4}
    assert checks["projection"]["status"] == "pass"
    assert checks["local_open_question"]["status"] == "unresolved"
    assert seen == ["replacement", "plain_callable"]


@pytest.mark.parametrize("case", load_cases("presets/inspection_cases.yaml")[1], ids=lambda case: case.name)
def test_saved_example_frames_preserve_symmetric_covariance_and_outputs(case):
    changed = orthogonal_frames(case.inputs)
    for reading in changed.readings:
        np.testing.assert_array_equal(reading.covariance, reading.covariance.T)
    report = inspect_case(case, masses=GRID, adapter=quoted_only_adapter, checks=selected_checks("frames"))
    check = report["checks"][0]
    assert check["status"] == "pass", check
