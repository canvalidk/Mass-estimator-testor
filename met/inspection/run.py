"""Thin orchestration and portable output for separately usable inspection functions."""
from collections import Counter
from collections.abc import Mapping
from datetime import datetime, timezone
import hashlib
from html import escape
from importlib import metadata
import json
import math
from pathlib import Path
import re

import numpy as np
import yaml

from .adapters import ours29_adapter
from .case import load_cases, parse_cases
from .checks import Context, DEFAULT_CHECKS, run_checks
from .geometry import profile


def validate(document):
    cases = parse_cases(document)
    return dict(kind="inspection", raw=document, cases=cases,
                cells=[dict(id=f"case_{index:03d}", name=case.name) for index, case in enumerate(cases)])


def inspect_case(case, *, masses=None, adapter=ours29_adapter, diagnostic=profile,
                 checks=DEFAULT_CHECKS, numerics=None):
    """Compose ordinary functions; call any stage separately or replace it here.

    Only case.inputs is passed to the adapter and checks. Truth stays in a separate
    output section for the investigator. The diagnostic window samples curves;
    native estimator grids still determine its reported probabilities and points.
    """
    settings = dict(log_half_span=10.0, nodes=201, check_tolerance=1e-7)
    settings.update(numerics or {})
    if masses is None:
        masses = settings.get("mass_grid")
    if masses is None:
        c = case.inputs.prior.get("center")
        if isinstance(c, bool) or not isinstance(c, (int, float)) or not math.isfinite(c) or c <= 0:
            raise ValueError("provide masses explicitly when the prior has no finite positive center")
        nodes, half = settings["nodes"], settings["log_half_span"]
        if type(nodes) is not int or nodes < 3 or not math.isfinite(half) or half <= 0:
            raise ValueError("diagnostic grid needs nodes >= 3 and a finite positive log_half_span")
        masses = np.exp(math.log(c) + np.linspace(-half, half, nodes))
    masses = np.asarray(masses, dtype=float)
    if (masses.ndim != 1 or masses.size < 3 or not np.all(np.isfinite(masses) & (masses > 0))
            or not np.all(np.diff(masses) > 0)):
        raise ValueError("inspection masses must be at least three increasing finite positive values")
    tolerance = settings["check_tolerance"]
    if isinstance(tolerance, bool) or not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("check_tolerance must be finite and positive")
    diagnostics = dict(mass=masses.copy(), readings=[])
    for reading in case.inputs.readings:
        diagnostics["readings"].append(dict(id=reading.id, **diagnostic(reading, masses)))
    for quantity in ("declared_cost", "residual_cost"):
        diagnostics["total_" + quantity] = np.sum([r[quantity] for r in diagnostics["readings"]], axis=0)
    evaluation = adapter(case.inputs, masses)
    context = Context(case.inputs, masses, evaluation, diagnostics, adapter, float(tolerance))
    return dict(format="met.inspection", format_version=1, case=case.to_dict(),
                investigator_truth=case.investigator_truth, evaluation=evaluation, diagnostics=diagnostics,
                checks=run_checks(context, checks), numerics=settings,
                assumptions=["Declared covariance is per reading; cross-reading errors are not represented.",
                             "Costs are squared corrections under declarations, not p-values.",
                             "Unconsumed case/reading fields are annotations and do not alter the calculation."],
                provenance={})


def _pack(value, arrays, path="result"):
    """Store arbitrary named arrays without making the report layout the data format."""
    if isinstance(value, np.ndarray):
        if value.dtype.hasobject:
            raise ValueError(f"object arrays cannot be exported safely: {path}")
        arrays[path] = value
        return {"npz": path, "shape": list(value.shape)}
    if isinstance(value, Mapping):
        if len({str(k) for k in value}) != len(value):
            raise ValueError(f"mapping keys collide after JSON string conversion: {path}")
        return {str(k): _pack(v, arrays, path + "/" + str(k).replace("~", "~0").replace("/", "~1"))
                for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_pack(v, arrays, path + "/" + str(i)) for i, v in enumerate(value)]
    if isinstance(value, np.generic):
        return _pack(value.item(), arrays, path)
    if isinstance(value, float) and not math.isfinite(value):
        return {"unavailable_scalar": str(value)}
    if isinstance(value, Path):
        return str(value)
    return value


def save_result(folder, result, case, *, reporter):
    folder.mkdir(parents=True, exist_ok=True)
    replay = dict(case.to_dict(include_truth=True), kind="inspection",
                  numerics=dict(result["numerics"], mass_grid=result["diagnostics"]["mass"].tolist()),
                  checks=[c["id"] for c in result["checks"]])
    (folder / "case.yaml").write_text(yaml.safe_dump(replay, sort_keys=False), encoding="utf-8")
    arrays = {}
    record = _pack(result, arrays)
    np.savez_compressed(folder / "arrays.npz", **arrays)
    record["figures_generated"] = bool(reporter(folder, result))
    (folder / "run.json").write_text(json.dumps(record, indent=2, allow_nan=False), encoding="utf-8")
    return record


def _callable_name(function):
    function = getattr(function, "func", function)
    return dict(module=getattr(function, "__module__", type(function).__module__),
                name=getattr(function, "__qualname__", type(function).__qualname__))


def _provenance(path, document, adapter, diagnostic, reporter, checks):
    from ..runner import code_version
    code = code_version()
    package = Path(__file__).parent
    code["inspection_file_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                      for p in sorted(package.glob("*.py"))}
    code["preset_sha256"] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    code["seed"] = document.get("seed")
    code["deterministic"] = True if adapter is ours29_adapter and diagnostic is profile else None
    code["determinism_scope"] = "Built-in numerical analysis; custom callable determinism is unknown. Figure identifiers may vary."
    code["callables"] = dict(adapter=_callable_name(adapter), diagnostic=_callable_name(diagnostic),
                             reporter=_callable_name(reporter),
                             checks=[_callable_name(getattr(c, "run", c)) for c in checks])
    code["dependencies"] = {}
    for name in ("numpy", "scipy", "PyYAML", "matplotlib"):
        try:
            code["dependencies"][name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            code["dependencies"][name] = None
    return code


def run(path, out_root="results", *, case_names=None, check_ids=None, progress=print,
        adapter=ours29_adapter, diagnostic=profile, reporter=None, checks=DEFAULT_CHECKS):
    """Inspect one case or a saved collection, serially; return output path and index.

    Case selection and check selection do not change the saved source. The index
    records the selections. A completed report can contain failed requirements;
    execution success does not mean all checks passed.
    """
    if reporter is None:
        from .report import write_report
        reporter = write_report
    document, cases = load_cases(path)
    if case_names is not None:
        requested = set(case_names)
        missing = requested - {case.name for case in cases}
        if missing:
            raise ValueError(f"unknown case names: {sorted(missing)}")
        cases = [case for case in cases if case.name in requested]
    if not cases:
        raise ValueError("select at least one case")
    checks = tuple(checks)
    identifier = lambda c: getattr(c, "id", getattr(c, "__name__", "custom"))
    if check_ids is None:
        check_ids = document.get("checks")
    if check_ids is not None:
        requested = set(check_ids)
        missing = requested - {identifier(check) for check in checks}
        if missing:
            raise ValueError(f"unknown check ids: {sorted(missing)}")
        checks = tuple(check for check in checks if identifier(check) in requested)
    name = re.sub(r"[^A-Za-z0-9_-]+", "-", str(document.get("name", "inspection"))).strip("-") or "inspection"
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    parent = Path(out_root) / name
    out = parent / stamp
    suffix = 1
    while out.exists():
        suffix += 1
        out = parent / f"{stamp}-{suffix}"
    code = _provenance(path, document, adapter, diagnostic, reporter, checks)
    out.mkdir(parents=True)
    (out / "preset.yaml").write_bytes(Path(path).read_bytes())
    bundle = dict(format="met.inspection.collection", format_version=1, provenance=code,
                  selected_cases=[c.name for c in cases], selected_checks=[identifier(c) for c in checks], cases=[])
    lines = [f"# {document.get('name', 'Case inspection')}", "", str(document.get("question", "")), "",
             "Execution can complete with failed or unresolved scientific checks.", "",
             "| Case | Pass | Fail | Unresolved | Inapplicable | Error | Report |",
             "|---|---:|---:|---:|---:|---:|---|"]
    links, finding_lines = [], []
    for index, case in enumerate(cases):
        result = inspect_case(case, adapter=adapter, diagnostic=diagnostic, checks=checks,
                              numerics=document.get("numerics"))
        result["provenance"] = code
        folder = out / f"case_{index:03d}"
        save_result(folder, result, case, reporter=reporter)
        counts = Counter(check["status"] for check in result["checks"])
        entry = dict(name=case.name, folder=folder.name, checks=dict(counts))
        bundle["cases"].append(entry)
        count_values = [counts[status] for status in ("pass", "fail", "unresolved", "inapplicable", "error")]
        lines.append("| " + " | ".join([case.name.replace("|", "\\|"), *map(str, count_values),
                     f"[Open]({folder.name}/report.html)"]) + " |")
        failures = [f'{c["id"]}: {c.get("observed", {})}' for c in result["checks"] if c["status"] in ("fail", "error")]
        links.append(f'<li><a href="{folder.name}/report.html">{escape(case.name)}</a> — '
                     + escape(", ".join(f"{n} {status}" for status, n in counts.items())) + "</li>")
        if progress:
            progress(f"{case.name}: " + ", ".join(f"{n} {status}" for status, n in counts.items()))
        if failures:
            finding_lines.extend(["", f"{case.name} findings:", *[f"- {failure}" for failure in failures], ""])
    lines.extend(finding_lines)
    (out / "run.json").write_text(json.dumps(bundle, indent=2, allow_nan=False), encoding="utf-8")
    (out / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (out / "report.html").write_text('<!doctype html><meta charset="utf-8"><title>Case inspection</title>'
        '<main style="font:16px/1.6 system-ui;max-width:1000px;margin:40px auto;padding:20px">'
        + f'<h1>{escape(str(document.get("name", "Case inspection")))}</h1>'
        + '<p>Inspect the individual checks and their stated authority. A completed run can contain failures and unresolved questions.</p><ul>'
        + "".join(links) + '</ul><p><a href="summary.md">Summary</a> · <a href="run.json">Provenance and index</a> · '
        '<a href="preset.yaml">Exact source preset</a></p></main>', encoding="utf-8")
    if progress:
        progress(f"wrote {out}")
    return out, bundle
