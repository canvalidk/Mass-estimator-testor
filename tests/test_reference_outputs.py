"""The reference scripts in reference/ print what they printed when they were carried over.

Each script under reference/ is a byte-for-byte copy of the one in the VD-docs records
(see reference/README.md for paths and SHA-256). Its checks print numbers from fixed
seeds. reference/expected/ holds the output each one gave when carried over (30 Sept).
This test runs the script again, from its own folder exactly as the records say to run it,
and compares its output with that: every piece of text must be identical, and every
number must agree to REL_TOL (relative) or ABS_TOL (absolute), which absorbs last-digit
floating-point differences between machines and nothing more.

This checks repeatability, not correctness: what each record claims about these numbers
is checked in tests/test_reference_records.py.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "reference"
EXPECTED = REFERENCE / "expected"

REL_TOL = 1e-6
ABS_TOL = 1e-9

# (script path under reference/, expected-output file, slow)
SCRIPTS = [
    ("ours27/ours28_reference_checks.py", "ours28_reference_checks.out", False),
    ("noise_vs_signal/ours29_reference_checks.py", "ours29_reference_checks.out", False),
    ("mass_estimator_equation/example.py", "example.out", False),
    # the ours29 record's §8 study (seed 30); about 20 minutes on an idle machine
    ("noise_vs_signal/ours29_reference_worlds.py", "ours29_reference_worlds.out", True),
    # ours28b's checks (seed 7); about 23 minutes on an idle machine
    ("noise_vs_signal/ours28b/ours28b_checks.py", "ours28b_checks.out", True),
]

_NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?|[-+]?(?:nan|inf)")


def _tokens(text):
    """Split output into alternating text and numbers: [(False, text), (True, number), ...]."""
    out, pos = [], 0
    for m in _NUMBER.finditer(text):
        out.append((False, text[pos:m.start()]))
        out.append((True, m.group()))
        pos = m.end()
    out.append((False, text[pos:]))
    return out


def _same_number(a, b):
    x, y = float(a), float(b)
    if x != x or y != y:
        return x != x and y != y
    if x == y:
        return True
    return abs(x - y) <= max(ABS_TOL, REL_TOL * max(abs(x), abs(y)))


def compare_outputs(got, expected):
    """None if the outputs agree, else a description of the first difference."""
    g, e = _tokens(got.replace("\r\n", "\n")), _tokens(expected.replace("\r\n", "\n"))
    if len(g) != len(e):
        return f"different structure: {len(g)} tokens against {len(e)} expected"
    for i, ((gn, gv), (en, ev)) in enumerate(zip(g, e)):
        if gn != en or (not gn and gv != ev):
            return f"text differs at token {i}: {gv!r} against expected {ev!r}"
        if gn and not _same_number(gv, ev):
            return f"number differs at token {i}: {gv} against expected {ev}"
    return None


def run_script(relpath):
    script = REFERENCE / relpath
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0")
    done = subprocess.run([sys.executable, "-u", script.name], cwd=script.parent, env=env,
                          capture_output=True, text=True, timeout=7200)
    assert done.returncode == 0, done.stderr[-2000:]
    return done.stdout


@pytest.mark.parametrize("relpath,expected,slow", SCRIPTS, ids=[s[0] for s in SCRIPTS])
def test_reference_script_output_is_unchanged(relpath, expected, slow, request):
    if slow and not os.environ.get("MET_SLOW"):
        pytest.skip("slow (over 5 minutes): set MET_SLOW=1 to run")
    got = run_script(relpath)
    want = (EXPECTED / expected).read_text(encoding="utf-8")
    problem = compare_outputs(got, want)
    assert problem is None, f"{relpath}: {problem}"


def test_reference_files_are_the_carried_over_bytes():
    """Every file listed in reference/SHA256SUMS still has the hash it was carried over with.
    The math in reference/ is never edited here; a change belongs in the records first."""
    import hashlib
    lines = (REFERENCE / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
    assert lines, "SHA256SUMS is empty"
    for line in lines:
        digest, name = line.split(maxsplit=1)
        path = REFERENCE / name.lstrip("*")
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, f"{name} differs from the carried-over copy"
    listed = {line.split(maxsplit=1)[1].lstrip("*") for line in lines}
    present = {p.relative_to(REFERENCE).as_posix() for p in REFERENCE.rglob("*.py")}
    assert present == listed, f"unlisted or missing reference scripts: {sorted(present ^ listed)}"


def test_the_comparison_is_strict_where_it_should_be():
    assert compare_outputs("readout 1.000000 x", "readout 1.000000 x") is None
    assert compare_outputs("readout 1.0000001 x", "readout 1.000000 x") is None
    assert compare_outputs("readout 1.00001 x", "readout 1.000000 x") is not None
    assert compare_outputs("median 1.0 x", "readout 1.0 x") is not None
    assert compare_outputs("a 1 2", "a 1 2 3") is not None
    assert compare_outputs("eta nan", "eta nan") is None
