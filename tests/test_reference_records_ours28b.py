r"""ours28b (28 Sept; withdrawn the same night by the no-borrowing ruling, kept as record): what its
records say reference/noise_vs_signal/ours28b/ours28b.py computes, checked against that script.

Records (VD-docs, _dscn_mass_estimator/noise_vs_signal/ours28b/):
- ours28b_one_scale_per_noise_pair_2026-09-28.md: §1 (the law of mass and its 1F1 factors, carried
  state, reductions), §2 (the grids,
  the checks table and the comparison with ours27), §3 (the strong zero-reading property given up),
  §5.1 (the jump at equal noise pairs).
- README.md: what ours28b() returns.

Inputs are built exactly as ours28b_checks.py builds them (its construction is copied below, not
imported, because it prints on import): rng = default_rng(7), then the six readings R1 (for each
reading F is drawn before a), then the six rot3() calls for "each reading rotated in its own frame",
then the three for the per-axis readings R2. The precise-null reading and the eps loop are the
checks script's too; §5.1's eps = 0.999 is that same construction with the value §5.1 states.

Every number is asserted to the precision the record prints it with: k decimals means within
0.5*10^-k (plus 1e-12). The exception is the rotation row's readout changes, 9e-14 and 4e-13:
those are floating-point roundoff, whose digits belong to the machine that ran them. There the
property the record states (the readout does not change under rotation) is asserted to 1e-11,
and the digits are kept in comments: record 9e-14 and 4e-13, this machine +4.4e-14 and -5.4e-13
(the checks script's own run printed "diff +4.4e-14" and "diff -5.4e-13"). Equalities the
record states without a number ("identical", "law = starting law", "law equal to the 1D law")
are asserted to floating-point rounding (1e-12 relative).

The 1F1 form of each law factor (§1) is stated as mathematics, without a precision. The code
computes the law by quadrature over the 1601-point beta grid, and reproduces the closed form to
3.0e-9 relative (one pair) and 2.8e-9 (per-axis), measured on 30 Sept; asserted to 1e-8.

Speed: the flat limit (flat=True, ours27 on this input type) takes well under a second. Every call
with the scale law (flat=False) took about 90-290 s per noise pair on this (loaded) machine (the A(m) readout
integrates over the 1601-point beta grid at each of 801 mass nodes), so every test that needs one
is marked slow and skipped unless MET_SLOW is set. Results are cached per process, so tests that
share a call (the six-reading example) pay for it once.

Not reproducible from what was kept (not asserted):
- ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 checks table, flat-limit row:
  "max relative error $1.4\times10^{-4}$ (far tails)"
  The two 3D readings are not recorded, and the comparison against the ours25b closed form for
  A(m) is not in ours28b_checks.py.
- ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 checks table, width-calibration row:
  "94–95% on the regular counterexample at 100–1600 repetitions; 94–98% one-pair"
  These are the derivation record's §5 coverage runs, and ours27_one_scale_per_noise_pair_2026-09-28.md
  §9 says: "Scripts in Claude's workspace, not kept; numpy and scipy."
- ours28b_one_scale_per_noise_pair_2026-09-28.md, §6:
  "median swap product 0.978 before correcting the"
  That was the code before the midpoint fix, which was not kept (the kept code is the fixed one).
- ours28b_no_borrowing_between_readings_2026-09-28.md, §2:
  "1.65 under the middle and 1.80 → 1.65 under ours28b. It does not move flat or"
  and the ours28b column of its coverage table, e.g.
  "| one pair, every push $r=1$, 200 readings | 68.4 | 84.0 | 96.2 | 96.4 | – |"
  The 20 weak readings and the simulated series are not recorded, the law there was on a
  "401-point grid", not ours28b.py's, and its §7 says: "Scripts in Claude's workspace, not kept; numpy and scipy."
- The derivation record's own tables (ours27_one_scale_per_noise_pair_2026-09-28.md §2, §3, §5) are
  not produced by ours28b.py; its §9 says their scripts were not kept.
"""

import functools
import os
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.special import hyp1f1

_REF = Path(__file__).resolve().parent.parent / "reference" / "noise_vs_signal" / "ours28b"
sys.path.insert(0, str(_REF))
_dont_write = sys.dont_write_bytecode
sys.dont_write_bytecode = True          # leave reference/ untouched (no __pycache__)
try:
    import ours28b as ours28b_module    # noqa: E402
    from ours28b import ours28b         # noqa: E402
finally:
    sys.dont_write_bytecode = _dont_write

SLOW_REASON = "slow (each call with the scale law takes minutes per noise pair): set MET_SLOW=1 to run"


def _printed(decimals):
    """Half a unit in the last printed decimal, plus 1e-12."""
    return 0.5 * 10.0 ** -decimals + 1e-12


# ---- construction copied verbatim from ours28b_checks.py (same seed, same order of draws) ----
rng = np.random.default_rng(7)
def diag(sF, sa): return np.diag(np.concatenate([np.atleast_1d(sF), np.atleast_1d(sa)]).astype(float)**2)  # noqa: E302,E704
def rot3(): return np.linalg.qr(rng.normal(size=(3, 3)))[0]  # noqa: E704
def rotated(F, a, S, R):  # noqa: E302
    B = np.kron(np.eye(2), R); return R@np.asarray(F), R@np.asarray(a), B@S@B.T  # noqa: E702
def swapped(R_):
    out = []
    for F, a, S in R_:
        d = len(F); P = np.block([[np.zeros((d, d)), np.eye(d)], [np.eye(d), np.zeros((d, d))]])  # noqa: E702
        out.append((a, F, P@S@P.T))
    return out


# weak-ish one-pair data, 3D, several readings
S1 = diag([1, 1, 1], [0.5, 0.5, 0.5])
R1 = [(rng.normal(size=3)*1.0 + np.array([1.0, 0.4, -0.3]), rng.normal(size=3)*0.5 + np.array([0.5, 0.2, -0.15]), S1) for _ in range(6)]
# checks: rt = ours28b([rotated(F, a, S, rot3()) for F, a, S in R1], c=2.0)  (the next six draws)
R1_ROTATED = [rotated(F, a, S, rot3()) for F, a, S in R1]
S2 = diag([1, 1, 1], [0.2, 0.5, 1.25])
R2 = [(np.array([2.0, 1.0, 0.5]), np.array([0.9, 0.6, 0.1]), S2), (np.array([-1.0, 2.0, 1.2]), np.array([-0.3, 1.1, 0.7]), S2),
      (np.array([0.4, -0.8, 1.9]), np.array([0.3, -0.2, 1.1]), S2)]
# checks: rt2 = ours28b([rotated(F, a, S, rot3()) for F, a, S in R2], c=2.0)  (the three draws after that)
R2_ROTATED = [rotated(F, a, S, rot3()) for F, a, S in R2]


def _null_reading(eps):
    """The checks script's precise-null reading: 1D signal F=3, a=1 plus two 0/0 axes of noise eps."""
    S3 = np.diag([1.0, eps**2, eps**2, 1.0, eps**2, eps**2])
    return [([3.0, 0.0, 0.0], [1.0, 0.0, 0.0], S3)]


_CALLS = {
    # one pair, six 3D readings, c = 2
    "base": lambda: ours28b(R1, c=2.0),
    "base_flat": lambda: ours28b(R1, c=2.0, flat=True),
    "swap": lambda: ours28b(swapped(R1), c=0.5),
    "units": lambda: ours28b([(3*np.asarray(F), a, np.diag(np.r_[9*np.ones(3), np.ones(3)])@S) for F, a, S in R1], c=6.0),
    "rotated": lambda: ours28b(R1_ROTATED, c=2.0),
    "split": lambda: ours28b([([F[k]], [a[k]], np.diag([S[k, k], S[3+k, 3+k]])) for F, a, S in R1 for k in range(3)], c=2.0),
    "zero_appended": lambda: ours28b(R1 + [(np.zeros(3), np.zeros(3), S1)], c=2.0),
    "zero_alone": lambda: ours28b([(np.zeros(3), np.zeros(3), S1)]*3, c=2.0),
    # per-axis, three 3D readings, c = 2
    "per_axis": lambda: ours28b(R2, c=2.0),
    "per_axis_flat": lambda: ours28b(R2, c=2.0, flat=True),
    "per_axis_rotated": lambda: ours28b(R2_ROTATED, c=2.0),
    # precise nulls, c = 1
    "one_1d": lambda: ours28b([([3.0], [1.0], np.diag([1.0, 1.0]))], c=1.0),
}


@functools.lru_cache(maxsize=None)
def _run(name, eps=None):
    if name == "null":
        return ours28b(_null_reading(eps), c=1.0)
    return _CALLS[name]()


# ---- fast: the flat limit, the grids, the classes ----

def test_returns_what_the_readme_lists_on_the_recorded_grid():
    r = _run("base_flat")
    # README.md, "Here now": "c, flat=False)` returns the readout, the median, the 95% interval, the law,"
    # README.md, "Here now": "$A(m)$ and the per-pair classes. `flat=True` gives ours27 on the same input"
    # (checked here on the fast flat=True call; the flat=False call the README describes is checked in
    # test_returns_what_the_readme_lists_with_the_scale_law, which is slow)
    assert {"readout", "median", "lo", "hi", "law", "A", "classes"} <= set(r)
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2: "The law is on an 801-point grid in $\log m$ (half-width 9), the scale on a"
    assert len(r["lm"]) == 801 and len(r["law"]) == 801 and len(r["A"]) == 801
    assert abs((r["lm"][-1] - r["lm"][0]) / 2 - 9.0) <= 1e-12
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2: "1601-point grid in $\operatorname{logit}\beta$, and $A(m)$ by a 721-point"
    assert len(ours28b_module.L) == 1601 and len(ours28b_module.BETA) == 1601


def test_one_pair_is_one_class_and_its_state_adds():
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §1 Reductions: "- One noise pair: a single class, law"
    cls = ours28b_module.classes(ours28b_module.to_components(R1))
    assert len(cls) == 1
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §1 "Components and classes": "of all components, over the whole series, with the same pair. For a class,"
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §1 "Components and classes": "H_g(m)=\sum_{c\in g}A_c^2,\qquad M_g=\lvert g\rvert ."
    # (six 3D readings with one pair: 18 components in the one class)
    assert cls[0]["M"] == 18
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §1: "**Carried state.** $(P_g,Q_g,D_g,M_g)$ per noise pair, with the pair itself."
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §1: "Combining series is addition (T11, T12 per pair)."
    first = ours28b_module.classes(ours28b_module.to_components(R1[:2]))
    second = ours28b_module.classes(ours28b_module.to_components(R1[2:]))
    assert len(first) == len(second) == 1
    assert (first[0]["sF"], first[0]["sa"]) == (second[0]["sF"], second[0]["sa"]) == (cls[0]["sF"], cls[0]["sa"])
    for key in ("P", "Q", "D"):
        assert abs(first[0][key] + second[0][key] - cls[0][key]) <= 1e-12 * max(1.0, abs(cls[0][key]))
    assert first[0]["M"] + second[0]["M"] == cls[0]["M"]


def test_ours27_column_one_pair():
    r = _run("base_flat")
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 "Against ours27 on the same data", ours27 column:
    # "| one pair, six weak 3D readings, $c=2$ | readout 4.194, median 4.276, 95% [1.61, 31.4] | readout 3.266, median 3.224, 95% [0.22, 51.9] |"
    assert abs(r["readout"] - 4.194) <= _printed(3)
    assert abs(r["median"] - 4.276) <= _printed(3)
    assert abs(r["lo"] - 1.61) <= _printed(2)
    assert abs(r["hi"] - 31.4) <= _printed(1)


def test_ours27_column_per_axis():
    r = _run("per_axis_flat")
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 "Against ours27 on the same data", ours27 column:
    # "| per-axis, three 3D readings, $c=2$ | readout 1.923, median 1.971, 95% [0.84, 3.93] | readout 1.909, median 1.945, 95% [0.58, 4.34] |"
    assert abs(r["readout"] - 1.923) <= _printed(3)
    assert abs(r["median"] - 1.971) <= _printed(3)
    assert abs(r["lo"] - 0.84) <= _printed(2)
    assert abs(r["hi"] - 3.93) <= _printed(2)


# ---- slow: every call with the arcsine scale law ----

@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_returns_what_the_readme_lists_with_the_scale_law():
    r = _run("base")
    # README.md, "Here now": "c, flat=False)` returns the readout, the median, the 95% interval, the law,"
    # README.md, "Here now": "$A(m)$ and the per-pair classes. `flat=True` gives ours27 on the same input"
    assert {"readout", "median", "lo", "hi", "law", "A", "classes"} <= set(r)


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
@pytest.mark.parametrize("name", ["base", "per_axis"])
def test_law_factors_are_1F1(name):
    r = _run(name)
    c = 2.0
    m = np.exp(r["lm"])
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §1 Types: "- A declared starting law $\pi$, currently the half-Cauchy $c/(c^2+m^2)$."
    closed = c / (c * c + m * m) * m          # as weights on the uniform log-m grid
    for g in r["classes"]:
        # ours28b_one_scale_per_noise_pair_2026-09-28.md, §1 "Components and classes": "\theta_g(m)=\arctan\frac{m}{s_g},\qquad"
        th = np.arctan(m / g["s"])
        # ours28b_one_scale_per_noise_pair_2026-09-28.md, §1 Reductions: "$H=P\sin^2\theta+Q\cos^2\theta+2D\sin\theta\cos\theta$."
        H = g["P"] * np.sin(th) ** 2 + g["Q"] * np.cos(th) ** 2 + 2 * g["D"] * np.sin(th) * np.cos(th)
        # ours28b_one_scale_per_noise_pair_2026-09-28.md, §1 "The law of mass":
        # "and with the arcsine law each factor is ${}_1F_1(\tfrac12;\,1+\tfrac{M_g}2;\,\tfrac{H_g}2)$"
        # "up to a constant."
        closed = closed * hyp1f1(0.5, 1 + g["M"] / 2, H / 2)
    if name == "base":
        # ours28b_one_scale_per_noise_pair_2026-09-28.md, §1 Reductions: "- One noise pair: a single class, law"
        # ours28b_one_scale_per_noise_pair_2026-09-28.md, §1 Reductions: "$\pi(m)\,{}_1F_1(\tfrac12;1+\tfrac M2;\tfrac{H(\theta)}2)$ with"
        assert len(r["classes"]) == 1
    # The record prints no precision for this identity; the code's beta-grid quadrature reproduces
    # it to 3.0e-9 (base) and 2.8e-9 (per_axis) relative, measured 30 Sept. Asserted to 1e-8.
    assert np.allclose(r["law"] / r["law"].sum(), closed / closed.sum(), rtol=1e-8, atol=0)


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_ours28b_column_one_pair():
    r = _run("base")
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 "Against ours27 on the same data", ours28b column:
    # "| one pair, six weak 3D readings, $c=2$ | readout 4.194, median 4.276, 95% [1.61, 31.4] | readout 3.266, median 3.224, 95% [0.22, 51.9] |"
    assert abs(r["readout"] - 3.266) <= _printed(3)
    assert abs(r["median"] - 3.224) <= _printed(3)
    assert abs(r["lo"] - 0.22) <= _printed(2)
    assert abs(r["hi"] - 51.9) <= _printed(1)


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_ours28b_column_per_axis():
    r = _run("per_axis")
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 "Against ours27 on the same data", ours28b column:
    # "| per-axis, three 3D readings, $c=2$ | readout 1.923, median 1.971, 95% [0.84, 3.93] | readout 1.909, median 1.945, 95% [0.58, 4.34] |"
    assert abs(r["readout"] - 1.909) <= _printed(3)
    assert abs(r["median"] - 1.945) <= _printed(3)
    assert abs(r["lo"] - 0.58) <= _printed(2)
    assert abs(r["hi"] - 4.34) <= _printed(2)


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_swap_products():
    base, sw = _run("base"), _run("swap")
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 checks table, swap row:
    # "| swap (channels exchanged, $c\to1/c$) | one pair, six 3D readings | readout product 1.000000, median product 1.000000 |"
    assert abs(base["readout"] * sw["readout"] - 1.0) <= 5e-7
    assert abs(base["median"] * sw["median"] - 1.0) <= 5e-7


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_units_ratio():
    base, un = _run("base"), _run("units")
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 checks table, units row:
    # "| units (force ×3) | same | readout ratio 3.000000 |"
    assert abs(un["readout"] / base["readout"] - 3.0) <= _printed(6)


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_rotation_one_pair():
    diff = _run("rotated")["readout"] - _run("base")["readout"]
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 checks table, rotation row (first number, one pair):
    # "| each reading rotated in its own frame, with its covariance | same; and per-axis $\sigma_a=(0.2,0.5,1.25)$ | readout changes by $9\times10^{-14}$ and $4\times10^{-13}$ |"
    # invariant under rotation, to roundoff. Record 9e-14; this machine +4.35e-14 (the checks script
    # printed "diff +4.4e-14"). Roundoff digits are not asserted.
    assert abs(diff) <= 1e-11


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_rotation_per_axis():
    diff = _run("per_axis_rotated")["readout"] - _run("per_axis")["readout"]
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 checks table, rotation row (second number, per-axis):
    # "| each reading rotated in its own frame, with its covariance | same; and per-axis $\sigma_a=(0.2,0.5,1.25)$ | readout changes by $9\times10^{-14}$ and $4\times10^{-13}$ |"
    # invariant under rotation, to roundoff. Record 4e-13; this machine -5.36e-13 (the checks script
    # printed "diff -5.4e-13"). Roundoff digits are not asserted.
    assert abs(diff) <= 1e-11


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_axis_splitting_is_identical():
    base, sp = _run("base"), _run("split")
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 checks table, axis-splitting row:
    # "| axis splitting | six 3D readings as eighteen 1D readings | identical |"
    assert sp["readout"] == base["readout"]


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_axis_splitting_law_is_identical():
    base, sp = _run("base"), _run("split")
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 checks table, axis-splitting row:
    # "| axis splitting | six 3D readings as eighteen 1D readings | identical |"
    # (ours28b_checks.py prints only the readout difference for this row; "identical" is read here
    # as the whole result too: the law, its median and its 95% interval, to 1e-12 relative)
    assert np.array_equal(sp["lm"], base["lm"])
    assert np.allclose(sp["law"], base["law"], rtol=1e-12, atol=0)
    for key in ("median", "lo", "hi"):
        assert abs(sp[key] - base[key]) <= 1e-12 * abs(base[key])


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_zero_experiment_alone_gives_the_starting_law():
    r = _run("zero_alone")
    c = 2.0
    m = np.exp(r["lm"])
    start = c / (c * c + m * m) * m          # half-Cauchy c/(c^2+m^2), as weights on the uniform log-m grid
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 checks table, zero-experiment row:
    # "| zero experiment alone | three 3D zero readings | law = starting law (median 2.000000 at $c=2$); readout exactly $s$ |"
    assert np.allclose(r["law"] / r["law"].sum(), start / start.sum(), rtol=1e-12, atol=0)
    assert abs(r["median"] - 2.0) <= _printed(6)
    # s = sigma_F / sigma_a = 1 / 0.5 = 2 for S1; "exactly" asserted to the 1e-12 floor.
    assert abs(r["readout"] - 2.0) <= 1e-12


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_strong_zero_reading_property_is_given_up():
    base, z = _run("base"), _run("zero_appended")
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §3 "Given up":
    # "reading to the six-reading example moves the median 3.224 → 2.970 and the"
    assert abs(base["median"] - 3.224) <= _printed(3)
    assert abs(z["median"] - 2.970) <= _printed(3)
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §3 "Given up": "readout 3.266 → 3.054"
    assert abs(base["readout"] - 3.266) <= _printed(3)
    assert abs(z["readout"] - 3.054) <= _printed(3)


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_precise_null_1d_readout():
    r = _run("one_1d")
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 checks table, precise-nulls row ("against 2.1140" is the 1D reading):
    # "| precise nulls tend to known-empty | 1D reading $F=3$, $a=1$ against the same with two 0/0 axes of noise $\varepsilon$ | $\varepsilon$ = 0.1, 0.01, 0.001: law equal to the 1D law; readout 2.1126, 2.1140, 2.1140 against 2.1140 |"
    assert abs(r["readout"] - 2.1140) <= _printed(4)


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
@pytest.mark.parametrize("eps,readout", [(0.1, 2.1126), (0.01, 2.1140), (0.001, 2.1140)])
def test_precise_nulls_tend_to_known_empty(eps, readout):
    one, r = _run("one_1d"), _run("null", eps)
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §2 checks table, precise-nulls row:
    # "| precise nulls tend to known-empty | 1D reading $F=3$, $a=1$ against the same with two 0/0 axes of noise $\varepsilon$ | $\varepsilon$ = 0.1, 0.01, 0.001: law equal to the 1D law; readout 2.1126, 2.1140, 2.1140 against 2.1140 |"
    assert np.array_equal(r["lm"], one["lm"])
    assert np.allclose(r["law"], one["law"], rtol=1e-12, atol=0)
    assert abs(r["readout"] - readout) <= _printed(4)


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_jump_at_equal_noise_pairs_just_below():
    one, r = _run("one_1d"), _run("null", 0.999)
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §5.1 "A jump at equal noise pairs":
    # "0.999 beside the 1D reading above, the law is the 1D law (median 2.125). At"
    # The 1D law's own median, which does not depend on the eps = 0.999 construction:
    assert abs(one["median"] - 2.125) <= _printed(3)
    assert np.array_equal(r["lm"], one["lm"])
    assert np.allclose(r["law"], one["law"], rtol=1e-12, atol=0)
    assert abs(r["median"] - 2.125) <= _printed(3)


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_jump_at_equal_noise_pairs_exactly_equal():
    r = _run("null", 1.0)
    # ours28b_one_scale_per_noise_pair_2026-09-28.md, §5.1 "A jump at equal noise pairs":
    # "exactly 1.0 they join the signal axis's class (median 1.822). Pooling by"
    assert len(r["classes"]) == 1 and r["classes"][0]["M"] == 3
    assert abs(r["median"] - 1.822) <= _printed(3)
