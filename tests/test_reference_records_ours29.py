r"""ours29 (30 Sept): what its record says reference/noise_vs_signal/ours29_reference.py computes, checked
against that script.

Record (VD-docs, _dscn_mass_estimator/noise_vs_signal/):
- ours29_misfit_weight_and_bartlett_power_2026-09-30.md: §3 (the readout weight A), §4 (the eigenvalue
  dictionary: the max H_i, off-line misfit and explained-energy rows), §7 (exact checks), §8 (worlds, and
  world (a)'s Q), §10 (the readout), §14 (the reference's grids and step).

Inputs for §7 are built exactly as ours29_reference_checks.py builds them (its construction is copied
below, not imported, because it prints on import): rng = default_rng(7), then the six readings R1, then
L for S2, then the forty readings R3, then check 2's F and a (40 each, F first), then the forty rot3()
calls of check 3 (one per reading of R3, in order), then the thirty readings Rw of check 8. Check 1
imports the ours28 reference from reference/ours27, as the checks script does.

Inputs for §8 are drawn exactly as ours29_reference_worlds.py draws them (seed 30, one generator for the
five worlds in the script's order; its functions are copied below). The script draws a series and then
evaluates it; the evaluation (variants) uses no random numbers, so drawing every series first and
evaluating afterwards gives the same series. Each world takes minutes (4 to 15 minutes per world on this
machine, under load), so those tests are marked slow and skipped unless MET_SLOW is set. With MET_SLOW=1
all five reproduce the §8 table. World (a)'s impossible reading is a constant (no random draws), so its
Q is checked in a fast test, on that reading as the worlds script builds it (c = 1, as rescue() returns).

§4's eigenvalue-dictionary rows are checked on check 4's misaligned reading yq (unit noise, so
Sigma = I_2 (x) I_3 and the whitened data are F and a; c = 1, so the record's angle theta is the code's phi).

Every number is asserted to the precision the record prints it with: k decimals means within
0.5*10^-k (plus 1e-12); "1.000000"-style products and ratios within 5e-7. A bound the record states
(to $10^{-10}$, to $3\times10^{-8}$) is asserted as that bound as written (plus 1e-12). Where the record
says "exactly" about an equality (the ten zero readings, the declared subspace), the test asserts bitwise
equality (it holds); where it says "exactly" of a printed number ("the tilt is exactly **1.500 nats"), the
number is asserted to its printed precision (within 5e-4). Equalities the record states symbolically
(Q = 1/sigma^2 for the impossible pair; Q = 0 and g = 1 under k = 1; median c for the all-zero experiment;
the eigenvalue-dictionary rows lambda_+ and lambda_- and the Rayleigh quotient) are asserted to
floating-point rounding (1e-12, relative for large values); eta = 1 (one reading, the all-zero experiment)
holds bitwise and is asserted so. The swap's "$\eta$ unchanged" has no number of its own, so the swapped
series' eta is asserted to round to the record's eta = 0.403.

Record statements the code does not reproduce (kept as strict xfail, not changed):
- §7, rotation: "unchanged to $3\times10^{-8}$". The code gives eta -3.8e-8, readout -6.8e-8,
  median -3.6e-8 (the same as the checks script prints).
- §10: "$\eta$ never moves the law's peak". On §7's weak 3D series (c = 2) the law, pi included, peaks at
  ln m = 1.1181 with eta = 0.403 and at ln m = 1.2181 with eta = 1 (as a density in m: 2.443 against
  2.910). The likelihood part alone, exp(eta * sum_i l_i), peaks at ln m = 1.3181 for both eta.

Not reproducible from what was kept (not asserted; each quote is on one line of its own):
- §7 (quoted line split in two here):
  "unchanged; the refusals. The tester and the reference agree to $10^{-8}$ on"
  "null series in 1D, 2D and 3D."
  Which null series (their sizes, noise, seed) and which outputs were compared are not recorded, and no
  script for the comparison was kept (tests/test_ours29.py does not compare with the reference). So no
  tester-against-reference test is written here.
- §4: "Checks: $\operatorname{tr}M_i=T_i$ and gap $=\lvert C_i\rvert$ to $10^{-12}$;"
  The readings this was checked on are not recorded, and the check is not in ours29_reference_checks.py.
- §4: "the median eccentricity of Wishart$_2(k,I)$ is 1.00, 0.87, 0.71, 0.38 for"
  §14 says: "§6.3 and §4's Wishart medians: sandbox, 2000 and 20 000 draws; not kept."
- §6.1, e.g. "| 1000 | 0.06 | **+0.440** | +0.031 |" and "kg, $N=20$, signal 3: law width 1.00 in"
  §14 says: "§6.1: sandbox scripts (numpy, scipy), 1D, seed 11; not kept."
- §6.3: "Rayleigh(1) under pure noise at every $N$ (medians 1.49, 1.20, 1.18 at"
  Not kept (§14, as above).
- §5 (the table cell, quoted in part):
  "falls back to the prior (3 of 40 pure-noise 3D series of 12 readings in one check)"
  The tester's ours25 on series whose seed is not recorded; the check was not kept.
- §2: "evidence of small signals and moves (ours28b, 1.80 → 1.65 on one zero"
  An ours28b number; the series is not given in this record.
- §7: "reached as soon as $Q\ge\nu$ (on ours28's base it was approached, 1.48 at"
  The misfit weight on ours28's base is not in the kept code (reference/ours27/ours28_reference.py has
  no misfit weight).
- §13 item 4 (quoted line split in two here):
  "bounded (about 2 nats; coverage 75–83% against a mass drawn from $\pi$ in"
  "2D and 3D, §8 (d), and near 80% in the morning's 1D check)."
  The 75-83% is §8 (d), checked below; the "about 2 nats" and the morning's 1D check have no recorded
  input or script here.

Not checked here (not output of the reference scripts):
- The tester study in §9 and §10's readout quartiles and accuracy figures are the tester's (met,
  presets/ours29_noise_vs_signal.yaml, seed 20260930).
- §7's tester paragraph, from "In the tester, `tests/test_ours29.py` (13 tests, all pass with the existing"
  on: a claim about the tester's own tests, not about the reference.
- Closed-form numbers with no input from the reference: §10's
  "normal reading exceeds $\nu$ often ($P(\chi^2_1>1)=32\%$ in 2D," and
  "$P(\chi^2_2>2)=37\%$ in 3D), so about a third of ordinary readings are", and §6.3's
  "$N=10$, 100, 1000; Rayleigh 1.18), and a zero reading leaves it alone. It is".
- §4's row "| energy $T_i$ | $\operatorname{tr}M_i$ |" on its own: the reference has no T_i (it computes
  E_i, §3). The lambda_+ and lambda_- rows are checked below on E_i - Q_i and Q_i, whose sum is E_i.
"""

import functools
import inspect
import os
import sys
from pathlib import Path

import numpy as np
import pytest

_REF = Path(__file__).resolve().parent.parent / "reference"
sys.path.insert(0, str(_REF / "noise_vs_signal"))
sys.path.insert(0, str(_REF / "ours27"))
_dont_write = sys.dont_write_bytecode
sys.dont_write_bytecode = True          # leave reference/ untouched (no __pycache__)
try:
    import ours29_reference                                                          # noqa: E402
    from ours29_reference import _span_quantities, bartlett_eta, misfit_weight, ours29, reading_state  # noqa: E402
    from ours28_reference import ours28                                              # noqa: E402
finally:
    sys.dont_write_bytecode = _dont_write

SLOW_REASON = "slow (each world of section 8 takes 4-15 minutes): set MET_SLOW=1 to run"


def _printed(decimals):
    """Half a unit in the last printed decimal, plus 1e-12."""
    return 0.5 * 10.0 ** -decimals + 1e-12


# ---- construction copied verbatim from ours29_reference_checks.py (same seed, same order of draws) ----
rng = np.random.default_rng(7)


def rot3():
    return np.linalg.qr(rng.normal(size=(3, 3)))[0]


def swap(rd):
    y, S = rd[0], rd[1]; d = len(y)//2; P = np.block([[np.zeros((d, d)), np.eye(d)], [np.eye(d), np.zeros((d, d))]])  # noqa: E702
    return (P@y, P@S@P.T) + tuple(rd[2:])


def rotate(rd, R):
    y, S = rd[0], rd[1]; B = np.kron(np.eye(2), R); out = (B@y, B@S@B.T)  # noqa: E702
    return out + ((R@rd[2],) if len(rd) > 2 else ())


def units(rd, lam):
    y, S = rd[0], rd[1]; d = len(y)//2; D = np.diag(np.r_[lam*np.ones(d), np.ones(d)]); return (D@y, D@S@D) + tuple(rd[2:])  # noqa: E702


S1 = np.diag([1.0]*3 + [0.25]*3)
R1 = [(np.r_[rng.normal(size=3) + [1.0, 0.4, -0.3], rng.normal(size=3)*0.5 + [0.5, 0.2, -0.15]], S1) for _ in range(6)]
L = rng.normal(size=(6, 6)); S2 = L@L.T/6 + 0.3*np.eye(6)          # full covariance, cross-channel correlation included  # noqa: E702
R2 = [(np.r_[[2.0, 1.0, 0.5], [0.9, 0.6, 0.1]], S2), (np.r_[[-1.0, 2.0, 1.2], [-0.3, 1.1, 0.7]], S2),
      (np.r_[[0.4, -0.8, 1.9], [0.3, -0.2, 1.1]], S2, np.array([[0.2], [-0.35], [0.915]])/np.linalg.norm([0.2, -0.35, 0.915]))]
# a weak series where both tempering steps act: 3D pure noise under the full covariance S2
Lc = np.linalg.cholesky(S2)
R3 = [(Lc@rng.normal(size=6), S2) for _ in range(40)]
# check 2 (1D): the next 80 draws
F, a = rng.normal(size=40) + .2, rng.normal(size=40) + .1
# check 3: rr = ours29([rotate(r, rot3()) for r in R3], c=2.0)  (the next forty rot3() calls)
R3_ROTATED = [rotate(r, rot3()) for r in R3]
# check 4: one 3D reading along x (F = 4, a = 2, unit noise), transverse axes read 0; a misaligned one
Su = np.eye(6); y = np.r_[[4.0, 0, 0], [2.0, 0, 0]]  # noqa: E702
yq = np.r_[[4.0, 1.5, 0], [2.0, -1.0, 0]]
ALONG_X = np.array([[1.0], [0], [0]])
# check 5: one reading (aligned, weak, zero) under S2
ONE_READINGS = (np.r_[[1.0, 0, 0], [0, 1.0, 0]], np.r_[[0.3, -0.2, 0.1], [0.1, 0.2, -0.4]], np.zeros(6))
# check 7: awkward single readings
AWKWARD = [("perpendicular", [2.0, 0, 0], [0, 1.0, 0]), ("anti-parallel", [2.0, 0, 0], [-1.0, 0, 0]),
           ("force only", [3.0, 0, 0], [0, 0, 0]), ("acceleration only", [0, 0, 0], [3.0, 0, 0])]
# check 8: axis splitting (the draws after check 3's rotations)
Rw = [(np.r_[rng.normal(size=3)*1.0 + [0.3, 0, 0], rng.normal(size=3)*1.0 + [0.2, 0, 0]], np.eye(6)) for _ in range(30)]
sp = [(np.r_[r[0][k], r[0][3+k]], np.eye(2)) for r in Rw for k in range(3)]


# ---- the checks script's calls (cached: several tests share one) ----
_CALLS = {
    "R1_ours28_flat": lambda: ours28(R1, c=2.0, flat=True),
    "R1_off": lambda: ours29(R1, c=2.0, misfit=False, calibrate=False),
    "R2_ours28_flat": lambda: ours28(R2, c=2.0, flat=True),
    "R2_off": lambda: ours29(R2, c=2.0, misfit=False, calibrate=False),
    "b3": lambda: ours29(R3, c=2.0),
    "b3_raw": lambda: ours29(R3, c=2.0, calibrate=False),
    "b3_ours27": lambda: ours29(R3, c=2.0, misfit=False, calibrate=False),
    "rr": lambda: ours29(R3_ROTATED, c=2.0),
    "s3": lambda: ours29([swap(r) for r in R3], c=0.5),
    "u3": lambda: ours29([units(r, 3.0) for r in R3], c=6.0),
    "z3": lambda: ours29(R3 + [(np.zeros(6), S2)]*10, c=2.0),
    "k1": lambda: ours29([(y, Su, ALONG_X)], c=1.0),
    "one_1d": lambda: ours29([(np.array([4.0, 2.0]), np.eye(2))], c=1.0),
    "misaligned": lambda: ours29([(yq, Su)], c=1.0, readout=False),
    "misaligned_k1": lambda: ours29([(yq, Su, ALONG_X)], c=1.0, readout=False),
    "all_zero": lambda: ours29([(np.zeros(6), np.eye(6))]*3, c=2.0),
    "a3": lambda: ours29(Rw, c=1.0, readout=False),
    "a1": lambda: ours29(sp, c=1.0, readout=False),
}


@functools.lru_cache(maxsize=None)
def _call(key):
    return _CALLS[key]()


# ---------------------------------------------------------------- §7 exact checks

@pytest.mark.parametrize("name", ["R1", "R2"])
def test_both_steps_off_is_ours27(name):
    base, off = _call(f"{name}_ours28_flat"), _call(f"{name}_off")
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", Both steps off is ours27:
    # "reference's flat branch to $10^{-10}$ on its two test series."
    assert abs(off["readout"] - base["readout"]) <= 1e-10
    assert abs(off["median"] - base["median"]) <= 1e-10


def test_1d_eta_against_the_carried_state_closed_form():
    # the checks script's check 2, verbatim
    Ci = (a**2 - F**2) + 2j*F*a; C = Ci.sum(); th = np.angle(C)/2  # noqa: E702
    Jc = (np.sum(abs(Ci)**2) - np.real(np.exp(4j*th)*np.conj(np.sum(Ci**2))))/8
    e = bartlett_eta([(np.array([F[i], a[i]]), np.eye(2)) for i in range(40)], 1.0, [1.0]*40)[0]
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks":
    # "- **1D $\eta$** against the carried-state closed form: $7\times10^{-9}$."
    assert abs(abs(e - min(1, abs(C)/Jc)) - 7e-9) <= 0.5e-9 + 1e-12


def test_weak_3d_series_eta_and_misfit_weights():
    b3 = _call("b3")
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", A weak 3D series:
    # "$\eta=0.403$; 4 of 40 readings weighted below 1 (min 0.62). Quoted 95%"
    assert abs(b3["eta"] - 0.403) <= _printed(3)
    assert sum(g < 1 for g in b3["g"]) == 4
    assert abs(min(b3["g"]) - 0.62) <= _printed(2)


def test_weak_3d_series_intervals():
    b3, raw, o27 = _call("b3"), _call("b3_raw"), _call("b3_ours27")
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", A weak 3D series:
    # "interval $[1.30,57.8]$ against the raw law's $[1.87,34.0]$ and ours27's"
    assert abs(b3["lo"] - 1.30) <= _printed(2)
    assert abs(b3["hi"] - 57.8) <= _printed(1)
    assert abs(raw["lo"] - 1.87) <= _printed(2)
    assert abs(raw["hi"] - 34.0) <= _printed(1)
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", A weak 3D series:
    # "$[1.63,15.0]$ ($c=2$)."
    assert abs(o27["lo"] - 1.63) <= _printed(2)
    assert abs(o27["hi"] - 15.0) <= _printed(1)


ROTATION_BOUND = 3e-8 + 1e-12       # the bound as the record writes it, "to $3\times10^{-8}$"


@pytest.mark.xfail(strict=True, reason="record states eta unchanged to 3x10^-8 under rotation, code gives eta diff -3.8e-8")
def test_weak_3d_series_rotation_eta():
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", A weak 3D series:
    # "- each reading rotated with its covariance: $\eta$, readout, median" / "unchanged to $3\times10^{-8}$;"
    # (quoted line split in two here)
    assert abs(_call("rr")["eta"] - _call("b3")["eta"]) <= ROTATION_BOUND


@pytest.mark.xfail(strict=True, reason="record states readout unchanged to 3x10^-8 under rotation, code gives readout diff -6.8e-8")
def test_weak_3d_series_rotation_readout():
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", A weak 3D series:
    # "- each reading rotated with its covariance: $\eta$, readout, median" / "unchanged to $3\times10^{-8}$;"
    # (quoted line split in two here)
    assert abs(_call("rr")["readout"] - _call("b3")["readout"]) <= ROTATION_BOUND


@pytest.mark.xfail(strict=True, reason="record states median unchanged to 3x10^-8 under rotation, code gives median diff -3.6e-8")
def test_weak_3d_series_rotation_median():
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", A weak 3D series:
    # "- each reading rotated with its covariance: $\eta$, readout, median" / "unchanged to $3\times10^{-8}$;"
    # (quoted line split in two here)
    assert abs(_call("rr")["median"] - _call("b3")["median"]) <= ROTATION_BOUND


def test_weak_3d_series_swap():
    b3, s3 = _call("b3"), _call("s3")
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", A weak 3D series:
    # "- swap: readout and median products 1.000000, $\eta$ unchanged;"
    assert abs(b3["readout"] * s3["readout"] - 1.0) <= 5e-7
    assert abs(b3["median"] * s3["median"] - 1.0) <= 5e-7
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", A weak 3D series:
    # "$\eta=0.403$; 4 of 40 readings weighted below 1 (min 0.62). Quoted 95%" (and "$\eta$ unchanged" above)
    assert abs(s3["eta"] - 0.403) <= _printed(3)


def test_weak_3d_series_units():
    b3, u3 = _call("b3"), _call("u3")
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", A weak 3D series:
    # "- units (force $\times3$): ratios 3.000000;"
    assert abs(u3["readout"] / b3["readout"] - 3.0) <= 5e-7
    assert abs(u3["median"] / b3["median"] - 3.0) <= 5e-7


def test_weak_3d_series_ten_zero_readings():
    b3, z3 = _call("b3"), _call("z3")
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", A weak 3D series:
    # "- ten exact zero readings appended: $\eta$ and the law unchanged exactly;"
    assert z3["eta"] == b3["eta"]
    assert np.array_equal(z3["law"], b3["law"])
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", A weak 3D series:
    # "**the readout moves, 3.5983 → 3.6017** (T29, below)."
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §10 "The readout", Z (O1) is untouched:
    # "3.5983 → 3.6017: the stacked cloud counts the zero readings (T29). ours29"
    assert abs(b3["readout"] - 3.5983) <= _printed(4)
    assert abs(z3["readout"] - 3.6017) <= _printed(4)


def test_declared_subspace_k1_is_the_1d_reading():
    k1, one = _call("k1"), _call("one_1d")
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", Declared subspace:
    # "- **Declared subspace:** a 3D reading declared $k=1$ gives the 1D reading's" / "law and readout exactly;"
    # (quoted line split in two here)
    assert np.array_equal(k1["law"], one["law"])
    assert k1["readout"] == one["readout"]


def test_declared_subspace_misaligned_reading():
    mis, mis_k1 = _call("misaligned"), _call("misaligned_k1")
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", Declared subspace:
    # "law and readout exactly; a misaligned 3D reading ($Q=2.34$, $g=0.853$)"
    assert abs(mis["Q"][0] - 2.34) <= _printed(2)
    assert abs(mis["g"][0] - 0.853) <= _printed(3)
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", Declared subspace:
    # "declared $k=1$ has $Q=0$, $g=1$."
    assert abs(mis_k1["Q"][0]) <= 1e-12
    assert abs(mis_k1["g"][0] - 1.0) <= 1e-12


def test_off_line_misfit_is_the_smaller_gram_eigenvalue():
    """On check 4's misaligned reading (unit noise, so Sigma = I_2 (x) I_3 and the whitened data are F and a)."""
    x, yy = yq[:3], yq[3:]
    M = np.array([[x @ x, x @ yy], [x @ yy, yy @ yy]])
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §4 "The eigenvalue dictionary", table:
    # "| off-line misfit $Q_i$ (L21's input) | $\lambda_-$ |"
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §4, Correction to the conversation:
    # "11:18. It is $\lambda_-=Q_{\min}=\tfrac12(T_i-\lvert C_i\rvert)$, as the"
    assert _call("misaligned")["Q"][0] == pytest.approx(np.linalg.eigvalsh(M)[0], rel=1e-12, abs=1e-12)


def _gram_yq():
    """Check 4's misaligned reading yq as the 2x2 Gram matrix M_i of §4 (unit noise: x = F, y = a)."""
    x, yy = yq[:3], yq[3:]
    return np.array([[x @ x, x @ yy], [x @ yy, yy @ yy]])


def test_energy_along_the_best_line_is_the_larger_gram_eigenvalue():
    """On check 4's misaligned reading: max H_i = E_i - Q_i (reading_state gives Q = E - max H, here > 0)."""
    st = reading_state(yq, Su, None, 1.0)
    assert st["Q"] > 0
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §4 "The eigenvalue dictionary", table:
    # "| energy along the best line, $\max H_i$ | $\lambda_+$ |"
    assert st["E"] - st["Q"] == pytest.approx(np.linalg.eigvalsh(_gram_yq())[1], rel=1e-12, abs=1e-12)


def test_explained_energy_is_the_rayleigh_quotient():
    """On check 4's misaligned reading (c = 1, so theta = phi), on the code's own 2001-point reading grid."""
    phi = np.linspace(-np.pi/2, np.pi/2, 2001, endpoint=False)
    H = _span_quantities(yq, Su, None, 1.0, phi)[0]
    u = np.stack([np.sin(phi), np.cos(phi)], axis=1)
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §4 "The eigenvalue dictionary", table:
    # "| explained energy $H_i(\theta)$ | the Rayleigh quotient $\mathbf u(\theta)^\top M_i\mathbf u(\theta)$, $\mathbf u=(\sin\theta,\cos\theta)$ |"
    assert H == pytest.approx(np.einsum("ni,ij,nj->n", u, _gram_yq(), u), rel=1e-12, abs=1e-12)


@pytest.mark.parametrize("which", range(3), ids=["aligned", "weak", "zero"])
def test_one_reading_gives_eta_one(which):
    got = ours29([(ONE_READINGS[which], S2)], c=2.0, readout=False)["eta"]
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", One reading:
    # "- **One reading:** $\eta=1$ for aligned, weak and zero readings."
    assert got == 1.0


@pytest.mark.parametrize("sg", [1.0, 0.5, 0.1, 0.005])
def test_impossible_pair(sg):
    # the checks script's check 6, verbatim
    r = ours29([(np.r_[[2.0, 0], [0, 1.0]], sg**2*np.eye(4))], c=2.0, readout=False, lm_half=14, nodes=4001)
    lpi = np.log(2.0/(4.0+np.exp(2*r['lm']))) + r['lm']
    lr = np.log(r['law']) - lpi
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", The impossible pair:
    # "2D): $Q=1/\sigma^2$, and the tilt is exactly **1.500 nats at every"
    assert r["Q"][0] == pytest.approx(1 / sg**2, rel=1e-12, abs=1e-12)
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", The impossible pair:
    # "2D): $Q=1/\sigma^2$, and the tilt is exactly **1.500 nats at every" / "$\sigma\le1$**, the limit $\nu a/(1-a)$."
    # (quoted line split in two here)
    assert abs(np.ptp(lr) - 1.500) <= _printed(3)


@pytest.mark.parametrize("label,F_,a_", AWKWARD, ids=[w[0] for w in AWKWARD])
def test_totality_awkward_readings_are_finite(label, F_, a_):
    r = ours29([(np.r_[F_, a_], np.eye(6))], c=1.0)
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", Totality:
    # "acceleration-only readings all finite; the all-zero experiment gives"
    assert np.isfinite(r["readout"]) and np.isfinite(r["median"])


def test_totality_all_zero_experiment():
    r0 = _call("all_zero")
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", Totality:
    # "median $c$ and $\eta=1$."   (c = 2 in the checks script)
    assert r0["median"] == pytest.approx(2.0, rel=1e-12, abs=1e-12)
    assert r0["eta"] == 1.0


def test_axis_splitting():
    a3, a1 = _call("a3"), _call("a1")
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", Axis splitting:
    # "- **Axis splitting:** 30 3D readings against their 90 axes, $\eta$ 0.458"
    assert abs(a3["eta"] - 0.458) <= _printed(3)
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", Axis splitting:
    # "against 0.445, intervals slightly different, as the rewritten L8 intends."
    assert abs(a1["eta"] - 0.445) <= _printed(3)
    assert (a3["lo"], a3["hi"]) != (a1["lo"], a1["hi"])


# ---------------------------------------------------------------- §3 and §10: the readout

def test_neither_misfit_nor_eta_touches_the_readout_weight():
    """A(m) of the §7 weak 3D series under ours29, its raw law and ours27, on the nodes all three evaluate."""
    b3, raw, o27 = _call("b3"), _call("b3_raw"), _call("b3_ours27")
    kept = (b3["A"] > 0) & (raw["A"] > 0) & (o27["A"] > 0)
    assert kept.any()
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §3 "The equation", Candidates and readout:
    # "$A(m)=E[\lVert\mathbf a^\ast_{\rm stack}\rVert\mid m]$. Neither $g_i$ nor" / "$\eta$ touches $A$."
    # (quoted line split in two here)
    assert np.array_equal(b3["A"][kept], raw["A"][kept])
    assert np.array_equal(b3["A"][kept], o27["A"][kept])


def test_eta_moves_the_readout():
    b3, raw = _call("b3"), _call("b3_raw")
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §10 "The readout":
    # "unchanged. $\eta$ never moves the law's peak, but it moves the readout,"
    assert b3["readout"] != raw["readout"]


@pytest.mark.xfail(strict=True, reason="record states eta never moves the law's peak, code gives the law's peak at "
                                       "ln m = 1.1181 with eta = 0.403 and at ln m = 1.2181 with eta = 1 "
                                       "(section 7's weak 3D series, c = 2; pi included; the likelihood part "
                                       "alone peaks at ln m = 1.3181 for both)")
def test_eta_never_moves_the_law_peak():
    b3, raw = _call("b3"), _call("b3_raw")
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §10 "The readout":
    # "unchanged. $\eta$ never moves the law's peak, but it moves the readout,"
    assert b3["lm"][np.argmax(b3["law"])] == raw["lm"][np.argmax(raw["law"])]


# ---------------------------------------------------------------- §14: the reference's grids and step

def test_law_grid_and_difference_step():
    b3 = _call("b3")
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §14 "Checks and reproduction":
    # "$h=10^{-4}$ in $\phi$; law on 721 points over $\ln c\pm9$."
    assert len(b3["lm"]) == 721
    assert b3["lm"][0] == pytest.approx(np.log(2.0) - 9, rel=1e-12, abs=1e-12)
    assert b3["lm"][-1] == pytest.approx(np.log(2.0) + 9, rel=1e-12, abs=1e-12)
    assert inspect.signature(bartlett_eta).parameters["h"].default == 1e-4


def test_projective_line_grids():
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §14 "Checks and reproduction":
    # "Projective line on 2001 (reading) and 4001 (series) points, refined by"
    assert "np.linspace(-np.pi/2, np.pi/2, 2001, endpoint=False)" in inspect.getsource(ours29_reference.reading_state)
    assert "np.linspace(-np.pi/2, np.pi/2, 4001, endpoint=False)" in inspect.getsource(ours29_reference.bartlett_eta)


# ---------------------------------------------------------------- §8 worlds (seed 30)

# ---- copied verbatim from ours29_reference_worlds.py ----
LM = np.linspace(-12, 12, 2401); M = np.exp(LM)  # noqa: E702
PHI = np.linspace(-np.pi/2, np.pi/2, 4001, endpoint=False)


def variants(readings, c):
    """Coverage-relevant CDFs for ours27, ours29_law, ours29 from one pass over the readings."""
    phi_law = np.arctan(M/c)
    Hl = np.array([_span_quantities(r[0], r[1], None, c, phi_law)[0] for r in readings])      # (N, G)
    Hf = np.array([_span_quantities(r[0], r[1], None, c, PHI)[0] for r in readings])          # (N, K)
    st = [reading_state(r[0], r[1], None, c) for r in readings]
    g = np.array([misfit_weight(s["Q"], s["nu"]) for s in st])
    logpi = np.log(c/(c*c + M*M)) + LM
    tot = (0.5*g[:, None]*Hf).sum(0); j = int(np.argmax(tot)); h = PHI[1] - PHI[0]  # noqa: E702
    ell = 0.5*g[:, None]*Hf[:, [j-1, j, (j+1) % len(PHI)]]
    psi = (ell[:, 2] - ell[:, 0])/(2*h); curv = -np.sum(ell[:, 2] - 2*ell[:, 1] + ell[:, 0])/h**2  # noqa: E702
    J = np.sum(psi**2); eta = 1.0 if J <= 1e-300 or curv <= 0 else min(1.0, curv/J)  # noqa: E702
    out = {}
    for name, gg, ee in (("ours27", np.ones_like(g), 1.0), ("ours29_law", g, 1.0), ("ours29", g, eta)):
        lp = logpi + ee*(0.5*gg[:, None]*Hl).sum(0); p = np.exp(lp - lp.max()); p /= p.sum()  # noqa: E702
        out[name] = np.cumsum(p)
    return out, eta


def covered(cdf, m0):
    k = np.searchsorted(M, m0); return 0.025 < cdf[min(k, len(M)-1)] < 0.975  # noqa: E702


@functools.lru_cache(maxsize=None)
def _world_series():
    """Every world's series, drawn as the worlds script draws them: one generator (seed 30), worlds in its order."""
    rng = np.random.default_rng(30)

    def reading(m0, snr, d, sf=1.0, sa=1.0):
        u = rng.normal(size=d); u /= np.linalg.norm(u)  # noqa: E702
        a_star = snr*sa*u
        return (np.r_[m0*a_star + sf*rng.normal(size=d), a_star + sa*rng.normal(size=d)], np.diag([sf**2]*d + [sa**2]*d))

    def rescue():
        m0 = np.exp(rng.uniform(-1, 1))
        R = [reading(m0, 3.0, 2) for _ in range(20)]
        R.append((np.r_[[2.0, 0.0], [0.0, 1.0]], 0.01*np.eye(4)))               # sigma = 0.1: Q = 100
        return R, m0, 1.0

    def strong_among_nulls():
        m0 = np.exp(rng.uniform(-1, 1))
        return [reading(m0, 6.0, 3)] + [reading(m0, 0.0, 3) for _ in range(100)], m0, 1.0

    def own_setups():
        m0 = np.exp(rng.uniform(-1, 1))
        R = [reading(m0, rng.choice([0, 0, 0.5, 1, 3, 10]), 3, 10**rng.uniform(-2, 2), 10**rng.uniform(-2, 2)) for _ in range(40)]
        return R, m0, 1.0

    def null_from_prior(d):
        def make():
            m0 = np.tan(rng.uniform(0, np.pi/2))                                   # half-Cauchy, c = 1
            return [reading(m0, 0.0, d) for _ in range(100)], m0, 1.0
        return make

    # the script's order: run(...) for (a) 150, (b) 100, (c) 100, (d) 2D 100, (d) 3D 100
    plan = [("a", rescue, 150), ("b", strong_among_nulls, 100), ("c", own_setups, 100),
            ("d2", null_from_prior(2), 100), ("d3", null_from_prior(3), 100)]
    return {key: [make() for _ in range(series)] for key, make, series in plan}


@functools.lru_cache(maxsize=None)
def _world(key):
    """Coverage (percent) of ours27, ours29_law, ours29 and the median eta, as the worlds script's run() computes them."""
    hits = {k: [] for k in ("ours27", "ours29_law", "ours29")}; etas = []  # noqa: E702
    for R, m0, c in _world_series()[key]:
        cdfs, eta = variants(R, c); etas.append(eta)  # noqa: E702
        for k in hits: hits[k].append(covered(cdfs[k], m0))  # noqa: E701
    return {k: 100*np.mean(v) for k, v in hits.items()}, float(np.median(etas))


def test_world_a_impossible_reading_q():
    """World (a)'s impossible reading, as rescue() appends it, read by reading_state as variants() reads it (c = 1)."""
    impossible = (np.r_[[2.0, 0.0], [0.0, 1.0]], 0.01*np.eye(4))               # sigma = 0.1: Q = 100 (the script's line)
    Q = reading_state(impossible[0], impossible[1], None, 1.0)["Q"]
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §8 "Worlds the tester cannot generate", table:
    # "| (a) rescue: 20 readings of signal 3 plus one impossible reading ($\sigma=0.1$, $Q=100$), 2D | 3.3 | 95.3 | 95.3 | 1.00 | 150 |"
    assert abs(Q - 100) <= _printed(0)
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §7 "Exact checks", The impossible pair:
    # "2D): $Q=1/\sigma^2$, and the tilt is exactly **1.500 nats at every"   (sigma = 0.1 here)
    assert Q == pytest.approx(1 / 0.1**2, rel=1e-12, abs=1e-12)


def _assert_world(key, ours27, law, quoted, eta):
    cov, med_eta = _world(key)
    assert abs(cov["ours27"] - ours27) <= _printed(1)
    assert abs(cov["ours29_law"] - law) <= _printed(1)
    assert abs(cov["ours29"] - quoted) <= _printed(1)
    assert abs(med_eta - eta) <= _printed(2)


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_world_a_rescue():
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §8 "Worlds the tester cannot generate", table:
    # "| (a) rescue: 20 readings of signal 3 plus one impossible reading ($\sigma=0.1$, $Q=100$), 2D | 3.3 | 95.3 | 95.3 | 1.00 | 150 |"
    _assert_world("a", 3.3, 95.3, 95.3, 1.00)
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §8, (a):
    # "already covers. One impossible reading otherwise takes ours27 to 3%."
    assert abs(_world("a")[0]["ours27"] - 3) <= _printed(0)


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_world_b_strong_among_nulls():
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §8 "Worlds the tester cannot generate", table:
    # "| (b) one reading of signal 6 among 100 pure-noise readings, 3D | 65.0 | 62.0 | 100.0 | 0.19 | 100 |"
    _assert_world("b", 65.0, 62.0, 100.0, 0.19)


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_world_c_own_setups():
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §8 "Worlds the tester cannot generate", table:
    # "| (c) every reading its own instrument ($\sigma_F,\sigma_a$ each $10^{U(-2,2)}$) and its own signal (0, 0, 0.5, 1, 3 or 10), 3D, $N=40$ | 92.0 | 93.0 | 94.0 | 1.00 | 100 |"
    _assert_world("c", 92.0, 93.0, 94.0, 1.00)


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_world_d_null_from_prior_2d():
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §8 "Worlds the tester cannot generate", table:
    # "| (d) pure noise, planted mass drawn from $\pi$ itself, 2D, $N=100$ | 53.0 | 50.0 | 83.0 | 0.18 | 100 |"
    _assert_world("d2", 53.0, 50.0, 83.0, 0.18)


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("MET_SLOW"), reason=SLOW_REASON)
def test_world_d_null_from_prior_3d():
    # ours29_misfit_weight_and_bartlett_power_2026-09-30.md, §8 "Worlds the tester cannot generate", table:
    # "| (d) the same, 3D | 42.0 | 40.0 | 75.0 | 0.17 | 100 |"
    _assert_world("d3", 42.0, 40.0, 75.0, 0.17)
