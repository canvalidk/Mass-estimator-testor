r"""ours28 (28 Sept, the live line under the no-borrowing ruling): what its records say
reference/ours27/ours28_reference.py computes, checked against that script.

Records (read-only, in the records folder Newton-analysis/_dscn_mass_estimator/):
- ours27/ours28_reference_implementation_and_readout_2026-09-28.md ("the reference record"):
  §1 (numerics), §2 (the checks table and the paragraphs under it), §4 (worked comparison).
- ours27/ours28_c3_narrow_form_and_the_reading_structure_2026-09-28.md ("the C3 record"):
  the exact readouts quoted in §3 and §4, the law under the two reading structures in §4,
  and the split readout in §1.
- mass_estimator_lemmas_2026-09-29.md: L8 in §2 (the split readout), T13 in §3 (totality, the
  all-zero experiment's median), T25 in §3 (exact zero components leave H_i and the law
  unchanged), with its inputs in §7.
Code: reference/ours27/ours28_reference.py (byte-for-byte copy), imported by putting
reference/ours27 on sys.path as ours28_reference_checks.py does.

Inputs are built exactly as ours28_reference_checks.py builds them (its construction is
copied below, not imported, because it prints on import): rng = default_rng(7), then the six
readings R1 (for each reading, F's three normals before a's three), then L for S2, then R2,
then the three rot3() calls the checks script makes for its rotated R2 (section 3 of its
printout). No other draw is made from rng. The zero readings, the split into eighteen 1D
readings, the swap, units and rotate maps, the declared-subspace reading and the totality
readings are the checks script's own. The C3 record, §7, on its §3 examples:
  "those of the reference checks (seed 7)." (the line before ends "The examples are");
its §4 uses the reference record's "six-reading example (one pair, diagonal noise).". Its
eighteen 1D readings (§4) are taken as the checks script's split, the only such construction
kept; the C3 record's own script was not kept.

Tolerance: a number printed with k decimals is asserted to within 0.5*10^-k (plus 1e-12);
products and ratios printed as 1.000000 or 3.000000 to 5e-7. A residual the record gives as
the size of a disagreement ("agree to $4\times10^{-7}$ in the log", "unchanged ($8\times10^{-17}$, $4\times10^{-16}$)",
"readout $3\times10^{-13}$, median $4\times10^{-16}$") is asserted two-sided at the precision it
is printed with: the size computed as the checks script computes it must lie within half a unit
of the stated digit (4e-7 means 3.5e-7 to 4.5e-7). The 1e-12 is not added there, since it would
pass anything below 1e-12.
Residuals at floating-point roundoff (below 1e-11: the rotation row's "readout
$3\times10^{-13}$, median $4\times10^{-16}$" and the zero-reading row's "unchanged
($8\times10^{-17}$, $4\times10^{-16}$)") are different: their digits are the rounding of the
machine that ran them, and another machine rounds differently. There the property the record
states in words (invariant under rotation; the law unchanged) is asserted, to 1e-11, and the
record's digits and this machine's are kept side by side in the test's comments, not asserted.
Equalities the record states without a number (identical, the same to every digit, exactly,
difference 0.0, "median $c$") are asserted to 1e-12.

Roundoff digits that differ from the record's (not a disagreement about the math; recorded here,
not asserted): the rotation row's readout change is +1.66e-12 on this machine against the
record's 3e-13 (the kept output reference/expected/ours28_reference_checks.out prints
"readout diff +1.7e-12"), its median change 4.4e-16 against 4e-16; the zero-reading row's law
changes are 6.6e-17 and 3.7e-16 against 8e-17 and 4e-16.

Speed: each ours28 call is computed once per process and shared between tests. The same set of
calls took 43 s in one script on this machine. Earlier versions of this file took 107 s to 177 s
(slowest test 31 s to 35 s); this version took 225 s in one run while many other jobs shared
the machine (slowest test 37 s; a test's time is mostly the calls it is the first to need).
Nothing here needs MET_SLOW.

Not reproducible from what was kept (not asserted):
- Reference record, §2 checks table, flat-limit row:
  "| flat limit against ours27 | one pair, two 3D readings; and per-axis | readout 1.983766 = 1.983766; 1.917844 = 1.917844 |"
  The two 3D readings and the per-axis data are not in ours28_reference_checks.py and are not
  stated in any record.
- Reference record, §2, "Gaussian noise is still forced":
  "- Gaussian noise: 1.0000 at every angle;" and
  "- Student-$t_5$ noise (unit variance): 1.0158, 1.0121, 1.0064, 1.0017,"
  §6 says "arcsine-spaced scales; that script was not kept." Both lines come from that script
  (scipy quadrature over the push), not from ours28_reference.py.
- Reference record, §1 Numerics:
  "changes the readout by $4\times10^{-5}$ relative (3.732955 at full"
  The settings of the full-resolution run are not recorded. (For information only, not asserted:
  bins=1 with every other argument at its default gives 3.7329546 here, 4.4e-5 relative to
  3.7327894, in 95 s in one run and 173 s in another. Taking the record's full resolution to
  mean bins=1 with tnodes=241 kept is a guess, so it is not asserted.)
  The coarse-grained 3.732789 is asserted below.
- Reference record, §3 table (all six rows), e.g.
  "| 5 | 1 | 200 | 98.0 | 0.913 | 0.898 | 0.531 | 1.130 |"
  and the bullets under it, e.g. "law's width** (0.42/0.46, 0.098/0.093, 0.169/0.166, 0.084/0.084,".
  §6: "- `ours28_reference_sim.py N omega series seed` reproduces a row of §3. Rows"; the seed
  is an argument and no row's seed is recorded.
- Lemma record, §3, T29: "ours28, 3D, +50 null readings: 2.2524 → 2.2406; law unchanged to $10^{-14}$."
  §7: "T29: 1, 5 and 50 exact zero readings appended to one 1D and one 3D reading" and
  "scale weights (3000 draws per mass node). Scripts not kept." The 3D reading is not recorded.
- Lemma record, §7, T25's law and readout at other transverse noise, as a reproduction of that run:
  "transverse SDs at 10, 1, 0.1, 0.01 and 0.001 (law identical in every case)."
  The starting-law centre c for that run is not recorded, nor is the 7-point mass grid of
  "T25: $H_i$ on a 7-point mass". The run's readouts are not stated. T25's statements about
  H_i and about the law are asserted below as properties, on ours28's own mass grid at c = 1
  (the checks script's c for this reading); that is not a reproduction of the record's run.
- C3 record, §3 (candidates sampled, 4000 per trial mass):
  "| $1$ (C3) | 3.7374 (exact 3.732789) | 2.2476 (exact 2.248458) |" (the sampled 3.7374 and
  2.2476; the exact values are asserted below), "| $1+q$ | +0.34% | +0.08% |" and the rows
  after it, and "| change | +5.7% / −6.6% | +10.2% / −9.5% | +20.1% / −15.6% |".
  §7: "Scripts in Claude's workspace, not kept." The sampling seed and draw order are not
  recorded.
- C3 record, §3 first table, last row:
  "| law's 95% interval | [0.49, 64.7] | [0.75, 25.2] |"
  §7 puts §3's law on another grid than the reference record's 721 points:
  "- **§3.** The law is on a 361-point grid in $\log m$. Per reading, at each"
  and does not record its half-width. (For information only, not asserted: with 361 points at
  half-width 9, a guess, the six-reading interval is [0.486, 64.75], whose upper end does not
  print as 64.7.) The same two intervals are stated in the reference record §4 and asserted
  there, on the reference record's own grid.
- C3 record, §4 (quoted line split in two here):
  "or as the same eighteen numbers as eighteen 1D readings, it gives the same law"
  "(to $5\times10^{-15}$) and the same candidates, draw for draw. But $q$"
  The measure behind 5e-15 and the grid of the law are not recorded (§7 says only
  "- **§4.** `ours28(..., flat=True)` gives the law."). On ours28's default grid the answer
  depends on the measure: the largest change relative to the law's peak (the checks script's
  measure for a law change) is 4.9e-15, the largest pointwise relative change 5.7e-15, and the
  largest absolute change 9.0e-17. So the number is not asserted; that the law is the same is
  asserted below to 1e-12 in the strictest of these measures. The candidates are drawn by the
  script that was not kept.
- C3 record, §4 table (sampled, 8000 draws per trial mass, seed 1), e.g.
  "| $1+q$ | 4.2015 | 4.2020 |". The sampling script is not kept (§7), so the order of the
  draws is not recorded. The exact 4.193908 and the law under the two structures are
  asserted below.
- null_series_tilt_and_property_N_2026-09-29.md (noise_vs_signal/), §2 tables, e.g.
  "| 100 | 9.03 | 4.49 | 0.50 |"; §3, "$\chi^2_{dN}$. So the tilt grows like $\sqrt{dN}$: measured $1.34\sqrt N$";
  §8's check of §1, "$9.5\times10^{-4}$ ($N=1$) to $1.5\times10^{-2}$ ($N=1000$).".
  §8: "unmodified (`flat=True` for ours27), `nodes` 241–361, `lm_half` 5–6," and
  "`tnodes=121`. Scripts not kept." The exact settings and draw order are not recorded.

Left to the ours29 family (not asserted here): ours29_misfit_weight_and_bartlett_power_2026-09-30.md
§7, "reference's flat branch to $10^{-10}$ on its two test series.", which compares
ours28_reference.py with ours29_reference.py.
"""

import functools
import inspect
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.special import hyp1f1

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "reference" / "ours27"))

from ours28_reference import BETA, L1B, LB, log_arcsine, lse, ours28, reading_quantities  # noqa: E402

REF = "ours28_reference_implementation_and_readout_2026-09-28.md"
C3 = "ours28_c3_narrow_form_and_the_reading_structure_2026-09-28.md"
LEMMAS = "mass_estimator_lemmas_2026-09-29.md"


def tol(k):
    """Half a unit of the k-th decimal, plus 1e-12."""
    return 0.5 * 10.0 ** -k + 1e-12


# --- construction copied from reference/ours27/ours28_reference_checks.py (same order of draws) ---
rng = np.random.default_rng(7)
def rot3(): return np.linalg.qr(rng.normal(size=(3, 3)))[0]
def swap(rd):
    y, S = rd[0], rd[1]; d = len(y)//2; P = np.block([[np.zeros((d, d)), np.eye(d)], [np.eye(d), np.zeros((d, d))]])
    return (P@y, P@S@P.T) + tuple(rd[2:])
def rotate(rd, R):
    y, S = rd[0], rd[1]; B = np.kron(np.eye(2), R); out = (B@y, B@S@B.T)
    return out + ((R@rd[2],) if len(rd) > 2 else ())
def units(rd, lam):
    y, S = rd[0], rd[1]; d = len(y)//2; D = np.diag(np.r_[lam*np.ones(d), np.ones(d)]); return (D@y, D@S@D) + tuple(rd[2:])

S1 = np.diag([1.0]*3 + [0.25]*3)
R1 = [(np.r_[rng.normal(size=3) + [1.0, 0.4, -0.3], rng.normal(size=3)*0.5 + [0.5, 0.2, -0.15]], S1) for _ in range(6)]
L = rng.normal(size=(6, 6)); S2 = L@L.T/6 + 0.3*np.eye(6)          # full covariance, cross-channel correlation included
R2 = [(np.r_[[2.0, 1.0, 0.5], [0.9, 0.6, 0.1]], S2), (np.r_[[-1.0, 2.0, 1.2], [-0.3, 1.1, 0.7]], S2),
      (np.r_[[0.4, -0.8, 1.9], [0.3, -0.2, 1.1]], S2, np.array([[0.2], [-0.35], [0.915]])/np.linalg.norm([0.2, -0.35, 0.915]))]
R2_ROTATED = [rotate(r, rot3()) for r in R2]                       # the checks script's three rot3() calls, in order
SPLIT = [(np.r_[r[0][k], r[0][3+k]], np.diag([S1[k, k], S1[3+k, 3+k]])) for r in R1 for k in range(3)]
Su = np.eye(6); Y_SUB = np.r_[[4.0, 0, 0], [2.0, 0, 0]]
TOTALITY = [("perpendicular", [2.0, 0, 0], [0, 1.0, 0]), ("anti-parallel", [2.0, 0, 0], [-1.0, 0, 0]),
            ("force only", [3.0, 0, 0], [0, 0, 0]), ("acceleration only", [0, 0, 0], [3.0, 0, 0])]
# --- end of the copied construction ---

CASES = {
    "R1": lambda: ours28(R1, c=2.0),
    "R1 flat": lambda: ours28(R1, c=2.0, flat=True),
    "R1 swap": lambda: ours28([swap(r) for r in R1], c=0.5),
    "R1 + zero": lambda: ours28(R1 + [(np.zeros(6), S1)], c=2.0),
    "zero alone": lambda: ours28([(np.zeros(6), S1)]*3, c=2.0),
    "R1 split": lambda: ours28(SPLIT, c=2.0),
    "R1 split flat": lambda: ours28(SPLIT, c=2.0, flat=True, readout=False),   # only its law is used
    "R2": lambda: ours28(R2, c=2.0),
    "R2 flat": lambda: ours28(R2, c=2.0, flat=True),
    "R2 rotated": lambda: ours28(R2_ROTATED, c=2.0),
    "R2 swap": lambda: ours28([swap(r) for r in R2], c=0.5),
    "R2 units": lambda: ours28([units(r, 3.0) for r in R2], c=6.0),
    "R2 + zero": lambda: ours28(R2 + [(np.zeros(6), S2)], c=2.0),
    "subspace k=3": lambda: ours28([(Y_SUB, Su)], c=1.0),
    "subspace k=1": lambda: ours28([(Y_SUB, Su, np.array([[1.0], [0], [0]]))], c=1.0),
    "subspace 1D": lambda: ours28([(np.array([4.0, 2.0]), np.eye(2))], c=1.0),
}
for _label, _F, _a in TOTALITY:
    CASES["totality " + _label] = (lambda F=_F, a=_a: ours28([(np.r_[F, a], np.eye(6))], c=1.0))


@functools.lru_cache(maxsize=None)
def run(name):
    return CASES[name]()


def law_change(new, old):
    """The checks script's measure: largest change of the law, relative to the law's peak."""
    return np.abs(new["law"] - old["law"]).max() / old["law"].max()


# ---------------------------------------------------------------- reference record, §1


def test_numerics_grids_are_as_stated():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §1 The equation as implemented, Numerics:
    # "- $\beta$ on a 481-point logit grid; the law on a 721-point grid in $\log m$"
    # "  (half-width 9)."
    assert len(LB) == 481 and len(BETA) == 481
    assert np.array_equal(BETA, 1 / (1 + np.exp(-LB)))
    lm = run("R1")["lm"]
    assert len(lm) == 721
    assert abs((lm[-1] - np.log(2.0)) - 9.0) <= 1e-12 and abs((np.log(2.0) - lm[0]) - 9.0) <= 1e-12
    # same section: "- For the readout, the $\beta$ weights are coarse-grained into bins of 8 and"
    # "  the Laplace integral uses 241 points."
    defaults = inspect.signature(ours28).parameters
    assert defaults["bins"].default == 8
    assert defaults["tnodes"].default == 241
    assert defaults["nodes"].default == 721 and defaults["lm_half"].default == 9.0


def test_six_weak_readings_readout_coarse_grained():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §1, Numerics:
    # "  resolution, 3.732789 coarse-grained)."
    assert abs(run("R1")["readout"] - 3.732789) <= tol(6)


# ---------------------------------------------------------------- reference record, §2 checks table


def test_per_reading_evidence_against_1F1():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §2 Checks, table row:
    # "| per-reading evidence against ${}_1F_1(\tfrac12;1+\tfrac k2;\tfrac H2)$ | $H$ = 0.5, 7, 40, 200; $k$ = 3, 3, 1, 2 | agree to $4\times10^{-7}$ in the log |"
    # (computed as the checks script computes it; the worst of the four disagreements, which this
    # machine gives as 3.77e-7; asserted two-sided at the printed digit)
    worst = 0.0
    for H, k in [(0.5, 3), (7.0, 3), (40.0, 1), (200.0, 2)]:
        lev = lse(log_arcsine() + 0.5*k*L1B + 0.5*BETA*H, 0) - lse(log_arcsine() + 0.5*k*L1B, 0)
        worst = max(worst, abs(lev - np.log(hyp1f1(0.5, 1 + k/2, H/2))))
    assert abs(worst - 4e-7) <= 0.5e-7


def test_swap_c_to_1_over_c_six_weak_readings():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §2 Checks, table row:
    # "| swap, $c\to1/c$ | one pair, six weak 3D readings | readout product 1.000000, median product 1.000000 |"
    b1, s1 = run("R1"), run("R1 swap")
    assert abs(b1["readout"] * s1["readout"] - 1.0) <= 5e-7
    assert abs(b1["median"] * s1["median"] - 1.0) <= 5e-7


def test_swap_and_units_full_covariance():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §2 Checks, table row:
    # "| swap, units (force ×3) | full covariance with cross-channel correlation, one reading with a declared 1D subspace | products 1.000000; ratios 3.000000 |"
    b2, s2, u2 = run("R2"), run("R2 swap"), run("R2 units")
    assert abs(b2["readout"] * s2["readout"] - 1.0) <= 5e-7
    assert abs(b2["median"] * s2["median"] - 1.0) <= 5e-7
    assert abs(u2["readout"] / b2["readout"] - 3.0) <= 5e-7
    assert abs(u2["median"] / b2["median"] - 3.0) <= 5e-7


def test_rotation_with_covariance_and_subspace_median():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §2 Checks, table row:
    # "| each reading rotated with its covariance and its subspace | same | readout $3\times10^{-13}$, median $4\times10^{-16}$ |"
    # (the median part: invariant under rotation, to roundoff. Record 4e-16; this machine +4.44e-16,
    # as the kept output prints "median diff +4.4e-16". Roundoff digits are not asserted.)
    assert abs(run("R2 rotated")["median"] - run("R2")["median"]) <= 1e-11


def test_rotation_with_covariance_and_subspace_readout():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §2 Checks, table row:
    # "| each reading rotated with its covariance and its subspace | same | readout $3\times10^{-13}$, median $4\times10^{-16}$ |"
    # (the readout part: invariant under rotation, to roundoff. Record 3e-13; this machine +1.66e-12,
    # as the checks script's own kept output prints "readout diff +1.7e-12". Roundoff digits are not asserted.)
    assert abs(run("R2 rotated")["readout"] - run("R2")["readout"]) <= 1e-11


@pytest.mark.parametrize("base, bound", [("R1", 1e-11), ("R2", 1e-11)],
                         ids=["six weak readings", "full covariance"])
def test_exact_zero_reading_leaves_the_law_unchanged(base, bound):
    # ours28_reference_implementation_and_readout_2026-09-28.md, §2 Checks, table row:
    # "| exact zero reading appended: the law | both examples | unchanged ($8\times10^{-17}$, $4\times10^{-16}$) |"
    # The property "unchanged": the law's change, measured as the checks script measures it, is at
    # roundoff (asserted: no more than 1e-11). Record 8e-17 and 4e-16; this machine 6.6e-17 and 3.7e-16
    # (the kept output prints "law max change 6.6e-17"). Roundoff digits are not asserted.
    assert law_change(run(base + " + zero"), run(base)) <= bound


def test_exact_zero_reading_moves_the_readout():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §2 Checks, table row:
    # "| exact zero reading appended: the readout | six weak readings | 3.732789 → 3.722316: the point moves, the law does not, as in ours27 |"
    assert abs(run("R1")["readout"] - 3.732789) <= tol(6)
    assert abs(run("R1 + zero")["readout"] - 3.722316) <= tol(6)


def test_zero_experiment_alone():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §2 Checks, table row:
    # "| zero experiment alone | three 3D zero readings, $c=s=2$ | median 2.000000, readout 2.000000 |"
    z0 = run("zero alone")
    assert abs(z0["median"] - 2.0) <= tol(6)
    assert abs(z0["readout"] - 2.0) <= tol(6)
    # mass_estimator_lemmas_2026-09-29.md, §3 What follows (theorems), T13 (Totality), under ours28:
    # "the all-zero experiment gives median $c$"
    # (a symbolic equality, so to 1e-12; T13 names no data, so the reference record's all-zero
    # experiment above is used, c = 2; this machine gives 1.999999999999999)
    assert abs(z0["median"] - 2.0) <= 1e-12


def test_declared_subspace_numbers():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §2 Checks, table row:
    # "| declared subspace | 3D reading $F=(4,0,0)$, $a=(2,0,0)$, unit noise | $k=3$: readout 1.8836, 95% [0.552, 16.5]; $k=1$ along $x$: 1.9174, [0.664, 14.6]; the same numbers as a 1D reading: 1.9174, [0.664, 14.6], identical |"
    k3, k1, one = run("subspace k=3"), run("subspace k=1"), run("subspace 1D")
    assert abs(k3["readout"] - 1.8836) <= tol(4)
    assert abs(k3["lo"] - 0.552) <= tol(3) and abs(k3["hi"] - 16.5) <= tol(1)
    for r in (k1, one):
        assert abs(r["readout"] - 1.9174) <= tol(4)
        assert abs(r["lo"] - 0.664) <= tol(3) and abs(r["hi"] - 14.6) <= tol(1)


def test_declared_k1_is_the_1d_reading_to_every_digit():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §2 Checks, paragraph
    # "**The declared subspace is known-empty, exactly.**":
    # "$k=1$ along an axis gives the same law and readout as the 1D reading, to every"
    # "digit."
    # Also mass_estimator_lemmas_2026-09-29.md, §3 What follows (theorems), T25:
    # "A 3D reading declared $k=1$ equals the 1D reading exactly."
    k1, one = run("subspace k=1"), run("subspace 1D")
    assert np.abs(k1["law"] - one["law"]).max() <= 1e-12
    for key in ("readout", "median", "lo", "hi"):
        assert abs(k1[key] - one[key]) <= 1e-12


def test_totality():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §2 Checks, table row:
    # "| totality | perpendicular, anti-parallel, force only, acceleration only | readouts 1.3015, 1.3300, 2.2914, 0.4364, all finite; the last two multiply to 1.000 (swap) |"
    got = {label: run("totality " + label)["readout"] for label, _, _ in TOTALITY}
    for label, want in [("perpendicular", 1.3015), ("anti-parallel", 1.3300),
                        ("force only", 2.2914), ("acceleration only", 0.4364)]:
        assert np.isfinite(got[label])
        assert abs(got[label] - want) <= tol(4)
    assert abs(got["force only"] * got["acceleration only"] - 1.000) <= tol(3)
    # mass_estimator_lemmas_2026-09-29.md, §3 What follows (theorems), T13 (Totality):
    # "proper $\pi$ $\Rightarrow$ $\widehat m\in(0,\infty)$; the totality checks pass under ours28 and ours29 (perpendicular, anti-parallel, force only, acceleration only;"
    for label in got:
        assert 0.0 < got[label] < np.inf


def test_axis_split_without_the_scale_carried():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §2 Checks, paragraph
    # "**Axis splitting needs the scale carried.**":
    # "eighteen separate 1D readings each with its own scale, give readout 4.0807"
    # "against 3.7328. A reading's axes share its scale. Splitting is the same law"
    # Also ours28_c3_narrow_form_and_the_reading_structure_2026-09-28.md, §1 What was open:
    # "eighteen 1D readings without the scale carried gives readout 4.0807 against"
    # "3.7328. The no-borrowing ruling keeps it that way, because a reading is one"
    # and mass_estimator_lemmas_2026-09-29.md, §2 The properties of mass (premises), L8:
    # "The old strong form failed under ours28 (4.0807 against 3.7328)"
    assert abs(run("R1 split")["readout"] - 4.0807) <= tol(4)
    assert abs(run("R1")["readout"] - 3.7328) <= tol(4)


# ---------------------------------------------------------------- reference record, §4


def test_worked_comparison_six_weak_readings():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §4 Worked comparison, table
    # "| data | ours27 | ours28 |", row:
    # "| one pair, six weak 3D readings, $c=2$ | 4.1939, median 4.276, [1.61, 31.4] | 3.7328, median 4.014, [0.49, 64.7] |"
    # (ours27 is ours28(..., flat=True), as the checks script prints it)
    f1, b1 = run("R1 flat"), run("R1")
    assert abs(f1["readout"] - 4.1939) <= tol(4)
    assert abs(f1["median"] - 4.276) <= tol(3)
    assert abs(f1["lo"] - 1.61) <= tol(2) and abs(f1["hi"] - 31.4) <= tol(1)
    assert abs(b1["readout"] - 3.7328) <= tol(4)
    assert abs(b1["median"] - 4.014) <= tol(3)
    assert abs(b1["lo"] - 0.49) <= tol(2) and abs(b1["hi"] - 64.7) <= tol(1)


def test_worked_comparison_full_covariance():
    # ours28_reference_implementation_and_readout_2026-09-28.md, §4 Worked comparison, table
    # "| data | ours27 | ours28 |", row:
    # "| full covariance with cross-channel correlation, three readings (one declared $k=1$), $c=2$ | 2.2066, median 2.238, [0.99, 12.4] | 2.2485, median 2.360, [0.75, 25.2] |"
    f2, b2 = run("R2 flat"), run("R2")
    assert abs(f2["readout"] - 2.2066) <= tol(4)
    assert abs(f2["median"] - 2.238) <= tol(3)
    assert abs(f2["lo"] - 0.99) <= tol(2) and abs(f2["hi"] - 12.4) <= tol(1)
    assert abs(b2["readout"] - 2.2485) <= tol(4)
    assert abs(b2["median"] - 2.360) <= tol(3)
    assert abs(b2["lo"] - 0.75) <= tol(2) and abs(b2["hi"] - 25.2) <= tol(1)


# ---------------------------------------------------------------- the C3 record


def test_c3_record_exact_readouts():
    # ours28_c3_narrow_form_and_the_reading_structure_2026-09-28.md, §3 How much it matters, table
    # "| $h$ | one pair, six weak 3D readings | full covariance, three readings (one $k=1$) |", row:
    # "| $1$ (C3) | 3.7374 (exact 3.732789) | 2.2476 (exact 2.248458) |"
    # (the exact values only; the sampled ones are not reproducible, see the docstring)
    # Input: the reference record's defaults (721-point grid, half-width 9), which give these
    # "exact" values to every printed digit; they are the reference record's §1/§2 numbers. C3's
    # §7 puts its §3 law on "a 361-point grid in $\log m$" for the sampled column and does not
    # record that grid's half-width. (For information only: 361 points at half-width 9, a guess,
    # give 3.732794 and 2.248459, which do not print as these values.)
    assert abs(run("R1")["readout"] - 3.732789) <= tol(6)
    assert abs(run("R2")["readout"] - 2.248458) <= tol(6)
    # same record, §4 The premise that closes it (ours27 on the six-reading example):
    # "readout is 4.193908. So every non-constant $h$ gives one candidate law two"
    assert abs(run("R1 flat")["readout"] - 4.193908) <= tol(6)


def test_c3_record_same_law_as_six_3d_or_eighteen_1d_readings_under_ours27():
    # ours28_c3_narrow_form_and_the_reading_structure_2026-09-28.md, §4 The premise that closes it:
    # "six-reading example (one pair, diagonal noise). Entered as six 3D readings,"
    # "or as the same eighteen numbers as eighteen 1D readings, it gives the same law"
    # "(to $5\times10^{-15}$) and the same candidates, draw for draw. But $q$"
    # §7: "- **§4.** `ours28(..., flat=True)` gives the law."
    # The record does not state the measure behind 5e-15 or the grid of the law, so the number is
    # not asserted (see "Not reproducible" in the docstring). The property "the same law" is
    # asserted to 1e-12 (the file's rule for an equality), in the largest pointwise relative
    # change, which bounds the change relative to the peak and the absolute change; this machine
    # gives 5.7e-15, 4.9e-15 and 9.0e-17 for the three. Grid: ours28's default (721 points,
    # half-width 9). The eighteen readings are the checks script's split.
    new, old = run("R1 split flat")["law"], run("R1 flat")["law"]
    assert old.min() > 0.0
    assert (np.abs(new - old) / old).max() <= 1e-12


# ---------------------------------------------------------------- the lemma record, T25


@pytest.mark.parametrize("transverse_sd", [1.0, 10.0, 0.1, 0.01, 0.001])
def test_T25_exact_zero_components_leave_H_unchanged(transverse_sd):
    # mass_estimator_lemmas_2026-09-29.md, §3 What follows (theorems), T25:
    # "measured exact zeros on some axes of a signal reading leave $H_i$ exactly as without those axes, at every $m$, whatever their $\sigma$."
    # §7 Sources:
    # "grid for the 3D reading $F=(4,0,0)$, $a=(2,0,0)$ against the 1D reading"
    # "$(4,2)$, unit noise (difference 0.0); then the law and readout with the"
    # "transverse SDs at 10, 1, 0.1, 0.01 and 0.001 (law identical in every case)."
    # The record's 7-point mass grid is not recorded; the statement is "at every $m$", so it is
    # checked on the 721-point grid ours28 uses for this reading at c = 1. The transverse SD is
    # given to every transverse component, force and acceleration alike (the statement is for
    # "whatever their $\sigma$").
    m = np.exp(np.linspace(-9.0, 9.0, 721))
    s2 = transverse_sd ** 2
    H3 = reading_quantities(Y_SUB, np.diag([1.0, s2, s2, 1.0, s2, s2]), None, m)[0]
    H1 = reading_quantities(np.array([4.0, 2.0]), np.eye(2), None, m)[0]
    assert np.abs(H3 - H1).max() <= 1e-12


@pytest.mark.parametrize("transverse_sd", [10.0, 0.1, 0.01, 0.001])
def test_T25_law_identical_whatever_the_transverse_sd(transverse_sd):
    # mass_estimator_lemmas_2026-09-29.md, §3 What follows (theorems), T25:
    # "measured exact zeros on some axes of a signal reading leave $H_i$ exactly as without those axes, at every $m$, whatever their $\sigma$. The law changes only through $k_i$."
    # §7 Sources:
    # "transverse SDs at 10, 1, 0.1, 0.01 and 0.001 (law identical in every case)."
    # A property check, not a reproduction of the record's run: that run's c is not recorded.
    # Here c = 1, the checks script's c for this reading (F=(4,0,0), a=(2,0,0)), on ours28's
    # default grid; k = 3 throughout (no subspace declared), so by the statement only k could
    # change the law. The law at each SD is compared with the law at SD 1 ("identical": 1e-12).
    # The readouts are not stated identical and are not asserted.
    def law(sd):
        s2 = sd ** 2
        return ours28([(Y_SUB, np.diag([1.0, s2, s2, 1.0, s2, s2]))], c=1.0, readout=False)["law"]
    assert np.abs(law(transverse_sd) - law(1.0)).max() <= 1e-12
