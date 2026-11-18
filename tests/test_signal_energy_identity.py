"""The signal-energy identity at the true mass (Claude, 18 Nov 2026; problems/O14.md, 18 Nov; tools/signal_energy_identity.py)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import signal_energy_identity as sei  # noqa: E402


@pytest.mark.parametrize("k", (1, 3))
@pytest.mark.parametrize("signal", (0.0, 1.0, 5.0))
def test_energy_minus_the_noise_share_is_unbiased_for_signal_energy_at_every_signal_size(k, signal):
    rng = np.random.default_rng(1 + k + int(10 * signal))
    f, al, phi0 = sei.series(rng, k, [signal], 200000)
    _, _, C = sei.parts(f, al, phi0)
    assert (C - 2 * k).mean() == pytest.approx(signal ** 2, abs=0.05 + 0.02 * signal ** 2)


@pytest.mark.parametrize("k", (1, 3))
@pytest.mark.parametrize("signals", ([0.0] * 5, [1.0] * 5, [5.0] * 5, [30.0] + [0.0] * 20, [10.0] * 2 + [0.5] * 20))
def test_the_weighted_leftover_has_its_declared_noise_mean_at_the_true_mass_whatever_the_pushes(k, signals):
    # sum (C - 2k) K and sum k (C - 2k + 2) have the same mean at the true mass: no push size enters, no reading is counted
    rng = np.random.default_rng(1811)
    assert sei.cell(rng, k, signals, 60000)["true"] == pytest.approx(1.0, abs=0.03)


def test_at_the_fitted_line_the_fit_removes_leftover_and_the_correction_is_not_exact():
    # recorded, not a requirement: at weak signal the fit takes far more than the guessed correction allows for
    rng = np.random.default_rng(7)
    weak = sei.cell(rng, 1, [0.0] * 5, 40000)
    strong = sei.cell(rng, 1, [5.0] * 20, 40000)
    assert weak["fitted"] < 0.5
    assert strong["fitted"] == pytest.approx(1.0, abs=0.03)


@pytest.mark.parametrize("k", (1, 3))
@pytest.mark.parametrize("signals", ([0.0] * 5, [5.0] * 5, [30.0] + [0.0] * 20, [5.0] * 2 + [0.0] * 200))
def test_covariance_form_exact_at_the_true_mass_and_never_above_it_at_its_own_least_line(k, signals):
    # E[(C - 2) K] = k E[C] at the true mass for every push and every k; and the weighted leftover's minimum over
    # lines is at most its value at the true line on every series, so its mean is at most the expected side's
    rng = np.random.default_rng(18 + k)
    c = sei.fit_cell(rng, k, signals, 20000 if len(signals) < 60 else 3000)
    assert c["true"] == pytest.approx(1.0, abs=0.04)
    assert c["never_above_true"] and c["own"] < 1.0


def test_read_at_the_gaussian_best_line_the_mean_can_exceed_the_expected_side():
    # two readings at signal 5 among 200 noise-only ones: the Gaussian best line is pulled by the noise-only readings
    rng = np.random.default_rng(5)
    c = sei.fit_cell(rng, 1, [5.0] * 2 + [0.0] * 200, 4000)
    assert c["gauss"] > 1.05 and c["own"] < 1.0


def test_the_first_order_fit_correction_at_strong_signal_and_what_it_leaves_at_weak_signal():
    rng = np.random.default_rng(9)
    assert sei.fit_cell(rng, 1, [5.0] * 20, 20000)["own_corrected"] == pytest.approx(1.0, abs=0.03)
    assert sei.fit_cell(rng, 3, [5.0] * 20, 20000)["own_corrected"] == pytest.approx(1.0, abs=0.03)
    assert sei.fit_cell(rng, 1, [5.0] * 2, 40000)["own_corrected"] == pytest.approx(0.87, abs=0.04)   # two readings
    assert sei.fit_cell(rng, 1, [0.0] * 20, 20000)["own_corrected"] < 0.7                              # no signal
