"""O14, 4 Nov: eta's target is the shared size law's shrinkage (met/size_law.py; tools/eta_vs_size_law.py).

What is pinned: the algebra of ours28b's profile shrinkage; that eta of the flat law and that shrinkage
agree in population on long series; that the shrinkage fails Z where eta does not; and the finite-N
defect of eta as a measurement (small N, weak signal). Nothing here changes the equation."""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from met.full import Reading, isotropic  # noqa: E402
from met.size_law import beta_bar, beta_hat, beta_profile, explained_energy, flat_eta, latent_count, oracle  # noqa: E402
from tools.eta_vs_size_law import cell, world  # noqa: E402


def test_beta_hat_is_the_profile_of_the_shared_size_law():
    # argmax_beta (M/2) log(1 - beta) + beta H^2/2 = 1 - M/H^2, the EB / James-Stein form
    for h2, m in ((40.0, 20.0), (300.0, 60.0), (7.0, 3.0)):
        assert beta_profile(h2, m) == pytest.approx(1.0 - m / h2, abs=1e-8)
    assert beta_profile(10.0, 20.0) == 0.0          # less explained energy than noise: no signal


def test_beta_bar_agrees_with_beta_hat_when_the_series_is_long():
    # the arcsine posterior mean (ours28b integrated) and the profile agree to O(1/M)
    assert beta_bar(600.0, 300.0) == pytest.approx(0.5, abs=0.01)
    assert beta_bar(1500.0, 300.0) == pytest.approx(0.8, abs=0.005)


def test_eta_tends_to_the_size_laws_shrinkage_in_population():
    # long series: mean eta within a few points of mean beta_hat, which tracks Delta/(Delta + M)
    for k, snr in ((1, 2.0), (3, 2.0), (3, 1.0)):
        _, _, _, eta, hat, bar, orc, _ = cell((900 + k, k, 100, snr, 40))
        assert abs(hat - orc) < 0.03, (k, snr, hat, orc)
        assert abs(eta - hat) < 0.07, (k, snr, eta, hat)
        assert abs(bar - hat) < 0.03


def test_eta_at_small_n_sits_far_above_the_target():
    # the finite-N defect, measured: with five weak readings eta barely tempers while the target is ~0.2
    _, _, _, eta, hat, _, orc, _ = cell((901, 1, 5, 0.5, 60))
    assert eta > 0.8 and orc < 0.35 and hat < 0.5


def test_the_size_laws_shrinkage_fails_z_where_eta_does_not():
    rng = np.random.default_rng(4)
    readings, _ = world(rng, 3, 6, 2.0)
    eta0, _, _, phi0 = flat_eta(readings, 1.0)
    hat0 = beta_hat(readings, 1.0, phi0)
    padded = readings + [isotropic(np.zeros(3), np.zeros(3), 1.0, 1.0)]
    eta1, _, _, phi1 = flat_eta(padded, 1.0)
    assert eta1 == eta0 and phi1 == phi0                       # Z for eta: exact-zero terms are 0
    assert latent_count(padded) == latent_count(readings) + 3
    assert explained_energy(padded, 1.0, phi1) == pytest.approx(explained_energy(readings, 1.0, phi0))
    assert beta_hat(padded, 1.0, phi1) < hat0                   # the count moved, the energy did not
