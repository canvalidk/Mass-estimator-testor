"""tools/cloud_readout_vs_median.py: its cloud readout is T10's, and its series are ours33_coverage.py's (O2, G6)."""
import math

import numpy as np
import pytest

from met.full import ours30
from met.ours33 import ours33
from tools import cloud_readout_vs_median as study
from tools import ours33_coverage


@pytest.mark.parametrize("signal", [0.0, 1.0, 5.0])
def test_the_cloud_helper_is_ours30s_own_readout_on_ours30s_law(signal):
    rng = np.random.default_rng(1)
    rd = study.series(rng, 2.0, 3, signal, 6)
    out = ours30(rd, 1.0, readout=True)
    assert study.cloud(out["lm"], out["law"], rd) == pytest.approx(out["readout"], rel=1e-12)


def test_the_cell_index_and_seeds_are_ours33_coverages():
    for k in study.KS:
        for n in study.NS:
            for r in study.SIGNALS:
                assert study.index("gaussian", k, n, r) == ours33_coverage.index("gaussian", k, n, r)


def test_at_the_all_zero_experiment_the_cloud_on_ours33_reads_s_when_c_is_s():
    """c = s = 1: the cloud's pull is toward s, so it agrees with the starting law's centre only here (O2, 5 Nov)."""
    rd = study.series(np.random.default_rng(0), 1.0, 3, 0.0, 1)
    rd = [type(r)(0 * r.force, 0 * r.acceleration, r.cov) for r in rd]
    out = ours33(rd, 1.0, study.NU)
    assert out["median"] == pytest.approx(1.0, rel=1e-6)
    assert study.cloud(out["lm"], out["law"], rd) == pytest.approx(1.0, rel=1e-3)


def test_one_small_cell_runs_and_reports_every_readout():
    k, r, n, m0, res = study.cell((3, 5.0, 20, 3, 1.0))
    assert set(res) == set(study.KEYS)
    assert all(math.isfinite(e) and math.isfinite(b) for e, b in res.values())
