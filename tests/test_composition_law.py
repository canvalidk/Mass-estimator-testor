"""L5 and L14 at the level of the law (restated 28 Oct, R-2026-10-28-1), on independent new
measurements: tools/composition_coverage.py. The full study (400 replicates per size, seed
2810) is recorded in problems/O2.md, 28 Oct; this runs a few replicates of it."""
import importlib.util
import os

import numpy as np
import pytest

_spec = importlib.util.spec_from_file_location(
    "composition_coverage", os.path.join(os.path.dirname(__file__), "..", "tools", "composition_coverage.py"))
cc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cc)


@pytest.fixture(scope="module")
def strong():
    return [cc._one((seed, 2.0)) for seed in range(12)]


def test_computed_and_measured_laws_agree_at_strong_signal(strong):
    """Clause (c): the computed law and the composite's own law describe one quantity."""
    for name in ("L5", "L14"):
        assert sum(r[name]["overlap"] for r in strong) >= 11
        assert sum(r[name]["own_covers"] for r in strong) >= 10
        assert sum(r[name]["computed_covers"] for r in strong) >= 10


def test_no_point_meets_composition_on_data(strong):
    """Neither the cloud readout nor the median meets L5 or L14 at finite noise; neither is
    systematically better (full study: L14 5.84% vs 5.85% strong, 91.4% vs 90.5% weak)."""
    for name in ("L5", "L14"):
        cloud = np.median([r[name]["defect_cloud"] for r in strong])
        median = np.median([r[name]["defect_median"] for r in strong])
        assert cloud > 1e-3 and median > 1e-3
        assert 0.8 < cloud / median < 1.25
