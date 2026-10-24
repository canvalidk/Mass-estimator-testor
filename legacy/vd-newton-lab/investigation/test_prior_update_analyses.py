"""Integration checks for analyses that rescore saved trial archives."""

from contextlib import redirect_stdout
import io
import math
from pathlib import Path
import unittest
from uuid import uuid4

import numpy as np

from archive_rescoring import failures, iter_cells
from experiment_config import run_config
from paired_dominance import analyze as dominance
from prior_weight_study import parameter_grid, summarize, sweep
from sequential_update import analyze_run
from switch_rules import analyze as switches


class PriorUpdateAnalysesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[3] / ".tools" / "mass_prior_update_tests" / uuid4().hex
        base = {"samples": 64, "batch_size": 32, "save_trials": True,
                "training_seed": 11, "heldout_seed": 12}
        with redirect_stdout(io.StringIO()):
            run_config({**base, "scenario": {"type": "pairs", "pairs": [[1, 1], [4, 1], [3, 3]]}},
                       root / "pairs")
            run_config({**base, "scenario": {"type": "mass_excitation", "masses": [1.],
                                             "acceleration_snrs": [0., 2.]}}, root / "ladder")
        cls.pairs, cls.ladder = root / "pairs", root / "ladder"

    def test_dominance_counts_and_pairing(self):
        result = dominance(self.pairs, factors=(1.5,))
        self.assertEqual(result["family_size"], 3 * 1 * len(result["compared_rules"]))
        self.assertAlmostEqual(result["bonferroni_z"], 3.0, delta=.5)
        entry, arrays = next(iter_cells(self.pairs))
        truth = entry["scenario"]["true_mass"]
        manual = (failures(arrays["point__flat_joint"], truth, 1.5)
                  - failures(arrays["point__norm_ratio"], truth, 1.5)).mean()
        cell = result["cells"][0]["factors"]["1.5"]["comparisons"]["norm_ratio"]
        self.assertAlmostEqual(cell["focal_minus_other_failure"], manual, places=12)
        excluded = dominance(self.pairs, factors=(1.5,), exclude=("tube_joint",))
        self.assertNotIn("tube_joint", excluded["compared_rules"])

    def test_switch_limits_and_oracle(self):
        result = switches(self.pairs, thresholds=(0., 1e9), fallbacks=("flat_median",), factors=(2.,))
        for cell in result["cells"]:
            block = cell["factors"]["2"]
            _, arrays = next((e, a) for e, a in iter_cells(self.pairs, [cell["cell"]]))
            truth = cell["true_mass"]
            norm = 1 - failures(arrays["point__norm_ratio"], truth, 2.).mean()
            median = 1 - failures(arrays["point__flat_median"], truth, 2.).mean()
            self.assertAlmostEqual(block["switch"]["flat_median"]["0"]["success"], norm, places=12)
            self.assertAlmostEqual(block["switch"]["flat_median"]["1e+09"]["success"], median, places=12)
            self.assertAlmostEqual(block["oracle"]["flat_median"]["success"], max(norm, median), places=12)

    def test_prior_sweep_reproduces_flat_and_scores_regret(self):
        parameters = parameter_grid("sech", lams=(0., .35, 8.))
        cells = sweep([self.pairs], parameters, factors=(1.5, 2.))
        for cell in cells:
            self.assertLess(cell["max_relative_deviation_lam0_vs_archived_flat"], 1e-10)
            self.assertAlmostEqual(cell["success"]["1.5"]["weighted"]["sech|centre=1|lam=0"],
                                   cell["success"]["1.5"]["archived"]["flat_joint"], places=12)
        table = summarize(cells, parameters, (1.5, 2.))
        flat = table["sech|centre=1|lam=0"]
        self.assertGreaterEqual(flat["1.5"]["max_regret"], 0.)
        self.assertEqual(flat["2"]["worst_cell_gain_over_flat"], 0.)
        oracle = parameter_grid("lognormal_factor", widths=(.01,), centre_factors=(1.,))
        exact = sweep([self.pairs], oracle, factors=(1.25,))
        for cell in exact:
            self.assertGreater(cell["success"]["1.25"]["weighted"]["lognormal_factor|centre_factor=1|width=0.01"], .95)
        halves = [sweep([self.pairs], parameters, factors=(2.,), half=h)[0]["readings"] for h in ("first", "second")]
        self.assertEqual(halves, [32, 32])

    def test_sequential_series(self):
        cells = analyze_run(self.ladder, series_length=8, prefixes=(1, 2, 8))
        null, signal = cells
        self.assertTrue(null["null_excitation"])
        self.assertNotIn("coverage_95", null["rules"]["symmetric"]["8"])
        self.assertIn("coverage_95", signal["rules"]["symmetric"]["8"])
        self.assertEqual(signal["series"], 8)
        # A single reading is the flat law under every rule.
        one = [signal["rules"][rule]["1"]["readout_ratio_quartiles"] for rule in ("symmetric", "mass_matching", "log_matching")]
        np.testing.assert_allclose(one[0], one[1], rtol=1e-12)
        np.testing.assert_allclose(one[0], one[2], rtol=1e-12)


if __name__ == "__main__":
    unittest.main()
