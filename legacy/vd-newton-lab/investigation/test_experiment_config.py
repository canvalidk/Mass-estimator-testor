"""Configuration validation and reproducibility through the shared engine."""

from contextlib import redirect_stdout
from copy import deepcopy
import io
import json
import math
from pathlib import Path
import unittest
from uuid import uuid4

import numpy as np

from comparison_estimators import baseline_points
from compare_estimators import scenarios
from experiment_config import (MAX_DIAGONAL_SAMPLES, MAX_DIAGONAL_TOTAL_TRIALS,
                               normalize_config, rms_sd, run_config, scenario_rows)


class ConfigurationTests(unittest.TestCase):
    def test_defaults_and_json_round_trip_do_not_modify_input(self):
        original = {"name": "  Boundary check  ", "scenario": {"type": "pairs", "pairs": [[8, .18], [8, 2.24]]}}
        before = deepcopy(original)
        config = normalize_config(original)
        self.assertEqual(original, before)
        self.assertEqual(config["name"], "Boundary check")
        self.assertEqual(config["samples"], 1024)
        self.assertEqual(config["order"], 96)
        self.assertEqual(config["noise_model"], "isotropic_unit_3d")
        self.assertEqual(normalize_config(json.loads(json.dumps(config, allow_nan=False))), config)
        config["scenario"]["pairs"][0][0] = 99
        self.assertEqual(original, before)
        self.assertEqual(len(scenario_rows({})), 16)

    def test_presets_match_existing_scenarios_exactly(self):
        for selection in ("all", "original", "fixed_mass"):
            self.assertEqual(scenario_rows({"scenario": {"type": "preset", "selection": selection}}),
                             scenarios(selection))

    def test_custom_paths_preserve_physical_identities_and_order(self):
        force_rows = scenario_rows({"scenario": {"type": "force_path", "force_snr": 8,
                                                 "acceleration_snrs": [.18, 1, 2.24]}})
        self.assertEqual([row["true_acceleration_snr"] for row in force_rows], [.18, 1, 2.24])
        self.assertEqual([row["true_force_snr"] for row in force_rows], [8, 8, 8])
        fixed_rows = scenario_rows({"scenario": {"type": "fixed_mass", "masses": [.25, 4],
                                                 "signal_snrs": [1, 8]}})
        np.testing.assert_allclose([row["true_mass"] for row in fixed_rows], [.25, .25, 4, 4])
        np.testing.assert_allclose([row["total_signal_snr"] for row in fixed_rows], [1, 8, 1, 8])
        pair_rows = scenario_rows({"scenario": {"type": "pairs", "pairs": [[16, 8], [4, .25]]}})
        self.assertEqual([row["true_mass"] for row in pair_rows], [2, 16])
        for rows in (force_rows, fixed_rows, pair_rows):
            self.assertEqual(len({row["id"] for row in rows}), len(rows))
            for row in rows:
                self.assertAlmostEqual(row["true_mass"] * row["true_acceleration_snr"], row["true_force_snr"])
                self.assertAlmostEqual(row["total_signal_snr"], math.hypot(row["true_force_snr"], row["true_acceleration_snr"]))

    def test_invalid_and_unknown_settings_are_rejected(self):
        invalid = [
            [], {"schema_version": 2}, {"schema_version": True}, {"name": ""}, {"name": " "},
            {"name": "x" * 101}, {"name": 1}, {"samples": 1}, {"samples": 1_048_577},
            {"samples": 1024.0}, {"samples": True}, {"order": 47}, {"order": 385},
            {"batch_size": 0}, {"batch_size": 4097}, {"training_seed": -1},
            {"heldout_seed": 2**32}, {"training_seed": 5, "heldout_seed": 5},
            {"save_trials": 1}, {"noise_model": "diagonal"}, {"noise_model": None},
            {"output": "anywhere"}, {"scenario": None}, {"scenario": {"type": "missing"}},
            {"scenario": {"type": "preset", "selection": "other"}},
            {"scenario": {"type": "preset", "ignored": 1}},
            {"scenario": {"type": "force_path", "force_snr": 8, "acceleration_snrs": []}},
            {"scenario": {"type": "fixed_mass", "masses": [1] * 9, "signal_snrs": [1] * 8}},
            {"scenario": {"type": "pairs", "pairs": [[1, 2, 3]]}},
            {"scenario": {"type": "pairs", "pairs": [[1, 2]] * 65}},
            {"samples": 1_048_576, "scenario": {"type": "pairs", "pairs": [[1, 2]] * 64}},
            {"scenario": {"type": "pairs", "pairs": [[1e308, 1e-308]]}},
        ]
        for value in (0, -1, math.inf, math.nan, "1", True, 10**1000):
            invalid.append({"scenario": {"type": "pairs", "pairs": [[1, value]]}})
        for value in invalid:
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    normalize_config(value)

    def test_mass_excitation_keeps_declared_mass_at_zero_and_round_trips(self):
        original = {"scenario": {"type": "mass_excitation", "masses": [.25, 4],
                                 "acceleration_snrs": [0, .5, 2]}}
        before = deepcopy(original)
        config = normalize_config(original)
        self.assertEqual(original, before)
        self.assertEqual(normalize_config(json.loads(json.dumps(config, allow_nan=False))), config)
        rows = scenario_rows(config)
        self.assertEqual([row["true_mass"] for row in rows], [.25, .25, .25, 4, 4, 4])
        self.assertEqual([row["true_acceleration_snr"] for row in rows], [0, .5, 2, 0, .5, 2])
        self.assertEqual(len({row["id"] for row in rows}), 6)
        for row in rows:
            self.assertEqual(row["true_acceleration_magnitude"], row["true_acceleration_snr"])
            self.assertEqual(row["true_force_magnitude"], row["true_mass"] * row["true_acceleration_magnitude"])
            self.assertEqual(row["true_force_snr"], row["true_force_magnitude"])
            self.assertEqual(row["total_signal_snr"], math.hypot(row["true_force_snr"], row["true_acceleration_snr"]))
        for row in (rows[0], rows[3]):
            self.assertEqual(row["true_force_magnitude"], 0.)
            self.assertEqual(row["true_acceleration_magnitude"], 0.)
            self.assertEqual(row["total_signal_snr"], 0.)
        config["scenario"]["acceleration_snrs"][0] = 10
        self.assertEqual(original, before)

    def test_mass_excitation_uses_true_acceleration_noise_and_physical_mass(self):
        config = {"noise_model": "diagonal_gaussian_3d", "direction": [0, 3, 4],
                  "noise": {"force_sd": [2, 4, 8], "acceleration_sd": [.5, 1, 2]},
                  "scenario": {"type": "mass_excitation", "masses": [.25, 4], "acceleration_snrs": [0, 2]}}
        rows = scenario_rows(config)
        self.assertEqual([row["true_mass"] for row in rows], [.25, .25, 4, 4])
        for row in rows:
            self.assertAlmostEqual(row["true_acceleration_magnitude"], row["true_acceleration_snr"] * rms_sd([.5, 1, 2]))
            self.assertAlmostEqual(row["true_force_magnitude"], row["true_mass"] * row["true_acceleration_magnitude"])
            self.assertAlmostEqual(row["true_force_snr"], row["true_force_magnitude"] / rms_sd([2, 4, 8]))
            np.testing.assert_allclose(row["direction"], [0, .6, .8])
        calibrated = {**config, "measurement": {"repeats": 10, "calibration": {
            "mode": "scaled", "force_scale": [.5] * 3, "acceleration_scale": [2] * 3}}}
        self.assertEqual(scenario_rows(calibrated), rows)

    def test_mass_excitation_changes_units_without_changing_signal_to_noise(self):
        original = {"noise_model": "diagonal_gaussian_3d",
                    "noise": {"force_sd": [.5, 1, 2], "acceleration_sd": [2, 1, .5]},
                    "scenario": {"type": "mass_excitation", "masses": [.25, 4], "acceleration_snrs": [0, .5, 2]}}
        changed = deepcopy(original)
        force_scale, acceleration_scale = 7., .4
        changed["noise"]["force_sd"] = [value * force_scale for value in original["noise"]["force_sd"]]
        changed["noise"]["acceleration_sd"] = [value * acceleration_scale for value in original["noise"]["acceleration_sd"]]
        changed["scenario"]["masses"] = [mass * force_scale / acceleration_scale for mass in original["scenario"]["masses"]]
        for before, after in zip(scenario_rows(original), scenario_rows(changed)):
            for key in ("true_force_snr", "true_acceleration_snr", "total_signal_snr"):
                self.assertAlmostEqual(after[key], before[key])
            self.assertAlmostEqual(after["true_force_magnitude"], before["true_force_magnitude"] * force_scale)
            self.assertAlmostEqual(after["true_acceleration_magnitude"], before["true_acceleration_magnitude"] * acceleration_scale)
            self.assertEqual(after["true_mass"], before["true_mass"] * force_scale / acceleration_scale)

    def test_mass_excitation_rejects_invalid_domain_and_excessive_product(self):
        scenario = {"type": "mass_excitation", "masses": [1], "acceleration_snrs": [0]}
        invalid = [{"masses": [0]}, {"masses": [-1]}, {"masses": [math.inf]},
                   {"masses": []}, {"acceleration_snrs": []}, {"acceleration_snrs": [0] * 65},
                   {"masses": [1] * 9, "acceleration_snrs": [0] * 8}, {"signal_snrs": [1]}]
        invalid += [{"acceleration_snrs": [value]} for value in (-1, math.inf, math.nan, "0", False, 10**1000)]
        for changes in invalid:
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                normalize_config({"scenario": {**scenario, **changes}})
        # Zero is still rejected by older signal designs.
        for old in ({"type": "pairs", "pairs": [[0, 0]]},
                    {"type": "force_path", "force_snr": 1, "acceleration_snrs": [0]},
                    {"type": "fixed_mass", "masses": [1], "signal_snrs": [0]}):
            with self.subTest(old=old), self.assertRaises(ValueError):
                normalize_config({"scenario": old})

    def test_mass_excitation_rejects_overflow_and_positive_underflow(self):
        cases = [(1e308, 2., None), (1e-300, 1e-300, None), (1., 1.7e308, None),
                 (1., 1e-300, 1e-300), (1e-300, 1e-300, 1e300), (1., 1e10, 1e300)]
        for mass, acceleration_snr, sd in cases:
            config = {"scenario": {"type": "mass_excitation", "masses": [mass],
                                   "acceleration_snrs": [acceleration_snr]}}
            if sd is not None:
                config.update(noise_model="diagonal_gaussian_3d",
                              noise={"force_sd": [sd] * 3, "acceleration_sd": [sd] * 3})
            with self.subTest(config=config), self.assertRaisesRegex(ValueError, "finite positive"):
                normalize_config(config)
            # The same supported scales/mass with an explicit zero are valid.
            config["scenario"]["acceleration_snrs"] = [0]
            row = scenario_rows(config)[0]
            self.assertEqual(row["true_mass"], mass)
            self.assertEqual(row["true_force_magnitude"], 0.)
            self.assertEqual(row["true_acceleration_magnitude"], 0.)

    def test_invalid_configuration_creates_no_output(self):
        output = Path(__file__).resolve().parents[3] / ".tools" / "mass_comparison_tests" / ("rejected-" + uuid4().hex)
        with self.assertRaisesRegex(ValueError, "noise_model"):
            run_config({"noise_model": "anisotropic"}, output)
        self.assertFalse(output.exists())

    def test_diagonal_defaults_and_direction_are_reproducible(self):
        original = {"noise_model": "diagonal_gaussian_3d", "direction": [3, -4, 2],
                    "noise": {"force_sd": [.5, 1, 2], "acceleration_sd": [2, 1, .5]}}
        before = deepcopy(original)
        config = normalize_config(original)
        self.assertEqual(original, before)
        self.assertAlmostEqual(math.hypot(*config["direction"]), 1.)
        np.testing.assert_allclose(config["direction"], np.array([3, -4, 2]) / math.sqrt(29))
        for _ in range(10):
            self.assertEqual(normalize_config(json.loads(json.dumps(config))), config)
        defaults = normalize_config({"noise_model": "diagonal_gaussian_3d"})
        self.assertEqual(defaults["noise"], {"force_sd": [1, 1, 1], "acceleration_sd": [1, 1, 1]})
        self.assertEqual(defaults["direction"], [1, 0, 0])
        self.assertNotIn("noise", normalize_config({}))
        self.assertNotIn("direction", normalize_config({}))
        config["noise"]["force_sd"][0] = 8
        self.assertEqual(original, before)

    def test_diagonal_pairs_use_rms_snr_and_physical_mass(self):
        config = {"noise_model": "diagonal_gaussian_3d",
                  "noise": {"force_sd": [2, 4, 8], "acceleration_sd": [.5, 1, 2]},
                  "direction": [0, 3, 4],
                  "scenario": {"type": "pairs", "pairs": [[8, 2], [4, 8]]}}
        rows = scenario_rows(config)
        np.testing.assert_allclose([row["true_mass"] for row in rows], [16, 2])
        for row in rows:
            self.assertAlmostEqual(row["true_force_magnitude"] / rms_sd([2, 4, 8]), row["true_force_snr"])
            self.assertAlmostEqual(row["true_acceleration_magnitude"] / rms_sd([.5, 1, 2]), row["true_acceleration_snr"])
            self.assertAlmostEqual(row["true_mass"] * row["true_acceleration_magnitude"], row["true_force_magnitude"])
            np.testing.assert_allclose(row["direction"], [0, .6, .8])
        self.assertEqual([row["true_force_snr"] for row in rows], [8, 4])
        self.assertEqual([row["true_acceleration_snr"] for row in rows], [2, 8])

    def test_diagonal_fixed_mass_and_preset_keep_physical_masses(self):
        common = {"noise_model": "diagonal_gaussian_3d",
                  "noise": {"force_sd": [2, 4, 8], "acceleration_sd": [.5, 1, 2]}}
        rows = scenario_rows({**common, "scenario": {"type": "fixed_mass", "masses": [.25, 4], "signal_snrs": [1, 8]}})
        np.testing.assert_allclose([row["true_mass"] for row in rows], [.25, .25, 4, 4])
        np.testing.assert_allclose([row["total_signal_snr"] for row in rows], [1, 8, 1, 8])
        for row in rows:
            self.assertAlmostEqual(row["true_force_snr"] / row["true_acceleration_snr"], row["true_mass"] / 4)
        for selection in ("all", "fixed_mass", "original"):
            actual = scenario_rows({**common, "scenario": {"type": "preset", "selection": selection}})
            expected = scenarios(selection)
            for row, old in zip(actual, expected):
                scale = 1 if old["family"] == "fixed_mass" else 4
                self.assertAlmostEqual(row["true_mass"], scale * old["true_mass"])
                self.assertAlmostEqual(math.hypot(row["true_force_snr"], row["true_acceleration_snr"]), old["total_signal_snr"])

    def test_diagonal_invalid_noise_and_truth_direction_are_rejected(self):
        invalid = [
            {"noise": None}, {"noise": {"ignored": [1, 1, 1]}},
            {"noise": {"force_sd": [1, 2]}}, {"noise": {"force_sd": [1, 1, 17]}},
            {"noise": {"acceleration_sd": [1, 1, 0]}},
            {"direction": [0, 0, 0]}, {"direction": [1, 0]},
            {"direction": [1, "0", 0]}, {"direction": [True, 0, 0]},
            {"direction": [math.inf, 0, 0]}, {"direction": [10**1000, 0, 0]},
            {"samples": MAX_DIAGONAL_SAMPLES + 1},
            {"samples": 8192, "scenario": {"type": "pairs", "pairs": [[1, 1]] * 33}},
            {"noise": {"force_sd": [1e300] * 3, "acceleration_sd": [1e-300] * 3}},
        ]
        for value in (0, -1, math.inf, math.nan, "1", True):
            invalid.append({"noise": {"force_sd": [1, 1, value]}})
        for settings in invalid:
            with self.subTest(settings=settings):
                with self.assertRaises(ValueError):
                    normalize_config({"noise_model": "diagonal_gaussian_3d", **settings})
        for settings in ({"noise": {}}, {"direction": [1, 0, 0]}):
            with self.assertRaisesRegex(ValueError, "require"):
                normalize_config(settings)

    def test_diagonal_normalization_handles_large_and_small_finite_scales(self):
        for scale in (1e300, 1e-300):
            self.assertAlmostEqual(rms_sd([scale] * 3) / scale, 1.)
            config = normalize_config({"noise_model": "diagonal_gaussian_3d", "direction": [scale, -scale, scale],
                                       "noise": {"force_sd": [scale] * 3, "acceleration_sd": [scale] * 3},
                                       "scenario": {"type": "pairs", "pairs": [[1, 1]]}})
            np.testing.assert_allclose(config["direction"], np.array([1, -1, 1]) / math.sqrt(3))
            self.assertAlmostEqual(scenario_rows(config)[0]["true_mass"], 1.)

    def test_measurement_contract_is_optional_canonical_and_preserves_truth(self):
        self.assertNotIn("measurement", normalize_config({}))
        original = {"measurement": {"repeats": 5},
                    "scenario": {"type": "pairs", "pairs": [[8, .25]]}}
        before = deepcopy(original)
        config = normalize_config(original)
        self.assertEqual(original, before)
        self.assertEqual(config["measurement"], {"repeats": 5, "calibration": {"mode": "known"}})
        self.assertEqual(normalize_config(json.loads(json.dumps(config))), config)
        self.assertEqual(scenario_rows(config), scenario_rows({"scenario": original["scenario"]}))
        for calibration in ({"mode": "scaled", "force_scale": [.9, 1, 1.1], "acceleration_scale": [1.5]*3},
                            {"mode": "estimated", "samples": 30, "seed": 42}):
            configured = normalize_config({**original, "measurement": {"repeats": 50, "calibration": calibration}})
            self.assertEqual(normalize_config(json.loads(json.dumps(configured))), configured)

    def test_measurement_rejects_unsupported_settings_before_output(self):
        invalid = [None, [], {"unknown": 1}, {"repeats": 0}, {"repeats": 1001}, {"repeats": True},
                   {"calibration": None}, {"calibration": {"mode": "missing"}},
                   {"calibration": {"mode": "known", "samples": 30}},
                   {"calibration": {"mode": "scaled", "force_scale": [1, 2], "acceleration_scale": [1]*3}},
                   {"calibration": {"mode": "scaled", "force_scale": [1, 1, 0], "acceleration_scale": [1]*3}},
                   {"calibration": {"mode": "scaled", "force_scale": [1, 1, 17], "acceleration_scale": [1]*3}},
                   {"calibration": {"mode": "estimated", "samples": 3, "seed": 42}},
                   {"calibration": {"mode": "estimated", "samples": 100001, "seed": 42}},
                   {"calibration": {"mode": "estimated", "samples": 30, "seed": -1}},
                   {"calibration": {"mode": "estimated", "samples": 30, "seed": 2026091601}},
                   {"calibration": {"mode": "estimated", "samples": 30, "seed": 2026091702}}]
        for value in invalid:
            with self.subTest(value=value), self.assertRaises(ValueError):
                normalize_config({"measurement": value})
        with self.assertRaisesRegex(ValueError, str(MAX_DIAGONAL_SAMPLES)):
            normalize_config({"samples": MAX_DIAGONAL_SAMPLES + 1, "measurement": {"repeats": 2}})
        # The bound concerns supplied SDs, not calibration factors alone.
        config = {"noise_model": "diagonal_gaussian_3d",
                  "noise": {"force_sd": [1, 4, 16]},
                  "measurement": {"calibration": {"mode": "scaled", "force_scale": [16, 4, 1], "acceleration_scale": [1]*3}}}
        self.assertEqual(normalize_config(config)["measurement"]["calibration"]["force_scale"], [16, 4, 1])
        config["measurement"]["calibration"]["force_scale"] = [1, 1, 2]
        with self.assertRaisesRegex(ValueError, "Supplied mean force"):
            normalize_config(config)

    def test_pooled_isotropic_calibration_contract_and_rejections(self):
        calibration = {"mode": "pooled_isotropic", "samples": 4, "seed": 42}
        settings = {"measurement": {"repeats": 4, "calibration": calibration},
                    "noise_model": "diagonal_gaussian_3d",
                    "noise": {"force_sd": [2]*3, "acceleration_sd": [.5]*3}}
        config = normalize_config(settings)
        self.assertEqual(config["measurement"]["calibration"], {**calibration, "calibration_order": 24})
        self.assertEqual(normalize_config(json.loads(json.dumps(config))), config)
        for changes in ({"samples": 3}, {"samples": True}, {"seed": 2026091601},
                        {"calibration_order": 7}, {"calibration_order": 97},
                        {"calibration_order": True}, {"calibration_order": 24.5},
                        {"force_scale": [1]*3}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                normalize_config({"measurement": {"calibration": {**calibration, **changes}}})
        for channel in ("force_sd", "acceleration_sd"):
            with self.subTest(channel=channel), self.assertRaisesRegex(ValueError, "requires isotropic true"):
                normalize_config({**settings, "noise": {**settings["noise"], channel: [1, 1, 2]}})
        with self.assertRaises(ValueError):
            normalize_config({"measurement": {"calibration": {"mode": "estimated", "samples": 30,
                                                               "seed": 42, "calibration_order": 24}}})

    def test_ten_thousand_trial_design_keeps_measurement_work_and_batch_guards(self):
        # The planned 12-cell study fits the old total-work budget even though
        # each cell now has enough replications for its prespecified precision.
        design = {"samples": 10_000, "batch_size": 128,
                  "scenario": {"type": "mass_excitation", "masses": [.25, 1, 4],
                               "acceleration_snrs": [0, .5, 2, 8]}}
        modes = [
            {"noise_model": "diagonal_gaussian_3d", "noise": {
                "force_sd": [.5, 1, 2], "acceleration_sd": [2, 1, .5]}},
            {"measurement": {"repeats": 4}},
            {"measurement": {"calibration": {"mode": "scaled", "force_scale": [.75] * 3,
                                              "acceleration_scale": [1.25] * 3}}},
            {"measurement": {"calibration": {"mode": "estimated", "samples": 30, "seed": 42}}},
        ]
        for mode in modes:
            with self.subTest(mode=mode):
                config = normalize_config({**design, **mode})
                self.assertEqual(config["samples"], 10_000)
                self.assertEqual(config["batch_size"], 128)
                self.assertEqual(len(scenario_rows(config)), 12)
                self.assertLessEqual(config["samples"] * 12, MAX_DIAGONAL_TOTAL_TRIALS)
                self.assertEqual(normalize_config(json.loads(json.dumps(config))), config)
                with self.assertRaisesRegex(ValueError, str(MAX_DIAGONAL_SAMPLES)):
                    normalize_config({**design, **mode, "samples": MAX_DIAGONAL_SAMPLES + 1})
        self.assertEqual(MAX_DIAGONAL_TOTAL_TRIALS, 262_144)

    def test_extended_capacity_rejects_total_work_before_creating_output(self):
        count = MAX_DIAGONAL_TOTAL_TRIALS // 10_000
        for mode in ({"noise_model": "diagonal_gaussian_3d"}, {"measurement": {"repeats": 2}}):
            config = {"samples": 10_000, "scenario": {"type": "mass_excitation", "masses": [1.],
                                                    "acceleration_snrs": list(range(count))}, **mode}
            self.assertEqual(len(scenario_rows(config)), count)
            config["scenario"]["acceleration_snrs"].append(count)
            output = Path(__file__).resolve().parents[3] / ".tools" / "mass_comparison_tests" / ("capacity-rejected-" + uuid4().hex)
            with self.subTest(mode=mode), self.assertRaisesRegex(ValueError, "262,144 total trials"):
                run_config(config, output)
            self.assertFalse(output.exists())


class ConfiguredRunTests(unittest.TestCase):
    def test_saved_configuration_reproduces_noise_and_paired_results(self):
        scratch = Path(__file__).resolve().parents[3] / ".tools" / "mass_comparison_tests" / ("configured-" + uuid4().hex)
        first, second = scratch / "first", scratch / "second"
        config = normalize_config({"name": "Seed replay", "samples": 8, "batch_size": 5,
                                   "training_seed": 101, "heldout_seed": 202, "save_trials": True,
                                   "scenario": {"type": "pairs", "pairs": [[8, .25], [4, 2]]}})
        events = []
        with redirect_stdout(io.StringIO()) as stdout:
            report = run_config(config, first, progress=events.append)
            replay = run_config(json.loads((first / "config.json").read_text(encoding="utf-8")), second)
        self.assertEqual(report["experiment_config"], config)
        self.assertEqual(report["training_seed"], 101)
        self.assertEqual(report["heldout_seed"], 202)
        self.assertIn("experiment_config.py", report["source_sha256"])
        self.assertEqual(report["rows"], replay["rows"])
        self.assertEqual([event["completed_count"] for event in events], [1, 2])
        self.assertTrue(all(event["total"] == 2 for event in events))
        self.assertEqual([event["completed"] for event in events], [row["id"] for row in report["rows"]])
        self.assertEqual(json.loads(stdout.getvalue().splitlines()[0]), events[0])
        self.assertEqual(json.loads((first / "results.json").read_text(encoding="utf-8")), report)

        noise = np.random.default_rng(101).standard_normal((8, 6))
        for row in report["rows"]:
            expected = baseline_points(noise[:, :3] + [row["true_force_snr"], 0, 0],
                                       noise[:, 3:] + [row["true_acceleration_snr"], 0, 0])
            with np.load(first / "trials" / (row["id"] + ".npz"), allow_pickle=False) as arrays:
                for method, points in expected.items():
                    np.testing.assert_array_equal(arrays["point__" + method], points)

        # Altering only the independent validation seed changes the prediction
        # task, while leaving fitted points and mass-accuracy scores unchanged.
        altered = {**config, "heldout_seed": 303, "save_trials": False}
        with redirect_stdout(io.StringIO()):
            changed = run_config(altered, scratch / "heldout_changed")
        old = report["rows"][0]["methods"]["flat_joint"]
        new = changed["rows"][0]["methods"]["flat_joint"]
        self.assertEqual(old["means"], new["means"])
        self.assertNotEqual(old["additional_tasks"]["heldout_known_acceleration_relative_squared_force_error"],
                            new["additional_tasks"]["heldout_known_acceleration_relative_squared_force_error"])


if __name__ == "__main__":
    unittest.main()
