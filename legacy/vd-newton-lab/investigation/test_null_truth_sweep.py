"""Truth changes require an unchanged null data law, outputs and verified bytes."""

import hashlib
import json
import math
from pathlib import Path
import unittest
from unittest.mock import patch
from uuid import uuid4

import numpy as np

from experiment_config import normalize_config, scenario_rows
import null_truth_sweep as sweep
from trial_archive import write_trial_archive


def arrays(extra_laws=()):
    point = np.array([.5, 1., 2., 4.])
    result = {"observation__force": np.zeros((4, 3)),
              "observation__acceleration": np.ones((4, 3)), "point__flat_joint": point}
    for law in sweep.LAWS + tuple(extra_laws):
        for level, factor in ((50, 1.25), (95, 2.)):
            result[f"posterior__{law}__log_lower_{level}"] = np.log(point / factor)
            result[f"posterior__{law}__log_upper_{level}"] = np.log(point * factor)
    return result


def write_fixture(root, alter_second=False, null=True, calibration=None, trial_factory=arrays):
    raw_config = {"samples": 4, "save_trials": True,
                  "scenario": {"type": "mass_excitation", "masses": [.25, 4],
                               "acceleration_snrs": [0 if null else 1]}}
    if calibration is not None:
        raw_config["measurement"] = {"repeats": 4, "calibration": calibration}
    config = normalize_config(raw_config)
    spec = {"experiments": [{"id": "n1", "config": config}]}
    relative = "experiments/n1/attempt-001"
    folder = root / relative
    folder.mkdir(parents=True)
    rows = scenario_rows(config)
    entries = []
    for i, row in enumerate(rows):
        trial = trial_factory()
        if alter_second and i == 1:
            trial["point__flat_joint"][0] += .1
        entries.append(write_trial_archive(folder, row, trial))
    manifest = {"schema_version": 1, "entries": entries}
    (folder / "trial_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    report = {"experiment_config": config, "samples_per_scenario": 4, "rows": rows}
    (folder / "results.json").write_text(json.dumps(report), encoding="utf-8")
    hashes = {p.relative_to(folder).as_posix(): sweep._hash(p) for p in folder.rglob("*") if p.is_file()}
    spec_hash = hashlib.sha256(json.dumps(spec, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
    state = {"schema_version": 1, "status": "completed", "specification": spec, "spec_sha256": spec_hash,
             "source_fingerprint": {"historical_sources": "deliberately not today's code"},
             "experiments": [{"id": "n1", "status": "completed", "attempt_directory": relative,
                              "output_sha256": hashes}]}
    (root / "study_checkpoint.json").write_text(json.dumps(state), encoding="utf-8")
    return state, folder


class ScoreTests(unittest.TestCase):
    def test_fixed_arrays_give_truth_specific_coverage_and_all_trial_losses(self):
        original = arrays()
        before = {k: v.copy() for k, v in original.items()}
        rows = sweep.score_null_trials(original, [1., 64.])
        first = rows[0]
        interval = first["intervals"]["flat_joint"]["95"]
        self.assertEqual(interval["coverage"]["mean"], .75)
        self.assertEqual(interval["truth_below_interval"]["mean"], .25)
        self.assertEqual(interval["truth_above_interval"]["mean"], 0.)
        self.assertAlmostEqual(interval["coverage"]["mcse"], .25)
        self.assertAlmostEqual(first["flat_point"]["capped_squared_log_factor_2"]["mean"], .75 * math.log(2) ** 2)
        self.assertEqual(first["flat_point"]["outside_factor_2"]["mean"], .25)
        self.assertEqual(rows[1]["intervals"]["flat_joint"]["95"]["coverage"]["mean"], 0.)
        for key in before:
            np.testing.assert_array_equal(original[key], before[key])

    def test_factor_two_endpoints_are_inclusive_even_at_large_powers_of_two(self):
        data = arrays()
        data["point__flat_joint"] = np.array([32., 64., 128., 129.])
        result = sweep.score_null_trials(data, [64.])[0]
        self.assertEqual(result["flat_point"]["outside_factor_2"]["mean"], .25)

    def test_invalid_outputs_remain_in_denominators_and_width_is_unavailable(self):
        data = arrays()
        data["point__flat_joint"][0] = np.nan
        data["posterior__flat_joint__log_lower_95"][0] = np.nan
        result = sweep.score_null_trials(data, [1.])[0]
        self.assertEqual(result["flat_point"]["invalid_estimate"]["mean"], .25)
        self.assertEqual(result["flat_point"]["outside_factor_2"]["mean"], .5)
        interval = result["intervals"]["flat_joint"]["95"]
        self.assertEqual(interval["invalid_interval"]["mean"], .25)
        self.assertEqual(sum(interval[k]["mean"] for k in ("coverage", "truth_below_interval", "truth_above_interval", "invalid_interval")), 1.)
        self.assertIsNone(interval["log_width"]["mean"])

    def test_rejects_invalid_truth_grid_and_missing_intervals(self):
        for masses in ([], [0], [True], [float("inf")], [1., 1.], [2., 1.]):
            with self.subTest(masses=masses), self.assertRaises(ValueError):
                sweep.score_null_trials(arrays(), masses)
        data = arrays()
        del data["posterior__tube_joint__log_lower_50"]
        with self.assertRaisesRegex(ValueError, "requires archived interval"):
            sweep.score_null_trials(data, [1.])

    def test_discovers_every_archived_law_and_retains_invalid_trials(self):
        data = arrays(("calibrated_joint", "oracle_joint", "future_law"))
        data["posterior__calibrated_joint__log_lower_95"][0] = np.nan
        rows = sweep.score_null_trials(data, [1., 64.])
        self.assertEqual(set(rows[0]["intervals"]),
                         set(sweep.LAWS) | {"calibrated_joint", "oracle_joint", "future_law"})
        interval = rows[0]["intervals"]["calibrated_joint"]["95"]
        self.assertEqual(interval["coverage"]["mean"], .5)
        self.assertEqual(interval["invalid_interval"]["mean"], .25)
        self.assertEqual(sum(interval[key]["mean"] for key in
                             ("coverage", "truth_below_interval", "truth_above_interval", "invalid_interval")), 1.)
        self.assertIsNone(interval["log_width"]["mean"])
        for law in ("oracle_joint", "future_law"):
            self.assertEqual(rows[0]["intervals"][law]["95"]["coverage"]["mean"], .75)
            self.assertEqual(rows[1]["intervals"][law]["95"]["coverage"]["mean"], 0.)
        self.assertEqual(rows[0]["flat_point"], sweep.score_null_trials(arrays(), [1.])[0]["flat_point"])

    def test_partial_extra_law_cannot_be_silently_omitted(self):
        for missing in ("log_lower_50", "log_upper_95"):
            data = arrays(("calibrated_joint",))
            del data[f"posterior__calibrated_joint__{missing}"]
            with self.subTest(missing=missing), self.assertRaisesRegex(ValueError, "calibrated_joint"):
                sweep.score_null_trials(data, [1.])
        data = arrays()
        data["posterior__oracle_joint__cdf_truth"] = np.ones(4)
        with self.assertRaisesRegex(ValueError, "oracle_joint"):
            sweep.score_null_trials(data, [1.])


class CompletedStudyTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(sweep.__file__).resolve().parents[3] / ".tools" / "null_truth_sweep_tests" / uuid4().hex
        self.study = self.root / "study"

    def test_valid_completed_study_checks_invariance_without_refitting(self):
        state, folder = write_fixture(self.study)
        before = {p: p.read_bytes() for p in self.study.rglob("*") if p.is_file()}
        with patch("comparison_posterior.infer_batch", side_effect=AssertionError("Must not fit")):
            report = sweep.build_null_sweep(self.study, [1., 64.])
        self.assertEqual(report["status"], "completed")
        self.assertEqual(report["diagnostic_runtime"]["numpy"], np.__version__)
        self.assertTrue(report["diagnostic_runtime"]["python"])
        self.assertNotIn("matplotlib", report["diagnostic_runtime"])
        experiment = report["experiments"][0]
        self.assertEqual(experiment["invariance"]["status"], "verified")
        self.assertEqual(len(experiment["source_null_scenarios"]), 2)
        self.assertEqual(experiment["rows"][0]["samples"], 4)
        self.assertEqual({p: p.read_bytes() for p in before}, before)

    def test_pooled_null_discovers_saved_laws_without_refitting(self):
        write_fixture(self.study, calibration={"mode": "pooled_isotropic", "samples": 4, "seed": 9827},
                      trial_factory=lambda: arrays(("calibrated_joint", "oracle_joint")))
        before = {path: path.read_bytes() for path in self.study.rglob("*") if path.is_file()}
        with patch("comparison_posterior.infer_batch", side_effect=AssertionError("Must not fit")), \
                patch("calibration_posterior.infer_batch", side_effect=AssertionError("Must not fit")):
            result = sweep.build_null_sweep(self.study, [1., 64.])
        experiment = result["experiments"][0]
        self.assertEqual(experiment["measurement"]["calibration"]["mode"], "pooled_isotropic")
        self.assertEqual(experiment["invariance"]["status"], "verified")
        for law in sweep.LAWS + ("calibrated_joint", "oracle_joint"):
            interval = experiment["rows"][0]["intervals"][law]["95"]
            self.assertEqual(interval["coverage"]["mean"], .75)
            self.assertAlmostEqual(interval["coverage"]["mcse"], .25)
        self.assertEqual({path: path.read_bytes() for path in before}, before)

    def test_completed_archive_with_incomplete_extra_law_is_rejected(self):
        def incomplete():
            data = arrays(("oracle_joint",))
            del data["posterior__oracle_joint__log_upper_95"]
            return data
        write_fixture(self.study, trial_factory=incomplete)
        with self.assertRaisesRegex(ValueError, "requires archived interval array posterior__oracle_joint"):
            sweep.build_null_sweep(self.study)

    def test_partial_law_only_in_later_null_archive_is_not_ignored(self):
        calls = 0
        def later_partial():
            nonlocal calls
            calls += 1
            data = arrays()
            if calls == 2:
                data["posterior__oracle_joint__cdf_truth"] = np.ones(4)
            return data
        write_fixture(self.study, trial_factory=later_partial)
        with self.assertRaisesRegex(ValueError, "requires archived interval array posterior__oracle_joint"):
            sweep.build_null_sweep(self.study)

    def test_changed_output_bytes_are_rejected_before_scoring(self):
        state, folder = write_fixture(self.study)
        archive = next((folder / "trials").glob("*.npz"))
        archive.write_bytes(archive.read_bytes() + b"changed")
        with self.assertRaisesRegex(ValueError, "output hash mismatch"):
            sweep.build_null_sweep(self.study)

    def test_truth_dependent_null_fit_is_rejected_even_when_hashes_match(self):
        write_fixture(self.study, alter_second=True)
        with self.assertRaisesRegex(ValueError, "vary across stipulated truth"):
            sweep.build_null_sweep(self.study)

    def test_non_null_or_incomplete_study_is_rejected(self):
        state, folder = write_fixture(self.study, null=False)
        with self.assertRaisesRegex(ValueError, "true-zero-excitation"):
            sweep.build_null_sweep(self.study)
        state["status"] = "failed"
        (self.study / "study_checkpoint.json").write_text(json.dumps(state), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "completed version"):
            sweep.build_null_sweep(self.study)

    def test_path_escape_and_existing_output_cannot_overwrite_source(self):
        state, folder = write_fixture(self.study)
        with self.assertRaisesRegex(ValueError, "must not overwrite"):
            sweep.write_null_sweep(self.study, folder / "diagnostic")
        with self.assertRaisesRegex(ValueError, "already exists"):
            sweep.write_null_sweep(self.study, self.root)
        state["experiments"][0]["attempt_directory"] = "../escape"
        (self.study / "study_checkpoint.json").write_text(json.dumps(state), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "inside the run directory"):
            sweep.build_null_sweep(self.study)


if __name__ == "__main__":
    unittest.main()
