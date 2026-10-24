"""Saved trials preserve pairing, failure penalties and independent rescoring."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import unittest
from uuid import uuid4

import numpy as np

from comparison_scores import paired_summary, point_losses, score_point
from trial_archive import (
    ARCHIVE_SCHEMA_VERSION, load_trial_archive, load_trial_manifest,
    rescore_run, write_trial_archive,
)


class TrialArchiveTests(unittest.TestCase):
    def setUp(self):
        self.output = (Path(__file__).resolve().parents[3] / ".tools" /
                       "mass_trial_archive_tests" / uuid4().hex)
        self.row = {"id": "00_fixture", "true_mass": 2.,
                    "true_force_magnitude": 0., "true_acceleration_magnitude": 0.,
                    "metric_applicability": {
                        "heldout_relative_force": {"status": "not_applicable", "reason": "zero_true_force"}}}
        self.arrays = {
            "observation__force": np.arange(18., dtype=float).reshape(6, 3),
            "observation__acceleration": np.ones((6, 3)),
            "point__flat_joint": np.array([2., 1., 4., 0., np.nan, np.inf]),
            "point__comparator": np.array([2., 2., 3., 8., 1., 4.]),
            "posterior__flat_joint__cdf_truth": np.linspace(.1, .9, 6),
            "noise__supplied_mean_force_sd": np.ones((6, 3)) * .5,
            "heldout__true_force": np.zeros((6, 3)),
        }

    def save(self):
        entry = write_trial_archive(self.output, self.row, self.arrays)
        manifest = {"schema_version": ARCHIVE_SCHEMA_VERSION, "entries": [entry],
                    "provenance": {"training_seed": 17, "source_sha256": {"runner": "fixture"}}}
        (self.output / "trial_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        return entry, manifest

    def test_round_trip_preserves_invalid_points_and_metadata(self):
        entry, manifest = self.save()
        self.assertNotIn("replicate_id", self.arrays)
        self.assertEqual(load_trial_manifest(self.output), manifest)
        arrays = load_trial_archive(self.output, entry)
        for name, expected in self.arrays.items():
            np.testing.assert_array_equal(arrays[name], expected)
        np.testing.assert_array_equal(arrays["replicate_id"], np.arange(6))
        self.assertEqual(entry["scenario"]["metric_applicability"], self.row["metric_applicability"])
        self.assertEqual(entry["point_methods"], ["comparator", "flat_joint"])
        with self.assertRaises(FileExistsError):
            write_trial_archive(self.output, self.row, self.arrays)

    def test_rescore_matches_point_engine_and_paired_custom_failures(self):
        entry, _ = self.save()
        before = {path: path.read_bytes() for path in self.output.rglob("*") if path.is_file()}
        result = rescore_run(self.output, self.output / "rescore.json", factor_tolerance=2., cap_factor=3.)
        row = result["rows"][0]
        self.assertEqual(row["archive_sha256"], entry["sha256"])
        for name in entry["point_methods"]:
            self.assertEqual(row["methods"][name]["original_point_scores"],
                             score_point(self.arrays["point__" + name], 2.))
        actual = row["paired_reference_minus_other"]["comparator"]
        for metric in point_losses(self.arrays["point__comparator"], 2.):
            # Truth below is deliberately non-default; point scores must use
            # scenario truth, not a unit-mass shortcut.
            expected = paired_summary(point_losses(self.arrays["point__flat_joint"], 2.)[metric],
                                      point_losses(self.arrays["point__comparator"], 2.)[metric])
            self.assertEqual(actual["original_point_losses"][metric], expected)
        outside_flat = np.array([0., 0., 0., 1., 1., 1.])
        outside_other = np.array([0., 0., 0., 1., 0., 0.])
        self.assertEqual(actual["configurable_losses"]["outside_factor"],
                         paired_summary(outside_flat, outside_other))
        cap = np.log(3.) ** 2
        custom_flat = np.array([0., np.log(2.)**2, np.log(2.)**2, cap, cap, cap])
        custom_other = np.minimum(point_losses(self.arrays["point__comparator"], 2.)["squared_log"], cap)
        expected = paired_summary(custom_flat, custom_other)
        self.assertAlmostEqual(actual["configurable_losses"]["capped_squared_log"]["mean_difference"],
                               expected["mean_difference"])
        self.assertEqual(row["methods"]["flat_joint"]["original_point_scores"]["invalid_count"], 3)
        self.assertIsNone(actual["original_point_losses"]["squared_log"]["mean_difference"])
        self.assertEqual(json.loads((self.output / "rescore.json").read_text()), result)
        for path, content in before.items():
            self.assertEqual(path.read_bytes(), content)

    def test_rejects_changed_archive_even_when_npz_is_readable(self):
        entry, _ = self.save()
        archive = self.output / entry["path"]
        changed = dict(self.arrays, replicate_id=np.arange(6))
        changed["point__comparator"] = changed["point__comparator"] + .1
        np.savez_compressed(archive, **changed)
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            load_trial_archive(self.output, entry)

    def test_known_risk_warnings_survive_finite_samples_and_reversed_pairing(self):
        self.arrays["point__flat_joint"] = np.array([2., 1., 4., 3., 2., 2.])
        note = "Known infinite population variance: empirical MCSE is descriptive only."
        self.row["methods"] = {"comparator": {
            "population_unbounded_risk": "infinite_under_declared_failure_policy",
            "population_risk_note": "Analytic population limitation persists in finite valid samples."}}
        self.row["paired_flat_minus_other"] = {"comparator": {
            "squared_mass": {"mcse_interpretation": note}}}
        entry, _ = self.save()
        self.assertEqual(entry["method_risk"]["comparator"]["source_paired_mcse_notes"]["squared_mass"], note)
        for reference in ("flat_joint", "comparator"):
            result = rescore_run(self.output, self.output / (reference + ".json"), reference_method=reference)
            row = result["rows"][0]
            summary = row["methods"]["comparator"]["original_point_scores"]
            self.assertEqual(summary["invalid_count"], 0)
            self.assertEqual(summary["population_unbounded_risk"], "infinite_under_declared_failure_policy")
            other = "comparator" if reference == "flat_joint" else "flat_joint"
            pair = row["paired_reference_minus_other"][other]["original_point_losses"]["squared_mass"]
            self.assertEqual(pair["mcse_interpretation"], note)
            self.assertIn("comparator", pair["source_method_risk"])
            self.assertIsNotNone(pair["mcse"])

    def test_writer_rejects_misalignment_objects_and_unsafe_ids_before_write(self):
        variants = []
        for key, value in (
            ("point__comparator", np.ones(5)),
            ("point__comparator", np.array([object()] * 6, dtype=object)),
            ("point__comparator", np.ones(6, dtype=complex)),
            ("replicate_id", np.array([0, 2, 1, 3, 4, 5])),
            ("observation__force", np.full((6, 3), np.nan)),
            ("noise__supplied_mean_force_sd", np.zeros((6, 3))),
        ):
            variants.append((self.row, dict(self.arrays, **{key: value})))
        variants.append((dict(self.row, id="../escape"), self.arrays))
        variants.append((dict(self.row, true_mass=0.), self.arrays))
        for row, arrays in variants:
            with self.subTest(row=row, keys=list(arrays)):
                with self.assertRaises(ValueError):
                    write_trial_archive(self.output, row, arrays)
        self.assertFalse(self.output.exists())

    def test_loader_rejects_unsafe_paths_schemas_shapes_and_duplicate_ids(self):
        entry, manifest = self.save()
        bad_entries = []
        for key, value in (("schema_version", 2), ("schema_version", True),
                           ("path", "../outside.npz"), ("samples", 5)):
            bad_entries.append(dict(entry, **{key: value}))
        descriptor = deepcopy(entry)
        descriptor["arrays"]["point__flat_joint"]["dtype"] = "|O"
        bad_entries.append(descriptor)
        for bad in bad_entries:
            with self.subTest(entry=bad):
                with self.assertRaises(ValueError):
                    load_trial_archive(self.output, bad)
        manifest["entries"].append(entry)
        (self.output / "trial_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            load_trial_manifest(self.output)

    def test_loader_checks_actual_shapes_after_matching_hash(self):
        entry, _ = self.save()
        archive = self.output / entry["path"]
        changed = dict(self.arrays, replicate_id=np.arange(6))
        changed["point__comparator"] = np.ones(5)
        np.savez_compressed(archive, **changed)
        entry["sha256"] = hashlib.sha256(archive.read_bytes()).hexdigest()
        with self.assertRaisesRegex(ValueError, "shape"):
            load_trial_archive(self.output, entry)

    def test_old_runs_and_existing_outputs_fail_without_overwriting(self):
        self.output.mkdir(parents=True)
        original = '{"rows": []}'
        (self.output / "results.json").write_text(original)
        with self.assertRaisesRegex(ValueError, "aggregate-only or legacy"):
            rescore_run(self.output, self.output / "new.json")
        self.save()
        with self.assertRaises(FileExistsError):
            rescore_run(self.output, self.output / "results.json")
        with self.assertRaisesRegex(ValueError, "config"):
            rescore_run(self.output, self.output / "config.json")
        self.assertEqual((self.output / "results.json").read_text(), original)
        self.assertFalse((self.output / "new.json").exists())
        for factor in (0., 1., float("nan"), float("inf"), True):
            with self.assertRaises(ValueError):
                rescore_run(self.output, self.output / "new.json", factor_tolerance=factor)

    def test_cli_rescores_without_importing_inference(self):
        self.save()
        invocation = [sys.executable, "-B", str(Path(__file__).with_name("trial_archive.py")),
                      "--run", str(self.output), "--output", str(self.output / "cli.json"),
                      "--factor-tolerance", "1.5", "--cap-factor", "4"]
        result = subprocess.run(invocation, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        report = json.loads((self.output / "cli.json").read_text())
        self.assertEqual(report["factor_tolerance"], 1.5)
        self.assertEqual(report["cap_factor"], 4.)
        code = ("import sys; import trial_archive; "
                "assert 'comparison_posterior' not in sys.modules; "
                "assert 'compare_estimators' not in sys.modules")
        imports = subprocess.run([sys.executable, "-B", "-c", code],
                                 cwd=Path(__file__).parent, text=True, capture_output=True)
        self.assertEqual(imports.returncode, 0, imports.stderr)


if __name__ == "__main__":
    unittest.main()
