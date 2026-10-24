"""HTTP boundary and asynchronous run lifecycle checks for the local browser."""
import io
import json
from pathlib import Path
import threading
import time
import unittest
from unittest.mock import Mock, patch
from urllib.error import HTTPError
from urllib.request import Request, urlopen
import uuid

from experiment_config import normalize_config
from results_app import ResultsServer, WORKSPACE


class ResultsAppTests(unittest.TestCase):
    def setUp(self):
        self.workspace = WORKSPACE / ".tools" / "mass_results_tests" / uuid.uuid4().hex
        self.workspace.mkdir(parents=True)
        self.server = ResultsServer(("127.0.0.1", 0), self.workspace)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.base = f"http://127.0.0.1:{self.server.server_port}"

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join(timeout=3)

    def request(self, path, data=None, headers=None):
        headers = dict(headers or {})
        if data is not None:
            headers.setdefault("Content-Type", "application/json")
        req = Request(self.base + path, data=json.dumps(data).encode() if data is not None else None,
                      headers=headers)
        try:
            with urlopen(req, timeout=5) as response:
                return response.status, response.read(), response.headers
        except HTTPError as exc:
            return exc.code, exc.read(), exc.headers

    def post(self, data):
        return self.request("/api/runs", data, {"X-Session-Token": self.server.token})

    def test_readonly_library_and_static_routes(self):
        status, body, headers = self.request("/api/runs")
        self.assertEqual(status, 200)
        self.assertEqual(json.loads(body)["runs"], [])
        self.assertEqual(headers["Cache-Control"], "no-store")
        self.assertEqual(self.request("/api/runs/not-a-run")[0], 404)
        self.assertEqual(self.request("/../compare_estimators.py")[0], 404)
        self.assertEqual(self.request("/api/runs/%2E%2E/raw")[0], 404)

    def test_other_websites_cannot_start_experiments(self):
        self.assertEqual(self.request("/api/runs", {})[0], 403)
        self.assertEqual(self.request("/api/session", headers={"Host": "evil.example"})[0], 403)
        self.assertEqual(self.request("/api/session", headers={"Origin": "https://evil.example"})[0], 403)
        self.assertEqual(self.request("/api/session")[0], 200)

    def test_invalid_configuration_never_creates_run(self):
        self.assertEqual(self.post({"noise_model": "diagonal"})[0], 400)
        self.assertEqual(self.post({"samples": True})[0], 400)
        self.assertEqual(self.post({"samples": 1})[0], 400)
        self.assertFalse((self.workspace / ".tools" / "mass_estimator_runs").exists())

    def test_one_worker_progress_and_saved_status(self):
        entered, release = threading.Event(), threading.Event()

        def fake_run(config, output, progress):
            entered.set()
            if not release.wait(5):
                raise RuntimeError("test worker timed out")
            progress({"completed": "00_custom", "completed_count": 1, "seconds": 0.1})
            return {"rows": [{}], "elapsed_seconds": 0.1}

        with patch("results_app.RunManager._run_worker", side_effect=fake_run):
            status, body, _ = self.post({"samples": 2, "scenario": {"type": "pairs", "pairs": [[8, 2]]}})
            self.assertEqual(status, 202)
            job_id = json.loads(body)["job"]["id"]
            self.assertTrue(entered.wait(3))
            try:
                self.assertEqual(self.post({"samples": 2})[0], 409)
                self.assertEqual(json.loads(self.request("/api/job")[1])["job"]["status"], "running")
                run_dir = self.workspace / ".tools" / "mass_estimator_runs" / job_id
                config = json.loads((run_dir / "config.json").read_text())
                self.assertEqual(config["scenario"]["pairs"], [[8, 2]])
            finally:
                release.set()
            self.wait_for_job("complete")
            saved = json.loads((run_dir / "status.json").read_text())
            self.assertEqual(saved["completed"], 1)

    def wait_for_job(self, expected):
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline:
            job = self.server.manager.snapshot()
            if job and job["status"] == expected:
                return job
            time.sleep(0.02)
        self.fail(f"Job did not reach {expected}")

    def test_worker_failure_is_explicit_and_persisted(self):
        with patch("results_app.RunManager._run_worker", side_effect=FloatingPointError("refinement failed")):
            status, body, _ = self.post({"samples": 2})
            self.assertEqual(status, 202)
            job = self.wait_for_job("failed")
            self.assertIn("refinement failed", job["error"])
            run_dir = self.workspace / ".tools" / "mass_estimator_runs" / job["id"]
            saved = json.loads((run_dir / "status.json").read_text())
            self.assertEqual(saved["status"], "failed")
            self.assertFalse((run_dir / "results.json").exists())

    def test_startup_write_failure_allows_retry(self):
        with patch("results_app.write_json", side_effect=OSError("disk unavailable")):
            status, body, _ = self.post({"samples": 2})
            self.assertEqual(status, 500)
            self.assertIn("disk unavailable", json.loads(body)["error"])
        failed = self.server.manager.snapshot()
        self.assertEqual(failed["status"], "failed")
        self.assertIn("status_save_error", failed)
        with patch("results_app.RunManager._run_worker", return_value={"rows": [{}], "elapsed_seconds": .1}):
            self.assertEqual(self.post({"samples": 2})[0], 202)
            retried = self.wait_for_job("complete")
        self.assertNotEqual(retried["id"], failed["id"])

    def test_thread_start_failure_is_persisted_and_does_not_block_retry(self):
        with patch("results_app.threading.Thread.start", side_effect=RuntimeError("cannot start thread")):
            with self.assertRaisesRegex(OSError, "cannot start thread"):
                self.server.manager.start({"samples": 2})
        job = self.server.manager.snapshot()
        self.assertEqual(job["status"], "failed")
        saved = self.workspace / ".tools" / "mass_estimator_runs" / job["id"] / "status.json"
        self.assertEqual(json.loads(saved.read_text())["status"], "failed")
        with patch("results_app.RunManager._run_worker", return_value={"rows": [{}], "elapsed_seconds": .1}):
            self.server.manager.start({"samples": 2})
            self.wait_for_job("complete")

    def test_real_worker_saves_current_config_and_source_provenance(self):
        status, body, _ = self.post({"name": "Actual child", "samples": 2,
                                    "scenario": {"type": "pairs", "pairs": [[8, 2]]}})
        self.assertEqual(status, 202)
        job = self.wait_for_job("complete")
        run_dir = self.workspace / ".tools" / "mass_estimator_runs" / job["id"]
        report = json.loads((run_dir / "results.json").read_text(encoding="utf-8"))
        self.assertEqual(report["experiment_config"], job["config"])
        self.assertEqual(len(report["rows"]), 1)
        self.assertEqual(job["completed"], 1)
        self.assertIn("compare_estimators.py", report["source_sha256"])
        self.assertTrue((run_dir / "worker.log").exists())
        self.assertIsNone(self.server.manager.process)

    def test_real_worker_publishes_repeated_readings_and_calibration_config(self):
        measurement = {"repeats": 4, "calibration": {"mode": "estimated", "samples": 20, "seed": 37}}
        status, _, _ = self.post({"name": "Calibrated child", "samples": 2,
                                 "scenario": {"type": "pairs", "pairs": [[8, 2]]},
                                 "measurement": measurement})
        self.assertEqual(status, 202)
        job = self.wait_for_job("complete")
        run_dir = self.workspace / ".tools" / "mass_estimator_runs" / job["id"]
        report = json.loads((run_dir / "results.json").read_text(encoding="utf-8"))
        self.assertEqual(report["experiment_config"]["measurement"], measurement)
        self.assertEqual(report["measurement"]["repeats"], 4)
        row = report["rows"][0]
        self.assertEqual((row["true_force_snr"], row["effective_force_snr"]), (8., 16.))
        self.assertEqual(len(row["methods"]), 16)
        self.assertIsNone(self.server.manager.process)

    def test_changed_provenance_is_retained_but_not_published_as_a_run(self):
        output = self.workspace / "unverified"
        output.mkdir()
        config = normalize_config({"samples": 2})
        (output / "config.json").write_text(json.dumps(config))
        (output / "results.json").write_text(json.dumps({"experiment_config": config, "source_sha256": {}}))
        process = Mock(stdout=io.StringIO(""))
        process.wait.return_value = 0
        process.poll.return_value = 0
        with patch("results_app.subprocess.Popen", return_value=process):
            with self.assertRaisesRegex(RuntimeError, "source changed"):
                self.server.manager._run_worker(config, output, lambda event: None)
        self.assertFalse((output / "results.json").exists())
        self.assertTrue((output / "unverified_results.json").exists())

    def test_shutdown_stops_child_and_preserves_interrupted_status(self):
        entered, release = threading.Event(), threading.Event()
        process = Mock()
        process.poll.return_value = None
        process.terminate.side_effect = release.set

        def blocked_worker(config, output, progress):
            with self.server.manager.lock:
                self.server.manager.process = process
            entered.set()
            release.wait(5)
            return {"rows": [{}], "elapsed_seconds": .1}

        with patch("results_app.RunManager._run_worker", side_effect=blocked_worker):
            self.server.manager.start({"samples": 2})
            self.assertTrue(entered.wait(3))
            self.server.manager.close()
        process.terminate.assert_called_once()
        job = self.server.manager.snapshot()
        self.assertEqual(job["status"], "failed")
        self.assertIn("InterruptedError", job["error"])
        status_path = self.workspace / ".tools" / "mass_estimator_runs" / job["id"] / "status.json"
        self.assertEqual(json.loads(status_path.read_text())["status"], "failed")
        self.assertFalse(self.server.manager.thread.is_alive())
        # Keep the shutdown fixture idempotent after this mock process ended.
        process.poll.return_value = 0
        with self.assertRaisesRegex(RuntimeError, "shutting down"):
            self.server.manager.start({"samples": 2})


if __name__ == "__main__":
    unittest.main()
