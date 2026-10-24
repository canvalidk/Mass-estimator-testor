"""Local results browser and configuration-driven experiment launcher.

Run from any directory: python -B results_app.py
Uses the existing comparison engine, Python's HTTP server, and NumPy. The
server listens on loopback only and keeps generated runs under .tools/.
"""
import argparse
from datetime import datetime, timezone
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import secrets
import subprocess
import sys
import threading
from urllib.parse import unquote, urlsplit
import uuid

from experiment_config import normalize_config, scenario_rows
from results_catalog import list_runs, load_run

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[2]
STATIC = HERE / "results_ui"
MAX_BODY = 65536
WORKER_SOURCES = ("compare_estimators.py", "comparison_scores.py", "comparison_posterior.py",
                  "comparison_estimators.py", "crossover_vector_checks.py", "experiment_config.py",
                  "comparison_diagonal_posterior.py", "comparison_diagonal_baselines.py", "trial_archive.py",
                  "calibration_posterior.py")


def json_bytes(value):
    return (json.dumps(value, allow_nan=False, ensure_ascii=False) + "\n").encode("utf-8")


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(json_bytes(value))
    temporary.replace(path)


class RunManager:
    """One calculation at a time; HTTP requests stay responsive during a run."""

    def __init__(self, workspace):
        self.workspace = Path(workspace).resolve()
        self.lock = threading.Lock()
        self.job = None
        self.process = None
        self.thread = None
        self.closed = False

    def snapshot(self):
        with self.lock:
            return dict(self.job) if self.job else None

    def start(self, data):
        config = normalize_config(data)
        total = len(scenario_rows(config))
        with self.lock:
            if self.closed:
                raise RuntimeError("The results server is shutting down")
            if self.job and self.job["status"] == "running":
                raise RuntimeError("An experiment is already running. Let it finish before starting another.")
            run_id = datetime.now(timezone.utc).strftime("run-%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:8]
            output = self.workspace / ".tools" / "mass_estimator_runs" / run_id
            job = {"id": run_id, "run_id": run_id, "status": "running",
                   "completed": 0, "total": total, "config": config,
                   "message": "Preparing the first signal condition",
                   "started_at": datetime.now(timezone.utc).isoformat()}
            self.job = job
            try:
                output.mkdir(parents=True, exist_ok=False)
                write_json(output / "config.json", config)
                write_json(output / "status.json", job)
                snapshot = dict(job)
                self.thread = threading.Thread(target=self._execute, args=(job, output),
                                               daemon=True, name="mass-comparison")
                self.thread.start()
            except Exception as exc:
                self._fail(job, output, exc, "Experiment could not start")
                raise OSError(f"Experiment could not start: {exc}") from exc
        return snapshot

    def _fail(self, job, output, error, message="Experiment failed; inspect worker.log for details"):
        """Called under the manager lock; preserve any earlier shutdown error."""
        if job["status"] != "running":
            return
        job.update(status="failed", error=f"{type(error).__name__}: {error}", message=message,
                   finished_at=datetime.now(timezone.utc).isoformat())
        try:
            write_json(output / "status.json", job)
        except OSError as exc:
            # A disk failure must not leave the in-memory manager stuck in
            # 'running', nor hide the original error behind another exception.
            job["status_save_error"] = str(exc)

    def _execute(self, job, output):
        def progress(event):
            with self.lock:
                if job["status"] == "running":
                    job.update(completed=event["completed_count"],
                               total=event.get("total", job["total"]),
                               message=f"Completed {event['completed']}", seconds=event["seconds"])
                    write_json(output / "status.json", job)

        try:
            report = self._run_worker(job["config"], output, progress)
            with self.lock:
                if job["status"] == "running":
                    completed = dict(job, status="complete", completed=len(report["rows"]), total=len(report["rows"]),
                                     message="Results saved", seconds=report["elapsed_seconds"],
                                     finished_at=datetime.now(timezone.utc).isoformat())
                    # Publish the completed in-memory state only after its
                    # durable status succeeds, keeping startup races excluded.
                    write_json(output / "status.json", completed)
                    job.update(completed)
        except Exception as exc:
            with self.lock:
                self._fail(job, output, exc)

    @staticmethod
    def _stop_process(process):
        if process is not None and process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=3)

    def _run_worker(self, config, output, progress):
        """Each experiment imports the current engine in a fresh interpreter."""
        source_hashes = {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                         for name in WORKER_SOURCES}
        process = None
        log_path = output / "worker.log"
        try:
            with log_path.open("w", encoding="utf-8") as log:
                with self.lock:
                    if self.closed:
                        raise InterruptedError("Results server stopped before the calculation began")
                    process = subprocess.Popen(
                        [sys.executable, "-B", str(HERE / "compare_estimators.py"),
                         "--config", str(output / "config.json"), "--output", str(output)],
                        cwd=str(HERE), stdout=subprocess.PIPE, stderr=log, text=True,
                        encoding="utf-8", errors="replace",
                        creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
                    self.process = process
                for line in process.stdout:
                    try:
                        event = json.loads(line)
                    except ValueError:
                        continue
                    if isinstance(event, dict) and "completed_count" in event:
                        progress(event)
                code = process.wait()
                if code:
                    log.flush()
                    with log_path.open("rb") as errors:
                        errors.seek(0, os.SEEK_END)
                        errors.seek(max(0, errors.tell() - 4096))
                        detail = errors.read().decode("utf-8", errors="replace").strip()
                    raise RuntimeError(f"Calculation exited with code {code}: {detail or 'see worker.log'}")
            result_path = output / "results.json"
            report = json.loads(result_path.read_text(encoding="utf-8"))
            if report.get("experiment_config") != config or report.get("source_sha256") != source_hashes:
                # Preserve the calculation for inspection, but do not publish
                # a run whose executed configuration/source cannot be verified.
                result_path.replace(output / "unverified_results.json")
                raise RuntimeError("Configuration or engine source changed during calculation; rerun the experiment")
            return report
        finally:
            self._stop_process(process)
            if process is not None and process.stdout is not None:
                process.stdout.close()
            with self.lock:
                if self.process is process:
                    self.process = None

    def close(self):
        """Persist interruption before stopping an active calculation."""
        with self.lock:
            self.closed = True
            job, process, thread = self.job, self.process, self.thread
            if job and job["status"] == "running":
                output = self.workspace / ".tools" / "mass_estimator_runs" / job["id"]
                self._fail(job, output, InterruptedError("Results server was stopped"),
                           "Experiment interrupted when the results server stopped; its configuration is retained")
        self._stop_process(process)
        if thread is not None and thread.ident is not None and thread is not threading.current_thread():
            thread.join(timeout=5)


class ResultsServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, address, workspace=WORKSPACE):
        super().__init__(address, ResultsHandler)
        self.workspace = Path(workspace).resolve()
        self.manager = RunManager(self.workspace)
        self.token = secrets.token_urlsafe(32)

    def server_close(self):
        if hasattr(self, "manager"):
            self.manager.close()
        super().server_close()


class ResultsHandler(BaseHTTPRequestHandler):
    server_version = "MassResults/1.0"

    def log_message(self, fmt, *args):
        # Keep the terminal useful: experiment progress and HTTP errors only.
        if len(args) >= 2 and str(args[1]).startswith(("4", "5")):
            super().log_message(fmt, *args)

    def _local_request(self):
        port = self.server.server_port
        allowed = {f"127.0.0.1:{port}", f"localhost:{port}"}
        if self.headers.get("Host", "").lower() not in allowed:
            self._json({"error": "Only this local server's address is accepted"}, 403)
            return False
        origin = self.headers.get("Origin")
        if origin and origin not in {f"http://{host}" for host in allowed}:
            self._json({"error": "Cross-origin requests are not accepted"}, 403)
            return False
        return True

    def _send(self, payload, content_type, status=200, filename=None):
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Content-Security-Policy",
                         "default-src 'self'; script-src 'self'; style-src 'self'; "
                         "img-src 'self' data:; connect-src 'self'; object-src 'none'; "
                         "base-uri 'none'; frame-ancestors 'none'")
        if filename:
            self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
        self.end_headers()
        self.wfile.write(payload)

    def _json(self, value, status=200, filename=None):
        self._send(json_bytes(value), "application/json; charset=utf-8", status, filename)

    def do_GET(self):
        if not self._local_request():
            return
        path = unquote(urlsplit(self.path).path)
        try:
            if path == "/api/session":
                self._json({"token": self.server.token})
            elif path == "/api/runs":
                self._json({"runs": list_runs(self.server.workspace),
                            "active": self.server.manager.snapshot()})
            elif path == "/api/job":
                self._json({"job": self.server.manager.snapshot()})
            elif path.startswith("/api/runs/"):
                parts = path.split("/")
                if len(parts) not in (4, 5) or not parts[3]:
                    raise KeyError(path)
                run_id = parts[3]
                if len(parts) == 4:
                    self._json(load_run(self.server.workspace, run_id))
                elif parts[4] == "config":
                    run = load_run(self.server.workspace, run_id)
                    if run.get("config") is None:
                        self._json({"error": "This historical study has no compatible saved configuration"}, 404)
                    else:
                        self._json(run["config"], filename=f"{run_id}-config.json")
                elif parts[4] == "raw":
                    entry = next((r for r in list_runs(self.server.workspace) if r["id"] == run_id), None)
                    if entry is None:
                        raise KeyError(run_id)
                    source = (self.server.workspace / entry["source"]).resolve()
                    if not source.is_relative_to((self.server.workspace / ".tools").resolve()):
                        raise ValueError("Result source must remain in the run library")
                    self._send(source.read_bytes(), "application/json", filename=f"{run_id}-results.json")
                else:
                    raise KeyError(path)
            elif path in ("/", "/index.html", "/app.js", "/style.css"):
                name = "index.html" if path == "/" else path[1:]
                mime = {"index.html": "text/html", "app.js": "text/javascript", "style.css": "text/css"}[name]
                self._send((STATIC / name).read_bytes(), mime + "; charset=utf-8")
            elif path == "/favicon.ico":
                self._send(b"", "image/x-icon", 204)
            else:
                raise KeyError(path)
        except (KeyError, FileNotFoundError):
            self._json({"error": "Not found"}, 404)
        except (ValueError, OSError) as exc:
            self._json({"error": str(exc)}, 400)

    def do_POST(self):
        # Consume bounded request data before an early rejection. Closing a
        # Windows socket with unread input can reset it before the error arrives.
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 < length <= MAX_BODY:
                raise ValueError("Configuration must contain between 1 and 65,536 bytes")
            payload = self.rfile.read(length)
        except ValueError as exc:
            self._json({"error": str(exc)}, 400)
            return
        if not self._local_request():
            return
        if urlsplit(self.path).path != "/api/runs":
            self._json({"error": "Not found"}, 404)
            return
        if not secrets.compare_digest(self.headers.get("X-Session-Token", ""), self.server.token):
            self._json({"error": "Refresh this local page before starting an experiment"}, 403)
            return
        if self.headers.get("Content-Type", "").split(";")[0].strip() != "application/json":
            self._json({"error": "Send an application/json configuration"}, 415)
            return
        try:
            def reject_constant(value):
                raise ValueError(f"Nonfinite JSON number: {value}")
            data = json.loads(payload, parse_constant=reject_constant)
            self._json({"job": self.server.manager.start(data)}, 202)
        except RuntimeError as exc:
            self._json({"error": str(exc)}, 409)
        except (ValueError, TypeError, UnicodeDecodeError) as exc:
            self._json({"error": str(exc)}, 400)
        except OSError as exc:
            self._json({"error": str(exc)}, 500)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8766)
    args = parser.parse_args()
    server = ResultsServer(("127.0.0.1", args.port))
    print(f"Mass estimator results: http://127.0.0.1:{server.server_port}", flush=True)
    print("Press Ctrl+C to stop. Runs and configurations are saved under .tools/mass_estimator_runs/.", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
