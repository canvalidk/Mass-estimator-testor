"""Run and record the behavioral suite, including code hashes and failures."""
from pathlib import Path
import hashlib
import json
import platform
import sys
import time
import unittest
import numpy as np

root = Path(__file__).resolve().parent
started = time.monotonic()
suite = unittest.defaultTestLoader.discover(str(root), pattern="test_*.py")
result = unittest.TextTestRunner(verbosity=2).run(suite)
record = dict(tests_run=result.testsRun, passed=result.wasSuccessful(),
              failures=[dict(test=str(test), traceback=tb) for test, tb in result.failures],
              errors=[dict(test=str(test), traceback=tb) for test, tb in result.errors],
              skipped=[dict(test=str(test), reason=reason) for test, reason in result.skipped],
              python=platform.python_version(), numpy=np.__version__,
              elapsed_seconds=time.monotonic()-started,
              source_sha256={name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                             for name in ("estimator.py","test_estimator.py","experiments.py",
                                          "uncertainty_walkthrough.py","test_uncertainty_walkthrough.py",
                                          "practical_formulas.py","test_practical_formulas.py",
                                          "practical_examples.py")})
out=root/"results"
out.mkdir(exist_ok=True)
(out/"tests.json").write_text(json.dumps(record,indent=2),encoding="utf-8")
sys.exit(0 if result.wasSuccessful() else 1)
