"""Minimal scientific inputs, with annotations and investigator truth kept separate.

Validation protects mathematical meaning, not an exhaustive document schema.
Unknown descriptive fields round-trip in extras; they are listed in the report
and have no effect on the calculation. Adapters explicitly accept or reject the
scientific declarations they understand, such as the prior family.
"""
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml


@dataclass
class Observation:
    id: str
    force: np.ndarray
    acceleration: np.ndarray
    covariance: np.ndarray
    subspace: np.ndarray = None
    window: tuple = None
    metadata: dict = field(default_factory=dict)
    extras: dict = field(default_factory=dict)

    def __post_init__(self):
        self.force = np.array(self.force, dtype=float, copy=True)
        self.acceleration = np.array(self.acceleration, dtype=float, copy=True)
        if (self.force.ndim != 1 or not self.force.size or self.acceleration.shape != self.force.shape
                or not np.all(np.isfinite(self.force)) or not np.all(np.isfinite(self.acceleration))):
            raise ValueError("force and acceleration must be finite nonempty vectors of equal length")
        d = self.force.size
        self.covariance = np.array(self.covariance, dtype=float, copy=True)
        if (self.covariance.shape != (2 * d, 2 * d) or not np.all(np.isfinite(self.covariance))
                or not np.allclose(self.covariance, self.covariance.T, rtol=1e-12, atol=0)):
            raise ValueError("covariance must be a finite symmetric 2d by 2d matrix")
        try:
            np.linalg.cholesky(self.covariance)
        except np.linalg.LinAlgError as error:
            raise ValueError("covariance must be positive definite") from error
        self.subspace = np.eye(d) if self.subspace is None else np.array(self.subspace, dtype=float, copy=True)
        if self.subspace.ndim == 1:
            self.subspace = self.subspace[:, None]
        n = self.subspace
        if (n.ndim != 2 or n.shape[0] != d or not 1 <= n.shape[1] <= d or not np.all(np.isfinite(n))
                or not np.allclose(n.T @ n, np.eye(n.shape[1]), rtol=1e-10, atol=1e-12)):
            raise ValueError("subspace must have orthonormal columns")
        if self.window is not None:
            window = np.asarray(self.window, dtype=float)
            if window.shape != (2,) or not np.all(np.isfinite(window)) or window[0] >= window[1]:
                raise ValueError("window must be a finite [start, end] with start < end")
            self.window = tuple(float(x) for x in window)
        self.id = str(self.id)
        self.metadata, self.extras = deepcopy(self.metadata), deepcopy(self.extras)

    def to_dict(self):
        return dict(deepcopy(self.extras), id=self.id, force=self.force.tolist(),
                    acceleration=self.acceleration.tolist(), covariance=self.covariance.tolist(),
                    subspace=self.subspace.tolist(), window=self.window, metadata=deepcopy(self.metadata))

    @classmethod
    def from_dict(cls, data, index=0):
        if not isinstance(data, dict):
            raise ValueError("each reading must be a mapping")
        known = {"id", "force", "acceleration", "covariance", "subspace", "window", "metadata"}
        missing = {"force", "acceleration", "covariance"} - data.keys()
        if missing:
            raise ValueError(f"reading is missing {sorted(missing)}")
        return cls(str(data.get("id", f"r{index:03d}")), data["force"], data["acceleration"],
                   data["covariance"], data.get("subspace"), data.get("window"),
                   data.get("metadata", {}), {k: deepcopy(v) for k, v in data.items() if k not in known})


@dataclass
class Inputs:
    """The entire estimator-facing boundary: observations and a declared prior."""
    readings: tuple
    prior: dict

    def to_dict(self):
        return dict(readings=[r.to_dict() for r in self.readings], prior=deepcopy(self.prior))


@dataclass
class Case:
    name: str
    question: str
    inputs: Inputs
    investigator_truth: dict = None
    metadata: dict = field(default_factory=dict)
    extras: dict = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data):
        if not isinstance(data, dict) or not {"name", "readings", "prior"} <= data.keys():
            raise ValueError("a case needs name, readings and an explicit prior")
        if not isinstance(data["name"], str) or not data["name"].strip():
            raise ValueError("case name must be nonempty")
        if not isinstance(data["readings"], list) or not data["readings"]:
            raise ValueError("a case needs a nonempty reading list")
        readings = tuple(Observation.from_dict(r, i) for i, r in enumerate(data["readings"]))
        if len({r.id for r in readings}) != len(readings):
            raise ValueError("reading ids must be unique within a case")
        if not isinstance(data["prior"], dict):
            raise ValueError("prior must be a declared mapping")
        known = {"name", "question", "readings", "prior", "truth", "metadata", "kind", "seed", "numerics", "checks"}
        return cls(data["name"], str(data.get("question", "")), Inputs(readings, deepcopy(data["prior"])),
                   deepcopy(data.get("truth")), deepcopy(data.get("metadata", {})),
                   {k: deepcopy(v) for k, v in data.items() if k not in known})

    def to_dict(self, *, include_truth=False):
        out = dict(deepcopy(self.extras), name=self.name, question=self.question,
                   **self.inputs.to_dict(), metadata=deepcopy(self.metadata))
        if include_truth and self.investigator_truth is not None:
            out["truth"] = deepcopy(self.investigator_truth)
        return out


def parse_cases(document):
    """Accept one case or a collection; the collection envelope stays intentionally small."""
    if not isinstance(document, dict):
        raise ValueError("inspection input must be a mapping")
    rows = document.get("cases", [document])
    if not isinstance(rows, list) or not rows:
        raise ValueError("cases must be a nonempty list")
    cases = [Case.from_dict(row) for row in rows]
    if len({c.name for c in cases}) != len(cases):
        raise ValueError("case names must be unique")
    return cases


def load_cases(path):
    with Path(path).open(encoding="utf-8") as handle:
        document = yaml.safe_load(handle)
    return document, parse_cases(document)
