"""Plot selected archived error curves; no estimator fitting."""
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "03_trace_and_evaluator/mass_estimator/investigation"))
import numpy as np
sys.path.append(str(ROOT / ".tools/assessment_plot_dependencies"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from trial_archive import load_trial_manifest, load_trial_archive, _custom_losses
from comparison_scores import paired_summary

output = Path(__file__).resolve().parent / "main-reviewed"
report = json.loads((output / "analysis.json").read_text(encoding="utf-8"))
destination = output / "tolerance-comparison"
if destination.with_suffix(".png").exists() or destination.with_suffix(".svg").exists():
    raise FileExistsError(destination)
specs = [(4, .25, 2., "Unequal signal strengths"),
         (4, 1., 2., "Balanced signal strengths"),
         (1, .25, .5, "Very weak signals")]
factors = np.geomspace(1.05, 10., 401)
fig, axes = plt.subplots(2, 3, figsize=(13.8, 8), sharex=True)
for col, (n, m, snr, title) in enumerate(specs):
    row = next(r for r in report["rows"] if r["scenario"]["repeats"] == n
               and r["scenario"]["true_mass"] == m and r["scenario"]["true_acceleration_snr"] == snr)
    run = Path(row["archive_path"]).parent.parent
    manifest = load_trial_manifest(run)
    entry = next(e for e in manifest["entries"] if e["id"] == row["id"])
    arrays = load_trial_archive(run, entry)
    assert entry["sha256"] == row["archive_sha256"]
    flat = arrays["point__flat_joint"]
    norm = arrays["point__norm_ratio"]
    fa, fb, differences, mcse = [], [], [], []
    for factor in factors:
        a = _custom_losses(flat, m, float(factor), 2.)["outside_factor"]
        b = _custom_losses(norm, m, float(factor), 2.)["outside_factor"]
        pair = paired_summary(a, b)
        fa.append(np.mean(a) * 100)
        fb.append(np.mean(b) * 100)
        differences.append(pair["mean_difference"] * 100)
        mcse.append(pair["mcse"] * 100)
    d, se = np.array(differences), np.array(mcse)
    width = row["uniform_difference_halfwidth"] * 100
    top, bottom = axes[:, col]
    top.plot(factors, fa, color="#195cab", lw=2.1, label="Our ratio of means")
    top.plot(factors, fb, color="#b76226", lw=2.1, label="Simple norm ratio")
    top.set_title(f"{title}\nn={n}, mass={m:g}, per-reading a-SNR={snr:g}", fontsize=11)
    top.set_ylim(-2, 102)
    bottom.fill_between(factors, d-width, d+width, color="#d4d8de", label="95% simultaneous bound")
    bottom.fill_between(factors, d-1.96*se, d+1.96*se, color="#80aeda", alpha=.85,
                        label="Pointwise ±1.96 paired MCSE")
    bottom.plot(factors, d, color="#183752", lw=1.6)
    bottom.axhline(0, color="#555555", lw=.8)
    bottom.set_ylim(-22, 14)
    bottom.set_xlabel("Allowed multiplicative error (factor K)")
    for ax in (top, bottom):
        ax.set_xscale("log")
        ax.set_xticks([1.1, 1.5, 2, 3, 5, 10], labels=["1.1", "1.5", "2", "3", "5", "10"])
        ax.grid(alpha=.18)
        ax.spines[["top", "right"]].set_visible(False)
axes[0, 0].set_ylabel("Outside tolerance (%)")
axes[1, 0].set_ylabel("Our failure rate − norm-ratio rate (pp)\nNegative favors our estimator")
axes[0, 0].legend(fontsize=9, loc="lower left")
axes[1, 0].legend(fontsize=8, loc="lower right")
fig.suptitle("The advantage depends on both signal conditions and error tolerance", fontsize=16, y=.98)
fig.text(.5, .035, "10,000 matched archived trials per condition; known isotropic Gaussian noise; no refitting.\n"
         "Success: true mass / K ≤ estimate ≤ true mass × K. Gray bounds cover all 24 cells, all archived methods and thresholds.\n"
         "Bounds concern Monte Carlo error only. Curves do not establish performance at other signal conditions.",
         ha="center", va="bottom", fontsize=9, color="#444444")
fig.tight_layout(rect=(0, .12, 1, .94))
fig.savefig(destination.with_suffix(".png"), dpi=160)
fig.savefig(destination.with_suffix(".svg"))
plt.close(fig)
provenance = {"analysis_sha256": hashlib.sha256((output/"analysis.json").read_bytes()).hexdigest(),
              "plot_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "numpy": np.__version__, "matplotlib": matplotlib.__version__,
              "selection": specs, "curve_factors": factors.tolist()}
(output/"plot_provenance.json").write_text(json.dumps(provenance, indent=2)+"\n", encoding="utf-8")
(output/"plot-source.py").write_bytes(Path(__file__).read_bytes())
print(destination.with_suffix(".png"))
