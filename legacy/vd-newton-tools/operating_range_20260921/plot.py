"""Render the verified discrete accuracy map; no inference or interpolation."""
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
source = HERE / "analysis/accuracy.json"
report = json.loads(source.read_text(encoding="utf-8"))
sys.path.append(str(ROOT / ".tools/assessment_plot_dependencies"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, TwoSlopeNorm
from matplotlib.patches import Patch

output = HERE / "figures"
# This renderer updates only its own named figure/provenance files; immutable
# simulation and accuracy inputs are never written.
rows = report["experiments"][0]["rows"]
force_grid = sorted({r["scenario"]["true_force_snr"] for r in rows})
accel_grid = sorted({r["scenario"]["true_acceleration_snr"] for r in rows})
assert len(rows) == len(force_grid) * len(accel_grid)
lookup = {(r["scenario"]["true_force_snr"],r["scenario"]["true_acceleration_snr"]): r for r in rows}
colors = ["#f1ddda", "#f4e7b5", "#cae6d6"]
classes = {"below_target": 0, "unresolved": 1, "supported": 2}
fig, axes = plt.subplots(2, 2, figsize=(13, 12))
for ax, factor in zip(axes.ravel()[:3], ["1.25", "1.5", "2"]):
    values = np.zeros((len(force_grid),len(accel_grid)))
    labels = values.copy()
    for i,f in enumerate(force_grid):
        for j,a in enumerate(accel_grid):
            metric = lookup[(f,a)]["methods"]["flat_joint"]["accuracy"][factor]
            assert metric["status"] == "available"
            values[i,j] = classes[metric["targets"]["0.95"]]
            labels[i,j] = metric["success_probability"] * 100
    ax.imshow(values, origin="lower", interpolation="nearest",
              cmap=ListedColormap(colors), norm=BoundaryNorm([-.5,.5,1.5,2.5],3))
    for (i,j),p in np.ndenumerate(labels):
        ax.text(j,i,f"{p:.1f}",ha="center",va="center",fontsize=9)
    ax.set_title(f"Our success within factor {factor}\nColor: assessment of a 95% success target",fontsize=11)

ax = axes.ravel()[3]
differences = np.zeros((len(force_grid),len(accel_grid)))
signs = {}
for i,f in enumerate(force_grid):
    for j,a in enumerate(accel_grid):
        pair=lookup[(f,a)]["paired_flat_minus_other"]["norm_ratio"]["1.5"]
        differences[i,j] = pair["mean_difference"] * 100
        band=pair["simultaneous_failure_difference_interval"]
        signs[(i,j)] = band["upper"] < 0 or band["lower"] > 0
norm=TwoSlopeNorm(vmin=-25,vcenter=0,vmax=25)
im=ax.imshow(differences,origin="lower",interpolation="nearest",cmap="RdBu_r",norm=norm)
for (i,j),d in np.ndenumerate(differences):
    ax.text(j,i,f"{d:+.1f}"+("*" if signs[(i,j)] else ""),ha="center",va="center",fontsize=8.5,
            color="white" if abs(d)>17 else "#202020")
ax.set_title("Failure difference at factor 1.5 (pp)\nOur rate minus norm ratio: negative favors ours",fontsize=11)
fig.colorbar(im,ax=ax,shrink=.78,pad=.02,label="Percentage points")
for ax in axes.ravel():
    ax.set_xticks(range(len(accel_grid)),[f"{x:g}" for x in accel_grid])
    ax.set_yticks(range(len(force_grid)),[f"{x:g}" for x in force_grid])
    ax.set_xlabel("Effective acceleration SNR")
    ax.set_ylabel("Effective force SNR")
    ax.set_xticks(np.arange(-.5,len(accel_grid),1),minor=True)
    ax.set_yticks(np.arange(-.5,len(force_grid),1),minor=True)
    ax.grid(which="minor",color="white",linewidth=1.5)
    ax.tick_params(which="minor",bottom=False,left=False)
fig.suptitle("Fresh operating-range study: absolute accuracy and relative advantage",fontsize=16,y=.99)
fig.legend(handles=[Patch(facecolor=colors[2],label="Target supported"),
                    Patch(facecolor=colors[1],label="Unresolved"),
                    Patch(facecolor=colors[0],label="Below target")],
           loc="lower center",bbox_to_anchor=(.5,.065),ncol=3,frameon=False)
fig.text(.5,.025,"Numbers in accuracy panels are observed success percentages; classification uses simultaneous bounds.\n"
         "10,000 fresh trials/cell. One 95% Monte Carlo event: marginal ±2.263 pp, differences ±4.527 pp; * excludes zero.\n"
         "Known isotropic Gaussian noise, fixed true cells. Discrete spacing; no interpolation or observed-data acceptance rule.",
         ha="center",va="bottom",fontsize=9,color="#444444")
fig.tight_layout(rect=(0,.105,1,.955),h_pad=2,w_pad=1.6)
output.mkdir(exist_ok=True)
for ext in ("png","svg"):
    fig.savefig(output/f"operating-range.{ext}",dpi=170,facecolor="white")
plt.close(fig)
(output/"provenance.json").write_text(json.dumps({"accuracy_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
    "plot_source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "matplotlib":matplotlib.__version__,"numpy":np.__version__,"interpolation":"none; displayed discrete grid cells"},indent=2)+"\n",encoding="utf-8")
print(output/"operating-range.png")
