"""Optional, local scientific plots for the deterministic misalignment study."""

from html import escape
from io import StringIO
from pathlib import Path
import json
import re

import numpy as np


def write_report(folder: Path, data: dict, curves: dict[str, np.ndarray]) -> bool:
    """Write PNG/SVG figures and a self-contained HTML report; False without matplotlib.

    Densities supplied by the investigator are densities per unit log mass. This
    renderer performs no estimation, inference, or pass/fail classification.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError:
        return False

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    centre = float(data["prior"]["center"])
    cells = data["cells"]
    groups = {}
    for cell in cells:
        groups.setdefault(cell["case"], []).append(cell)
    for group in groups.values():
        group.sort(key=lambda row: row["noise_scale"], reverse=True)

    def fmt(value, digits=5):
        if value is None:
            return "—"
        return format(float(value), f".{digits}g")

    def boundary(cell):
        location = cell["positive_min_location"]
        return {
            "mass_to_zero": "m → 0⁺",
            "mass_to_infinity": "m → ∞",
            "both_boundaries": "m → 0⁺ or ∞",
            "all_masses": "every positive mass",
        }.get(location, "m = " + fmt(cell["positive_min_mass"]))

    def arr(cell, name):
        return np.asarray(curves[f'{cell["id"]}/{name}'])

    def save(fig, name):
        fig.savefig(folder / f"{name}.png", dpi=170, facecolor="white")
        fig.savefig(folder / f"{name}.svg", facecolor="white")
        stream = StringIO()
        fig.savefig(stream, format="svg", facecolor="white")
        plt.close(fig)
        svg = stream.getvalue()
        # Each figure has its own clip paths. Prefix IDs before embedding several
        # SVG documents in the same HTML page, where IDs otherwise collide.
        svg = svg[svg.index("<svg"):]
        ids = re.findall(r'\bid="([^"]+)"', svg)
        for ident in sorted(set(ids), key=len, reverse=True):
            svg = svg.replace(f'id="{ident}"', f'id="{name}_{ident}"')
            svg = svg.replace(f'#{ident})', f'#{name}_{ident})')
            svg = svg.replace(f'="#{ident}"', f'="#{name}_{ident}"')
        return svg

    style = {
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.titleweight": "bold", "axes.labelcolor": "#334155",
        "axes.edgecolor": "#94a3b8", "grid.alpha": .2,
        "svg.fonttype": "none", "savefig.bbox": "tight",
    }
    fragments = []
    captions = {}
    with plt.rc_context(style):
        fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), layout="constrained")
        colors = plt.get_cmap("tab10")
        for i, (name, group) in enumerate(groups.items()):
            color = colors(i % 10)
            noise = np.array([row["noise_scale"] for row in group])
            raw = np.array([row["raw"]["point"] / centre for row in group])
            quoted = np.array([row["quoted"]["point"] / centre for row in group])
            intervals = np.array([row["quoted"]["interval"] for row in group]) / centre
            axes[0, 0].fill_between(noise, intervals[:, 0], intervals[:, 1], color=color, alpha=.1)
            axes[0, 0].plot(noise, quoted, "o-", color=color, label=name, markersize=4)
            axes[0, 0].plot(noise, raw, "--", color=color, alpha=.8)
            axes[0, 1].plot(noise, [row["positive_min_cost"] for row in group], "o-", color=color)
            axes[0, 1].plot(noise, [row["projective_min_cost"] for row in group], "--", color=color)
            weights = np.array([row["misfit_weight"] for row in group])
            axes[1, 0].plot(noise, np.ma.masked_less_equal(weights, 0), "o-", color=color)
            axes[1, 1].plot(noise, [row["quoted"]["point"] / row["quoted"]["median"]
                                  for row in group], "o-", color=color)
        for ax in axes.flat:
            ax.set_xscale("log")
            ax.invert_xaxis()
            ax.set_xlabel("Declared noise scale (precision increases →)")
            ax.grid(True, which="major")
        axes[0, 0].set(title="Reported mass and quoted 95% interval", ylabel="Mass / prior centre c", yscale="log")
        case_legend = axes[0, 0].legend(fontsize=9, loc="best")
        axes[0, 0].add_artist(case_legend)
        axes[0, 0].legend(handles=[Line2D([], [], color="#475569", label="Quoted", marker="o"),
                                 Line2D([], [], color="#475569", label="Unweighted", linestyle="--")],
                         fontsize=8, loc="lower left")
        axes[0, 1].set(title="Least squared correction cost", ylabel="Squared correction in noise units")
        axes[0, 1].set_yscale("symlog", linthresh=1)
        axes[0, 1].set_ylim(bottom=0)
        axes[0, 1].legend(handles=[Line2D([], [], color="#475569", label="Positive-mass infimum"),
                                 Line2D([], [], color="#475569", linestyle="--", label="Whole projective line")],
                         fontsize=8)
        positive_weights = [row["misfit_weight"] for row in cells if row["misfit_weight"] > 0]
        axes[1, 0].set(title="Per-reading misfit weight", ylabel="g(Q)", yscale="log")
        if positive_weights:
            axes[1, 0].set_ylim(min(positive_weights) * .75, max(positive_weights) * 1.2)
        axes[1, 1].set(title="Readout relative to the law's median", ylabel="Quoted point / quoted median", yscale="log")
        axes[1, 1].axhline(1, color="#94a3b8", linestyle=":", linewidth=1)
        fig.suptitle("Fixed measured vectors; shrinking declared uncertainty", fontsize=15)
        fragments.append(("Noise sweep", "overview", save(fig, "overview")))

        for index, (name, group) in enumerate(groups.items()):
            cell = min(group, key=lambda row: row["noise_scale"])
            fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), layout="constrained")
            lm = arr(cell, "profile_lm") - np.log(centre)
            prior_99 = np.log(np.tan(np.pi * np.array([.005, .995]) / 2))
            density_bounds = list(prior_99)
            for law, color in (("raw", "#2563eb"), ("quoted", "#dc2626")):
                law_lm = arr(cell, f"{law}_lm") - np.log(centre)
                density = arr(cell, f"{law}_density")
                visible = law_lm[density >= np.max(density) * 1e-7]
                density_bounds.extend(visible[[0, -1]])
                density_bounds.extend(np.log(np.array(cell[law]["interval"]) / centre))
                density_bounds.append(np.log(cell[law]["point"] / centre))
                label = "Unweighted" if law == "raw" else "Current ours29"
                axes[0, 0].plot(law_lm, np.ma.masked_less(density, np.max(density) * 1e-7), label=label, color=color)
                axes[0, 0].axvline(np.log(cell[law]["point"] / centre), color=color, linestyle=":", alpha=.8)
                interval_lm = np.log(np.array(cell[law]["interval"]) / centre)
                axes[0, 0].axvspan(*interval_lm, color=color, alpha=.06)
                axes[1, 1].plot(law_lm, arr(cell, f"{law}_A"), label=label, color=color,
                                linestyle="--" if law == "raw" else "-")
            prior = arr(cell, "prior_density")
            axes[0, 0].plot(lm, np.ma.masked_less_equal(prior, 0), color="#64748b", linestyle="--", label="Prior")
            axes[0, 0].set(title="Normalized laws; dotted lines mark readouts", ylabel="Density per unit ln m", yscale="log")
            axes[0, 0].set_xlim(min(density_bounds), max(density_bounds))
            axes[0, 0].legend(fontsize=9)
            axes[0, 1].plot(lm, arr(cell, "correction_cost"), color="#334155", label="Compatible-candidate cost")
            axes[0, 1].plot(lm, arr(cell, "positive_excess_cost"), color="#ea580c", linestyle="--",
                              label="Cost above positive infimum")
            axes[0, 1].axhline(cell["positive_min_cost"], color="#16a34a", linestyle="--",
                              label="Positive infimum: " + boundary(cell))
            axes[0, 1].axhline(cell["projective_min_cost"], color="#9333ea", linestyle=":", label="Projective minimum")
            axes[0, 1].set(title="Absolute correction cost across positive masses", ylabel="Squared correction in noise units")
            axes[0, 1].set_yscale("symlog", linthresh=1)
            axes[0, 1].set_ylim(bottom=.5 * min(cell["positive_min_cost"], cell["projective_min_cost"]))
            axes[0, 1].legend(fontsize=8)
            axes[1, 1].set(title="Candidate-cloud size entering the readout", ylabel="E[|a*| | m]")
            axes[1, 1].legend(fontsize=9)
            for ax in (axes[0, 0], axes[0, 1], axes[1, 1]):
                ax.set_xlabel("ln(m/c)")
                ax.grid(True)

            # F and m*a share force units. Plotting unscaled a on these axes
            # would incorrectly compare vectors with different dimensions.
            ax = axes[1, 0]
            force = np.asarray(cell["force"], dtype=float)
            point = cell["quoted"]["point"]
            scaled_a = point * np.asarray(cell["acceleration"], dtype=float)
            star = force - np.asarray(cell["quoted"]["force_correction_at_point"])
            for vector, label, color in ((force, "Observed F", "#2563eb"),
                                          (scaled_a, "m × observed a", "#ea580c"),
                                          (star, "Closest compatible F* = m a*", "#16a34a")):
                ax.annotate("", xy=vector, xytext=(0, 0),
                            arrowprops={"arrowstyle": "->", "color": color, "lw": 2})
                ax.plot(*vector, "o", color=color, label=label, markersize=4)
            ax.annotate("", xy=force, xytext=scaled_a,
                        arrowprops={"arrowstyle": "->", "color": "#dc2626", "lw": 1.5})
            ax.plot([], [], color="#dc2626", label="Residual r = F − m a")
            ax.plot([star[0], force[0]], [star[1], force[1]], ":", color="#2563eb", alpha=.7)
            ax.plot([star[0], scaled_a[0]], [star[1], scaled_a[1]], ":", color="#ea580c", alpha=.7)
            ax.update_datalim(np.array([[0, 0], force, scaled_a, star]))
            ax.autoscale_view()
            ax.margins(.2)
            ax.set_aspect("equal", adjustable="datalim")
            ax.set(title=f"Force-space geometry at quoted m = {fmt(point)}", xlabel="Force coordinate 1", ylabel="Force coordinate 2")
            ax.grid(True)
            ax.legend(fontsize=8, loc="best")
            angle = "angle undefined (a zero vector)" if cell["angle_degrees"] is None else f'angle {fmt(cell["angle_degrees"])}°'
            fig.suptitle(f'{name}: {angle}; noise scale {fmt(cell["noise_scale"])}', fontsize=14)
            stem = f"case_{index:02d}_{cell['id']}"
            correction_rows = []
            for label in ("raw", "quoted"):
                result = cell[label]
                vector = lambda value: "(" + ", ".join(fmt(v) for v in value) + ")"
                residual = np.asarray(cell["force"]) - result["point"] * np.asarray(cell["acceleration"])
                correction_rows.append("<tr>" + "".join("<td>" + escape(str(value)) + "</td>" for value in
                    (label.capitalize(), vector(residual), vector(result["force_correction_at_point"]),
                     vector(result["acceleration_correction_at_point"]), fmt(result["cost_at_point"]))) + "</tr>")
            captions[stem] = ('<div class="table-scroll"><table><thead><tr><th>At readout</th><th>r = F − m a</th>'
                              '<th>δF (force units)</th><th>δa (acceleration units)</th><th>Correction cost</th></tr></thead>'
                              '<tbody>' + "".join(correction_rows) + '</tbody></table></div>')
            fragments.append((name + " at smallest noise", stem, save(fig, stem)))

    def table(rows):
        headings = ["Case", "Noise scale", "Raw point", "Raw 95% interval", "Quoted point", "Quoted median",
                    "Quoted 95% interval", "Positive cost infimum", "Where", "Q (projective)",
                    "Extra cost from m > 0", "Quoted log-evidence span", "g(Q)", "η"]
        parts = ['<div class="table-scroll"><table><thead><tr>']
        parts += ["<th>" + escape(text) + "</th>" for text in headings]
        parts.append("</tr></thead><tbody>")
        for row in rows:
            quoted = row["quoted"]
            values = [row["case"], fmt(row["noise_scale"]), fmt(row["raw"]["point"]),
                      "[" + ", ".join(fmt(v) for v in row["raw"]["interval"]) + "]",
                      fmt(quoted["point"]), fmt(quoted["median"]),
                      "[" + ", ".join(fmt(v) for v in quoted["interval"]) + "]",
                      fmt(row["positive_min_cost"]), boundary(row), fmt(row["equation_Q"]),
                      fmt(row["positive_domain_penalty"]), fmt(quoted["log_evidence_span"]),
                      fmt(row["misfit_weight"]), fmt(row["bartlett_eta"])]
            parts.append("<tr>" + "".join("<td>" + escape(str(value)) + "</td>" for value in values) + "</tr>")
        parts.append("</tbody></table></div>")
        return "".join(parts)

    endpoint_rows = []
    for group in groups.values():
        endpoint_rows.append(group[0])
        if len(group) > 1:
            endpoint_rows.append(group[-1])
    figures = "".join(
        f'<section><h2>{escape(title)}</h2>{svg}{captions.get(stem, "")}<p class="downloads">'
        f'<a href="{stem}.png">PNG</a> · <a href="{stem}.svg">SVG</a></p></section>'
        for title, stem, svg in fragments)
    warned = [(cell, label) for cell in cells for label in ("raw", "quoted")
              if cell[label].get("support_edge_warning", False)]
    warning = ""
    if warned:
        labels = "; ".join(f'{cell["case"]}, noise {fmt(cell["noise_scale"])}, {label} '
                           f'(edge/peak density {fmt(cell[label]["support_edge_density_fraction"])})'
                           for cell, label in warned)
        warning = ('<section class="warning"><h2>Numerical support warning</h2><p>'
                   'Some computed densities remain above 10<sup>−7</sup> of their peak at a finite '
                   'estimator grid edge. Tail probabilities, intervals and readouts may be truncated. '
                   'These edges are numerical limits, not physical mass boundaries.</p><p>' + escape(labels) + '</p></section>')
    metadata = escape(json.dumps({key: data.get(key) for key in
                                 ("seed", "deterministic", "prior", "code_version")}, indent=2, default=str))
    html = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{escape(data["name"])}</title><style>
*{{box-sizing:border-box}}body{{margin:0;background:#f3f6fa;color:#1e293b;font:16px/1.6 system-ui,sans-serif}}
main{{max-width:1240px;margin:auto;padding:32px 24px 64px}}h1{{font-size:30px;line-height:1.2;margin:0 0 16px}}
h2{{font-size:21px;margin:0 0 16px}}p{{max-width:100ch}}section{{background:white;padding:26px;margin:24px 0;border:1px solid #dce3ed;border-radius:10px}}
.lead{{font-size:18px}}.note{{border-left:4px solid #2563eb;padding-left:16px}}a{{color:#1d4ed8}}
.warning{{border:2px solid #b45309;background:#fffbeb}}
svg{{width:100%;height:auto;display:block}}.downloads{{font-size:14px;margin:10px 0 0}}.table-scroll{{overflow:auto}}
table{{border-collapse:collapse;width:100%;font-size:13px;font-variant-numeric:tabular-nums}}th,td{{padding:9px 10px;text-align:left;border-bottom:1px solid #e2e8f0;white-space:nowrap}}
th{{color:#475569;background:#f8fafc}}pre{{overflow:auto;font-size:12px;background:#f8fafc;padding:14px}}details{{margin-top:18px}}summary{{cursor:pointer;font-weight:600}}
@media print{{body{{background:white}}main{{padding:0}}section{{break-inside:avoid;border:0;padding:12px 0}}}}
</style></head><body><main>
<h1>{escape(data["name"])}</h1><p class="lead">{escape(data["question"])}</p>
<p class="note">A deterministic stress study of fixed measured vectors. No true mass was planted;
these figures make no interval-coverage claim and impose no automatic pass/fail rule on the equation.</p>
{warning}
<section><h2>How to read this investigation</h2>
<p>The measured force and acceleration stay fixed while their declared standard deviations shrink together.
The prior is a fixed half-Cauchy with centre c = {fmt(centre)}. “Unweighted” (the raw export key)
is the ours27 baseline, with both misfit weighting and Bartlett calibration disabled. “Quoted”
is current ours29, with both enabled. The estimator's native raw_interval disables only Bartlett;
that separate quantity is not the unweighted baseline in this comparison. Every displayed mass density is normalized
and expressed per unit ln m. The density plots use a logarithmic vertical axis;
dotted vertical lines mark readouts and pale shaded bands mark each law's 95% interval.</p>
<p>For readability, density curves omit samples below 10<sup>−7</sup> of their own peak.
The horizontal view includes the remaining samples, both readouts and 95% intervals,
and the prior's central 99%. Full sampled supports are retained in the data files.</p>
<p>At each candidate mass, r(m) = F − m a. For independent isotropic channel errors, its covariance is
(σ<sub>F</sub>² + m²σ<sub>a</sub>²)I and the displayed squared correction cost is
|r(m)|²/(σ<sub>F</sub>² + m²σ<sub>a</sub>²). This is the cost of the closest compatible pair under the declared noise geometry.</p>
<p>The positive-mass infimum is evaluated analytically, including limits m → 0⁺ and m → ∞.
A boundary infimum is never presented as a finite-grid optimum. The projective minimum also permits negative slopes;
it is the quantity used by the current equation's misfit weight. Cost plots use a symmetric logarithmic scale,
linear between −1 and 1, so zero is displayed exactly. The g(Q) axis is logarithmic and includes g = 1 exactly;
any nonpositive weights would be masked rather than replaced by an artificial positive number.</p>
<p>Corrections are observed minus closest-compatible-candidate vectors. They are not recovered actual instrument errors.
The geometry panel puts F, m a and F* = m a* in the same force units, at the quoted point.
The conditional size E[|a*| | m] is a candidate-cloud mean, rather than the length of the closest acceleration vector.
Its two curves evaluate the same conditional mean on each law's own mass grid, so overlapping segments coincide.</p>
<p><strong>Absolute correction cost is not a marginal likelihood, Bayes factor or p-value.</strong>
Normalized mass laws compare candidates within the assumed construction; normalization does not supply an absolute compatibility test.</p>
<p>Write K(m) for the correction cost and Q<sub>+</sub> for its infimum over positive mass.
Since H(m) = E − K(m), the normalized law can be written as prior(m) ×
exp[−η g (K(m) − Q<sub>+</sub>)/2], followed by normalization. The constant Q<sub>+</sub>
cancels from this expression. The orange cost curve shows the remaining excess K(m) − Q<sub>+</sub>.
The table separates Q<sub>+</sub> − Q from the projective misfit Q; this extra cost comes from
forbidding negative masses and is not substituted into the estimator's weight.</p>
<p>The log-evidence span is the range of log(law/prior) over positive mass, in nats.
It describes relative preference, independently of a constant cost offset. A stable point or median
need not mean a stable law: reciprocal anti-parallel readings can keep both centres fixed while
probability moves toward zero and infinite mass. The full intervals and density curves must be read with the centres.</p>
<p>In the anti-parallel case with equal whitened channel energies and c = σ<sub>F</sub>/σ<sub>a</sub>, reciprocal symmetry
fixes the mathematical median at c. At small noise the numerical cumulative probability can plateau at one-half,
so its exported median can drift. This is an identified numerical discrepancy, recorded in the tests;
it does not establish a physical asymmetry between the channels.</p>
</section>
{figures}
<section><h2>Endpoints of each noise sweep</h2><p>Values use the supplied mass units; plots scale mass by c.
Numbers here are rounded for reading. JSON/CSV retain full numerical precision.</p>{table(endpoint_rows)}
<details><summary>All study cells</summary>{table(cells)}</details></section>
<section><h2>Reproducibility</h2><p>Generated locally with no remote assets. The HTML embeds every plot;
PNG and SVG copies are stored beside it. The seed is recorded for a consistent study interface;
this deterministic sweep draws no random samples.</p>
<p><a href="run.json">Full-precision run JSON</a> · <a href="summary.csv">Cell CSV</a> · <a href="values.npz">All sampled curves (NumPy)</a> · <a href="preset.yaml">Exact preset</a></p>
<details><summary>Recorded inputs and code version</summary><pre>{metadata}</pre></details></section>
</main></body></html>'''
    (folder / "report.html").write_text(html, encoding="utf-8")
    return True
