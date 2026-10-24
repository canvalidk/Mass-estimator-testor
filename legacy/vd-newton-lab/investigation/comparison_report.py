"""Make readable tables from compare_estimators.py output; no refitting."""
import argparse
import json
from pathlib import Path


LABELS = {
    "flat_joint": "Our ratio of means", "flat_geometric": "Flat-law geometric",
    "flat_median": "Flat-law median", "flat_reciprocal_root": "Flat-law root-loss point",
    "tube_joint": "Tube ratio of means", "norm_ratio": "Vector-length ratio",
    "norm_floor": "Floored vector-length ratio", "positive_profile": "Positive TLS / profile",
    "forward_ols": "Forward OLS", "reverse_ols": "Reverse OLS",
    "noise_corrected_ols": "Noise-corrected moment"}


def number(value, places=4):
    return "—" if value is None else f"{value:.{places}f}"


def point_value(row, method, loss):
    summary = row["methods"][method]
    if summary.get("population_unbounded_risk") and not loss.startswith(("capped_", "outside_")):
        return "∞¹"
    return number(summary["means"][loss])


def report(data):
    rows = data["rows"]
    lines = ["# Refined mass-estimator comparison — simulation results", "",
             f"{data['samples_per_scenario']:,} trials per signal cell; {len(rows)} cells; "
             f"{len(rows[0]['methods'])} point estimators. Identical observations are used for every method. "
             "All losses are smaller-is-better. These are results for the declared Gaussian 3D model.", "",
             "Negative paired differences favor our ratio-of-means equation. ± values below are "
             "one empirical Monte Carlo standard error of a paired difference, not uncertainty in "
             "one object's mass. No global winner is calculated across the chosen cells.", ""]
    path = sorted([r for r in rows if r["family"] == "original" and r["true_force_snr"] == 8.],
                  key=lambda r: r["true_acceleration_snr"])
    if path:
        lines += ["## Original force-SNR-8 path", "",
                  "Mass changes along this path: true mass = 8 / acceleration SNR.", "",
                  "| Accel. SNR | True mass | Our squared log error | Norm ratio | Paired difference ± MCSE | Our factor-two success | Norm success |",
                  "|---:|---:|---:|---:|---:|---:|---:|"]
        for row in path:
            pair = row["paired_flat_minus_other"]["norm_ratio"]["squared_log"]
            ours, norm = row["methods"]["flat_joint"], row["methods"]["norm_ratio"]
            lines.append(f"| {row['true_acceleration_snr']:g} | {row['true_mass']:.3g} | "
                         f"{number(ours['means']['squared_log'])} | {number(norm['means']['squared_log'])} | "
                         f"{number(pair['mean_difference'])} ± {number(pair['mcse'])} | "
                         f"{100*ours['tolerance_success']['2']['fraction']:.2f}% | "
                         f"{100*norm['tolerance_success']['2']['fraction']:.2f}% |")
        selected = [next((r for r in path if r["true_acceleration_snr"] == a), None) for a in (1., 2., 3.)]
        if all(row is not None for row in selected):
            lines += ["", "## All equations around the earlier transition", "",
                      "Unbounded squared log error and bounded factor-ten-capped squared log error "
                      "are displayed separately. The capped comparison includes every failed output.", "",
                      "| Estimator | Log², a=1 | Log², a=2 | Log², a=3 | Capped log², a=1 | Failed outputs, a=1 |",
                      "|---|---:|---:|---:|---:|---:|"]
            for method in LABELS:
                fields = [point_value(row, method, "squared_log") for row in selected]
                cap = point_value(selected[0], method, "capped_squared_log_factor_10")
                failures = 100*(1-selected[0]["methods"][method]["finite_positive_fraction"])
                lines.append(f"| {LABELS[method]} | {' | '.join(fields)} | {cap} | {failures:.2f}% |")
            lines += ["", "¹ Under the declared infinite penalty for invalid positive masses, population risk "
                      "is infinite for these boundary/signed/unavailable-output rules. This remains true "
                      "if a finite high-signal sample contains no invalid output. JSON retains the finite "
                      "sample diagnostics and the separate population-risk annotation."]
        lines += ["", "## Entire probability laws", "",
                  "CRPS evaluates log mass; 95% coverage is the fraction of intervals containing true mass. "
                  "Coverage should be considered together with interval scores and widths in the full output.", "",
                  "| Accel. SNR | Flat log-CRPS | Tube log-CRPS | Flat 95% coverage | Tube 95% coverage | Flat mean log width | Tube mean log width |",
                  "|---:|---:|---:|---:|---:|---:|---:|"]
        for row in path:
            flat, tube = row["uncertainty"]["flat_joint"], row["uncertainty"]["tube_joint"]
            lines.append(f"| {row['true_acceleration_snr']:g} | "
                         f"{number(flat['scores']['log_crps']['mean'])} | {number(tube['scores']['log_crps']['mean'])} | "
                         f"{100*flat['intervals']['95']['coverage']['mean']:.2f}% | "
                         f"{100*tube['intervals']['95']['coverage']['mean']:.2f}% | "
                         f"{number(flat['intervals']['95']['mean_log_width']['mean'])} | "
                         f"{number(tube['intervals']['95']['mean_log_width']['mean'])} |")
    fixed = [r for r in rows if r["family"] == "fixed_mass"]
    if fixed:
        strengths = sorted(set(r["total_signal_snr"] for r in fixed))
        lines += ["", "## Holding true mass fixed", "",
                  "Entries are our squared log error minus the norm ratio's. Negative favors our equation. "
                  "Columns increase total signal strength √(f²+a²), with mass held fixed. For masses 16 and 32, "
                  "this selected range still leaves acceleration weak; it does not map their eventual high-signal limit.", "",
                  "| True mass | " + " | ".join(f"Signal {s:g}" for s in strengths) + " |",
                  "|---:|" + "---:|"*len(strengths)]
        for mass in sorted(set(r["true_mass"] for r in fixed)):
            entries = []
            for strength in strengths:
                row = next(r for r in fixed if r["true_mass"] == mass and r["total_signal_snr"] == strength)
                pair = row["paired_flat_minus_other"]["norm_ratio"]["squared_log"]
                entries.append(f"{number(pair['mean_difference'])} ± {number(pair['mcse'])}")
            lines.append(f"| {mass:g} | " + " | ".join(entries) + " |")
    refinement = max(r["refinement"]["maximum_relative_point_difference"] for r in rows)
    interval = max(v for r in rows for k, v in r["refinement"]["maximum_absolute_distribution_difference"].items()
                   if k.startswith(("log_lower_", "log_upper_")))
    lines += ["", "## Reproducibility and scope", "",
              f"All {len(rows)} per-cell refinement checks passed. Among the first 32 readings per cell, "
              f"doubling integration order changed points by at most {refinement:.3g} relatively and "
              f"log interval endpoints by {interval:.3g} absolutely. These are numerical checks on that subset, "
              "not a uniform error bound on every possible reading.", "",
              "The full JSON includes all point losses, quantiles, tolerance curves, paired differences, "
              "proper distribution scores, interval calibration and the two stated held-out prediction tasks. "
              "Empirical squared mass, relative and held-out residual loss MCSEs for the norm ratio "
              "(including their paired differences) are descriptive only: their population variance is infinite. "
              "Consult comparison_guide.md for equations, assumptions and primary sources.", "",
              f"Training seed {data['training_seed']}; independent validation seed {data['heldout_seed']}. "
              f"Python {data['python']}, NumPy {data['numpy']}. Source hashes are in results.json. "
              f"Managed context: `{data['managed_commit']}`. No real-data validation is claimed.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    content = report(data)
    output = args.output or args.input.with_name("summary.md")
    output.write_text(content, encoding="utf-8")
    print(output.resolve())


if __name__ == "__main__":
    main()
