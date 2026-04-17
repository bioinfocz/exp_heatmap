#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import pandas as pd


MODE_ORDER = ["load", "regions", "static", "interactive", "focus", "compare"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create manuscript-ready benchmark assets from a benchmark_results.tsv file."
    )
    parser.add_argument("results_path", type=Path, help="Path to benchmark_results.tsv")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("manuscript_assets"),
        help="Directory for the generated figure and summary tables",
    )
    return parser.parse_args()


def format_window(window_size: int) -> str:
    if window_size >= 1_000_000:
        return f"{window_size // 1_000_000} Mb"
    if window_size >= 1_000:
        return f"{window_size // 1_000} kb"
    return str(window_size)


def format_seconds(value: float | None) -> str:
    if value is None or pd.isna(value):
        return "-"
    if value < 0.01:
        return f"{value:.3f}"
    if value < 0.1:
        return f"{value:.2f}"
    return f"{value:.2f}"


def format_size(value: float | None) -> str:
    if value is None or pd.isna(value):
        return "-"
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f} MB"
    if value >= 1_000:
        return f"{value / 1_000:.1f} KB"
    return f"{int(round(value))} B"


def build_scenario_table(results: pd.DataFrame) -> tuple[pd.DataFrame, list[tuple[int, int]]]:
    scenarios = sorted(
        {(int(row.window_size), int(row.population_count)) for row in results.itertuples()},
        key=lambda item: (item[0], item[1]),
    )
    scenario_labels = [f"{format_window(window)} / {pop_count} pops" for window, pop_count in scenarios]

    runtime_matrix = pd.DataFrame(index=MODE_ORDER, columns=scenario_labels, dtype=float)
    size_matrix = pd.DataFrame(index=MODE_ORDER, columns=scenario_labels, dtype=float)

    for window_size, population_count in scenarios:
        label = f"{format_window(window_size)} / {population_count} pops"
        subset = results[
            (results["window_size"] == window_size)
            & (results["population_count"] == population_count)
        ].set_index("mode")
        for mode in MODE_ORDER:
            if mode not in subset.index:
                continue
            runtime_matrix.loc[mode, label] = subset.loc[mode, "seconds"]
            size_matrix.loc[mode, label] = subset.loc[mode, "artifact_size_bytes"]

    return runtime_matrix, size_matrix


def draw_heatmap(ax: plt.Axes, data: pd.DataFrame, title: str, formatter, *, log_scale: bool) -> None:
    matrix = data.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(matrix)
    cmap = plt.cm.YlGnBu.copy()
    cmap.set_bad("#f1f3f5")

    finite_values = matrix[np.isfinite(matrix)]
    if finite_values.size == 0:
        raise ValueError(f"No finite values available for {title}")

    if log_scale:
        positive_values = finite_values[finite_values > 0]
        norm = colors.LogNorm(vmin=positive_values.min(), vmax=positive_values.max())
    else:
        norm = colors.Normalize(vmin=finite_values.min(), vmax=finite_values.max())

    image = ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto")
    ax.set_title(title, fontsize=12, weight="bold")
    ax.set_xticks(np.arange(data.shape[1]), labels=data.columns, rotation=20, ha="right")
    ax.set_yticks(np.arange(data.shape[0]), labels=[label.title() for label in data.index])
    ax.tick_params(axis="both", labelsize=9)

    for row_idx in range(data.shape[0]):
        for col_idx in range(data.shape[1]):
            value = matrix[row_idx, col_idx]
            text = formatter(value)
            if pd.isna(value):
                text_color = "#52606d"
            else:
                text_color = "white" if norm(value) > 0.55 else "#102a43"
            ax.text(
                col_idx,
                row_idx,
                text,
                ha="center",
                va="center",
                fontsize=8,
                color=text_color,
            )

    colorbar = plt.colorbar(image, ax=ax, fraction=0.045, pad=0.03)
    colorbar.ax.tick_params(labelsize=8)


def dataframe_to_markdown(dataframe: pd.DataFrame) -> str:
    headers = list(dataframe.columns)
    divider = ["---"] * len(headers)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(divider) + " |",
    ]
    for row in dataframe.itertuples(index=False):
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def make_summary_tables(results: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    scenario_rows = []
    for (window_size, population_count), subset in results.groupby(
        ["window_size", "population_count"], sort=True
    ):
        by_mode = subset.set_index("mode")
        scenario_rows.append(
            {
                "window": format_window(int(window_size)),
                "population_count": int(population_count),
                "variant_columns": int(by_mode["columns"].max()),
                "matrix_rows": int(by_mode["rows"].max()),
                "load_s": by_mode.at["load", "seconds"],
                "regions_s": by_mode.at["regions", "seconds"],
                "static_s": by_mode.at["static", "seconds"],
                "interactive_s": by_mode.at["interactive", "seconds"],
                "focus_s": by_mode.at["focus", "seconds"],
                "compare_s": by_mode.at["compare", "seconds"],
                "static_png": format_size(by_mode.at["static", "artifact_size_bytes"]),
                "interactive_html": format_size(by_mode.at["interactive", "artifact_size_bytes"]),
                "focus_html": format_size(by_mode.at["focus", "artifact_size_bytes"]),
                "compare_html": format_size(by_mode.at["compare", "artifact_size_bytes"]),
            }
        )

    summary = pd.DataFrame(scenario_rows)
    numeric_cols = ["load_s", "regions_s", "static_s", "interactive_s", "focus_s", "compare_s"]
    summary[numeric_cols] = summary[numeric_cols].round(3)

    markdown_summary = dataframe_to_markdown(summary)
    return summary, markdown_summary


def main() -> None:
    args = parse_args()
    results_path = args.results_path.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    results = pd.read_csv(results_path, sep="\t")
    failed_rows = results[results["status"] != "ok"]
    if not failed_rows.empty:
        raise ValueError("Benchmark results contain non-ok rows; aborting asset generation.")

    runtime_matrix, size_matrix = build_scenario_table(results)

    figure_path = out_dir / "ggvp_chr21_benchmark_overview.png"
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)
    draw_heatmap(
        axes[0],
        runtime_matrix,
        "Runtime per mode (seconds)",
        format_seconds,
        log_scale=True,
    )
    draw_heatmap(
        axes[1],
        size_matrix,
        "Artifact size per mode",
        format_size,
        log_scale=True,
    )
    fig.suptitle(
        "ExP Heatmap local GGVP chr21 benchmark\n50 kb and 1 Mb windows, 3- and 5-population subsets",
        fontsize=13,
        weight="bold",
    )
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    summary_table, markdown_summary = make_summary_tables(results)
    tsv_path = out_dir / "ggvp_chr21_benchmark_summary.tsv"
    md_path = out_dir / "ggvp_chr21_benchmark_summary.md"
    summary_table.to_csv(tsv_path, sep="\t", index=False)

    min_runtime = results["seconds"].min()
    max_runtime = results["seconds"].max()
    html_sizes = results[results["mode"].isin(["interactive", "focus", "compare"])]["artifact_size_bytes"]
    html_min_mb = html_sizes.min() / 1_000_000
    html_max_mb = html_sizes.max() / 1_000_000

    md_path.write_text(
        "\n".join(
            [
                "# GGVP chr21 benchmark summary",
                "",
                f"- Runtime range across measured modes: {format_seconds(min_runtime)}-{format_seconds(max_runtime)} s",
                f"- Interactive HTML artifact-size range: {html_min_mb:.2f}-{html_max_mb:.2f} MB",
                "- Source benchmark: local GGVP chr21 smoke-test compute output",
                "",
                markdown_summary,
                "",
            ]
        ),
        encoding="utf-8",
    )

    print(f"Wrote {figure_path}")
    print(f"Wrote {tsv_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
