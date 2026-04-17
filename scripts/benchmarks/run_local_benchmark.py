#!/usr/bin/env python3

import argparse
import json
import os
import platform
import sys
from pathlib import Path
from time import perf_counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from exp_heatmap.benchmark_utils import (
    build_comparison_regions,
    centered_window,
    infer_available_populations,
    normalize_population_counts,
    read_variant_range,
    resolve_population_selection,
)
from exp_heatmap.interactive import (
    create_comparison_view,
    create_population_focus_view,
    plot_interactive_heatmap,
)
from exp_heatmap.logging import setup_logging
from exp_heatmap.plot import (
    create_plot_input,
    extract_top_regions,
    plot_exp_heatmap,
    summarize_by_superpopulation,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run local ExP Heatmap benchmarks across region widths, population counts, and output modes."
    )
    parser.add_argument("input_dir", help="Directory with TSV files from 'exp_heatmap compute'")
    parser.add_argument(
        "--out-dir",
        default="benchmark_output",
        help="Directory where benchmark outputs and reports will be written.",
    )
    parser.add_argument(
        "--window-sizes",
        nargs="+",
        type=int,
        default=[50000, 250000, 1000000],
        help="Region widths in base pairs to benchmark.",
    )
    parser.add_argument(
        "--population-counts",
        nargs="+",
        type=int,
        help="Population counts to benchmark. Defaults to a small and full-panel run.",
    )
    parser.add_argument(
        "--modes",
        nargs="+",
        default=["load", "regions", "static", "interactive", "focus", "compare", "summary"],
        choices=["load", "regions", "static", "interactive", "focus", "compare", "summary"],
        help="Benchmark modes to execute.",
    )
    parser.add_argument(
        "--rank-scores",
        default="directional",
        choices=["directional", "2-tailed", "ascending", "descending"],
        help="Rank-score mode passed into plotting helpers.",
    )
    parser.add_argument(
        "--max-columns",
        type=int,
        default=3000,
        help="Column budget for static and interactive benchmark renders.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="DPI for static benchmark renders.",
    )
    return parser.parse_args()


def _output_size(path):
    return os.path.getsize(path) if os.path.exists(path) else None


def _record_result(results, **kwargs):
    results.append(kwargs)


def _run_mode(mode, plot_input, population_selection, start, end, artifacts_dir, window_size, pop_count, max_columns, dpi):
    artifact_base = artifacts_dir / f"{mode}_w{window_size}_p{pop_count}"
    artifact_path = None

    if mode == "regions":
        artifact_path = artifact_base.with_suffix(".tsv")
        top_regions = extract_top_regions(
            plot_input,
            n_top=10,
            window_size=max(10000, window_size // 10),
            min_gap=max(5000, window_size // 20),
        )
        top_regions.to_csv(artifact_path, sep="\t", index=False)

    elif mode == "static":
        plot_exp_heatmap(
            plot_input,
            start=start,
            end=end,
            title=f"Benchmark static {window_size}bp / {pop_count} pops",
            output=str(artifact_base),
            populations=population_selection,
            max_columns=max_columns,
            dpi=dpi,
        )
        artifact_path = artifact_base.with_suffix(".png")

    elif mode == "interactive":
        plot_interactive_heatmap(
            plot_input,
            start=start,
            end=end,
            title=f"Benchmark interactive {window_size}bp / {pop_count} pops",
            output=str(artifact_base),
            populations=population_selection,
            max_columns=max_columns,
            show_superpop_annotations=population_selection == "1000Genomes",
        )
        artifact_path = artifact_base.with_suffix(".html")

    elif mode == "focus":
        focus_population = population_selection[0] if population_selection != "1000Genomes" else "ACB"
        create_population_focus_view(
            plot_input,
            focus_population=focus_population,
            start=start,
            end=end,
            title=f"Benchmark focus {focus_population}",
            output=str(artifact_base),
            max_columns=max_columns,
        )
        artifact_path = artifact_base.with_suffix(".html")

    elif mode == "compare":
        region1, region2 = build_comparison_regions(start, end)
        create_comparison_view(
            plot_input,
            region1=region1,
            region2=region2,
            title=f"Benchmark compare {window_size}bp / {pop_count} pops",
            output=str(artifact_base),
            max_columns=max_columns,
        )
        artifact_path = artifact_base.with_suffix(".html")

    elif mode == "summary":
        if population_selection != "1000Genomes":
            raise ValueError("summary benchmarking is only available for the full canonical 1000 Genomes panel")

        summary_df = summarize_by_superpopulation(plot_input, populations="1000Genomes", agg_func="mean")
        summary_path = artifact_base.with_suffix(".tsv")
        summary_df.to_csv(summary_path, sep="\t")
        plot_interactive_heatmap(
            summary_df,
            start=int(summary_df.columns[0]),
            end=int(summary_df.columns[-1]),
            title="Benchmark summary",
            output=str(artifact_base),
            populations="1000Genomes",
            max_columns=max_columns,
            show_superpop_annotations=False,
        )
        artifact_path = artifact_base.with_suffix(".html")

    else:
        raise ValueError(f"Unsupported benchmark mode '{mode}'")

    return artifact_path


def render_runtime_plot(results_df, output_path):
    completed = results_df[results_df["status"] == "ok"].copy()
    if completed.empty:
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    grouped = completed.groupby(["mode", "population_count"])
    for (mode, population_count), group in grouped:
        group = group.sort_values("window_size")
        ax.plot(
            group["window_size"],
            group["seconds"],
            marker="o",
            label=f"{mode} / {population_count} pops",
        )

    ax.set_xlabel("Window size (bp)")
    ax.set_ylabel("Runtime (seconds)")
    ax.set_title("ExP Heatmap local benchmark runtimes")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main():
    args = parse_args()
    setup_logging("benchmark", log_to_file=False, verbose=False)

    input_dir = Path(args.input_dir)
    output_dir = Path(args.out_dir)
    artifacts_dir = output_dir / "artifacts"
    output_dir.mkdir(parents=True, exist_ok=True)
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    range_start, range_end = read_variant_range(str(input_dir))
    available_populations, full_population_mode = infer_available_populations(str(input_dir))
    population_counts = normalize_population_counts(available_populations, args.population_counts)

    results = []

    for window_size in args.window_sizes:
        start, end = centered_window(range_start, range_end, window_size)
        for pop_count in population_counts:
            population_selection = resolve_population_selection(
                available_populations,
                full_population_mode,
                pop_count,
            )

            load_start = perf_counter()
            try:
                plot_input = create_plot_input(
                    str(input_dir),
                    start=start,
                    end=end,
                    populations=population_selection,
                    rank_scores=args.rank_scores,
                )
                load_seconds = perf_counter() - load_start
                _record_result(
                    results,
                    mode="load",
                    window_size=window_size,
                    population_count=pop_count,
                    start=start,
                    end=end,
                    seconds=load_seconds,
                    rows=int(plot_input.shape[0]),
                    columns=int(plot_input.shape[1]),
                    artifact_path=None,
                    artifact_size_bytes=None,
                    status="ok",
                    note="",
                )
            except Exception as exc:
                _record_result(
                    results,
                    mode="load",
                    window_size=window_size,
                    population_count=pop_count,
                    start=start,
                    end=end,
                    seconds=perf_counter() - load_start,
                    rows=0,
                    columns=0,
                    artifact_path=None,
                    artifact_size_bytes=None,
                    status="error",
                    note=str(exc),
                )
                continue

            for mode in args.modes:
                if mode == "load":
                    continue

                started = perf_counter()
                try:
                    artifact_path = _run_mode(
                        mode=mode,
                        plot_input=plot_input,
                        population_selection=population_selection,
                        start=int(plot_input.columns[0]),
                        end=int(plot_input.columns[-1]),
                        artifacts_dir=artifacts_dir,
                        window_size=window_size,
                        pop_count=pop_count,
                        max_columns=args.max_columns,
                        dpi=args.dpi,
                    )
                    _record_result(
                        results,
                        mode=mode,
                        window_size=window_size,
                        population_count=pop_count,
                        start=int(plot_input.columns[0]),
                        end=int(plot_input.columns[-1]),
                        seconds=perf_counter() - started,
                        rows=int(plot_input.shape[0]),
                        columns=int(plot_input.shape[1]),
                        artifact_path=str(artifact_path) if artifact_path else None,
                        artifact_size_bytes=_output_size(artifact_path) if artifact_path else None,
                        status="ok",
                        note="",
                    )
                except Exception as exc:
                    _record_result(
                        results,
                        mode=mode,
                        window_size=window_size,
                        population_count=pop_count,
                        start=int(plot_input.columns[0]),
                        end=int(plot_input.columns[-1]),
                        seconds=perf_counter() - started,
                        rows=int(plot_input.shape[0]),
                        columns=int(plot_input.shape[1]),
                        artifact_path=None,
                        artifact_size_bytes=None,
                        status="error",
                        note=str(exc),
                    )

    results_df = pd.DataFrame(results)
    results_path = output_dir / "benchmark_results.tsv"
    results_df.to_csv(results_path, sep="\t", index=False)

    summary = {
        "input_dir": str(input_dir),
        "range_start": range_start,
        "range_end": range_end,
        "available_populations": list(available_populations),
        "full_population_mode": full_population_mode if isinstance(full_population_mode, str) else list(full_population_mode),
        "window_sizes": args.window_sizes,
        "population_counts": population_counts,
        "modes": args.modes,
        "python_version": sys.version,
        "platform": platform.platform(),
        "results_path": str(results_path),
    }
    summary_path = output_dir / "benchmark_summary.json"
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    runtime_plot_path = output_dir / "benchmark_runtime.png"
    render_runtime_plot(results_df, runtime_plot_path)

    print(f"Wrote benchmark results to {results_path}")
    print(f"Wrote benchmark summary to {summary_path}")
    if runtime_plot_path.exists():
        print(f"Wrote runtime plot to {runtime_plot_path}")


if __name__ == "__main__":
    main()
