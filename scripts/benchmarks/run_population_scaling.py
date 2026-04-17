#!/usr/bin/env python3

import argparse
import itertools
import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import psutil

from exp_heatmap.benchmark_utils import read_variant_range, centered_window


def dataframe_to_markdown(df):
    headers = [str(column) for column in df.columns]
    rows = [[str(value) for value in row] for row in df.itertuples(index=False, name=None)]
    widths = [len(header) for header in headers]

    for row in rows:
        for index, value in enumerate(row):
            widths[index] = max(widths[index], len(value))

    def format_row(values):
        return "| " + " | ".join(
            value.ljust(widths[index]) for index, value in enumerate(values)
        ) + " |"

    separator = "| " + " | ".join("-" * widths[index] for index in range(len(widths))) + " |"
    lines = [format_row(headers), separator]
    lines.extend(format_row(row) for row in rows)
    return "\n".join(lines)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run a compact ordered-pair display scaling benchmark by synthesizing larger "
            "population panels from an existing compute directory."
        )
    )
    parser.add_argument(
        "template_dir",
        help="Directory with benchmark template TSV files from 'exp_heatmap compute'.",
    )
    parser.add_argument(
        "--out-dir",
        default="local_data/benchmarks/population_scaling",
        help="Directory where synthetic inputs, logs, and summary files will be written.",
    )
    parser.add_argument(
        "--population-counts",
        nargs="+",
        type=int,
        default=[5, 10, 20, 30, 50],
        help="Population counts to benchmark.",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=1_000_000,
        help="Fixed genomic window width used for all scaling measurements.",
    )
    parser.add_argument(
        "--rank-scores",
        default="directional",
        choices=["directional", "2-tailed", "ascending", "descending"],
        help="Rank-score mode used for display scaling.",
    )
    parser.add_argument(
        "--max-columns",
        type=int,
        default=3000,
        help="Column budget for static and interactive outputs.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="DPI for static outputs.",
    )
    parser.add_argument(
        "--worker-stage",
        choices=["load", "static", "interactive"],
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--input-dir", help=argparse.SUPPRESS)
    parser.add_argument("--start", type=int, help=argparse.SUPPRESS)
    parser.add_argument("--end", type=int, help=argparse.SUPPRESS)
    parser.add_argument("--out-prefix", help=argparse.SUPPRESS)
    parser.add_argument("--populations-json", help=argparse.SUPPRESS)
    return parser.parse_args()


def _rss_for_process_tree(process):
    try:
        processes = [process] + process.children(recursive=True)
    except (psutil.Error, PermissionError):
        processes = [process]

    total = 0
    for proc in processes:
        try:
            total += proc.memory_info().rss
        except psutil.Error:
            continue
    return total


def run_monitored_command(cmd, cwd, env, log_path):
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write(f"$ {shlex.join(cmd)}\n\n")
        log_handle.flush()

        start_time = time.perf_counter()
        process = subprocess.Popen(
            cmd,
            cwd=cwd,
            env=env,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
        ps_process = psutil.Process(process.pid)
        peak_rss = 0

        while True:
            peak_rss = max(peak_rss, _rss_for_process_tree(ps_process))
            if process.poll() is not None:
                break
            time.sleep(0.05)

        try:
            peak_rss = max(peak_rss, _rss_for_process_tree(ps_process))
        except psutil.Error:
            pass

        seconds = time.perf_counter() - start_time

    return {
        "returncode": int(process.returncode),
        "seconds": seconds,
        "peak_rss_bytes": int(peak_rss),
    }


def create_synthetic_compute_dir(template_dir, output_dir, population_count):
    output_dir.mkdir(parents=True, exist_ok=True)
    template_files = sorted(Path(template_dir).glob("*.tsv"))
    if not template_files:
        raise ValueError(f"No TSV files found in template directory '{template_dir}'")

    population_names = tuple(f"P{index:02d}" for index in range(1, population_count + 1))

    pair_iter = itertools.combinations(population_names, 2)
    for index, (pop1, pop2) in enumerate(pair_iter):
        template = template_files[index % len(template_files)]
        target = output_dir / f"{pop1}_{pop2}.tsv"
        if target.exists():
            continue
        os.symlink(template.resolve(), target)

    return population_names


def run_worker(args):
    from exp_heatmap.interactive import plot_interactive_heatmap
    from exp_heatmap.plot import create_plot_input, plot_exp_heatmap

    populations = tuple(json.loads(args.populations_json))
    plot_input = create_plot_input(
        args.input_dir,
        start=args.start,
        end=args.end,
        populations=populations,
        rank_scores=args.rank_scores,
    )

    if args.worker_stage == "load":
        stats = {
            "rows": int(plot_input.shape[0]),
            "columns": int(plot_input.shape[1]),
        }
        Path(f"{args.out_prefix}.json").write_text(json.dumps(stats), encoding="utf-8")
        return

    if args.worker_stage == "static":
        plot_exp_heatmap(
            plot_input,
            start=int(plot_input.columns[0]),
            end=int(plot_input.columns[-1]),
            title=f"Synthetic scaling static ({len(populations)} populations)",
            output=args.out_prefix,
            populations=populations,
            max_columns=args.max_columns,
            dpi=args.dpi,
        )
        return

    if args.worker_stage == "interactive":
        plot_interactive_heatmap(
            plot_input,
            start=int(plot_input.columns[0]),
            end=int(plot_input.columns[-1]),
            title=f"Synthetic scaling interactive ({len(populations)} populations)",
            output=args.out_prefix,
            populations=populations,
            max_columns=args.max_columns,
            show_superpop_annotations=False,
        )
        return

    raise ValueError(f"Unsupported worker stage '{args.worker_stage}'")


def main():
    args = parse_args()
    if args.worker_stage:
        run_worker(args)
        return

    template_dir = Path(args.template_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    synthetic_root = out_dir / "synthetic_inputs"
    logs_dir = out_dir / "logs"
    artifacts_dir = out_dir / "artifacts"
    synthetic_root.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(exist_ok=True)
    artifacts_dir.mkdir(exist_ok=True)

    range_start, range_end = read_variant_range(str(template_dir))
    start, end = centered_window(range_start, range_end, args.window_size)

    env = os.environ.copy()
    env["PYTHONPATH"] = "src"
    env.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

    results = []
    repo_root = Path.cwd().resolve()

    for population_count in args.population_counts:
        synthetic_dir = synthetic_root / f"pops_{population_count}"
        populations = create_synthetic_compute_dir(
            template_dir,
            synthetic_dir,
            population_count,
        )

        for stage in ("load", "static", "interactive"):
            out_prefix = artifacts_dir / f"{stage}_p{population_count}"
            log_path = logs_dir / f"{stage}_p{population_count}.log"
            command = [
                sys.executable,
                str(Path(__file__).resolve()),
                "--worker-stage",
                stage,
                "--input-dir",
                str(synthetic_dir),
                "--start",
                str(start),
                "--end",
                str(end),
                "--out-prefix",
                str(out_prefix),
                "--populations-json",
                json.dumps(populations),
                "--rank-scores",
                args.rank_scores,
                "--max-columns",
                str(args.max_columns),
                "--dpi",
                str(args.dpi),
                str(template_dir),
            ]

            stats = run_monitored_command(command, cwd=repo_root, env=env, log_path=log_path)
            if stats["returncode"] != 0:
                raise RuntimeError(
                    f"Population scaling stage '{stage}' failed for population count {population_count}. "
                    f"See log: {log_path}"
                )

            rows = None
            columns = None
            artifact_path = None
            artifact_size = None

            if stage == "load":
                load_json = Path(f"{out_prefix}.json")
                payload = json.loads(load_json.read_text(encoding="utf-8"))
                rows = int(payload["rows"])
                columns = int(payload["columns"])
                artifact_path = str(load_json)
                artifact_size = load_json.stat().st_size
            elif stage == "static":
                png_path = Path(f"{out_prefix}.png")
                artifact_path = str(png_path)
                artifact_size = png_path.stat().st_size
            else:
                html_path = Path(f"{out_prefix}.html")
                artifact_path = str(html_path)
                artifact_size = html_path.stat().st_size

            results.append(
                {
                    "population_count": int(population_count),
                    "ordered_pair_rows": int(population_count * (population_count - 1)),
                    "window_start": int(start),
                    "window_end": int(end),
                    "stage": stage,
                    "rows": rows,
                    "columns": columns,
                    "seconds": round(stats["seconds"], 6),
                    "peak_rss_bytes": int(stats["peak_rss_bytes"]),
                    "peak_rss_gb": round(stats["peak_rss_bytes"] / (1024 ** 3), 4),
                    "artifact_path": artifact_path,
                    "artifact_size_bytes": int(artifact_size),
                    "artifact_size_mb": round(artifact_size / (1024 ** 2), 4),
                    "log_path": str(log_path),
                    "command": shlex.join(command),
                }
            )

    results_df = pd.DataFrame(results)
    results_tsv = out_dir / "population_scaling_results.tsv"
    results_df.to_csv(results_tsv, sep="\t", index=False)

    summary = results_df.pivot(
        index=["population_count", "ordered_pair_rows"],
        columns="stage",
        values=["seconds", "peak_rss_gb", "artifact_size_mb", "rows", "columns"],
    ).reset_index()
    summary.columns = [
        "_".join([str(part) for part in column if part]).rstrip("_")
        for column in summary.columns.to_flat_index()
    ]
    summary_tsv = out_dir / "population_scaling_summary.tsv"
    summary.to_csv(summary_tsv, sep="\t", index=False)

    note_lines = [
        "# Population scaling summary",
        "",
        f"- Template compute directory: `{template_dir}`",
        f"- Fixed window: `{start}-{end}`",
        f"- Rank-score mode: `{args.rank_scores}`",
        f"- Output directory: `{out_dir}`",
        "",
        dataframe_to_markdown(summary),
        "",
        "## Interpretation",
        "",
        "This synthetic display-scaling benchmark increases the number of population labels while keeping the score profiles and genomic window structure fixed. It therefore benchmarks the ordered-pair display and reporting layer rather than genotype-level statistic computation. The ordered-pair row count grows quadratically with the number of populations, and the outputs remain computationally tractable in this local test, but the display becomes progressively denser as population count rises. That is precisely the regime where focused views, grouping, or summary layers become more important for interpretation.",
    ]
    (out_dir / "population_scaling_summary.md").write_text(
        "\n".join(note_lines) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
