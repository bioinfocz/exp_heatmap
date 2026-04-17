#!/usr/bin/env python3

import argparse
import json
import os
import platform
import shlex
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import psutil
import pysam


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
            "Run a controlled full-pipeline benchmark covering prepare, compute, and plot "
            "for a smaller and a larger input derived from the same dataset."
        )
    )
    parser.add_argument("input_vcf", help="Filtered input VCF used as the large benchmark case.")
    parser.add_argument("panel_file", help="Panel file matching the VCF sample order.")
    parser.add_argument(
        "--out-dir",
        default="local_data/benchmarks/pipeline_benchmark",
        help="Output directory for benchmark results, logs, and generated artifacts.",
    )
    parser.add_argument(
        "--small-start",
        type=int,
        help="Start position for the smaller subset benchmark case. Defaults to a centered 1 Mb window.",
    )
    parser.add_argument(
        "--small-end",
        type=int,
        help="End position for the smaller subset benchmark case. Defaults to a centered 1 Mb window.",
    )
    parser.add_argument(
        "--test",
        default="xpehh",
        choices=["xpehh", "xpnsl", "delta_tajima_d", "hudson_fst"],
        help="Statistic to compute during benchmarking.",
    )
    parser.add_argument(
        "--rank-scores",
        default="directional",
        choices=["directional", "2-tailed", "ascending", "descending"],
        help="Rank-score mode used for the plot stage.",
    )
    parser.add_argument(
        "--max-columns",
        type=int,
        default=3000,
        help="Column budget for the plot stage.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="DPI for the static plot stage.",
    )
    return parser.parse_args()


def scan_vcf_range_and_count(vcf_path):
    first = None
    last = None
    count = 0
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf:
            pos = int(record.pos)
            if first is None:
                first = pos
            last = pos
            count += 1
    if first is None or last is None:
        raise ValueError(f"No records found in VCF '{vcf_path}'")
    return first, last, count


def centered_window(range_start, range_end, width):
    width = int(width)
    available = int(range_end) - int(range_start)
    if width >= available:
        return int(range_start), int(range_end)
    center = (int(range_start) + int(range_end)) // 2
    start = center - width // 2
    end = start + width
    if start < range_start:
        start = int(range_start)
        end = start + width
    if end > range_end:
        end = int(range_end)
        start = end - width
    return int(start), int(end)


def write_subset_vcf(input_vcf, output_vcf, start, end):
    kept = 0
    first = None
    last = None
    with pysam.VariantFile(input_vcf) as src, pysam.VariantFile(output_vcf, "w", header=src.header) as dst:
        for record in src:
            pos = int(record.pos)
            if start <= pos <= end:
                dst.write(record)
                kept += 1
                if first is None:
                    first = pos
                last = pos
    if kept == 0:
        raise ValueError(
            f"Subset region {start}-{end} produced no records from '{input_vcf}'"
        )
    return {
        "records": kept,
        "start": int(first),
        "end": int(last),
    }


def directory_size_bytes(path):
    total = 0
    path = Path(path)
    if not path.exists():
        return None
    if path.is_file():
        return path.stat().st_size
    for entry in path.rglob("*"):
        if entry.is_file():
            total += entry.stat().st_size
    return total


def _rss_for_process_tree(process):
    try:
        processes = [process] + process.children(recursive=True)
    except psutil.Error:
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


def _sysctl_value(name):
    try:
        result = subprocess.run(
            ["sysctl", "-n", name],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()
    except Exception:
        return None


def machine_specs():
    total_memory = psutil.virtual_memory().total
    return {
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "machine": platform.machine(),
        "processor": _sysctl_value("machdep.cpu.brand_string") or platform.processor(),
        "cpu_count_logical": os.cpu_count(),
        "total_memory_bytes": int(total_memory),
        "total_memory_gb": round(total_memory / (1024 ** 3), 2),
    }


def main():
    args = parse_args()

    input_vcf = Path(args.input_vcf).resolve()
    panel_file = Path(args.panel_file).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    run_label = time.strftime("%Y%m%d_%H%M%S")
    cases_dir = out_dir / "cases"
    logs_dir = out_dir / "logs" / run_label
    runs_dir = out_dir / "runs" / run_label
    cases_dir.mkdir(exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    runs_dir.mkdir(parents=True, exist_ok=True)

    range_start, range_end, full_records = scan_vcf_range_and_count(input_vcf)

    if args.small_start is not None and args.small_end is not None:
        small_start = int(args.small_start)
        small_end = int(args.small_end)
    elif args.small_start is None and args.small_end is None:
        small_start, small_end = centered_window(range_start, range_end, 1_000_000)
    else:
        raise ValueError("Provide both --small-start and --small-end, or neither.")

    small_vcf = cases_dir / "ggvp_chr21_small_1mb.vcf"
    small_case = write_subset_vcf(input_vcf, small_vcf, small_start, small_end)

    cases = [
        {
            "case": "small_1mb_region",
            "input_vcf": small_vcf,
            "input_records": int(small_case["records"]),
            "input_start": int(small_case["start"]),
            "input_end": int(small_case["end"]),
            "plot_start": int(small_case["start"]),
            "plot_end": int(small_case["end"]),
        },
        {
            "case": "full_chr21",
            "input_vcf": input_vcf,
            "input_records": int(full_records),
            "input_start": int(range_start),
            "input_end": int(range_end),
            "plot_start": int(range_start),
            "plot_end": int(range_end),
        },
    ]

    benchmark_env = os.environ.copy()
    benchmark_env["PYTHONPATH"] = "src"
    benchmark_env.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

    results = []
    workspace = input_vcf.parent.parent
    repo_root = Path.cwd().resolve()

    for case in cases:
        case_name = case["case"]
        case_root = runs_dir / case_name
        case_root.mkdir(parents=True, exist_ok=True)

        zarr_dir = case_root / "zarr"
        compute_dir = case_root / "compute"
        plot_prefix = case_root / "plot"

        commands = {
            "prepare": [
                sys.executable,
                "-m",
                "exp_heatmap.cli",
                "prepare",
                str(case["input_vcf"]),
                "-o",
                str(zarr_dir),
                "--no-log",
            ],
            "compute": [
                sys.executable,
                "-m",
                "exp_heatmap.cli",
                "compute",
                str(zarr_dir),
                str(panel_file),
                "-o",
                str(compute_dir),
                "-t",
                args.test,
                "-c",
                "--no-log",
            ],
            "plot": [
                sys.executable,
                "-m",
                "exp_heatmap.cli",
                "plot",
                str(compute_dir),
                "-s",
                str(case["plot_start"]),
                "-e",
                str(case["plot_end"]),
                "-o",
                str(plot_prefix),
                "--rank-scores",
                args.rank_scores,
                "--max-columns",
                str(args.max_columns),
                "--dpi",
                str(args.dpi),
                "--no-log",
            ],
        }

        for stage_name in ("prepare", "compute", "plot"):
            log_path = logs_dir / f"{case_name}_{stage_name}.log"
            command = commands[stage_name]
            stats = run_monitored_command(command, cwd=repo_root, env=benchmark_env, log_path=log_path)

            if stage_name == "prepare":
                artifact_path = zarr_dir
            elif stage_name == "compute":
                artifact_path = compute_dir
            else:
                artifact_path = plot_prefix.with_suffix(".png")

            artifact_size = directory_size_bytes(artifact_path)
            results.append(
                {
                    "case": case_name,
                    "stage": stage_name,
                    "input_vcf": str(case["input_vcf"]),
                    "input_records": int(case["input_records"]),
                    "input_start": int(case["input_start"]),
                    "input_end": int(case["input_end"]),
                    "plot_start": int(case["plot_start"]),
                    "plot_end": int(case["plot_end"]),
                    "seconds": round(stats["seconds"], 6),
                    "peak_rss_bytes": int(stats["peak_rss_bytes"]),
                    "peak_rss_gb": round(stats["peak_rss_bytes"] / (1024 ** 3), 4),
                    "returncode": int(stats["returncode"]),
                    "artifact_path": str(artifact_path),
                    "artifact_size_bytes": artifact_size,
                    "artifact_size_mb": round((artifact_size or 0) / (1024 ** 2), 4),
                    "command": shlex.join(command),
                    "log_path": str(log_path),
                }
            )

            if stats["returncode"] != 0:
                raise RuntimeError(
                    f"Benchmark stage '{stage_name}' failed for case '{case_name}'. "
                    f"See log: {log_path}"
                )

    results_df = pd.DataFrame(results)
    results_tsv = out_dir / "pipeline_stage_results.tsv"
    results_df.to_csv(results_tsv, sep="\t", index=False)

    machine_json = out_dir / "machine_specs.json"
    machine_json.write_text(json.dumps(machine_specs(), indent=2), encoding="utf-8")

    serialized_cases = []
    for case in cases:
        serialized_cases.append(
            {
                key: str(value) if isinstance(value, Path) else value
                for key, value in case.items()
            }
        )

    metadata = {
        "input_vcf": str(input_vcf),
        "panel_file": str(panel_file),
        "test": args.test,
        "rank_scores": args.rank_scores,
        "max_columns": int(args.max_columns),
        "dpi": int(args.dpi),
        "run_label": run_label,
        "cases": serialized_cases,
        "workspace": str(workspace),
        "repo_root": str(repo_root),
    }
    (out_dir / "benchmark_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    summary = (
        results_df.pivot(index="case", columns="stage", values=["seconds", "peak_rss_gb", "artifact_size_mb"])
        .sort_index(axis=1)
        .reset_index()
    )
    summary.columns = [
        "case" if column == ("case", "") else f"{column[1]}_{column[0]}"
        for column in summary.columns.to_flat_index()
    ]
    summary_tsv = out_dir / "pipeline_stage_summary.tsv"
    summary.to_csv(summary_tsv, sep="\t", index=False)

    specs = machine_specs()
    lines = [
        "# Full-pipeline benchmark summary",
        "",
        f"- Input VCF: `{input_vcf}`",
        f"- Panel file: `{panel_file}`",
        f"- Test: `{args.test}`",
        f"- Rank-score mode for plotting: `{args.rank_scores}`",
        f"- Output directory: `{out_dir}`",
        "",
        dataframe_to_markdown(summary),
        "",
        "## Machine context",
        "",
        f"- Platform: `{specs['platform']}`",
        f"- Processor: `{specs['processor']}`",
        f"- Logical CPU count: `{specs['cpu_count_logical']}`",
        f"- Total RAM (GB): `{specs['total_memory_gb']}`",
        "",
        "## Notes",
        "",
        "- The smaller case uses a 1 Mb subset derived from the same filtered GGVP chromosome 21 VCF.",
        "- The larger case uses the full filtered GGVP chromosome 21 VCF.",
        "- Peak RAM is the observed maximum resident set size across the monitored process tree.",
        "- Raw command logs are stored in the `logs/` subdirectory.",
    ]
    (out_dir / "pipeline_stage_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
