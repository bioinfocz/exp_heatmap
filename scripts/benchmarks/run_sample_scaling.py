#!/usr/bin/env python3
"""
Measure how compute cost scales with the number of samples.

The other benchmark scripts vary the genomic window and the number of population labels.
Both exercise the display layer, and the population-scaling script does so on synthetic
panels. This one varies sample count on real genotypes and runs the statistic itself, so it
measures the stage that dominates whole-chromosome runtime.

Method: hold the variant set and the population set fixed, subsample each population down
to a target number of samples per population, and run `compute` at each size. Sampling is
deterministic for a given seed so a run can be repeated exactly.

compute requires the panel to list exactly the samples in the Zarr store, in the same
order, so a smaller panel alone is rejected. Each sample size therefore gets its own Zarr
store, written by slicing calldata/GT along the sample axis and copying every variants/*
array unchanged. Because variants/AF is copied rather than recomputed, the allele-frequency
filter keeps the same loci at every sample size, which is what makes the sizes comparable.
Building those stores is untimed setup; only compute is measured.
"""

import argparse
import json
import shlex
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import psutil
import zarr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bench_stats import (
    dataframe_to_markdown,
    machine_specs,
    repeat_measurements,
    summarize_runs,
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("zarr_dir", help="Zarr store from 'exp_heatmap prepare'")
    parser.add_argument("panel_file", help="Panel file matching the Zarr sample order")
    parser.add_argument(
        "--out-dir",
        default="local_data/benchmarks/sample_scaling",
        help="Directory for panels, logs, compute outputs and summaries.",
    )
    parser.add_argument(
        "--samples-per-population",
        nargs="+",
        type=int,
        default=[10, 25, 50, 100],
        help="Sample counts per population to benchmark. Counts larger than the smallest "
             "available population are skipped and reported.",
    )
    parser.add_argument(
        "--populations",
        nargs="+",
        help="Population codes to include. Defaults to every population in the panel.",
    )
    parser.add_argument("--test", default="xpehh",
                        choices=["xpehh", "xpnsl", "delta_tajima_d", "hudson_fst"])
    parser.add_argument("--chunked", action="store_true",
                        help="Pass -c to compute, for lower memory on large inputs.")
    parser.add_argument("--repeats", type=int, default=1,
                        help="Measured compute runs per sample size.")
    parser.add_argument("--warmup", type=int, default=0,
                        help="Unmeasured runs executed first and discarded.")
    parser.add_argument("--seed", type=int, default=0,
                        help="Seed for deterministic subsampling.")
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


def run_monitored_command(cmd, log_path):
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write(f"$ {shlex.join(cmd)}\n\n")
        log_handle.flush()
        start_time = time.perf_counter()
        process = subprocess.Popen(cmd, stdout=log_handle, stderr=subprocess.STDOUT, text=True)
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
    return {"returncode": int(process.returncode), "seconds": seconds,
            "peak_rss_bytes": int(peak_rss)}


def write_sample_subset_zarr(source_dir, target_dir, sample_indices):
    """
    Copy a Zarr callset keeping only the given samples.

    variants/* arrays are copied unchanged so the allele-frequency filter selects the same
    loci at every sample size; only calldata/GT is sliced along the sample axis.
    """
    if target_dir.exists():
        return target_dir

    import numcodecs

    def create(group, name, data, **kwargs):
        """Create an array, supplying a codec when the source array holds Python objects."""
        data = np.asarray(data)
        if data.dtype == object:
            kwargs["object_codec"] = numcodecs.VLenUTF8()
        return group.create_dataset(name, data=data, **kwargs)

    source = zarr.open_group(str(source_dir), mode="r")
    staging = target_dir.with_name(target_dir.name + ".partial")
    if staging.exists():
        import shutil

        shutil.rmtree(staging)
    target = zarr.open_group(str(staging), mode="w")

    indices = np.asarray(sample_indices, dtype=int)
    create(target, "samples", np.asarray(source["samples"][:])[indices])

    variants = source["variants"]
    variant_group = target.create_group("variants")
    for name in variants.array_keys():
        create(variant_group, name, variants[name][:])

    genotypes = source["calldata/GT"]
    calldata = target.create_group("calldata")
    subset = np.asarray(genotypes[:, indices, :])
    create(calldata, "GT", subset, chunks=genotypes.chunks)

    staging.rename(target_dir)
    return target_dir


def subsample_panel(panel, populations, per_population, seed):
    """
    Take the first `per_population` samples of each population, preserving panel order.

    compute requires the panel sample order to match the Zarr sample order exactly, so the
    subset must stay in the original order rather than being shuffled into a new one. The
    seed selects which samples are taken while the surviving rows keep panel order.
    """
    selected_indices = []
    for population in populations:
        rows = panel.index[panel["pop"] == population]
        if len(rows) < per_population:
            return None
        chosen = pd.Series(list(rows)).sample(
            n=per_population, random_state=seed
        ).tolist()
        selected_indices.extend(chosen)
    return panel.loc[sorted(selected_indices)].copy()


def main():
    args = parse_args()
    zarr_dir = Path(args.zarr_dir).resolve()
    panel_file = Path(args.panel_file).resolve()
    out_dir = Path(args.out_dir).resolve()
    panels_dir, logs_dir, runs_dir = out_dir / "panels", out_dir / "logs", out_dir / "runs"
    zarr_dir_root = out_dir / "subset_zarr"
    for directory in (panels_dir, logs_dir, runs_dir, zarr_dir_root):
        directory.mkdir(parents=True, exist_ok=True)

    panel = pd.read_csv(panel_file, sep="\t")
    populations = args.populations or sorted(panel["pop"].unique())
    available = {pop: int((panel["pop"] == pop).sum()) for pop in populations}
    smallest = min(available.values())
    print(f"populations: {len(populations)} | smallest has {smallest} samples")

    results, skipped = [], []
    for per_population in sorted(args.samples_per_population):
        if per_population > smallest:
            skipped.append(per_population)
            print(f"  skipping {per_population}/population: exceeds the smallest population "
                  f"({smallest})")
            continue

        subset = subsample_panel(panel, populations, per_population, args.seed)
        subset_path = panels_dir / f"panel_n{per_population}.tsv"
        subset.to_csv(subset_path, sep="\t", index=False)
        total_samples = int(subset.shape[0])

        # Untimed setup: a Zarr store holding exactly this panel's samples.
        print(f"    building subset store for {total_samples} samples...", flush=True)
        subset_zarr = write_sample_subset_zarr(
            zarr_dir, zarr_dir_root / f"subset_n{per_population}.zarr", list(subset.index)
        )

        compute_dir = runs_dir / f"n{per_population}"
        command = [
            sys.executable, "-m", "exp_heatmap.cli", "compute",
            str(subset_zarr), str(subset_path), "-o", str(compute_dir),
            "-t", args.test, "--no-log",
        ]
        if args.chunked:
            command.append("-c")

        def run_once(run_index, per_population=per_population, command=command):
            run_log = logs_dir / f"compute_n{per_population}_run{run_index}.log"
            stats = run_monitored_command(command, run_log)
            if stats["returncode"] != 0:
                raise RuntimeError(
                    f"compute failed for {per_population} samples/population on run "
                    f"{run_index}. See log: {run_log}"
                )
            stats["log_path"] = str(run_log)
            return stats

        print(f"  {per_population}/population ({total_samples} samples): "
              f"{args.warmup} warmup + {args.repeats} measured", flush=True)
        runs = repeat_measurements(run_once, repeats=args.repeats, warmup=args.warmup)
        seconds = summarize_runs([run["seconds"] for run in runs])
        peak_rss = summarize_runs([run["peak_rss_bytes"] for run in runs])

        pair_files = sorted(compute_dir.glob("*.tsv"))
        loci = 0
        if pair_files:
            loci = int(pd.read_csv(pair_files[0], sep="\t", usecols=["variant_pos"]).shape[0])

        results.append({
            "samples_per_population": per_population,
            "total_samples": total_samples,
            "populations": len(populations),
            "population_pairs": len(populations) * (len(populations) - 1) // 2,
            "loci_per_pair": loci,
            "test": args.test,
            "chunked": bool(args.chunked),
            "replicates": seconds["n"],
            "warmup_runs": int(args.warmup),
            "seconds": round(seconds["mean"], 6),
            "seconds_std": round(seconds["std"], 6),
            "seconds_cv_percent": round(seconds["cv_percent"], 4),
            "seconds_ci95_lower": round(seconds["ci95_lower"], 6),
            "seconds_ci95_upper": round(seconds["ci95_upper"], 6),
            "seconds_runs": ";".join(f"{value:.6f}" for value in seconds["values"]),
            "peak_rss_gb": round(peak_rss["mean"] / (1024 ** 3), 4),
            "peak_rss_gb_std": round(peak_rss["std"] / (1024 ** 3), 4),
            "seconds_per_pair": round(
                seconds["mean"] / max(1, len(populations) * (len(populations) - 1) // 2), 6
            ),
            "panel_path": str(subset_path),
            "subset_zarr": str(subset_zarr),
            "log_path": ";".join(run["log_path"] for run in runs),
            "command": shlex.join(command),
        })

    if not results:
        raise SystemExit(
            f"No sample size was small enough to run. The smallest population has "
            f"{smallest} samples; pass --samples-per-population values at or below that."
        )

    results_df = pd.DataFrame(results)
    results_df.to_csv(out_dir / "sample_scaling_results.tsv", sep="\t", index=False)

    specs = machine_specs()
    (out_dir / "machine_specs.json").write_text(json.dumps(specs, indent=2), encoding="utf-8")

    lines = [
        "# Sample-size scaling of the compute stage",
        "",
        f"- Zarr store: `{zarr_dir}`",
        f"- Panel: `{panel_file}`",
        f"- Populations: {len(populations)} "
        f"({len(populations) * (len(populations) - 1) // 2} unordered pairs)",
        f"- Statistic: `{args.test}` | chunked: `{bool(args.chunked)}`",
        f"- Replicates: `{args.repeats}` measured, `{args.warmup}` warmup | seed `{args.seed}`",
        f"- Machine: {specs['processor']} | "
        f"{specs['cpu_count_physical']} physical / {specs['cpu_count_logical']} logical cores | "
        f"{specs['total_memory_gb']} GB RAM",
        f"- Platform: {specs['platform']} | Python {specs['python']}",
    ]
    if skipped:
        lines.append(
            f"- Skipped sizes (larger than the smallest population, {smallest} samples): "
            f"{', '.join(str(value) for value in skipped)}"
        )
    lines += ["", dataframe_to_markdown(results_df), ""]
    (out_dir / "sample_scaling_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"\nwrote {out_dir}/sample_scaling_results.tsv")
    print(results_df[["samples_per_population", "total_samples", "loci_per_pair",
                      "seconds", "seconds_std", "peak_rss_gb"]].to_string(index=False))


if __name__ == "__main__":
    main()
