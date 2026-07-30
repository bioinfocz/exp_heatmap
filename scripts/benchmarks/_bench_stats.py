"""
Shared replicate handling for the benchmark scripts.

Reviewer feedback on the manuscript asked for repeated benchmark runs rather than single
measurements. The formulas here match the earlier `benchmark_stats.py` harness kept on the
`benchmarks` branch: sample standard deviation (ddof=1), a Student-t confidence interval on
the mean, and the coefficient of variation as a percentage.
"""

import platform
import subprocess
import sys

import numpy as np
from scipy import stats


def summarize_runs(values, confidence_level=0.95):
    """
    Summarize repeated measurements of one quantity.

    Returns mean, sample standard deviation, coefficient of variation and a Student-t
    confidence interval on the mean. A single measurement yields zero spread and a
    degenerate interval, which is reported rather than hidden so that single-run rows stay
    distinguishable from replicated ones.
    """
    array = np.asarray([float(value) for value in values], dtype=float)
    n = int(array.size)
    if n == 0:
        raise ValueError("summarize_runs requires at least one measurement")

    mean = float(np.mean(array))
    std = float(np.std(array, ddof=1)) if n > 1 else 0.0

    if n > 1 and std > 0:
        alpha = 1.0 - confidence_level
        t_critical = float(stats.t.ppf(1.0 - alpha / 2.0, df=n - 1))
        half_width = t_critical * std / np.sqrt(n)
        ci_lower, ci_upper = mean - half_width, mean + half_width
    else:
        ci_lower = ci_upper = mean

    return {
        "n": n,
        "mean": mean,
        "std": std,
        "cv_percent": (std / mean * 100.0) if mean else 0.0,
        "ci95_lower": float(ci_lower),
        "ci95_upper": float(ci_upper),
        "min": float(np.min(array)),
        "max": float(np.max(array)),
        "values": [float(value) for value in array],
    }


def format_summary(summary, unit="s", decimals=2):
    """Render a summary as 'mean ± std unit (CV x.x%, n=k)', or plain mean for n=1."""
    if summary["n"] < 2:
        return f"{summary['mean']:.{decimals}f} {unit}"
    return (
        f"{summary['mean']:.{decimals}f} ± {summary['std']:.{decimals}f} {unit} "
        f"(CV {summary['cv_percent']:.1f}%, n={summary['n']})"
    )


def repeat_measurements(run_once, repeats=1, warmup=0, on_result=None):
    """
    Run a measured operation warmup+repeats times and return only the measured results.

    `run_once` takes the 1-based run index and returns the per-run stats dict. Warmup runs
    are executed and discarded, which matters here because the first run of a stage pays
    filesystem cache and import costs that later runs do not.
    """
    if repeats < 1:
        raise ValueError("repeats must be at least 1")
    if warmup < 0:
        raise ValueError("warmup cannot be negative")

    measured = []
    for index in range(1, warmup + repeats + 1):
        result = run_once(index)
        is_warmup = index <= warmup
        if on_result is not None:
            on_result(index, result, is_warmup)
        if not is_warmup:
            measured.append(result)
    return measured


def _linux_cpu_model():
    """Read the CPU model name from /proc/cpuinfo."""
    try:
        with open("/proc/cpuinfo", "r", encoding="utf-8") as handle:
            for line in handle:
                if line.lower().startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return None


def _macos_cpu_model():
    """Read the CPU brand string via sysctl."""
    try:
        result = subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip() or None
    except Exception:
        return None


def cpu_model():
    """
    Best available CPU model name for the current platform.

    platform.processor() returns an empty string or a bare architecture on most Linux
    distributions, so the benchmark provenance records were previously uninformative there.
    """
    if sys.platform == "darwin":
        return _macos_cpu_model() or platform.processor() or platform.machine()
    return _linux_cpu_model() or platform.processor() or platform.machine()


def machine_specs(extra=None):
    """Machine and interpreter context recorded alongside every benchmark result."""
    import psutil

    total_memory = psutil.virtual_memory().total
    specs = {
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "machine": platform.machine(),
        "processor": cpu_model(),
        "cpu_count_physical": psutil.cpu_count(logical=False),
        "cpu_count_logical": psutil.cpu_count(logical=True),
        "total_memory_bytes": int(total_memory),
        "total_memory_gb": round(total_memory / (1024 ** 3), 2),
    }
    if extra:
        specs.update(extra)
    return specs


def dataframe_to_markdown(df):
    """
    Render a DataFrame as a Markdown table without pandas' optional tabulate dependency.

    Mirrors the local helper the other benchmark scripts already use, so summaries render
    identically across all of them.
    """
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

    separator = "| " + " | ".join("-" * width for width in widths) + " |"
    return "\n".join([format_row(headers), separator, *(format_row(row) for row in rows)])
