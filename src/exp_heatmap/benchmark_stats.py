"""
Statistical analysis module for benchmarking with multiple replicates.

This module provides tools for analyzing multiple benchmark runs,
computing statistical summaries (mean, std, confidence intervals),
and running benchmarks with replicates for research paper standards.
"""

from typing import List, Tuple, Callable, Any, Dict, Optional
from dataclasses import dataclass, field
import pandas as pd
import numpy as np
from scipy import stats

from exp_heatmap.logging import get_logger

logger = get_logger(__name__)


@dataclass
class BenchmarkStatistics:
    """
    Container for statistical analysis of multiple benchmark runs.
    
    Provides mean, standard deviation, confidence intervals, and other
    statistics required for research paper reporting.
    """
    operation: str
    n_replicates: int
    
    # Runtime statistics
    runtime_mean: float
    runtime_std: float
    runtime_min: float
    runtime_max: float
    runtime_ci_lower: float  # 95% CI lower bound
    runtime_ci_upper: float  # 95% CI upper bound
    
    # CPU time statistics
    cpu_time_mean: float
    cpu_time_std: float
    
    # Memory statistics
    peak_memory_mean: float
    peak_memory_std: float
    peak_memory_ci_lower: float
    peak_memory_ci_upper: float
    
    # Throughput statistics (if available)
    throughput_mean: Optional[float] = None
    throughput_std: Optional[float] = None
    throughput_unit: Optional[str] = None
    
    # Efficiency statistics (if available)
    memory_efficiency_mean: Optional[float] = None
    memory_efficiency_std: Optional[float] = None
    memory_efficiency_unit: Optional[str] = None
    
    # Coefficient of variation (for assessing measurement stability)
    runtime_cv: Optional[float] = None  # CV = std/mean * 100
    memory_cv: Optional[float] = None
    
    # Additional parameters from the benchmark
    parameters: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for DataFrame creation."""
        result = {
            'operation': self.operation,
            'n_replicates': self.n_replicates,
            'runtime_mean': round(self.runtime_mean, 3),
            'runtime_std': round(self.runtime_std, 3),
            'runtime_min': round(self.runtime_min, 3),
            'runtime_max': round(self.runtime_max, 3),
            'runtime_ci_lower': round(self.runtime_ci_lower, 3),
            'runtime_ci_upper': round(self.runtime_ci_upper, 3),
            'cpu_time_mean': round(self.cpu_time_mean, 3),
            'cpu_time_std': round(self.cpu_time_std, 3),
            'peak_memory_mean': round(self.peak_memory_mean, 2),
            'peak_memory_std': round(self.peak_memory_std, 2),
            'peak_memory_ci_lower': round(self.peak_memory_ci_lower, 2),
            'peak_memory_ci_upper': round(self.peak_memory_ci_upper, 2),
        }
        
        if self.throughput_mean is not None:
            result['throughput_mean'] = round(self.throughput_mean, 2)
            result['throughput_std'] = round(self.throughput_std, 2) if self.throughput_std else None
            result['throughput_unit'] = self.throughput_unit
            
        if self.memory_efficiency_mean is not None:
            result['memory_efficiency_mean'] = round(self.memory_efficiency_mean, 2)
            result['memory_efficiency_std'] = round(self.memory_efficiency_std, 2) if self.memory_efficiency_std else None
            result['memory_efficiency_unit'] = self.memory_efficiency_unit
            
        if self.runtime_cv is not None:
            result['runtime_cv_percent'] = round(self.runtime_cv, 2)
        if self.memory_cv is not None:
            result['memory_cv_percent'] = round(self.memory_cv, 2)
            
        result.update(self.parameters)
        return result


def compute_statistics(results: List['BenchmarkResult'], 
                       confidence_level: float = 0.95) -> BenchmarkStatistics:
    """
    Compute statistical summary from multiple benchmark runs.
    
    Parameters
    ----------
    results : List[BenchmarkResult]
        List of benchmark results from multiple runs
    confidence_level : float
        Confidence level for confidence intervals (default: 0.95 for 95% CI)
        
    Returns
    -------
    BenchmarkStatistics
        Statistical summary including mean, std, and confidence intervals
    """
    if not results:
        raise ValueError("No results provided for statistical analysis")
    
    n = len(results)
    operation = results[0].operation
    
    # Extract arrays
    runtimes = np.array([r.runtime_seconds for r in results])
    cpu_times = np.array([r.cpu_time_seconds for r in results])
    peak_memories = np.array([r.peak_memory_mb for r in results])
    
    # Compute basic statistics
    runtime_mean = np.mean(runtimes)
    runtime_std = np.std(runtimes, ddof=1) if n > 1 else 0.0
    cpu_time_mean = np.mean(cpu_times)
    cpu_time_std = np.std(cpu_times, ddof=1) if n > 1 else 0.0
    peak_memory_mean = np.mean(peak_memories)
    peak_memory_std = np.std(peak_memories, ddof=1) if n > 1 else 0.0
    
    # Compute confidence intervals
    if n > 1:
        alpha = 1 - confidence_level
        t_critical = stats.t.ppf(1 - alpha/2, df=n-1)
        
        runtime_se = runtime_std / np.sqrt(n)
        runtime_ci_lower = runtime_mean - t_critical * runtime_se
        runtime_ci_upper = runtime_mean + t_critical * runtime_se
        
        memory_se = peak_memory_std / np.sqrt(n)
        peak_memory_ci_lower = peak_memory_mean - t_critical * memory_se
        peak_memory_ci_upper = peak_memory_mean + t_critical * memory_se
    else:
        runtime_ci_lower = runtime_ci_upper = runtime_mean
        peak_memory_ci_lower = peak_memory_ci_upper = peak_memory_mean
    
    # Compute coefficient of variation
    runtime_cv = (runtime_std / runtime_mean * 100) if runtime_mean > 0 else None
    memory_cv = (peak_memory_std / peak_memory_mean * 100) if peak_memory_mean > 0 else None
    
    # Throughput statistics (if available)
    throughputs = [r.throughput for r in results if r.throughput is not None]
    throughput_mean = np.mean(throughputs) if throughputs else None
    throughput_std = np.std(throughputs, ddof=1) if len(throughputs) > 1 else None
    throughput_unit = results[0].throughput_unit
    
    # Memory efficiency statistics (if available)
    efficiencies = [r.memory_efficiency for r in results if r.memory_efficiency is not None]
    efficiency_mean = np.mean(efficiencies) if efficiencies else None
    efficiency_std = np.std(efficiencies, ddof=1) if len(efficiencies) > 1 else None
    efficiency_unit = results[0].memory_efficiency_unit
    
    # Collect common parameters (from first result)
    parameters = results[0].parameters.copy()
    
    return BenchmarkStatistics(
        operation=operation,
        n_replicates=n,
        runtime_mean=runtime_mean,
        runtime_std=runtime_std,
        runtime_min=np.min(runtimes),
        runtime_max=np.max(runtimes),
        runtime_ci_lower=runtime_ci_lower,
        runtime_ci_upper=runtime_ci_upper,
        cpu_time_mean=cpu_time_mean,
        cpu_time_std=cpu_time_std,
        peak_memory_mean=peak_memory_mean,
        peak_memory_std=peak_memory_std,
        peak_memory_ci_lower=peak_memory_ci_lower,
        peak_memory_ci_upper=peak_memory_ci_upper,
        throughput_mean=throughput_mean,
        throughput_std=throughput_std,
        throughput_unit=throughput_unit,
        memory_efficiency_mean=efficiency_mean,
        memory_efficiency_std=efficiency_std,
        memory_efficiency_unit=efficiency_unit,
        runtime_cv=runtime_cv,
        memory_cv=memory_cv,
        parameters=parameters,
    )


def run_benchmark_with_replicates(
    benchmark_func: Callable[..., 'BenchmarkResult'],
    n_replicates: int = 3,
    warmup_runs: int = 0,
    **benchmark_kwargs
) -> Tuple[BenchmarkStatistics, List['BenchmarkResult']]:
    """
    Run a benchmark function multiple times and compute statistics.
    
    This function is designed for research paper standards where multiple
    runs are needed to establish statistical significance.
    
    Parameters
    ----------
    benchmark_func : Callable
        The benchmark function to run (e.g., benchmark_prepare, benchmark_compute)
    n_replicates : int
        Number of replicate runs (default: 3, recommended: 3-10)
    warmup_runs : int
        Number of warmup runs to discard (useful for JIT compilation, caching)
    **benchmark_kwargs : dict
        Arguments to pass to the benchmark function
        
    Returns
    -------
    Tuple[BenchmarkStatistics, List[BenchmarkResult]]
        Statistical summary and list of individual results
        
    Example
    -------
    >>> from exp_heatmap.benchmark import benchmark_compute
    >>> stats, results = run_benchmark_with_replicates(
    ...     benchmark_compute,
    ...     n_replicates=5,
    ...     zarr_dir="data.zarr",
    ...     panel_file="panel.tsv",
    ...     output_dir="output"
    ... )
    >>> print(f"Runtime: {stats.runtime_mean:.2f} ± {stats.runtime_std:.2f} s")
    """
    all_results = []
    
    # Warmup runs (discarded)
    for i in range(warmup_runs):
        logger.debug(f"  Warmup run {i+1}/{warmup_runs}...")
        try:
            benchmark_func(**benchmark_kwargs)
        except Exception as e:
            logger.debug(f"    Warmup failed: {e}")
    
    # Actual measurement runs
    for i in range(n_replicates):
        logger.debug(f"  Replicate {i+1}/{n_replicates}...")
        try:
            result = benchmark_func(**benchmark_kwargs)
            all_results.append(result)
            logger.debug(f"    Runtime: {result.runtime_seconds:.2f}s, Memory: {result.peak_memory_mb:.1f}MB")
        except Exception as e:
            logger.debug(f"    FAILED: {e}")
    
    if not all_results:
        raise RuntimeError("All benchmark runs failed")
    
    statistics = compute_statistics(all_results)
    return statistics, all_results


def run_full_benchmark_with_replicates(
    vcf_file: str, 
    panel_file: str, 
    start: int, 
    end: int,
    output_prefix: str = "benchmark",
    test: str = "xpehh",
    chunked: bool = False,
    n_replicates: int = 3,
    warmup_runs: int = 0,
    skip_prepare: bool = False
) -> Tuple[pd.DataFrame, 'SystemInfo']:
    """
    Run a complete benchmark of the ExP Heatmap pipeline with multiple replicates.
    
    This is the recommended function for generating publication-quality benchmarks.
    
    Parameters
    ----------
    vcf_file : str
        Path to input VCF file
    panel_file : str
        Path to population panel file
    start : int
        Start genomic position for plot
    end : int
        End genomic position for plot
    output_prefix : str
        Prefix for output files
    test : str
        Statistical test to compute
    chunked : bool
        Whether to use chunked arrays
    n_replicates : int
        Number of replicate runs for each step
    warmup_runs : int
        Number of warmup runs to discard
    skip_prepare : bool
        If True, skip prepare step (useful if ZARR already exists)
        
    Returns
    -------
    Tuple[pd.DataFrame, SystemInfo]
        DataFrame with statistical results and system information
    """
    from exp_heatmap.benchmark import (
        SystemInfo, benchmark_prepare, benchmark_compute, benchmark_plot
    )
    
    # Capture system information
    system_info = SystemInfo.capture()
    
    results = []
    zarr_dir = f"{output_prefix}_zarr"
    compute_dir = f"{output_prefix}_compute"
    plot_output = f"{output_prefix}_plot"
    
    logger.info("=" * 60)
    logger.info("ExP Heatmap Full Pipeline Benchmark (with replicates)")
    logger.info("=" * 60)
    logger.debug(f"System: {system_info.cpu_model}")
    logger.debug(f"Replicates: {n_replicates}, Warmup: {warmup_runs}")
    
    # Step 1: Prepare (usually run once since it's deterministic)
    if not skip_prepare:
        logger.info("[1/3] Benchmarking PREPARE step...")
        try:
            # Prepare typically only needs one run as it's I/O bound and deterministic
            result = benchmark_prepare(vcf_file, zarr_dir)
            stats = compute_statistics([result])
            results.append(stats.to_dict())
            logger.debug(f"      Runtime: {stats.runtime_mean:.2f}s, Memory: {stats.peak_memory_mean:.1f}MB")
        except Exception as e:
            logger.error(f"      FAILED: {e}")
    
    # Step 2: Compute (run with replicates)
    logger.info(f"[2/3] Benchmarking COMPUTE step ({n_replicates} replicates)...")
    try:
        stats, _ = run_benchmark_with_replicates(
            benchmark_compute,
            n_replicates=n_replicates,
            warmup_runs=warmup_runs,
            zarr_dir=zarr_dir,
            panel_file=panel_file,
            output_dir=compute_dir,
            test=test,
            chunked=chunked
        )
        results.append(stats.to_dict())
        logger.debug(f"      Runtime: {stats.runtime_mean:.2f} ± {stats.runtime_std:.2f}s")
        logger.debug(f"      Memory: {stats.peak_memory_mean:.1f} ± {stats.peak_memory_std:.1f}MB")
    except Exception as e:
        logger.error(f"      FAILED: {e}")
    
    # Step 3: Plot (run with replicates)
    logger.info(f"[3/3] Benchmarking PLOT step ({n_replicates} replicates)...")
    try:
        stats, _ = run_benchmark_with_replicates(
            benchmark_plot,
            n_replicates=n_replicates,
            warmup_runs=warmup_runs,
            input_dir=compute_dir,
            start=start,
            end=end,
            output=plot_output
        )
        results.append(stats.to_dict())
        logger.debug(f"      Runtime: {stats.runtime_mean:.2f} ± {stats.runtime_std:.2f}s")
        logger.debug(f"      Memory: {stats.peak_memory_mean:.1f} ± {stats.peak_memory_std:.1f}MB")
    except Exception as e:
        logger.error(f"      FAILED: {e}")
    
    logger.info("=" * 60)
    logger.info("Benchmark Complete")
    logger.info("=" * 60)
    
    return pd.DataFrame(results), system_info

