"""
Benchmarking module for ExP Heatmap performance evaluation.

This module provides tools to measure runtime and memory usage across
the prepare, compute, and plot steps of the ExP Heatmap pipeline.
Results can be used to evaluate scalability and identify bottlenecks.

Features aligned with research paper benchmarking standards:
- Multi-run statistical analysis (mean, std, confidence intervals)
- CPU time and wall time measurement
- System information capture for reproducibility
- Throughput and efficiency metrics
- JSON/Markdown export for publications
"""

import os
import sys
import time
import json
import platform
import resource
import hashlib
import tracemalloc
import functools
from datetime import datetime
from typing import Dict, List, Optional, Callable, Any, Union, Tuple
from dataclasses import dataclass, field, asdict
import pandas as pd
import numpy as np
from scipy import stats


def _get_package_version(package_name: str) -> str:
    """Safely get package version, return 'unknown' if not found."""
    try:
        import importlib.metadata
        return importlib.metadata.version(package_name)
    except Exception:
        return "unknown"


@dataclass
class SystemInfo:
    """
    Container for system and environment information.
    
    Captures all relevant details for reproducibility in research papers.
    """
    timestamp: str
    python_version: str
    os_name: str
    os_version: str
    architecture: str
    cpu_model: str
    cpu_cores_physical: int
    cpu_cores_logical: int
    total_ram_gb: float
    exp_heatmap_version: str
    scikit_allel_version: str
    zarr_version: str
    pandas_version: str
    numpy_version: str
    scipy_version: str
    
    @classmethod
    def capture(cls) -> 'SystemInfo':
        """Capture current system information."""
        import multiprocessing
        
        # Get CPU model
        cpu_model = "unknown"
        try:
            if platform.system() == "Linux":
                with open("/proc/cpuinfo", "r") as f:
                    for line in f:
                        if "model name" in line:
                            cpu_model = line.split(":")[1].strip()
                            break
            elif platform.system() == "Darwin":
                import subprocess
                cpu_model = subprocess.check_output(
                    ["sysctl", "-n", "machdep.cpu.brand_string"]
                ).decode().strip()
            elif platform.system() == "Windows":
                cpu_model = platform.processor()
        except Exception:
            pass
        
        # Get RAM
        total_ram_gb = 0.0
        try:
            if platform.system() == "Linux":
                with open("/proc/meminfo", "r") as f:
                    for line in f:
                        if "MemTotal" in line:
                            # Value is in kB
                            total_ram_gb = int(line.split()[1]) / (1024 * 1024)
                            break
            elif platform.system() == "Darwin":
                import subprocess
                mem_bytes = int(subprocess.check_output(["sysctl", "-n", "hw.memsize"]).decode().strip())
                total_ram_gb = mem_bytes / (1024 ** 3)
        except Exception:
            pass
        
        # Get physical cores
        try:
            physical_cores = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else multiprocessing.cpu_count()
        except Exception:
            physical_cores = multiprocessing.cpu_count()
        
        return cls(
            timestamp=datetime.now().isoformat(),
            python_version=platform.python_version(),
            os_name=platform.system(),
            os_version=platform.release(),
            architecture=platform.machine(),
            cpu_model=cpu_model,
            cpu_cores_physical=physical_cores,
            cpu_cores_logical=multiprocessing.cpu_count(),
            total_ram_gb=round(total_ram_gb, 1),
            exp_heatmap_version=_get_package_version("exp_heatmap"),
            scikit_allel_version=_get_package_version("scikit-allel"),
            zarr_version=_get_package_version("zarr"),
            pandas_version=_get_package_version("pandas"),
            numpy_version=_get_package_version("numpy"),
            scipy_version=_get_package_version("scipy"),
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)
    
    def to_report_lines(self) -> List[str]:
        """Generate formatted report lines."""
        return [
            f"Timestamp: {self.timestamp}",
            f"Python: {self.python_version}",
            f"OS: {self.os_name} {self.os_version} ({self.architecture})",
            f"CPU: {self.cpu_model}",
            f"CPU Cores: {self.cpu_cores_physical} physical, {self.cpu_cores_logical} logical",
            f"RAM: {self.total_ram_gb} GB",
            f"exp_heatmap: {self.exp_heatmap_version}",
            f"scikit-allel: {self.scikit_allel_version}",
            f"zarr: {self.zarr_version}",
            f"pandas: {self.pandas_version}",
            f"numpy: {self.numpy_version}",
            f"scipy: {self.scipy_version}",
        ]


@dataclass
class BenchmarkResult:
    """
    Container for benchmark results from a single operation.
    
    Extended to include CPU time, throughput, and efficiency metrics
    for research paper standards.
    """
    operation: str
    runtime_seconds: float  # Wall clock time
    cpu_time_seconds: float  # CPU time (user + system)
    peak_memory_mb: float
    current_memory_mb: float
    # Throughput metrics (operation-specific)
    throughput: Optional[float] = None
    throughput_unit: Optional[str] = None
    # Efficiency metrics
    memory_efficiency: Optional[float] = None  # MB per unit of work
    memory_efficiency_unit: Optional[str] = None
    # Additional parameters
    parameters: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for DataFrame creation."""
        result = {
            'operation': self.operation,
            'runtime_seconds': round(self.runtime_seconds, 3),
            'cpu_time_seconds': round(self.cpu_time_seconds, 3),
            'peak_memory_mb': round(self.peak_memory_mb, 2),
            'current_memory_mb': round(self.current_memory_mb, 2),
        }
        if self.throughput is not None:
            result['throughput'] = round(self.throughput, 2)
            result['throughput_unit'] = self.throughput_unit
        if self.memory_efficiency is not None:
            result['memory_efficiency'] = round(self.memory_efficiency, 2)
            result['memory_efficiency_unit'] = self.memory_efficiency_unit
        result.update(self.parameters)
        return result


class Benchmarker:
    """
    Context manager and utility class for benchmarking ExP Heatmap operations.
    
    Captures both wall clock time and CPU time for research paper standards.
    
    Usage:
        # As context manager
        with Benchmarker("prepare") as bench:
            prepare(vcf_file, zarr_dir)
        result = bench.result
        
        # Or using the decorator
        @Benchmarker.measure("compute")
        def my_compute_function(...):
            ...
    """
    
    def __init__(self, operation_name: str, **parameters):
        """
        Initialize benchmarker.
        
        Parameters
        ----------
        operation_name : str
            Name of the operation being benchmarked (e.g., 'prepare', 'compute', 'plot')
        **parameters : dict
            Additional parameters to record with the benchmark (e.g., n_variants, n_populations)
        """
        self.operation_name = operation_name
        self.parameters = parameters
        self.start_time: Optional[float] = None
        self.start_cpu_time: Optional[float] = None
        self.start_resources: Optional[resource.struct_rusage] = None
        self.result: Optional[BenchmarkResult] = None
        
    def __enter__(self):
        """Start benchmarking."""
        tracemalloc.start()
        self.start_time = time.perf_counter()
        self.start_cpu_time = time.process_time()
        self.start_resources = resource.getrusage(resource.RUSAGE_SELF)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Stop benchmarking and record results."""
        end_time = time.perf_counter()
        end_cpu_time = time.process_time()
        end_resources = resource.getrusage(resource.RUSAGE_SELF)
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        # Calculate CPU time (user + system)
        wall_time = end_time - self.start_time
        cpu_time = end_cpu_time - self.start_cpu_time
        
        # Alternative: use resource module for more accurate CPU time
        user_time = end_resources.ru_utime - self.start_resources.ru_utime
        sys_time = end_resources.ru_stime - self.start_resources.ru_stime
        total_cpu_time = user_time + sys_time
        
        # Use the more accurate of the two
        cpu_time = max(cpu_time, total_cpu_time)
        
        self.result = BenchmarkResult(
            operation=self.operation_name,
            runtime_seconds=wall_time,
            cpu_time_seconds=cpu_time,
            peak_memory_mb=peak / (1024 * 1024),
            current_memory_mb=current / (1024 * 1024),
            parameters=self.parameters
        )
        return False  # Don't suppress exceptions
    
    @staticmethod
    def measure(operation_name: str, **default_params):
        """
        Decorator to benchmark a function.
        
        Parameters
        ----------
        operation_name : str
            Name of the operation
        **default_params : dict
            Default parameters to record
            
        Returns
        -------
        Callable
            Decorated function that records benchmark results
        """
        def decorator(func: Callable) -> Callable:
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                params = {**default_params}
                with Benchmarker(operation_name, **params) as bench:
                    result = func(*args, **kwargs)
                # Store benchmark result as function attribute
                wrapper.last_benchmark = bench.result
                return result
            wrapper.last_benchmark = None
            return wrapper
        return decorator


def _count_variants_in_zarr(zarr_dir: str) -> int:
    """Count the number of variants in a ZARR directory."""
    try:
        import zarr
        callset = zarr.open_group(zarr_dir, mode="r")
        return callset['variants/POS'].shape[0]
    except Exception:
        return 0


def _get_file_checksum(filepath: str, algorithm: str = 'md5') -> str:
    """Calculate checksum of a file for reproducibility tracking."""
    hash_func = hashlib.new(algorithm)
    try:
        with open(filepath, 'rb') as f:
            for chunk in iter(lambda: f.read(8192), b''):
                hash_func.update(chunk)
        return hash_func.hexdigest()
    except Exception:
        return "unknown"


def benchmark_prepare(vcf_file: str, zarr_dir: str, 
                      compute_checksum: bool = False) -> BenchmarkResult:
    """
    Benchmark the prepare step (VCF to ZARR conversion).
    
    Parameters
    ----------
    vcf_file : str
        Path to input VCF file
    zarr_dir : str
        Path for output ZARR directory
    compute_checksum : bool
        If True, compute MD5 checksum of input file (slower but aids reproducibility)
        
    Returns
    -------
    BenchmarkResult
        Benchmark results including runtime, memory, throughput and efficiency metrics
    """
    from exp_heatmap import prepare as prepare_func
    
    # Get file size for context
    file_size_mb = os.path.getsize(vcf_file) / (1024 * 1024)
    
    params = {'vcf_file_size_mb': file_size_mb}
    if compute_checksum:
        params['vcf_checksum_md5'] = _get_file_checksum(vcf_file)
    
    with Benchmarker("prepare", **params) as bench:
        prepare_func(vcf_file, zarr_dir)
    
    # Calculate throughput metrics
    n_variants = _count_variants_in_zarr(zarr_dir)
    if n_variants > 0 and bench.result.runtime_seconds > 0:
        bench.result.throughput = n_variants / bench.result.runtime_seconds
        bench.result.throughput_unit = "variants/second"
        bench.result.parameters['n_variants'] = n_variants
        
        # Memory efficiency: MB per million variants
        if n_variants > 0:
            bench.result.memory_efficiency = (bench.result.peak_memory_mb / n_variants) * 1_000_000
            bench.result.memory_efficiency_unit = "MB/million_variants"
    
    return bench.result


def benchmark_compute(zarr_dir: str, panel_file: str, output_dir: str, 
                      test: str = "xpehh", chunked: bool = False) -> BenchmarkResult:
    """
    Benchmark the compute step (statistical test computation).
    
    Parameters
    ----------
    zarr_dir : str
        Path to ZARR directory
    panel_file : str
        Path to population panel file
    output_dir : str
        Path for output TSV files
    test : str
        Statistical test to compute
    chunked : bool
        Whether to use chunked arrays
        
    Returns
    -------
    BenchmarkResult
        Benchmark results including runtime, memory, throughput and efficiency metrics
    """
    from exp_heatmap import compute as compute_func
    
    # Count populations from panel file
    panel = pd.read_csv(panel_file, sep="\t")
    n_populations = panel['pop'].nunique()
    n_samples = len(panel)
    n_pairs = n_populations * (n_populations - 1) // 2
    
    # Count variants
    n_variants = _count_variants_in_zarr(zarr_dir)
    
    with Benchmarker("compute", 
                     test=test,
                     n_populations=n_populations,
                     n_samples=n_samples,
                     n_pairs=n_pairs,
                     n_variants=n_variants,
                     chunked=chunked) as bench:
        compute_func(zarr_dir, panel_file, output_dir, test, chunked)
    
    # Calculate throughput metrics
    if bench.result.runtime_seconds > 0:
        # Pairs per second
        bench.result.throughput = n_pairs / bench.result.runtime_seconds
        bench.result.throughput_unit = "pairs/second"
        
        # Memory efficiency: MB per million variants
        if n_variants > 0:
            bench.result.memory_efficiency = (bench.result.peak_memory_mb / n_variants) * 1_000_000
            bench.result.memory_efficiency_unit = "MB/million_variants"
        
        # Additional derived metrics
        if n_variants > 0 and n_pairs > 0:
            bench.result.parameters['seconds_per_pair'] = round(bench.result.runtime_seconds / n_pairs, 3)
            bench.result.parameters['variants_x_pairs'] = n_variants * n_pairs
    
    return bench.result


def benchmark_plot(input_dir: str, start: int, end: int, 
                   output: str = "benchmark_plot", **plot_kwargs) -> BenchmarkResult:
    """
    Benchmark the plot step (heatmap generation).
    
    Parameters
    ----------
    input_dir : str
        Path to directory with TSV files
    start : int
        Start genomic position
    end : int
        End genomic position
    output : str
        Output filename
    **plot_kwargs : dict
        Additional arguments for plot function
        
    Returns
    -------
    BenchmarkResult
        Benchmark results including runtime, memory, throughput and efficiency metrics
    """
    from exp_heatmap.plot import plot
    import glob
    
    # Count files and estimate data size
    tsv_files = glob.glob(os.path.join(input_dir, "*.tsv"))
    n_files = len(tsv_files)
    region_size = end - start
    
    with Benchmarker("plot",
                     n_files=n_files,
                     region_size_bp=region_size,
                     start=start,
                     end=end) as bench:
        plot(input_dir, start=start, end=end, title="Benchmark", output=output, **plot_kwargs)
    
    # Calculate throughput metrics
    if bench.result.runtime_seconds > 0:
        # Region size per second (bp/s)
        bench.result.throughput = region_size / bench.result.runtime_seconds
        bench.result.throughput_unit = "bp/second"
        
        # Memory efficiency: MB per megabase
        if region_size > 0:
            bench.result.memory_efficiency = (bench.result.peak_memory_mb / region_size) * 1_000_000
            bench.result.memory_efficiency_unit = "MB/megabase"
    
    return bench.result


# =============================================================================
# Statistical Analysis for Multi-Run Benchmarks
# =============================================================================

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


def compute_statistics(results: List[BenchmarkResult], 
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
    benchmark_func: Callable[..., BenchmarkResult],
    n_replicates: int = 3,
    warmup_runs: int = 0,
    **benchmark_kwargs
) -> Tuple[BenchmarkStatistics, List[BenchmarkResult]]:
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
        print(f"  Warmup run {i+1}/{warmup_runs}...")
        try:
            benchmark_func(**benchmark_kwargs)
        except Exception as e:
            print(f"    Warmup failed: {e}")
    
    # Actual measurement runs
    for i in range(n_replicates):
        print(f"  Replicate {i+1}/{n_replicates}...")
        try:
            result = benchmark_func(**benchmark_kwargs)
            all_results.append(result)
            print(f"    Runtime: {result.runtime_seconds:.2f}s, Memory: {result.peak_memory_mb:.1f}MB")
        except Exception as e:
            print(f"    FAILED: {e}")
    
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
) -> Tuple[pd.DataFrame, SystemInfo]:
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
    # Capture system information
    system_info = SystemInfo.capture()
    
    results = []
    zarr_dir = f"{output_prefix}_zarr"
    compute_dir = f"{output_prefix}_compute"
    plot_output = f"{output_prefix}_plot"
    
    print("=" * 60)
    print("ExP Heatmap Full Pipeline Benchmark (with replicates)")
    print("=" * 60)
    print(f"\nSystem: {system_info.cpu_model}")
    print(f"Replicates: {n_replicates}, Warmup: {warmup_runs}")
    
    # Step 1: Prepare (usually run once since it's deterministic)
    if not skip_prepare:
        print("\n[1/3] Benchmarking PREPARE step...")
        try:
            # Prepare typically only needs one run as it's I/O bound and deterministic
            result = benchmark_prepare(vcf_file, zarr_dir)
            stats = compute_statistics([result])
            results.append(stats.to_dict())
            print(f"      Runtime: {stats.runtime_mean:.2f}s, Memory: {stats.peak_memory_mean:.1f}MB")
        except Exception as e:
            print(f"      FAILED: {e}")
    
    # Step 2: Compute (run with replicates)
    print(f"\n[2/3] Benchmarking COMPUTE step ({n_replicates} replicates)...")
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
        print(f"      Runtime: {stats.runtime_mean:.2f} ± {stats.runtime_std:.2f}s")
        print(f"      Memory: {stats.peak_memory_mean:.1f} ± {stats.peak_memory_std:.1f}MB")
    except Exception as e:
        print(f"      FAILED: {e}")
    
    # Step 3: Plot (run with replicates)
    print(f"\n[3/3] Benchmarking PLOT step ({n_replicates} replicates)...")
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
        print(f"      Runtime: {stats.runtime_mean:.2f} ± {stats.runtime_std:.2f}s")
        print(f"      Memory: {stats.peak_memory_mean:.1f} ± {stats.peak_memory_std:.1f}MB")
    except Exception as e:
        print(f"      FAILED: {e}")
    
    print("\n" + "=" * 60)
    print("Benchmark Complete")
    print("=" * 60)
    
    return pd.DataFrame(results), system_info


def run_full_benchmark(vcf_file: str, panel_file: str, 
                       start: int, end: int,
                       output_prefix: str = "benchmark",
                       test: str = "xpehh",
                       chunked: bool = False) -> pd.DataFrame:
    """
    Run a complete benchmark of the ExP Heatmap pipeline.
    
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
        
    Returns
    -------
    pd.DataFrame
        DataFrame with benchmark results for all steps
    """
    results = []
    
    zarr_dir = f"{output_prefix}_zarr"
    compute_dir = f"{output_prefix}_compute"
    plot_output = f"{output_prefix}_plot"
    
    print("=" * 60)
    print("ExP Heatmap Full Pipeline Benchmark")
    print("=" * 60)
    
    # Step 1: Prepare
    print("\n[1/3] Benchmarking PREPARE step...")
    try:
        result = benchmark_prepare(vcf_file, zarr_dir)
        results.append(result.to_dict())
        print(f"      Runtime: {result.runtime_seconds:.2f}s, Peak Memory: {result.peak_memory_mb:.1f}MB")
    except Exception as e:
        print(f"      FAILED: {e}")
    
    # Step 2: Compute
    print("\n[2/3] Benchmarking COMPUTE step...")
    try:
        result = benchmark_compute(zarr_dir, panel_file, compute_dir, test, chunked)
        results.append(result.to_dict())
        print(f"      Runtime: {result.runtime_seconds:.2f}s, Peak Memory: {result.peak_memory_mb:.1f}MB")
    except Exception as e:
        print(f"      FAILED: {e}")
    
    # Step 3: Plot
    print("\n[3/3] Benchmarking PLOT step...")
    try:
        result = benchmark_plot(compute_dir, start, end, plot_output)
        results.append(result.to_dict())
        print(f"      Runtime: {result.runtime_seconds:.2f}s, Peak Memory: {result.peak_memory_mb:.1f}MB")
    except Exception as e:
        print(f"      FAILED: {e}")
    
    print("\n" + "=" * 60)
    print("Benchmark Complete")
    print("=" * 60)
    
    return pd.DataFrame(results)


def generate_benchmark_report(
    results: pd.DataFrame, 
    system_info: Optional[SystemInfo] = None,
    output_file: str = None,
    include_raw_data: bool = False
) -> str:
    """
    Generate a formatted benchmark report suitable for research papers.
    
    Supports both single-run results (with runtime_seconds column) and
    multi-run statistical results (with runtime_mean column).
    
    Parameters
    ----------
    results : pd.DataFrame
        DataFrame with benchmark results (single-run or statistical)
    system_info : SystemInfo, optional
        System information to include in report
    output_file : str, optional
        If provided, save report to this file
    include_raw_data : bool
        If True, include raw data table at the end
        
    Returns
    -------
    str
        Formatted benchmark report
    """
    report_lines = [
        "ExP Heatmap Benchmark Report",
        "=" * 50,
        "",
    ]
    
    # Add system information section if provided
    if system_info:
        report_lines.extend([
            "System Information",
            "-" * 30,
        ])
        report_lines.extend(system_info.to_report_lines())
        report_lines.append("")
    
    # Detect if this is a statistical report (has runtime_mean) or single-run (has runtime_seconds)
    is_statistical = 'runtime_mean' in results.columns
    
    report_lines.extend([
        "Summary",
        "-" * 30,
    ])
    
    if is_statistical:
        # Statistical report
        n_replicates = results['n_replicates'].iloc[0] if 'n_replicates' in results.columns else 1
        total_runtime = results['runtime_mean'].sum()
        total_runtime_std = np.sqrt((results['runtime_std'] ** 2).sum()) if 'runtime_std' in results.columns else 0
        max_memory = results['peak_memory_mean'].max()
        max_memory_std = results.loc[results['peak_memory_mean'].idxmax(), 'peak_memory_std'] if 'peak_memory_std' in results.columns else 0
        
        report_lines.extend([
            f"Number of Replicates: {n_replicates}",
            f"Total Runtime: {total_runtime:.2f} ± {total_runtime_std:.2f} seconds (mean ± std)",
            f"Peak Memory Usage: {max_memory:.1f} ± {max_memory_std:.1f} MB",
            "",
            "Step-by-Step Results",
            "-" * 30,
        ])
        
        for _, row in results.iterrows():
            op_name = row['operation'].upper()
            report_lines.append(f"\n{op_name}:")
            
            # Runtime with statistics
            runtime_mean = row.get('runtime_mean', row.get('runtime_seconds', 0))
            runtime_std = row.get('runtime_std', 0)
            report_lines.append(f"  Runtime: {runtime_mean:.3f} ± {runtime_std:.3f} seconds")
            
            # CPU time if available
            if 'cpu_time_mean' in row and pd.notna(row['cpu_time_mean']):
                cpu_mean = row['cpu_time_mean']
                cpu_std = row.get('cpu_time_std', 0)
                report_lines.append(f"  CPU Time: {cpu_mean:.3f} ± {cpu_std:.3f} seconds")
            
            # Memory with statistics
            mem_mean = row.get('peak_memory_mean', row.get('peak_memory_mb', 0))
            mem_std = row.get('peak_memory_std', 0)
            report_lines.append(f"  Peak Memory: {mem_mean:.1f} ± {mem_std:.1f} MB")
            
            # Confidence intervals if available
            if 'runtime_ci_lower' in row and pd.notna(row['runtime_ci_lower']):
                report_lines.append(f"  Runtime 95% CI: [{row['runtime_ci_lower']:.3f}, {row['runtime_ci_upper']:.3f}] seconds")
            
            # Throughput if available
            if 'throughput_mean' in row and pd.notna(row['throughput_mean']):
                tp_mean = row['throughput_mean']
                tp_std = row.get('throughput_std', 0) or 0
                tp_unit = row.get('throughput_unit', '')
                report_lines.append(f"  Throughput: {tp_mean:.2f} ± {tp_std:.2f} {tp_unit}")
            
            # Memory efficiency if available
            if 'memory_efficiency_mean' in row and pd.notna(row['memory_efficiency_mean']):
                eff_mean = row['memory_efficiency_mean']
                eff_std = row.get('memory_efficiency_std', 0) or 0
                eff_unit = row.get('memory_efficiency_unit', '')
                report_lines.append(f"  Memory Efficiency: {eff_mean:.2f} ± {eff_std:.2f} {eff_unit}")
            
            # Coefficient of variation if available
            if 'runtime_cv_percent' in row and pd.notna(row['runtime_cv_percent']):
                report_lines.append(f"  Runtime CV: {row['runtime_cv_percent']:.1f}%")
            
            # Additional parameters
            skip_keys = {
                'operation', 'n_replicates', 
                'runtime_mean', 'runtime_std', 'runtime_min', 'runtime_max',
                'runtime_ci_lower', 'runtime_ci_upper', 'runtime_cv_percent',
                'cpu_time_mean', 'cpu_time_std',
                'peak_memory_mean', 'peak_memory_std', 
                'peak_memory_ci_lower', 'peak_memory_ci_upper', 'memory_cv_percent',
                'throughput_mean', 'throughput_std', 'throughput_unit',
                'memory_efficiency_mean', 'memory_efficiency_std', 'memory_efficiency_unit',
                'current_memory_mb'
            }
            for key, value in row.items():
                if key not in skip_keys and pd.notna(value):
                    report_lines.append(f"  {key}: {value}")
    else:
        # Single-run report (backward compatible)
        total_runtime = results['runtime_seconds'].sum()
        max_memory = results['peak_memory_mb'].max()
        
        report_lines.extend([
            f"Total Runtime: {total_runtime:.2f} seconds",
            f"Peak Memory Usage: {max_memory:.1f} MB",
            "",
            "Step-by-Step Results",
            "-" * 30,
        ])
        
        for _, row in results.iterrows():
            report_lines.extend([
                f"\n{row['operation'].upper()}:",
                f"  Runtime: {row['runtime_seconds']:.3f} seconds",
            ])
            
            # CPU time if available
            if 'cpu_time_seconds' in row and pd.notna(row['cpu_time_seconds']):
                report_lines.append(f"  CPU Time: {row['cpu_time_seconds']:.3f} seconds")
            
            report_lines.append(f"  Peak Memory: {row['peak_memory_mb']:.1f} MB")
            
            # Throughput if available
            if 'throughput' in row and pd.notna(row['throughput']):
                tp = row['throughput']
                tp_unit = row.get('throughput_unit', '')
                report_lines.append(f"  Throughput: {tp:.2f} {tp_unit}")
            
            # Memory efficiency if available
            if 'memory_efficiency' in row and pd.notna(row['memory_efficiency']):
                eff = row['memory_efficiency']
                eff_unit = row.get('memory_efficiency_unit', '')
                report_lines.append(f"  Memory Efficiency: {eff:.2f} {eff_unit}")
            
            # Add operation-specific parameters
            skip_keys = {
                'operation', 'runtime_seconds', 'cpu_time_seconds',
                'peak_memory_mb', 'current_memory_mb',
                'throughput', 'throughput_unit',
                'memory_efficiency', 'memory_efficiency_unit'
            }
            for key, value in row.items():
                if key not in skip_keys and pd.notna(value):
                    report_lines.append(f"  {key}: {value}")
    
    # Add raw data table if requested
    if include_raw_data:
        report_lines.extend([
            "",
            "",
            "Raw Data",
            "-" * 30,
            results.to_string(),
        ])
    
    report = "\n".join(report_lines)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report)
        print(f"Report saved to: {output_file}")
    
    return report


def export_benchmark_json(
    results: pd.DataFrame,
    system_info: Optional[SystemInfo] = None,
    output_file: str = "benchmark_results.json",
    metadata: Optional[Dict[str, Any]] = None
) -> str:
    """
    Export benchmark results to JSON format for machine-readable output.
    
    This format is suitable for automated analysis, comparison with other tools,
    and integration into CI/CD pipelines.
    
    Parameters
    ----------
    results : pd.DataFrame
        DataFrame with benchmark results
    system_info : SystemInfo, optional
        System information to include
    output_file : str
        Output JSON file path
    metadata : dict, optional
        Additional metadata to include (e.g., experiment description, git commit)
        
    Returns
    -------
    str
        JSON string of the benchmark results
    """
    export_data = {
        "benchmark_version": "2.0",
        "export_timestamp": datetime.now().isoformat(),
    }
    
    # Add system information
    if system_info:
        export_data["system_info"] = system_info.to_dict()
    
    # Add metadata
    if metadata:
        export_data["metadata"] = metadata
    
    # Detect if statistical or single-run results
    is_statistical = 'runtime_mean' in results.columns
    export_data["result_type"] = "statistical" if is_statistical else "single_run"
    
    # Convert results to records
    export_data["results"] = results.to_dict(orient='records')
    
    # Add summary statistics
    if is_statistical:
        export_data["summary"] = {
            "total_runtime_mean": float(results['runtime_mean'].sum()),
            "total_runtime_std": float(np.sqrt((results['runtime_std'] ** 2).sum())),
            "peak_memory_max": float(results['peak_memory_mean'].max()),
            "n_operations": len(results),
        }
    else:
        export_data["summary"] = {
            "total_runtime": float(results['runtime_seconds'].sum()),
            "peak_memory_max": float(results['peak_memory_mb'].max()),
            "n_operations": len(results),
        }
    
    # Write to file
    json_str = json.dumps(export_data, indent=2, default=str)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(json_str)
        print(f"JSON export saved to: {output_file}")
    
    return json_str


def export_benchmark_markdown(
    results: pd.DataFrame,
    system_info: Optional[SystemInfo] = None,
    output_file: str = "benchmark_results.md",
    title: str = "ExP Heatmap Benchmark Results",
    include_methodology: bool = True
) -> str:
    """
    Export benchmark results to Markdown format for research papers and documentation.
    
    Generates publication-ready tables and formatted text suitable for inclusion
    in papers, README files, and documentation.
    
    Parameters
    ----------
    results : pd.DataFrame
        DataFrame with benchmark results
    system_info : SystemInfo, optional
        System information to include
    output_file : str
        Output Markdown file path
    title : str
        Title for the report
    include_methodology : bool
        If True, include a methodology section
        
    Returns
    -------
    str
        Markdown string of the benchmark results
    """
    lines = [
        f"# {title}",
        "",
    ]
    
    # System information section
    if system_info:
        lines.extend([
            "## System Configuration",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Date | {system_info.timestamp} |",
            f"| Python | {system_info.python_version} |",
            f"| OS | {system_info.os_name} {system_info.os_version} |",
            f"| CPU | {system_info.cpu_model} |",
            f"| CPU Cores | {system_info.cpu_cores_physical} physical, {system_info.cpu_cores_logical} logical |",
            f"| RAM | {system_info.total_ram_gb} GB |",
            f"| exp_heatmap | {system_info.exp_heatmap_version} |",
            f"| scikit-allel | {system_info.scikit_allel_version} |",
            f"| zarr | {system_info.zarr_version} |",
            "",
        ])
    
    # Methodology section
    if include_methodology:
        is_statistical = 'runtime_mean' in results.columns
        n_replicates = results['n_replicates'].iloc[0] if 'n_replicates' in results.columns else 1
        
        lines.extend([
            "## Methodology",
            "",
        ])
        
        if is_statistical:
            lines.extend([
                f"- **Replicates**: {n_replicates} runs per operation",
                "- **Metrics**: Wall clock time, CPU time, peak memory usage",
                "- **Statistics**: Mean ± standard deviation, 95% confidence intervals",
                "- **Memory tracking**: Python `tracemalloc` module",
                "- **Timing**: `time.perf_counter()` for wall time, `resource.getrusage()` for CPU time",
                "",
            ])
        else:
            lines.extend([
                "- **Runs**: Single run per operation",
                "- **Metrics**: Wall clock time, CPU time, peak memory usage",
                "- **Memory tracking**: Python `tracemalloc` module",
                "",
            ])
    
    # Results section
    lines.extend([
        "## Results",
        "",
    ])
    
    # Detect result type and create appropriate table
    is_statistical = 'runtime_mean' in results.columns
    
    if is_statistical:
        lines.extend([
            "### Performance Summary",
            "",
            "| Operation | Runtime (s) | CPU Time (s) | Peak Memory (MB) | Throughput |",
            "|-----------|-------------|--------------|------------------|------------|",
        ])
        
        for _, row in results.iterrows():
            op = row['operation'].upper()
            runtime = f"{row['runtime_mean']:.2f} ± {row.get('runtime_std', 0):.2f}"
            cpu_time = f"{row.get('cpu_time_mean', 0):.2f} ± {row.get('cpu_time_std', 0):.2f}"
            memory = f"{row['peak_memory_mean']:.1f} ± {row.get('peak_memory_std', 0):.1f}"
            
            if pd.notna(row.get('throughput_mean')):
                throughput = f"{row['throughput_mean']:.2f} {row.get('throughput_unit', '')}"
            else:
                throughput = "N/A"
            
            lines.append(f"| {op} | {runtime} | {cpu_time} | {memory} | {throughput} |")
        
        # Add totals row
        total_runtime = results['runtime_mean'].sum()
        total_runtime_std = np.sqrt((results['runtime_std'] ** 2).sum())
        max_memory = results['peak_memory_mean'].max()
        
        lines.extend([
            "",
            f"**Total Runtime**: {total_runtime:.2f} ± {total_runtime_std:.2f} seconds",
            f"**Peak Memory**: {max_memory:.1f} MB",
            "",
        ])
        
        # Confidence intervals table
        if 'runtime_ci_lower' in results.columns:
            lines.extend([
                "### 95% Confidence Intervals",
                "",
                "| Operation | Runtime CI (s) | Memory CI (MB) |",
                "|-----------|----------------|----------------|",
            ])
            
            for _, row in results.iterrows():
                op = row['operation'].upper()
                runtime_ci = f"[{row['runtime_ci_lower']:.2f}, {row['runtime_ci_upper']:.2f}]"
                memory_ci = f"[{row['peak_memory_ci_lower']:.1f}, {row['peak_memory_ci_upper']:.1f}]"
                lines.append(f"| {op} | {runtime_ci} | {memory_ci} |")
            
            lines.append("")
    
    else:
        # Single-run results
        lines.extend([
            "### Performance Summary",
            "",
            "| Operation | Runtime (s) | CPU Time (s) | Peak Memory (MB) | Throughput |",
            "|-----------|-------------|--------------|------------------|------------|",
        ])
        
        for _, row in results.iterrows():
            op = row['operation'].upper()
            runtime = f"{row['runtime_seconds']:.2f}"
            cpu_time = f"{row.get('cpu_time_seconds', 0):.2f}"
            memory = f"{row['peak_memory_mb']:.1f}"
            
            if pd.notna(row.get('throughput')):
                throughput = f"{row['throughput']:.2f} {row.get('throughput_unit', '')}"
            else:
                throughput = "N/A"
            
            lines.append(f"| {op} | {runtime} | {cpu_time} | {memory} | {throughput} |")
        
        lines.extend([
            "",
            f"**Total Runtime**: {results['runtime_seconds'].sum():.2f} seconds",
            f"**Peak Memory**: {results['peak_memory_mb'].max():.1f} MB",
            "",
        ])
    
    # Additional parameters section
    lines.extend([
        "### Benchmark Parameters",
        "",
    ])
    
    # Collect unique parameters across all operations
    skip_keys = {
        'operation', 'n_replicates',
        'runtime_mean', 'runtime_std', 'runtime_min', 'runtime_max',
        'runtime_ci_lower', 'runtime_ci_upper', 'runtime_cv_percent',
        'runtime_seconds', 'cpu_time_seconds', 'cpu_time_mean', 'cpu_time_std',
        'peak_memory_mean', 'peak_memory_std', 'peak_memory_mb',
        'peak_memory_ci_lower', 'peak_memory_ci_upper', 'memory_cv_percent',
        'current_memory_mb',
        'throughput', 'throughput_mean', 'throughput_std', 'throughput_unit',
        'memory_efficiency', 'memory_efficiency_mean', 'memory_efficiency_std', 'memory_efficiency_unit',
    }
    
    for _, row in results.iterrows():
        op = row['operation'].upper()
        params = []
        for key, value in row.items():
            if key not in skip_keys and pd.notna(value):
                params.append(f"  - {key}: {value}")
        
        if params:
            lines.append(f"**{op}**:")
            lines.extend(params)
            lines.append("")
    
    md_str = "\n".join(lines)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(md_str)
        print(f"Markdown export saved to: {output_file}")
    
    return md_str


def export_benchmark_csv(
    results: pd.DataFrame,
    output_file: str = "benchmark_results.csv"
) -> None:
    """
    Export benchmark results to CSV format.
    
    Parameters
    ----------
    results : pd.DataFrame
        DataFrame with benchmark results
    output_file : str
        Output CSV file path
    """
    results.to_csv(output_file, index=False)
    print(f"CSV export saved to: {output_file}")


def scalability_test(base_data_generator: Callable,
                     population_counts: List[int] = [10, 20, 30, 50],
                     variant_counts: List[int] = [10000, 50000, 100000],
                     output_dir: str = "scalability_results") -> pd.DataFrame:
    """
    Run scalability tests with varying population and variant counts.
    
    This function is designed for systematic evaluation of how runtime
    and memory scale with data size. Requires a data generator function
    that can create test datasets of specified sizes.
    
    Parameters
    ----------
    base_data_generator : Callable
        Function that generates test data. Should accept (n_populations, n_variants)
        and return paths to generated VCF and panel files.
    population_counts : List[int]
        List of population counts to test
    variant_counts : List[int]
        List of variant counts to test
    output_dir : str
        Directory for output files
        
    Returns
    -------
    pd.DataFrame
        DataFrame with scalability test results
    """
    os.makedirs(output_dir, exist_ok=True)
    all_results = []
    
    for n_pops in population_counts:
        for n_vars in variant_counts:
            print(f"\nTesting: {n_pops} populations, {n_vars} variants")
            
            try:
                # Generate test data
                vcf_file, panel_file = base_data_generator(n_pops, n_vars)
                
                # Run benchmark
                results = run_full_benchmark(
                    vcf_file=vcf_file,
                    panel_file=panel_file,
                    start=0,
                    end=n_vars * 100,  # Approximate genomic range
                    output_prefix=os.path.join(output_dir, f"test_{n_pops}pop_{n_vars}var")
                )
                
                # Add configuration info
                results['n_populations_config'] = n_pops
                results['n_variants_config'] = n_vars
                results['n_pairs_config'] = n_pops * (n_pops - 1)
                
                all_results.append(results)
                
            except Exception as e:
                print(f"  FAILED: {e}")
                continue
    
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(os.path.join(output_dir, "scalability_results.csv"), index=False)
        return combined
    
    return pd.DataFrame()


