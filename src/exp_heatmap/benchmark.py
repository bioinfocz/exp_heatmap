"""
Core benchmarking module for ExP Heatmap performance evaluation.

This module provides the essential tools to measure runtime and memory usage
across the prepare, compute, and plot steps of the ExP Heatmap pipeline.

For statistical analysis with replicates, see `benchmark_stats.py`
For export/reporting functions, see `benchmark_export.py`
For scalability testing, see `benchmark_scalability.py`

Basic Usage
-----------
    from exp_heatmap.benchmark import (
        benchmark_prepare, benchmark_compute, benchmark_plot,
        run_full_benchmark
    )
    
    # Run full pipeline benchmark
    results = run_full_benchmark(
        vcf_file="data.vcf",
        panel_file="panel.tsv",
        start=1000000,
        end=2000000
    )
    
    # Or benchmark individual steps
    result = benchmark_prepare("data.vcf", "data.zarr")
    print(f"Prepare took {result.runtime_seconds:.2f}s")

For Advanced Usage
------------------
    # Multi-run statistical benchmarking
    from exp_heatmap.benchmark_stats import run_full_benchmark_with_replicates
    
    results_df, system_info = run_full_benchmark_with_replicates(
        vcf_file="data.vcf",
        panel_file="panel.tsv",
        start=1000000,
        end=2000000,
        n_replicates=5
    )
    
    # Export results
    from exp_heatmap.benchmark_export import (
        generate_benchmark_report,
        export_benchmark_markdown
    )
    
    report = generate_benchmark_report(results_df, system_info)
    export_benchmark_markdown(results_df, system_info, "benchmark_report.md")
"""

import os
import time
import platform
import resource
import hashlib
import tracemalloc
import functools
from datetime import datetime
from typing import Dict, List, Optional, Callable, Any
from dataclasses import dataclass, field, asdict
import pandas as pd

from exp_heatmap.logging import get_logger

logger = get_logger(__name__)


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


def run_full_benchmark(vcf_file: str, panel_file: str, 
                       start: int, end: int,
                       output_prefix: str = "benchmark",
                       test: str = "xpehh",
                       chunked: bool = False) -> pd.DataFrame:
    """
    Run a complete benchmark of the ExP Heatmap pipeline.
    
    This is a simple single-run benchmark. For multiple replicates and
    statistical analysis, use `run_full_benchmark_with_replicates` from
    the `benchmark_stats` module.
    
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
    
    logger.info("=" * 60)
    logger.info("ExP Heatmap Full Pipeline Benchmark")
    logger.info("=" * 60)
    
    # Step 1: Prepare
    logger.info("[1/3] Benchmarking PREPARE step...")
    try:
        result = benchmark_prepare(vcf_file, zarr_dir)
        results.append(result.to_dict())
        logger.debug(f"      Runtime: {result.runtime_seconds:.2f}s, Peak Memory: {result.peak_memory_mb:.1f}MB")
    except Exception as e:
        logger.error(f"      FAILED: {e}")
    
    # Step 2: Compute
    logger.info("[2/3] Benchmarking COMPUTE step...")
    try:
        result = benchmark_compute(zarr_dir, panel_file, compute_dir, test, chunked)
        results.append(result.to_dict())
        logger.debug(f"      Runtime: {result.runtime_seconds:.2f}s, Peak Memory: {result.peak_memory_mb:.1f}MB")
    except Exception as e:
        logger.error(f"      FAILED: {e}")
    
    # Step 3: Plot
    logger.info("[3/3] Benchmarking PLOT step...")
    try:
        result = benchmark_plot(compute_dir, start, end, plot_output)
        results.append(result.to_dict())
        logger.debug(f"      Runtime: {result.runtime_seconds:.2f}s, Peak Memory: {result.peak_memory_mb:.1f}MB")
    except Exception as e:
        logger.error(f"      FAILED: {e}")
    
    logger.info("=" * 60)
    logger.info("Benchmark Complete")
    logger.info("=" * 60)
    
    return pd.DataFrame(results)


# =============================================================================
# Re-export advanced functionality from submodules for backward compatibility
# =============================================================================

# Import statistical analysis functions
from exp_heatmap.benchmark_stats import (
    BenchmarkStatistics,
    compute_statistics,
    run_benchmark_with_replicates,
    run_full_benchmark_with_replicates,
)

# Import export/report functions
from exp_heatmap.benchmark_export import (
    generate_benchmark_report,
    export_benchmark_json,
    export_benchmark_markdown,
    export_benchmark_csv,
)

# Import scalability testing
from exp_heatmap.benchmark_scalability import (
    scalability_test,
)

__all__ = [
    # Core classes
    'SystemInfo',
    'BenchmarkResult',
    'Benchmarker',
    # Core benchmark functions
    'benchmark_prepare',
    'benchmark_compute',
    'benchmark_plot',
    'run_full_benchmark',
    # Statistical analysis (from benchmark_stats)
    'BenchmarkStatistics',
    'compute_statistics',
    'run_benchmark_with_replicates',
    'run_full_benchmark_with_replicates',
    # Export/reporting (from benchmark_export)
    'generate_benchmark_report',
    'export_benchmark_json',
    'export_benchmark_markdown',
    'export_benchmark_csv',
    # Scalability testing (from benchmark_scalability)
    'scalability_test',
]
