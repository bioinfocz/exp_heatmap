"""
Benchmarking module for ExP Heatmap performance evaluation.

This module provides tools to measure runtime and memory usage across
the prepare, compute, and plot steps of the ExP Heatmap pipeline.
Results can be used to evaluate scalability and identify bottlenecks.
"""

import os
import time
import tracemalloc
import functools
from typing import Dict, List, Optional, Callable, Any
from dataclasses import dataclass, field
import pandas as pd
import numpy as np


@dataclass
class BenchmarkResult:
    """Container for benchmark results from a single operation."""
    operation: str
    runtime_seconds: float
    peak_memory_mb: float
    current_memory_mb: float
    parameters: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for DataFrame creation."""
        return {
            'operation': self.operation,
            'runtime_seconds': round(self.runtime_seconds, 3),
            'peak_memory_mb': round(self.peak_memory_mb, 2),
            'current_memory_mb': round(self.current_memory_mb, 2),
            **self.parameters
        }


class Benchmarker:
    """
    Context manager and utility class for benchmarking ExP Heatmap operations.
    
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
        self.result: Optional[BenchmarkResult] = None
        
    def __enter__(self):
        """Start benchmarking."""
        tracemalloc.start()
        self.start_time = time.perf_counter()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Stop benchmarking and record results."""
        end_time = time.perf_counter()
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        self.result = BenchmarkResult(
            operation=self.operation_name,
            runtime_seconds=end_time - self.start_time,
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


def benchmark_prepare(vcf_file: str, zarr_dir: str) -> BenchmarkResult:
    """
    Benchmark the prepare step (VCF to ZARR conversion).
    
    Parameters
    ----------
    vcf_file : str
        Path to input VCF file
    zarr_dir : str
        Path for output ZARR directory
        
    Returns
    -------
    BenchmarkResult
        Benchmark results including runtime and memory usage
    """
    from exp_heatmap import prepare as prepare_func
    
    # Get file size for context
    file_size_mb = os.path.getsize(vcf_file) / (1024 * 1024)
    
    with Benchmarker("prepare", vcf_file_size_mb=file_size_mb) as bench:
        prepare_func(vcf_file, zarr_dir)
    
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
        Benchmark results including runtime and memory usage
    """
    from exp_heatmap import compute as compute_func
    
    # Count populations from panel file
    panel = pd.read_csv(panel_file, sep="\t")
    n_populations = panel['pop'].nunique()
    n_samples = len(panel)
    n_pairs = n_populations * (n_populations - 1) // 2
    
    with Benchmarker("compute", 
                     test=test,
                     n_populations=n_populations,
                     n_samples=n_samples,
                     n_pairs=n_pairs,
                     chunked=chunked) as bench:
        compute_func(zarr_dir, panel_file, output_dir, test, chunked)
    
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
        Benchmark results including runtime and memory usage
    """
    from exp_heatmap.plot import plot
    
    # Count files and estimate data size
    import glob
    tsv_files = glob.glob(os.path.join(input_dir, "*.tsv"))
    n_files = len(tsv_files)
    region_size = end - start
    
    with Benchmarker("plot",
                     n_files=n_files,
                     region_size_bp=region_size,
                     start=start,
                     end=end) as bench:
        plot(input_dir, start=start, end=end, title="Benchmark", output=output, **plot_kwargs)
    
    return bench.result


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


def generate_benchmark_report(results: pd.DataFrame, output_file: str = None) -> str:
    """
    Generate a formatted benchmark report.
    
    Parameters
    ----------
    results : pd.DataFrame
        DataFrame with benchmark results
    output_file : str, optional
        If provided, save report to this file
        
    Returns
    -------
    str
        Formatted benchmark report
    """
    report_lines = [
        "ExP Heatmap Benchmark Report",
        "=" * 50,
        "",
        "Summary",
        "-" * 30,
    ]
    
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
            f"  Peak Memory: {row['peak_memory_mb']:.1f} MB",
        ])
        
        # Add operation-specific parameters
        for key, value in row.items():
            if key not in ['operation', 'runtime_seconds', 'peak_memory_mb', 'current_memory_mb']:
                if pd.notna(value):
                    report_lines.append(f"  {key}: {value}")
    
    report = "\n".join(report_lines)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report)
        print(f"Report saved to: {output_file}")
    
    return report


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


