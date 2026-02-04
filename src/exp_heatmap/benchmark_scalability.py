"""
Scalability testing module for benchmarking.

This module provides tools for systematic evaluation of how runtime
and memory scale with data size (population counts, variant counts).
"""

import os
from typing import Callable, List
import pandas as pd

from exp_heatmap.logging import get_logger

logger = get_logger(__name__)


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
        
    Example
    -------
    >>> def generate_test_data(n_pops, n_vars):
    ...     # Generate synthetic VCF and panel files
    ...     vcf_file = f"test_{n_pops}pop_{n_vars}var.vcf"
    ...     panel_file = f"test_{n_pops}pop_panel.tsv"
    ...     # ... generation logic ...
    ...     return vcf_file, panel_file
    >>> 
    >>> results = scalability_test(
    ...     generate_test_data,
    ...     population_counts=[5, 10, 20],
    ...     variant_counts=[1000, 5000, 10000]
    ... )
    """
    from exp_heatmap.benchmark import run_full_benchmark
    
    os.makedirs(output_dir, exist_ok=True)
    all_results = []
    
    for n_pops in population_counts:
        for n_vars in variant_counts:
            logger.debug(f"Testing: {n_pops} populations, {n_vars} variants")
            
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
                logger.error(f"  FAILED: {e}")
                continue
    
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(os.path.join(output_dir, "scalability_results.csv"), index=False)
        return combined
    
    return pd.DataFrame()

