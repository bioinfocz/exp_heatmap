"""
Export and reporting module for benchmark results.

This module provides functions to export benchmark results in various formats
suitable for research papers, documentation, and automated analysis:
- Human-readable text reports
- JSON for machine-readable output
- Markdown for documentation
- CSV for data analysis
"""

import sys
import json
from datetime import datetime
from typing import Dict, List, Optional, Any
import pandas as pd
import numpy as np

from exp_heatmap.logging import get_logger

logger = get_logger(__name__)


def generate_benchmark_report(
    results: pd.DataFrame, 
    system_info: Optional['SystemInfo'] = None,
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
    command = " ".join(sys.argv) if hasattr(sys, "argv") else ""
    report_lines = [
        "ExP Heatmap Benchmark Report",
        "=" * 50,
        "",
        f"Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Command: {command}",
        
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
        total_runtime = results['runtime_mean'].sum()
        total_runtime_std = np.sqrt((results['runtime_std'] ** 2).sum()) if 'runtime_std' in results.columns else 0
        max_memory = results['peak_memory_mean'].max()
        max_memory_std = results.loc[results['peak_memory_mean'].idxmax(), 'peak_memory_std'] if 'peak_memory_std' in results.columns else 0
        
        # Get replicate counts per step
        if 'n_replicates' in results.columns:
            replicate_counts = results.set_index('operation')['n_replicates'].to_dict()
            unique_counts = set(replicate_counts.values())
            if len(unique_counts) == 1:
                replicate_str = f"{list(unique_counts)[0]}"
            else:
                replicate_str = ", ".join(f"{op.upper()}: {n}" for op, n in replicate_counts.items())
        else:
            replicate_str = "1"
        
        report_lines.extend([
            f"Number of Replicates: {replicate_str}",
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
        logger.info(f"Report saved to: {output_file}")
    
    return report


def export_benchmark_json(
    results: pd.DataFrame,
    system_info: Optional['SystemInfo'] = None,
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
        logger.info(f"JSON export saved to: {output_file}")
    
    return json_str


def export_benchmark_markdown(
    results: pd.DataFrame,
    system_info: Optional['SystemInfo'] = None,
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
        logger.info(f"Markdown export saved to: {output_file}")
    
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
    logger.info(f"CSV export saved to: {output_file}")

