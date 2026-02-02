import sys
import click
from exp_heatmap import prepare, compute, plot, __version__
from exp_heatmap.logging import setup_logging, get_logger, get_log_file_path

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(help="ExP Heatmap - Population genetics visualization tool", epilog="For more information, see the documentation at:\nhttps://github.com/bioinfocz/exp_heatmap/", context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__, '-v', '--version', prog_name='exp_heatmap')
def cli():
    pass

# prepare command
@cli.command(name='prepare', short_help='Convert VCF file to ZARR format', context_settings=CONTEXT_SETTINGS)
@click.argument('vcf_file', type=click.Path(exists=True, readable=True, dir_okay=False), required=True, metavar='<vcf_file>')
@click.option('-o', '--out', 'output', type=click.Path(), default='zarr_output', show_default=True, help='Directory for ZARR files')
@click.option('--no-log', is_flag=True, help='Disable logging to file')
@click.option('--verbose', is_flag=True, help='Show detailed debug output in console')
def prepare_cmd(vcf_file, output, no_log, verbose):
    """
    <vcf_file>  PATH  Recoded VCF file
    """
    setup_logging("prepare", log_to_file=not no_log, verbose=verbose)
    logger = get_logger(__name__)
    logger.info("Starting prepare command...")
    prepare(vcf_file, output)
    log_file = get_log_file_path()
    if log_file:
        logger.info(f"Log file: {log_file}")

# compute command
@cli.command(name='compute', short_help='Compute population genetics statistics', context_settings=CONTEXT_SETTINGS)
@click.argument('zarr_dir', type=click.Path(exists=True, readable=True, file_okay=False), required=True, metavar='<zarr_dir>')
@click.argument('panel_file', type=click.Path(exists=True, readable=True, dir_okay=False), required=True, metavar='<panel_file>')
@click.option('-o', '--out', 'output', type=click.Path(), default='output', show_default=True, help='Directory for output files')
@click.option('-t', '--test', type=click.Choice(['xpehh', 'xpnsl', 'delta_tajima_d', 'hudson_fst']), default='xpehh', show_default=True, help='Statistical test to compute')
@click.option('-c', '--chunked', is_flag=True, help='Use chunked array to avoid memory exhaustion')
@click.option('--no-log', is_flag=True, help='Disable logging to file')
@click.option('--verbose', is_flag=True, help='Show detailed debug output in console')
def compute_cmd(zarr_dir, panel_file, output, test, chunked, no_log, verbose):
    """
    <zarr_dir>  PATH  Directory with ZARR files from 'exp_heatmap prepare'
    <panel_file>  PATH  Population panel file
    """
    setup_logging("compute", log_to_file=not no_log, verbose=verbose)
    logger = get_logger(__name__)
    logger.info("Starting compute command...")
    compute(zarr_dir, panel_file, output, test, chunked)
    log_file = get_log_file_path()
    if log_file:
        logger.info(f"Log file: {log_file}")

# plot command
@cli.command(name='plot', short_help='Generate heatmap visualization', context_settings=CONTEXT_SETTINGS)
@click.argument('input_dir', type=click.Path(exists=True, file_okay=False), required=True, metavar='<input_dir>')
@click.option('-s', '--start', type=int, help='Start position for the displayed region.')
@click.option('-e', '--end', type=int, help='End position for the displayed region.')
@click.option('-m', '--mid', type=int, help='Middle of the displayed area. The start and end positions will be calculated (mid ± 500 kb)')
@click.option('-t', '--title', type=str, help='Title of the heatmap')
@click.option('-o', '--out', 'output', type=click.Path(), default='ExP_heatmap', show_default=True, help='The output heatmap')
@click.option('-c', '--cmap', type=str, default='Blues', show_default=True, help='Matplotlib colormap for heatmap visualization')
@click.option('--interactive', is_flag=True, help='Generate interactive HTML visualization')
@click.option('--no-log', is_flag=True, help='Disable logging to file')
@click.option('--verbose', is_flag=True, help='Show detailed debug output in console')

def plot_cmd(input_dir, start, end, mid, title, output, cmap, interactive, no_log, verbose):
    """
    <input_dir>  PATH  Directory with TSV files from 'exp_heatmap compute'
    
    For positional arguments, use either [-s with -e] or [-m].
    """
    setup_logging("plot", log_to_file=not no_log, verbose=verbose)
    logger = get_logger(__name__)
    logger.info("Starting plot command...")
    
    start_end_provided = start is not None and end is not None
    mid_provided = mid is not None
    
    # Check for invalid combinations of arguments
    if (start is None) != (end is None):
        raise click.UsageError("--start and --end must be used together")
    if not start_end_provided and not mid_provided:
        raise click.UsageError("Either (--start and --end) or --mid must be provided")
    if start_end_provided and mid_provided:
        raise click.UsageError("Cannot use both (--start and --end) and --mid at the same time")
    
    start_pos, end_pos = (mid - 500000, mid + 500000) if mid_provided else (start, end)
    
    if interactive:
        # Use interactive Plotly visualization
        from exp_heatmap.interactive import plot_interactive
        plot_interactive(
            input_dir,
            start=start_pos,
            end=end_pos,
            title=title,
            output=output,
            colorscale=cmap if cmap != 'expheatmap' else 'Blues'
        )
    else:
        # Use standard matplotlib visualization
        try:
            plot(
                input_dir, 
                start=start_pos, 
                end=end_pos, 
                title=title, 
                output=output, 
                cmap=cmap
            )
        except ValueError as e:
            raise click.ClickException(str(e))
    
    log_file = get_log_file_path()
    if log_file:
        logger.info(f"Log file: {log_file}")

# benchmark command
@cli.command(name='benchmark', short_help='Run performance benchmarks', context_settings=CONTEXT_SETTINGS)
@click.argument('vcf_file', type=click.Path(exists=True, readable=True, dir_okay=False), required=True, metavar='<vcf_file>')
@click.argument('panel_file', type=click.Path(exists=True, readable=True, dir_okay=False), required=True, metavar='<panel_file>')
@click.option('-s', '--start', type=int, required=True, help='Start position for plot benchmark')
@click.option('-e', '--end', type=int, required=True, help='End position for plot benchmark')
@click.option('-o', '--out', 'output', type=click.Path(), default='benchmark', show_default=True, help='Output prefix for benchmark files')
@click.option('-t', '--test', type=click.Choice(['xpehh', 'xpnsl', 'delta_tajima_d', 'hudson_fst']), default='xpehh', show_default=True, help='Statistical test to benchmark')
@click.option('-c', '--chunked', is_flag=True, help='Use chunked array during compute')
@click.option('-r', '--replicates', type=int, default=1, show_default=True, help='Number of replicate runs for statistical analysis')
@click.option('-w', '--warmup', type=int, default=0, show_default=True, help='Number of warmup runs to discard')
@click.option('--report', type=click.Path(), default=None, help='Save benchmark report to file')
@click.option('--no-log', is_flag=True, help='Disable logging to file')
@click.option('--verbose', is_flag=True, help='Show detailed debug output in console')
def benchmark_cmd(vcf_file, panel_file, start, end, output, test, chunked, replicates, warmup, report, no_log, verbose):
    """
    <vcf_file>  PATH  VCF file for benchmarking
    <panel_file>  PATH  Population panel file
    
    Run performance benchmarks on the full ExP Heatmap pipeline.
    
    Use --replicates to run multiple times and get statistical analysis
    (mean, std, 95% confidence intervals).
    """
    setup_logging("benchmark", log_to_file=not no_log, verbose=verbose)
    logger = get_logger(__name__)
    logger.info("Starting benchmark command...")
    
    from exp_heatmap.benchmark import (
        run_full_benchmark, 
        run_full_benchmark_with_replicates,
        generate_benchmark_report
    )
    
    if replicates > 1:
        results, system_info = run_full_benchmark_with_replicates(
            vcf_file=vcf_file,
            panel_file=panel_file,
            start=start,
            end=end,
            output_prefix=output,
            test=test,
            chunked=chunked,
            n_replicates=replicates,
            warmup_runs=warmup
        )
        
        if report:
            generate_benchmark_report(results, system_info=system_info, output_file=report)
        else:
            logger.info("\n" + generate_benchmark_report(results, system_info=system_info))
    else:
        results = run_full_benchmark(
            vcf_file=vcf_file,
            panel_file=panel_file,
            start=start,
            end=end,
            output_prefix=output,
            test=test,
            chunked=chunked
        )
        
        if report:
            generate_benchmark_report(results, output_file=report)
        else:
            logger.info("\n" + generate_benchmark_report(results))
    
    log_file = get_log_file_path()
    if log_file:
        logger.info(f"Log file: {log_file}")


if __name__ == "__main__":
    cli()
