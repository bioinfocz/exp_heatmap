import click
from exp_heatmap import __version__
from exp_heatmap.logging import setup_logging, get_logger, get_log_file_path

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
RANK_SCORE_CHOICES = ['directional', '2-tailed', 'ascending', 'descending']


def _finalize_log_file(logger):
    log_file = get_log_file_path()
    if log_file:
        logger.info(f"Log file: {log_file}")

@click.group(help="ExP Heatmap - Population genetics visualization tool", epilog="For more information, see the documentation at:\nhttps://github.com/bioinfocz/exp_heatmap/", context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__, '-v', '--version', prog_name='exp_heatmap')
def cli():
    pass


@cli.command(name='filter-vcf', short_help='Keep only biallelic SNP records', context_settings=CONTEXT_SETTINGS)
@click.argument('input_vcf', type=click.Path(exists=True, readable=True, dir_okay=False), required=True, metavar='<input_vcf>')
@click.option('-o', '--out', 'output', type=click.Path(), required=True, help='Output VCF path (.vcf or .vcf.gz)')
@click.option('--chrom', type=str, help='Optional chromosome name for region-scoped filtering')
@click.option('--start', type=int, help='Optional region start position (inclusive)')
@click.option('--end', type=int, help='Optional region end position (inclusive)')
@click.option('--no-log', is_flag=True, help='Disable logging to file')
@click.option('--verbose', is_flag=True, help='Show detailed debug output in console')
def filter_vcf_cmd(input_vcf, output, chrom, start, end, no_log, verbose):
    """
    <input_vcf>  PATH  Input VCF or VCF.GZ file
    """
    from exp_heatmap.vcf_utils import filter_biallelic_snp_vcf

    setup_logging("filter_vcf", log_to_file=not no_log, verbose=verbose)
    logger = get_logger(__name__)
    logger.info("Starting filter-vcf command...")

    if (start is None) != (end is None):
        raise click.UsageError("--start and --end must be provided together")
    if start is not None and end is not None and start > end:
        raise click.UsageError("--start must be less than or equal to --end")

    stats = filter_biallelic_snp_vcf(
        input_vcf,
        output,
        chromosome=chrom,
        start=start,
        end=end,
    )

    region_text = ""
    if chrom is not None or start is not None or end is not None:
        region_bits = []
        if chrom is not None:
            region_bits.append(str(chrom))
        if start is not None and end is not None:
            region_bits.append(f"{start}-{end}")
        region_text = f" for region {'/'.join(region_bits)}"

    logger.info(
        "Filtered VCF written to %s%s (%s kept / %s total records)",
        output,
        region_text,
        stats["kept_records"],
        stats["total_records"],
    )
    _finalize_log_file(logger)

@cli.command(name='full', short_help='Run the full pipeline', context_settings=CONTEXT_SETTINGS)
@click.argument('vcf_file', type=click.Path(exists=True, readable=True, dir_okay=False), required=True, metavar='<vcf_file>')
@click.argument('panel_file', type=click.Path(exists=True, readable=True, dir_okay=False), required=True, metavar='<panel_file>')
@click.option('-o', '--out', 'output', type=click.Path(), default='exp_heatmap', show_default=True, help='Prefix for all output files')
@click.option('-s', '--start', type=int, required=True, help='Start position for the displayed region.')
@click.option('-e', '--end', type=int, required=True, help='End position for the displayed region.')
@click.option('-t', '--test', type=click.Choice(['xpehh', 'xpnsl', 'delta_tajima_d', 'hudson_fst']), default='xpehh', show_default=True, help='Statistical test to compute')
@click.option('-c', '--chunked', is_flag=True, help='Use chunked array to avoid memory exhaustion')
@click.option('--title', type=str, help='Title of the heatmap')
@click.option('--cmap', type=str, default='Blues', show_default=True, help='Matplotlib colormap for heatmap visualization')
@click.option('--interactive', is_flag=True, help='Generate interactive HTML visualization')
@click.option('--rank-scores', type=click.Choice(RANK_SCORE_CHOICES), default='directional', show_default=True, help='Rank-score mode to visualize')
@click.option('--max-columns', type=int, help='Maximum rendered columns for wide regions; enables explicit downsampling for static or interactive output')
@click.option('--column-aggregation', type=click.Choice(['max', 'mean', 'median']), default='max', show_default=True, help='Aggregation used when static output is downsampled')
@click.option('--dpi', type=int, default=400, show_default=True, help='DPI for saved static figures')
@click.option('--no-log', is_flag=True, help='Disable logging to file')
@click.option('--verbose', is_flag=True, help='Show detailed debug output in console')
def full_cmd(vcf_file, panel_file, output, start, end, test, chunked, title, cmap, interactive, rank_scores, max_columns, column_aggregation, dpi, no_log, verbose):
    from exp_heatmap.compute import compute
    from exp_heatmap.interactive import plot_interactive
    from exp_heatmap.plot import plot
    from exp_heatmap.prepare import prepare

    setup_logging("full", log_to_file=not no_log, verbose=verbose)
    logger = get_logger(__name__)
    logger.info("Starting full command...")
    
    # Prepare
    zarr_dir = f"{output}_zarr"
    prepare(vcf_file, zarr_dir)
    
    # Compute
    compute_dir = f"{output}_compute"
    compute(zarr_dir, panel_file, compute_dir, test, chunked)
    
    # Plot
    plot_output = f"{output}_plot"
    if interactive:
        # Use interactive Plotly visualization
        colorscale=cmap if cmap != 'expheatmap' else 'Blues'
        plot_interactive(
            compute_dir,
            start,
            end,
            title,
            plot_output,
            rank_scores=rank_scores,
            colorscale=colorscale,
            max_columns=max_columns,
        )
    else:
        # Use static matplotlib visualization
        plot(
            compute_dir,
            start,
            end,
            title,
            plot_output,
            cmap,
            max_columns=max_columns,
            column_aggregation=column_aggregation,
            dpi=dpi,
            rank_scores=rank_scores,
        )
    
    _finalize_log_file(logger)

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
    from exp_heatmap.prepare import prepare

    setup_logging("prepare", log_to_file=not no_log, verbose=verbose)
    logger = get_logger(__name__)
    logger.info("Starting prepare command...")
    
    prepare(vcf_file, output)
    _finalize_log_file(logger)

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
    from exp_heatmap.compute import compute

    setup_logging("compute", log_to_file=not no_log, verbose=verbose)
    logger = get_logger(__name__)
    logger.info("Starting compute command...")
    
    compute(zarr_dir, panel_file, output, test, chunked)
    _finalize_log_file(logger)

# plot command
@cli.command(name='plot', short_help='Generate heatmap visualization', context_settings=CONTEXT_SETTINGS)
@click.argument('input_dir', type=click.Path(exists=True, file_okay=False), required=True, metavar='<input_dir>')
@click.option('-s', '--start', type=int, required=True, help='Start position for the displayed region.')
@click.option('-e', '--end', type=int, required=True, help='End position for the displayed region.')
@click.option('-t', '--title', type=str, help='Title of the heatmap')
@click.option('-o', '--out', 'output', type=click.Path(), default='ExP_heatmap', show_default=True, help='The output heatmap')
@click.option('-c', '--cmap', type=str, default='Blues', show_default=True, help='Matplotlib colormap for heatmap visualization')
@click.option('--interactive', is_flag=True, help='Generate interactive HTML visualization')
@click.option('--rank-scores', type=click.Choice(RANK_SCORE_CHOICES), default='directional', show_default=True, help='Rank-score mode to visualize')
@click.option('--max-columns', type=int, help='Maximum rendered columns for wide regions; enables explicit downsampling for static or interactive output')
@click.option('--column-aggregation', type=click.Choice(['max', 'mean', 'median']), default='max', show_default=True, help='Aggregation used when static output is downsampled')
@click.option('--dpi', type=int, default=400, show_default=True, help='DPI for saved static figures')
@click.option('--no-log', is_flag=True, help='Disable logging to file')
@click.option('--verbose', is_flag=True, help='Show detailed debug output in console')

def plot_cmd(input_dir, start, end, title, output, cmap, interactive, rank_scores, max_columns, column_aggregation, dpi, no_log, verbose):
    from exp_heatmap.interactive import plot_interactive
    from exp_heatmap.plot import plot

    """
    <input_dir>  PATH  Directory with TSV files from 'exp_heatmap compute'
    """
    setup_logging("plot", log_to_file=not no_log, verbose=verbose)
    logger = get_logger(__name__)
    logger.info("Starting plot command...")
    
    if interactive:
        # Use interactive Plotly visualization
        colorscale=cmap if cmap != 'expheatmap' else 'Blues'
        plot_interactive(
            input_dir,
            start,
            end,
            title,
            output,
            rank_scores=rank_scores,
            colorscale=colorscale,
            max_columns=max_columns,
        )
    else:
        # Use static matplotlib visualization
        plot(
            input_dir,
            start,
            end,
            title,
            output,
            cmap,
            max_columns=max_columns,
            column_aggregation=column_aggregation,
            dpi=dpi,
            rank_scores=rank_scores,
        )
    
    _finalize_log_file(logger)


@cli.command(name='summary', short_help='Summarize to superpopulation pairs', context_settings=CONTEXT_SETTINGS)
@click.argument('input_dir', type=click.Path(exists=True, file_okay=False), required=True, metavar='<input_dir>')
@click.option('-s', '--start', type=int, required=True, help='Start position for the displayed region.')
@click.option('-e', '--end', type=int, required=True, help='End position for the displayed region.')
@click.option('-o', '--out', 'output', type=click.Path(), default='ExP_superpopulation_summary', show_default=True, help='Output prefix for summary files')
@click.option('-t', '--title', type=str, help='Title for the summary view')
@click.option('-c', '--cmap', type=str, default='Blues', show_default=True, help='Plotly colorscale for the interactive summary heatmap')
@click.option('--rank-scores', type=click.Choice(RANK_SCORE_CHOICES), default='directional', show_default=True, help='Rank-score mode to visualize')
@click.option('--agg-func', type=click.Choice(['mean', 'median', 'max', 'min']), default='mean', show_default=True, help='Aggregation used to collapse population pairs to superpopulation pairs')
@click.option('--max-columns', type=int, help='Maximum rendered columns for the interactive summary heatmap')
@click.option('--no-plot', is_flag=True, help='Only write the collapsed TSV without generating an HTML view')
@click.option('--no-log', is_flag=True, help='Disable logging to file')
@click.option('--verbose', is_flag=True, help='Show detailed debug output in console')
def summary_cmd(input_dir, start, end, output, title, cmap, rank_scores, agg_func, max_columns, no_plot, no_log, verbose):
    from exp_heatmap.plot import create_plot_input, summarize_by_superpopulation
    from exp_heatmap.interactive import plot_interactive_heatmap

    """
    <input_dir>  PATH  Directory with TSV files from 'exp_heatmap compute'
    """
    setup_logging("summary", log_to_file=not no_log, verbose=verbose)
    logger = get_logger(__name__)
    logger.info("Starting summary command...")

    plot_input = create_plot_input(input_dir, start=start, end=end, rank_scores=rank_scores)
    summary_df = summarize_by_superpopulation(plot_input, agg_func=agg_func)
    summary_path = f"{output}.tsv"
    summary_df.to_csv(summary_path, sep="\t")
    logger.info(f"Superpopulation summary written to {summary_path}")

    if not no_plot:
        plot_interactive_heatmap(
            summary_df,
            start=int(summary_df.columns[0]),
            end=int(summary_df.columns[-1]),
            title=title or f"Superpopulation summary ({agg_func})",
            output=output,
            colorscale=cmap,
            max_columns=max_columns,
            show_superpop_annotations=False,
        )

    _finalize_log_file(logger)


@cli.command(name='regions', short_help='Extract top-scoring genomic regions', context_settings=CONTEXT_SETTINGS)
@click.argument('input_dir', type=click.Path(exists=True, file_okay=False), required=True, metavar='<input_dir>')
@click.option('-s', '--start', type=int, required=True, help='Start position for the scanned region.')
@click.option('-e', '--end', type=int, required=True, help='End position for the scanned region.')
@click.option('-o', '--out', 'output', type=click.Path(), default='ExP_top_regions.tsv', show_default=True, help='Output TSV path')
@click.option('--rank-scores', type=click.Choice(RANK_SCORE_CHOICES), default='directional', show_default=True, help='Rank-score mode to use before region extraction')
@click.option('--n-top', type=int, default=100, show_default=True, help='Number of top regions to report')
@click.option('--window-size', type=int, default=10000, show_default=True, help='Window size in base pairs')
@click.option('--min-gap', type=int, default=5000, show_default=True, help='Minimum spacing between reported centers')
@click.option('--no-log', is_flag=True, help='Disable logging to file')
@click.option('--verbose', is_flag=True, help='Show detailed debug output in console')
def regions_cmd(input_dir, start, end, output, rank_scores, n_top, window_size, min_gap, no_log, verbose):
    from exp_heatmap.plot import create_plot_input, extract_top_regions

    """
    <input_dir>  PATH  Directory with TSV files from 'exp_heatmap compute'
    """
    setup_logging("regions", log_to_file=not no_log, verbose=verbose)
    logger = get_logger(__name__)
    logger.info("Starting regions command...")

    plot_input = create_plot_input(input_dir, start=start, end=end, rank_scores=rank_scores)
    top_regions = extract_top_regions(
        plot_input,
        n_top=n_top,
        window_size=window_size,
        min_gap=min_gap,
    )
    top_regions.to_csv(output, sep="\t", index=False)
    logger.info(f"Top regions written to {output}")
    _finalize_log_file(logger)


@cli.command(name='focus', short_help='Generate a population-focused interactive view', context_settings=CONTEXT_SETTINGS)
@click.argument('input_dir', type=click.Path(exists=True, file_okay=False), required=True, metavar='<input_dir>')
@click.option('-s', '--start', type=int, required=True, help='Start position for the displayed region.')
@click.option('-e', '--end', type=int, required=True, help='End position for the displayed region.')
@click.option('-p', '--population', 'focus_population', type=str, required=True, help='Population code to focus on')
@click.option('-o', '--out', 'output', type=click.Path(), default='ExP_population_focus', show_default=True, help='Output HTML prefix')
@click.option('-t', '--title', type=str, help='Title of the focus view')
@click.option('-c', '--cmap', type=str, default='Blues', show_default=True, help='Plotly colorscale for the interactive heatmap')
@click.option('--rank-scores', type=click.Choice(RANK_SCORE_CHOICES), default='directional', show_default=True, help='Rank-score mode to visualize')
@click.option('--max-columns', type=int, help='Maximum rendered columns for wide focus views')
@click.option('--no-log', is_flag=True, help='Disable logging to file')
@click.option('--verbose', is_flag=True, help='Show detailed debug output in console')
def focus_cmd(input_dir, start, end, focus_population, output, title, cmap, rank_scores, max_columns, no_log, verbose):
    from exp_heatmap.interactive import create_population_focus_view
    from exp_heatmap.plot import create_plot_input

    """
    <input_dir>  PATH  Directory with TSV files from 'exp_heatmap compute'
    """
    setup_logging("focus", log_to_file=not no_log, verbose=verbose)
    logger = get_logger(__name__)
    logger.info("Starting focus command...")

    plot_input = create_plot_input(input_dir, start=start, end=end, rank_scores=rank_scores)
    create_population_focus_view(
        plot_input,
        focus_population=focus_population,
        start=int(plot_input.columns[0]),
        end=int(plot_input.columns[-1]),
        title=title,
        output=output,
        colorscale=cmap,
        max_columns=max_columns,
    )

    _finalize_log_file(logger)


@cli.command(name='compare', short_help='Generate an interactive side-by-side region comparison', context_settings=CONTEXT_SETTINGS)
@click.argument('input_dir', type=click.Path(exists=True, file_okay=False), required=True, metavar='<input_dir>')
@click.option('--start-1', type=int, required=True, help='Start position of the first region')
@click.option('--end-1', type=int, required=True, help='End position of the first region')
@click.option('--start-2', type=int, required=True, help='Start position of the second region')
@click.option('--end-2', type=int, required=True, help='End position of the second region')
@click.option('-o', '--out', 'output', type=click.Path(), default='ExP_comparison', show_default=True, help='Output HTML prefix')
@click.option('-t', '--title', type=str, help='Title of the comparison view')
@click.option('-c', '--cmap', type=str, default='Blues', show_default=True, help='Plotly colorscale for the interactive heatmap')
@click.option('--rank-scores', type=click.Choice(RANK_SCORE_CHOICES), default='directional', show_default=True, help='Rank-score mode to visualize')
@click.option('--max-columns', type=int, help='Maximum rendered columns per region in the comparison view')
@click.option('--no-log', is_flag=True, help='Disable logging to file')
@click.option('--verbose', is_flag=True, help='Show detailed debug output in console')
def compare_cmd(input_dir, start_1, end_1, start_2, end_2, output, title, cmap, rank_scores, max_columns, no_log, verbose):
    from exp_heatmap.interactive import create_comparison_view
    from exp_heatmap.plot import create_plot_input

    """
    <input_dir>  PATH  Directory with TSV files from 'exp_heatmap compute'
    """
    setup_logging("compare", log_to_file=not no_log, verbose=verbose)
    logger = get_logger(__name__)
    logger.info("Starting compare command...")

    start = min(start_1, start_2)
    end = max(end_1, end_2)
    plot_input = create_plot_input(input_dir, start=start, end=end, rank_scores=rank_scores)
    create_comparison_view(
        plot_input,
        region1=(start_1, end_1),
        region2=(start_2, end_2),
        title=title,
        output=output,
        colorscale=cmap,
        max_columns=max_columns,
    )

    _finalize_log_file(logger)

if __name__ == "__main__":
    cli()
