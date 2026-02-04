import click
from exp_heatmap import prepare, compute, plot, __version__
from exp_heatmap.logging import setup_logging, get_logger, get_log_file_path
from exp_heatmap.interactive import plot_interactive

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(help="ExP Heatmap - Population genetics visualization tool", epilog="For more information, see the documentation at:\nhttps://github.com/bioinfocz/exp_heatmap/", context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__, '-v', '--version', prog_name='exp_heatmap')
def cli():
    pass

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
@click.option('--no-log', is_flag=True, help='Disable logging to file')
@click.option('--verbose', is_flag=True, help='Show detailed debug output in console')
def full_cmd(vcf_file, panel_file, output, start, end, test, chunked, title, cmap, interactive, no_log, verbose):
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
        plot_interactive(compute_dir,start,end,title,plot_output,colorscale)
    else:
        # Use static matplotlib visualization
        plot(compute_dir, start, end, title, plot_output, cmap)
    
    log_file = get_log_file_path()
    if log_file:
        logger.info(f"Log file: {log_file}")

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
@click.option('-s', '--start', type=int, required=True, help='Start position for the displayed region.')
@click.option('-e', '--end', type=int, required=True, help='End position for the displayed region.')
@click.option('-t', '--title', type=str, help='Title of the heatmap')
@click.option('-o', '--out', 'output', type=click.Path(), default='ExP_heatmap', show_default=True, help='The output heatmap')
@click.option('-c', '--cmap', type=str, default='Blues', show_default=True, help='Matplotlib colormap for heatmap visualization')
@click.option('--interactive', is_flag=True, help='Generate interactive HTML visualization')
@click.option('--no-log', is_flag=True, help='Disable logging to file')
@click.option('--verbose', is_flag=True, help='Show detailed debug output in console')

def plot_cmd(input_dir, start, end, title, output, cmap, interactive, no_log, verbose):
    """
    <input_dir>  PATH  Directory with TSV files from 'exp_heatmap compute'
    """
    setup_logging("plot", log_to_file=not no_log, verbose=verbose)
    logger = get_logger(__name__)
    logger.info("Starting plot command...")
    
    if interactive:
        # Use interactive Plotly visualization
        colorscale=cmap if cmap != 'expheatmap' else 'Blues'
        plot_interactive(input_dir,start,end,title,output,colorscale)
    else:
        # Use static matplotlib visualization
        plot(input_dir, start, end, title, output, cmap)
    
    log_file = get_log_file_path()
    if log_file:
        logger.info(f"Log file: {log_file}")

if __name__ == "__main__":
    cli()
