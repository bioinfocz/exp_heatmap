import sys
import click
from exp_heatmap import prepare, compute, plot, __version__

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(help="ExP Heatmap - Population genetics visualization tool", epilog="For more information, see the documentation at:\nhttps://github.com/bioinfocz/exp_heatmap/", context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__, '-v', '--version', prog_name='exp_heatmap')
def cli():
    pass

# prepare command
@cli.command(name='prepare', context_settings=CONTEXT_SETTINGS)
@click.argument('vcf_file', type=click.Path(exists=True, readable=True, dir_okay=False), required=True, metavar='<vcf_file>')
@click.option('-o', '--output', type=click.Path(), default='zarr_output', help='Directory for ZARR files')
def prepare_cmd(vcf_file, output):
    """
    Convert VCF file to ZARR format.
    
    <vcf_file>  PATH  Recoded VCF file
    """
    prepare(vcf_file, output)

# compute command
@cli.command(name='compute',context_settings=CONTEXT_SETTINGS)
@click.argument('zarr_dir', type=click.Path(exists=True, readable=True, file_okay=False), required=True, metavar='<zarr_dir>')
@click.argument('panel_file', type=click.Path(exists=True, readable=True, dir_okay=False), required=True, metavar='<panel_file>')
@click.option('-o', '--output', type=click.Path(), default='output', help='Directory for output files')
@click.option('-t', '--test', type=click.Choice(['xpehh', 'xpnsl', 'delta_tajima_d', 'hudson_fst'], case_sensitive=False), default='xpehh', show_default=True, help='Statistical test to compute')
@click.option('-c', '--chunked', is_flag=True, help='Use chunked array to avoid memory exhaustion')
def compute_cmd(zarr_dir, panel_file, output, test, chunked):
    """
    Compute population genetics statistics.
    
    \b
    <zarr_dir>  PATH  Directory with ZARR files from 'exp_heatmap prepare'
    <panel_file>  PATH  Population panel file
    """
    compute(zarr_dir, panel_file, output, test, chunked)

# plot command
@cli.command(name='plot', context_settings=CONTEXT_SETTINGS)
@click.argument('input_dir', type=click.Path(exists=True, file_okay=False), required=True, metavar='<input_dir>')
@click.option('-s', '--start', type=int, help='Start position for the displayed region. Uses nearest available position if exact match not found in the input data [USE WITH -e]')
@click.option('-e', '--end', type=int, help='End position for the displayed region. Uses nearest available position if exact match not found in the input data [USE WITH -s]')
@click.option('-m', '--mid', type=int, help='Middle of the displayed area. The start and end positions will be calculated (mid ± 500 kb)')
@click.option('-t', '--title', type=str, help='Title of the heatmap')
@click.option('-o', '--output', type=click.Path(), default='ExP_heatmap',show_default=True, help='The output heatmap')
@click.option('-c', '--cmap', type=str, default='Blues', show_default=True, help='Matplotlib colormap for heatmap visualization')
def plot_cmd(input_dir, start, end, mid, title, output, cmap):
    """
    Generate heatmap visualization.
    
    <input_dir>  PATH  Directory with TSV files from 'exp_heatmap compute'
    """
    # TODO: Reduce the size of the validation step
    has_start_end = start is not None and end is not None
    has_mid = mid is not None
    
    if not has_start_end and not has_mid:
        raise click.UsageError("Either (--start and --end) or --mid must be provided")
    elif has_start_end and has_mid:
        raise click.UsageError("Cannot use both (--start and --end) and --mid at the same time")
    elif (start is None) != (end is None):
        raise click.UsageError("--start and --end must be used together")
    
    # Calculate start and end positions
    if has_mid:
        start_pos, end_pos = mid - 500000, mid + 500000
    else:
        start_pos, end_pos = start, end
    
    plot(input_dir, start=start_pos, end=end_pos, title=title, output=output, cmap=cmap)

if __name__ == "__main__":
    cli()