import sys
import click
from exp_heatmap import prepare, compute, plot, __version__

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(help="ExP Heatmap - Population genetics visualization tool", epilog="For more information, see the documentation at:\nhttps://github.com/bioinfocz/exp_heatmap/", context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__, '-v', '--version', prog_name='exp_heatmap')
def cli():
    pass

# prepare command
@cli.command(name='prepare', short_help='Convert VCF file to ZARR format', context_settings=CONTEXT_SETTINGS)
@click.argument('vcf_file', type=click.Path(exists=True, readable=True, dir_okay=False), required=True, metavar='<vcf_file>')
@click.option('-o', '--output', type=click.Path(), default='zarr_output', show_default=True, help='Directory for ZARR files')
def prepare_cmd(vcf_file, output):
    """
    <vcf_file>  PATH  Recoded VCF file
    """
    prepare(vcf_file, output)

# compute command
@cli.command(name='compute', short_help='Compute population genetics statistics', context_settings=CONTEXT_SETTINGS)
@click.argument('zarr_dir', type=click.Path(exists=True, readable=True, file_okay=False), required=True, metavar='<zarr_dir>')
@click.argument('panel_file', type=click.Path(exists=True, readable=True, dir_okay=False), required=True, metavar='<panel_file>')
@click.option('-o', '--output', type=click.Path(), default='output', show_default=True, help='Directory for output files')
@click.option('-t', '--test', type=click.Choice(['xpehh', 'xpnsl', 'delta_tajima_d', 'hudson_fst']), default='xpehh', show_default=True, help='Statistical test to compute')
@click.option('-c', '--chunked', is_flag=True, help='Use chunked array to avoid memory exhaustion')
def compute_cmd(zarr_dir, panel_file, output, test, chunked):
    """
    <zarr_dir>  PATH  Directory with ZARR files from 'exp_heatmap prepare'
    <panel_file>  PATH  Population panel file
    """
    compute(zarr_dir, panel_file, output, test, chunked)

# plot command
@cli.command(name='plot', short_help='Generate heatmap visualization', context_settings=CONTEXT_SETTINGS)
@click.argument('input_dir', type=click.Path(exists=True, file_okay=False), required=True, metavar='<input_dir>')
@click.option('-s', '--start', type=int, help='Start position for the displayed region.')
@click.option('-e', '--end', type=int, help='End position for the displayed region.')
@click.option('-m', '--mid', type=int, help='Middle of the displayed area. The start and end positions will be calculated (mid ± 500 kb)')
@click.option('-t', '--title', type=str, help='Title of the heatmap')
@click.option('-o', '--output', type=click.Path(), default='ExP_heatmap',show_default=True, help='The output heatmap')
@click.option('-c', '--cmap', type=str, default='Blues', show_default=True, help='Matplotlib colormap for heatmap visualization')
def plot_cmd(input_dir, start, end, mid, title, output, cmap):
    """
    <input_dir>  PATH  Directory with TSV files from 'exp_heatmap compute'
    
    For positional arguments, use either [-s with -e] or [-m].
    """
    start_end_provided = start is not None and end is not None
    mid_provided = mid is not None
    
    #Check for invalid combinations of arguments
    if (start is None) != (end is None):
        raise click.UsageError("--start and --end must be used together")
    if not start_end_provided and not mid_provided:
        raise click.UsageError("Either (--start and --end) or --mid must be provided")
    if start_end_provided and mid_provided:
        raise click.UsageError("Cannot use both (--start and --end) and --mid at the same time")
    
    start_pos, end_pos = (mid - 500000, mid + 500000) if mid_provided else (start, end)
    plot(input_dir, start=start_pos, end=end_pos, title=title, output=output, cmap=cmap)

if __name__ == "__main__":
    cli()