import argparse
import sys

from exp_heatmap import prepare, compute, plot, __version__

def main():
    parser = argparse.ArgumentParser(
        prog="exp_heatmap",
        description="Generate heatmaps from population genetics data using VCF files",
        epilog="For more information, see the documentation at https://github.com/bioinfocz/exp_heatmap/"
    )
    parser.add_argument(
        '-v', '--version', action='version', version=f'%(prog)s {__version__}'
    )
    subparser = parser.add_subparsers()

    parser.set_defaults(func=lambda _: parser.print_help())

    
    # prepare
    prepare_parser = subparser.add_parser("prepare", help="Convert VCF file to ZARR format")
    prepare_parser.add_argument(
        "recode_file", metavar="RECODE_FILE", help="Path to recoded VCF file"
    )
    prepare_parser.add_argument(
        "zarr_dir", metavar="ZARR_DIR", help="Path where ZARR directory will be saved"
    )
    prepare_parser.set_defaults(
        func=lambda args: prepare(args.recode_file, args.zarr_dir)
    )


    # compute
    compute_parser = subparser.add_parser("compute", help="Compute population genetics statistics")
    compute_parser.add_argument(
        "zarr_dir", metavar="ZARR_DIR", help="Path to ZARR directory"
    )
    compute_parser.add_argument(
        "panel_file", metavar="PANEL_FILE", help="Path to panel file"
    )
    compute_parser.add_argument(
        "output_dir", metavar="OUTPUT_DIR", help="Directory for output files"
    )
    compute_parser.add_argument(
        "-t", "--test", choices=['xpehh', 'xpnsl', 'delta_tajima_d', 'hudson_fst'], 
        default="xpehh", help="Statistical test to compute (default: %(default)s)"
    )
    compute_parser.add_argument(
        "-c", "--chunked", action="store_true", help="Use chunked array to avoid memory exhaustion"
    )
    compute_parser.set_defaults(
        func=lambda args: compute(args.zarr_dir, args.panel_file, args.output_dir, args.test, args.chunked)
    )


    # plot
    plot_parser = subparser.add_parser("plot", help="Generate heatmap visualization")
    plot_parser.add_argument(
        "input_dir", metavar="INPUT_DIR", help="Directory containing TSV files from 'exp_heatmap compute'"
    )
    plot_parser.add_argument(
        "--start", type=int, default=None, help="Start position for the displayed region. Uses nearest available position if exact match not found in the input data"
    )
    plot_parser.add_argument(
        "--end", type=int, default=None, help="End position for the displayed region. Uses nearest available position if exact match not found in the input data"
    )
    plot_parser.add_argument(
        "--mid", type=int, default=None, help="Middle of displayed area. The start and end positions will be calculated"
    )
    plot_parser.add_argument(
        "--title", default=None, help="Title of the figure"
    )
    plot_parser.add_argument(
        "--output", default="ExP_heatmap", help="The figure will be saved as OUTPUT + .png (default: %(default)s)"
    )
    plot_parser.add_argument(
        "--cmap", default="Blues", help="Matplotlib colormap for heatmap visualization. Common options: Blues, Reds, Greens, viridis, plasma, inferno, magma, coolwarm, RdYlBu, Spectral, jet (default: %(default)s)"
    )

    def plot_wrapper(args):
        """
        Validate the correct usage of position arguments and call plot function with calculated start/end coordinates.
        """
        has_start_end = args.start is not None and args.end is not None
        has_mid = args.mid is not None
        
        if not has_start_end and not has_mid:
            plot_parser.error("Either (--start and --end) or --mid must be provided")
        elif has_start_end and has_mid:
            plot_parser.error("Cannot use both (--start and --end) and --mid at the same time")
        elif (args.start is None) != (args.end is None):
            plot_parser.error("--start and --end must be used together")
        
        if has_mid:
            start, end = args.mid - 500000, args.mid + 500000
        else:
            start, end = args.start, args.end
            
        plot(args.input_dir, start=start, end=end, title=args.title, output=args.output, cmap=args.cmap)

    plot_parser.set_defaults(func=plot_wrapper)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("")
        sys.exit(1)
