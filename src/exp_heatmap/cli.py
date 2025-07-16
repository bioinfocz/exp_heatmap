import argparse

from exp_heatmap.prepare import prepare
from exp_heatmap.compute import compute
from exp_heatmap.plot import plot


def main():
    parser = argparse.ArgumentParser(
        prog="exp_heatmap",
        description="Generate heatmaps from population genetics data using VCF files",
        epilog="For more information, see the documentation at https://github.com/bioinfocz/exp_heatmap/"
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
    compute_parser.set_defaults(
        func=lambda args: compute(args.zarr_dir, args.panel_file, args.output_dir, args.test)
    )


    # plot
    plot_parser = subparser.add_parser("plot", help="Generate heatmap visualization")
    plot_parser.add_argument(
        "input_dir", metavar="INPUT_DIR", help="Directory containing TSV files from 'exp_heatmap compute'"
    )
    plot_parser.add_argument("--begin", type=int, required=True,
        help="Beginning of displayed area (if this position is not in your input data, will start from the next closest position instead)")
    plot_parser.add_argument("--end", type=int, required=True,
        help="End of displayed area (if this position is not in your input data, will finish on the previous closest position instead)")
    plot_parser.add_argument("--title", default=None, help="Title of the figure")
    plot_parser.add_argument("--output", default="ExP_heatmap", help="The figure will be saved as OUTPUT + .png (default: %(default)s)")
    plot_parser.add_argument("--cmap", default="Blues", help="Matplotlib colormap for heatmap visualization. Common options: Blues, Reds, Greens, viridis, plasma, inferno, magma, coolwarm, RdYlBu, Spectral, jet (default: %(default)s)")

    plot_parser.set_defaults(
        func=lambda args: plot(
            args.input_dir,
            begin=args.begin,
            end=args.end,
            title=args.title,
            output=args.output,
            cmap=args.cmap,
        )
    )

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
