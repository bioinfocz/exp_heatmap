import argparse

from exp_selection.prepare import prepare
from exp_selection.compute import compute
from exp_selection.plot import plot


def main():
    parser = argparse.ArgumentParser(
        prog="exp-selection",
    )
    subparser = parser.add_subparsers()

    parser.set_defaults(func=lambda _: parser.print_help())

    prepare_parser = subparser.add_parser("prepare")
    prepare_parser.add_argument(
        "recode_file", metavar="RECODE_FILE", help="Where is recoded VCF"
    )
    prepare_parser.add_argument(
        "zarr_dir", metavar="ZARR_DIR", help="Where ZARR dir will be saved"
    )
    prepare_parser.set_defaults(
        func=lambda args: prepare(args.recode_file, args.zarr_dir)
    )

    compute_parser = subparser.add_parser("compute")
    compute_parser.add_argument(
        "zarr_dir", metavar="ZARR_DIR", help="Where ZARR dir is located"
    )
    compute_parser.add_argument(
        "panel_file", metavar="PANEL_FILE", help="Where panel file will be saved"
    )
    compute_parser.add_argument(
        "output_dir", metavar="OUTPUT_DIR", help="Where output dir will be saved"
    )
    compute_parser.set_defaults(
        func=lambda args: compute(args.zarr_dir, args.panel_file, args.output_dir)
    )

    plot_parser = subparser.add_parser("plot")
    plot_parser.add_argument(
        "output_dir", metavar="OUTPUT_DIR", help="Where output dir is located"
    )
    plot_parser.add_argument("--begin", type=int, required=True)
    plot_parser.add_argument("--end", type=int, required=True)
    plot_parser.add_argument("--cmap", default=None)
    plot_parser.add_argument("--title", default=None)
    plot_parser.add_argument("--output", default=None)
    plot_parser.set_defaults(
        func=lambda args: plot(
            args.output_dir,
            begin=args.begin,
            end=args.end,
            title=args.title,
            cmap=args.cmap,
            output=args.output,
        )
    )

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
