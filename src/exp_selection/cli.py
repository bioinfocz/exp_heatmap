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

    prepare_parse = subparser.add_parser("prepare")
    prepare_parse.add_argument(
        "recode_file", metavar="RECODE_FILE", help="Where is recoded VCF"
    )
    prepare_parse.add_argument(
        "zarr_dir", metavar="ZARR_DIR", help="Where ZARR dir will be saved"
    )
    prepare_parse.set_defaults(
        func=lambda args: prepare(args.recode_file, args.zarr_dir)
    )

    compute_parse = subparser.add_parser("compute")
    compute_parse.add_argument(
        "zarr_dir", metavar="ZARR_DIR", help="Where ZARR dir is located"
    )
    compute_parse.add_argument(
        "panel_file", metavar="PANEL_FILE", help="Where panel file will be saved"
    )
    compute_parse.add_argument(
        "xpehh_dir", metavar="XPEHH_DIR", help="Where xpehh dir will be saved"
    )
    compute_parse.set_defaults(
        func=lambda args: compute(args.zarr_dir, args.panel_file, args.xpehh_dir)
    )

    compute_parse = subparser.add_parser("plot")
    compute_parse.add_argument(
        "xpehh_dir", metavar="XPEHH_DIR", help="Where xpehh dir is located"
    )
    parser.add_argument("--begin", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--cmap", default=None)
    parser.add_argument("--title", default=None)
    compute_parse.set_defaults(
        func=lambda args: plot(
            args.xpehh_dir,
            begin=args.begin,
            end=args.end,
            title=args.title,
            cmap=args.cmap,
        )
    )

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
