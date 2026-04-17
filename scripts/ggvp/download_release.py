#!/usr/bin/env python3

import argparse
from pathlib import Path

from exp_heatmap.ggvp import (
    GGVP_MANIFEST_URL,
    GGVP_RELEASE_BASE_URL,
    download_file,
    filter_manifest_for_chromosome,
    parse_manifest,
    read_text_url,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download GGVP release files for selected chromosomes and verify checksums."
    )
    parser.add_argument(
        "--chromosomes",
        nargs="+",
        required=True,
        help="Chromosomes to download, e.g. 2 15 X",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Destination directory for downloaded files.",
    )
    parser.add_argument(
        "--manifest-url",
        default=GGVP_MANIFEST_URL,
        help="Manifest URL used for file discovery and checksum validation.",
    )
    parser.add_argument(
        "--base-url",
        default=GGVP_RELEASE_BASE_URL,
        help="Base GGVP release URL.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest_entries = parse_manifest(read_text_url(args.manifest_url))
    selected_entries = []
    for chromosome in args.chromosomes:
        selected_entries.extend(filter_manifest_for_chromosome(manifest_entries, chromosome))

    if not selected_entries:
        raise SystemExit("No manifest entries matched the requested chromosomes.")

    for entry in selected_entries:
        source_url = f"{args.base_url}/{entry.filename}"
        destination = output_dir / entry.filename
        print(f"Downloading {entry.filename} -> {destination}")
        download_file(source_url, destination, expected_md5=entry.md5)

    print(f"Downloaded {len(selected_entries)} files into {output_dir}")


if __name__ == "__main__":
    main()
