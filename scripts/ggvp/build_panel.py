#!/usr/bin/env python3

import argparse
from pathlib import Path

from exp_heatmap.ggvp import (
    GGVP_ALIGNMENT_INDEX_URL,
    ONEKG_PANEL_URL,
    ONEKG_PHASE3_SEQUENCE_INDEX_URL,
    build_integrated_panel,
    fetch_vcf_sample_names,
    load_1kg_panel_mapping,
    load_ggvp_sample_mapping,
    load_population_superpop_mapping,
    load_sequence_index_population_mapping,
    official_ggvp_vcf_url,
    read_text_url,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a GGVP-compatible panel file in VCF sample order."
    )
    parser.add_argument(
        "--vcf-url",
        help="GGVP VCF URL to read the sample header from. Defaults to chr2.",
    )
    parser.add_argument(
        "--chromosome",
        default="2",
        help="Chromosome used when constructing the default GGVP VCF URL.",
    )
    parser.add_argument(
        "--alignment-index-url",
        default=GGVP_ALIGNMENT_INDEX_URL,
        help="Official GGVP alignment index URL used to map SC_* samples to populations.",
    )
    parser.add_argument(
        "--onekg-panel-url",
        default=ONEKG_PANEL_URL,
        help="Official 1000 Genomes panel URL used to map HG* samples already present in the integrated GGVP VCF.",
    )
    parser.add_argument(
        "--onekg-sequence-index-url",
        default=ONEKG_PHASE3_SEQUENCE_INDEX_URL,
        help="Official 1000 Genomes phase-3 sequence index used to recover HG* samples missing from the legacy ALL.panel file.",
    )
    parser.add_argument(
        "--ggvp-super-pop",
        default="AFR",
        help="super_pop value assigned to GGVP-specific populations.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV path for the integrated panel file.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    vcf_url = args.vcf_url or official_ggvp_vcf_url(args.chromosome)

    sample_names = fetch_vcf_sample_names(vcf_url)
    ggvp_mapping = load_ggvp_sample_mapping(
        read_text_url(args.alignment_index_url),
        super_pop=args.ggvp_super_pop,
    )
    onekg_panel_text = read_text_url(args.onekg_panel_url)
    onekg_mapping = load_1kg_panel_mapping(onekg_panel_text)
    population_to_super_pop = load_population_superpop_mapping(onekg_panel_text)
    onekg_sequence_mapping = load_sequence_index_population_mapping(
        read_text_url(args.onekg_sequence_index_url),
        population_to_super_pop,
    )
    onekg_sequence_mapping.update(onekg_mapping)

    panel_df = build_integrated_panel(sample_names, ggvp_mapping, onekg_sequence_mapping)
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    panel_df.to_csv(output_path, sep="\t", index=False)

    print(f"Wrote panel with {len(panel_df)} rows to {output_path}")
    print(panel_df["pop"].value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
