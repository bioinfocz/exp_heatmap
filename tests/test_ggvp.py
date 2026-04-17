from exp_heatmap.ggvp import (
    build_integrated_panel,
    filter_manifest_for_chromosome,
    load_1kg_panel_mapping,
    load_ggvp_sample_mapping,
    load_population_superpop_mapping,
    load_sequence_index_population_mapping,
    parse_manifest,
)


def test_load_ggvp_sample_mapping_parses_population_from_alignment_paths():
    index_text = "\n".join(
        [
            "#CRAM_FILE\tCRAM_MD5",
            "ftp://example/ggvp/data/FULA/SC_GMFUL1/alignment/file.cram\tabc",
            "ftp://example/ggvp/data/JOLA/SC_GMJOL1/alignment/file.cram\tdef",
        ]
    )

    mapping = load_ggvp_sample_mapping(index_text, super_pop="AFR")

    assert mapping["SC_GMFUL1"] == {"pop": "FULA", "super_pop": "AFR"}
    assert mapping["SC_GMJOL1"] == {"pop": "JOLA", "super_pop": "AFR"}


def test_build_integrated_panel_respects_vcf_sample_order():
    ggvp_mapping = {"SC_GMFUL1": {"pop": "FULA", "super_pop": "AFR"}}
    onekg_mapping = {"HG00096": {"pop": "GBR", "super_pop": "EUR"}}

    panel = build_integrated_panel(["HG00096", "SC_GMFUL1"], ggvp_mapping, onekg_mapping)

    assert panel["sample"].tolist() == ["HG00096", "SC_GMFUL1"]
    assert panel["pop"].tolist() == ["GBR", "FULA"]
    assert panel["super_pop"].tolist() == ["EUR", "AFR"]


def test_load_1kg_panel_mapping_reads_official_panel_shape():
    panel_text = "sample\tpop\tsuper_pop\tgender\nHG00096\tGBR\tEUR\tmale\n"
    mapping = load_1kg_panel_mapping(panel_text)
    assert mapping["HG00096"] == {"pop": "GBR", "super_pop": "EUR"}


def test_load_sequence_index_population_mapping_uses_population_lookup():
    panel_text = "sample\tpop\tsuper_pop\tgender\nHG00096\tGBR\tEUR\tmale\nHG02461\tGWD\tAFR\tfemale\n"
    index_text = "\n".join(
        [
            "FASTQ_FILE\tMD5\tRUN_ID\tSTUDY_ID\tSTUDY_NAME\tCENTER_NAME\tSUBMISSION_ID\tSUBMISSION_DATE\tSAMPLE_ID\tSAMPLE_NAME\tPOPULATION",
            "data/example.fastq.gz\t.\tRUN1\tSTUDY\tStudy\tSC\tSUB\t2020-01-01\tS1\tHG03033\tGWD",
        ]
    )

    mapping = load_sequence_index_population_mapping(
        index_text,
        load_population_superpop_mapping(panel_text),
    )

    assert mapping["HG03033"] == {"pop": "GWD", "super_pop": "AFR"}


def test_manifest_filter_selects_vcf_and_tbi_for_requested_chromosome():
    manifest_text = "\n".join(
        [
            "./ALL_GGVP.chr2.shapeit2_integrated_snvindels_v1b_20200120.GRCh38.phased.vcf.gz\t10\tabc",
            "./ALL_GGVP.chr2.shapeit2_integrated_snvindels_v1b_20200120.GRCh38.phased.vcf.gz.tbi\t20\tdef",
            "./ALL_GGVP.chr3.shapeit2_integrated_snvindels_v1b_20200120.GRCh38.phased.vcf.gz\t30\tghi",
        ]
    )
    entries = parse_manifest(manifest_text)
    selected = filter_manifest_for_chromosome(entries, "2")
    filenames = [entry.filename for entry in selected]
    assert filenames == [
        "ALL_GGVP.chr2.shapeit2_integrated_snvindels_v1b_20200120.GRCh38.phased.vcf.gz",
        "ALL_GGVP.chr2.shapeit2_integrated_snvindels_v1b_20200120.GRCh38.phased.vcf.gz.tbi",
    ]
