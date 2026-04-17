import gzip

from exp_heatmap.vcf_utils import filter_biallelic_snp_vcf, is_biallelic_snp


def test_is_biallelic_snp_accepts_single_base_substitution():
    assert is_biallelic_snp("A", "G")
    assert not is_biallelic_snp("AT", "G")
    assert not is_biallelic_snp("A", "GA")
    assert not is_biallelic_snp("A", "G,T")


def test_filter_biallelic_snp_vcf_filters_indels_and_multiallelic_records(tmp_path):
    input_path = tmp_path / "input.vcf.gz"
    output_path = tmp_path / "output.vcf"
    with gzip.open(input_path, "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        handle.write("21\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\n")
        handle.write("21\t110\trs2\tAT\tG\t.\tPASS\t.\tGT\t0|1\n")
        handle.write("21\t120\trs3\tA\tGA\t.\tPASS\t.\tGT\t0|1\n")
        handle.write("21\t130\trs4\tC\tT,G\t.\tPASS\t.\tGT\t1|2\n")
        handle.write("21\t140\trs5\tG\tA\t.\tPASS\t.\tGT\t0|0\n")

    stats = filter_biallelic_snp_vcf(str(input_path), str(output_path))

    assert stats["total_records"] == 5
    assert stats["kept_records"] == 2
    assert output_path.read_text().splitlines() == [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        "21\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1",
        "21\t140\trs5\tG\tA\t.\tPASS\t.\tGT\t0|0",
    ]


def test_filter_biallelic_snp_vcf_can_limit_to_region(tmp_path):
    input_path = tmp_path / "input.vcf.gz"
    output_path = tmp_path / "output.vcf"
    with gzip.open(input_path, "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        handle.write("chr2\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\n")
        handle.write("chr2\t150\trs2\tC\tT\t.\tPASS\t.\tGT\t0|1\n")
        handle.write("chr2\t250\trs3\tG\tA\t.\tPASS\t.\tGT\t0|0\n")
        handle.write("chr3\t150\trs4\tT\tC\t.\tPASS\t.\tGT\t0|1\n")

    stats = filter_biallelic_snp_vcf(
        str(input_path),
        str(output_path),
        chromosome="2",
        start=120,
        end=200,
    )

    assert stats["total_records"] == 4
    assert stats["kept_records"] == 1
    assert output_path.read_text().splitlines() == [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        "chr2\t150\trs2\tC\tT\t.\tPASS\t.\tGT\t0|1",
    ]
