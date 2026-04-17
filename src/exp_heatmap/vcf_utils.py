import gzip


def is_biallelic_snp(ref, alt):
    if len(ref) != 1 or not alt or alt == ".":
        return False
    alleles = alt.split(",")
    return len(alleles) == 1 and len(alleles[0]) == 1


def _open_text(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def _normalize_chromosome_name(chromosome):
    value = str(chromosome).strip()
    if value.lower().startswith("chr"):
        return value[3:]
    return value


def _record_is_within_region(record_chromosome, record_position, chromosome=None, start=None, end=None):
    if chromosome is not None and _normalize_chromosome_name(record_chromosome) != _normalize_chromosome_name(chromosome):
        return False
    if start is not None and record_position < start:
        return False
    if end is not None and record_position > end:
        return False
    return True


def filter_biallelic_snp_vcf(input_path, output_path, chromosome=None, start=None, end=None):
    header_lines = 0
    total_records = 0
    kept_records = 0

    with _open_text(input_path, "rt") as in_handle, _open_text(output_path, "wt") as out_handle:
        for line in in_handle:
            if line.startswith("#"):
                out_handle.write(line)
                header_lines += 1
                continue

            total_records += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue

            record_chromosome = fields[0]
            try:
                record_position = int(fields[1])
            except ValueError:
                continue

            if not _record_is_within_region(
                record_chromosome,
                record_position,
                chromosome=chromosome,
                start=start,
                end=end,
            ):
                continue

            ref = fields[3]
            alt = fields[4]
            if is_biallelic_snp(ref, alt):
                out_handle.write(line)
                kept_records += 1

    return {
        "header_lines": header_lines,
        "total_records": total_records,
        "kept_records": kept_records,
        "dropped_records": total_records - kept_records,
    }
