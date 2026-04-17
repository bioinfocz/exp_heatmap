import csv
import gzip
import hashlib
import os
import urllib.request
from dataclasses import dataclass

import pandas as pd


GGVP_RELEASE_BASE_URL = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "gambian_genome_variation_project/release/20200217_biallelic_SNV"
)
GGVP_MANIFEST_URL = f"{GGVP_RELEASE_BASE_URL}/20200217_GGVP_vc_manifest.txt"
GGVP_ALIGNMENT_INDEX_URL = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "gambian_genome_variation_project/gambian_genome_variation_project.GRCh38DH.alignment.index"
)
ONEKG_PANEL_URL = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)
ONEKG_PHASE3_SEQUENCE_INDEX_URL = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.sequence.index"
)


@dataclass(frozen=True)
class ManifestEntry:
    path: str
    size: int
    md5: str

    @property
    def filename(self):
        return os.path.basename(self.path)


def read_text_url(url):
    with urllib.request.urlopen(url) as response:
        return response.read().decode("utf-8")


def parse_manifest(manifest_text):
    entries = []
    for line in manifest_text.splitlines():
        line = line.strip()
        if not line:
            continue
        path, size, md5 = line.split("\t")
        entries.append(ManifestEntry(path=path, size=int(size), md5=md5))
    return entries


def filter_manifest_for_chromosome(entries, chromosome):
    chromosome = str(chromosome)
    needle = f".chr{chromosome}."
    selected = [entry for entry in entries if needle in entry.filename]
    selected.extend(
        entry for entry in entries if needle in entry.filename.replace(".vcf.gz.tbi", ".vcf.gz")
        and entry.filename.endswith(".vcf.gz.tbi")
    )
    deduped = {}
    for entry in selected:
        deduped[entry.filename] = entry
    return [deduped[name] for name in sorted(deduped)]


def iter_gzip_text_lines_from_url(url):
    with urllib.request.urlopen(url) as response:
        with gzip.GzipFile(fileobj=response) as gz:
            for raw_line in gz:
                yield raw_line.decode("utf-8")


def fetch_vcf_sample_names(vcf_url):
    for line in iter_gzip_text_lines_from_url(vcf_url):
        if line.startswith("#CHROM"):
            fields = line.rstrip("\n").split("\t")
            return fields[9:]
    raise ValueError(f"No #CHROM header line found in {vcf_url}")


def load_ggvp_sample_mapping(index_text, super_pop="AFR"):
    sample_mapping = {}
    for line in index_text.splitlines():
        if not line or line.startswith("#"):
            continue
        path = line.split("\t", 1)[0]
        parts = path.split("/")
        population = parts[-4]
        sample = parts[-3]
        previous = sample_mapping.get(sample)
        if previous and previous["pop"] != population:
            raise ValueError(f"Sample {sample} appears with conflicting GGVP populations")
        sample_mapping[sample] = {"pop": population, "super_pop": super_pop}
    return sample_mapping


def load_1kg_panel_mapping(panel_text):
    reader = csv.DictReader(panel_text.splitlines(), delimiter="\t")
    mapping = {}
    for row in reader:
        mapping[row["sample"]] = {
            "pop": row["pop"],
            "super_pop": row["super_pop"],
        }
    return mapping


def load_population_superpop_mapping(panel_text):
    reader = csv.DictReader(panel_text.splitlines(), delimiter="\t")
    mapping = {}
    for row in reader:
        population = row["pop"]
        super_pop = row["super_pop"]
        previous = mapping.get(population)
        if previous and previous != super_pop:
            raise ValueError(
                f"Population {population} appears with conflicting super-population labels"
            )
        mapping[population] = super_pop
    return mapping


def load_sequence_index_population_mapping(index_text, population_to_super_pop):
    reader = csv.DictReader(index_text.splitlines(), delimiter="\t")
    mapping = {}
    for row in reader:
        sample = row["SAMPLE_NAME"]
        population = row["POPULATION"]
        if not sample or not population:
            continue
        super_pop = population_to_super_pop.get(population)
        if super_pop is None:
            continue

        previous = mapping.get(sample)
        if previous and previous["pop"] != population:
            raise ValueError(
                f"Sample {sample} appears with conflicting 1000 Genomes populations"
            )
        mapping[sample] = {"pop": population, "super_pop": super_pop}
    return mapping


def build_integrated_panel(sample_names, ggvp_mapping, onekg_mapping):
    rows = []
    missing = []
    for sample in sample_names:
        if sample in ggvp_mapping:
            entry = ggvp_mapping[sample]
        elif sample in onekg_mapping:
            entry = onekg_mapping[sample]
        else:
            missing.append(sample)
            continue
        rows.append(
            {
                "sample": sample,
                "pop": entry["pop"],
                "super_pop": entry["super_pop"],
            }
        )

    if missing:
        missing_preview = ", ".join(missing[:10])
        raise ValueError(
            f"Unable to map {len(missing)} samples from the GGVP VCF header. "
            f"Examples: {missing_preview}"
        )

    return pd.DataFrame(rows, columns=["sample", "pop", "super_pop"])


def official_ggvp_vcf_url(chromosome):
    chromosome = str(chromosome)
    return (
        f"{GGVP_RELEASE_BASE_URL}/"
        f"ALL_GGVP.chr{chromosome}.shapeit2_integrated_snvindels_v1b_20200120.GRCh38.phased.vcf.gz"
    )


def official_ggvp_tbi_url(chromosome):
    return official_ggvp_vcf_url(chromosome) + ".tbi"


def download_file(url, destination_path, expected_md5=None):
    with urllib.request.urlopen(url) as response, open(destination_path, "wb") as out_handle:
        digest = hashlib.md5()
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out_handle.write(chunk)
            digest.update(chunk)

    if expected_md5 and digest.hexdigest() != expected_md5:
        raise ValueError(
            f"Checksum mismatch for {destination_path}: expected {expected_md5}, got {digest.hexdigest()}"
        )

    return destination_path
