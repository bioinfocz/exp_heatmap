#!/bin/bash
set -euo pipefail

# Author-endorsed full 1000 Genomes raw-data validation workflow for SLC24A5.
# This mirrors the logic of the README example while using exp_heatmap filter-vcf
# instead of vcftools so the run stays self-contained in the current local environment.

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
WORK_DIR="${ROOT_DIR}/local_data/1kg_chr15_slc24a5"
SHARED_DIR="${ROOT_DIR}/local_data/1kg_shared"

mkdir -p "${WORK_DIR}" "${SHARED_DIR}"

export PYTHONPATH="${ROOT_DIR}/src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib}"

CHR15_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
PANEL_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

curl -L "${CHR15_URL}" -o "${WORK_DIR}/chr15.vcf.gz"
curl -L "${PANEL_URL}" -o "${SHARED_DIR}/genotypes.panel"

python3 -m exp_heatmap.cli filter-vcf "${WORK_DIR}/chr15.vcf.gz" \
  -o "${WORK_DIR}/chr15_snps.vcf" \
  --no-log

python3 -m exp_heatmap.cli prepare "${WORK_DIR}/chr15_snps.vcf" \
  -o "${WORK_DIR}/chr15.zarr" \
  --no-log

python3 -m exp_heatmap.cli compute "${WORK_DIR}/chr15.zarr" "${SHARED_DIR}/genotypes.panel" \
  -o "${WORK_DIR}/chr15_output" \
  -t xpehh \
  -c \
  --no-log

python3 -m exp_heatmap.cli plot "${WORK_DIR}/chr15_output" \
  --start 47924019 \
  --end 48924019 \
  --title "SLC24A5" \
  --cmap gist_heat \
  --out "${WORK_DIR}/SLC24A5_heatmap" \
  --no-log

python3 -m exp_heatmap.cli plot "${WORK_DIR}/chr15_output" \
  --start 47924019 \
  --end 48924019 \
  --interactive \
  --title "SLC24A5" \
  --out "${WORK_DIR}/SLC24A5_heatmap_interactive" \
  --no-log
