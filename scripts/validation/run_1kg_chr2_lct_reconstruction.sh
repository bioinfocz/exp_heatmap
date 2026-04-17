#!/bin/bash
set -euo pipefail

# Local public-data reconstruction of the 1000 Genomes chr2 LCT workflow.
# This is not the archived preprepared LCT tutorial bundle; it is a transparent
# reconstruction from public inputs using the same pipeline logic.
#
# The run is region-scoped on purpose. We keep the public raw-data -> filter ->
# prepare -> compute -> plot pipeline, but filter chromosome 2 to the LCT
# plotting window with 1 Mb of flanking sequence on each side before prepare.
# This keeps the reconstruction practical while still leaving enough context
# around the displayed locus for local figure generation and consistency checks.

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
WORK_DIR="${ROOT_DIR}/local_data/1kg_chr2_lct_region"
SHARED_DIR="${ROOT_DIR}/local_data/1kg_shared"

mkdir -p "${WORK_DIR}" "${SHARED_DIR}"

export PYTHONPATH="${ROOT_DIR}/src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib}"

CHR2_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
PANEL_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

PLOT_START=136108646
PLOT_END=137108646
FLANK_BP=1000000
FILTER_CHROM=2
FILTER_START=$((PLOT_START - FLANK_BP))
FILTER_END=$((PLOT_END + FLANK_BP))

curl -L "${CHR2_URL}" -o "${WORK_DIR}/chr2.vcf.gz"
curl -L "${PANEL_URL}" -o "${SHARED_DIR}/genotypes.panel"

python3 -m exp_heatmap.cli filter-vcf "${WORK_DIR}/chr2.vcf.gz" \
  -o "${WORK_DIR}/chr2_lct_region_snps.vcf" \
  --chrom "${FILTER_CHROM}" \
  --start "${FILTER_START}" \
  --end "${FILTER_END}" \
  --no-log

python3 -m exp_heatmap.cli prepare "${WORK_DIR}/chr2_lct_region_snps.vcf" \
  -o "${WORK_DIR}/chr2_lct_region.zarr" \
  --no-log

python3 -m exp_heatmap.cli compute "${WORK_DIR}/chr2_lct_region.zarr" "${SHARED_DIR}/genotypes.panel" \
  -o "${WORK_DIR}/chr2_lct_region_output" \
  -t xpehh \
  -c \
  --no-log

python3 -m exp_heatmap.cli plot "${WORK_DIR}/chr2_lct_region_output" \
  --start "${PLOT_START}" \
  --end "${PLOT_END}" \
  --title "LCT gene" \
  --out "${WORK_DIR}/LCT_xpehh" \
  --no-log

python3 -m exp_heatmap.cli plot "${WORK_DIR}/chr2_lct_region_output" \
  --start "${PLOT_START}" \
  --end "${PLOT_END}" \
  --interactive \
  --title "LCT gene" \
  --out "${WORK_DIR}/LCT_xpehh_interactive" \
  --no-log
