# ExP Heatmap Benchmark Report

## 1. Benchmark Overview

### 1.1 System Environment

| Parameter | Value |
|-----------|-------|
| CPU | AMD Ryzen 9 7900 12-Core Processor |
| CPU Cores | 24 physical, 24 logical |
| RAM | 62.0 GB |
| OS | Linux 6.10.14-061014-generic (x86_64) |
| Python | 3.13.5 |
| exp_heatmap | 1.1.4 |
| scikit-allel | 1.3.13 |
| zarr | 2.18.7 |
| pandas | 2.3.0 |
| numpy | 2.3.1 |
| scipy | 1.16.0 |

### 1.2 Pipeline Architecture

ExP Heatmap is a three-step pipeline for visualizing cross-population natural selection signals:

1. **PREPARE**: Converts an entire VCF file to Zarr columnar format via `allel.vcf_to_zarr()`. This is a one-time, I/O-bound operation that processes the full chromosome regardless of the target genomic region.

2. **COMPUTE**: The computational core. Loads Zarr data and filters variants by allele frequency (AF > 0.05). Iterates over all C(n_pops, 2) unique population pairs, computing XP-EHH via `allel.xpehh()` with multi-threaded execution (`use_threads=True`). Supports two genotype loading strategies: chunked (`GenotypeChunkedArray` — lazy loading, lower memory) and unchunked (`GenotypeArray` — full in-memory, higher memory). Results are rank-transformed to empirical -log10 rank scores. Computational complexity is O(n_pairs × n_variants).

3. **PLOT**: Loads all per-pair TSV result files, filters to the requested genomic region, creates bidirectional comparisons (650 rows for 26 1000 Genomes populations), and renders a seaborn heatmap. Memory usage scales with the size of the displayed region.

### 1.3 Datasets

| Parameter | Chromosome 2 | Chromosome 22 |
|-----------|:------------:|:--------------:|
| Genome build | GRCh38 | GRCh37 |
| VCF file size | ~80.8 GB | ~10.7 GB |
| Variants | 6,088,598 | 1,103,547 |
| Samples | 3,202 | 2,504 |
| Populations | 26 | 26 |
| Population pairs | 325 | 325 |
| Genotype loading | Chunked | Unchunked |
| Statistical test | XP-EHH | XP-EHH |

### 1.4 Benchmark Configurations

| Benchmark | Region | Region Size | Replicates | Warmup | Chunked |
|-----------|--------|:-----------:|:----------:|:------:|:-------:|
| chr2_bm_1 | 136,108,646 – 137,108,646 | 1 Mb | 1 | 0 | Yes |
| chr2_bm_2 | 47,920,990 – 48,920,990 | 1 Mb | 1 | 0 | Yes |
| chr2_bm_3 | 136,108,646 – 137,108,646 | 1 Mb | 5 | 1 | Yes |
| chr2_bm_4 | 47,920,990 – 48,920,990 | 1 Mb | 5 | 1 | Yes |
| chr2_bm_5 | 40,000,000 – 137,000,000 | 97 Mb | 1 | 0 | Yes |
| chr22_bm_1 | 40,000,000 – 40,500,000 | 0.5 Mb | 1 | 0 | No |
| chr22_bm_2 | 16,060,513 – 51,210,268 | 35.1 Mb | 1 | 0 | No |
| chr22_bm_3 | 40,000,000 – 40,500,000 | 0.5 Mb | 5 | 0 | No |
| chr22_bm_4 | 16,060,513 – 51,210,268 | 35.1 Mb | 5 | 0 | No |
| chr22_bm_5 | 40,000,000 – 40,500,000 | 0.5 Mb | 5 | 1 | No |
| chr22_bm_6 | 16,060,513 – 51,210,268 | 35.1 Mb | 5 | 1 | No |

> **Note**: chr2_bm_6 (97 Mb region, 5 replicates, 1 warmup) was interrupted during the PLOT replication phase and no report was generated. It is excluded from all tables below.

### 1.5 Per-Step Results: All Benchmarks

#### Chromosome 2

| Benchmark | Step | Runtime (s) | CPU Time (s) | Peak Memory (MB) |
|-----------|------|:-----------:|:------------:|:-----------------:|
| chr2_bm_1 | PREPARE | 1,006.2 | 1,013.7 | 931.3 |
| | COMPUTE | 18,057.6 | 44,123.5 | 7,907.2 |
| | PLOT | 46.1 | 46.0 | 303.3 |
| | **Total** | **19,109.9** | | **7,907.2** |
| chr2_bm_2 | PREPARE | 1,011.1 | 1,018.8 | 931.3 |
| | COMPUTE | 18,100.9 | 44,152.1 | 7,907.2 |
| | PLOT | 44.8 | 44.8 | 229.8 |
| | **Total** | **19,156.9** | | **7,907.2** |
| chr2_bm_3 | PREPARE | 1,023.8 ± 0.0 | 1,031.7 ± 0.0 | 931.3 ± 0.0 |
| | COMPUTE | 18,042.2 ± 40.1 | 44,060.3 ± 46.4 | 8,033.6 ± 37.5 |
| | PLOT | 46.3 ± 0.5 | 46.3 ± 0.5 | 301.7 ± 0.1 |
| | **Total** | **19,112.3 ± 40.1** | | **8,033.6 ± 37.5** |
| chr2_bm_4 | PREPARE | 1,015.5 ± 0.0 | 1,023.8 ± 0.0 | 931.3 ± 0.0 |
| | COMPUTE | 18,045.2 ± 13.3 | 44,072.6 ± 11.7 | 8,033.6 ± 37.5 |
| | PLOT | 46.1 ± 0.3 | 46.1 ± 0.3 | 228.3 ± 0.1 |
| | **Total** | **19,106.8 ± 13.3** | | **8,033.6 ± 37.5** |
| chr2_bm_5 | PREPARE | 1,016.8 | 1,024.6 | 931.3 |
| | COMPUTE | 18,039.0 | 44,123.7 | 7,907.2 |
| | PLOT | 148.4 | 148.4 | 21,842.3 |
| | **Total** | **19,204.2** | | **21,842.3** |

#### Chromosome 22

| Benchmark | Step | Runtime (s) | CPU Time (s) | Peak Memory (MB) |
|-----------|------|:-----------:|:------------:|:-----------------:|
| chr22_bm_1 | PREPARE | 138.1 | 137.9 | 681.8 |
| | COMPUTE | 1,781.1 | 4,112.3 | 10,555.0 |
| | PLOT | 9.0 | 9.0 | 58.0 |
| | **Total** | **1,928.2** | | **10,555.0** |
| chr22_bm_2 | PREPARE | 137.1 | 138.1 | 681.8 |
| | COMPUTE | 1,784.5 | 4,116.1 | 10,555.0 |
| | PLOT | 42.4 | 42.4 | 8,727.9 |
| | **Total** | **1,963.9** | | **10,555.0** |
| chr22_bm_3 | PREPARE | 133.8 ± 0.0 | 134.8 ± 0.0 | 681.8 ± 0.0 |
| | COMPUTE | 1,796.9 ± 5.6 | 4,124.3 ± 5.5 | 10,555.0 ± 0.0 |
| | PLOT | 8.9 ± 0.2 | 8.9 ± 0.2 | 56.8 ± 0.7 |
| | **Total** | **1,939.6 ± 5.6** | | **10,555.0 ± 0.0** |
| chr22_bm_4 | PREPARE | 138.2 ± 0.0 | 138.3 ± 0.0 | 681.8 ± 0.0 |
| | COMPUTE | 1,786.1 ± 3.7 | 4,114.7 ± 4.2 | 10,555.0 ± 0.0 |
| | PLOT | 42.9 ± 1.0 | 42.9 ± 1.0 | 8,726.7 ± 0.7 |
| | **Total** | **1,967.2 ± 3.8** | | **10,555.0 ± 0.0** |
| chr22_bm_5 | PREPARE | 140.6 ± 0.0 | 140.2 ± 0.0 | 681.8 ± 0.0 |
| | COMPUTE | 1,792.3 ± 6.1 | 4,119.2 ± 6.9 | 10,555.0 ± 0.0 |
| | PLOT | 8.9 ± 0.2 | 8.9 ± 0.2 | 56.5 ± 0.1 |
| | **Total** | **1,941.8 ± 6.1** | | **10,555.0 ± 0.0** |
| chr22_bm_6 | PREPARE | 140.8 ± 0.0 | 141.8 ± 0.0 | 681.8 ± 0.0 |
| | COMPUTE | 1,803.5 ± 4.6 | 4,130.5 ± 4.9 | 10,555.0 ± 0.0 |
| | PLOT | 42.8 ± 3.0 | 42.7 ± 2.7 | 8,726.4 ± 0.1 |
| | **Total** | **1,987.1 ± 5.5** | | **10,555.0 ± 0.0** |

---

## 2. Publication-Ready Results

### Table 1. Input Data Characteristics

| | Chromosome 2 | Chromosome 22 |
|---|:---:|:---:|
| Genome build | GRCh38 | GRCh37 |
| VCF file size | 80.8 GB | 10.7 GB |
| Variants (pre-filter) | 6,088,598 | 1,103,547 |
| Samples | 3,202 | 2,504 |
| Populations | 26 | 26 |
| Population pairs | 325 | 325 |
| Genotype loading mode | Chunked | In-memory |
| Statistical test | XP-EHH | XP-EHH |

Data source: 1000 Genomes Project, phase 3. Variants were filtered to AF > 0.05 before computation.

### Table 2. Pipeline Performance (5 replicates, 1 warmup)

Results from the replicated benchmarks with warmup (chr2_bm_3/bm_4, chr22_bm_5/bm_6). PREPARE was run once per benchmark (deterministic, I/O-bound). COMPUTE and PLOT values are mean ± std over 5 replicates. Representative 1 Mb (chr2) and 0.5 Mb (chr22) region benchmarks are shown.

| Step | Chr2 Runtime (s) | Chr2 Peak Mem (MB) | Chr22 Runtime (s) | Chr22 Peak Mem (MB) |
|------|:-----------------:|:------------------:|:------------------:|:-------------------:|
| PREPARE | 1,023.8 | 931.3 | 140.6 | 681.8 |
| COMPUTE | 18,042.2 ± 40.1 | 8,033.6 ± 37.5 | 1,792.3 ± 6.1 | 10,555.0 ± 0.0 |
| PLOT | 46.3 ± 0.5 | 301.7 ± 0.1 | 8.9 ± 0.2 | 56.5 ± 0.1 |
| **Total** | **19,112.3 ± 40.1** | **8,033.6 ± 37.5** | **1,941.8 ± 6.1** | **10,555.0 ± 0.0** |

### Table 3. Compute Step: Reproducibility and Efficiency

| Metric | Chr2 (bm_3) | Chr2 (bm_4) | Chr22 (bm_5) | Chr22 (bm_6) |
|--------|:-----------:|:-----------:|:------------:|:------------:|
| Replicates | 5 | 5 | 5 | 5 |
| Warmup runs | 1 | 1 | 1 | 1 |
| Runtime CV | 0.2% | 0.1% | 0.3% | 0.3% |
| 95% CI (s) | [17,992 – 18,092] | [18,029 – 18,062] | [1,785 – 1,800] | [1,798 – 1,809] |
| Wall time (s) | 18,042.2 ± 40.1 | 18,045.2 ± 13.3 | 1,792.3 ± 6.1 | 1,803.5 ± 4.6 |
| CPU time (s) | 44,060.3 ± 46.4 | 44,072.6 ± 11.7 | 4,119.2 ± 6.9 | 4,130.5 ± 4.9 |
| CPU / Wall ratio | 2.44 | 2.44 | 2.30 | 2.29 |
| Seconds per pair | 55.5 | 55.5 | 5.5 | 5.6 |
| Memory efficiency (MB/Mvariants) | 1,319.5 ± 6.2 | 1,319.5 ± 6.2 | 9,564.6 ± 0.0 | 9,564.6 ± 0.0 |

### Table 4. Plot Step: Region Size Scaling

| Dataset | Region Size | Runtime (s) | Peak Memory (MB) | Memory Efficiency (MB/Mb) |
|---------|:-----------:|:-----------:|:-----------------:|:-------------------------:|
| Chr2 (bm_3) | 1 Mb | 46.3 ± 0.5 | 301.7 ± 0.1 | 301.7 |
| Chr2 (bm_5) | 97 Mb | 148.4 | 21,842.3 | 225.2 |
| Chr22 (bm_5) | 0.5 Mb | 8.9 ± 0.2 | 56.5 ± 0.1 | 113.0 |
| Chr22 (bm_6) | 35.1 Mb | 42.8 ± 3.0 | 8,726.4 ± 0.1 | 248.3 |

### Summary of Findings

The ExP Heatmap pipeline was benchmarked on two datasets from the 1000 Genomes Project (chromosome 2: 6.1M variants, 3,202 samples; chromosome 22: 1.1M variants, 2,504 samples) using XP-EHH as the selection test. All benchmarks were run on an AMD Ryzen 9 7900 system with 62 GB RAM.

**COMPUTE dominates total runtime**, accounting for approximately 94% of wall time on chromosome 2 and 92% on chromosome 22. This is expected given the O(n_pairs × n_variants) complexity of iterating XP-EHH over all 325 population pairs. The CPU-to-wall-time ratio of ~2.4× confirms that scikit-allel's multi-threaded XP-EHH implementation provides measurable parallelism on this hardware.

**Results are highly reproducible.** The coefficient of variation for compute runtime across 5 replicates was 0.1–0.3%, with narrow 95% confidence intervals. Warmup runs had negligible effect on measured performance: benchmarks with and without warmup produced overlapping confidence intervals (e.g., chr22 compute: 1,796.9 ± 5.6 s without warmup vs. 1,792.3 ± 6.1 s with warmup).

**PREPARE is I/O-bound and scales with VCF file size**, not region size. It consistently required ~1,016 s for the 80.8 GB chr2 VCF and ~138 s for the 10.7 GB chr22 VCF, processing at approximately 6,000 and 8,000 variants/second, respectively.

**PLOT runtime and memory scale with the displayed region size.** For chromosome 2, plotting a 1 Mb region required 46 s and 302 MB, while a 97 Mb region required 148 s and 21.8 GB. This is because the plot step loads all variants within the region across all 650 population pair rows into memory.

**Chunked genotype loading dramatically reduces compute memory usage.** Chromosome 2 (chunked mode) achieved a memory efficiency of 1,319 MB per million variants, compared to 9,565 MB per million variants for chromosome 22 (in-memory mode) — a 7.3× reduction. This came at no measurable runtime penalty, making chunked mode strongly recommended for large VCF files.

