# v1.3.1 benchmark results

Produced with the replicate-capable harness on `main` (commit `73ad651`, branch
`an-benchmark-replicates`). Every directory contains the results TSV, a Markdown summary
with full provenance, `machine_specs.json`, and one log per individual run.

Machine for all runs: AMD Ryzen 9 7900, 12 physical / 24 logical cores, 61.96 GB RAM,
Linux 6.10.14, Python 3.12.2.

| Directory | Benchmark | Replicates |
|---|---|---|
| `population_scaling/` | Display and reporting cost from 5 to 50 populations, 600 dpi, 1 Mb window, seeded from the 1000 Genomes chr2 LCT window | 5 + 1 warmup |
| `sample_scaling_xpehh/` | compute cost versus samples per population, XP-EHH, 1000 Genomes chr22 | 2 + 1 warmup |
| `sample_scaling_hudson_fst/` | The same for Hudson's FST | 2 + 1 warmup |
| `pipeline_ggvp_chr21/` | Full prepare -> compute -> plot on GGVP chr21, downloaded fresh from the official release. Replaces manuscript Table 3 | 5 + 1 warmup (prepare once) |

Statistics are mean, sample standard deviation (ddof=1), coefficient of variation and a
Student-t 95% confidence interval, matching `src/exp_heatmap/benchmark_stats.py` on this
branch. Individual per-run timings are preserved in the `seconds_runs` column.

Population scaling grows the panel by symlinking template per-pair outputs under invented
labels, so it measures the display and reporting layer only and never invokes `compute`.
Sample scaling varies real genotypes and does invoke `compute`.

The earlier `chr2/` and `chr22/` results in the parent directory were measured on v1.1.4
and are retained for the variability result; their absolute timings do not describe v1.3.1.

## GGVP chr21 provenance

VCF and index from the 20200217_biallelic_SNV release, MD5 verified against the release
manifest. Panel built by `scripts/ggvp/build_panel.py`: 505 samples across FULA (92),
GWD (118), JOLA (98), MANDINKA (99), WOLOFF (98). `filter-vcf` retained 509,924 biallelic
SNP records from 553,906 variant records, matching the counts in the manuscript. The
allele-frequency filter then retains 163,102 of those 509,924 SNPs.

Run at 600 dpi with a 4,000-column budget, so the plot-stage figures describe the cost of
producing a publication figure rather than a screen-resolution render.
