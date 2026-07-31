# Full-pipeline benchmark summary

- Input VCF: `/home/adam/Documents/PhD/Cursor_projects/exp_heatmap/local_data/ggvp/ggvp_chr21_snps.vcf`
- Panel file: `/home/adam/Documents/PhD/Cursor_projects/exp_heatmap/local_data/ggvp/ggvp_chr21.panel`
- Test: `xpehh`
- Rank-score mode for plotting: `directional`
- Output directory: `/tmp/claude-1002/-home-adam-Documents-PhD-Cursor-projects-exp-heatmap/81ba78c2-d814-4cf3-8ccc-ac8eeb745d51/scratchpad/pipeline_ggvp`

| case             | compute_artifact_size_mb | plot_artifact_size_mb | prepare_artifact_size_mb | compute_peak_rss_gb | plot_peak_rss_gb | prepare_peak_rss_gb | compute_seconds | plot_seconds | prepare_seconds |
| ---------------- | ------------------------ | --------------------- | ------------------------ | ------------------- | ---------------- | ------------------- | --------------- | ------------ | --------------- |
| full_chr21       | 71.5035                  | 1.0084                | 79.2802                  | 0.5708              | 2.7633           | 0.237               | 54.835455       | 6.313262     | 13.971018       |
| small_1mb_region | 1.4911                   | 0.2579                | 2.0261                   | 0.1249              | 0.4369           | 0.1614              | 1.720863        | 1.517092     | 0.880199        |

## Machine context

- Platform: `Linux-6.10.14-061014-generic-x86_64-with-glibc2.39`
- Processor: `AMD Ryzen 9 7900 12-Core Processor`
- Logical CPU count: `24`
- Total RAM (GB): `61.96`

## Notes

- The smaller case uses a 1 Mb subset derived from the same filtered GGVP chromosome 21 VCF.
- The larger case uses the full filtered GGVP chromosome 21 VCF.
- Peak RAM is the observed maximum resident set size across the monitored process tree.
- Raw command logs are stored in the `logs/` subdirectory.
