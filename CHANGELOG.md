# Changelog

All notable changes to ExP Heatmap are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.2] - 2026-07-31

### Added

#### Declared population panels on the command line (#9)

- New `--populations` option on `exp_heatmap plot` and `exp_heatmap full`, taking
  comma-separated population codes, for example `--populations "GWD,MSL,ESN"`. It is the
  command-line equivalent of the existing `populations=` API parameter.
- Declaring the panel makes an incomplete compute directory fail loudly. Previously,
  compute outputs covering only a subset of populations formed a self-consistent smaller
  set and were rendered without warning, so an interrupted `compute` run silently
  produced a smaller heatmap. `verify_declared_panel_is_complete` now raises with the
  number of missing pairs out of the expected total, a preview of up to five missing
  pairs, and the names of any populations that have no output at all.
- Input is validated before any work starts: fewer than two codes, or duplicate codes,
  raise a `UsageError` naming the problem.
- A declared panel that exactly matches the canonical 26-population 1000 Genomes tuple is
  mapped back to `"1000Genomes"` mode, so it keeps the standard 1000G row order and
  super-population annotations. Any other declared order is used verbatim.
- `ValueError` raised by the plotting layer is surfaced as a `ClickException` rather than
  a traceback.

#### Warning when data exceeds the color-scale maximum

- `warn_if_values_are_clipped` logs a warning giving the number of clipped cells, the
  total cell count, the color-scale maximum, and the observed data maximum. It is wired
  into both the static (matplotlib) and interactive (Plotly) renderers.
- Only the upper bound is reported. Values below the lower bound are clipped by design,
  since the 1000 Genomes lower bound is the 5% tail cutoff.
- Motivation: empirical rank scores reach log10(n) for n ranked loci, so runs with more
  than roughly 68,000 loci per pair exceed the fixed 4.833 ceiling. Those cells all render
  at the top color and cannot be told apart, which previously gave no visible signal.
- The bounds `1.301` and `4.833`, previously duplicated as literals across `plot.py` and
  `interactive.py`, are now the named constants `CBAR_VMIN_1000G` and `CBAR_VMAX_1000G`.

#### Allele-frequency filter reporting

- `filter_by_AF` now reports, at info level, how many variants were kept, how many were
  removed, out of how many, and at what threshold. This filter determines which loci reach
  the output and was not otherwise visible to the user.
- The "using precomputed alternate allele frequencies from `variants/AF`" message was
  promoted from debug to info. `INFO/AF` is cohort-wide whereas the genotype-derived
  fallback covers only the samples present in the file, so the two can select different
  variants and users need to know which source was used.

#### Benchmark harness: replicates, warmup, and sample scaling

- New `scripts/benchmarks/_bench_stats.py`: shared replicate handling providing the sample
  standard deviation (ddof=1), a Student-t confidence interval, and the coefficient of
  variation, matching the formulas of the `benchmark_stats.py` harness on the `benchmarks`
  branch.
- `--repeats` and `--warmup` added to the pipeline, population-scaling, and display
  benchmarks. Every run keeps its own log, and the individual per-run timings are recorded
  in a `seconds_runs` column alongside the summary statistics.
- New `scripts/benchmarks/run_sample_scaling.py`: varies samples per population on real
  genotypes and measures `compute`, the stage that dominates whole-chromosome runtime.
  Because `compute` requires the panel to match the store exactly, each sample size gets
  its own Zarr store; `variants/AF` is copied rather than recomputed so the
  allele-frequency filter keeps the same loci at every size. Store construction is untimed
  setup.

#### Tests

- New `tests/test_compute_fidelity.py`: asserts that the per-pair `xpehh` and `hudson_fst`
  values equal a direct scikit-allel call on the same haplotypes, and that two runs are
  byte-identical. This tests the wrapper rather than re-testing scikit-allel.
- New coverage in `tests/test_cli.py` (population forwarding and validation),
  `tests/test_plot.py` (panel resolution, completeness checks, clipping warnings, numeric
  column index, region-crop fallbacks), and `tests/test_xp_utils.py` (filter counts,
  frequency source, genotype fallback warning). 24 test functions added in total.

### Fixed

- **Coordinate slicing after plot-input construction.** Adding the `first_pop` and
  `second_pop` sort keys promoted the genomic-position column index to object dtype, and
  dropping them did not restore it, so callers could not slice a region by coordinate.
  `create_plot_input` now restores a numeric column index.
- **Region cropping on externally built frames.** `plot_exp_heatmap` caught only
  `KeyError` when a requested bound was absent, but a non-numeric column index raises
  `TypeError`. Both are now caught and fall back to the nearest available positions.
- **Machine specs on Linux.** `machine_specs` read a macOS-only sysctl key and fell back
  to `platform.processor()`, which records nothing useful on most Linux distributions. It
  now reads `/proc/cpuinfo` on Linux, and reports physical and logical core counts
  separately.
- **Zenodo badge rendering in the README.**

### Changed

- `plot_exp_heatmap` no longer duplicates the `seaborn.heatmap` call across a 1000 Genomes
  branch and a non-1000 Genomes branch. The color bounds are resolved first into
  `effective_vmin` and `effective_vmax`, and there is a single call site.
- `plot_interactive` accepts a `populations` argument and forwards it to
  `create_plot_input`. The value recovered from the data is now held separately as
  `resolved_populations`.

### Dependencies

- `scipy>=1.10` added to the `benchmarks` optional-dependency group, for the Student-t
  confidence interval.

### Documentation

- Zenodo DOI badge added to the README.
- New "Color Scale Range" section, including "Values above the upper bound are clipped".
- Panel resolution documented as inferred (default) versus declared, with a worked GGVP
  example of a `compute` run interrupted before it reached the `MANDINKA` pairs.
- Two new troubleshooting entries: "Strongest Signals All Render as One Flat Block" and
  "Heatmap Has Fewer Populations Than Expected".
- Benchmark commands corrected. The documented invocations used flags that do not exist;
  the arguments are positional. `--repeats`, `--warmup`, and the `--dpi` caveat are now
  documented, along with a new "Compute scaling with sample count" section.
- Guidance added on direction-aware reciprocal ranking versus one-sided ranking from the
  same compute output, and on passing the region as rendered rather than as requested,
  since `create_plot_input` snaps the window to the nearest variant positions present in
  the data.

## [1.3.1] - 2026-05-19

### Fixed

- `extract_top_regions`: reported regions are now separated by a minimum edge-to-edge gap
  in base pairs rather than by center-to-center distance, so adjacent windows are no
  longer reported as distinct hits. `start` and `end` are now the ±`window_size`/2
  boundaries around `center`, which keeps them well defined when only one column position
  falls inside the window, as happens with pre-summarized 25 kb windows. A `n_positions`
  column was added; the legacy `n_variants` name is retained as a duplicate for backwards
  compatibility.
- Interactive downsampling: bins are now formed with `np.array_split`, so all columns are
  covered when the column count is not an exact multiple of `max_columns`. The previous
  integer-division binning dropped the trailing remainder, and the representative position
  is now taken from the bin's own indices.
- `cbar_kws` is normalized across input types: colorbar ticks are only added when ticks
  were actually supplied.
- Bare `except:` in the region-cropping path narrowed to `except KeyError`.
- VCF-to-Zarr conversion errors are logged with `exc_info=True` instead of the message
  alone.
- A rank/value length mismatch now raises a `ValueError` reporting both lengths, instead
  of only logging an error and continuing.

### Documentation

- README workflow graph corrected.

## [1.3.0] - 2026-04-17

The release used for the analyses in the ExP Heatmap manuscript (commit `996c053`).

### Added

- Ordered-pair directional ranking and the reciprocal-direction-aware empirical
  rank-score framing.
- User-defined population panels through the `populations=` API parameter, for
  non-canonical population layouts. (The `--populations` CLI flag arrived later, in
  1.3.2.)
- Pair-specific missing-data handling.
- Coordinate-aware, max-preserving column reducer for dense regions.
- MIT license file; heatmap gallery assets.

## Earlier releases

No changelog was kept before 1.3.2; the entries above were reconstructed from the git
history. Earlier releases are recorded on PyPI:

| Version | Released on PyPI |
| ------- | ---------------- |
| 1.2.0   | 2026-02-25       |
| 1.1.5   | 2025-10-21       |
| 1.1.2a0 | 2024-04-02       |

Their contents are not reconstructed here, because the git tags do not delimit them:

- `v1.2.0` points at commit `6f95e2c` (2026-05-19), two commits *after* `v1.3.1` and three
  months after the 1.2.0 release went to PyPI. The tag was created against the wrong
  commit and cannot be used to recover the 1.2.0 range.
- `Version` points at commit `c068517` (2026-02-04) and carries no version number.

A prototype, v0.9.0, was disseminated separately as a public software resource and is not
covered here.

[1.3.2]: https://github.com/bioinfocz/exp_heatmap/compare/v1.3.1...v1.3.2
[1.3.1]: https://github.com/bioinfocz/exp_heatmap/compare/v1.3.0...v1.3.1
[1.3.0]: https://github.com/bioinfocz/exp_heatmap/releases/tag/v1.3.0
