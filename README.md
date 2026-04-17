# ExP Heatmap

[![PyPI version](https://badge.fury.io/py/exp-heatmap.svg)](https://pypi.org/project/exp_heatmap/)
![Python](https://img.shields.io/badge/python-≥3.10-blue)
![Version](https://img.shields.io/badge/version-1.3.0-green)
![License](https://img.shields.io/badge/license-MIT-green)
[![Tests](https://github.com/bioinfocz/exp_heatmap/actions/workflows/tests.yml/badge.svg)](https://github.com/bioinfocz/exp_heatmap/actions/workflows/tests.yml)

> An ordered-pair workflow for regional cross-population genomic visualization. ExP Heatmap takes a VCF through filter → prepare → compute → plot and renders ordered population-pair matrices as empirical rank-score heatmaps.

ExP Heatmap is designed for visualizing cross-population selection signals, differentiation, and other statistical summaries across many populations at once. It is most useful when a regional view (e.g. a gene locus or an ordered-pair matrix over a candidate window) would be hard to interpret as dozens of separate per-pair tracks. The package supports the canonical 1000 Genomes Project panel out of the box and now also handles arbitrary custom panels (e.g. the Gambian Genome Variation Project).

<img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/heatmap_gallery/LCT_gene.png" width="800" alt="ExP heatmap of LCT gene">

*ExP heatmap of the human lactose (LCT) gene showing population differences between 26 populations from the 1000 Genomes Project, displaying empirical rank scores for cross-population extended haplotype homozygosity (XPEHH) selection test. Create your own LCT heatmap with the [Quick Start](#quick-start) Guide*

**Developed by the [Laboratory of Genomics and Bioinformatics](https://www.img.cas.cz/group/michal-kolar/), Institute of Molecular Genetics of the Academy of Sciences of the Czech Republic**

## Features

- **End-to-end workflow**: Built-in `filter-vcf` → `prepare` → `compute` → `plot` pipeline with no external preprocessing dependencies required
- **Multiple statistical tests**: XP-EHH, XP-nSL, Delta Tajima's D, and Hudson's Fst
- **Flexible input**: Work from raw VCF, pre-computed statistics, or ready-to-plot TSVs
- **Custom population panels**: Automatic detection of non-1000G panels (e.g. GGVP) with fallback to the canonical 1000G layout when the full 26-population set is present
- **Correct pair-specific missingness**: A NaN in one population pair no longer silently removes that locus from all other pairs
- **Wide-region static rendering**: Explicit, documented column downsampling (`max`/`mean`/`median`) instead of implicit raster compression
- **Interactive HTML**: Plotly-based zoom, pan, hover, region comparison, and population-focused views
- **Reproducibility scaffolding**: Conda environment, Docker image, and CI-tested pytest suite
- **Benchmark harness**: Scripts for local display, full-pipeline, and population-scaling benchmarks

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Pipeline Overview](#pipeline-overview)
  - [Command-Line Interface](#command-line-interface)
  - [Python Package](#python-package)
- [Reproducibility](#reproducibility)
- [Input File Formats](#input-file-formats)
- [Output File Formats](#output-file-formats)
- [Advanced Features](#advanced-features)
- [Workflow Examples](#workflow-examples)
- [Gallery](#gallery)
- [Statistical Methodology](#statistical-methodology)
- [Benchmarks](#benchmarks)
- [Validation on Public Data](#validation-on-public-data)
- [1000 Genomes Population Reference](#1000-genomes-population-reference)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Requirements

- Python ≥ 3.10 (tested on 3.10, 3.11, 3.12)
- `vcftools` (optional; only needed if you prefer external SNP-only preprocessing)

ExP Heatmap now includes a built-in `filter-vcf` command for preparing biallelic SNP-only VCFs locally, so `vcftools` is no longer required for the standard preprocessing path.

### Python Dependencies

ExP Heatmap requires the following Python packages (automatically installed):

| Package | Version | Purpose |
|---------|---------|---------|
| scikit-allel | latest | Population genetics computations |
| zarr | **< 3.0.0** | Efficient array storage (v3+ not supported) |
| numpy | latest | Numerical operations |
| pandas | latest | Data manipulation |
| matplotlib | latest | Static visualizations |
| seaborn | latest | Heatmap rendering |
| click | latest | Command-line interface |
| plotly | latest | Interactive visualizations |
| tqdm | latest | Progress bars |

> **Important**: ExP Heatmap requires `zarr < 3.0.0`. The package will automatically install a compatible version, but if you encounter issues, ensure you have the correct version:
> ```bash
> pip install 'zarr<3.0.0'
> ```

### Install from PyPI

```bash
pip install exp_heatmap
```

### Install from GitHub (latest version)

```bash
pip install git+https://github.com/bioinfocz/exp_heatmap.git
```

## Quick Start

Get started with ExP Heatmap in three simple steps:

**Step 1**: Download the prepared results of the extended haplotype homozygosity (XPEHH) selection test for the part of human chromosome 2, 1000 Genomes Project data either directly via [Zenodo](https://zenodo.org/records/16364351) or via command:
```bash
wget "https://zenodo.org/records/16364351/files/chr2_output.tar.gz"
```
**Step 2**: Decompress the downloaded folder in your working directory: 
```bash
tar -xzf chr2_output.tar.gz
```
**Step 3**: Run the exp_heatmap plot command:
```bash
exp_heatmap plot chr2_output/ --start 136108646 --end 137108646 --title "LCT gene" --out LCT_xpehh
```
The `exp_heatmap` package will read the files from `chr2_output/` folder and create the ExP heatmap and save it as `LCT_xpehh.png` file.

<br/>

## Usage

### Pipeline Overview

ExP Heatmap follows a simple three-step workflow: **prepare** → **compute** → **plot**. If your input VCF still contains indels or multiallelic sites, use **filter-vcf** first. Each step can be used independently depending on your data format.

```
┌─────────────┐      ┌─────────────┐      ┌─────────────┐      ┌─────────────┐
│   VCF File  │ ──── │ filter-vcf  │ ──── │ SNP-only VCF│ ──── │   prepare   │
└─────────────┘      └─────────────┘      └─────────────┘      └──────┬──────┘
                                                                        │
                                                                        ▼
                                                                 ┌─────────────┐
                                                                 │   ZARR Dir  │
                                                                 └──────┬──────┘
                                                                        │
                                                                        ▼
┌─────────────┐      ┌─────────────┐      ┌─────────────┐      ┌─────────────┐
│ Panel File  │ ──── │   compute   │ ──── │  TSV Files  │ ──── │   Heatmap   │
└─────────────┘      └─────────────┘      └──────┬──────┘      │ PNG / HTML  │
                                                 │             └─────────────┘
                                                 ▼
                                          ┌─────────────┐
                                          │    plot     │
                                          └─────────────┘
```

**Starting Points:**
- **From raw VCF**: Use `filter-vcf` first if needed, then `prepare` → `compute` → `plot`
- **From VCF**: Use all three steps (`prepare` → `compute` → `plot`)
- **From ZARR**: Skip `prepare`, use `compute` → `plot`
- **From TSV results**: Skip to `plot` directly

### Command-Line Interface

#### 1. VCF Filtering - `filter-vcf`

>Keep only biallelic SNP records before `prepare`.

```bash
exp_heatmap filter-vcf [OPTIONS] <input_vcf>
```

| Argument/Option | Type | Default | Description |
|-----------------|------|---------|-------------|
| `<input_vcf>` | PATH | required | Input VCF or VCF.GZ file |
| `-o, --out` | PATH | required | Output VCF path (`.vcf` or `.vcf.gz`) |
| `--no-log` | flag | - | Disable logging to file |
| `--verbose` | flag | - | Show detailed debug output in console |

**Example:**
```bash
exp_heatmap filter-vcf chr21.vcf.gz -o chr21_snps.vcf
```

#### 2. Full Pipeline - `full`

>Run the complete pipeline (prepare → compute → plot) in a single command.

```bash
exp_heatmap full [OPTIONS] <vcf_file> <panel_file>
```

| Argument/Option | Type | Default | Description |
|-----------------|------|---------|-------------|
| `<vcf_file>` | PATH | required | Recoded VCF file (SNPs only recommended) |
| `<panel_file>` | PATH | required | Population panel file |
| `-o, --out` | PATH | `exp_heatmap` | Prefix for all output files |
| `-s, --start` | INT | required | Start position for displayed region |
| `-e, --end` | INT | required | End position for displayed region |
| `-t, --test` | choice | `xpehh` | Statistical test: `xpehh`, `xpnsl`, `delta_tajima_d`, `hudson_fst` |
| `-c, --chunked` | flag | - | Use chunked array to avoid memory exhaustion |
| `--title` | STR | - | Title of the heatmap |
| `--cmap` | STR | `Blues` | Colormap for visualization |
| `--interactive` | flag | - | Generate interactive HTML visualization |
| `--max-columns` | INT | auto/static, `30000` interactive | Explicit column budget for wide regions |
| `--column-aggregation` | choice | `max` | Reducer for static downsampling: `max`, `mean`, `median` |
| `--dpi` | INT | `400` | DPI for saved static figures |
| `--no-log` | flag | - | Disable logging to file |
| `--verbose` | flag | - | Show detailed debug output |

**Example:**
```bash
exp_heatmap full chr15_snps.recode.vcf genotypes.panel -s 47924019 -e 48924019 -o slc24a5_analysis --title "SLC24A5"
```

This creates: `slc24a5_analysis_zarr/`, `slc24a5_analysis_compute/`, and `slc24a5_analysis_plot.png`

#### 3. Data Preparation - `prepare`

>Convert VCF files to efficient Zarr format for faster computation.

```bash
exp_heatmap prepare [OPTIONS] <vcf_file>
```

| Argument/Option | Type | Default | Description |
|-----------------|------|---------|-------------|
| `<vcf_file>` | PATH | required | VCF file (biallelic SNP-only recommended) |
| `-o, --out` | PATH | `zarr_output` | Directory for ZARR output files |
| `--no-log` | flag | - | Disable logging to file |
| `--verbose` | flag | - | Show detailed debug output in console |

**Example:**
```bash
exp_heatmap prepare chr15_snps.recode.vcf -o chr15.zarr
```

#### 4. Statistical Analysis - `compute`

>Calculate population genetic statistics across all genomic positions.

```bash
exp_heatmap compute [OPTIONS] <zarr_dir> <panel_file>
```

| Argument/Option | Type | Default | Description |
|-----------------|------|---------|-------------|
| `<zarr_dir>` | PATH | required | Directory with ZARR files from `prepare` step |
| `<panel_file>` | PATH | required | Population panel file (see [Input File Formats](#input-file-formats)) |
| `-o, --out` | PATH | `output` | Directory for output TSV files |
| `-t, --test` | choice | `xpehh` | Statistical test: `xpehh`, `xpnsl`, `delta_tajima_d`, `hudson_fst` |
| `-c, --chunked` | flag | - | Use chunked array to avoid memory exhaustion |
| `--no-log` | flag | - | Disable logging to file |
| `--verbose` | flag | - | Show detailed debug output in console |

**Statistical Tests:**
- `xpehh`: Cross-population Extended Haplotype Homozygosity - detects recent positive selection
- `xpnsl`: Cross-population Number of Segregating sites by Length - robust to variation in recombination rate
- `delta_tajima_d`: Delta Tajima's D - measures difference in allele frequency spectrum
- `hudson_fst`: Hudson's Fst - genetic differentiation between populations

> **Note**: The `-t` flag has different meanings for different commands: for `compute` it specifies the statistical test (`--test`), while for `plot` it specifies the heatmap title (`--title`).

**Example:**
```bash
exp_heatmap compute chr15.zarr genotypes.panel -o chr15_results -t xpehh
```

#### 5. Visualization - `plot`

>Generate heatmap visualizations from computed statistics.

```bash
exp_heatmap plot [OPTIONS] <input_dir>
```

| Argument/Option | Type | Default | Description |
|-----------------|------|---------|-------------|
| `<input_dir>` | PATH | required | Directory with TSV files from `compute` step |
| `-s, --start` | INT | required | Start genomic position |
| `-e, --end` | INT | required | End genomic position |
| `-t, --title` | STR | - | Title of the heatmap |
| `-o, --out` | PATH | `ExP_heatmap` | Output filename (without extension) |
| `-c, --cmap` | STR | `Blues` | Colormap name (see [Colormap Options](#colormap-options)) |
| `--interactive` | flag | - | Generate interactive HTML visualization |
| `--max-columns` | INT | auto/static, `30000` interactive | Explicit column budget for wide regions |
| `--column-aggregation` | choice | `max` | Reducer for static downsampling: `max`, `mean`, `median` |
| `--dpi` | INT | `400` | DPI for saved static figures |
| `--no-log` | flag | - | Disable logging to file |
| `--verbose` | flag | - | Show detailed debug output in console |

When the input TSVs represent the full canonical 1000 Genomes panel, ExP Heatmap uses the standard 1000G row ordering and annotations. For other datasets, the plotter now infers the population set directly from the TSV filenames and renders a custom-population heatmap automatically.

**Example:**
```bash
# Basic usage
exp_heatmap plot chr15_results/ --start 47924019 --end 48924019 --title "SLC24A5" --out slc24a5

# Interactive HTML output
exp_heatmap plot chr15_results/ --start 47924019 --end 48924019 --interactive --out slc24a5_interactive

# Static output with explicit wide-region downsampling
exp_heatmap plot chr15_results/ --start 47924019 --end 48924019 --max-columns 6000 --column-aggregation max --out slc24a5_binned
```

#### 6. Superpopulation Summary - `summary`

>Collapse population-pair rows to superpopulation pairs, write the collapsed TSV, and optionally create an interactive summary heatmap.

```bash
exp_heatmap summary [OPTIONS] <input_dir>
```

**Example:**
```bash
exp_heatmap summary chr15_results/ --start 47924019 --end 48924019 --agg-func mean --out slc24a5_summary
```

This creates:
- `slc24a5_summary.tsv`
- `slc24a5_summary.html` unless `--no-plot` is used

#### 7. Population-Focused View - `focus`

>Generate an interactive view containing only rows involving one chosen population.

```bash
exp_heatmap focus [OPTIONS] <input_dir>
```

**Example:**
```bash
exp_heatmap focus chr15_results/ --start 47924019 --end 48924019 --population CEU --out slc24a5_ceu_focus
```

#### 8. Region Comparison - `compare`

>Generate an interactive side-by-side comparison of two genomic regions from the same computed result set.

```bash
exp_heatmap compare [OPTIONS] <input_dir>
```

**Example:**
```bash
exp_heatmap compare chr15_results/ --start-1 47924019 --end-1 48924019 --start-2 55000000 --end-2 57000000 --out region_comparison
```

#### 9. Top-Region Extraction - `regions`

>Extract top-scoring windows into a TSV.

```bash
exp_heatmap regions [OPTIONS] <input_dir>
```

**Example:**
```bash
exp_heatmap regions chr15_results/ --start 47924019 --end 48924019 --n-top 25 --out slc24a5_top_regions.tsv
```

---

### Python Package

The Python API offers more flexibility and customization options. Choose the appropriate scenario based on your data format:

#### Scenario A: Ready-to-Plot Data

**Use when:** You have pre-computed rank scores in a TSV file.

**Data format:** TSV file with columns: `CHROM`, `POS`, followed by pairwise columns for population comparisons.

```python
from exp_heatmap.plot import plot_exp_heatmap
import pandas as pd

# Load your data
data = pd.read_csv("rank_scores.tsv", sep="\t")

# Create heatmap
plot_exp_heatmap(
    data,
    start=135287850,
    end=136287850,
    title="Population Differences in LCT Gene",
    cmap="Blues",
    output="lct_analysis",
    populations="1000Genomes"  # Predefined population set
)
```

#### Scenario B: Statistical Results to Rank Scores

**Use when:** You have computed statistical test results that need conversion to rank scores.

```python
from exp_heatmap.plot import plot_exp_heatmap, create_plot_input

# Convert statistical results to empirical rank scores
data_to_plot = create_plot_input(
    "results_directory/",      # Directory with test results
    start=135287850, 
    end=136287850, 
    populations="1000Genomes",
    rank_scores="directional"  # Options: "directional" (legacy alias: "2-tailed"), "ascending", "descending"
)

# Create heatmap
plot_exp_heatmap(
    data_to_plot,
    start=135287850,
    end=136287850,
    title="XP-NSL Test Results",
    cmap="expheatmap",         # Custom ExP colormap
    output="xpnsl_results"
)
```

#### Scenario C: Complete VCF Workflow

**Use when:** Starting from raw VCF files. Combine CLI commands with Python plotting:

```python
import subprocess
from exp_heatmap.plot import plot_exp_heatmap, create_plot_input

# 1. Prepare data (using CLI)
subprocess.run(["exp_heatmap", "prepare", "data_snps.recode.vcf", "-o", "data.zarr"])

# 2. Compute statistics (using CLI) 
subprocess.run(["exp_heatmap", "compute", "data.zarr", "populations.panel", "-o", "results/"])

# 3. Create custom plots (using Python)
data_to_plot = create_plot_input("results/", start=47000000, end=49000000)
plot_exp_heatmap(data_to_plot, start=47000000, end=49000000, 
                 title="Custom Analysis", output="custom_plot")
```

#### Additional Python API Functions

**Summarize by Superpopulation:**
```python
from exp_heatmap.plot import create_plot_input, summarize_by_superpopulation

# Load data
data = create_plot_input("results/", start=47000000, end=49000000)

# Aggregate to superpopulation level (AFR, EUR, EAS, SAS, AMR)
superpop_data = summarize_by_superpopulation(data, agg_func='mean')
# Result has 20 rows (5×4 superpopulation pairs) instead of 650
```

**Extract Top Regions:**
```python
from exp_heatmap.plot import create_plot_input, extract_top_regions

# Load data
data = create_plot_input("results/", start=47000000, end=49000000)

# Find genomic windows with highest selection signals
top_regions = extract_top_regions(data, n_top=50, window_size=10000)
print(top_regions[['center', 'mean_score', 'top_population_pair']])
```

**Custom Colorbar Parameters:**
```python
from exp_heatmap.plot import prepare_cbar_params

# Calculate optimal colorbar settings based on data range
cmin, cmax, cbar_ticks = prepare_cbar_params(data_to_plot, n_cbar_ticks=6)
```

---

## Reproducibility

Two local reproducibility entry points are included:

- `environment.yml`: conda environment with Python 3.11, optional `vcftools`, and the package installed in editable mode with the `dev` extra
- `Dockerfile`: minimal container image that installs the CLI and optional external preprocessing tools

Quick start:

```bash
conda env create -f environment.yml
conda activate exp-heatmap-dev
```

Or:

```bash
docker build -t exp-heatmap-local .
docker run --rm exp-heatmap-local --help
```

---

## Input File Formats

### VCF File

Standard VCF format. For best results:
- Filter to biallelic SNPs only before `prepare`
- Use the built-in CLI filter or your preferred external normalization tool

```bash
exp_heatmap filter-vcf input.vcf.gz -o input_snps.vcf
```

### Panel File

Tab-separated file defining population membership for each sample. **Required columns:**

| Column | Description |
|--------|-------------|
| `sample` | Sample identifier (must match VCF sample names exactly) |
| `pop` | Population code (e.g., "CEU", "YRI") |
| `super_pop` | Superpopulation code (e.g., "EUR", "AFR") |

**Example panel file:**
```
sample	pop	super_pop	gender
HG00096	GBR	EUR	male
HG00097	GBR	EUR	female
HG00099	GBR	EUR	female
NA18486	YRI	AFR	male
NA18487	YRI	AFR	female
NA18488	YRI	AFR	male
```

> **Important**: The sample order in the panel file must match the sample order in the VCF/ZARR file exactly.

**1000 Genomes Panel File:**

Download the official panel file:
```bash
wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" -O genotypes.panel
```

---

## Output File Formats

### Compute Output (TSV Files)

The `compute` step generates one TSV file per population pair, named `POP1_POP2.tsv`.

Rows remain position-aligned across all pair files. If one population pair yields no valid value at a locus, that locus is preserved as `NaN` for that pair rather than being globally removed from every output file.

**Columns:**

| Column | Description |
|--------|-------------|
| `name` | Dataset identifier (derived from input filename) |
| `variant_pos` | Genomic position of the variant |
| `{test}` | Raw test statistic value (e.g., `xpehh`, `xpnsl`) |
| `-log10_p_value_ascending` | Empirical rank score (ascending sort) |
| `-log10_p_value_descending` | Empirical rank score (descending sort) |

**Example output (`CEU_YRI.tsv`):**
```
name	variant_pos	xpehh	-log10_p_value_ascending	-log10_p_value_descending
chr15_snps	48120990	0.523	1.234	2.456
chr15_snps	48120995	-0.891	2.891	1.123
chr15_snps	48121003	1.245	0.891	3.012
```

> **Note**: Column names contain "p_value" for backward compatibility, but these are **empirical rank scores**, not classical p-values. See [Statistical Methodology](#statistical-methodology) for details.

---

## Advanced Features

### Interactive Visualizations

Generate HTML-based interactive heatmaps with zoom, pan, and hover tooltips:

```python
from exp_heatmap.interactive import plot_interactive

# Create interactive HTML visualization
plot_interactive(
    "results_directory/",
    start=135287850,
    end=136287850,
    title="Interactive LCT Analysis",
    output="lct_interactive"  # Saves as lct_interactive.html
)
```

Or via CLI:
```bash
exp_heatmap plot results/ --start 135287850 --end 136287850 --interactive --out lct_interactive
```

For static PNG output, ExP Heatmap now supports explicit wide-region downsampling. When the genomic window is wider than the effective export width, the static renderer can aggregate neighboring SNP columns using `max`, `mean`, or `median` before drawing the heatmap. The default static behavior uses an automatic column budget derived from figure width and DPI.

**Additional Interactive Functions:**

```python
from exp_heatmap.interactive import create_comparison_view, create_population_focus_view
from exp_heatmap.plot import create_plot_input

data = create_plot_input("results/", start=40000000, end=60000000)

# Compare two genomic regions side-by-side
create_comparison_view(
    data,
    region1=(47000000, 49000000),
    region2=(55000000, 57000000),
    title="Region Comparison",
    output="comparison"
)

# Focus on comparisons involving a specific population
create_population_focus_view(
    data,
    focus_population="CEU",
    start=47000000,
    end=49000000,
    title="CEU Selection Signals",
    output="ceu_focus"
)
```

### Advanced Customization

Fine-tune your visualizations with advanced options:

```python
from exp_heatmap.plot import plot_exp_heatmap, prepare_cbar_params, superpopulations

# Custom colorbar settings
cmin, cmax, cbar_ticks = prepare_cbar_params(data_to_plot, n_cbar_ticks=6)

# Advanced plot with multiple customizations
plot_exp_heatmap(
    data_to_plot,
    start=135000000,
    end=137000000,
    title="Selection Signals in African Populations",
    
    # Population filtering
    populations=superpopulations["AFR"],  # Focus on African populations
    # Available: ["AFR", "AMR", "EAS", "EUR", "SAS"] or custom list
    
    # Visual customizations
    cmap="expheatmap",                    # Custom ExP colormap
    display_limit=1.60,                   # Filter noise (values below limit = white)
    display_values="higher",              # Show values above display_limit
    
    # Annotations
    vertical_line=[                       # Mark important SNPs
        [135851073, "rs41525747"],        # [position, label]
        [135851081, "rs41380347"]
    ],
    
    # Colorbar customization
    cbar_vmin=cmin,
    cbar_vmax=cmax,
    cbar_ticks=cbar_ticks,
    
    # Output
    output="african_populations_analysis",
    xlabel="Custom region description"
)
```

### Colormap Options

ExP Heatmap supports all matplotlib colormaps plus a custom colormap:

| Colormap | Description | Best For |
|----------|-------------|----------|
| `Blues` | Sequential blue gradient (default) | General use |
| `expheatmap` | Custom colormap based on `gist_ncar_r` with white background | Highlighting strong signals |
| `gist_heat` | Heat-style gradient | High contrast visualization |
| `Reds` | Sequential red gradient | Alternative sequential |
| `viridis` | Perceptually uniform | Color-blind friendly |
| `RdBu` | Diverging red-blue | Bidirectional data |

See the [full list of matplotlib colormaps](https://matplotlib.org/stable/users/explain/colors/colormaps.html).

**Custom `expheatmap` colormap:**
- Based on `gist_ncar_r` 
- White background for low values (better noise filtering)
- Dark blue for highest values
- Optimized for selection signal visualization

## Workflow Examples

### Complete Analysis: SLC24A5 Gene

<img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/heatmap_gallery/SLC24A5_gene.png" width="800" alt="ExP heatmap of SLC24A5 gene">

This example demonstrates a full workflow analyzing the [SLC24A5](https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000188467;r=15:48120990-48142672) gene, known for its role in human skin pigmentation using 1000 Genomes Project data. SLC24A5 is also known to show strong selection signals, which makes it a suitable example.

```bash
#!/bin/bash

# Download 1000 Genomes data
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" -O chr15.vcf.gz
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" -O genotypes.panel

# Filter to biallelic SNPs only
exp_heatmap filter-vcf chr15.vcf.gz -o chr15_snps.vcf

# Prepare data
exp_heatmap prepare chr15_snps.vcf -o chr15_snps.zarr

# Compute statistics
exp_heatmap compute chr15_snps.zarr genotypes.panel -o chr15_snps_output

# Generate heatmap for SLC24A5 region
exp_heatmap plot chr15_snps_output \
    --start 47924019 \
    --end 48924019 \
    --title "SLC24A5" \
    --cmap gist_heat \
    --out SLC24A5_heatmap
```

## Gallery

### Different Rank Score Computations

The same XP-EHH test data for the ADM2 gene region, showing different rank score calculation methods:

**Directional rank scores** *(reciprocal ordered-pair ranking; legacy alias: "2-tailed")*:
<img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/heatmap_gallery/ADM2_2tailed.png" width="800" alt="Directional rank scores">

**Ascending rank scores** *(one-sided, highlights lowest test values)*:
<img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/heatmap_gallery/ADM2_ascending.png" width="800" alt="Ascending rank scores">

**Descending rank scores** *(one-sided, highlights highest test values)*:
<img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/heatmap_gallery/ADM2_descending.png" width="800" alt="Descending rank scores">

### Noise Filtering

Using `display_limit` and `display_values` parameters to filter noisy data and highlight the most extreme regions:

<img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/heatmap_gallery/ADM2_display_limit.png" width="800" alt="Filtered display">

*Same data as above, but with display_limit=1.60 to filter noise and highlight the strongest signals.*

## Statistical Methodology

### Empirical Rank Scores

ExP Heatmap uses **empirical rank scores** to visualize selection signals. Despite column names containing "p_value" for backward compatibility, these values are **not inferential p-values**.

**How rank scores are computed:**

1. **Rank all variants**: For each population pair, sort all genome-wide test statistics.
2. **Convert the rank to a fraction**: `rank_fraction = rank / total_variants`
3. **Transform to the display scale**: `empirical_rank_score = -log10(rank_fraction)`

**Interpretation:**
- Higher values indicate more extreme variants within that empirical ranking.
- A value of 3.0 means the variant falls in the top 0.1% of the ranked distribution.
- A value of 2.0 means the variant falls in the top 1% of the ranked distribution.

### Rank Score Options

The `rank_scores` parameter in `create_plot_input()` controls how rank scores are calculated:

| Option | Description | Use Case |
|--------|-------------|----------|
| `"directional"` | For POP1_POP2: use descending; for POP2_POP1: use ascending | **Default** - Captures reciprocal directionality |
| `"2-tailed"` | Legacy alias of `"directional"` | Backward-compatible legacy name |
| `"ascending"` | Lowest test values ranked first | Detect negative selection signals |
| `"descending"` | Highest test values ranked first | Detect positive selection signals |

**Example:**
```python
# Direction-aware reciprocal ranking (recommended for most analyses)
data = create_plot_input("results/", start=47000000, end=49000000, rank_scores="directional")

# One-sided ranking for a specific hypothesis
data = create_plot_input("results/", start=47000000, end=49000000, rank_scores="descending")
```

### Statistical Tests Explained

| Test | Detects | Interpretation |
|------|---------|----------------|
| **XP-EHH** | Recent positive selection | Positive values: selection in pop1; Negative: selection in pop2 |
| **XP-nSL** | Selection (robust to recombination rate variation) | Similar to XP-EHH but more robust |
| **Delta Tajima's D** | Difference in allele frequency spectrum | Positive: excess rare variants in pop1 |
| **Hudson's Fst** | Population differentiation | Higher values: greater genetic distance |

## Benchmarks

All benchmarks were run locally on a macOS 14 / arm64 (Apple M4 Pro) workstation with 14 logical CPUs and 48 GB RAM. Python 3.12, zarr 2.x. Full benchmark scripts are provided under [`scripts/benchmarks/`](scripts/benchmarks).

### Full pipeline on GGVP chr21

Input: Gambian Genome Variation Project integrated chromosome 21 VCF (553,906 input records), filtered to 509,924 biallelic SNPs across 505 samples and 5 population labels (`FULA`, `GWD`, `JOLA`, `MANDINKA`, `WOLOFF`). Statistic: XP-EHH. Default settings otherwise.

| Stage | Small case (1 Mb, 11,620 SNPs) | Full chromosome (509,924 SNPs) |
|-------|-------------------------------:|-------------------------------:|
| `prepare` | 0.99 s  /  0.28 GB peak RSS | 12.64 s  /  0.38 GB peak RSS |
| `compute` | 2.76 s  /  0.21 GB peak RSS | 74.73 s  /  0.99 GB peak RSS |
| `plot` (static overview) | 2.34 s  /  0.26 GB peak RSS | 2.27 s  /  0.74 GB peak RSS |

`compute` dominates wall-clock cost at full-chromosome scale, while static overview plotting remains cheap once per-pair TSVs exist. Peak RSS stayed below 1 GB for every stage in both cases.

### Display scaling with population count

Synthetic scaling test that holds the genomic window fixed and varies the number of ordered population-pair rows by symlinking existing per-pair outputs. All runs use the GGVP chr21 compute output as the seed.

| Populations | Ordered rows | Static plot | Interactive HTML | Peak RSS (interactive) |
|-------------|-------------:|------------:|-----------------:|-----------------------:|
| 5 | 20 | 1.80 s, 0.10 MB PNG | 1.83 s, 5.4 MB HTML | 0.43 GB |
| 10 | 90 | 3.25 s | 3.57 s, 10.3 MB HTML | 0.62 GB |
| 20 | 380 | 9.61 s | 9.24 s, 32.4 MB HTML | 1.08 GB |
| 30 | 870 | 16.87 s | 16.34 s, 54.1 MB HTML | 1.48 GB |
| 50 | 2,450 | 38.36 s, 2.47 MB PNG | 36.09 s, 99.5 MB HTML | 2.17 GB |

The ordered-pair display stays computationally tractable up to 50 populations, but visual density and HTML size rise quickly. At 50 populations we recommend using the `focus`, `summary`, or `regions` workflows in addition to the full matrix view.

### Reproducing the benchmarks

```bash
pip install -e .[benchmarks]

# Full pipeline benchmark on your own VCF
python scripts/benchmarks/run_pipeline_benchmark.py \
    --vcf path/to/input.vcf.gz \
    --panel path/to/panel.tsv \
    --out-dir local_data/benchmarks/pipeline

# Population scaling benchmark (seeds from an existing compute output directory)
python scripts/benchmarks/run_population_scaling.py \
    --compute-dir path/to/compute_output \
    --out-dir local_data/benchmarks/population_scaling
```

## Validation on Public Data

ExP Heatmap ships with two reproducible 1000 Genomes Project validation pipelines under [`scripts/validation/`](scripts/validation). Both scripts start from public Phase 3 release URLs and run end-to-end through `filter-vcf` → `prepare` → `compute` → `plot`.

### chr15 / SLC24A5 (full-pipeline validation)

The `run_1kg_chr15_slc24a5.sh` script reproduces the pigmentation-locus showcase using raw public 1000 Genomes chromosome 15 data:

| Stage | Wall-clock |
|-------|-----------:|
| `filter-vcf` | 99.59 s (biallelic-SNP filter; retains 6,456,568 / 6,477,157 records) |
| `prepare` | 252.98 s (VCF → Zarr) |
| `compute` | 3034.42 s (XP-EHH, all 650 ordered population pairs from 26 populations) |
| `plot` (static) | 28.43 s |
| `plot` (interactive) | 14.04 s |

Output artifacts include a 21.95 GB filtered VCF, 499.69 MB Zarr store, 2.94 GB of pairwise TSV results, and a 442 KB PNG / 21.54 MB HTML heatmap of the SLC24A5 window (chr15:47,924,019-48,924,019).

### chr2 / LCT locus (region-scoped reconstruction)

The `run_1kg_chr2_lct_reconstruction.sh` script reconstructs the canonical lactase-persistence locus as a region-scoped public-data run.The script does not try to re-run whole-chromosome compute; instead it filters chromosome 2 to the plotted LCT window plus 1 Mb of flanking sequence on each side before preparing and computing.

| Stage | Wall-clock |
|-------|-----------:|
| `filter-vcf` (chr2:135,108,646-138,108,646) | 268.18 s (retains 81,595 / 7,081,600 records) |
| `prepare` | 11.18 s |
| `compute` | 117.41 s |
| `plot` (static) | 4.61 s |
| `plot` (interactive) | 2.11 s |

This scales the whole-locus reconstruction to ~7 minutes of wall time end-to-end while keeping the plotted window and interpretation identical. See [`scripts/validation/summarize_validation_run.py`](scripts/validation/summarize_validation_run.py) to regenerate the JSON/Markdown summaries from log files.

## 1000 Genomes Population Reference

ExP Heatmap is optimized for the 1000 Genomes Project Phase 3 data (26 populations).

### Population Codes

| Code | Population | Superpopulation |
|------|------------|-----------------|
| **ACB** | African Caribbean in Barbados | AFR |
| **ASW** | African Ancestry in SW USA | AFR |
| **ESN** | Esan in Nigeria | AFR |
| **GWD** | Gambian in Western Division | AFR |
| **LWK** | Luhya in Webuye, Kenya | AFR |
| **MSL** | Mende in Sierra Leone | AFR |
| **YRI** | Yoruba in Ibadan, Nigeria | AFR |
| **BEB** | Bengali in Bangladesh | SAS |
| **GIH** | Gujarati Indians in Houston | SAS |
| **ITU** | Indian Telugu in the UK | SAS |
| **PJL** | Punjabi in Lahore, Pakistan | SAS |
| **STU** | Sri Lankan Tamil in the UK | SAS |
| **CDX** | Chinese Dai in Xishuangbanna | EAS |
| **CHB** | Han Chinese in Beijing | EAS |
| **CHS** | Han Chinese South | EAS |
| **JPT** | Japanese in Tokyo | EAS |
| **KHV** | Kinh in Ho Chi Minh City, Vietnam | EAS |
| **CEU** | Utah residents (CEPH) with European ancestry | EUR |
| **FIN** | Finnish in Finland | EUR |
| **GBR** | British in England and Scotland | EUR |
| **IBS** | Iberian populations in Spain | EUR |
| **TSI** | Toscani in Italy | EUR |
| **CLM** | Colombian in Medellín | AMR |
| **MXL** | Mexican Ancestry in Los Angeles | AMR |
| **PEL** | Peruvian in Lima | AMR |
| **PUR** | Puerto Rican in Puerto Rico | AMR |

### Superpopulations

| Code | Name | Populations |
|------|------|-------------|
| **AFR** | African | ACB, ASW, ESN, GWD, LWK, MSL, YRI |
| **SAS** | South Asian | BEB, GIH, ITU, PJL, STU |
| **EAS** | East Asian | CDX, CHB, CHS, JPT, KHV |
| **EUR** | European | CEU, FIN, GBR, IBS, TSI |
| **AMR** | American (Admixed) | CLM, MXL, PEL, PUR |

**Reference:** [1000 Genomes Project](https://www.internationalgenome.org/category/population/)

## Troubleshooting

### Common Errors and Solutions

#### Zarr Version Error

**Error:**
```
Unsupported zarr version: 3.x.x
Please downgrade to zarr version < 3.0.0
```

**Solution:**
```bash
pip install 'zarr<3.0.0'
```

#### No Data Found in Genomic Region

**Error:**
```
ValueError: No data found in the requested genomic region (X - Y). 
The data contains positions from A to B.
```

**Cause:** The specified `--start`/`--end` coordinates don't overlap with the data.

**Solution:** 
- Check that your coordinates match the chromosome in your data
- Use coordinates within the available range shown in the error message
- Verify you're using the correct output directory

#### Sample Order Mismatch

**Error:**
```
Sample order differs! Found X mismatches
```

**Cause:** The sample order in the panel file doesn't match the VCF/ZARR file.

**Solution:**
- Ensure you're using the correct panel file for your VCF
- Check that both files are from the same data release/phase
- Verify sample IDs match exactly (case-sensitive)

#### Memory Exhaustion During Compute

**Symptoms:** Process killed, system becomes unresponsive, or "MemoryError"

**Solution:** Use the `--chunked` flag to process data in smaller chunks:
```bash
exp_heatmap compute data.zarr panel.tsv -o output/ --chunked
```

#### Empty or All-NaN Results

**Symptoms:** One pair-specific TSV file contains only `NaN` values in the statistic and rank-score columns.

**Cause:** The statistical test produced no valid results for that specific population pair, often due to:
- Too few variants
- Insufficient allele frequency variation
- For `delta_tajima_d`: window size too large

**Solution:**
- Check your data has sufficient variants
- For `delta_tajima_d`, the default window size is 13 SNPs; ensure your data has enough variants
- Inspect logs to identify which population pair failed; other pairs will still be preserved in the output directory

### Getting Help

- Check the [GitHub Issues](https://github.com/bioinfocz/exp_heatmap/issues) for known problems
- Enable verbose logging with `--verbose` flag for detailed debug output
- Log files are saved automatically (disable with `--no-log`)

## Contributing

We welcome contributions! Feel free to contact us or submit issues or pull requests.

### Development Setup

```bash
git clone https://github.com/bioinfocz/exp_heatmap.git
cd exp_heatmap
pip install -e ".[dev]"
```

### Running the test suite

```bash
python -m pytest
```

The suite currently contains 28 tests covering rank-score generation with ties and missing data, pair-specific output preservation, AF fallback behavior, static downsampling, VCF filtering, custom-panel inference, GGVP metadata handling, and CLI wiring. GitHub Actions runs the suite on Python 3.10, 3.11, and 3.12 for every push and pull request.

### Building documentation-facing assets

Reproducibility scripts, benchmark drivers, and validation pipelines live under [`scripts/`](scripts).

## License

This project is licensed under the [MIT License](https://github.com/bioinfocz/exp_heatmap/blob/main/LICENSE).

## Contributors

- **Edvard Ehler** ([@EdaEhler](https://github.com/EdaEhler)) - Lead Developer & Corresponding Author
- **Adam Nógell** ([@AdamNogell](https://github.com/AdamNogell)) - Maintainer
- **Jan Pačes** ([@hpaces](https://github.com/hpaces)) - Developer
- **Mariana Komárková** ([@satrovam](https://github.com/satrovam)) - Developer
- **Ondřej Moravčík** ([@ondra-m](https://github.com/ondra-m)) - Original Developer

## Acknowledgments

<div align="center">

<a href="http://genomat.img.cas.cz">
  <img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/funding/genomat.png" width="100" alt="GenoMat">
</a>
&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://www.img.cas.cz/en">
  <img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/funding/img.png" width="100" alt="IMG CAS">
</a>
&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://www.elixir-czech.cz">
  <img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/funding/elixir.png" width="100" alt="ELIXIR">
</a>

</div>

---

*If you use ExP Heatmap in your research, please cite our paper [citation details will be added upon publication].*
