# ExP Heatmap

[![PyPI version](https://badge.fury.io/py/exp-heatmap.svg)](https://pypi.org/project/exp_heatmap/)
![Python](https://img.shields.io/badge/python-≥3.8-blue)
![Version](https://img.shields.io/badge/version-1.2.0-green)

> A powerful Python package and command-line tool for visualizing multidimensional population genetics data through intuitive heatmaps.

ExP Heatmap specializes in displaying cross-population data, including differences, similarities, p-values, and other statistical parameters between multiple groups or populations. This tool enables efficient evaluation of millions of statistical values in a single, comprehensive visualization.

<img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/LCT_gene.png" width="800" alt="ExP heatmap of LCT gene">

*ExP heatmap of the human lactose (LCT) gene showing population differences between 26 populations from the 1000 Genomes Project, displaying empirical rank scores for cross-population extended haplotype homozygosity (XPEHH) selection test. Create your own LCT heatmap with the [Quick Start](#quick-start) Guide*

**Developed by the [Laboratory of Genomics and Bioinformatics](https://www.img.cas.cz/group/michal-kolar/), Institute of Molecular Genetics of the Academy of Sciences of the Czech Republic**

## Features

- **Multiple Statistical Tests**: Support for XPEHH, XP-NSL, Delta Tajima's D, and Hudson's Fst
- **Flexible Input Formats**: Work with VCF files, pre-computed statistics, or ready-to-plot data
- **Command-Line Interface**: Easy-to-use CLI for standard workflows
- **Python API**: Full programmatic control for custom analyses
- **Efficient Processing**: Zarr-based data storage for fast computation
- **Customizable Visualization**: Multiple color schemes, resolution options, and display settings
- **Interactive Mode**: Plotly-based HTML visualizations with zoom, pan, and hover tooltips

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Pipeline Overview](#pipeline-overview)
  - [Command-Line Interface](#command-line-interface)
  - [Python Package](#python-package)
- [Input File Formats](#input-file-formats)
- [Output File Formats](#output-file-formats)
- [Advanced Features](#advanced-features)
- [Workflow Examples](#workflow-examples)
- [Gallery](#gallery)
- [Statistical Methodology](#statistical-methodology)
- [1000 Genomes Population Reference](#1000-genomes-population-reference)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Requirements

- Python ≥ 3.8
- `vcftools` (for genomic data preparation - optional if using preprocessed data)

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

ExP Heatmap follows a simple three-step workflow: **prepare** → **compute** → **plot**. Each step can be used independently depending on your data format.

```
┌─────────────┐      ┌─────────────┐      ┌─────────────┐      ┌─────────────┐
│   VCF File  │ ──── │   prepare   │ ──── │   ZARR Dir  │      │             │
└─────────────┘      └─────────────┘      └──────┬──────┘      │             │
                                                 │             │             │
                                                 ▼             │             │
┌─────────────┐      ┌─────────────┐      ┌─────────────┐      │   Heatmap   │
│ Panel File  │ ──── │   compute   │ ──── │  TSV Files  │ ──── │    (PNG)    │
└─────────────┘      └─────────────┘      └──────┬──────┘      │             │
                                                 │             │             │
                                                 ▼             │             │
                                          ┌─────────────┐      │             │
                                          │    plot     │ ──── │             │
                                          └─────────────┘      └─────────────┘
```

**Starting Points:**
- **From VCF**: Use all three steps (`prepare` → `compute` → `plot`)
- **From ZARR**: Skip `prepare`, use `compute` → `plot`
- **From TSV results**: Skip to `plot` directly

### Command-Line Interface

#### 1. Full Pipeline - `full`

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
| `--no-log` | flag | - | Disable logging to file |
| `--verbose` | flag | - | Show detailed debug output |

**Example:**
```bash
exp_heatmap full chr15_snps.recode.vcf genotypes.panel -s 47924019 -e 48924019 -o slc24a5_analysis --title "SLC24A5"
```

This creates: `slc24a5_analysis_zarr/`, `slc24a5_analysis_compute/`, and `slc24a5_analysis_plot.png`

#### 2. Data Preparation - `prepare`

>Convert VCF files to efficient Zarr format for faster computation.

```bash
exp_heatmap prepare [OPTIONS] <vcf_file>
```

| Argument/Option | Type | Default | Description |
|-----------------|------|---------|-------------|
| `<vcf_file>` | PATH | required | Recoded VCF file (SNPs only recommended) |
| `-o, --out` | PATH | `zarr_output` | Directory for ZARR output files |
| `--no-log` | flag | - | Disable logging to file |
| `--verbose` | flag | - | Show detailed debug output in console |

**Example:**
```bash
exp_heatmap prepare chr15_snps.recode.vcf -o chr15.zarr
```

#### 3. Statistical Analysis - `compute`

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

#### 4. Visualization - `plot`

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
| `--no-log` | flag | - | Disable logging to file |
| `--verbose` | flag | - | Show detailed debug output in console |

**Example:**
```bash
# Basic usage
exp_heatmap plot chr15_results/ --start 47924019 --end 48924019 --title "SLC24A5" --out slc24a5

# Interactive HTML output
exp_heatmap plot chr15_results/ --start 47924019 --end 48924019 --interactive --out slc24a5_interactive
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
    rank_scores="2-tailed"    # Options: "2-tailed", "ascending", "descending"
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

## Input File Formats

### VCF File

Standard VCF format. For best results:
- Filter to SNPs only (remove indels)
- Use recoded VCF from vcftools

```bash
vcftools --gzvcf input.vcf.gz --remove-indels --recode --recode-INFO-all --out snps_only
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

<img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/SLC24A5_gene.png" width="800" alt="ExP heatmap of SLC24A5 gene">

This example demonstrates a full workflow analyzing the [SLC24A5](https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000188467;r=15:48120990-48142672) gene, known for its role in human skin pigmentation using 1000 Genomes Project data. SLC24A5 is also known to show strong selection signals, which makes it a suitable example.

```bash
#!/bin/bash

# Download 1000 Genomes data
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" -O chr15.vcf.gz
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" -O genotypes.panel

# Filter to SNPs only
vcftools --gzvcf chr15.vcf.gz \
    --remove-indels \
    --recode \
    --recode-INFO-all \
    --out chr15_snps

# Prepare data
exp_heatmap prepare chr15_snps.recode.vcf -o chr15_snps.recode.zarr

# Compute statistics
exp_heatmap compute chr15_snps.recode.zarr genotypes.panel -o chr15_snps_output

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

**Two-tailed rank scores:**
<img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/ADM2_2tailed.png" width="800" alt="Two-tailed rank scores">

**Ascending rank scores:**
<img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/ADM2_ascending.png" width="800" alt="Ascending rank scores">

**Descending rank scores:**
<img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/ADM2_descending.png" width="800" alt="Descending rank scores">

### Noise Filtering

Using `display_limit` and `display_values` parameters to filter noisy data and highlight significant regions:

<img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/ADM2_display_limit.png" width="800" alt="Filtered display">

*Same data as above, but with display_limit=1.60 to filter noise and highlight significant signals.*

## Statistical Methodology

### Empirical Rank Scores (Not P-Values)

ExP Heatmap uses **empirical rank scores** to visualize selection signals. Despite column names containing "p_value" (kept for backward compatibility), these are **not classical statistical p-values**.

**How rank scores are computed:**

1. **Rank all variants**: For each population pair, sort all genome-wide test statistics
2. **Calculate percentile**: `rank_score = rank / total_variants`
3. **Transform to -log10 scale**: `display_value = -log10(rank_score)`

**Interpretation:**
- Higher values indicate more extreme (potentially selected) variants
- A value of 3.0 means the variant is in the top 0.1% genome-wide
- A value of 2.0 means the variant is in the top 1% genome-wide

### Rank Score Options

The `rank_scores` parameter in `create_plot_input()` controls how rank scores are calculated:

| Option | Description | Use Case |
|--------|-------------|----------|
| `"2-tailed"` | For POP1_POP2: use descending; for POP2_POP1: use ascending | **Default** - Captures selection in either direction |
| `"ascending"` | Lowest test values ranked first | Detect negative selection signals |
| `"descending"` | Highest test values ranked first | Detect positive selection signals |

**Example:**
```python
# Two-tailed (recommended for most analyses)
data = create_plot_input("results/", start=47000000, end=49000000, rank_scores="2-tailed")

# One-tailed for specific hypothesis
data = create_plot_input("results/", start=47000000, end=49000000, rank_scores="descending")
```

### Statistical Tests Explained

| Test | Detects | Interpretation |
|------|---------|----------------|
| **XP-EHH** | Recent positive selection | Positive values: selection in pop1; Negative: selection in pop2 |
| **XP-nSL** | Selection (robust to recombination rate variation) | Similar to XP-EHH but more robust |
| **Delta Tajima's D** | Difference in allele frequency spectrum | Positive: excess rare variants in pop1 |
| **Hudson's Fst** | Population differentiation | Higher values: greater genetic distance |

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

**Error:**
```
All positions have NaN results. No output will be generated.
```

**Cause:** The statistical test produced no valid results, often due to:
- Too few variants
- Insufficient allele frequency variation
- For `delta_tajima_d`: window size too large

**Solution:**
- Check your data has sufficient variants
- For `delta_tajima_d`, the default window size is 13 SNPs; ensure your data has enough variants

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
pip install -e .
```

## License

This project is licensed under a Custom Non-Commercial License based on the MIT License - see the [LICENSE](https://github.com/bioinfocz/exp_heatmap/blob/main/LICENSE) file for details.

**Note**: This is NOT the standard MIT License. Commercial use requires explicit written permission from the authors.

For commercial licensing under different terms, please contact: edvard.ehler@img.cas.cz

## Contributors

- **Edvard Ehler** ([@EdaEhler](https://github.com/EdaEhler)) - Lead Developer
- **Adam Nógell** ([@AdamNogell](https://github.com/AdamNogell)) - Developer
- **Jan Pačes** ([@hpaces](https://github.com/hpaces)) - Developer
- **Mariana Šatrová** ([@satrovam](https://github.com/satrovam)) - Developer  
- **Ondřej Moravčík** ([@ondra-m](https://github.com/ondra-m)) - Developer

## Acknowledgments

<div align="center">

<a href="http://genomat.img.cas.cz">
  <img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/genomat.png" width="100" alt="GenoMat">
</a>
&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://www.img.cas.cz/en">
  <img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/img.png" width="100" alt="IMG CAS">
</a>
&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://www.elixir-czech.cz">
  <img src="https://github.com/bioinfocz/exp_heatmap/raw/main/assets/elixir.png" width="100" alt="ELIXIR">
</a>

</div>

---

*If you use ExP Heatmap in your research, please cite our paper [citation details will be added upon publication].*
