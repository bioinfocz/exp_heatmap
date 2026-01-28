# ExP Heatmap

[![PyPI version](https://badge.fury.io/py/exp-heatmap.svg)](https://pypi.org/project/exp_heatmap/)

> A powerful Python package and command-line tool for visualizing multidimensional population genetics data through intuitive heatmaps.

ExP Heatmap specializes in displaying cross-population data, including differences, similarities, p-values, and other statistical parameters between multiple groups or populations. This tool enables efficient evaluation of millions of statistical values in a single, comprehensive visualization.

<img src="https://github.com/bioinfocz/exp_heatmap/raw/master/assets/LCT_gene.png" width="800" alt="ExP heatmap of LCT gene">

*ExP heatmap of the human lactose (LCT) gene showing population differences between 26 populations from the 1000 Genomes Project, displaying empirical rank scores for cross-population extended haplotype homozygosity (XPEHH) selection test. Create your own LCT heatmap with the [Quick Start](#quick-start) Guide*

**Developed by the [Laboratory of Genomics and Bioinformatics](https://www.img.cas.cz/group/michal-kolar/), Institute of Molecular Genetics of the Academy of Sciences of the Czech Republic**

## Features

- **Multiple Statistical Tests**: Support for XPEHH, XP-NSL, Delta Tajima's D, and Hudson's Fst
- **Flexible Input Formats**: Work with VCF files, pre-computed statistics, or ready-to-plot data
- **Command-Line Interface**: Easy-to-use CLI for standard workflows
- **Python API**: Full programmatic control for custom analyses
- **Efficient Processing**: Zarr-based data storage for fast computation
- **Customizable Visualization**: Multiple color schemes, resolution options, and display settings
TODO: ADD THE POINTS BELOW ONLY IF IMPLEMENTED INTO PROD
- **Interactive Mode**: Plotly-based HTML visualizations with zoom, pan, and hover tooltips
- **Row Clustering**: Hierarchical clustering to reveal patterns across population pairs
- **Superpopulation Annotations**: Color-coded bars indicating continental ancestry groups
- **Performance Benchmarking**: Built-in tools for measuring runtime and memory usage

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Command-Line Interface](#command-line-interface)
  - [Python Package](#python-package)
- [Advanced Features](#advanced-features)
- [Workflow Examples](#workflow-examples)
- [Gallery](#gallery)
- [Statistical Methodology](#statistical-methodology)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Requirements

- Python ≥ 3.8
- `vcftools` (for genomic data preparation - optional if using preprocessed data)

### Install from PyPI

```bash
pip install exp_heatmap
```

### Install from GitHub (latest version)

```bash
pip install git+https://github.com/bioinfocz/exp_heatmap.git
```
TODO: ADD THE SECTION BELOW ONLY IF IMPLEMENTED INTO PROD
### Optional Dependencies

For interactive visualizations:
```bash
pip install plotly
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
exp_heatmap plot chr2_output/ --start 136108646 --end 137108646 --title "LCT gene" --output LCT_xpehh
```
The `exp_heatmap` package will read the files from `chr2_output/` folder and create the ExP heatmap and save it as `LCT_xpehh.png` file.

<br/>

## Usage

ExP Heatmap follows a simple three-step workflow: **prepare** → **compute** → **plot**. Each step can be used independently depending on your data format.

### Command-Line Interface

#### 1.  Data Preparation - `prepare`

>Convert VCF files to efficient Zarr format for faster computation.

```bash
exp_heatmap prepare [OPTIONS] <vcf_file>
```

- `<vcf_file> [PATH]`: Recoded VCF file
- `-o, --output [PATH]`: Directory for output files

#### 2. Statistical Analysis - `compute`

>Calculate population genetic statistics across all genomic positions.

```bash
exp_heatmap compute [OPTIONS] <zarr_dir> <panel_file>
```

`<zarr_dir> [PATH]`: Directory with ZARR files from `prepare` step
`<panel_file>[PATH]`: Population panel file
- `-o, --output`: Directory for output files
- `-t, --test`: Statistical test to compute
  - `xpehh`: Cross-population Extended Haplotype Homozygosity (default)
  - `xpnsl`: Cross-population Number of Segregating sites by Length  
  - `delta_tajima_d`: Delta Tajima's D
  - `hudson_fst`: Hudson's Fst genetic distance
- `-c, --chunked`: Use chunked array to avoid memory exhaustion

#### 3. Visualization - `plot`

>Generate heatmap visualizations from computed statistics.

```bash
exp_heatmap plot [OPTIONS] <input_dir>
```

- `<input_dir>`: Directory with TSV files from `compute` step
- `-s, --start & -e, --end`: Genomic coordinates for the region to display
- `-m, --mid`: Alternative way to specify region (mid ± 500 kb)
- `-t, --title`: Title of the heatmap
- `-o, --output`: Output filename (without extension)
- `-c, --cmap`: Matplotlib colormap - [list of colormaps](https://matplotlib.org/stable/users/explain/colors/colormaps.html)
TODO: ADD THE POINTS BELOW ONLY IF IMPLEMENTED INTO PROD
- `--dpi`: Resolution of output image (default: 400)
- `--figsize`: Figure size as "WIDTH,HEIGHT" in inches
- `--cluster-rows`: Cluster rows by similarity for pattern discovery
- `--no-superpop-colors`: Disable superpopulation color annotation bar
- `--interactive`: Generate interactive HTML visualization (requires plotly)

TODO: ADD THE SECTION BELOW ONLY IF IMPLEMENTED INTO PROD
#### 4. Performance Benchmarking - `benchmark`

>Measure runtime and memory usage across the pipeline.

```bash
exp_heatmap benchmark <vcf_file> <panel_file> -s <start> -e <end> [OPTIONS]
```

- `-o, --output`: Output prefix for benchmark files
- `-t, --test`: Statistical test to benchmark
- `--report`: Save benchmark report to file

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
    rank_pvalues="2-tailed"    # Options: "2-tailed", "ascending", "descending"
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
subprocess.run(["exp_heatmap", "prepare", "data_snps.recode.vcf", "data.zarr"])

# 2. Compute statistics (using CLI) 
subprocess.run(["exp_heatmap", "compute", "data.zarr", "populations.panel", "results/"])

# 3. Create custom plots (using Python)
data_to_plot = create_plot_input("results/", start=47000000, end=49000000)
plot_exp_heatmap(data_to_plot, start=47000000, end=49000000, 
                 title="Custom Analysis", output="custom_plot")
```
TODOSTART: ADD THE SECTION BELOW ONLY IF IMPLEMENTED INTO PROD
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
exp_heatmap plot results/ --start 135287850 --end 136287850 --interactive --output lct_interactive
```

### Row Clustering

Cluster population pairs by similarity to reveal hidden patterns:

```python
from exp_heatmap.plot import plot_exp_heatmap, cluster_rows

# Cluster rows before plotting
clustered_data = cluster_rows(data_to_plot, method='average', metric='euclidean')
plot_exp_heatmap(clustered_data, ...)

# Or use the built-in parameter
plot_exp_heatmap(data_to_plot, cluster_rows_by='euclidean', ...)
```

### Superpopulation Summaries

Collapse data to superpopulation-level summaries:

```python
from exp_heatmap.plot import summarize_by_superpopulation

# Aggregate by superpopulation
summary_data = summarize_by_superpopulation(data_to_plot, agg_func='mean')
```

### Top Regions Extraction

Automatically identify top-scoring genomic regions:

```python
from exp_heatmap.plot import extract_top_regions

# Find top 50 selection signal regions
top_regions = extract_top_regions(data_to_plot, n_top=50, window_size=10000)
print(top_regions)
```

### Performance Benchmarking

Measure and report pipeline performance:

```python
from exp_heatmap.benchmark import run_full_benchmark, generate_benchmark_report

# Run complete benchmark
results = run_full_benchmark(
    vcf_file="data.vcf",
    panel_file="populations.panel",
    start=135000000,
    end=137000000
)

# Generate report
report = generate_benchmark_report(results, output_file="benchmark_report.txt")
```
TODOEND: ADD THE SECTION BELOW ONLY IF IMPLEMENTED INTO PROD

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
    TODO: ADD THE VALUES BELOW ONLY IF IMPLEMENTED INTO PROD
    dpi=600,                              # High resolution output
    figsize=(20, 8),                      # Custom figure size
    
    # Annotations
    vertical_line=[                       # Mark important SNPs
        [135851073, "rs41525747"],        # [position, label]
        [135851081, "rs41380347"]
    ],
    TODO: ADD THE VALUE BELOW ONLY IF IMPLEMENTED INTO PROD
    show_superpop_colors=True,            # Show superpopulation color bar
    
    # Colorbar customization
    cbar_vmin=cmin,
    cbar_vmax=cmax,
    cbar_ticks=cbar_ticks,
    
    # Output
    output="african_populations_analysis",
    xlabel="Custom region description"
)
```

## Workflow Examples

### Complete Analysis: SLC24A5 Gene

<img src="https://github.com/bioinfocz/exp_heatmap/raw/master/assets/SLC24A5_gene.png" width="800" alt="ExP heatmap of SLC24A5 gene">

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
exp_heatmap prepare chr15_snps.recode.vcf chr15_snps.recode.zarr

# Compute statistics
exp_heatmap compute chr15_snps.recode.zarr genotypes.panel chr15_snps_output

# Generate heatmap for SLC24A5 region
exp_heatmap plot chr15_snps_output \
    --start 47924019 \
    --end 48924019 \
    --title "SLC24A5" \
    --cmap gist_heat \
    --output SLC24A5_heatmap
```

## Gallery

### Different Rank Score Computations

The same XP-EHH test data for the ADM2 gene region, showing different rank score calculation methods:

**Two-tailed rank scores:**
<img src="https://github.com/bioinfocz/exp_heatmap/raw/master/assets/ADM2, chr22, XP-EHH, pvals: 2-tailed.png" width="800" alt="Two-tailed rank scores">

**Ascending rank scores:**
<img src="https://github.com/bioinfocz/exp_heatmap/raw/master/assets/ADM2, chr22, XP-EHH, pvals: ascending.png" width="800" alt="Ascending rank scores">

**Descending rank scores:**
<img src="https://github.com/bioinfocz/exp_heatmap/raw/master/assets/ADM2, chr22, XP-EHH, pvals: descending.png" width="800" alt="Descending rank scores">

### Noise Filtering

Using `display_limit` and `display_values` parameters to filter noisy data and highlight significant regions:

<img src="https://github.com/bioinfocz/exp_heatmap/raw/master/assets/ADM2_XP-EHH_display_limit.png" width="800" alt="Filtered display">

*Same data as above, but with display_limit=1.60 to filter noise and highlight significant signals.*

## Contributing

We welcome contributions! Feel free to contact us or submit issues or pull requests.

### Development Setup

```bash
git clone https://github.com/bioinfocz/exp_heatmap.git
cd exp_heatmap
pip install -e .
```

## License

This project is licensed under a Custom Non-Commercial License based on the MIT License - see the [LICENSE](https://github.com/bioinfocz/exp_heatmap/blob/master/LICENSE) file for details.

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
  <img src="https://github.com/bioinfocz/exp_heatmap/raw/master/assets/genomat.png" width="100" alt="GenoMat">
</a>
&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://www.img.cas.cz/en">
  <img src="https://github.com/bioinfocz/exp_heatmap/raw/master/assets/img.png" width="100" alt="IMG CAS">
</a>
&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://www.elixir-czech.cz">
  <img src="https://github.com/bioinfocz/exp_heatmap/raw/master/assets/elixir.png" width="100" alt="ELIXIR">
</a>

</div>

---

*If you use ExP Heatmap in your research, please cite our paper [citation details will be added upon publication].*
