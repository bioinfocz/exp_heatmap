# ExP Heatmap

Welcome to the ExP Heatmap `python` package and command-line tool. Our software is focused on displaying multidimensional data, expecially the so-called cross-population data - differences/similarities/p-values/or any other parameters of your choice between several groups/populations. Our method allows the user to quickly and efficiently evaluate millions of p-values or test statistics in one figure.

This tool is being developed in the [Laboratory of Genomics and Bioinformatics](https://www.img.cas.cz/group/michal-kolar/), Institute of Molecular Genetics of the Academy of Sciences of the Czech Republic, v. v. i.


The ExP Heatmap manual is divided into following sections:
1. [**Requirements and install**](#1-requirements-and-install)

2. [**Simple example**](#2-simple-example)

3. [**Workflow**](#3-workflow)

4. [**Command-line tool**](#4-command-line-tool)

5. [**Python package**](#5-python-package)

6. [**Galery**](#6-galery)

7. [**Licence and final remarks**](#7-licence-and-final-remarks)

<br/>

#### ExP heatmap example - LCT gene

<img src="https://github.com/bioinfocz/exp_heatmap/raw/master/assets/LCT_gene.png" width=800>

This is the ExP heatmap of human lactose (LCT) gene on chromosome 2 and its surrounding genomic region displaying population differences between 26 populations of 1000 Genomes Project, phase 3. Displayed values are the adjusted rank p-values for cross-population extended haplotype homozygosity (XPEHH) selection test.


## 1. Requirements and install

### Requirements

- `python` >= 3.8
- `vcftools` for genomic data preparation (not needed if you just want just plot your data) ([repository](https://github.com/vcftools/vcftools))
- space on disk (.vcf files are usually quite large)

### Install

PyPI repository link ([exp_heatmap](https://pypi.org/project/exp_heatmap/))

```bash
pip install exp_heatmap
```


Install the latest version directly from this GitHub
```bash
pip install git+https://github.com/bioinfocz/exp_heatmap.git
```
<br/>

## 2. Simple example

After installing the package, try to construct ExP heatmap in **three simple steps:**
1. **Download** the prepared results of the extended haplotype homozygosity (XPEHH) selection test for the part of human chromosome 2, 1000 Genomes Project data either directly via [Zenodo](https://zenodo.org/records/16364351) or via command:
```bash
wget "https://zenodo.org/records/16364351/files/chr2_output.tar.gz"
```
2. **Decompress** the downloaded folder in your working directory: `tar -xzf chr2_output.tar.gz`
3. **Run** the following command:
```bash
exp_heatmap plot chr2.xpehh.example/ --begin 136070087 --end 137070087 --title "LCT gene" --output LCT_xpehh
```
The `exp_heatmap` package will read the files from `chr2.xpehh.example/` folder and create the ExP heatmap and save it as `LCT_xpehh.png` file.

<br/>


## 3. Workflow
<img src="https://github.com/bioinfocz/exp_heatmap/raw/master/assets/SLC24A5_gene.png" width=800>

As a workflow example, we present an analysis of 1000 Genomes Project, phase 3 data of chromosome 15. It is focused on the [SLC24A5](https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000188467;r=15:48120990-48142672) gene, known for its role in human skin pigmentation. It is also known to show strong selection signals, which makes it a suitable example.

The bash script below presents a complete functional workflow analysis, including high-level comments that describe the entire process.

```bash
#!/bin/bash

#=============================================================================
# exp_heatmap Example Workflow
#=============================================================================
# This script demonstrates the complete workflow for using exp_heatmap, a tool
# for generating population genetics heatmaps from VCF (Variant Call Format) data.
# 
# The workflow covers:
# 1. Data acquisition from the 1000 Genomes Project
# 2. Data preprocessing and quality filtering
# 3. Data preparation and format conversion
# 4. Computation of population genetic statistics
# 5. Visualization of genetic variation patterns
#
# Final output: A heatmap visualization showing genetic variation patterns
# across different populations for the SLC24A5 gene region on chromosome 15.
#=============================================================================


#-----------------------------------------------------------------------------
# STEP 1: Download genomic data from 1000 Genomes Project
#-----------------------------------------------------------------------------
# Download the VCF file for chromosome 15 from the 1000 Genomes Project
# (Phase 3 data, GRCh37 reference genome build)
echo "Downloading chromosome 15 VCF file from 1000 Genomes Project..."
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" -O chr15.vcf.gz

# Download the panel file containing sample-to-population mappings
echo "Downloading population panel file..."
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" -O genotypes.panel

#-----------------------------------------------------------------------------
# STEP 2: Filter and preprocess genetic variants
#-----------------------------------------------------------------------------
# Use VCFtools to filter the raw VCF data:
# - Remove insertions/deletions (indels) to focus only on SNPs
# - Keep all INFO field annotations for downstream analysis
# - Output: chr15_snps.recode.vcf (SNPs-only VCF file)
echo "Filtering VCF to retain only SNPs (removing indels)..."
vcftools --gzvcf chr15.vcf.gz --remove-indels --recode --recode-INFO-all --out chr15_snps

#-----------------------------------------------------------------------------
# STEP 3: Prepare data for exp_heatmap analysis
#-----------------------------------------------------------------------------
# Convert the VCF file to Zarr format for efficient data access
echo "Converting VCF to Zarr format for efficient processing..."
exp_heatmap prepare chr15_snps.recode.vcf chr15_snps.recode.zarr

#-----------------------------------------------------------------------------
# STEP 4: Compute population genetic statistics
#-----------------------------------------------------------------------------
# Calculate population genetics metrics across the chromosome
# 
# Inputs: Zarr-formatted genotype data + population panel
# Output: Processed data ready for plotting
echo "Computing population genetic statistics..."
exp_heatmap compute chr15_snps.recode.zarr genotypes.panel chr15_snps_output

#-----------------------------------------------------------------------------
# STEP 5: Generate heatmap visualization
#-----------------------------------------------------------------------------
# Create a heatmap for the SLC24A5 gene region
#
# Parameters:
# --begin 47924019 --end 48024019: SLC24A5 gene region (±500kb)
# --title "SLC24A5": Label for the output plot
# --out SLC24A5_heatmap: Output filename prefix
echo "Generating heatmap for SLC24A5 gene region..."
exp_heatmap plot chr15_snps_output --begin 47924019 --end 48924019 --title "SLC24A5" --out SLC24A5_heatmap

echo "Workflow completed!"
echo ""
echo "Expected outputs:"
echo "- SLC24A5_heatmap.png: Main heatmap visualization"
echo "- chr15_snps_output/: Directory with computed statistics"
```
<br/>

## 4. Command-line tool

After installing the `exp_heatmap` using `pip` as described above, you can use its basic functionality directly from the command line interface.

### Get the data

- VCF files (e.g. [1000 Genomes Project](https://www.internationalgenome.org/data) and [Phase 3, chr22](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz))
- Panel file (e.g. [1000 Genomes Project](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel))


### Prepare the data

&emsp;  **Extract only SNP**

You can use a .vcf or .vcf.gz file

```bash
# we will use SNPs only, so we remove insertion/deletion polymorphisms
# another option would be to use only biallelic SNPs (--min-alleles 2 --max-alleles 2),
# probably with minor allele frequency above 5% (--maf 0.05)
# ouput VCF will be named DATA.recode.vcf
DATA="ALL.chr22_GRCh38.genotypes.20170504"


# Gziped VCF
vcftools --gzvcf $DATA.vcf.gz --remove-indels --recode --recode-INFO-all --out $DATA

# Plain VCF
vcftools --vcf $DATA.vcf --remove-indels --recode --recode-INFO-all --out $DATA
```

&emsp;  **Prepare data for computing**

```bash
# DATA.recode.vcf a vcf from previous step
# DATA.zarr is path (folder) where zarr representation of the VCF input will be saved
# this will vastly increase the speed of follow-up computations
exp_heatmap prepare DATA.recode.vcf DATA.zarr
```


### Compute pairwise values

```bash
# DATA.zarr a zarr data from previous step
# DATA.output a path (folder) where the results will be saved
# in this step, by default Cross-population extended haplotype homozygosity (XPEHH) score will be computed for all positions, together with their -log10 rank p-values.
exp_heatmap compute DATA.zarr genotypes.panel DATA.output
```
Besides the default cross-population extended haplotype homozygosity (XPEHH) test, you can use this `exp_heatmap compute` with optional parameter `-t` and one of the keywords:
- `xpehh` - computes cross-population extended haplotype homozygosity (XPEHH) test (default),
- `xpnsl` - computes cross-population number of segregating sites by length (NSL) test,
- `delta_tajima_d` - computes delta Tajima's D,
- `hudson_fst` - computes pairwise genetic distance Fst (using the method od Hudson (1992)).

```bash
# computing the XP-NSL test
exp_heatmap compute DATA.zarr genotypes.panel DATA.output -t xpnsl
```


### Display ExP heatmap

- `--begin`, `--end` (required)
  - plot boundaries
- `--title` (optional)
  - name of the image
- `--cmap` (optional)
  - color schema
  - [more informations at seaborn package](http://seaborn.pydata.org/tutorial/color_palettes.html)
- `--output` (optional)
  - png output path

```bash
exp_heatmap plot DATA.output --begin BEING --end END --title TITLE --output NAME
```
<br/>

## 5. Python package

Besides using ExP Heatmap as a standalone command-line tool, more options and user-defined parameters' changes are available when ExP Heatmap is imported directly into your Python script.

Test files used in these examples (p-values, test results, VCF files etc.) can be downloaded [HERE](http://genomat.img.cas.cz/). They are based on results of cross-population selection tests of the lactase ([LCT](https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000115850;r=2:135787850-135837184)) gene area (
chr2:135,787,850-135,837,184).

Here we outline a solution to 3 possible and most common scenarios where the ExP is being used.
Possible model scenarios:
* **a) you have values ready to display**
* **b) you have some kind of parameters/test results, need to compute the p-values and display them**
* **c) you only have the input data (VCF), need to compute the parameters/tests, turn them into p-values and display them as ExP heatmap**


### a) you have values ready to display 
Your data are in a \*.tsv file, tab-delimited text file (table), where the results or p-values are stored in columns, first column is 'CHROM', second column 'POS', followed by X columns of pairwise parameters (i.e. rank p-values). For 1000 Genomes data, that would mean 325 columns of pair-wise p-values for 26 populations.

```python
from exp_heatmap.plot import plot_exp_heatmap
import pandas as pd

# input data in the form of pandas DataFrame, expected shape (x, 327)
# 327 columns consisting of CHROM, POS and 325 columns of pairwise p-values
# x represents the number of SNPs to display
# column names are expected to include the 1000 Genomes population abbreviations
data_to_plot = pd.read_csv("LCT_xpnsl_pvalues.csv", sep="\t")


plot_exp_heatmap(data_to_plot,
                 begin=135287850,
                 end=136287850,
                 title="XP-NSL test on LCT gene in 1000 Genomes Project (phase 3) populations",
                 cmap="Blues",
                 output=False,  # enter the save file name here
                 populations="1000Genomes",
                 xlabel="LCT gene, 1 Mbp window, 2:135,287,850-136,287,850, GRCh38")
```

<br/>

### b) you have some kind of parameters/test results, need to compute the p-values and display them

Here, you will need to compute the p-values using a prepared function in `exp_heatmap` python package.

```python
from exp_heatmap.plot import plot_exp_heatmap, create_plot_input, superpopulations, prepare_cbar_params

# input data are in the form of pairwise population results per file
# here, the results of XP-NSL test for populations of 1000 Genome Project dataset
results_directory = "chr2_xpnsl_1000Genomes.test/"

# compute ranked p-values and prepare data for ExP heatmap
data_to_plot = create_plot_input("chr2_xpnsl_1000Genomes.test/", begin=135287850, end=136287850, populations="1000Genomes")


plot_exp_heatmap(data_to_plot,
                 begin=135287850,
                 end=136287850,
                 title="XP-NSL test on LCT gene in 1000 Genomes Project (phase 3) populations",
                 cmap="Blues",
                 output=False,  # enter the save file name here
                 populations="1000Genomes",
                 xlabel="LCT gene, 1 Mbp window, 2:135,287,850-136,287,850, GRCh38")


#######################################################################################
# you can tweak different paramaters in the ExP heatmap plot
# prepare custom colorbar parameters

cmin, cmax, cbar_ticks = prepare_cbar_params(data_to_plot, n_cbar_ticks=4)

# display custom population set
plot_exp_heatmap(data_to_plot,
                 begin=135000000,
                 end=137000000,
                 title="XP-NSL test results in African populations",
                 cmap="expheatmap",  # custom heatmap
                 output="xpnsl_Africa",  # save results
                 vertical_line=([135851073, "rs41525747"], [135851081, "rs41380347"], [135851176, "rs145946881"]), # 3 vertical lines marking SNPs with described selection pressure (https://doi.org/10.1093/gbe/evab065)
                 populations=superpopulations["AFR"], # custom population set
                 xlabel="LCT gene region, 2:135,000,000-137,000,000, GRCh38")

```

<br/>

### c) you only have the input data (vcf)...

...and need to compute the parameters/tests, turn them into p-values and display them as ExP heatmap.
Here the process will differ depending on what kind test you want to run. Below we give different examples using common tools (`VCFtools`)
and pythonic library `scikit-allel`.

```python
XX VCF to zarr
XX Compute test (different!)
XX Display!
XX
XX
```

<br/>

### d) custom population set
xxxx


## 6. Galery


## 7. Licence and final remarks

The ExP Heatmap package is available under the MIT License. ([link](https://github.com/bioinfocz/exp_heatmap?tab=MIT-1-ov-file "ExP Heatmap MIT licence"))

If you are interested in using this method in your commercial software under another licence, please, contact us at edvard.ehler@img.cas.cz.



<br/>

# Contributors

- Eda Ehler ([@EdaEhler](https://github.com/EdaEhler))
- Jan Pačes ([@hpaces](https://github.com/hpaces))
- Mariana Šatrová ([@satrovam](https://github.com/satrovam))
- Ondřej Moravčík ([@ondra-m](https://github.com/ondra-m))

# Acknowledgement

<a href="http://genomat.img.cas.cz">
  <img src="https://github.com/bioinfocz/exp_heatmap/raw/master/assets/genomat.png" width=100>
</a>

---

<a href="https://www.img.cas.cz/en">
  <img src="https://github.com/bioinfocz/exp_heatmap/raw/master/assets/img.png" width=100>
</a>

---

<a href="https://www.elixir-czech.cz">
  <img src="https://github.com/bioinfocz/exp_heatmap/raw/master/assets/elixir.png" width=100>
</a>
