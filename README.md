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
1. **Download** the prepared results of the extended haplotype homozygosity (XPEHH) selection test for the part of human chromosome 2, 1000 Genomes Project data: [example results](https://github.com/bioinfocz/exp_heatmap/raw/master/assets/chr2.xpehh.example.zip)
2. **Unpack** the zipped folder `chr2.xpehh.example/` in your working directory: `unzip chr2.xpehh.example.zip`
3. **Run** the following command:
```bash
exp_heatmap plot chr2.xpehh.example/ --begin 135287850 --end 136287850 --title "LCT gene" --output LCT_xpehh
```
The `exp_heatmap` package will read the files from `chr2.xpehh.example/` folder and create the ExP heatmap and save it as `LCT_xpehh.png` file.

<br/>


## 3. Workflow

<img src="https://github.com/bioinfocz/exp_heatmap/blob/master/assets/ExP_process_schema.png" width=1100>

As a workflow example we present an analysis of 1000 Genomes Project, phase 3 data of chromosome 22, chosen especially for its small size and thus reasonable fast computations. It is focused on ADM2 gene ([link](https://www.ensembl.org/Homo_sapiens/Gene/Phenotype?db=core;g=ENSG00000128165;r=22:50481543-50486440)), which is active especially in reproductive system, and angiogenesis and cardiovascular system in general.

```bash
################
# GET THE DATA #
################
# Download chromosome 22 from 1000genomes ftp
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz" -O chr22.genotypes.vcf.gz
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" -O genotypes.panel

# OR

# The 1000 Genomes Project alternative ftp mirror (GRCh37 version);
wget "https://ddbj.nig.ac.jp/public/mirror_database/1000genomes/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -O chr22.genotypes.vcf.gz
wget "https://ddbj.nig.ac.jp/public/mirror_database/1000genomes/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" -O genotypes.panel

####################
# PREPARE THE DATA #
####################
# Filter the VCF
vcftools --gzvcf chr22.genotypes.vcf.gz \
         --remove-indels \
         --recode \
         --recode-INFO-all \
         --out chr22.genotypes

exp_heatmap prepare chr22.genotypes.recode.vcf chr22.genotypes.recode.zarr

###########################
# COMPUTE PAIRWISE VALUES #
###########################
exp_heatmap compute chr22.genotypes.recode.zarr genotypes.panel chr22.genotypes.output

#######################
# DISPLAY ExP HEATMAP #
#######################
# Plot heatmap
exp_heatmap plot chr22.genotypes.output --begin 50481556 --end 50486440 --title ADM2 --output adm2_GRCh38

# OR

# use this if you used the GRCh37 version of the VCF input files.
exp_heatmap plot chr22.genotypes.output --begin 50910000 --end 50950000 --title ADM2 --output adm2_GRCh37

# The heatmap is saved as adm2_GRCh38.png or adm2_GRCh37.png, depending on which version of the plot function you are using.
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
