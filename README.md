# ExP Selection

**ExP heatmap example - LCT gene**

<img src="https://github.com/ondra-m/exp-selection/raw/master/assets/LCT_gene.png" width=800>

This is the ExP heatmap of human lactose (LCT) gene on chromosome 2 and its surrounding genomic region displaying population differences between 26 populations of 1000 Genomes Project, phase 3. Displayed values are the adjusted rank p-values for cross-population extended haplotype homozygosity (XPEHH) selection test.

## Requirements

- python >= 3.8
- vcftools ([repository](https://github.com/vcftools/vcftools))
- space on disk (.vcf files are usually quite large)

## Install

```bash
pip install exp-selection
```

## Workflow

<img src="https://github.com/ondra-m/exp-selection/blob/master/assets/ExP_process_schema.png" width=1100>

## Usage

### Get the data

- VCF files (e.g. [1000 Genomes Project](https://www.internationalgenome.org/data) and [Phase 3, chr22](ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz))
- Panel file (e.g. [1000 Genomes Project](ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel))


### Prepare the data

&emsp;  **Extract only SNP**

You can give an .vcf or .vcf.gz file

```bash
# we will use SNPs only, so we remove insertion/deletion polymorphisms
# another option would be to use only biallelic SNPs (--min-alleles 2 --max-alleles 2),
# probably with minor allele frequency above 5% (--maf 0.05)
# ouput VCF will be named DATA.recode.vcf

# Gziped VCF
vcftools --gzvcf DATA.vcf.gz --remove-indels --recode --recode-INFO-all --out DATA

# Plain VCF
vcftools --vcf DATA.vcf --remove-indels --recode --recode-INFO-all --out DATA
```

&emsp;  **Prepare data for computing**

```bash
# DATA.recode.vcf a vcf from previous step
# DATA.zarr is path (folder) where zarr representation of the VCF input will be saved
exp-selection prepare DATA.recode.vcf DATA.zarr
```


### Compute pairwise values

```bash
# DATA.zarr a zarr data from previous step
# DATA.output a path (folder) where the results will be saved
# in this step, by default Cross-population extended haplotype homozygosity (XPEHH) score will be computed for all positions provided together with their -log10 rank p-values.
exp-selection compute DATA.zarr genotypes.panel DATA.output
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
exp-selection plot DATA.output --begin BEING --end END --title TITLE --output NAME
```


## Example

This example shows an analysis of 1000 Genomes Project, phase 3 data of chromosome 22, chosen especially for its small size and thus reasonable fast computations. It is focused on ADM2 gene ([link](https://www.ensembl.org/Homo_sapiens/Gene/Phenotype?db=core;g=ENSG00000128165;r=22:50481543-50486440)), which is active especially in reproductive system, and angiogenesis and cardiovascular system in general.

```bash
# Download datasets
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz" -O chr22.genotypes.vcf.gz
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" -O genotypes.panel

# The 1000 Genomes Project ftp seems not working, you can get the VCF files (GRCh37 version) at this mirror
wget "https://ddbj.nig.ac.jp/public/mirror_database/1000genomes/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" -O chr22.genotypes.vcf.gz
wget "https://ddbj.nig.ac.jp/public/mirror_database/1000genomes/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" -O genotypes.panel


# Compute files for graph
vcftools --gzvcf chr22.genotypes.vcf.gz \
         --remove-indels \
         --recode \
         --recode-INFO-all \
         --out chr22.genotypes

exp-selection prepare chr22.genotypes.recode.vcf chr22.genotypes.recode.zarr
exp-selection compute chr22.genotypes.recode.zarr genotypes.panel chr22.genotypes.output

# Plot heatmap
exp-selection plot chr22.genotypes.output --begin 50481556 --end 50486440 --title ADM2 --output adm2_GRCh38
exp-selection plot chr22.genotypes.output --begin 50910000 --end 50950000 --title ADM2 --output adm2_GRCh37 # use this plotting if you use GRCh37 version of the VCF input files.

# The heatmap is saved as adm2_GRCh38.png or adm2_GRCh37.png, depending on which version of plot function are you using.
```

# Contributors

- Eda Ehler ([@EdaEhler](https://github.com/EdaEhler))
- Jan Pačes ([@hpaces](https://github.com/hpaces))
- Mariana Šatrová ([@satrovam](https://github.com/satrovam))
- Ondřej Moravčík ([@ondra-m](https://github.com/ondra-m))

# Acknowledgement

<a href="http://genomat.img.cas.cz">
  <img src="https://github.com/ondra-m/exp-selection/raw/master/assets/genomat.png" width=100>
</a>

---

<a href="https://www.img.cas.cz/en">
  <img src="https://github.com/ondra-m/exp-selection/raw/master/assets/img.png" width=100>
</a>

---

<a href="https://www.elixir-czech.cz">
  <img src="https://github.com/ondra-m/exp-selection/raw/master/assets/elixir.png" width=100>
</a>
