# ExP Selection

**LCT gene**

<img src="https://github.com/ondra-m/exp-selection/raw/master/assets/LCT_gene.png" width=400>

## Requirements

- python >= 3.8
- vcftools ([repository](https://github.com/vcftools/vcftools))
- space on disk (.vcf files are usually quite large)

## Install

```bash
pip install exp-selection
```

## Usage

**Get data**

- VCF files (e.g. [1000 Genomes Project](https://www.internationalgenome.org/data) and [Phase 3, chr22](ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz))
- Panel file (e.g. [1000 Genomes Project](ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel))

**Extract only SNP**

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

**Prepare data for computing**

```bash
# DATA.recode.vcf a vcf from previous step
# DATA.zarr is path (folder) where zarr representation of the VCF input will be saved
exp-selection prepare DATA.recode.vcf DATA.zarr
```

**Compute**

```bash
# DATA.zarr a zarr data from previous step
# DATA.xpehh a path (folder) where the results will be saved
# in this step, by default Cross-population extended haplotype homozygosity (XPEHH) score will be computed for all positions provided together with their -log10 rank p-values.
exp-selection compute DATA.zarr genotypes.panel DATA.xpehh
```

**Plot graph**

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
exp-selection plot DATA.xpehh --begin BEING --end END --title TITLE --output NAME
```

## Example

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
exp-selection compute chr22.genotypes.recode.zarr genotypes.panel chr22.genotypes.recode.xpehh

# Plot heatmap
exp-selection plot chr22.genotypes.recode.xpehh --begin 50481556 --end 50486440 --title ADM2 --output adm2_GRCh38
exp-selection plot chr22.genotypes.recode.xpehh --begin 50910000 --end 50950000 --title ADM2_test --output adm2_GRCh37 # use this plotting if you use GRCh37 version of the VCF input files.

# A heatmap is saved at adm2.png
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
