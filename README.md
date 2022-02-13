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
# Gziped VCF
vcftools --gzvcf DATA.vcf.gz --remove-indels --recode --recode-INFO-all --out DATA.recode.vcf

# Plan VCF
vcftools --vcf DATA.vcf --remove-indels --recode --recode-INFO-all --out DATA.recode.vcf
```

**Prepare data for computing**

```bash
# DATA.recode.vcf a vcf from previous step
# DATA.zarr is path where zarr will be saved
exp-selection prepare DATA.recode.vcf DATA.zarr
```

**Compute**

```bash
# DATA.zarr a zarr data from previous step
# DATA.xpehh a path where xpehh will be saved
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
exp-selection DATA.recode.xpehh --begin BEING --end END --title TITLE --output NAME
```

## Example

```bash
# Download datasets
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz" -O chr22.genotypes.vcf.gz
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" -O genotypes.panel


# Compute files for graph
vcftools --gzvcf chr22.genotypes.vcf.gz \
         --remove-indels \
         --recode \
         --recode-INFO-all \
         --out chr22.genotypes.recode.vcf

exp-selection prepare chr22.genotypes.recode.vcf
exp-selection compute chr22.genotypes.recode.zarr genotypes.panel

# Plot heatmap
exp-selection plot chr22.genotypes.recode.xpehh --begin 50481556 --end 50486440 --title ADM2 --output adm2

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
