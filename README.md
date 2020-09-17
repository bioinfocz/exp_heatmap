# ExP Selection

## Requirements

- python3
- vcftools ([repository](https://github.com/vcftools/vcftools))
- space on disk (.vcf files are usually quite large)

## Install

Get the repository

```bash
git clone git@github.com:ondra-m/exp-selection.git
```

Install python packages

```bash
pip3 install -r requirements.txt
```

## Usage

**Get data**

- VCF files (e.g. [1000 Genomes Project](https://www.internationalgenome.org/data) and [Phase 3, chr22](ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz))
- Panel file (e.g. [1000 Genomes Project](ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel))

**Extract only SNP**

You can give an .vcf or .vcf.gz file

```bash
bin/recode FILE.vcf
# => FILE.recode.vcf
```

**Prepare data for computing**

```bash
bin/prepare FILE.recode.vcf
# => FILE.recode.zarr
```

**Compute**

```bash
bin/compute FILE.recode.zarr PANEL_FILE
# => FILE.recode.xpehh
```

**Plot graph**

- `begin`, `end` (required)
   - plot boundaries
- `title` (optional)
   - name of the image
- `cmap` (optional)
   - color schema
   - [more informations at seaborn package](http://seaborn.pydata.org/tutorial/color_palettes.html)

```bash
bin/plot FILE.recode.xpehh --begin BEING --end END --title TITLE
# => TITLE.png
```

## Example

```bash
# Download datasets
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz" -O chr22.genotypes.vcf.gz
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" -O genotypes.panel


# Compute files for graph
bin/recode chr22.genotypes.vcf.gz
bin/prepare chr22.genotypes.recode.vcf
bin/compute chr22.genotypes.recode.zarr genotypes.panel

# Plot heatmap
bin/plot chr22.genotypes.recode.xpehh --begin 50481556 --end 50486440 --title ADM2

# A heatmap is saved at ADM2.png
```

# Contributors

- Eda Ehler ([@EdaEhler](https://github.com/EdaEhler))
- Jan Pačes ([@hpaces](https://github.com/hpaces))
- Mariana Šatrová ([@satrovam](https://github.com/satrovam))
- Ondřej Moravčík ([@ondra-m](https://github.com/ondra-m))
