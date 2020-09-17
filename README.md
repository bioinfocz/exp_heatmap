# ExP Selection

## Requirements

- python3
- vcftools ([github repo](https://github.com/vcftools/vcftools))
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

- VCF files (e.g. [1000 Genomes Project](https://www.internationalgenome.org/data))
- Panel file (e.g. [1000 Genomes Project](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel))

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
```

# Contributors

- Eda Ehler ([@EdaEhler](https://github.com/EdaEhler))
- Jan Pačes ([@hpaces](https://github.com/hpaces))
- Mariana Šatrová ([@satrovam](https://github.com/satrovam))
- Ondřej Moravčík ([@ondra-m](https://github.com/ondra-m))
