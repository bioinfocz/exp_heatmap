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

### Get data

- VCF files ([1000 Genomes Project](https://www.internationalgenome.org/data))
- Panel file ([1000 Genomes Project](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel))

### Extract only SNP

```bash
# Automatic way
bin/recode VCF_FILE

# Or manual way
vcftools --gzvcf VCF_FILE --remove-indels --recode --recode-INFO-all --out OUTPUT
```

### Prepare data for computing

```bash
bin/prepare RECODE_FILE
```

### Compute

```bash
bin/compute ZARR_DIR PANEL_FILE
```

# Contributors

- Eda Ehler
- Jan Pačes
- Mariana Šatrová
- Ondřej Moravčík ([@ondra-m](https://github.com/ondra-m))
