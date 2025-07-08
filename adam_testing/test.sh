#! /bin/bash

################
# GET THE DATA #
################
# Download chromosome 22 from 1000genomes ftp (GRch38 version)
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" -O chr22.genotypes.vcf.gz
wget "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" -O genotypes.panel

####################
# PREPARE THE DATA #
####################
# Filter the VCF
vcftools --gzvcf chr22.genotypes.vcf.gz --remove-indels --recode --recode-INFO-all --out chr22.genotypes

exp_heatmap prepare chr22.genotypes.recode.vcf chr22.genotypes.recode.zarr

#################################################
# COMPUTE PAIRWISE VALUES (default XP-EHH test) #
#################################################
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