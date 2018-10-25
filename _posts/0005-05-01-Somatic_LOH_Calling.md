---
feature_text: |
  ## Precision Medicine
title: Somatic LOH Calling
categories:
    - Module-05-Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-05-01
---

### Objectives
-Calculate VAFs in normal sample, find heterozygyous germline SNPs

-Run bam-readcount on tumor sample, use readcounts to calculate tumor VAFs at heterozygous positions

-Determine regions of LOH in tumor sample

Note: We are only going to analyze LOH in chromosomes 6 and 17

Get SNPs from vcf file

We will use gatk SelectVariants to extract SNPs from VQSR filtered file create in Germline SNV and Indel Calling section

```bash
#Need to put where to go in file system

#Extract SNPs
gatk --java-options '-Xmx64g' SelectVariants -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated.PASS.vcf -select-type SNP -O WGS_Norm_HC_calls.recalibrated.PASS.snps.vcf

```

Next, we will we use gatk VariantsToTable to extract columns with data relevant for calculating germline VAFs

```bash

#Create table that has columns for chromosome, position, genotype, allele depth, total depth
gatk --java-options '-Xmx64g' VariantsToTable -V WGS_Norm_HC_calls.recalibrated.PASS.snps.vcf -F CHROM -F POS -GF GT -GF AD -GF DP -O WGS_Norm_HC_calls.recalibrated.PASS.snps.table

```

The code below will calcule germline VAFs, filtering down to heterozygous positions for chromosome 6 and chromosome 17

```bash
# Start R
R

# Set working directory, see what that should be
setwd("")
 
# Load libraries
library(dplyr)
library(tidyr)

# Read in table you made with VariantsToTable command
normal_calls <- read.delim("WGS_Norm_HC_calls.recalibrated.PASS.snps.table", header = TRUE, col.names = c("CHROM", "POS", "GT", "GT_AD", "DP"))

# Filter down to data for chromosomes 6 and 17, ";" separates commands so we can run multiple commands in one line
chr6_normal_calls <- subset.data.frame(normal_calls, CHROM == "chr6"); chr17_normal_calls <- subset.data.frame(normal_calls, CHROM == "chr17"); normal_calls <- merge(chr6_normal_calls, chr17_normal_calls, all = TRUE)

# Before starting calculations, filter for positions with at least 20x coverage
normal_calls <- subset.data.frame(normal_calls, DP >= 20)

# Split GT_AD column so each allele depth has it's own entry in rows, will make three new entries because some positions have multiple alternate alleles
normal_calls <- normal_calls %>% separate(GT_AD, c(AD_1, AD_2, AD_3), remove = FALSE)

# Make sure allele depth columns are numeric so they can be used to calculate VAFs
normal_calls$AD_1 <- as.numeric(normal_calls$AD_1); normal_calls$AD_2 <- as.numeric(normal_calls$AD_2); normal_calls$AD_3 <- as.numeric(normal_calls$AD_3) 

# Calculate VAFs, allele depth over total dpeth
normal_calls$VAF_1 <- normal_calls$AD_1/normal_calls$DP; normal_calls$VAF_2 <- normal_calls$AD_2/normal_calls$DP; normal_calls$VAF_3 <- normal_calls$AD_3/normal_calls$DP

# Split into separate files for each VAF, remerge so we have one column with all VAFs to filter for heterozygous variants
VAF_1 <- select(normal_calls, CHROM, POS, GT, AD_1, DP, VAF_1); VAF_2 <- select(normal_calls, CHROM, POS, GT, AD_2, DP, VAF_2); VAF_3 <- select(normal_calls, CHROM, POS, GT, AD_3, DP, VAF_3)
colnames(VAF_1) <- c("CHROM", "POS", "GT", "AD", "DP", "VAF"); colnames(VAF_2) <- c("CHROM", "POS", "GT", "AD", "DP", "VAF"); colnames(VAF_3) <- c("CHROM", "POS", "GT", "AD", "DP", "VAF") 
normal_calls <- normal_calls <- merge(VAF_1, VAF_2, all = TRUE); normal_calls <- merge(normal_calls, VAF_3, all = TRUE)

# Filter for heterozygous posistions, i.e. positions with VAF between .4 and .6
normal_calls <- subset.data.frame(normal_calls, VAF >= .4); normal_calls <- subset.data.frame(normal_calls, VAF >= .6)

# Create table containing info about heterozygous SNPs, will use later to plot results
write.table(normal_calls, file="heterozygous.snps.table", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Get list of heterozygous positions to use with bam-readcount
heterozygous.positions <- unique(normal_calls[c("CHROM", "POS")])
write.table(heterozygous.positions[ ,c("CHROM", "POS", "POS")], file = "heterozygous.positions.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
 
```
