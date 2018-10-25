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

###Get SNPs from vcf file

We will use gatk SelectVariants to extract SNPs from VQSR filtered file create in Germline SNV and Indel Calling section

```bash
#Need to put where to go in file system

#Extract SNPs
gatk --java-options '-Xmx64g' SelectVariants -R ~/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V ~/data/germline_variants/WGS_Norm_HC_calls_recalibrated.PASS.vcf -select-type SNP -O WGS_Norm_HC_calls.recalibrated.PASS.snps.vcf

```

Next, we will we use gatk VariantsToTable to extract columns with data relevant for calculating germline VAFs

```bash

#Create table that has columns for chromosome, position, genotype, allele depth, total depth
gatk --java-options '-Xmx64g' VariantsToTable -V WGS_Norm_HC_calls.recalibrated.PASS.snps.vcf -F CHROM -F POS -GF GT -GF AD -GF DP -O WGS_Norm_HC_calls.recalibrated.PASS.snps.table

```

To calculate VAFs, start R, follow scripts 
