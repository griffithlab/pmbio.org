---
feature_text: |
  ## Precision Medicine
title: Germline SNV/Indel Filtering/Annotation/Review
categories:
    - Module-04-Germline
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0004-02-02
---

### Module objectives

- Perform quality filtering of germline SNVs and indels 

The raw output of GATK HaplotypeCaller will include many variants with varying degrees of quality. For various reasons we might wish to further filter these to a higher confidence set of variants. The recommended approach is to use GATK VQSR. However, this requires a large number of variants. Sufficient numbers may be available for even a single sample of whole genome sequencing data. However for targeted sequencing (e.g., exome data) it is recommended to include a larger number (i.e., at least 30-50), preferably platform-matched (i.e., similar sequencing strategy), samples with variant calls. For our purposes, we will first demonstrate a less optimal hard filtering strategy. Then we will demonstrate VQSR filtering.  

It is strongly recommended to read the following documentation from GATK:
- [How to run VQSR](https://software.broadinstitute.org/gatk/documentation/article?id=2805)
- [How to apply hard filters](https://software.broadinstitute.org/gatk/documentation/article?id=2806) 

### Perform hard-filtering on Exome data

```
# Extract the SNPs from the call set
cd /workspace/germline/
gatk --java-options '-Xmx60g' SelectVariants -R /workspace/inputs/references/genome/ref_genome.fa -V /workspace/germline/Exome_Norm_HC_calls.vcf -select-type SNP -O /workspace/germline/Exome_Norm_HC_calls.snps.vcf

# Extract the indels from the call set
gatk --java-options '-Xmx60g' SelectVariants -R /workspace/inputs/references/genome/ref_genome.fa -V /workspace/germline/Exome_Norm_HC_calls.vcf -select-type INDEL -O /workspace/germline/Exome_Norm_HC_calls.indels.vcf

# Apply basic filters to the SNP call set
gatk --java-options '-Xmx64g' VariantFiltration -R /workspace/inputs/references/genome/ref_genome.fa -V /workspace/germline/Exome_Norm_HC_calls.snps.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name "basic_snp_filter" -O /workspace/germline/Exome_Norm_HC_calls.snps.filtered.vcf

# Apply basic filters to the INDEL call set
gatk --java-options '-Xmx64g' VariantFiltration -R /workspace/inputs/references/genome/ref_genome.fa -V /workspace/germline/Exome_Norm_HC_calls.indels.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name "basic_indel_filter" -O /workspace/germline/Exome_Norm_HC_calls.indels.filtered.vcf

# Merge filtered SNP and INDEL vcfs back together
gatk --java-options '-Xmx64g' MergeVcfs -I /workspace/germline/Exome_Norm_HC_calls.snps.filtered.vcf -I /workspace/germline/Exome_Norm_HC_calls.indels.filtered.vcf -O /workspace/germline/Exome_Norm_HC_calls.filtered.vcf

# Extract PASS variants only
gatk --java-options '-Xmx64g' SelectVariants -R /workspace/inputs/references/genome/ref_genome.fa -V /workspace/germline/Exome_Norm_HC_calls.filtered.vcf -O /workspace/germline/Exome_Norm_HC_calls.filtered.PASS.vcf --exclude-filtered

```









### Perform hard-filtering on WGS data

```
# Extract the SNPs from the call set
gatk --java-options '-Xmx64g' SelectVariants -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.vcf -select-type SNP -O /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.snps.vcf

# Extract the indels from the call set
gatk --java-options '-Xmx64g' SelectVariants -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.vcf -select-type INDEL -O /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.indels.vcf

# Apply basic filters to the SNP call set
gatk --java-options '-Xmx64g' VariantFiltration -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.snps.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "basic_snp_filter" -O /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.snps.filtered.vcf

# Apply basic filters to the INDEL call set
gatk --java-options '-Xmx64g' VariantFiltration -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.indels.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "basic_indel_filter" -O /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.indels.filtered.vcf

# Merge filtered SNP and INDEL vcfs back together
gatk --java-options '-Xmx64g' MergeVcfs -I /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.snps.filtered.vcf -I /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.indels.filtered.vcf -O /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.filtered.vcf

# Extract PASS variants only
gatk --java-options '-Xmx64g' SelectVariants -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.filtered.vcf -O /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.filtered.PASS.vcf --exclude-filtered

```

NOTE: At the time of writing (and as of gatk v4.0.10.1) it seems as if gatk VariantFiltration is broken for filtering on SOR or DP (site/INFO level). Filtering for these properties (e.g., SOR > 3.0 or DP < 20) works when applied as individual filters but not when included as multi-filter expressions as in the above example commands. I have reported this issue in the [GATK forums](https://gatkforums.broadinstitute.org/gatk/discussion/10928/variantfiltration-not-filtering-correctly).

NOTE: There are a number of MIXED type variants (multi-allelic with both SNP and INDEL alleles) that are currently dropped by the above workflow (selecting for only SNPs and INDELs). In the future we could consider converting these with gatk LeftAlignAndTrimVariants --split-multi-allelics.


### Perform VQSR-filtering of Exome variants that were joint called with 1KG exomes

When using VQSR filtering on exome data, there will be a smaller number of variants per sample compared to WGS. These are typically insufficient to build a robust recalibration model. If running on only a few samples, GATK recommends that you analyze samples jointly in cohorts of at least 30 samples. If necessary, add exomes from 1000G Project or comparable. Ideally, these should be processed with similar technical generation (technology, capture, read length, depth).

See also: https://drive.google.com/drive/folders/0BzI1CyccGsZiWlU5SXNvNnFGbVE (GATK Workshops -> 1709 -> GATKwr19-05-Variant_filtering.pdf)

```

#Build SNP recalibration model
#Note: Any annotations specified ("-an XX") below must actually be present in your VCF. They should have been added at an earlier step in the GATK workflow. If not, they could be added with VariantAnnotator
#Note: For exome data, exclude "-an DP" as this coverage metric should only be used if extreme deviations in coverage are not expected and indicative of errors
#Note: The "-an InbreedingCoeff" option is for a population level statistic that requires at least 10 samples in order to be computed (When?). For projects with fewer samples, or that includes many closely related samples (such as a family) please omit this annotation from the command line.

gatk --java-options '-Xmx64g' VariantRecalibrator -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/Exome_GGVCFs_jointcalls.vcf --resource hapmap,known=false,training=true,truth=true,prior=15.0:/home/ubuntu/data/reference/hapmap_3.3.hg38.vcf.gz --resource omni,known=false,training=true,truth=true,prior=12.0:/home/ubuntu/data/reference/1000G_omni2.5.hg38.vcf.gz --resource 1000G,known=false,training=true,truth=false,prior=10.0:/home/ubuntu/data/reference/1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource dbsnp,known=true,training=false,truth=false,prior=2.0:/home/ubuntu/data/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum --mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -O recalibrate_SNP.recal --tranches-file recalibrate_SNP.tranches --rscript-file recalibrate_SNP_plots.R 

#Apply recalibration to SNPs 
gatk --java-options '-Xmx64g' ApplyVQSR -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/Exome_GGVCFs_jointcalls.vcf --mode SNP --truth-sensitivity-filter-level 99.0 --recal-file /home/ubuntu/data/germline_variants/recalibrate_SNP.recal --tranches-file /home/ubuntu/data/germline_variants/recalibrate_SNP.tranches -O /home/ubuntu/data/germline_variants/Exome_GGVCFs_jointcalls_recalibrated_snps_raw_indels.vcf

#Build Indel recalibration model
gatk --java-options '-Xmx64g' VariantRecalibrator -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/Exome_GGVCFs_jointcalls_recalibrated_snps_raw_indels.vcf --resource mills,known=false,training=true,truth=true,prior=12.0:/home/ubuntu/data/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --resource dbsnp,known=true,training=false,truth=false,prior=2.0:/home/ubuntu/data/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum --mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --max-gaussians 4 -O recalibrate_INDEL.recal --tranches-file recalibrate_INDEL.tranches --rscript-file recalibrate_INDEL_plots.R  

#Apply recalibration to Indels
gatk --java-options '-Xmx64g' ApplyVQSR -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/Exome_GGVCFs_jointcalls_recalibrated_snps_raw_indels.vcf --mode INDEL --truth-sensitivity-filter-level 99.0 --recal-file /home/ubuntu/data/germline_variants/recalibrate_INDEL.recal --tranches-file /home/ubuntu/data/germline_variants/recalibrate_INDEL.tranches -O /home/ubuntu/data/germline_variants/Exome_GGVCFs_jointcalls_recalibrated.vcf

#Extract PASS variants only and only variants actually called and non-reference in our sample of interest
gatk --java-options '-Xmx64g' SelectVariants -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/Exome_GGVCFs_jointcalls_recalibrated.vcf -O /home/ubuntu/data/germline_variants/Exome_Norm_GGVCFs_jointcalls_recalibrated.PASS.vcf --exclude-filtered --exclude-non-variants --remove-unused-alternates --sample-name HCC1395BL_DNA

```


### Perform VQSR-filtering of WGS variants

```

#Build SNP recalibration model
#Note: Any annotations specified ("-an XX") below must actually be present in your VCF. They should have been added at an earlier step in the GATK workflow. If not, they could be added with VariantAnnotator
#Note: For exome data, exclude "-an DP" as this coverage metric should only be used if extreme deviations in coverage are not expected and indicative of errors
#Note: The "-an InbreedingCoeff" option is for a population level statistic that requires at least 10 samples in order to be computed (When?). For projects with fewer samples, or that includes many closely related samples (such as a family) please omit this annotation from the command line.

gatk --java-options '-Xmx64g' VariantRecalibrator -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.vcf --resource hapmap,known=false,training=true,truth=true,prior=15.0:/home/ubuntu/data/reference/hapmap_3.3.hg38.vcf.gz --resource omni,known=false,training=true,truth=true,prior=12.0:/home/ubuntu/data/reference/1000G_omni2.5.hg38.vcf.gz --resource 1000G,known=false,training=true,truth=false,prior=10.0:/home/ubuntu/data/reference/1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource dbsnp,known=true,training=false,truth=false,prior=2.0:/home/ubuntu/data/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an DP --mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -O recalibrate_SNP.WGS.recal --tranches-file recalibrate_SNP.WGS.tranches --rscript-file recalibrate_SNP_plots.R

#Apply recalibration to SNPs
gatk --java-options '-Xmx64g' ApplyVQSR -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.vcf --mode SNP --truth-sensitivity-filter-level 99.0 --recal-file /home/ubuntu/data/germline_variants/recalibrate_SNP.WGS.recal --tranches-file /home/ubuntu/data/germline_variants/recalibrate_SNP.WGS.tranches -O /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated_snps_raw_indels.vcf

#Build Indel recalibration model
gatk --java-options '-Xmx64g' VariantRecalibrator -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated_snps_raw_indels.vcf --resource mills,known=false,training=true,truth=true,prior=12.0:/home/ubuntu/data/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --resource dbsnp,known=true,training=false,truth=false,prior=2.0:/home/ubuntu/data/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an DP --mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --max-gaussians 4 -O recalibrate_INDEL.WGS.recal --tranches-file recalibrate_INDEL.WGS.tranches --rscript-file recalibrate_INDEL_plots.R

#Apply recalibration to Indels
gatk --java-options '-Xmx64g' ApplyVQSR -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated_snps_raw_indels.vcf --mode INDEL --truth-sensitivity-filter-level 99.0 --recal-file /home/ubuntu/data/germline_variants/recalibrate_INDEL.WGS.recal --tranches-file /home/ubuntu/data/germline_variants/recalibrate_INDEL.WGS.tranches -O /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated.vcf

#Extract PASS variants only and only variants actually called and non-reference in our sample of interest
gatk --java-options '-Xmx64g' SelectVariants -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated.vcf -O /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated.PASS.vcf --exclude-filtered --exclude-non-variants --remove-unused-alternates --sample-name HCC1395BL_DNA

```


### Perform VEP annotation of hard-filtered results

```
#VEP annotate hard-filtered exome results
#Output VEP VCF
~/bin/ensembl-vep/vep --cache --dir_cache /home/ubuntu/data/vep_cache --dir_plugins /home/ubuntu/data/vep_cache/Plugins --fasta /home/ubuntu/data/vep_cache/homo_sapiens/91_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --fork 8 --assembly=GRCh38 --offline --vcf --plugin Downstream --everything --terms SO --pick --coding_only --transcript_version -i /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.filtered.PASS.vcf -o /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.filtered.PASS.vep.vcf --force_overwrite 

#Output tabular VEP
~/bin/ensembl-vep/vep --cache --dir_cache /home/ubuntu/data/vep_cache --dir_plugins /home/ubuntu/data/vep_cache/Plugins --fasta /home/ubuntu/data/vep_cache/homo_sapiens/91_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --fork 8 --assembly=GRCh38 --offline --tab --plugin Downstream --everything --terms SO --pick --coding_only --transcript_version -i /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.filtered.PASS.vcf -o /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.filtered.PASS.vep.tsv --force_overwrite

#VEP annotate hard-filtered WGS results
#Output VEP VCF
~/bin/ensembl-vep/vep --cache --dir_cache /home/ubuntu/data/vep_cache --dir_plugins /home/ubuntu/data/vep_cache/Plugins --fasta /home/ubuntu/data/vep_cache/homo_sapiens/91_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --fork 8 --assembly=GRCh38 --offline --vcf --plugin Downstream --everything --terms SO --pick --coding_only --transcript_version -i /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.filtered.PASS.vcf -o /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.filtered.PASS.vep.vcf --force_overwrite

#Output tabular VEP
~/bin/ensembl-vep/vep --cache --dir_cache /home/ubuntu/data/vep_cache --dir_plugins /home/ubuntu/data/vep_cache/Plugins --fasta /home/ubuntu/data/vep_cache/homo_sapiens/91_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --fork 8 --assembly=GRCh38 --offline --tab --plugin Downstream --everything --terms SO --pick --coding_only --transcript_version -i /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.filtered.PASS.vcf -o /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.filtered.PASS.vep.tsv --force_overwrite


#VEP annotate VQSR-filtered exome results
#Output VEP VCF
~/bin/ensembl-vep/vep --cache --dir_cache /home/ubuntu/data/vep_cache --dir_plugins /home/ubuntu/data/vep_cache/Plugins --fasta /home/ubuntu/data/vep_cache/homo_sapiens/91_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --fork 8 --assembly=GRCh38 --offline --vcf --plugin Downstream --everything --terms SO --pick --coding_only --transcript_version -i /home/ubuntu/data/germline_variants/Exome_Norm_GGVCFs_jointcalls_recalibrated.PASS.vcf -o /home/ubuntu/data/germline_variants/Exome_Norm_GGVCFs_jointcalls_recalibrated.PASS.vep.vcf --force_overwrite

#Output tabular VEP
~/bin/ensembl-vep/vep --cache --dir_cache /home/ubuntu/data/vep_cache --dir_plugins /home/ubuntu/data/vep_cache/Plugins --fasta /home/ubuntu/data/vep_cache/homo_sapiens/91_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --fork 8 --assembly=GRCh38 --offline --tab --plugin Downstream --everything --terms SO --pick --coding_only --transcript_version -i /home/ubuntu/data/germline_variants/Exome_Norm_GGVCFs_jointcalls_recalibrated.PASS.vcf -o /home/ubuntu/data/germline_variants/Exome_Norm_GGVCFs_jointcalls_recalibrated.PASS.vep.tsv --force_overwrite

#VEP annotate VQSR-filtered WGS results
#Output VEP VCF
~/bin/ensembl-vep/vep --cache --dir_cache /home/ubuntu/data/vep_cache --dir_plugins /home/ubuntu/data/vep_cache/Plugins --fasta /home/ubuntu/data/vep_cache/homo_sapiens/91_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --fork 8 --assembly=GRCh38 --offline --vcf --plugin Downstream --everything --terms SO --pick --coding_only --transcript_version -i /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated.PASS.vcf -o /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated.PASS.vep.vcf --force_overwrite

#Output tabular VEP
~/bin/ensembl-vep/vep --cache --dir_cache /home/ubuntu/data/vep_cache --dir_plugins /home/ubuntu/data/vep_cache/Plugins --fasta /home/ubuntu/data/vep_cache/homo_sapiens/91_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --fork 8 --assembly=GRCh38 --offline --tab --plugin Downstream --everything --terms SO --pick --coding_only --transcript_version -i /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated.PASS.vcf -o /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls_recalibrated.PASS.vep.tsv --force_overwrite

```

