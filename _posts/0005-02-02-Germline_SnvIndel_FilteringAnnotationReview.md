---
feature_text: |
  ## Precision Medicine
title: Germline SNV/Indel Filtering/Annotation/Review
categories:
    - Module-05-Germline
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-02-02
---

### Module objectives

- Perform filtering of germline SNVs and indels 

The raw output of GATK HaplotypeCaller will include many variants with varying degrees of quality. For various reasons we might wish to further filter these to a higher confident set of variants. The recommended approach is to use GATK VQSR. However, this requires a large (i.e., at least 30-50), preferably platform-matched (i.e., similar sequencing strategy), set of samples with variant calls. For our purposes, we will first demonstrate the less optimal hard filtering strategy. 

It is strongly recommended to read the following documentation from GATK:
- [How to run VQSR](https://software.broadinstitute.org/gatk/documentation/article?id=2805)
- [How to apply hard filters](https://software.broadinstitute.org/gatk/documentation/article?id=2806) 

### Perform hard filtering

```
# Extract the SNPs from the call set
gatk --java-options '-Xmx64g' SelectVariants -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.vcf -select-type SNP -O /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.snps.vcf

# Extract the indels from the call set
gatk --java-options '-Xmx64g' SelectVariants -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.vcf -select-type INDEL -O /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.indels.vcf

# Apply basic filters to the SNP call set
gatk --java-options '-Xmx64g' VariantFiltration -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.snps.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "basic_snp_filter" -O /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.snps.filtered.vcf

# Apply basic filters to the INDEL call set
gatk --java-options '-Xmx64g' VariantFiltration -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.indels.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "basic_indel_filter" -O /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.indels.filtered.vcf

# Merge filtered SNP and INDEL vcfs back together
gatk --java-options '-Xmx64g' MergeVcfs -I /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.snps.filtered.vcf -I /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.indels.filtered.vcf -O /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.filtered.vcf

```

### Perform VEP annotation of filtered results
~/bin/ensembl-vep/vep --cache --dir_cache /home/ubuntu/data/vep_cache --vcf -i /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.filtered.vcf -o /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.filtered.vep.vcf 



Note, for running VQSR on exome data. There will be a smaller number of variants per sample compared to WGS. These are typically insufficient to build a robust recalibration model. If running on only a few samples, GATK recommends that you analyze samples jointly in cohorts of at least 30 samples. If necessary, add exomes from 1000G Project or comparable. These should be processed with similar technical generation (technology, capture, read length, depth). 

See also: https://drive.google.com/drive/folders/0BzI1CyccGsZiWlU5SXNvNnFGbVE (GATK Workshops -> 1709 -> GATKwr19-05-Variant_filtering.pdf)




