---
feature_text: |
  ## Precision Medicine
title: Germline SNV and Indel Calling
categories:
    - Module 05. Germline
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0004-02-01
---

Analysis based on suggestions here:
https://gatkforums.broadinstitute.org/gatk/discussion/7869/howto-discover-variants-with-gatk-a-gatk-workshop-tutorial

### Run GATK HaplotypeCaller

Include option to generate bam output from haplotype caller so that local reassembly/ralignment around called variants can be visualized.

Runtimes: Exome, 160min; 

TO DO: Consider whether the following variant calling steps should be run on alt contigs as well

```bash
cd ~/data
mkdir germline_variants
cd germline_variants

gatk --java-options '-Xmx64g' HaplotypeCaller -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/Exome_Norm_sorted_mrkdup_bqsr.bam -O /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.vcf --bam-output /home/ubuntu/data/germline_variants/Exome_Norm_HC_out.bam -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

gatk --java-options '-Xmx64g' HaplotypeCaller -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/WGS_Norm_merged_sorted_mrkdup_bqsr.bam -O /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.vcf --bam-output /home/ubuntu/data/germline_variants/WGS_Norm_HC_out.bam -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM                                      
```

### Run HaplotypeCaller in GVCF mode with single sample calling, followed by joint calling

TO DO: Consider implementing the following options from CCDG workflow: -ERC GVCF --max_alternate_alleles 3 -variant_index_type LINEAR -variant_index_parameter 128000 -L $CHR -o $chr.g.vcf.gz -contamination $freemix --read_filter OverclippedRead
See: https://confluence.ris.wustl.edu/display/BIO/Proposed+CCDG+Analysis+Workflow+-+2017.7.14

```bash
gatk --java-options '-Xmx64g' HaplotypeCaller -ERC GVCF -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/Exome_Norm_sorted_mrkdup_bqsr.bam -O /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.g.vcf --bam-output /home/ubuntu/data/germline_variants/Exome_Norm_HC_GVCF_out.bam -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

gatk --java-options '-Xmx64g' HaplotypeCaller -ERC GVCF -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/WGS_Norm_merged_sorted_mrkdup_bqsr.bam -O /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.g.vcf --bam-output /home/ubuntu/data/germline_variants/WGS_Norm_HC_GVCF_out.bam -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

```

### Perform joint genotype calling

TO DO: Note, this is a somewhat artificial (and possibly ill-advised?) example where the two gvcfs that we are joint genotyping are wgs and exome calls from the same individual. This is not the normal use of joint genotyping and is provided just as an illustration of how the command should be formulated. For a more realistic example, we should download a few 1000 genome exomes or something.

gatk --java-options '-Xmx64g' GenotypeGVCFs -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.g.vcf -V /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.g.vcf -O Exome_WGS_GGVCFs_jointcalls.vcf -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

