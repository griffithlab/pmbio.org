---
feature_text: |
  ## Precision Medicine
title: Germline SNV and Indel Calling
categories:
    - Module 4
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0004-02-01
---

Analysis based on suggestions here:
https://gatkforums.broadinstitute.org/gatk/discussion/7869/howto-discover-variants-with-gatk-a-gatk-workshop-tutorial

Run GATK HaplotypeCaller
Include option to generate bam output from haplotype caller so that local reassembly/ralignment around called variants can be visualized.

Runtimes: Exome, 160min; 

```bash
cd ~/data
mkdir germline_variants
cd germline_variants

gatk --java-options '-Xmx64g' HaplotypeCaller -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/Exome_Norm_sorted_mrkdup_bqsr.bam -O /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.vcf --bam-output /home/ubuntu/data/germline_variants/Exome_Norm_HC_out.bam -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 

gatk --java-options '-Xmx64g' HaplotypeCaller -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/WGS_Norm_merged_sorted_mrkdup_bqsr.bam -O /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.vcf --bam-output /home/ubuntu/data/germline_variants/WGS_Norm_HC_out.bam -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22                                      
```

######
#Note - above commands should probably be rerun with chrX and chrY?
# Or maybe just all chromosomes???
######

Also, consider switching to GVCF method with single calling followed by joint calling

Here are options from CCDG workflow: -ERC GVCF --max_alternate_alleles 3 -variant_index_type LINEAR -variant_index_parameter 128000 -L $CHR -o $chr.g.vcf.gz -contamination $freemix --read_filter OverclippedRead


