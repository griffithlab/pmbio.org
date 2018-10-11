---
feature_text: |
  ## Precision Medicine
title: Germline SNV and Indel Calling
categories:
    - Module-05-Germline
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-02-01
---

### Acknowledgements and citations

- GATK

### Module objectives

- Perform single sample germline variant calling with GATK on WGS and Exome data
- Perform single sample germline variant calling GATK GVCF workflow on Exome data for HCC1395 and 1000genome exomes

In this module we will use the GATK HaplotypeCaller to call variants from our aligned bams. Since we are only interested in germline variants in this module, we will only call variants in the normal samples (WGS and exome). The following tutorial and example commands are based on suggestions from following [GATK tutorial](https://gatkforums.broadinstitute.org/gatk/discussion/7869/howto-discover-variants-with-gatk-a-gatk-workshop-tutorial), provided by the Broad Institute. 

### Run GATK HaplotypeCaller

Include option to generate bam output from haplotype caller so that local reassembly/realignment around called variants can be visualized.

Runtimes: Exome, 160min; 

NOTE: In the following command we have limited to calling variants on chr1-22, X, Y and MT. You may wish to also run on alt contigs.

```bash
cd ~/data
mkdir germline_variants
cd germline_variants

#Call variants for exome data
gatk --java-options '-Xmx64g' HaplotypeCaller -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/Exome_Norm_sorted_mrkdup_bqsr.bam -O /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.vcf --bam-output /home/ubuntu/data/germline_variants/Exome_Norm_HC_out.bam -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

#Call variants for WGS data
gatk --java-options '-Xmx64g' HaplotypeCaller -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/WGS_Norm_merged_sorted_mrkdup_bqsr.bam -O /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.vcf --bam-output /home/ubuntu/data/germline_variants/WGS_Norm_HC_out.bam -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM                                      

```

### Run HaplotypeCaller in GVCF mode with single sample calling, followed by joint calling (for exomes)

TO DO: Consider implementing the following options from CCDG workflow: -ERC GVCF --max_alternate_alleles 3 -variant_index_type LINEAR -variant_index_parameter 128000 -L $CHR -o $chr.g.vcf.gz -contamination $freemix --read_filter OverclippedRead
See: https://confluence.ris.wustl.edu/display/BIO/Proposed+CCDG+Analysis+Workflow+-+2017.7.14

```bash

#Call variants in GVCF mode for exome data
gatk --java-options '-Xmx64g' HaplotypeCaller -ERC GVCF -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/Exome_Norm_sorted_mrkdup_bqsr.bam -O /home/ubuntu/data/germline_variants/Exome_Norm_HC_calls.g.vcf --bam-output /home/ubuntu/data/germline_variants/Exome_Norm_HC_GVCF_out.bam -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

#Call variants in GVCF mode for WGS data
gatk --java-options '-Xmx64g' HaplotypeCaller -ERC GVCF -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/WGS_Norm_merged_sorted_mrkdup_bqsr.bam -O /home/ubuntu/data/germline_variants/WGS_Norm_HC_calls.g.vcf --bam-output /home/ubuntu/data/germline_variants/WGS_Norm_HC_GVCF_out.bam -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

```


### Obtain 1000 genomes exome bam files for joint genotyping and VQSR steps

Recall that the cell line being used in this tutorial ([HCC1395](https://www.atcc.org/Products/All/CRL-2324.aspx)) was derived from a 43 year old caucasian female (by a research group in Dallas Texas). Therefore, for a "matched" set of reference alignments we might limit to those of European descent.    

#### Get a list of 1000 genome GBR samples

- See summary of samples/ethnicities available: http://www.internationalgenome.org/data-portal/population
- Now, we are guessing, but perhaps GBR (British in England and Scotland) would be the best match
- Go to: http://www.internationalgenome.org/data-portal/sample
- Filter by population -> GBR
- Filter by analysis group -> Exome
- Filter by data collection -> 1000 genomes on GRCh38
- Download the list to get sample details (igsr_samples_GBR.tsv)

#### Get a list of data files for GBR samples

- Go to Populations (http://www.internationalgenome.org/data-portal/population)
- Select GBR
- Scroll down to 'Data collections for GBR'
- Choose '1000 Genomes on GRCh38' tab
- Select 'Data types' -> 'Alignment'
- Select 'Analysis groups' -> 'Exome'
- 'Download the list' (igsr_GBR_GRCh38_exome_alignment.tsv). 

#### Download exome cram and crai files for 30 female GBR cases

```

cd /home/ubuntu/data/reference/
mkdir 1000genomes
cd 1000genomes
wget http://genomedata.org/pmbio-workshop/references/igsr_samples_GBR.tsv
wget http://genomedata.org/pmbio-workshop/references/igsr_GBR_GRCh38_exome_alignment.tsv 
grep "female" igsr_samples_GBR.tsv | cut -f 1 | head -30 > GBR_female_30_samples.txt
grep -f GBR_female_30_samples.txt igsr_GBR_GRCh38_exome_alignment.tsv | cut -f 1 | grep ".cra" | perl -ne 'chomp; print "wget $_\n";' > wget_files.sh
bash wget_files.sh

```

Note: These cram files were created in a generally compatible way with the anlysis done so far in this tutorial. The were aligned with BWA-mem to GRCh38 with alternative sequences, plus decoys and HLA. This was followed by GATK BAM improvement steps as in the 1000 Genomes phase 3 pipeline (GATK IndelRealigner, BQSR, and Picard MarkDuplicates). Of course, there will be batch effects related to potentially different sample preparation, library construction, exome capture reagent and protocol, sequencing pipeline, etc. For more details see:

- http://www.internationalgenome.org/analysis
- https://media.nature.com/original/nature-assets/nature/journal/v526/n7571/extref/nature15393-s1.pdf

#### Perform germline variant calling on 1KG exomes with GATK in GVCF mode

```

gatk --java-options '-Xmx64g' HaplotypeCaller -ERC GVCF -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /data/reference/1000genomes/HG00099.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram -O /data/germline_variants/HG00099_HC_calls.g.vcf -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

gatk --java-options '-Xmx64g' HaplotypeCaller -ERC GVCF -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /data/reference/1000genomes/HG00102.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram -O /data/germline_variants/HG00102_HC_calls.g.vcf -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

gatk --java-options '-Xmx64g' HaplotypeCaller -ERC GVCF -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /data/reference/1000genomes/HG00104.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram -O /data/germline_variants/HG00104_HC_calls.g.vcf -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

gatk --java-options '-Xmx64g' HaplotypeCaller -ERC GVCF -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /data/reference/1000genomes/HG00106.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram -O /data/germline_variants/HG00106_HC_calls.g.vcf -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

gatk --java-options '-Xmx64g' HaplotypeCaller -ERC GVCF -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /data/reference/1000genomes/HG00118.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram -O /data/germline_variants/HG00118_HC_calls.g.vcf -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

```

### Perform joint genotype calling

Create joint genotype vcf, combining HCC1395 Exome normal together with a set of 1KG exomes, for use in VQSR filtering.

```
#Combine gvcfs into a single vcf for use with GenotypeGVCFs
gatk --java-options '-Xmx64g' CombineGVCFs -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /data/germline_variants/Exome_Norm_HC_calls.g.vcf -V /data/germline_variants/HG00099_HC_calls.g.vcf -V /data/germline_variants/HG00102_HC_calls.g.vcf -V /data/germline_variants/HG00104_HC_calls.g.vcf -V /data/germline_variants/HG00106_HC_calls.g.vcf -V /data/germline_variants/HG00118_HC_calls.g.vcf -O /data/germline_variants/Exome_Norm_1KG_HC_calls_combined.g.vcf -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

#Perform joint genotyping
gatk --java-options '-Xmx64g' GenotypeGVCFs -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -V /data/germline_variants/Exome_Norm_1KG_HC_calls_combined.g.vcf -O /data/germline_variants/Exome_GGVCFs_jointcalls.vcf -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

```
