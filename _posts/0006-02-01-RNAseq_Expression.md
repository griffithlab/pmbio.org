---
feature_text: |
  ## Precision Medicine
title: RNAseq Expression Estimation
categories:
    - Module-06-RNAseq
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0006-02-01
---


### Adapter Trimming FASTQ files
```bash
cd ~/workspace/inputs/references
wget -c http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa
cd ~/workspace/inputs/data/fastq/RNAseq_Tumor

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters ~/workspace/inputs/references/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads RNAseq_Tumor_Lane1_R1.fastq.gz --reads2 RNAseq_Tumor_Lane1_R2.fastq.gz --target ~/workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane1
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters ~/workspace/inputs/references/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads RNAseq_Tumor_Lane2_R1.fastq.gz --reads2 RNAseq_Tumor_Lane2_R2.fastq.gz --target ~/workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane2
```

#### Alignment
First, we will assign a path for temporary directories:
```bash
mkdir -p /workspace/rnaseq/alignments
cd /workspace/rnaseq/alignments
TUMOR_DATA_1_TEMP=`mktemp -d /workspace/rnaseq/alignments/2895626107.XXXXXXXXXXXX`
TUMOR_DATA_2_TEMP=`mktemp -d /workspace/rnaseq/alignments/2895626112.XXXXXXXXXXXX`
NORMAL_DATA_1_TEMP=`mktemp -d /workspace/rnaseq/alignments/2895625992.XXXXXXXXXXXX`
NORMAL_DATA_2_TEMP=`mktemp -d /workspace/rnaseq/alignments/2895626097.XXXXXXXXXXXX`

hisat2 -p 4 --dta -x /workspace/inputs/references/transcriptome/ref_genome --rg-id 2895626107 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-GCCAAT.4 --rg LB:rna_tumor_lib1 --rg SM:HCC1395_RNA --rna-strandness RF -1 /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane1_R1.fastq.gz -2  /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane1_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 4 -m 8G --tmpdir $TUMOR_DATA_1_TEMP -o /workspace/rnaseq/alignments/RNAseq_Tumor_Lane1.bam /dev/stdin

rmdir $TUMOR_DATA_1_TEMP/* $TUMOR_DATA_1_TEMP

hisat2 -p 4 --dta -x /workspace/inputs/references/transcriptome/ref_genome --rg-id 2895626112 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-GCCAAT.5 --rg LB:rna_tumor_lib1 --rg SM:HCC1395_RNA --rna-strandness RF -1 /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane2_R1.fastq.gz -2  /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane2_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | ~/bin/sambamba sort -t 4 -m 8G --tmpdir $TUMOR_DATA_2_TEMP -o /workspace/rnaseq/alignments/RNAseq_Tumor_Lane2.bam /dev/stdin

rmdir $TUMOR_DATA_2_TEMP/* $TUMOR_DATA_2_TEMP

hisat2 -p 4 --dta -x /workspace/inputs/references/transcriptome/ref_genome --rg-id 2895625992 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-CTTGTA.4 --rg LB:rna_norm_lib1 --rg SM:HCC1395BL_RNA --rna-strandness RF -1 /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Norm_Lane1_R1.fastq.gz -2  /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Norm_Lane1_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 4 -m 8G --tmpdir $NORMAL_DATA_1_TEMP -o /workspace/rnaseq/alignments/RNAseq_Norm_Lane1.bam /dev/stdin

rmdir $NORMAL_DATA_1_TEMP/* $NORMAL_DATA_1_TEMP

hisat2 -p 4 --dta -x /workspace/inputs/references/transcriptome/ref_genome --rg-id 2895626097 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-CTTGTA.5 --rg LB:rna_norm_lib1 --rg SM:HCC1395BL_RNA --rna-strandness RF -1 /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Norm_Lane2_R1.fastq.gz -2  /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Norm_Lane2_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 4 -m 8G --tmpdir $NORMAL_DATA_1_TEMP -o /workspace/rnaseq/alignments/RNAseq_Norm_Lane2.bam /dev/stdin

rmdir $NORMAL_DATA_2_TEMP/* $NORMAL_DATA_2_TEMP
```
#### **Merging Bams**

```bash
sambamba merge -t 4 /workspace/rnaseq/alignments/RNAseq_Norm.bam /workspace/rnaseq/alignments/RNAseq_Norm_Lane1.bam /workspace/rnaseq/alignments/RNAseq_Norm_Lane2.bam

sambamba merge -t 4 /workspace/rnaseq/alignments/RNAseq_Tumor.bam /workspace/rnaseq/alignments/RNAseq_Tumor_Lane1.bam /workspace/rnaseq/alignments/RNAseq_Tumor_Lane2.bam
```
#### **Assembling transcript from merged bams**

```bash
~/bin/stringtie -G /workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o /workspace/inputs/references/transcriptome/RNAseq_Tumor.gtf -p 4 -l RNAseq_Tumor /workspace/rnaseq/alignments/RNAseq_Tumor.bam

~/bin/stringtie -G /workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o /workspace/inputs/references/transcriptome/RNAseq_Norm.gtf -p 4 -l RNAseq_Norm /workspace/rnaseq/alignments/RNAseq_Norm.bam
```
#### **Merging Transcripts from merged bams**

```bash
~/bin/stringtie --merge -p 4 -G /workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o /workspace/inputs/references/transcriptome/stringtie_merged_bams.gtf /workspace/rnaseq/alignments/RNAseq_Tumor.gtf /workspace/inputs/references/transcriptome/RNAseq_Norm.gtf
```
#### **Comparing transcripts**

```bash
mkdir -p /workspace/rnaseq/transcripts/
cd /workspace/rnaseq/transcripts/
gffcompare -r /workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o /workspace/rnaseq/transcripts/gffcmp /workspace/inputs/references/transcriptome/stringtie_merged_bams.gtf
```
#### **Estimate Abundance**

```bash
mkdir -p /workspace/rnaseq/ballgown/RNAseq_Tumor_Lane1
mkdir -p /workspace/rnaseq/ballgown/RNAseq_Tumor_Lane2
mkdir -p /workspace/rnaseq/ballgown/RNAseq_Norm_Lane1
mkdir -p /workspace/rnaseq/ballgown/RNAseq_Norm_Lane2

stringtie -e -B -G /workspace/rnaseq/transcripts/gffcmp.annotated.gtf -o /workspace/rnaseq/ballgown/RNAseq_Tumor_Lane1/RNAseq_Tumor_Lane1.gtf -p 4 /workspace/rnaseq/alignments/RNAseq_Tumor_Lane1.bam

stringtie -e -B -G /workspace/rnaseq/transcripts/gffcmp.annotated.gtf -o /workspace/rnaseq/ballgown/RNAseq_Tumor_Lane2/RNAseq_Tumor_Lane2.gtf -p 4 /workspace/rnaseq/alignments/RNAseq_Tumor_Lane2.bam

stringtie -e -B -G /workspace/rnaseq/transcripts/gffcmp.annotated.gtf -o /workspace/rnaseq/ballgown/RNAseq_Norm_Lane1/RNAseq_Norm_Lane1.gtf -p 4 /workspace/rnaseq/alignments/RNAseq_Norm_Lane1.bam

stringtie -e -B -G /workspace/rnaseq/transcripts/gffcmp.annotated.gtf -o /workspace/rnaseq/ballgown/RNAseq_Norm_Lane2/RNAseq_Tumor_Lane2.gtf -p 4 /workspace/rnaseq/alignments/RNAseq_Norm_Lane2.bam
```
