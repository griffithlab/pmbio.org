---
feature_text: |
  ## Precision Medicine
title: RNAseq Expression Estimation
categories:
    - Module-06-RNAseq
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0006-02-01
---

### Annotation
To continue with the annotation step, you will need to first download the proper gtf file named `converted_Homo_sapiens.GRCh38.92.gtf` from genomedata.org. You may need to create a folder in `/data` for storing annotation files by running `mkdir -p /data/annotation`.
```bash
~/bin/hisat2-2.0.4/hisat2_extract_splice_sites.py /data/refseq/converted_Homo_sapiens.GRCh38.92.gtf > /data/annotation/GRCh38_ss.tsv
(/workspace/inputs/references/transcriptome/splicesites.tsv)
~/bin/hisat2-2.0.4/hisat2_extract_exons.py /data/refseq/converted_Homo_sapiens.GRCh38.92.gtf > /data/annotation/GRCh38_exons.tsv
(/workspace/inputs/references/transcriptome/exons.tsv)
```

### Indexing

\* Note that this step may require up to 200 GB of RAM.
```bash
~/bin/hisat2-2.0.4/hisat2-build -p 1 --ss /data/annotation/GRCh38_ss.tsv --exon /data/annotation/GRCh38_exons.tsv /data/refseq/GRCh38_full_analysis_set_plus_decoy_hla.fa /data/refseq/GRCh38_tran
(/workspace/inputs/references/transcriptome/ref_genome)
```
### Adapter Trimming FASTQ files


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
