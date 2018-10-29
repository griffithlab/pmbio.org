---
feature_text: |
  ## Precision Medicine
title: RNAseq Expression Estimation
categories:
    - Module-06-RNAseq
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0006-02-01
---

#### **Downloading RNA BAMs**

First, if not done so already, make a separate folder named `/data/RNA_seq` and a subfolder called `fastqs_RNA` and download the RNA_seq data from `genomedata.org` to your instance.
In order to prevent storage space from running out, you may want to unzip the files sequentially and delete the original zipped file once the unzipped file have been obtained.

You will need to have the following software installed, including `HISAT`, `Sambamba`, `StringTie`, `Gffcompare`, `R`. If you are missing any of the following software, or you run into problems with running commands using your currently installed version, please refer to the installation page.

#### **Annotation**
To continue with the annotation step, you will need to first download the proper gtf file named `converted_Homo_sapiens.GRCh38.92.gtf` from genomedata.org. You may need to create a folder in `/data` for storing annotation files by running `mkdir -p /data/annotation`.
```bash
~/bin/hisat2-2.0.4/hisat2_extract_splice_sites.py /data/refseq/converted_Homo_sapiens.GRCh38.92.gtf > /data/annotation/GRCh38_ss.tsv
(/workspace/inputs/references/transcriptome/splicesites.tsv)
~/bin/hisat2-2.0.4/hisat2_extract_exons.py /data/refseq/converted_Homo_sapiens.GRCh38.92.gtf > /data/annotation/GRCh38_exons.tsv
(/workspace/inputs/references/transcriptome/exons.tsv)
```

#### **Indexing**
\* Note that this step may require up to 200 GB of RAM.
```bash
~/bin/hisat2-2.0.4/hisat2-build -p 1 --ss /data/annotation/GRCh38_ss.tsv --exon /data/annotation/GRCh38_exons.tsv /data/refseq/GRCh38_full_analysis_set_plus_decoy_hla.fa /data/refseq/GRCh38_tran
(/workspace/inputs/references/transcriptome/ref_genome)
```

#### **Alignment**
First, we will assign a path for temporary directories:
```bash
mkdir -p /workspace/rnaseq/alignments
cd /workspace/rnaseq/alignments
TUMOR_DATA_1_TEMP=`mktemp -d /workspace/rnaseq/alignments/2895626107.XXXXXXXXXXXX`
TUMOR_DATA_2_TEMP=`mktemp -d /workspace/rnaseq/alignments/2895626112.XXXXXXXXXXXX`
NORMAL_DATA_1_TEMP=`mktemp -d /workspace/rnaseq/alignments/2895625992.XXXXXXXXXXXX`
NORMAL_DATA_2_TEMP=`mktemp -d /workspace/rnaseq/alignments2895626097.XXXXXXXXXXXX`

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
sambamba merge -t 4 /workspace/rnaseq/alignments/HCC1395BL_RNA.bam /workspace/rnaseq/alignments/RNAseq_Norm_Lane1.bam /workspace/rnaseq/alignments/RNAseq_Norm_Lane2.bam

sambamba merge -t 4 /workspace/rnaseq/alignments/HCC1395_RNA.bam /workspace/rnaseq/alignments/RNAseq_Tumor_Lane1.bam /workspace/rnaseq/alignments/RNAseq_Tumor_Lane2.bam
```
#### **Assembling transcript from merged bams**

```bash
~/bin/stringtie -G /data/RNA_seq/refseq/converted_Homo_sapiens.GRCh38.92.gtf -o /data/RNA_seq/transcripts/HCC1395_RNA.gtf -p 4 -l HCC1395_RNA /data/RNA_seq/alignments/HCC1395_RNA.bam

~/bin/stringtie -G /data/RNA_seq/refseq/converted_Homo_sapiens.GRCh38.92.gtf -o /data/RNA_seq/transcripts/HCC1395BL_RNA.gtf -p 4 -l HCC1395BL_RNA /data/RNA_seq/alignments/HCC1395BL_RNA.bam
```
#### **Merging Transcripts from merged bams**

```bash
~/bin/stringtie --merge -p 4 -G /data/RNA_seq/refseq/converted_Homo_sapiens.GRCh38.92.gtf -o /data/RNA_seq/transcripts/stringtie_merged_bams.gtf /data/RNA_seq/transcripts/HCC1395_RNA.gtf $RNA_HOME/transcripts/HCC1395BL_RNA.gtf
```
#### **Comparing transcripts**

```bash
~/bin/gffcompare -r /data/RNA_seq/refseq/converted_Homo_sapiens.GRCh38.92.gtf -o /data/RNA_seq/transcripts/gffcmp /data/RNA_seq/transcripts/stringtie_merged_bams.gtf
```
#### **Estimate Abundance**

```bash
mkdir -p /data/RNA_seq/ballgown/HCC1395_RNA_H3MYFBBXX_4_GCCAAT
mkdir -p /data/RNA_seq/ballgown/HCC1395_RNA_H3MYFBBXX_5_GCCAAT
mkdir -p /data/RNA_seq/ballgown/HCC1395BL_RNA_H3MYFBBXX_4_CTTGTA
mkdir -p /data/RNA_seq/ballgown/HCC1395BL_RNA_H3MYFBBXX_5_CTTGTA

~/bin/stringtie -e -B -G /data/RNA_seq/transcripts/gffcmp.annotated.gtf -o /data/RNA_seq/ballgown/HCC1395_RNA_H3MYFBBXX_4_GCCAAT/HCC1395_RNA_H3MYFBBXX_4_GCCAAT.gtf -p 4 /data/RNA_seq/alignments/HCC1395_RNA_H3MYFBBXX_4_GCCAAT.bam

~/bin/stringtie -e -B -G /data/RNA_seq/transcripts/gffcmp.annotated.gtf -o /data/RNA_seq/ballgown/HCC1395_RNA_H3MYFBBXX_5_GCCAAT/HCC1395_RNA_H3MYFBBXX_5_GCCAAT.gtf -p 4 /data/RNA_seq/alignments/HCC1395_RNA_H3MYFBBXX_5_GCCAAT.bam

~/bin/stringtie -e -B -G /data/RNA_seq/transcripts/gffcmp.annotated.gtf -o /data/RNA_seq/ballgown/HCC1395BL_RNA_H3MYFBBXX_4_CTTGTA/HCC1395BL_RNA_H3MYFBBXX_4_CTTGTA.gtf -p 4 /data/RNA_seq/alignments/HCC1395BL_RNA_H3MYFBBXX_4_CTTGTA.bam

~/bin/stringtie -e -B -G /data/RNA_seq/transcripts/gffcmp.annotated.gtf -o /data/RNA_seq/ballgown/HCC1395BL_RNA_H3MYFBBXX_5_CTTGTA/HCC1395BL_RNA_H3MYFBBXX_5_CTTGTA.gtf -p 4 /data/RNA_seq/alignments/HCC1395BL_RNA_H3MYFBBXX_5_CTTGTA.bam
```
