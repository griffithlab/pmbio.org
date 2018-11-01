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
The purpose of adapter trimming is to remove sequences in our data that correspond to the Illumina sequence adapters.  The most common adapter trimming scenario is the removal of adapter sequences that occur at the end of read sequences. This happens when a DNA (or cDNA) fragment is shorter than the read length.  For example if we sequence our RNA-seq fragments to 150 base length and a fragment is only 140 bases long the read will end with 10 bases of adapter sequence. Since this adapter sequence does not correspond to the genome, it will not align. Too much adapter sequence can actually prevent reads from aligning at all. Adapter trimming may therefore sometime improve the overall alignment success rate for an RNA-seq data set.  Adapter trimming involves a simplistic alignment itself and therefore can be computationally expensive.
```bash
cd ~/workspace/inputs/references
wget -c http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa
cd ~/workspace/inputs/data/fastq/RNAseq_Tumor
#runtime: ~25min each flexbar command
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

# Runtime: ~15 min each run
hisat2 -p 8 --dta -x /workspace/inputs/references/transcriptome/ref_genome --rg-id 2895626107 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-GCCAAT.4 --rg LB:rna_tumor_lib1 --rg SM:HCC1395_RNA --rna-strandness RF -1 /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane1_R1.fastq.gz -2  /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane1_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 8 -m 32G --tmpdir $TUMOR_DATA_1_TEMP -o /workspace/rnaseq/alignments/RNAseq_Tumor_Lane1.bam /dev/stdin

rmdir $TUMOR_DATA_1_TEMP/* $TUMOR_DATA_1_TEMP

hisat2 -p 8 --dta -x /workspace/inputs/references/transcriptome/ref_genome --rg-id 2895626112 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-GCCAAT.5 --rg LB:rna_tumor_lib1 --rg SM:HCC1395_RNA --rna-strandness RF -1 /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane2_R1.fastq.gz -2  /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane2_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 8 -m 32G --tmpdir $TUMOR_DATA_2_TEMP -o /workspace/rnaseq/alignments/RNAseq_Tumor_Lane2.bam /dev/stdin

rmdir $TUMOR_DATA_2_TEMP/* $TUMOR_DATA_2_TEMP

hisat2 -p 8 --dta -x /workspace/inputs/references/transcriptome/ref_genome --rg-id 2895625992 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-CTTGTA.4 --rg LB:rna_norm_lib1 --rg SM:HCC1395BL_RNA --rna-strandness RF -1 /workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane1_R1.fastq.gz -2  /workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane1_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 8 -m 32G --tmpdir $NORMAL_DATA_1_TEMP -o /workspace/rnaseq/alignments/RNAseq_Norm_Lane1.bam /dev/stdin

rmdir $NORMAL_DATA_1_TEMP/* $NORMAL_DATA_1_TEMP

hisat2 -p 8 --dta -x /workspace/inputs/references/transcriptome/ref_genome --rg-id 2895626097 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-CTTGTA.5 --rg LB:rna_norm_lib1 --rg SM:HCC1395BL_RNA --rna-strandness RF -1 /workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane2_R1.fastq.gz -2  /workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane2_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 8 -m 32G --tmpdir $NORMAL_DATA_1_TEMP -o /workspace/rnaseq/alignments/RNAseq_Norm_Lane2.bam /dev/stdin

rmdir $NORMAL_DATA_2_TEMP/* $NORMAL_DATA_2_TEMP

```

#### Merging Bams

```bash
# Runtime: ~ 8m each merging command
sambamba merge -t 8 /workspace/rnaseq/alignments/RNAseq_Norm.bam /workspace/rnaseq/alignments/RNAseq_Norm_Lane1.bam /workspace/rnaseq/alignments/RNAseq_Norm_Lane2.bam

sambamba merge -t 8 /workspace/rnaseq/alignments/RNAseq_Tumor.bam /workspace/rnaseq/alignments/RNAseq_Tumor_Lane1.bam /workspace/rnaseq/alignments/RNAseq_Tumor_Lane2.bam
```
#### Assembling transcript from merged bams

```bash
# Runtime: ~8min
stringtie -G /workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o /workspace/inputs/references/transcriptome/RNAseq_Tumor.gtf -p 8 -l RNAseq_Tumor /workspace/rnaseq/alignments/RNAseq_Tumor.bam

stringtie -G /workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o /workspace/inputs/references/transcriptome/RNAseq_Norm.gtf -p 8 -l RNAseq_Norm /workspace/rnaseq/alignments/RNAseq_Norm.bam

```

#### Merging Transcripts from merged bams

```bash
~/bin/stringtie --merge -p 8 -G /workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o /workspace/inputs/references/transcriptome/stringtie_merged_bams.gtf /workspace/rnaseq/alignments/RNAseq_Tumor.gtf /workspace/inputs/references/transcriptome/RNAseq_Norm.gtf

```

#### Comparing transcripts

```bash
mkdir -p /workspace/rnaseq/transcripts/
cd /workspace/rnaseq/transcripts/
gffcompare -r /workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o /workspace/rnaseq/transcripts/gffcmp /workspace/inputs/references/transcriptome/stringtie_merged_bams.gtf

```

#### Estimate Abundance

```bash
mkdir -p /workspace/rnaseq/ballgown/RNAseq_Tumor_Lane1
mkdir -p /workspace/rnaseq/ballgown/RNAseq_Tumor_Lane2
mkdir -p /workspace/rnaseq/ballgown/RNAseq_Norm_Lane1
mkdir -p /workspace/rnaseq/ballgown/RNAseq_Norm_Lane2

stringtie -e -B -G /workspace/rnaseq/transcripts/gffcmp.annotated.gtf -o /workspace/rnaseq/ballgown/RNAseq_Tumor_Lane1/RNAseq_Tumor_Lane1.gtf -p 8 /workspace/rnaseq/alignments/RNAseq_Tumor_Lane1.bam

stringtie -e -B -G /workspace/rnaseq/transcripts/gffcmp.annotated.gtf -o /workspace/rnaseq/ballgown/RNAseq_Tumor_Lane2/RNAseq_Tumor_Lane2.gtf -p 8 /workspace/rnaseq/alignments/RNAseq_Tumor_Lane2.bam

stringtie -e -B -G /workspace/rnaseq/transcripts/gffcmp.annotated.gtf -o /workspace/rnaseq/ballgown/RNAseq_Norm_Lane1/RNAseq_Norm_Lane1.gtf -p 8 /workspace/rnaseq/alignments/RNAseq_Norm_Lane1.bam

stringtie -e -B -G /workspace/rnaseq/transcripts/gffcmp.annotated.gtf -o /workspace/rnaseq/ballgown/RNAseq_Norm_Lane2/RNAseq_Tumor_Lane2.gtf -p 8 /workspace/rnaseq/alignments/RNAseq_Norm_Lane2.bam

```
