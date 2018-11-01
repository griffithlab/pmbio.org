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

#### Merging BAMss
Since we have multiple BAMs of each sample that just represent additional data for the same sequence library, we should combine them into a single BAM for convenience before proceeding.

```bash
# Runtime: ~ 8m each merging command
sambamba merge -t 8 /workspace/rnaseq/alignments/RNAseq_Norm.bam /workspace/rnaseq/alignments/RNAseq_Norm_Lane1.bam /workspace/rnaseq/alignments/RNAseq_Norm_Lane2.bam

sambamba merge -t 8 /workspace/rnaseq/alignments/RNAseq_Tumor.bam /workspace/rnaseq/alignments/RNAseq_Tumor_Lane1.bam /workspace/rnaseq/alignments/RNAseq_Tumor_Lane2.bam
```

#### Indexing BAMs
In order to be able to view our BAM files in IGV, as usual we need to index them
```bash
cd  /workspace/rnaseq/alignments/
samtools index RNAseq_Norm.bam
samtools index RNAseq_Tumor.bam

```

#### IGV exercise
Since we have RNA alignments now, we should compare these to the DNA alignments we generated previously. Open IGV and load six BAM files:
* Normal Exome BAM: http://s#.pmbio.org/align/Exome_Norm_sorted_mrkdup_bqsr.bam
* Tumor Exome BAM: http://s#.pmbio.org/align/Exome_Tumor_sorted_mrkdup_bqsr.bam
* Normal WGS BAM: http://s#.pmbio.org/align/WGS_Norm_merged_sorted_mrkdup_bqsr.bam
* Tumor WGS BAM: http://s#.pmbio.org/align/WGS_Tumor_merged_sorted_mrkdup_bqsr.bam
* Normal RNAseq BAM: http://s#.pmbio.org/rnaseq/alignments/RNAseq_Norm.bam
* Tumor RNAseq BAM: http://s#.pmbio.org/rnaseq/alignments/RNAseq_Tumor.bam

Specific exercises:
* Can you find the BRCA1 germline variant we discussed previously? Hint: search the annotated germline VCF for BRCA1 variants that have "HIGH" impact if you have to.
* Is the BRCA germline variant expressed?
* Can you find the TP53 somatic mutation we discussed previously? Hint: load our final merged somatic VCF and visually compare the normal and tumor exome BAMs. 
* Is the TP53 variant expressed?

#### Assembling transcripts from merged bams
We are now going to use `stringtie` to perform a reference guided transcriptome assembly and then determine transcript abundance estimates for those transcript. The so called "reference guided" mode is specified with `-G ref_transcriptome.gtf`.  If you would like to constrain StringTie to just calculate abundance estimates for those transcripts we already konw about (in the GTF) you would also add the `-e` option. This simplifies the output, makes it easier to integrate expression values with variant data for known genes and is faster.  

```bash
# Runtime: ~8min
stringtie -G /workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o /workspace/inputs/references/transcriptome/RNAseq_Tumor.gtf -p 8 -l RNAseq_Tumor /workspace/rnaseq/alignments/RNAseq_Tumor.bam

# Runtime: ~5min
stringtie -G /workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o /workspace/inputs/references/transcriptome/RNAseq_Norm.gtf -p 8 -l RNAseq_Norm /workspace/rnaseq/alignments/RNAseq_Norm.bam

```

#### Merging Transcripts from merged bams
Since we performed transcript discovery on the tumor and normal sample independently, we want to create a unified version of the transcriptome before proceeding to comparing expression between samples. The StringTie `--merge` option is used to combine multiple GTFs into a single GTF.  If we supply our reference transcriptome at the same time (using `-G` ref_transcriptome.gtf) it will also take this information into account. 

```bash
stringtie --merge -p 8 -G /workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o /workspace/inputs/references/transcriptome/stringtie_merged.gtf /workspace/inputs/references/transcriptome/RNAseq_Tumor.gtf /workspace/inputs/references/transcriptome/RNAseq_Norm.gtf

```

#### Comparing transcripts
The `gffcompare` tool can be used to give a more detailed comparison of the transcripts assembled by StringTie and the reference transcriptome.

```bash
mkdir -p /workspace/rnaseq/transcripts/
cd /workspace/rnaseq/transcripts/
gffcompare -r /workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o /workspace/rnaseq/transcripts/gffcmp /workspace/inputs/references/transcriptome/stringtie_merged.gtf

```

#### Estimate Abundance
Below we will estimate abundance using our final combined GTF.  We will do this for each of the individual lanes of RNA-seq data, but also for our merged normal and tumor BAMs.  

```bash
mkdir -p /workspace/rnaseq/ballgown/RNAseq_Tumor_Lane1
mkdir -p /workspace/rnaseq/ballgown/RNAseq_Tumor_Lane2
mkdir -p /workspace/rnaseq/ballgown/RNAseq_Tumor

mkdir -p /workspace/rnaseq/ballgown/RNAseq_Norm_Lane1
mkdir -p /workspace/rnaseq/ballgown/RNAseq_Norm_Lane2
mkdir -p /workspace/rnaseq/ballgown/RNAseq_Norm

cd /workspace/rnaseq/ballgown
# Runtime: ~3min each stringtie command below
stringtie -e -B -G /workspace/rnaseq/transcripts/gffcmp.annotated.gtf -o /workspace/rnaseq/ballgown/RNAseq_Tumor_Lane1/RNAseq_Tumor_Lane1.gtf -p 8 /workspace/rnaseq/alignments/RNAseq_Tumor_Lane1.bam

stringtie -e -B -G /workspace/rnaseq/transcripts/gffcmp.annotated.gtf -o /workspace/rnaseq/ballgown/RNAseq_Tumor_Lane2/RNAseq_Tumor_Lane2.gtf -p 8 /workspace/rnaseq/alignments/RNAseq_Tumor_Lane2.bam

stringtie -e -B -G /workspace/rnaseq/transcripts/gffcmp.annotated.gtf -o /workspace/rnaseq/ballgown/RNAseq_Tumor/RNAseq_Tumor.gtf -p 8 /workspace/rnaseq/alignments/RNAseq_Tumor.bam

stringtie -e -B -G /workspace/rnaseq/transcripts/gffcmp.annotated.gtf -o /workspace/rnaseq/ballgown/RNAseq_Norm_Lane1/RNAseq_Norm_Lane1.gtf -p 8 /workspace/rnaseq/alignments/RNAseq_Norm_Lane1.bam

stringtie -e -B -G /workspace/rnaseq/transcripts/gffcmp.annotated.gtf -o /workspace/rnaseq/ballgown/RNAseq_Norm_Lane2/RNAseq_Norm_Lane2.gtf -p 8 /workspace/rnaseq/alignments/RNAseq_Norm_Lane2.bam

stringtie -e -B -G /workspace/rnaseq/transcripts/gffcmp.annotated.gtf -o /workspace/rnaseq/ballgown/RNAseq_Norm/RNAseq_Norm.gtf -p 8 /workspace/rnaseq/alignments/RNAseq_Norm.bam

```
