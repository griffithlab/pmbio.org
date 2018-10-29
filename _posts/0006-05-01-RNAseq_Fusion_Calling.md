---
feature_text: |
  ## Precision Medicine
title: RNAseq Fusion Calling
categories:
    - Module-06-RNAseq
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0006-05-01
---

DEV NOTES
- bam -> fastq process not documented, assuming using output from https://pmbio.org/module-10-appendix/0010/02/01/Developer_Notes/ 
- need to update paths the chr 6 and 17 alignment BAMs

---

# Introduction

In addition to providing information about gene expression, RNA-seq data can be used to discover transcripts which result from chromosomal translocations. Translocations and their resultant chimeric (AKA fusion) transcripts are important driver mutations in many cancers. A [variety of specialized alignment and filtering strategies](https://www.ncbi.nlm.nih.gov/pubmed/27485475) have been developed to identify fusion transcripts from RNA, but these programs suffer from low specificity (i.e. many false-positives) and poor correlation across methods.

This tutorial uses the [kallisto](https://pachterlab.github.io/kallisto/about) and [pizzly](https://github.com/pmelsted/pizzly) tools for fusion detection from RNA-seq data. Kallisto quantifies transcript abundance through [pseudoalignment](http://tinyheero.github.io/2015/09/02/pseudoalignments-kallisto.html). Pizzly aligns reads which kallisto has flagged as potentially spanning fusion junctions. Running the tutorial requires RNA fastq files, a reference transcriptome, and a gene annotation file- see below.

# Setup

**Prerequisites- This module assumes you have completed prior modules
including:**
- From [Installation](http://pmbio.org/module-01-setup/0001/04/01/Software_Installation/), installed samtools, Conda package manager, sambamba, R, pizzly, and kallisto.
- From [RNAseq Expression Estimation](http://pmbio.org/module-06-rnaseq/0006/02/01/RNAseq_Expression/), have alignments of normal and tumor RNA as .bam files

###CHECK###
at /workspace/inputs/data/fastqs/chr6_and_chr17/. 

**Additional setup:**

**_Important: pizzly will not work with the most recent Ensembl human GTF annotation file. Download the version 87 GTF as shown in the below code block. Fasta splitting programs which do no preserve the full Ensembl header such as gffread or gtf_to_fasta will not work with pizzly._**

- Download Ensembl GTF and fasta and parse to include only chromosomes 6 and 17 (12 min): 

```bash
# Get files from source
cd /workspace/inputs/reference
mkdir -p ./fusion
cd /workspace/inputs/reference/fusion/
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz       
wget ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz
gunzip -k Homo_sapiens.GRCh38.87.gtf.gz

# Get annotations for only chromosomes 6 and 17
cat Homo_sapiens.GRCh38.87.gtf | grep --color=never -w "^6" >> chr617.gtf
cat Homo_sapiens.GRCh38.87.gtf | grep --color=never -w "^17">> chr617.gtf 

# Parse transcriptome fasta, preserving full ensembl headers
mkdir per-feature
cd per-feature
csplit -s -z ../Homo_sapiens.GRCh38.cdna.all.fa '/>/' '{*}'
for f in xx*; do awk -F ":" 'NR==1{print $2 "." $3 "." $4 "." $5 }' $f | xargs -I{} mv $f {}.fa; done
cd ../
cat ./per-feature/GRCh38.17.*.fa ./per-feature/GRCh38.17.*.fa > chr617.fa
rm -rf per-feature
```

- Merge the chr6_and_chr17 only RNA-seq fastqs:
```bash
mkdir -p /workspace/inputs/data/fastq/chr6_and_chr17/fusion
cd /workspace/inputs/data/fastq/chr6_and_chr17/fusion
cat ../RNAseq_Norm/RNAseq_Norm_Lane1_R1.fastq.gz ../RNAseq_Norm/RNAseq_Norm_Lane2_R1.fastq.gz > RNAseq_Norm_R1.fastq.gz
cat ../RNAseq_Norm/RNAseq_Norm_Lane1_R2.fastq.gz ../RNAseq_Norm/RNAseq_Norm_Lane2_R2.fastq.gz > RNAseq_Norm_R2.fastq.gz
cat ../RNAseq_Tumor/RNAseq_Tumor_Lane1_R1.fastq.gz ../RNAseq_Tumor/RNAseq_Tumor_Lane2_R1.fastq.gz > RNAseq_Tumor_R1.fastq.gz
cat ../RNAseq_Tumor/RNAseq_Tumor_Lane1_R2.fastq.gz ../RNAseq_Tumor/RNAseq_Tumor_Lane2_R2.fastq.gz > RNAseq_Tumor_R2.fastq.gz
```

- Subsample fastqs to allow fusion alignment to run quickly (10 min):
```bash
# Unzip fastqs
gunzip *.fastq.gz

# Use seqtk to take subsamples of the 10% of the fastq read pairs  
seqtk sample -s100 RNAseq_Norm_R1.fastq 0.1 > subRNAseq_Norm_R1.fastq
seqtk sample -s100 RNAseq_Norm_R2.fastq 0.1 > subRNAseq_Norm_R2.fastq
seqtk sample -s100 RNAseq_Tumor_R1.fastq 0.1 > subRNAseq_Tumor_R1.fastq
seqtk sample -s100 RNAseq_Tumor_R2.fastq 0.1 > subRNAseq_Tumor_R2.fastq

```

# Run Fusion Alignment 
- Create kallisto index (4 min):

```bash
mkdir -p /workspace/rnaseq/fusion
cd /workspace/rnaseq/fusion
kallisto index -i index.617.idx -k 31 --make-unique /workspace/inputs/reference/fusion/chr617.fa
```

- Quantify potential fusions (_min):

```bash
kallisto quant -i index.617.idx --fusion -o kquant-norm617 /workspace/inputs/data/fastq/chr6_and_chr17/fusion/subRNAseq_Norm_R1.fastq /workspace/inputs/data/fastq/chr6_and_chr17/fusion/subRNAseq_Norm_R2.fastq

kallisto quant -i index.617.idx --fusion -o kquant-tumor617 /workspace/inputs/data/fastq/chr6_and_chr17/fusion/subRNAseq_Tumor_R1.fastq /workspace/inputs/data/fastq/chr6_and_chr17/fusion/subRNAseq_Tumor_R2.fastq
```

- Call fusions with pizzly (min):

```bash
# (Output is created in the directory where command is run)
cd /workspace/rnaseq/fusion

pizzly -k 31 --gtf /workspace/inputs/reference/fusion/chr617.gtf --cache index-norm617.cache.txt --align-score 2 --insert-size 400 --fasta /workspace/inputs/reference/fusion/chr617.fa --output norm-fuse617 kquant-norm617/fusion.txt

pizzly -k 31 --gtf /workspace/inputs/reference/fusion/chr617.gtf --cache index-tumor617.cache.txt --align-score 2 --insert-size 400 --fasta /workspace/inputs/reference/fusion/chr617.fa --output tumor-fuse617 kquant-tumor617/fusion.txt

```





##########MOVE TO DEV##################
# To generate RNA fusion alignments from full chromosome 6 and 17 RNA-seq data 
**_Machine resources are insufficient to run fusion for the entire transcriptome_**

## Setup

**_Important: pizzly will not work with the most recent Ensembl human GTF annotation file. Download the version 87 GTF as shown in the below code block._** Pizzly also requires a corresponding fasta file with full Ensembl headers. 

- Download Ensembl GTF and fasta:

```bash
mkdir -p /workspace/inputs/references/fusion
cd /workspace/inputs/references/fusion
wget ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
```

- If not already done, download and unpack chr6_and_chr17 RNA fastqs to /workspace/inputs/data/fastqs as described for [Module 2- Data](https://pmbio.org/module-02-inputs/0002/05/01/Data/).

- Get 6-17 sepcific...

wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz
cat

## Run Fusion Alignment 
- Create kallisto index (_min):

```bash
mkdir -p /workspace/rnaseq/fusion
cd /workspace/rnaseq/fusion
kallisto index -i index.idx -k 31 /workspace/inputs/reference/fusion/Homo_sapiens.GRCh38.cdna.all.fa.gz
```
## GET INDEX##
## RESUME HERE ##
- Quantify potential fusions (_min):

```bash
# UPDATE FULL PATH 
kallisto quant -i index.idx --fusion -o kquant-normal-full 
/workspace/inputs/data/fastqs/
/workspace/inputs/data/fastqs/


kallisto quant -i index.idx --fusion -o kquant-tumor-full /workspace/inputs/data/fastqs/ /workspace/inputs/data/fastqs/
```

- Call fusions with pizzly (_min):

```bash
cd /workspace/rnaseq/fusion
# (Output is created in the directory where command is run)
pizzly -k 31 --gtf /workspace/inputs/reference/fusion/Homo_sapiens.GRCh38.87.gtf.gz --cache index-normal-full.cache.txt --align-score 2 --insert-size 400 --fasta /workspace/inputs/reference/fusion/Homo_sapiens.GRCh38.cnda.all.fa.gz --output normal-full kquant-normal-full/fusion.txt
pizzly -k 31 --gtf /workspace/inputs/reference/fusion/Homo_sapiens.GRCh38.87.gtf.gz --cache index-tumor-full.cache.txt --align-score 2 --insert-size 400 --fasta /workspace/inputs/reference/fusion/Homo_sapiens.GRCh38.cnda.all.fa.gz --output tumor-full kquant-normal-full/fusion.txt
```

#### NOTES TO DEL 


create a chr6-17-specific fasta with gtf617:

or
 gtf_to_fasta ./chr6_and_chr17/chr617.gtf Homo_sapiens.GRCh38.dna.toplevel.fa testgtftofasta.fa

samtools faidx Homo_sapiens.GRCh38.dna.toplevel.fa

gffread -w test.fa -g Homo_sapiens.GRCh38.dna.toplevel.fa ./chr6_and_chr17/chr617.gtf


```bash
# Create kallisto index
cd /workspace/rnaseq/fusion
kallisto index -i index..full.idx -k 31 /workspace/inputs/reference/fusion/Homo_sapiens.GRCh38.cnda.all.fa

# Create kallisto index for chromosomes 6 and 17

cd /workspace/inputs/reference/fusion


# Quantify potential fusions
~/workspace/bin/kallisto quant -i /workspace/rnaseq/fusion/index.idx --fusion -o /workspace/rnaseq/fusion/normal-quant /data/RNA_seq/fastqs_RNA/normal-RNA-R1.fastq.gz /data/RNA_seq/fastqs_RNA/normal-RNA-R2.fastq.gz
~/workspace/bin/kallisto quant -i /workspace/rnaseq/fusion/index.idx --fusion -o /workspace/rnaseq/fusion/tumor-quant /data/RNA_seq/fastqs_RNA/tumor-RNA-R1.fastq.gz /data/RNA_seq/fastqs_RNA/tumor-RNA-R2.fastq.gz

# Call fusions with pizzly
# (use unzipped GTF file)
pizzly -k 31 -gtf /data/reference/Homo_sapiens.GRCh38.87.gtf --cache /workspace/rnaseq/fusion/normal-index.cache.txt --align-score 2 --insert-size 400 --fasta /data/reference/Homo_sapiens.GRCh38.cdna.all.fa.gz --output /workspace/rnaseq/fusion/normal-fusion /workspace/rnaseq/fusion/normal-quant/fusion.txt 

pizzly -k 31 -gtf /data/reference/Homo_sapiens.GRCh38.87.gtf --cache /workspace/rnaseq/fusion/tumor-index.cache.txt --align-score 2 --insert-size 400 --fasta /data/reference/Homo_sapiens.GRCh38.cdna.all.fa.gz --output /workspace/rnaseq/fusion/tumor-fusion /workspace/rnaseq/fusion/tumor-quant/fusion.txt 
```



```bash
# Slice merged alignment bams from HISAT2 alignment: 
mkdir /data/RNA_seq/alignments/sliced
~/workspace/bin/sambamba slice -o /data/RNA_seq/alignments/sliced/HCC1395BL-6p.bam /data/RNA_seq/alignments/HCC1395BL.bam chr6:1-57200000
~/workspace/bin/sambamba slice -o /data/RNA_seq/alignments/sliced/HCC1395BL-17p.bam /data/RNA_seq/alignments/HCC1395BL.bam chr17:1-22700000
~/workspace/bin/sambamba slice -o /data/RNA_seq/alignments/sliced/HCC1395-6p.bam /data/RNA_seq/alignments/HCC1395.bam chr6:1-57200000
~/workspace/bin/sambamba slice -o /data/RNA_seq/alignments/sliced/HCC1395-17p.bam /data/RNA_seq/alignments/HCC1395.bam chr17:1-22700000

# Sort
mkdir /data/tmp
for f in /data/RNA_seq/alignments/sliced/*; do unsorted+=("$f"); done
for f in "${unsorted[@]}"; do ~/workspace/bin/sambamba sort --tmpdir /data/tmp -t 4 $f; done

# Merge
~/workspace/bin/sambamba merge -t 4 /data/RNA_seq/alignments/normal-6-17.bam /data/RNA_seq/alignments/sliced/HCC1395BL-6p.sorted.bam /data/RNA_seq/alignments/sliced/HCC1395BL-17p.sorted.bam
~/workspace/bin/sambamba merge -t 4 /data/RNA_seq/alignments/tumnor-6-17.bam /data/RNA_seq/alignments/sliced/HCC1395-6p.sorted.bam /data/RNA_seq/alignments/sliced/HCC1395-17p.sorted.bam


### ON NEW INSTANCE 
### move to dev
### update paths
# Bam to fastq (10 min)
samtools bam2fq /workspace/normal-6-17.bam > /workspace/normal-R-chr6-17.fastq.gz
samtools bam2fq /workspace/tumor-6-17.bam > /workspace/tumor-R-chr6-17.fastq.gz

# Split fastqs (8 min)
cat /workspace/normal-R-chr6-17.fastq.gz | grep '^@.*/1$' -A 3 --no-group-separator > /workspace/normal-R1-chr6-17.fastq.gz
cat /workspace/normal-R-chr6-17.fastq.gz | grep '^@.*/2$' -A 3 --no-group-separator > /workspace/normal-R2-chr6-17.fastq.gz
cat /workspace/tumor-R-chr6-17.fastq.gz | grep '^@.*/1$' -A 3 --no-group-separator > /workspace/tumor-R1-chr6-17.fastq.gz
cat /workspace/tumor-R-chr6-17.fastq.gz | grep '^@.*/2$' -A 3 --no-group-separator > /workspace/tumor-R2-chr6-17.fastq.gz

```
