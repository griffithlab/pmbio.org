---
feature_text: |
  ## Precision Medicine
title: RNAseq Fusion Calling
categories:
    - Module-06-RNAseq
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0006-05-01
---

# Introduction

In addition to providing information about gene expression, RNA-seq data can be used to discover transcripts which result from chromosomal translocations. Translocations and their resultant chimeric (AKA fusion) transcripts are important driver mutations in many cancers. A variety of specialized RNA alignment and filtering strategies have been developed to identify fusion transcripts, but these programs suffer from low specificity (many false-positives) and poor correlation across programs.

This tutorial uses the [kallisto](https://pachterlab.github.io/kallisto/about) and [pizzly](https://github.com/pmelsted/pizzly) tools for fusion detection. kallisto quantifies transcript abundance through pseudoalignment. pizzly aligns reads which kallisto has flagged as potentially spanning fusion junctions. In addition to RNA fastq files, fusion calling requires a reference transcriptome and gene annotation file- see below.

# Setup

**Prerequisites- This module assumes you have completed prior modules
including:**
- From [Installation](http://pmbio.org/module-01-setup/0001/04/01/Software_Installation/), installed Conda package manager, conda-forge and bioconda channels, pizzly, and kallisto.
- From [Data](http://pmbio.org/module-02-inputs/0002/05/01/Data/), downloaded and unzipped RNA-seq fastq .tar files from <http://genomedata.org/pmbio-workshop/fastqs/all/> to \~/data/fastqs

**Additional setup:**
- Download Ensembl GTF and cDNA files:
```bash
cd ~/data/reference
wget ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
```

- IF NECESSARY merge rna seq files:
```bash
cd ~/data/fastqs/RNAseq_Tumor
cat RNAseq_Tumor_Lane1_R1.fastq.gz RNAseq_Tumor_Lane2_R1.fastq.gz > RNAseq_TumorR1.fastq.gz
cat RNAseq_Tumor_Lane1_R2.fastq.gz RNAseq_Tumor_Lane2_R2.fastq.gz > RNAseq_TumorR2.fastq.gz
cd ~/data/fastqs/RNAseq_Norm
cat RNAseq_Norm_Lane1_R1.fastq.gz RNAseq_Norm_Lane2_R1.fastq.gz > RNAseq_NormalR1.fastq.gz
cat RNAseq_Norm_Lane1_R2.fastq.gz RNAseq_Norm_Lane2_R2.fastq.gz > RNAseq_NormalR2.fastq.gz
```
- To save space, you may wish to remove the original zipped fastq files as well as the unzipped but unmerged fastq files.


# Run Fusion Alignment and Filtering
- Create kallisto index:

```bash
cd /data
mkdir rna-fusion
cd rna-fusion
kallisto index -i index.idx -k 31 ~/data/reference/Homo_sapiens.GRCh38.cdna.all.fa.gz
```

- Quantify potential fusions:

```bash
cd /data/rna-fusion
kallisto quant -i index.idx --fusion -o output-tumor ~/data/fastqs/RNAseq_Tumor/RNAseq_TumorR1.fastq.gz ~/data/fastqs/RNAseq_Tumor/RNAseq_TumorR2.fastq.gz
kallisto quant -i index.idx --fusion -o output-normal ~/data/fastqs/RNAseq_Norm/RNAseq_NormalR1.fastq.gz ~/data/fastqs/RNAseq_Norm/RNAseq_NormalR2.fastq.gz
```

- Filter fusions with pizzly:

```bash
cd /data/rna-fusion
pizzly -k 31 --gtf /data/annotation/Homo_sapiens.GRCh38.93.gtf --cache index.cache.txt --align-score 2 -ignore-protein --insert-size 400 --fasta /data/annotation/Homo_sapiens.GRCh38.cdna.all.fa.gz --output fusTumor output-tumor/fusion.txt

pizzly -k 31 --gtf /data/annotation/Homo_sapiens.GRCh38.93.gtf --cache index.cache.txt --align-score 2 -ignore-protein --insert-size 400 --fasta /data/annotation/Homo_sapiens.GRCh38.cdna.all.fa.gz --output fusNormal output-normal/fusion.txt
```
