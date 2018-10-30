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

In addition to providing information about gene expression, RNA-seq data can be used to discover transcripts which result from chromosomal translocations. Translocations and their resultant chimeric (AKA fusion) transcripts are important driver mutations in many cancers. A [variety of specialized alignment and filtering strategies](https://www.ncbi.nlm.nih.gov/pubmed/27485475) have been developed to identify fusion transcripts from RNA, but these programs suffer from low specificity (i.e. many false-positives) and poor correlation across methods.

This tutorial uses the [kallisto](https://pachterlab.github.io/kallisto/about) and [pizzly](https://github.com/pmelsted/pizzly) tools for fusion detection from RNA-seq data. Kallisto quantifies transcript abundance through [pseudoalignment](http://tinyheero.github.io/2015/09/02/pseudoalignments-kallisto.html). Pizzly aligns reads which kallisto has flagged as potentially spanning fusion junctions. Running the tutorial requires RNA fastq files, a reference transcriptome, and a gene annotation file- see below.

# Setup

**Prerequisites- This module assumes you have completed prior modules
including:**
- From [Installation](http://pmbio.org/module-01-setup/0001/04/01/Software_Installation/), installed Conda package manager, R, pizzly, and kallisto.
- From [Data](http://127.0.0.1:4000/module-02-inputs/0002/05/01/Data/), have fastq sequence files of normal and tumor RNA at /workspace/inputs/data/fastq/chr6_and_chr17/. 

**Additional setup:**

**_Important: pizzly will not work with the most recent Ensembl human GTF annotation file. Download the version 87 GTF as shown in the below code block. We will subset the reference transcriptome fasta file. Fasta splitting programs which do no preserve the full Ensembl header such as gffread or gtf_to_fasta will not work with pizzly._**

- Download Ensembl GTF and fasta and parse to include only chromosomes 6 and 17 (12 min): 

```bash
# Get files from source
mkdir -p workspace/inputs/reference/fusion
cd /workspace/inputs/reference/fusion/
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz       
wget ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz
gunzip -k Homo_sapiens.GRCh38.87.gtf.gz

# Get annotations for only chromosomes 6 and 17
cat Homo_sapiens.GRCh38.87.gtf | grep --color=never -w "^6" >> chr617.gtf
cat Homo_sapiens.GRCh38.87.gtf | grep --color=never -w "^17">> chr617.gtf 

# Parse transcriptome fasta, preserving full ensembl headers
#  Setup directory
mkdir /workspace/inputs/reference/fusion/per-feature
cd /workspace/inputs/reference/fusion/per-feature
#  Split fasta at each instance of a sequence header and write to new file
csplit -s -z ../Homo_sapiens.GRCh38.cdna.all.fa '/>/' '{*}'
#  If from chromosomes 6 or 17, rename files using the columns of the original ensemble header
#  (This step takes about 12 minutes. You can proceed with the next section in /workspace/inputs/data/fastq/chr6_and_chr17)
for f in xx*; do awk -F ":" 'NR==1 && $3=="6" || $3=="17"{print $2 "." $3 "." $4 "." $5}' $f | xargs -I{} mv $f {}.fa; done

#  Concatenate features from chromsomes 6 and 17 to a new reference fasta  
cd /workspace/inputs/reference/fusion
cat ./per-feature/GRCh38.17.*.fa ./per-feature/GRCh38.17.*.fa > chr617.fa
rm -rf per-feature
```

- To get one read pair each for normal and tumor, merge the chr6_and_chr17 only RNA-seq fastqs (2 min):
```bash
mkdir -p /workspace/inputs/data/fastq/chr6_and_chr17/fusion
cd /workspace/inputs/data/fastq/chr6_and_chr17/fusion
cat ../RNAseq_Norm/RNAseq_Norm_Lane1_R1.fastq.gz ../RNAseq_Norm/RNAseq_Norm_Lane2_R1.fastq.gz > RNAseq_Norm_R1.fastq.gz
cat ../RNAseq_Norm/RNAseq_Norm_Lane1_R2.fastq.gz ../RNAseq_Norm/RNAseq_Norm_Lane2_R2.fastq.gz > RNAseq_Norm_R2.fastq.gz
cat ../RNAseq_Tumor/RNAseq_Tumor_Lane1_R1.fastq.gz ../RNAseq_Tumor/RNAseq_Tumor_Lane2_R1.fastq.gz > RNAseq_Tumor_R1.fastq.gz
cat ../RNAseq_Tumor/RNAseq_Tumor_Lane1_R2.fastq.gz ../RNAseq_Tumor/RNAseq_Tumor_Lane2_R2.fastq.gz > RNAseq_Tumor_R2.fastq.gz
```

- Subsample fastqs to allow fusion alignment to run quickly (5 min):
```bash
# Use seqtk to take subsamples of the 10% of the fastq read pairs  
seqtk sample -s100 RNAseq_Norm_R1.fastq.gz 0.1 > subRNAseq_Norm_R1.fastq.gz
seqtk sample -s100 RNAseq_Norm_R2.fastq.gz 0.1 > subRNAseq_Norm_R2.fastq.gz
seqtk sample -s100 RNAseq_Tumor_R1.fastq.gz 0.1 > subRNAseq_Tumor_R1.fastq.gz
seqtk sample -s100 RNAseq_Tumor_R2.fastq.gz 0.1 > subRNAseq_Tumor_R2.fastq.gz
```

# Run Fusion Alignment 
- Create kallisto index:

```bash
mkdir -p /workspace/rnaseq/fusion
cd /workspace/rnaseq/fusion
kallisto index -i index.617.idx -k 31 --make-unique /workspace/inputs/reference/fusion/chr617.fa
```

- Quantify potential fusions (4 min):

```bash
kallisto quant -i index.617.idx --fusion -o kquant-norm617 /workspace/inputs/data/fastq/chr6_and_chr17/fusion/subRNAseq_Norm_R1.fastq.gz /workspace/inputs/data/fastq/chr6_and_chr17/fusion/subRNAseq_Norm_R2.fastq.gz
kallisto quant -i index.617.idx --fusion -o kquant-tumor617 /workspace/inputs/data/fastq/chr6_and_chr17/fusion/subRNAseq_Tumor_R1.fastq.gz /workspace/inputs/data/fastq/chr6_and_chr17/fusion/subRNAseq_Tumor_R2.fastq.gz
```

- Call fusions with pizzly:

```bash
# (Output is created in the directory where command is run)
cd /workspace/rnaseq/fusion
pizzly -k 31 --gtf /workspace/inputs/reference/fusion/chr617.gtf --cache index-norm617.cache.txt --align-score 2 --insert-size 400 --fasta /workspace/inputs/reference/fusion/chr617.fa --output norm-fuse617 kquant-norm617/fusion.txt
pizzly -k 31 --gtf /workspace/inputs/reference/fusion/chr617.gtf --cache index-tumor617.cache.txt --align-score 2 --insert-size 400 --fasta /workspace/inputs/reference/fusion/chr617.fa --output tumor-fuse617 kquant-tumor617/fusion.txt
```

- See next section to investigate the output of the pizzly fusion calling

# Possible Additional Analysis
- Indexing full transcriptome in kallisto is not recommended for a machine instance with course specifications, and will likely result in a memory fault

- Fusion calling for full RNA-seq data set is possible, but will take about _minutes. After following dirctions above through the indexing step:

```bash
# Quantify potential fusions
cd /workspace/rnaseq/fusion
kallisto quant -i index.617.idx --fusion -o kquant-norm617-full /workspace/inputs/data/fastq/chr6_and_chr17/fusion/RNAseq_Norm_R1.fastq.gz /workspace/inputs/data/fastq/chr6_and_chr17/fusion/RNAseq_Norm_R2.fastq.gz
kallisto quant -i index.617.idx --fusion -o kquant-tumor617-full /workspace/inputs/data/fastq/chr6_and_chr17/fusion/RNAseq_Tumor_R1.fastq.gz /workspace/inputs/data/fastq/chr6_and_chr17/fusion/RNAseq_Tumor_R2.fastq.gz

# Call fusions
pizzly -k 31 --gtf /workspace/inputs/reference/fusion/chr617.gtf --cache index-norm617.cache.txt --align-score 2 --insert-size 400 --fasta /workspace/inputs/reference/fusion/chr617.fa --output norm-fuse617 kquant-norm617/fusion.txt
pizzly -k 31 --gtf /workspace/inputs/reference/fusion/chr617.gtf --cache index-tumor617.cache.txt --align-score 2 --insert-size 400 --fasta /workspace/inputs/reference/fusion/chr617.fa --output tumor-fuse617 kquant-tumor617/fusion.txt
```
```


