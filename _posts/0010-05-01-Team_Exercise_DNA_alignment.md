---
feature_text: |
  ## Precision Medicine Bioinformatics
  Introduction to bioinformatics for DNA and RNA sequence analysis
title: Team Exercise: DNA Alignment & IGV
categories:
    - Module-10-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0010-04-01
---

The goal of this DNA alignment lab is to familiarize students with the steps of DNA alignment as well as visualizing their data in IGV.
Here, we have prepared 5 individual datasets with hidden identities. Utilizing the knowledge learnt, students will attempt to uncover the organism (Human or Mouse) and type of data (WGS, WES or RNA) they were assigned.


### Obtain Fastq Data

First, you will need to download fastq data files for the upcoming analysis.
Before downloading the data, please create a directory specifically for this team exercise.
```bash
cd ~/workspace
mkdir dna_alignment_exercise
cd dna_alignment_exercise
```


According to your assigned team number, please download using the corresponding set of commands below.

```bash
# TEAM A
wget -c http://genomedata.org/seq-tec-workshop/read_data/dna_alignment_exercise/dataset_A/dataset_A_R1.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/dna_alignment_exercise/dataset_A/dataset_A_R2.fastq.gz

# TEAM B
wget -c http://genomedata.org/seq-tec-workshop/read_data/dna_alignment_exercise/dataset_B/dataset_B_R1.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/dna_alignment_exercise/dataset_B/dataset_B_R2.fastq.gz

# TEAM C
wget -c http://genomedata.org/seq-tec-workshop/read_data/dna_alignment_exercise/dataset_C/dataset_C_R1.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/dna_alignment_exercise/dataset_C/dataset_C_R2.fastq.gz

# TEAM D
wget -c http://genomedata.org/seq-tec-workshop/read_data/dna_alignment_exercise/dataset_D/dataset_D_R1.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/dna_alignment_exercise/dataset_D/dataset_D_R2.fastq.gz

# TEAM E
wget -c http://genomedata.org/seq-tec-workshop/read_data/dna_alignment_exercise/dataset_E/dataset_E_R1.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/dna_alignment_exercise/dataset_E/dataset_E_R2.fastq.gz
```
### Obtain Reference Files

In order to perform alignment, teams will also need to identify the correct reference file. The human and mouse reference files can be downloaded as following.

```bash
cd ~/workspace/dna_alignment_exercise
mkdir references human mouse

# Human reference files
cd ~/workspace/dna_alignment_exercise/references/human
wget -c http://genomedata.org/seq-tec-workshop/references/human/chr21/*

# Mouse reference files
cd ~/workspace/dna_alignment_exercise/references/mouse
wget -c http://genomedata.org/seq-tec-workshop/references/mouse/chr19/*

```

### Assignment

After obtaining the necessary fastq files as well as references, students now have all necessary materials to perform alignment and subsequent visualization analysis. As mentioned previously, teams now need to figure out which of the following choices corresponds to their given dataset.

1. Mouse WGS
2. Human WGS
3. Mouse WES
4. Human WES
5. Mouse RNA
6. Human RNA

Note: These datasets have been subsetted for optimized runtimes. Human data subsetted to `chr22:19,591,397-27,525,431` and mouse data subsetted to `chr19:26,436,477-35,060,791`.

When teams are confident with the result they have found, please send the instructors/TAs the ip address of the instance where the aligned bam and index files are located.

At the end of this exercise, teams will reveal their answer and present the pieces of supporting evidence (e.g. IGV) to the class.

###
