---
feature_text: |
  ## Precision Medicine Bioinformatics
  Introduction to bioinformatics for DNA and RNA sequence analysis
title: DNA alignment lab
categories:
    - Module-10-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0010-04-01
---

This lab will introduce students to a basic, and very common commandline bioinformatics task. We will obtain raw DNA sequence data, a reference sequence, perform any necessary indexing, align the raw data to reference, mark duplicate reads, and then prepare the alignment results for visualization in a genome viewer. 

### Prepare directories for analysis

First, let's create a directory in your workspace for all your results (e.g. dna_alignment_lab). Within that directory we will also make three subdirectories to organize our analysis: 

- /workspace/dna_alignment_lab/alignment_results
- /workspace/dna_alignment_lab/fastq_files
- /workspace/dna_alignment_lab/reference_sequences

```bash
cd /workspace
mkdir dna_alignment_lab
cd dna_alignment_lab
mkdir alignment_results
mkdir fastq_files
mkdir reference_sequences
```

### Obtain sequence data for alignment (fastq files) 

We have provided a sample dataset for this exercise. This data represents paired (2 x 100bp) DNA exome sequence data in fastq format. DNA was isolated from a breast cancer cell line (HCC1395). We have provided sequence reads that have been limited to chromosome 21. 

{% include question.html question="How could raw sequence data be limited to chr21?" answer='We previously aligned the complete dataset, extracted only chr21 alignments, then reverted back to the raw (unaligned) reads.' %}

```bash
#download seq data 
cd /workspace/dna_alignment_lab/fastq_files 

wget http://genomedata.org/seq-tec-workshop/read_data/dna_alignment_exercise/dataset_lab/HCC1395_Exome_chr21_R1.fastq.gz
wget http://genomedata.org/seq-tec-workshop/read_data/dna_alignment_exercise/dataset_lab/HCC1395_Exome_chr21_R2.fastq.gz
```
### Obtain a reference sequence

For convenience, we also provide a reference genome file, limited to chromosome 21.

```bash 
#download ref data 
cd /workspace/dna_alignment_lab/reference_sequences 

wget http://genomedata.org/seq-tec-workshop/references/human/chr21/chr21_references.fa
```

### Index reference file with bwa 

Now that we have our data, we need to create the files necessary to run an alignment. The first thing we need to do is index our reference sequence. Indexing your reference sequence allows the aligner to rapidly narrow down the potential origin of the query sequence within the genome, saving both time and memory.

```bash
cd /workspace/dna_alignment_lab/reference_sequences 

bwa index chr21_references.fa
```

### Run bwa mem to create an alignment 

OK, now that we have an indexed reference sequence we are ready to create an alignment.

```bash
cd /workspace/dna_alignment_lab/alignment_results

bwa mem -t 8 -o /workspace/dna_alignment_lab/alignment_results/HCC1395_Exome_chr21.sam /workspace/dna_alignment_lab/reference_sequences/chr21_references.fa /workspace/dna_alignment_lab/fastq_files/HCC1395_Exome_chr21_R1.fastq.gz /workspace/dna_alignment_lab/fastq_files/HCC1395_Exome_chr21_R2.fastq.gz
```

### Convert sam to bam format

We've got an alignment! But we have several post processing steps to complete. The first step is to convert your large .sam file to compressed .bam file.

```bash
cd /workspace/dna_alignment_lab/alignment_results

samtools view -@ 8 -h -b -o HCC1395_Exome_chr21.bam HCC1395_Exome_chr21.sam
```

### Query name sort bam files

Next we need to sort the reads by their read names. This is expected for the next step (duplicate marking).

```bash
cd /workspace/dna_alignment_lab/alignment_results

java -Xmx60g -jar /home/ubuntu/bin/picard.jar SortSam I=HCC1395_Exome_chr21.bam O=HCC1395_Exome_chr21_namesorted_picard.bam SO=queryname
```

### Perform duplicate marking

Next we need to mark the duplicate reads within our data. Duplicate reads are typically defined as reads with identical start and stop alignment positions. These reads are likely to be exact copies of a single DNA molecule. These reads introduce bias in your variant calling. If you did not mark duplicates, you risk over-representing areas that are preferentially amplified during PCR. 

```bash
cd /workspace/dna_alignment_lab/alignment_results

java -Xmx60g -jar /home/ubuntu/bin/picard.jar MarkDuplicates I=HCC1395_Exome_chr21_namesorted_picard.bam  O=HCC1395_Exome_chr21_namesorted_picard_mrkdup.bam ASSUME_SORT_ORDER=queryname METRICS_FILE=HCC1395_Exome_chr21_mrk_dup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT
# This command will also print out a txt file that gives you some metrics about the number of duplicates identified 
```

### Position sort bam files
Next we need to position sort the bam file. Indexing requires a position-sorted bam. 

```bash
cd /workspace/dna_alignment_lab/alignment_results

java -Xmx60g -jar /home/ubuntu/bin/picard.jar SortSam I=HCC1395_Exome_chr21_namesorted_picard_mrkdup.bam O=HCC1395_Exome_chr21_pos_sorted_mrkdup_picard.bam SO=coordinate
```

### Index bam file 

In order to efficiently load and search a bam file, downstream applications typically require an index. This is very similar to the index you created of your reference file.

```bash
cd /workspace/dna_alignment_lab/alignment_results

java -Xmx60g -jar /home/ubuntu/bin/picard.jar BuildBamIndex I=HCC1395_Exome_chr21_pos_sorted_mrkdup_picard.bam
```

### Generate a flagstat report for bam file

There are many different ways to assess the quality of an alignment. We are going to use one method to get a fast summary of some common alignment statistics.

```bash 
cd /workspace/dna_alignment_lab/alignment_results

samtools flagstat HCC1395_Exome_chr21_pos_sorted_mrkdup_picard.bam > HCC1395_Exome_chr21_pos_sorted_mrkdup_picard_flagstat.flagstat
```

### Clean up un-needed sam/bam files

Keep final sorted, duplicated marked, bam/bai files and mrkdup.txt files. Delete everything else.

```bash
cd /workspace/dna_alignment_lab/alignment_results

mkdir final
mv HCC1395_Exome_chr21_pos_sorted_mrkdup_picard.* final/
mv *.txt final/
mv *.flagstat

rm *.sam
rm *.bam
```

### View your alignment  

Finally, we will visualize the alignments in IGV, a popular genome viewer. We will provide instructions for this. This is the URL you will use to view your bam. Note, you will need to replace `student.i.p.address` with the IP address of your personal AWS instance.

```bash
http://student.i.p.address/workspace/dna_alignment_lab/alignment_results/final/HCC1395_Exome_chr21_pos_sorted_mrkdup_picard.bam
```
