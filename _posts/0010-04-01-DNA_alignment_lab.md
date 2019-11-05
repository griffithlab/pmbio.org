---
feature_text: |
  ## Precision Medicine Bioinformatics Introduction to bioinformatics for DNA and RNA sequence analysis
title: DNA alignment lab
categories:
    - Module-10-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0010-04-01
---
First you want to create a directory in your workspace for all your results of this exercise  (e.g. alignment_exercise).  Within that directory make three subdirectories: 

- /workspace/username/align_lab/alignment_results
- /workspace/username/align_lab/fastq_files 
-/workspace/username/align_lab/reference_sequences 

### Obtain sequence data for alignment (fastq files) 
```bash
#download seq data 
cd /workspace/username/align_lab/fastq_files 

wget  http://genomedata.org/seq-tec-workshop/read_data/dna_alignment_exercise/dataset_lab/2891351068_1.fastq.gz
wget  http://genomedata.org/seq-tec-workshop/read_data/dna_alignment_exercise/dataset_lab/2891351068_2.fastq.gz
````
### Obtain a reference for seq data
````bash 
#download ref data 
cd /workspace/username/align_lab/reference_sequences 

wget -r -l1 -np -nd -A dict,fa,fai http://genomedata.org/seq-tec-workshop/references/human/chr21
````
Now that we have our data we need to create the files necessary to create an alignment 
The first thing we need to do is index our reference sequence. Indexing your reference sequence allows the alignment program to narrow down the potential origin of a query sequence within the genome, saving both time and memory.

### Index reference file with bwa 
````bash
cd cd /workspace/username/align_lab/reference_sequences 

bwa index chr21_references.fa
````
OK, now that we have an indexed reference sequence we are ready to create an alignment 

### Run bwa mem to create an alignment 
````bash
cd /workspace/username/align_lab/alignment_results

bwa mem -t 8 -o /workspace/username/align_lab/alignment_results/2891351068.sam /workspace/username/align_lab/reference_sequences/chr21_references.fa /workspace/username/align_lab/fastq_files/2891351068_1.fastq.gz /workspace/username/align_lab/fastq_files/2891351068_2.fastq.gz
````
We've got an alignment, but we have several post processing steps to complete. 
The first step is to convert your large .sam file, to compressed binary.bam file

### Convert sam to bam format
```bash
cd /workspace/username/align_lab/alignment_results

samtools view -@ 8 -h -b -o 2891351068.bam 2891351068.sam
````
Next we need to sort the reads by their read names 

### Query name sort bam files
```bash
cd /workspace/username/align_lab/alignment_results

java -Xmx60g -jar /home/ubuntu/bin/picard.jar SortSam I=2891351068.bam O=2891351068_namesorted_picard.bam SO=queryname
````
Next we need to mark the duplicate reads within our data. Duplicate reads are typically defined as reads with identical start and stop alignment positions. These reads are likely to be exact copies of a single DNA molecule. These reads introduce bias in your variant calling. If you did not mark duplicates, you risk over-representing areas that preferentially amplified during PCR. 
```bash
cd /workspace/username/align_lab/alignment_results

java -Xmx60g -jar /home/ubuntu/bin/picard.jar MarkDuplicates I=2891351068_namesorted_picard.bam  O=2891351068_namesorted_picard_mrkdup.bam ASSUME_SORT_ORDER=queryname METRICS_FILE=2891351068_mrk_dup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT
## This  command will also print out a txt file that gives you some metrics about the number of duplicates identified 
````
Next we need to position sort the bam file. For indexing and other possible subsequent steps a position-sorted bam is required. 

```bash
cd /workspace/username/align_lab/alignment_results

java -Xmx60g -jar /home/ubuntu/bin/picard.jar SortSam I=2891351068_namesorted_picard_mrkdup.bam O=2891351068_pos_sorted_mrkdup_picard.bam SO=coordinate
````
In order to efficiently load and search a bam file, downstream applications typically require an index.  This is very similar to the index you created of your reference file

###index bam file 
java -Xmx60g -jar /home/ubuntu/bin/picard.jar BuildBamIndex I=2891351068_pos_sorted_mrkdup_picard.bam

View your alignment!!! 
Go to IGV - we will provide instructions for this.  This is the URL you will use to view your bam
````bash
http://34.239.1.158/workspace/username/align_lab/alignment_results/2891351068_pos_sorted_mrkdup_picard.bam
````
