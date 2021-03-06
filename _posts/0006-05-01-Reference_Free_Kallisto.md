---
feature_text: |
  ## Precision Medicine Bioinformatics
  Introduction to bioinformatics for DNA and RNA sequence analysis
title: Reference Free Expression Analysis with Kallisto
categories:
    - Module-06-RNAseq
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0006-05-01
---

### Transcript annotation and sequences
Remember that in previous sections we have been using reference **genome** fasta sequences for the reference for alignment and subsequent steps. However, Kallisto works directly on target cDNA/transcript sequences. Remember also that in the [Annotation section of the Inputs module](/module-02-inputs/0002/03/01/Annotation/) we have previously obtained transcript annotations for genes on our subset of chromosomes (i.e. chr6 and chr17). The transcript models were downloaded from Ensembl in GTF format. This GTF contains a description of the coordinates of exons that make up each transcript but it does not contain the transcript sequences themselves. We used the `gtf_to_fasta` tool and those exon coordinates to construct the transcript sequences needed by Kallisto. There are many places we could obtain such transcript sequences. For example, we could have download them directly in Fasta format from the Ensembl FTP site (or from UCSC or NCBI).

Also remember that in the [Indexing section of the Inputs module](/module-02-inputs/0002/04/01/Indexing/) we created a Kallisto index using the sequences in the transcripts Fasta file using the default k-mer size (which is 31).

Make sure you still have the transcript GTF, Fasta, and Kallisto indexes:
```bash
cd /workspace/rnaseq/
mkdir kallisto
cd kallisto

# first check that the GTF and Fasta file are present
head /workspace/inputs/references/transcriptome/ref_transcriptome.gtf
head /workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_clean.fa

# now check for the kallisto index that we previously created
ls /workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_kallisto_index

# create a list of all transcript IDs for later use:
cd /workspace/rnaseq/kallisto/
cat /workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_clean.fa | grep ">" | perl -ne '$_ =~ s/\>//; print $_' | sort | uniq > transcript_id_list.txt

```

### Generate abundance estimates for all samples using Kallisto
As we did with StringTie we will generate transcript abundances for each of our demonstration samples using Kallisto. Here we are treating the two lanes for each sample as if they were independent samples.

```bash
cd /workspace/rnaseq/kallisto/
mkdir quants
cd quants

kallisto quant --index=/workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_kallisto_index --output-dir=RNAseq_Norm_Lane1 --threads=8 --plaintext /workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane1_R1.fastq.gz /workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane1_R2.fastq.gz
kallisto quant --index=/workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_kallisto_index --output-dir=RNAseq_Norm_Lane2 --threads=8 --plaintext /workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane2_R1.fastq.gz /workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane2_R2.fastq.gz

kallisto quant --index=/workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_kallisto_index --output-dir=RNAseq_Tumor_Lane1 --threads=8 --plaintext /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane1_R1.fastq.gz /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane1_R2.fastq.gz
kallisto quant --index=/workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_kallisto_index --output-dir=RNAseq_Tumor_Lane2 --threads=8 --plaintext /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane2_R1.fastq.gz /workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane2_R2.fastq.gz

```

Create a single TSV file that has the TPM abundance estimates for all six samples.
```bash
cd /workspace/rnaseq/kallisto/quants/
paste */abundance.tsv | cut -f 1,2,5,10,15,20 > transcript_tpms_all_samples.tsv
ls -1 */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > header.tsv
cat header.tsv transcript_tpms_all_samples.tsv | grep -v "tpm" > transcript_tpms_all_samples.tsv2
mv transcript_tpms_all_samples.tsv2 transcript_tpms_all_samples.tsv
rm -f header.tsv
head transcript_tpms_all_samples.tsv

```

### Compare transcript and gene abundance estimates from Kallisto to isoform abundance estimates from StringTie
How similar are the results we obtained from each approach?

We can compare the expression value for each Ensembl transcript.

To do this comparison, we need to gather the expression estimates for each of our replicates from each approach. The Kallisto transcript results were neatly organized into a single file above. For Kallisto gene expression estimates, we will simply sum the TPM values for transcripts of the same gene. The following R code will pull together the various expression matrix files we created in previous steps and create some visualizations to compare them (for both transcript and gene estimates).

First create the gene version of the Kallisto TPM matrix
```bash
cd /workspace/rnaseq/kallisto/quants
wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/kallisto_gene_matrix.pl
chmod +x kallisto_gene_matrix.pl
./kallisto_gene_matrix.pl --gtf_file=/workspace/inputs/references/transcriptome/ref_transcriptome.gtf --kallisto_transcript_matrix_in=transcript_tpms_all_samples.tsv --kallisto_transcript_matrix_out=gene_tpms_all_samples.tsv

```

Now load files and summarize results from each approach in R
```R
# example code here:
# https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_comparisons.R

# start R
R

# set proper working directory and load libraries
setwd("/workspace/rnaseq/kallisto/quants")
library(ggplot2)
library(reshape2)

# read in data
kallisto_transcript_tpm <- read.delim("transcript_tpms_all_samples.tsv")
stringtie_transcript_tpm <- read.delim("~/workspace/rnaseq/ref-only-expression/transcript_tpm_all_samples.tsv")

# minor reformatting
kallisto_transcript_tpm <- kallisto_transcript_tpm[,-2]
kallisto_transcript_tpm <- melt(kallisto_transcript_tpm, id.vars=c("target_id"))
stringtie_transcript_tpm <- melt(stringtie_transcript_tpm, id.vars=c("Transcript_ID"))

# merge the data
kallisto_stringtie_tpm <- merge(kallisto_transcript_tpm, stringtie_transcript_tpm, by.x=c("target_id", "variable"), by.y=c("Transcript_ID", "variable"), suffixes=c(".kallisto", ".stringtie"))

########### plot the result ######################

pdf(file="transcript_stringtie_v_kallisto.pdf", height=8, width=8)
ggplot(data=kallisto_stringtie_tpm, aes(x=log2(value.kallisto + .1), y=log2(value.stringtie + .1))) + geom_point(alpha=.1) + facet_wrap(~variable) + theme_bw()
dev.off()

pdf(file="transcript_stringtie_v_kallisto_density_hex.pdf", height=8, width=8)
ggplot(data=kallisto_stringtie_tpm, aes(x=log2(value.kallisto + .1), y=log2(value.stringtie + .1))) + geom_point(alpha=.1) + facet_wrap(~variable) + theme_bw() + stat_binhex() + scale_fill_continuous(trans="sqrt")
dev.off()

```

### Exercise: Create a custom Kallisto index for a selected subset of gene. Calculate abundances for this set only
For example, you could create a custom GTF of rRNA genes, or cancer genes. Then create a custom transcriptome fasta file for this set. Then create a kallisto index.  Then get kallisto quants using this custom index.

### Exercise: Use Sleuth to perform differential expression analysis with Kallisto output
You may need to refer to the Sleuth documentation to determine if kallisto needs to be rerun in some specific way...
