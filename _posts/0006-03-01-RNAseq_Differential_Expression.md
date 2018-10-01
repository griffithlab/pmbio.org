---
feature_text: |
  ## Precision Medicine
title: RNAseq Differential Expression
categories:
    - Module-06-RNAseq
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0006-03-01
---

#### **Differential Expression**

```bash
cd /data/RNA_seq
# Make CSV/TSV for phenotype data
printf "\"ids\",\"type\"\n\"HCC1395BL_RNA_H3MYFBBXX_4_CTTGTA\",\"normal\"\n\"HCC1395BL_RNA_H3MYFBBXX_5_CTTGTA\",\"normal\"\n\"HCC1395_RNA_H3MYFBBXX_4_GCCAAT\",\"tumor\"\n\"HCC1395_RNA_H3MYFBBXX_5_GCCAAT\",\"tumor\"\n" > HCC1395.csv
#Run R
~/bin/R
# In R, run the following commands
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
# Load the phenotype data for each sample
pheno_data = read.csv("HCC1395.csv")
# Load ballgown data structures for each sample
bg = ballgown(dataDir = "ballgown", samplePattern = "HCC1395", pData=pheno_data)
# Filter low-abundance genes
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)
# Identify signficant differently expressed Transcripts
results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
# Identify significant differently expressed Genes
results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
# Add gene names and gene IDs to the retuls_transcripts data frame
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filt),geneIDs=ballgown::geneIDs(bg_filt), results_transcripts)
# Sort from the smallest P value to largest
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)
# Output as CSV
write.csv(results_transcripts,"HCC1395_transcript_results.csv",row.names=FALSE)
write.csv(results_genes,"HCC1395_genes_results.csv",row.names=FALSE)
# Output as TSV
write.table(results_transcripts,"HCC1395_transcript_results.tsv",sep="\t")
write.table(results_genes,"HCC1395_gene_results.tsv",sep="\t")
# Identify genes with p value < 0.05
subset(results_transcripts,results_transcripts$pval<0.05)
subset(results_genes,results_genes$pval<0.05)
```
