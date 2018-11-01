---
feature_text: |
  ## Precision Medicine
title: RNAseq Differential Expression
categories:
    - Module-06-RNAseq
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0006-03-01
---


#### Running Ballgown for differential expression
We will use `Ballgown` to perform differential expression analysis on the HISAT/StringTie abundance estimates for each of our samples "replicates"?

```bash
# change working directory
cd /workspace/rnaseq/

# start R
R

# In R, run the following commands
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

# create the phenotype data for each sample
#pheno_data = read.csv("RNA_data.csv")
pheno_data <- data.frame(ids=c("RNAseq_Norm", "RNAseq_Norm_Lane1", "RNAseq_Norm_Lane2", "RNAseq_Tumor", "RNAseq_Tumor_Lane1", "RNAseq_Tumor_Lane2"), type=c("normal", "normal", "normal", "tumor", "tumor", "tumor"))

# Load ballgown data structures for each sample
bg = ballgown(dataDir = "ballgown", samplePattern = "RNAseq", pData=pheno_data)

# Filter low-abundance genes
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Identify signficant differently expressed Transcripts
results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")

# Identify significant differently expressed Genes
results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")

# Add transcript/gene names and transcript/gene IDs to the results_transcripts data frame
results_transcripts = data.frame(transcriptNames=ballgown::transcriptNames(bg_filt),transcriptIDs=ballgown::transcriptIDs(bg_filt),geneNames=ballgown::geneNames(bg_filt),geneIDs=ballgown::geneIDs(bg_filt),results_transcripts)

# Add common gene names on
tmp <- unique(results_transcripts[,c("geneNames", "geneIDs")])
results_genes <- merge(results_genes, tmp, by.x=c("id"), by.y=c("geneIDs"), all.x=TRUE)
results_genes = data.frame(geneNames=ballgown::geneNames(bg_filt),geneIDs=ballgown::geneIDs(bg_filt), results_genes)

# Sort from the smallest P value to largest
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)

# Output as CSV
write.csv(results_transcripts,"RNAseq_transcript_results.csv",row.names=FALSE)
write.csv(results_genes,"RNAseq_genes_results.csv",row.names=FALSE)

# Output as TSV
write.table(results_transcripts,"RNAseq_transcript_results.tsv",sep="\t", quote=FALSE, row.names=FALSE)
write.table(results_genes,"RNAseq_gene_results.tsv",sep="\t", quote=FALSE, row.names=FALSE)

# Identify genes with p value < 0.05
subset(results_transcripts,results_transcripts$pval<0.05)
subset(results_genes,results_genes$pval<0.05)

# quit R
q()
```
