---
feature_text: |
  ## Precision Medicine
title: RNAseq Fusion Filtering/Annotation/Review
categories:
    - Module-06-RNAseq
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0006-05-02
---

# Introduction
Pizzly generates outputs in `.fasta` and `.json` formats. Some initial filtering is performed automatically in pizzly, for example removing alignments where the distance of the breakpoint to exon boundaries is 10 or more base pairs. These automatically filtered reads are included in the outputs with `unfiltered.` prefix.

# Setup
JSON data are name/value pairs separated by a colon. Pairs are organized into objects withing curly braces adn arrays within square brackets. JSON data can reorganized into a delimited text file using many tools and programming languages. Below, we will use the R programming langauge and code from the [grolar](https://github.com/MattBashton/grolar/blob/master/grolar.R) GitHub repository to annotate pizzly output and create a tabular output. 

# Filtering Example
```bash
# Open the R environment 
R

# Install and load reuqired packages
install.packages("jsonlite")
library(jsonlite)
install.packages("dplyr")
library(dplyr)
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("ensembldb")
biocLite("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
edb = EnsDb.Hsapiens.v86

# Load data
setwd("/data/RNA_seq/fusion/")
suffix = "fusion.json"
JSON_files = list.files(path = "/data/RNA_seq/fusion", pattern = paste0("*",suffix))
Ids <- gsub(suffix, "", JSON_files)

 GetFusionz_and_namez
function(sample, suffix) {
  JSON_file <- paste0(sample, suffix)
  JSON <- fromJSON(JSON_file, flatten = TRUE)
  JSON_level1 <- JSON$genes
  output <- JSON_level1
  output <- output[,c(-3,-4)]
  output$geneA.id <- gsub(".\\d+$", "", output$geneA.id, perl = TRUE)
  output$geneB.id <- gsub(".\\d+$", "", output$geneB.id, perl = TRUE)
  output <- cbind(rownames(output), output)
  tmp <- colnames(output)[-1]
  colnames(output) <- c("ID", tmp)
  geneAs <- output$geneA.id
  geneAs_info <- genes(edb,
        columns = c("gene_id","seq_name","gene_seq_start","gene_seq_end","seq_strand"),
        filter = GeneIdFilter(geneAs),
        return.type = "DataFrame")
  geneBs <- output$geneB.id
  geneBs_info <- genes(edb,
                       columns = c("gene_id","seq_name","gene_seq_start","gene_seq_end","seq_strand"),
                       filter = GeneIdFilter(geneBs),
                       return.type = "DataFrame")
  colnames(geneAs_info) <- paste0("geneA.", colnames(geneAs_info))
  colnames(geneBs_info) <- paste0("geneB.", colnames(geneBs_info))
  tmp1 <- merge(output, geneAs_info, by.x = "geneA.id", by.y = "geneA.gene_id", all.x = TRUE)
  tmp2 <- merge(output, geneBs_info, by.x = "geneB.id", by.y = "geneB.gene_id", all.x = TRUE)
  output <- merge(tmp1, tmp2[,c(2,8:11)], by = "ID")[,c(1,3,4,2,5,8,9,10,11,6,7,12,13,14,15)]
  super_identical <- Vectorize(identical, c("x", "y"))
  output <- as.data.frame(output)
  output <- mutate(output, same_chr = super_identical(output$geneA.seq_name,output$geneB.seq_name))
  identical_idx <- which(output$same_chr == TRUE)
  idx <- order(output$splitcount, output$paircount, decreasing = TRUE)
  output <- output[idx,]
  write.table(output, file = paste0(sample, "_fusions_filt_sorted.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
```