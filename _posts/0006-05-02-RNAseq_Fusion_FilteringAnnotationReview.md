---
feature_text: |
  ## Precision Medicine
title: RNAseq Fusion Filtering/Annotation/Review
categories:
    - Module-06-RNAseq
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0006-05-02
---

DEV NOTES
- consider moving package install to prebuilt
- can just wget to script https://raw.githubusercontent.com/MattBashton/grolar/master/grolar.R

---

# Introduction
Pizzly generates outputs in `.fasta` and `.json` formats. Some initial filtering is performed automatically in pizzly, for example removing alignments where the distance of the breakpoint to exon boundaries is 10 or more base pairs. These automatically filtered reads are included in the outputs with `unfiltered.` prefix. In this module we will perform additional annotation, filtering and visualization of the `.json` output.

# Annotation
JSON data are name/value pairs separated by a colon. Pairs are organized into objects within curly braces and arrays within square brackets. JSON data can reorganized into a delimited text file using many tools and programming languages. Below, we will use R and a modified script from the [grolar](https://github.com/MattBashton/grolar/blob/master/grolar.R) GitHub repository to annotate pizzly output and create a tabular, annotated output. 

```bash
# Open the R environment and set working directory 
R
setwd("/workspace/rnaseq/fusion/")

# Install packages if necessary
#  If asked about updating old packages (e.g. Old packages: 'MASS', 'devtools', 'foreign', 'ggplot2', 'mgcv', 'rlang',
  'survival' Update all/some/none? [a/s/n]:), select n
# install.packages("jsonlite")
# install.packages("dplyr")
# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("ensembldb")
# biocLite("EnsDb.Hsapiens.v86")

# Load required packages
library(jsonlite)
library(dplyr)
library(EnsDb.Hsapiens.v86)
edb = EnsDb.Hsapiens.v86

# Load data
suffix = "fusion.json"
JSON_files = list.files(path = "/workspace/rnaseq/fusion", pattern = paste0("*",suffix))
Ids = gsub(suffix, "", JSON_files)

# Load grolar.R script and run the GetFusionz_and_namez function to annotate
source("/usr/local/bin/mod-grolar.R)

# The function is too long to copy out line by line, but we can view it by calling it without variables
GetFusionz_and_namez()

# Run annotation script on each file in /workspace/rna/fusion suffixed with fusion.json
lapply(Ids, function(x) GetFusionz_and_namez(x, suffix=suffix)))
```