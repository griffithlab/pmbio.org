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
Pizzly generates outputs in `.fasta` and `.json` formats. Some initial filtering is performed automatically in pizzly, for example removing alignments where the distance of the breakpoint to exon boundaries is 10 or more base pairs. These automatically filtered reads are included in the outputs with `unfiltered.` suffix. In this module we will perform additional annotation, filtering and visualization of the `.json` output.

# Annotation
JSON data are name/value pairs separated by a colon. Pairs are organized into objects within curly braces and arrays within square brackets. JSON data can be reorganized into a delimited text file using many tools and programming languages. Below, we will use R and a modified script from the [grolar](https://github.com/MattBashton/grolar/blob/master/grolar.R) GitHub repository to annotate pizzly output and create a tabular, annotated output. 

```bash
# Get R scripts for later use
cd /workspace/rnaseq/fusion
wget https://raw.githubusercontent.com/griffithlab/pmbio.org/master/assets/scripts/mod.grolar.R
wget https://raw.githubusercontent.com/griffithlab/pmbio.org/master/assets/scripts/importPizzly.R
```

```R
# Open the R environment and set working directory 
R
setwd("/workspace/rnaseq/fusion/")

# Install packages if necessary
#  If asked about updating old packages (e.g. Old packages: 'MASS', 'devtools'... Update all/some/none? [a/s/n]:), select n
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
suffix = "617.json"
JSON_files = list.files(path = "/workspace/rnaseq/fusion", pattern = paste0("*",suffix))
Ids = gsub(suffix, "", JSON_files)

# Load grolar.R script and run the GetFusionz_and_namez function to annotate
source("./mod.grolar.R")

# The function is too long to copy out line by line, but we can view it by calling it without variables
GetFusionz_and_namez

# Run annotation script on each file in /workspace/rna/fusion suffixed with fusion.json
lapply(Ids, function(x) GetFusionz_and_namez(x, suffix=suffix))
```

# Filtering

- The GetFusionz_and_namez script wrote two new files into /workspace/rnaseq/fusion: ```norm-fuse_fusions_filt_sorted.txt``` and ```tumor-fuse_fusions_filt_sorted.txt```. Let's read those back into R and filter further:

```R
normal=read.table("./norm-fuse_fusions_filt_sorted.txt", header=T)
tumor=read.table("./tumor-fuse_fusions_filt_sorted.txt", header=T)
names(normal)
```

- In addition to values present in the orignial ```.json`` files, each fusion now has a unique identifier, sequence positions, and distance values for genes from the same chromosome: 

```R
> names(normal)
 [1] "ID"                   "paircount"            "splitcount"
 [4] "geneA.id"             "geneA.name"           "geneA.seq_name"
 [7] "geneA.gene_seq_start" "geneA.gene_seq_end"   "geneA.seq_strand"
[10] "geneB.id"             "geneB.name"           "geneB.seq_name"
[13] "geneB.gene_seq_start" "geneB.gene_seq_end"   "geneB.seq_strand"
[16] "same_chr"             "gene_distance"
>
```

- Common filtering tasks for fusion output include removing fusions from the tumor sample which are present in the normal, and removing fusions for which the is little support by pair and split read counts: 

```R
normal$sample="normal"
tumor$sample="tumor"
allfusions=rbind(normal, tumor)
normal$genepr=paste0(normal$geneA.name,".",normal$geneB.name)
tumor$genepr=paste0(tumor$geneA.name,".",tumor$geneB.name)
uniqueTumor=subset(tumor, !(tumor$genepr %in% normal$genepr))
nrow(uniqueTumor)==nrow(tumor)
[1] TRUE
# All fusions from the tumor sample are unique to the tumor (and therefore all fusions in the normal sample are also unique from the tumor). 
```

- Split read and paired read counts are generally low, consistent with our heavily downsampled fastqs. To see the highest counts for normal fusions:

```R
highNormal=subset(normal, normal$paircount + normal$splitcount >1)
nrow(normal)
[1] 21
nrow(highNormal)
[1] 3
```

# Visualization: 
```
# Install and load chimeraviz 
source("https://bioconductor.org/biocLite.R")
biocLite("chimeraviz")
library(chimeraviz)

# Use the pizzly importer script to import fusion data
source("./importPizzly.R")
fusions = importPizzly(allfuion,"hg38")
```