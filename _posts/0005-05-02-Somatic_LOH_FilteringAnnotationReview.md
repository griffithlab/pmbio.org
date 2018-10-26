---
feature_text: |
  ## Precision Medicine
title: Somatic LOH Filtering/Annotation/Review
categories:
    - Module-05-Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-05-02
---

### Objectives
-Determine segments of LOH

To find regions of LOH, will use R, DNAcopy package

```bash
# Start R
R

# Set working directory
setwd("")

# Load libraries
library(dplyr); library(tidyr); library(ggplot2); library(DNAcopy)

# Read tables with germline VAFs, tumor VAFs into R
germline <- read.delim("heterozygous.snps.table", header = TRUE, col.names("CHROM", "POS", "GT", "AD", "DP", "VAF"))
tumor <- read.delim("tumor_VAFs.table", header = TRUE, col.names("CHROM", "POS", "TUMOR_DP", "AD", "VAF"))

# Make sure tumor VAFs are numeric for calculations
tumor$VAF <- as.numeric(tumor$VAF)

# Create column for abs(.5 - tumor VAF)
tumor$ABS <- abs(.5 - tumor$VAF)

# Create CNA.object for DNAcopy segmentation analysis
LOH_CNA.object <- CNA(genomdat = tumor$ABS, chrom = tumor$CHROM, maploc = tumor$POS, data.type = 'binary')
LOH.segmentation <- segment(LOH_CNA.object)
LOH.segmentation.output <- LOH.segmentation$output

# Extract columns needed to plot segmentation results, change column names
LOH_segments <- select(LOH.segmentation.output, chrom, loc.start, loc.end, seg.mean)
colnames(LOH_segments) <- c("CHROM", "START", "END", "MEAN")

# Create columns for .5 + segment means, .5 - segment means
LOH_segments$TOP <- .5 + LOH_segments$MEAN; LOH_segments$BOTTOM <- .5 - LOH_segments$MEAN

# Plot VAFs, LOH segments for chr6

# Plot VAFs, LOH segments for chr17
```
