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
germline <- read.delim("heterozygous.snps.table", header = TRUE, col.names = c("CHROM", "POS", "GT", "AD", "DP", "VAF"))
tumor <- read.delim("tumor.VAFs.table", header = TRUE, col.names = c("CHROM", "POS", "TUMOR_DP", "AD", "VAF"))

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

# Split LOH_segments file by chromosome for plotting
chr6_LOH_segments <- subset.data.frame(LOH_segments, CHROM == "chr6"); chr17_LOH_segments <- subset.data.frame(LOH_segments, CHROM == "chr17")

# Plot VAFs, LOH segments for chr6
png("chr6_segmentation.png", width = 1340, height = 300)
chr6_segment_plot <- ggplot() + geom_point(data = germline[germline$CHROM == "chr6", ], aes(POS,VAF), color="blue") + geom_point(data = tumor[tumor$CHROM == "chr6", ], aes(POS,VAF), color="green") + geom_segment(data = chr6_LOH_segments, aes(x = chr6_LOH_segments$START,y = chr6_LOH_segments$TOP,xend = chr6_LOH_segments$END,yend = chr6_LOH_segments$TOP), size = 1.5) + geom_segment(data = chr6_LOH_segments, aes(x = chr6_LOH_segments$START,y = chr6_LOH_segments$BOTTOM,xend = chr6_LOH_segments$END,yend = chr6_LOH_segments$BOTTOM), size = 1.5) + xlab("Chr6 Position") + ylab("VAF")
plot(chr6_segment_plot)
dev.off()

# Plot VAFs, LOH segments for chr17
png("chr17_segmentation.png", width = 1340, height = 300)
chr17_segment_plot <- ggplot() + geom_point(data = germline[germline$CHROM == "chr17", ], aes(POS,VAF), color="blue") + geom_point(data = tumor[tumor$CHROM == "chr17", ], aes(POS,VAF), color="green") + geom_segment(data = chr17_LOH_segments, aes(x = chr17_LOH_segments$START,y = chr17_LOH_segments$TOP,xend = chr17_LOH_segments$END,yend = chr17_LOH_segments$TOP), size = 1.5) + geom_segment(data = chr17_LOH_segments, aes(x = chr17_LOH_segments$START,y = chr17_LOH_segments$BOTTOM,xend = chr17_LOH_segments$END,yend = chr17_LOH_segments$BOTTOM), size = 1.5) + xlab("Chr17 Position") + ylab("VAF")
plot(chr17_segment_plot)
dev.off()
```
