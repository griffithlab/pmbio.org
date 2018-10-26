---
feature_text: |
  ## Precision Medicine
title: Somatic LOH Filtering/Annotation/Review
categories:
    - Module-05-Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-05-02
---

### Module objectives
- Find segments of chromosomes 6 and 17 that show LOH

### Run DNAcopy
To find regions of LOH, we will use R. In particular, the package that we will use to determine segments of LOH is the R package DNAcopy.

```bash
# Start R
R

# Load libraries
library(dplyr); library(tidyr); library(ggplot2); library(DNAcopy)

# Read table with tumor VAFs into R
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

# Read in table with germline VAFs for plotting
germline <- read.delim("heterozygous.snps.table", header = TRUE, col.names = c("CHROM", "POS", "GT", "AD", "DP", "VAF"))

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

# Exit R, no need to save workspace
q()
```
In your loh directory, you should now see two new files
1. chr6_segmentation.png
2. chr17_segmentation.png

The plots created will look similar to these:
{% include figure.html image="/assets/module_5/chr6_segmentation.png" position="left" width="1340" %}

{% include figure.html image="/assets/module_5/chr17_segmentation.png" position="left" width="1340" %}
Blue is our normal heterozygous VAFs, green is our tumor VAFs, and black is our segmentation results. Black lines close to 0 or 1 indicate segments of complete LOH, while lines between .4 to .6 indicate regions that are still heterozygous.
