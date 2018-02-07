################################################################################
######################### liftOver Tutorial ####################################

# install and load rtracklayer
source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer")
library("rtracklayer")

# specify coordinates to liftover
grObject <- GRanges(seqnames=c("chr1"), ranges=IRanges(start=226061851, end=226071523))

# import the chain file
chainObject <- import.chain("~/Desktop/hg38ToCanFam3.over.chain")

# run liftOver
results <- as.data.frame(liftOver(grObject, chainObject))