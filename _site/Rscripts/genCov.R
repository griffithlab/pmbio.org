################################################################################
###################### Tutorial ################################################
################################################################################

# read the coverage data into R
covData <- read.delim("STAT1_mm9_coverage.tsv")

# rename the columns
colnames(covData) <- c("chromosome", "start", "end", "TAC245", "TAC265")

# create a function to split the data frame into lists of data frames
samples <- c("TAC245", "TAC265")
a <- function(x, y){
    col_names <- c("chromosome", "end", x)
    y <- y[,col_names]
    colnames(y) <- c("chromosome", "end", "cov")
    return(y)
}
covData <- lapply(samples, a, covData)

names(covData) <- samples

# install and load the BSgenome package and list available genomes
# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome")
library("BSgenome")
available.genomes()

# install and load the mm9 BSgenome from UCSC
# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Mmusculus.UCSC.mm9")
library("BSgenome.Mmusculus.UCSC.mm9")
genomeObject <- BSgenome.Mmusculus.UCSC.mm9

# Install and load a TxDb object for the mm9 genome
# source("https://bioconductor.org/biocLite.R")
# biocLite("TxDb.Mmusculus.UCSC.mm9.knownGene")
library("TxDb.Mmusculus.UCSC.mm9.knownGene")
TxDbObject <- TxDb.Mmusculus.UCSC.mm9.knownGene

# get the chromosome, and the start and end of our coverage data
chromosome <- as.character(unique(covData[[1]]$chromosome))
start <- as.numeric(min(covData[[1]]$end))
end <- as.numeric(max(covData[[1]]$end))

# define the genomic range
grObject <- GRanges(seqnames=c("chr1"), ranges=IRanges(start=start, end=end))

# construct an initial plot
genCov(x=covData, txdb=TxDbObject, gr=grObject, genome=genomeObject, cov_plotType="line")

# add a flank and create another plot
grObject <- GRanges(seqnames=c("chr1"), ranges=IRanges(start=start-500, end=end+500))
genCov(x=covData, txdb=TxDbObject, gr=grObject, genome=genomeObject, cov_plotType="line", gene_isoformSel="uc007aya.1")

# adjust compression of genomic features
genCov(x=covData, txdb=TxDbObject, gr=grObject, genome=genomeObject, cov_plotType="line", gene_isoformSel="uc007aya.1", base=NA, transform=NA)

################################################################################
############### Exercises ######################################################
################################################################################

genCov(x=covData, txdb=TxDbObject, gr=grObject, genome=genomeObject, cov_plotType="point", gene_isoformSel="uc007aya.1", label_bgFill="darkorchid4", label_txtFill="grey80", cov_colour="tomato2", gene_colour="tomato3")
