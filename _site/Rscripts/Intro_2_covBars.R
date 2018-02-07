################################################################################
############### Tutorial #######################################################
################################################################################

# load the ouput of samtools depth into R
seqData <- read.delim(url("http://genomedata.org/gen-viz-workshop/GenVisR/ALL1_CaptureDepth.tsv"), header=F)
colnames(seqData) <- c("chr", "position", "samp1", "samp2", "samp3", "samp4")
seqData <- seqData[,c(3:6)]

# Count the occurrences of each coverage value
# install.packages("plyr")
library(plyr)
seqCovList <- apply(seqData, 2, count)

# rename the columns for each dataframe
renameCol <- function(x){
    colnames(x) <- c("coverage", "freq")
    return(x)
}
seqCovList <- lapply(seqCovList, renameCol)

# create framework data frame with entries for the min to the max coverage
maximum <- max(unlist(lapply(seqCovList, function(x) max(x$coverage))))
minimum <- min(unlist(lapply(seqCovList, function(x) min(x$coverage))))
covFramework <- data.frame("coverage"=minimum:maximum)

# Merge the framework data frame with the coverage
# seqCovList <- lapply(seqCovList, function(x, y) merge(x, y, by="coverage", all=TRUE), covFramework)
seqCovList <- lapply(seqCovList, merge, covFramework, by="coverage", all=TRUE)

# merge all data frames together
seqCovDataframe <- Reduce(function(...) merge(..., by="coverage", all=T), seqCovList)

# set all NA values to 0
seqCovDataframe[is.na(seqCovDataframe)] <- 0

# set the rownames, remove the extra column, and convert to a matrix
rownames(seqCovDataframe) <- seqCovDataframe$coverage
seqCovDataframe$coverage <- NULL
seqCovMatrix <- as.matrix(seqCovDataframe)

# rename columns
colnames(seqCovMatrix) <- c("Skin_d42_I", "SL_d3072_I", "SB_d3072_A", "M_d3068_A")

# run covBars
#pdf(file="../assets/GenVisR/Coverage_Summary_v2.pdf", height=6, width=8)
covBars(seqCovMatrix)
#dev.off()

# ceiling pileups to 1200
column_sums <- colSums(seqCovMatrix[1200:nrow(seqCovMatrix),])
column_sums <- t(as.matrix(column_sums))
rownames(column_sums) <- 1200
seqCovMatrix2 <- seqCovMatrix[1:1199,]
seqCovMatrix2 <- rbind(seqCovMatrix2, column_sums)

# run covBars
#pdf(file="../assets/GenVisR/Coverage_Summary_v3.pdf", height=6, width=8)
covBars(seqCovMatrix2)
#dev.off()

# run covBars
rainbowColorRamp <- rainbow(1200)[1:1050]
#pdf(file="../assets/GenVisR/Coverage_Summary_v4.pdf", height=6, width=8)
covBars(seqCovMatrix2, colour=rainbowColorRamp)
#dev.off()

################################################################################
################# Exercise 1 ###################################################
################################################################################

# load ggplot2
library(ggplot2)

############### Make additional layers to pass to covBars ######################

# alter the y-axis facet labels
colnames(seqCovMatrix2) <- c("Skin_d42_I\n(skin normal)", "SL_d3072_I\n(sorted lymphs)", "SB_d3072_A\n(relapse 2 sorted blasts)", "M_d3068_A\n(relapse 2 bluk marrow)")

# remove the grey facet boxes
layer1 <- theme(strip.background=element_rect(fill="white"))

# position the facet text on the right side
layer2 <- facet_grid(sample ~ ., switch="y")
layer3 <- theme(strip.text.y=element_text(angle=180))

# add a plot title
layer4 <- ggtitle("Custom Capture Data")

# change the x-axis title
layer5 <- xlab("Cumulative Coverageof Region\nTargeted by and Successfully\nTiled by Nimblegen")

# change the legend
rainbowColorRamp <- rainbow(1200)[1:1050]
layer6 <- scale_fill_gradientn(name="Seq.\nDepth", colours=rainbowColorRamp, labels=c("0", "300", "600", "900", "1200+"))

# run covBars with all of these layers
covBars(seqCovMatrix2, plotLayer = list(layer1, layer2, layer3, layer4, layer5, layer6))


