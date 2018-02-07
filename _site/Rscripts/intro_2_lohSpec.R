################################################################################
############################ Tutorial ##########################################
################################################################################
# read in the varscan data 
lohData <- read.delim("~/Desktop/HCC1395.varscan.tsv", header=FALSE)

# grab only those columns which are required and name them
lohData <- lohData[,c("V1", "V2", "V7", "V11")]
colnames(lohData) <- c("chromosome", "position", "n_vaf", "t_vaf")

# convert the normal and tumor vaf columns to fractions
lohData$n_vaf <- as.numeric(gsub("%", "", lohData$n_vaf))/100
lohData$t_vaf <- as.numeric(gsub("%", "", lohData$t_vaf))/100

# add a sample column
lohData$sample <- "HCC1395"

# limit to just the genome of interest
lohData <- lohData[grepl("^\\d|X|Y", lohData$chromosome),]

# create an inital plot
lohSpec(lohData)

# install and load a benchmarking package
# install.packages("microbenchmark")
library(microbenchmark)

# run benchmark tests
microbenchmark(lohSpec(x=lohData, window_size = 2500000, step = 1000000), lohSpec(x=lohData, window_size = 2500000, step = 1500000), times = 5L)

################################################################################
############################ Exercises #########################################
################################################################################

# exercise 1
# create a custom genome
chr10 <- data.frame("chromosome"="chr10", "start"=0, "end"=135534747)

# create custom layer to highlight q23.1
library(ggplot2)
layer1 <- geom_vline(xintercept=c(89500000, 92900000), colour="chartreuse", linetype=2, size=1)

# create the plot
lohSpec(lohData[lohData$chromosome==10,], y=chr10, plotLayer = layer1)

# exercise2
lohSpec(lohData, colourScheme = "viridis", plotLayer=list(facet_grid(.~chromosome, scales="free"), ggtitle("Loss of Heterozygosity")))
