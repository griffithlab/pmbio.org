################################################################################
###################### Tutorial ################################################
################################################################################

# load extra libraries
# install.packages("stringr")
library("stringr")

# get locations of all cn files
files <- Sys.glob("~/Desktop/*cc2.tsv")

# create function to read in and format data
a <- function(x){
    # read data and set column names
    data <- read.delim(x, header=FALSE)
    colnames(data) <- c("chromosome", "start", "end", "probes", "segmean")
    
    # get the sample name from the file path
    sampleName <- str_extract(x, "H_OM.+cc2")
    sampleName <- gsub(".cc2", "", sampleName)
    data$sample <- sampleName
    
    # return the data
    return(data)
}

# run the anonymous function defined above
cnData <- lapply(files, a)

# turn the list of data frames into a single data frame
cnData <- do.call("rbind", cnData)

# plot a typical plot
cnFreq(cnData, genome="hg19")

# change the CN cutoffs
cnFreq(cnData, genome="hg19", CN_low_cutoff = 0, CN_high_cutoff = .1)

# highlight erbb2
library(ggplot2)
layer1 <- geom_vline(xintercept=c(39709170))
cnFreq(cnData, genome="hg19", plotChr="chr17", plotLayer=layer1)

################################################################################
################################# exercises ####################################

# reset scales 
layer1 <-facet_grid(.~chromosome, scales="free")
cnFreq(cnData, genome="hg19", plotLayer = layer1 )
