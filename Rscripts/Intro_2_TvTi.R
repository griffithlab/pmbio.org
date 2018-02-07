# load the data into R
mutationData <- read.delim("~/Desktop/ggplot2ExampleData.tsv")

# change the "Simple_name" column to "sample"
colnames(mutationData)[colnames(mutationData) %in% "sample"] <- "sample_1"
colnames(mutationData)[colnames(mutationData) %in% "Simple_name"] <- "sample"

# subset the data to just the discovery cohort
mutationData <- mutationData[mutationData$dataset == "discovery",]
mutationData$sample <- factor(mutationData$sample, levels=unique(mutationData$sample))

# run TvTi
TvTi(mutationData, fileType="MGI")

# install and load gridExtra
# install.packages("gridExtra")
library(gridExtra)
library(grid)

# make a frequency and proportion plot
tvtiFreqGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="frequency")
tvtiPropGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="proportion")

# combine the two plots
finalGrob <- arrangeGrob(tvtiFreqGrob, tvtiPropGrob, ncol=1)

# draw the plot
grid.draw(finalGrob)

# determine the exisitng order of samples
levels(mutationData$sample)

# create a custom sample order to match the manuscript
sampleOrder <- c("LYM145-Naive", "LYM023-Naive", "LYM045-Naive", "LYM238-Sorted", "LYM238-Naive", "LYM058-Naive", "LYM125-Naive", "LYM100-Naive", "LYM088-Naive", "LYM167-Naive", "LYM153-Naive", "LYM002-Naive", "LYM005-Naive", "LYM120-Naive", "LYM120-Treated", "LYM215-Sorted", "LYM215-Treated", "LYM139-Treated", "LYM177-Treated", "LYM202-Sorted", "LYM202-Treated", "LYM062-Treated", "LYM057-tNHL", "LYM001-tNHL", "LYM175-tNHL", "LYM193-tNHL", "LYM237-tNHL", "LYM013-tNHL")

# make sure all samples have been specified
all(mutationData$sample %in% sampleOrder)

# alter the sample levels of the input data
mutationData$sample <- factor(mutationData$sample, levels=sampleOrder)

# recreate the plot
tvtiFreqGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="frequency")
tvtiPropGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="proportion")
finalGrob <- arrangeGrob(tvtiFreqGrob, tvtiPropGrob, ncol=1)
grid.draw(finalGrob)

# read in the clinical data
clinicalData <- read.delim("~/Desktop/FL_ClinicalData.tsv")

# subset just the discovery cohort
clinicalData <- clinicalData[clinicalData$Simple_name %in% mutationData$sample,]
all(sort(unique(as.character(clinicalData$Simple_name))) == sort(unique(as.character(mutationData$sample))))

# convert to long format
library(reshape2)
clinicalData <- melt(data=clinicalData, id.vars="Simple_name")
colnames(clinicalData) <- c("sample", "variable", "value")

# define a few clinical parameter inputs
clin_colors <- c('0'="lightblue",'1'='dodgerblue2','2'='blue1','3'='darkblue','4'='black','FL'='green3','t-NHL'='darkgreen','Treated'='darkred','TreatmentNaive'='red','NA'='lightgrey')

# add the clinical data to the plots
tvtiFreqGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="frequency")
tvtiPropGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="proportion", clinData=clinicalData, clinLegCol = 2, clinVarCol = clin_colors)
finalGrob <- arrangeGrob(tvtiFreqGrob, tvtiPropGrob, ncol=1, heights=c(2,5))
grid.draw(finalGrob)

# grab the widths
proportionPlotWidth <- finalGrob$grobs[[2]]$grobs[[1]]$widths
clinicalPlotWidth <- finalGrob$grobs[[2]]$grobs[[2]]$widths
transitionPlotWidth <- finalGrob$grobs[[1]]$grobs[[1]]$widths

# find the max width of each of these
maxWidth <- unit.pmax(proportionPlotWidth, clinicalPlotWidth, transitionPlotWidth)

# re-assign the max widths
finalGrob$grobs[[2]]$grobs[[1]]$widths <- as.list(maxWidth)
finalGrob$grobs[[2]]$grobs[[2]]$widths <- as.list(maxWidth)
finalGrob$grobs[[1]]$grobs[[1]]$widths <- as.list(maxWidth)

# plot the plot
grid.draw(finalGrob)

# load ggplot2
library(ggplot2)

# create a layer for vertical lines
layer1 <- geom_vline(xintercept=c(14.5, 22.5))

# create a layer to remove redundant x-axis labels
layer2 <- theme(axis.text.x=element_blank(), axis.title.x=element_blank())

# create the plot 
tvtiFreqGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="frequency", layers=list(layer1, layer2))
tvtiPropGrob <- TvTi(mutationData, fileType="MGI", out="grob", type="proportion", clinData=clinicalData, clinLegCol = 2, clinVarCol = clin_colors, layers=layer1, clinLayer=layer1)
finalGrob <- arrangeGrob(tvtiFreqGrob, tvtiPropGrob, ncol=1, heights=c(2, 5))
grid.draw(finalGrob)
