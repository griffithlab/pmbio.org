################################################################################
################## Tutorial ####################################################

# install gage
# source("https://bioconductor.org/biocLite.R")
# biocLite("gage")
library(gage)

# load the differential expression results fro the previous section
load(url("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/deseq2Data_v1.RData"))

# extract the results from the deseq2 data
library(DESeq2)
tumor_v_normal_DE <- results(deseq2Data, contrast=c("tissueType", "primary colorectal cancer", "normal-looking surrounding colonic epithelium"))

# set up kegg database
kg.hsa <- kegg.gsets("hsa")
kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]

# set up go database
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]

# load in libraries to annotate data
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("AnnotationDbi","org.Hs.eg.db"))
library(AnnotationDbi)
library(org.Hs.eg.db)

# annotate the deseq2 results with additional gene identifiers
tumor_v_normal_DE$symbol <- mapIds(org.Hs.eg.db, keys=row.names(tumor_v_normal_DE), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
tumor_v_normal_DE$entrez <- mapIds(org.Hs.eg.db, keys=row.names(tumor_v_normal_DE), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
tumor_v_normal_DE$name <- mapIds(org.Hs.eg.db, keys=row.names(tumor_v_normal_DE), column="GENENAME", keytype="ENSEMBL", multiVals="first")

# grab the log fold changes for everything
tumor_v_normal_DE.fc <- tumor_v_normal_DE$log2FoldChange
names(tumor_v_normal_DE.fc) <- tumor_v_normal_DE$entrez

# Run enrichment analysis on all log fc
fc.kegg.sigmet.p <- gage(tumor_v_normal_DE.fc, gsets = kegg.sigmet.gs)
fc.kegg.dise.p <- gage(tumor_v_normal_DE.fc, gsets = kegg.dise.gs)
fc.go.bp.p <- gage(tumor_v_normal_DE.fc, gsets = go.bp.gs)
fc.go.mf.p <- gage(tumor_v_normal_DE.fc, gsets = go.mf.gs)
fc.go.cc.p <- gage(tumor_v_normal_DE.fc, gsets = go.cc.gs)

# covert the kegg results to data frames
fc.kegg.sigmet.p.up <- as.data.frame(fc.kegg.sigmet.p$greater)
fc.kegg.dise.p.up <- as.data.frame(fc.kegg.dise.p$greater)

fc.kegg.sigmet.p.down <- as.data.frame(fc.kegg.sigmet.p$less)
fc.kegg.dise.p.down <- as.data.frame(fc.kegg.dise.p$less)

# convert the go results to data frames
fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
fc.go.mf.p.up <- as.data.frame(fc.go.mf.p$greater)
fc.go.cc.p.up <- as.data.frame(fc.go.cc.p$greater)

fc.go.bp.p.down <- as.data.frame(fc.go.bp.p$less)
fc.go.mf.p.down <- as.data.frame(fc.go.mf.p$less)
fc.go.cc.p.down <- as.data.frame(fc.go.cc.p$less)

################################################################################
################ Exercise 1 ####################################################

#!!!!!!!!!! "fc.go.bp.p.up" and "go.bp.gs" should already be in your R environment, see tutorial above !!!!!!!!!#

# get genes associated with pathways
a <- function(sigPathways, geneSetDB){

    # grab the names of all significant pathways
    sigPathwaysID <- rownames(sigPathways)
    
    # subset the geneSet list to only those in the significant pathways
    sigPathwaysGenes <- geneSetDB[which(names(geneSetDB) %in% sigPathwaysID)]
    numberSigPathways <- length(sigPathwaysGenes)
    sigPathwaysGenes <- unlist(sigPathwaysGenes)
    
    # count the number of times a gene occurs in these significant pathways
    sigPathwaysTable <- plyr::count(sigPathwaysGenes)
    
    # annotate these final results with the gene symbol and some extra information
    sigPathwaysTable$symbol <- mapIds(org.Hs.eg.db, keys=as.character(sigPathwaysTable$x), column="SYMBOL", keytype="ENTREZID", multiVals="first")
    sigPathwaysTable$proportion <- sigPathwaysTable$freq/numberSigPathways
    
    # sort and return the results
    sigPathwaysTable <- sigPathwaysTable[order(-sigPathwaysTable$proportion),]
    
    return(sigPathwaysTable)
}

# run the function defined above
geneTable <- a(na.omit(fc.go.bp.p.up[fc.go.bp.p.up$q.val <= .05,]), go.bp.gs)