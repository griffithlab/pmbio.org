#' Import results from a pizzly run into a list of Fusion objects.
#'
#' A function that imports the results from a pizzly run, typically from a
#' grolar .tsv file, into a list of Fusion objects.
#' 
#' This import function was contributed by jeszyman
#'
#' @param filename Filename for the oncofuse results .tsv file.
#' @param genomeVersion Which genome was used in mapping- (must be hg19, hg38, or mm10).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' ADD EXAMPLE
#'
#' @export

importPizzly <- function (filename, genomeVersion, limit) {

  # Is the genome version valid?
  valid_genomes <- c("hg19", "hg38", "mm10")
  if (is.na(match(tolower(genomeVersion), tolower(valid_genomes)))) {
    stop("Invalid genome version given")
  }

  # If the limit is set, is the value valid?
  if (missing(limit) == FALSE) {
    if (is.numeric(limit) == FALSE || limit <= 0) {
      stop("limit must be a numeric value bigger than 0")
    }
  }

  # Try to read the fusion report
  report <- withCallingHandlers({
      col_types <- c(
          "ID" = "integer",
          "geneB.seq_name"="character",
          "geneA.seq_name"="character",
          "geneA.gene_seq_end"="numeric",
          "geneB.gene_seq_start"="numeric",
          "geneA.id"="character",
          "geneB.id"="character",
          "paircount" = "integer"
      )
      if (missing(limit)) {
        # Read all lines
        data.table::fread(
          input = filename,
          colClasses = col_types,
          showProgress = FALSE
        )
      } else {
        # Only read up to the limit
        data.table::fread(
          input = filename,
          colClasses = col_types,
          showProgress = FALSE,
          nrows = limit
        )
      }
    },
    error = function(cond) {
      message(paste0("Reading ", filename, " caused an error: ", cond[[1]]))
      stop(cond)
    },
    warning = function(cond) {
      # Begin Exclude Linting
      message(paste0("Reading ", filename, " caused a warning: ", cond[[1]]))
      # End Exclude Linting
    }
    )

  # Set variables
    id <- NA
    inframe <- NA
    fusionTool <- "pizzly"
    spanningReadsCount <- NA
    splitReadsCount <- NA
    junctionSequence <- NA
    fusionReadsAlignment <- NA
    # List to hold all Fusion objects
    fusionList <- vector("list", dim(report)[1])
    for (i in 1:dim(report)[1]) {

        fusionToolSpecificData = list()
        fusionToolSpecificData[["distance"]] = report[[i, "gene_distance"]]
    # ID
    id = as.character(report[[i, "ID"]])
    # Chromosome names
    chromosomeA <- paste("chr", report[[i, "geneA.seq_name"]], sep = "")
    chromosomeB <- paste("chr", report[[i, "geneB.seq_name"]], sep = "")

 # Breakpoints
        breakpointA <- tryCatch(breakpointA <- as.numeric(report[[i, 
            "geneA.gene_seq_end"]]), warning = function(w) {
            breakpointA <- -1
        })
        breakpointB <- tryCatch(breakpointB <- as.numeric(report[[i, 
            "geneB.gene_seq_start"]]), warning = function(w) {
            breakpointB <- -1
        })

    
    # Strand
    strandA <- as.character(report[[i, "geneA.seq_strand"]])
    strandB <- as.character(report[[i, "geneB.seq_strand"]])

    # Number of supporting reads
    splitReadsCount <- report[[i, "splitcount"]]
 
        spanningReadsCount <- tryCatch(spanningReadsCount <- as.numeric(report[[i, 
            "paircount"]]), warning = function(w) {
            spanningReadsCount <- 0
        })

    
     # SOAPfuse doesn't provide the fusion sequence
    junctionSequenceA <- Biostrings::DNAString()
    junctionSequenceB <- Biostrings::DNAString()

        # Gene names
    nameA <- report[[i, "geneA.name"]]
    nameB <- report[[i, "geneB.name"]]

    # Ensembl ids
        ensemblIdA <- NA_character_
        ensemblIdB <- NA_character_
    # PartnerGene objects
    geneA=new(Class = "PartnerGene",
              name = nameA,
              ensemblId=ensemblIdA,
              chromosome = chromosomeA,
              breakpoint= breakpointA,
              strand=strandA,
              junctionSequence=junctionSequenceA,
              transcripts = GenomicRanges::GRangesList()
              )
    # PartnerGene objects
    geneB=new(Class = "PartnerGene",
              name = nameB,
              ensemblId=ensemblIdB,
              chromosome = chromosomeB,
              breakpoint= breakpointB,
              strand=strandB,
              junctionSequence=junctionSequenceB,
              transcripts = GenomicRanges::GRangesList()
              )

          fusionList[[i]] <- new(Class = "Fusion", id = id, fusionTool = fusionTool, 
            genomeVersion = genomeVersion, spanningReadsCount = spanningReadsCount, 
            splitReadsCount = splitReadsCount, fusionReadsAlignment = Gviz::AlignmentsTrack(), 
            geneA = geneA, geneB = geneB, inframe = inframe, 
            fusionToolSpecificData = fusionToolSpecificData)
    }
    fusionList
}
    
    
   

