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

importPizzly <- function (filename, genome_version, limit) {
 # Is the genome version valid?
  valid_genomes <- c("hg19", "hg38", "mm10")
  if (is.na(match(tolower(genome_version), tolower(valid_genomes)))) {
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
      "ID" = "character",
      "paircount" = "integer",
      "splitcount" = "integer",
      "geneA.id" = "character",
      "geneA.name" = "character",
      "geneA.seq_name" = "integer",
      "geneA.gene_seq_start" = "integer",
      "geneA.gene_seq_end" = "integer",
      "geneA.seq_strand" = "character",
      "geneB.id" = "character",
      "geneB.name" = "character",
      "geneB.seq_name" = "integer",
      "geneB.gene_seq_start" = "integer",
      "geneB.gene_seq_end" = "integer",
      "geneB.seq_strand" = "character",
      "same_chr" = "character",
      "gene_distance" = "character")
      
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
  id                   <- NA
  inframe              <- NA
  fusion_tool          <- "pizzly"
  spanning_reads_count <- NA
  split_reads_count    <- NA

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import oncofuse-specific fields
    fusion_tool_specific_data <- list()
    fusion_tool_specific_data[["gene_distance"]] <- report[[i, "gene_distance"]]

    # Cluster id
    id <- report[[i, "ID"]]

    # Is the downstream fusion partner in-frame?
#    if (report[[i, "FPG_FRAME_DIFFERENCE"]] == "0") {
#      inframe <- TRUE
#    } else {
#      inframe <- FALSE
#    }

    # Number of supporting reads
    split_reads_count <- report[[i, "splitcount"]]
    spanning_reads_count <- report[[i, "paircount"]]

    chromosome_upstream <- paste("chr", report[[i, "geneA.seq_name"]],
            sep = "")
    chromosome_downstream <- paste("chr", report[[i, "geneB.seq_name"]],
            sep = "")
     breakpoint_upstream <- tryCatch(as.numeric(report[[i,
            "geneA.gene_seq_end"]]), warning = function(w) {
            -1
        })
     breakpoint_downstream <- tryCatch(as.numeric(report[[i,
            "geneB.gene_seq_start"]]), warning = function(w) {
            -1
        })
     strand_upstream <- report[[i, "geneA.seq_strand"]]
     strand_downstream <- report[[i, "geneB.seq_strand"]]

     junction_sequence_upstream <- Biostrings::DNAString()
     junction_sequence_downstream <- Biostrings::DNAString()
     name_upstream <- report[[i, "geneA.id"]]
     name_downstream <- report[[i, "geneB.id"]]
     ensembl_id_upstream <- report[[i, "geneA.name"]]
     ensembl_id_downstream <- report[[i, "geneB.name"]]
     gene_upstream <- new(Class = "PartnerGene", name = name_upstream,
            ensembl_id = ensembl_id_upstream, chromosome = chromosome_upstream,
            breakpoint = breakpoint_upstream, strand = strand_upstream,
            junction_sequence = junction_sequence_upstream, transcripts = GenomicRanges::GRangesList())
        gene_downstream <- new(Class = "PartnerGene", name = name_downstream,
            ensembl_id = ensembl_id_downstream, chromosome = chromosome_downstream,
            breakpoint = breakpoint_downstream, strand = strand_downstream,
            junction_sequence = junction_sequence_downstream,
            transcripts = GenomicRanges::GRangesList())
        fusion_list[[i]] <- new(Class = "Fusion", id = id, fusion_tool = fusion_tool,
            genome_version = genome_version, spanning_reads_count = spanning_reads_count,
            split_reads_count = split_reads_count, fusion_reads_alignment = Gviz::AlignmentsTrack(),
            gene_upstream = gene_upstream, gene_downstream = gene_downstream,
            inframe = inframe, fusion_tool_specific_data = fusion_tool_specific_data)
    }
    fusion_list
}
