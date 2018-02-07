# Load the BioMart library
library("biomaRt")

# Output all available databases for use from BioMart
databases <- listEnsembl()
databases

# Output all available datasets for the "Ensembl Genes" database
ensembl <- useEnsembl(biomart="ensembl")
datasets <- listDatasets(ensembl)
datasets

# Connect to the live Ensembl Gene Human dataset 
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# Output all attributes to see which are available and how they are named
attributes <- listAttributes(ensembl)
attributes

# Output all filters to see which are available and how they are named
filters <- listFilters(ensembl)
filters

# Get Ensembl gene records and relevant attributes for a list of HGNC symbols
hgnc_symbols <- c("NEUROD2", "PPP1R1B", "STARD3", "TCAP", "PNMT", "PGAP3", "ERBB2", "MIR4728", "MIEN1", "GRB7", "IKZF3")
annotations_ENSG <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","external_gene_name","hgnc_symbol"), filter="hgnc_symbol", values=hgnc_symbols, mart=ensembl)
annotations_ENSG
