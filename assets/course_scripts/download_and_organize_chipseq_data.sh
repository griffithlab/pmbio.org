#!/usr/bin/env

# Download chipseq tar
wget -O chipseq_data.tar.gz "https://app.box.com/index.php?rm=box_download_shared_file&shared_name=pu626zopurawjb16fmy97auwi2a5kq0m&file_id=f_1138767053169"

# Untar chipseq_data.tar.gz  -x extracts files ; -v lists files as they are extracted/ processed ; -z unzips the files ; -f use archive file
# This will make a new folder called chipseq_data with the BAM files in the tar
tar -xzvf chipseq_data.tar.gz

# Move BAM files inside chipseq_data to the current directory for ease of access
mv chipseq_data/* .

# Delete the now empty chipseq_data folder
rm -r chipseq_data

# Make a directory to hold the reference
mkdir references

# Download the hg38 reference that corresponds to this data to the reference folder
wget -O references/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz "https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz"









