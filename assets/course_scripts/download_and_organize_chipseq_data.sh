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

# The BAM files downloaded will be restricted to the permissions set by the original file owner
# We can reset these permissions using `chmod` so that we can download them and work with them easily.
chmod ugo+r *.bam








