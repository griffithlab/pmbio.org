---
feature_text: |
  ## Precision Medicine Bioinformatics
  Introduction to bioinformatics for DNA and RNA sequence analysis
title: Directory Setup
categories:
    - Module-01-Setup
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-06-01
---

In this section we will set up top level directories for the course. 

All exercises for this course will take place in the `/workspace` directory. We have created a symlink so that you can also specify `/home/ubuntu/workspace` or `~/workspace`.  All three of these are equivalent.

On each student cloud instance we have mounted `/workspace` to a large volume (2 TB) so that it can hold large input and results files.

Within `/workspace` we will now create a series of empty directories. These are named according to each section of the workshop. As we move through each module, we will place results relating to each analysis topic in the appropriate directory (e.g. germline variant calls in `/workspace/germline`, rna-seq expression estimates in `/workspace/rnaseq`, etc.).

### Create and view our starting directory structure as follows: 
```bash

# change to the workspace directory where all results will be stored
cd /workspace

# make a directory named for each module of the course
mkdir inputs align germline somatic rnaseq clinical immune cwl

# view the directory structure
tree

# try omitting the bin directory
tree -I bin

```
