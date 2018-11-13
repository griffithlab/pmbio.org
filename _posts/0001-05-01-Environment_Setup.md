---
feature_text: |
  ## Precision Medicine Bioinformatics
  Introduction to bioinformatics for DNA and RNA sequence analysis
title: Environment
categories:
    - Module-01-Setup
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-05-01
---

In this section we will set up a number of environment variables and perform other operating system configurations for various reasons of convenience. This might include creating shortcuts for long file paths, setting default values for commonly used parameters, and other configurations. All of these variables will be defined in the `.bashrc` file located in your home directory. This file gets sourced whenever a terminal is opened so the variables that are put in here will automatically be set each time you open a session. You can edit your `.bashrc` file with a text editor of your choice (e.g. `nano` or `vi`, `vim`, etc.).

After you change your `.bashrc` you need to either source it (`source ~/.bashrc`) or exit your current session and log in again.

### Set up environment/path variables for convenience
```bash

# create shortcuts to some long files paths that we will use frequently (e.g. reference file location, etc.)

# create shortcuts to some frequently used JAR files
export PICARD='/usr/local/bin/picard.jar'

# Depending on the instance configuration (i.e., small root volume) it may be necessary to specify a custom temp dir for java processes
# We might also Consider changing this to separate volume for performance
# export _JAVA_OPTIONS=-Djava.io.tmpdir=/workspace/tmp/

# If student wish to use versions of tools the installed, add to path 
# export PATH=/home/ubuntu/workspace/bin:$PATH

# the data subset we will use in this course
export CHRS='chr6_and_chr17'

# GATK regions string for many gatk commands that require it
# export GATK_REGIONS='-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM'
#export GATK_REGIONS='-L chr6 -L chr17'
export GATK_REGIONS='-L chr17' 

```
