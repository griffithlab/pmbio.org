---
feature_text: |
  ## Precision Medicine
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
export PICARD='/home/ubuntu/bin/picard.jar'

# specify a custom temp dir for java processes
export _JAVA_OPTIONS=-Djava.io.tmpdir=/workspace/tmp/ #Consider changing this to separate volume for performance

# add some locally installed tools to our path
export PATH=/home/ubuntu/bin/samtools-1.7/bin:/home/ubuntu/bin/bwa-0.7.17:/home/ubuntu/bin/gatk-4.0.2.1:$PATH

# the data subset we will use in this course
export CHRS='chr6_and_chr17'

```
