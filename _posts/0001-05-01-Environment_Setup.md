---
feature_text: |
  ## Precision Medicine
title: Environment
categories:
    - Module 01. Setup
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-05-01
---

We will set up a number of environment variables and perform other operating system configurations for various reasons of convenience. This might include creating shortcuts for longer file paths, setting default values for commonly used parameters, and other configurations.

### Set up environment/path variables for convenience

```bash
export PICARD='/home/ubuntu/bin/picard.jar'
export _JAVA_OPTIONS=-Djava.io.tmpdir=/data/tmp/ #Consider changing this to separate volume for performance
export PATH=/home/ubuntu/bin/samtools-1.7/bin:/home/ubuntu/bin/bwa-0.7.17:/home/ubuntu/bin/gatk-4.0.2.1:$PATH
```

