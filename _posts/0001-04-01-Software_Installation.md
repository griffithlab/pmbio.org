---
feature_text: |
  ## Precision Medicine
title: Installation
categories:
    - Module 1
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-04-01
---

This workshop requires a large number of different bioinformatics tools. The instructions for installing these tools exist here. Note that depending on the operating system and environment, some additional dependencies would likely be needed. See our [AWS Setup](/module )

### Prepare for installation
First, choose a single directory for installing tools

```bash
cd ~
mkdir bin
```

### Install Picard

```bash
cd ~/bin
wget https://github.com/broadinstitute/picard/releases/download/2.17.6/picard.jar
export PICARD='/home/ubuntu/bin/picard.jar'
java -jar $PICARD -h
```

### Install Samtools

```bash
cd ~/bin
wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2
tar --bzip2 -xvf samtools-1.7.tar.bz2
cd samtools-1.7/
./configure --prefix=/home/ubuntu/bin/samtools-1.7
make
make install
./samtools
```

### Install BWA
```bash
cd ~/bin
wget https://cytranet.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
tar --bzip2 -xvf bwa-0.7.17.tar.bz2
cd  bwa-0.7.17
make
./bwa
```

### Environment setup

```bash
export PICARD='/home/ubuntu/bin/picard.jar'
export _JAVA_OPTIONS=-Djava.io.tmpdir=/data/tmp/ #Consider changing this to separate volume for performance
export PATH=/home/ubuntu/bin/samtools-1.7/bin:/home/ubuntu/bin/bwa-0.7.17:$PATH
```

