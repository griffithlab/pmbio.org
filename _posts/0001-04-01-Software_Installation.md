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

### Install GATK 4

```bash
cd ~/bin
wget https://github.com/broadinstitute/gatk/releases/download/4.0.2.1/gatk-4.0.2.1.zip
unzip gatk-4.0.2.1.zip
cd gatk-4.0.2.1/
./gatk
```

### Install gsutil - should these apt-get commands be moved to basic set up of instance?

```bash
sudo apt-get install gcc python-dev python-setuptools libffi-dev
sudo apt-get install python-pip
sudo pip install gsutil
```

### Install VEP

During the installation make sure to accept (y) when asked whether you'd like to install cache files, fastas, and plugins.
Install all homo sapiens build 38 cache files (For VEP v91, this corresponded to options 172, 174 and 176), the homo sapiens fasta (option 40), and the Downstream plugin (option 29) or all plugins (option 0).
If you choose to download all plugins, some will not work without installation (dbNSFP, Carol, Condel, PolyPhen_SIFT, LoF, dbscSNV, GeneSplicer, MaxEntScan) or downloading data (dbNSFP, CADD, FATHMM_MKL, Gwava, LoF, LoFtool, ExACpLI, MPC, MTR, dbscSNV, AncestralAllele, ExAC)

Note, the cache, fasta and plugin files can be quite large, therefore the default (root volume) location may be insufficient. Specify another path with the CACHEDIR option.

Note, VEP natively supports gnomad allele frequencies but it is unclear if this works for all variants or only for dbSNP subset of variants.
See: http://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html#assembly

To access the full gnomAD data set, it is possible to use VEP's custom annotation feature to retrieve the frequency data directly from the gnomAD VCF files. See: http://useast.ensembl.org/info/docs/tools/vep/script/vep_example.html#gnomad

```bash
cd ~/data
mkdir vep_cache
cd ~/bin
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl --CACHEDIR /home/ubuntu/data/vep_cache
```

Download additional data files need for various VEP plugins - CADD, gnomAD, 

```bash
cd ~/data/vep_cache
mkdir data
cd ~/data/vep_cache/data
wget http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz
wget http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz.tbi
wget http://krishna.gs.washington.edu/download/CADD/v1.3/InDels.tsv.gz
wget http://krishna.gs.washington.edu/download/CADD/v1.3/InDels.tsv.gz.tbi

wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz
wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz.tbi
```

