---
feature_text: |
  ## Precision Medicine
title: Installation
categories:
    - Module 01. Setup
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-04-01
---

This workshop requires a large number of different bioinformatics tools. The instructions for installing these tools exist here. Note that depending on the operating system and environment, some additional dependencies would likely be needed. If you are using the AWS instance built for this course these dependencies have already been installed. However if you are interested in the underlying dependencies and how they were installed see the [AWS Setup](http://pmbio.org/module%2010.%20appendix/0010/02/28/AWS_AMI_Setup/) page. The remainder of this section will assume that you are on the AWS instance, however these instructions should work on any xenial ubuntu distribution.

### Prepare for installation
First we must choose a single directory for installing tools, typically in linux user compiled tools are installed in `/usr/local/bin` however it doesn't really matter where our tools are installed as long as we can find them. In this tutorial we will install tools in `~/bin`, the AMI has copies of all these tools in `/usr/local/bin` as well. For this workshop we will be using the workspace folder to store all of our results. Lets go ahead and make a `bin` directory in there as well for the tools we will be installing.
```bash
cd ~/workspace
mkdir bin
```

### Install PICARD
PICARD is a set of java based tools developed by the Broad institute. It is usefull for manipulating next generation sequencing data and is available under an open source MIT license. Go ahead and follow the instructions below to download the pre compiled jar file.
```bash
cd ~/workspace/bin
wget https://github.com/broadinstitute/picard/releases/download/2.18.14/picard.jar
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

### Install Varscan
```bash
cd ~/bin
curl -L -k -o VarScan.v2.4.2.jar https://github.com/dkoboldt/varscan/releases/download/2.4.2/VarScan.v2.4.2.jar
java -jar ~/bin/VarScan.v2.4.2.jar
```


### Install BCFtools
```bash
cd ~/bin
curl -L -k -o bcftools-1.3.1.tar.bz2                https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2
tar --bzip2 -xvf bcftools-1.3.1.tar.bz2
cd bcftools-1.3.1
make -j
make prefix=~/bin/software install

./bcftools -h
```

### Install Strelka
```bash
cd ~/bin
curl -L -k -o strelka-2.7.1.centos5_x86_64.tar.bz2 https://github.com/Illumina/strelka/releases/download/v2.7.1/strelka-2.7.1.centos5_x86_64.tar.bz2
tar --bzip2 -xvf strelka-2.7.1.centos5_x86_64.tar.bz2

.strelka-2.7.1.centos5_x86_64/bin/configureStrelkaWorkflow.py -h
```

### Install Sambamba

```bash
cd ~/bin
curl -L -k -o sambamba_v0.6.4_linux.tar.bz2 https://github.com/lomereiter/sambamba/releases/download/v0.6.4/sambamba_v0.6.4_linux.tar.bz2
tar --bzip2 -xvf sambamba_v0.6.4_linux.tar.bz2
ln -s ~/bin/sambamba_v0.6.4 ~/bin/sambamba

./sambamba
```

### Install HISAT

```bash
cd ~/bin
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-Linux_x86_64.zip
unzip hisat2-2.0.4-Linux_x86_64.zip
ln -s ~/bin/sambamba_v0.6.4 ~/bin/sambamba

./hisat2-2.0.4/hisat2
```

### Install StringTie

```bash
cd ~/bin
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.0.Linux_x86_64.tar.gz
tar -xzvf stringtie-1.3.0.Linux_x86_64.tar.gz
ln -s ~/bin/stringtie-1.3.0.Linux_x86_64/stringtie ~/bin/stringtie

./stringtie -h
```
### Install Gffcompare

```bash
cd ~/bin
wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.9.8.Linux_x86_64.tar.gz
tar -xzvf gffcompare-0.9.8.Linux_x86_64.tar.gz
ln -s ~/bin/gffcompare-0.9.8.Linux_x86_64/gffcompare ~/bin/gffcompare

./gffcompare
```
### Install R

```bash
cd ~/bin

export R_LIBS=
wget https://cran.r-project.org/src/base/R-3/R-3.5.1.tar.gz
tar -zxvf R-3.5.1.tar.gz
cd R-3.5.1
./configure --prefix=/home/ubuntu/bin/ --with-x=no
make
make install

cd ~/bin
./R-3.5.1/bin/Rscript --version

R-3.5.1/bin/R --vanilla -e 'install.packages(c("devtools", "BiocManager"), repos="http://cran.us.r-project.org")'
```

### Install copyCat
```
R-3.5.1/bin/R --vanilla -e 'devtools::install_github("chrisamiller/copycat")'
```

### Install CNVnator
```bash
wget https://github.com/abyzovlab/CNVnator/releases/download/v0.3.3/CNVnator_v0.3.3.zip
unzip CNVnator_v0.3.3.zip
cd CNVnator_v0.3.3/src/samtools
make
cd ../
make
```

### Install cnvkit
```bash
cd ~/bin
wget https://github.com/etal/cnvkit/archive/v0.9.5.zip
unzip v0.9.5.zip
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n cnvkit cnvkit
```
