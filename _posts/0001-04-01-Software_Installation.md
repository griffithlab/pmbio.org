---
feature_text: |
  ## Precision Medicine
title: Installation
categories:
    - Module-01-Setup
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-04-01
---

This workshop requires a large number of different bioinformatics tools. The instructions for installing these tools exist here. Note that depending on the operating system and environment, some additional dependencies would likely be needed. If you are using the AWS instance built for this course these dependencies have already been installed. However if you are interested in the underlying dependencies and how they were installed see the [AWS AMI Setup](http://pmbio.org/module%2010.%20appendix/0010/02/28/AWS_AMI_Setup/) page. The remainder of this section will assume that you are on the AWS instance, however these instructions should work on any xenial ubuntu distribution with the required dependencies.

### Prepare for installation
For this workshop we will be using the workspace folder to store results, executables, and input files. To start we must choose a single directory for installing tools, typically in linux, user compiled tools are installed in `/usr/local/bin` however backups of the tools we will be using have already been installed there. In this tutorial we will install tools in `~/workspace/bin`. Lets go ahead and make a `bin` directory in `~/workspace` to get started.
```bash
# make a bin directory
mkdir -p ~/workspace/bin
```

### Install Samtools
[Samtools](http://www.htslib.org/) is a software package based in C which provies utilities for manipulating alignment files (SAM/BAM/CRAM). It is open source, available on github, and is under an [MIT license](https://opensource.org/licenses/MIT). Let's go ahead and download the source code from github to our bin directory and extract it with `tar`. Next we need to `cd` into our extracted samtools source code and configure the software. Running `./configure` will make sure all dependencies are available and will also let the software know where it should install to. After that we will need to run `make` to actually build the software. Finally we can run `make install` which will copy the built software and the underlying libraries, documentation, etc. to their final locations. We can check the installation and print out the help message by providing the full path to the executable.
```bash
# change to bin directory
cd ~/workspace/bin

# download and extract the source code
wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2
tar --bzip2 -xvf samtools-1.7.tar.bz2

# configure and compile
cd samtools-1.7/
./configure --prefix=/home/ubuntu/workspace/
make
make install

# check instalation
~/workspace/bin/samtools --help
```

### Install PICARD
[PICARD](https://broadinstitute.github.io/picard/) is a set of java based tools developed by the Broad institute. It is usefull for manipulating next generation sequencing data and is available under an open source [MIT license](https://opensource.org/licenses/MIT). The version of Picard we will be using requires java 8 which has already been installed. All we need to do is download the jar file which is a package file used to distribute java code. We can do this with `wget` from there, to run the software, we simply need to call java with the -jar option and provide the jar file.
```bash
# change to the bin and download the jar file
cd ~/workspace/bin
wget https://github.com/broadinstitute/picard/releases/download/2.18.14/picard.jar

# check the installation
java -jar ~/workspace/bin/picard.jar -h
```

### Install BWA
[BWA](http://bio-bwa.sourceforge.net/) is a popular DNA alignment tool used for mapping sequences to a reference genome. It is available under an open source [GPLv3 license](https://opensource.org/licenses/GPL-3.0). To install BWA, we first need to download and extract the source code. Unlike with samtools theres no `./configure` file so we can just run `make` to build the software. We can then make a symlink with `ln -s` which is just a reference to another file. In this case we will make a symlink so the executable in `~/workspace/bin/bwa-0.7.17/bwa` and be found in `~/workspace/bin/bwa`.
```bash
# change to the bin folder, download, and extract the source code
cd ~/workspace/bin
wget https://cytranet.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
tar --bzip2 -xvf bwa-0.7.17.tar.bz2

# build the software
cd  bwa-0.7.17
make

# make symlink
ln -s ~/workspace/bin/bwa-0.7.17/bwa ~/workspace/bin/bwa

# check the installation
~/workspace/bin/bwa
```

### Install GATK 4
[GATK](https://software.broadinstitute.org/gatk/) is a toolkit developed by the broad institute focused primarily on variant discovery and genotyping. It is open source, hosted on github, and available under a [BSD 3-clause license](https://opensource.org/licenses/BSD-3-Clause). First let's download and unzip GATK from github. The creators of GATK recommend running GATK through [conda](https://conda.io/docs/) which is a package, environment, and dependency management software, in essence conda basically creates a virtual environment from which to run software. The next step then is to tell conda to create a virtual environment for GATK by using the yaml file included within GATK as the instructions for creating the virtual environment. We do this with the command `conda env create`, we also use the `-p` option to specify where this environment should be stored. We will also make a symlink so the executable downloaded is available directly from our `bin` folder. To run GATK we must first start up the virtual environment with the command `source activate`, we can then run the program by providing the path to the executable. To exit the virutal environment run the command `source deactivate`.
```bash
# download and unzip
cd ~/workspace/bin
wget https://github.com/broadinstitute/gatk/releases/download/4.0.2.1/gatk-4.0.2.1.zip
unzip gatk-4.0.2.1.zip

# create conda environment for gatk
cd gatk-4.0.2.1/
conda env create -f gatkcondaenv.yml -p ~/workspace/bin/conda/gatk_4021

# make symlink
ln -s ~/workspace/bin/gatk-4.0.2.1/gatk ~/workspace/bin/gatk

# test installation
source activate ~/workspace/bin/conda/gatk_4021
~/workspace/bin/gatk

# to exit the virtual environment
source deactivate
```

### Install VEP 93.4
[VEP](https://ensembl.org/info/docs/tools/vep/index.html) is a variant annotation tool developed by ensembl and written in [perl](https://www.perl.org/). By default VEP will perform annotations by making web-based API queries however it is much faster to have a local copy of cache and fasta files. The AWS AMI image we're using already has these files as they can take a bit of time to download so to start we'll make a symlink to the directory where the vep_cache containing the HG38 annotation database which already exists. Next we need to download vep from github using `wget` and unzip VEP. From there we can use the INSTALL.pl script vep provides to install the software, the program will ask a series of questions listed below to get a feel for what it's like installing cache and fasta files we'll tell the program to install these for petromyzon_marinus, a much smaller genome.

1. Do you wish to exit so you can get updates (y) or continue (n): n [ENTER]
2. Do you want to continue installing the API (y/n)? y [ENTER]
3. Do you want to install any cache files (y/n)? y [ENTER] 147 [ENTER]
4. Do you want to install any FASTA files (y/n)? y [ENTER] 96 [ENTER]
5. Do you want to install any plugins (y/n)? y [ENTER] 0 [ENTER]

```bash
# symlink the vep_cache directory
mkdir -p ~/workspace/data
ln -s ~/workspace/instructor/data/vep_cache/ ~/workspace/data/vep_cache

# download and unzip vep
cd ~/workspace/bin
wget https://github.com/Ensembl/ensembl-vep/archive/release/93.5.zip
unzip 93.5.zip

# run the INSTALL.pl script provided by VEP
cd ensembl-vep-release-93.5/
perl INSTALL.pl --CACHEDIR ~/workspace/data/vep_cache
#1. Do you wish to exit so you can get updates (y) or continue (n): n [ENTER]
#2. Do you want to continue installing the API (y/n)? y [ENTER]
#3. Do you want to install any cache files (y/n)? y [ENTER] 147 [ENTER]
#4. Do you want to install any FASTA files (y/n)? y [ENTER] 96 [ENTER]
#5. Do you want to install any plugins (y/n)? y [ENTER] 0 [ENTER]

# make a symlink
ln -s ~/workspace/bin/ensembl-vep-release-93.5/vep ~/workspace/bin/vep

# test the Installation
~/workspace/bin/vep --help

# download additional required files for plugins
cd ~/workspace/data/vep_cache/data
wget -c http://krishna.gs.washington.edu/download/CADD/v1.4/GRCh38/InDels.tsv.gz
wget -c http://krishna.gs.washington.edu/download/CADD/v1.4/GRCh38/InDels.tsv.gz.tbi
```

### Install Varscan
[Varscan](http://dkoboldt.github.io/varscan/) is a java program designed to call variants in sequencing data. It was developed at the Genome Institute at Washington University and is hosted on [github](https://github.com/dkoboldt/varscan/). To use Varscan we simply need to download the distributed jar file into our `~/workspace/bin`. As with the other java programs which have already been installed in this section we can envoke Varscan via `java -jar`.
```bash
# Install Varscan
cd ~/workspace/bin
curl -L -k -o VarScan.v2.4.2.jar https://github.com/dkoboldt/varscan/releases/download/2.4.2/VarScan.v2.4.2.jar
java -jar ~/workspace/bin/VarScan.v2.4.2.jar
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
