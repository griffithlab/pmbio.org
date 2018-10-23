---
feature_text: |
  ## Precision Medicine
title: Installation
categories:
    - Module-01-Setup
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-04-01
---

This workshop requires a large number of different bioinformatics tools. The instructions for installing these tools exist here. Note that depending on the operating system and environment, some additional dependencies would likely be needed. If you are using the AWS instance built for this course these dependencies have already been installed. However if you are interested in the underlying dependencies and how they were installed see the [AWS AMI Setup]({{ site.baseurl }}{% link _posts/0010-03-01-AWS_AMI_Setup.md %}) page. The remainder of this section will assume that you are on the AWS instance, however these instructions should work on any xenial ubuntu distribution with the required dependencies.

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
conda env create -f gatkcondaenv.yml -p ~/workspace/bin/conda/gatk

# make symlink
ln -s ~/workspace/bin/gatk-4.0.2.1/gatk ~/workspace/bin/gatk

# test installation
source activate ~/workspace/bin/conda/gatk
~/workspace/bin/gatk

# to exit the virtual environment
source deactivate
```

### Install VEP 93.4
[VEP](https://ensembl.org/info/docs/tools/vep/index.html) is a variant annotation tool developed by ensembl and written in [perl](https://www.perl.org/). By default VEP will perform annotations by making web-based API queries however it is much faster to have a local copy of cache and fasta files. The AWS AMI image we're using already has these files for hg38 in the directory `/opt/vep_cache` as they can take a bit of time to download. To get an idea of what it's like to install these we will install a vep_cache for petromyzon_marinus, a much smaller genome. so to start we need to download vep from github using `wget` and unzip VEP. From there we can use the INSTALL.pl script vep provides to install the software which will ask a series of questions listed below. We also make a symlink when the installer completes.

1. Do you wish to exit so you can get updates (y) or continue (n): n [ENTER]
2. Do you want to continue installing the API (y/n)? y [ENTER]
3. Do you want to install any cache files (y/n)? y [ENTER] 147 [ENTER]
4. Do you want to install any FASTA files (y/n)? y [ENTER] 96 [ENTER]
5. Do you want to install any plugins (y/n)? y [ENTER] 0 [ENTER]

```bash
# download and unzip vep
cd ~/workspace/bin
wget https://github.com/Ensembl/ensembl-vep/archive/release/93.5.zip
unzip 93.5.zip

# run the INSTALL.pl script provided by VEP
cd ensembl-vep-release-93.5/
perl INSTALL.pl --CACHEDIR /opt/vep_cache
#1. Do you wish to exit so you can get updates (y) or continue (n): n [ENTER]
#2. Do you want to continue installing the API (y/n)? y [ENTER] (if asked)
#3. Do you want to install any cache files (y/n)? y [ENTER] 147 [ENTER]
#4. Do you want to install any FASTA files (y/n)? y [ENTER] 96 [ENTER]
#5. Do you want to install any plugins (y/n)? n [ENTER]

# make a symlink
ln -s ~/workspace/bin/ensembl-vep-release-93.5/vep ~/workspace/bin/vep

# test the Installation
~/workspace/bin/vep --help
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
[BCFtools](https://samtools.github.io/bcftools/) is an open source program for variant calling and manipulating files in Variant Call Format (VCF) or Binary Variant Call Format (BCF). To install we first need to download and extract the source code with `curl` and `tar` respectively. We can then call `make` to build the program and `make install` to copy the program to the desired directory.
```bash
# download and extract
cd ~/workspace/bin
curl -L -k -o bcftools-1.3.1.tar.bz2                https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2
tar --bzip2 -xvf bcftools-1.3.1.tar.bz2

# install the software
cd bcftools-1.3.1
make -j
make prefix=~/workspace/ install

# test installation
~/workspace/bin/bcftools -h
```

### Install Strelka
[Strekla](https://github.com/Illumina/strelka) is a germline and somatic variant caller developed by illumina and available under an open source [GPLv3 license](https://opensource.org/licenses/GPL-3.0). The binary distribution for strelka is already built and hosted on github so to install all we have to do is download and extract the software. It is important to note that strelka is built on `python 2` and won't work for python 3. The AMI we're using contains both python versions so we just have to make sure we invoke strelka with python2, you can view the python versions on the AMI with `python2 --version` and `python3 --version`.
```bash
# download and extract
cd ~/workspace/bin
curl -L -k -o strelka-2.7.1.centos5_x86_64.tar.bz2 https://github.com/Illumina/strelka/releases/download/v2.7.1/strelka-2.7.1.centos5_x86_64.tar.bz2
tar --bzip2 -xvf strelka-2.7.1.centos5_x86_64.tar.bz2

# test installation
python2 ~/workspace/bin/strelka-2.7.1.centos5_x86_64/bin/configureStrelkaWorkflow.py -h
```

### Install Sambamba
[Sambamba](http://lomereiter.github.io/sambamba/) is a high performance alternative to samtools and provides a subset of samtools functionality. It is up to 6x faster for duplicate read marking and 4x faster for viewing alignment files. To install sambamba we can just download the binary distribution and extract it. From there we just make a symlink to make using it a bit more intuitive.
```bash
# download and extract
cd ~/workspace/bin
curl -L -k -o sambamba_v0.6.4_linux.tar.bz2 https://github.com/lomereiter/sambamba/releases/download/v0.6.4/sambamba_v0.6.4_linux.tar.bz2
tar --bzip2 -xvf sambamba_v0.6.4_linux.tar.bz2

# create symlink
ln -s ~/workspace/bin/sambamba_v0.6.4 ~/workspace/bin/sambamba

# test installation
~/workspace/bin/sambamba
```

### Install HISAT2
[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) is a graph based alignment algorithm devoloped at Johns Hopkins University. It is heavily used in the bioinformatics community for RNAseq based alignments. To Install we will need to download and extract the binary executable. We then make a symlink to put it with the other executables we've installed.
```bash
# download and extract
cd ~/workspace/bin
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-Linux_x86_64.zip
unzip hisat2-2.0.4-Linux_x86_64.zip

# create symlink
ln -s ~/workspace/bin/hisat2-2.0.4/hisat2 ~/workspace/bin/hisat2

# test installation
~/workspace/bin/hisat2 --help
```

### Install StringTie
[StringTie](https://ccb.jhu.edu/software/stringtie/) is a software program to perform transcript assembly and quantification of RNAseq data. The binary distributions are available so to install we can just download this distribution and extract it. Like with our other programs we also make a symlink to make it easier to find.
```bash
# download and extract
cd ~/workspace/bin
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.0.Linux_x86_64.tar.gz
tar -xzvf stringtie-1.3.0.Linux_x86_64.tar.gz

# make symlink
ln -s ~/workspace/bin/stringtie-1.3.0.Linux_x86_64/stringtie ~/workspace/bin/stringtie

# test installation
~/workspace/bin/stringtie -h
```
### Install Gffcompare
[Gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) is a prgram that is used to perform operations on general feature format (GFF) and general transfer format (GTF) files. It has a binary distribution compatible with the linux we're using so we will just download, extract, and make a symlink.
```bash
# download and extract
cd ~/workspace/bin
wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.9.8.Linux_x86_64.tar.gz
tar -xzvf gffcompare-0.9.8.Linux_x86_64.tar.gz

# make symlink
ln -s ~/workspace/bin/gffcompare-0.9.8.Linux_x86_64/gffcompare ~/workspace/bin/gffcompare

# check Installation
~/workspace/bin/gffcompare
```
### Install R
[R](https://www.r-project.org/) is an feature rich interpretive programming language originally released in 1995. It is heavily used in the bioinformatics community largely due to numerous R libraries available on [bioconductor](https://www.bioconductor.org/). It takes a couple minutes to compile so we'll use one which has already been setup if we were to install R we first would need to download and extract the source code. Next we'd configure the installation with `--with-x=no` which tells R to install without X11, a windowing system for displays. We'd also specify `--prefix` which is where the R.framework will go, this includes the additional R libraries we'll download later. From there we'd do make and make install to build the software and copy the files to their proper location and create symlinks for the executables. Finally we'd install the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) and [Biocmanager](https://cran.r-project.org/web/packages/BiocManager/index.html) packages from the command line to make installing additional packages easier. We've commented out the code below however it is exactly what was run to set up the R we will be using, except the installation location.
```bash
## download and extract
#cd ~/workspace/bin
#wget https://cran.r-project.org/src/base/R-3/R-3.5.1.tar.gz
#tar -zxvf R-3.5.1.tar.gz

## configure the installation, build the code
#cd R-3.5.1
#./configure --prefix=/home/ubuntu/workspace/bin --with-x=no
#make
#make install

## make symlinks
#ln -s ~/workspace/bin/R-3.5.1/bin/Rscript ~/workspace/bin/Rscript
#ln -s ~/workspace/bin/R-3.5.1/bin/R ~/workspace/bin/R

## test installation
#cd ~/workspace/bin
#~/workspace/bin/Rscript --version

## install additional packages
#~/workspace/bin/R --vanilla -e 'install.packages(c("devtools", "BiocManager"), repos="http://cran.us.r-project.org")'
```

### Install copyCat
[copyCat](https://github.com/chrisamiller/copyCat) is an R library for detecting copy number aberrations in sequencing data. The library is only available on github so we will have to use the [BiocManager](https://cran.r-project.org/web/packages/BiocManager/index.html) library to install a few of the underlying package dependencies. If installing a package from cran or bioconductor these dependencies would be automatically installed. After these dependencies are installed we can use the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package to install copycat directory from its github repository.
```bash
# Install R Library dependencies
R --vanilla -e 'BiocManager::install(c("IRanges", "DNAcopy"))'

# install copyCat
R --vanilla -e 'devtools::install_github("chrisamiller/copycat")'
```

### Install CNVnator
[CNVnator](https://github.com/abyzovlab/CNVnator) is a depth based copy number caller. It is open source and available on github under a creative common public license (CCPL). To install we first download and extract the source code. CNVnator relies on a specific version of samtools which is distributed with CNVnator, so our first step is to run `make` on that samtools. To finnish the installation process we can then run `make` in CNVnator's main source directory.
```bash
# download and decompress
cd ~/workspace/bin
wget https://github.com/abyzovlab/CNVnator/releases/download/v0.3.3/CNVnator_v0.3.3.zip
unzip CNVnator_v0.3.3.zip

# make the samtools dependency distributed with CNVnator
cd CNVnator_v0.3.3/src/samtools
make

# make CNVnator
cd ../
make

# make a symlink
ln -s ~/workspace/bin/CNVnator_v0.3.3/src/cnvnator ~/workspace/bin/cnvnator

# test installation
~/workspace/bin/cnvnator
```

### Install CNVkit
[CNVkit](https://cnvkit.readthedocs.io/en/stable/) is a python based copy number caller designed for use with hybrid capture. To install we can download and extract the package. We then must use [conda](https://conda.io/docs/) to set up the environment to run cnvkit. This process, while straight forward, takes some time so we've commented out the installation instructions for this tool and will use the conda environment that has already been set up.
```bash
## download and unzip
cd ~/workspace/bin
wget https://github.com/etal/cnvkit/archive/v0.9.5.zip
unzip v0.9.5.zip

## add conda channels
#conda config --add channels defaults
#conda config --add channels conda-forge
#conda config --add channels bioconda

## create conda environment
#conda create -n cnvkit cnvkit

# test installation
source activate cnvkit
python ~/workspace/bin/cnvkit-0.9.5/cnvkit.py --help

# to exit the virtual environment
source deactivate
```

### Install Kallisto
[Kallisto](https://pachterlab.github.io/kallisto/) is a kmer-based alignment algorithm used for quantifying transcripts in RNAseq data. Kallisto has a binary distribution available so to use the program we only have to download and extract the software from github.
```bash
# download and extract
cd ~/workspace/bin
wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz
tar -zxvf kallisto_linux-v0.44.0.tar.gz

# make symlink
ln -s ~/workspace/bin/kallisto_linux-v0.44.0/kallisto ~/workspace/bin/kallisto

# test installation
~/workspace/bin/kallisto --help
```

### Install Pizzly
[Pizzly](https://github.com/pmelsted/pizzly) is a fusion detection algorithm which uses output from Kallisto. Pizzly is has a binary distribution so we can download and extract that from github to get started.
```bash
# download and extract
cd ~/workspace/bin
mkdir pizzly-v0.37.3
cd pizzly-v0.37.3
wget https://github.com/pmelsted/pizzly/releases/download/v0.37.3/pizzly_linux.tar.gz
tar -zxvf pizzly_linux.tar.gz

# test executable
~/workspace/bin/pizzly --help
```

### Manta
[Manta](https://github.com/Illumina/manta) is a structural variant caller developed by Illumina and available on gitub under the GPL_v3 license. It uses paired-end sequencing reads to build a breakend association graph to identify structural varaints.
```bash
# download and extract
cd ~/workspace/bin
wget https://github.com/Illumina/manta/releases/download/v1.4.0/manta-1.4.0.centos6_x86_64.tar.bz2
tar --bzip2 -xvf manta-1.4.0.centos6_x86_64.tar.bz2

# test installation
python2 ~/workspace/bin/manta-1.4.0.centos6_x86_64/bin/configManta.py --help
```

### mosdepth
[mosdepth](https://github.com/brentp/mosdepth) is a program for determining depth in sequencing data. The easiest way to install mosdepth is through `bioconda` a channel for the `conda` package manager. The AMI already has conda setup to install to `/usr/local/bin/miniconda` and so we've already installed mosdepth for you. However below are the commands used during the installation.
```bash
# add the bioconda channel
conda config --add channels bioconda

# install mosdepth with the conda package manager
conda install mosdepth
```

### bam-readcount
[bam-readcount](https://github.com/genome/bam-readcount) is a program for determing read support for individual variants (SNVs and Indels only). We are going to point this local install of bam-readcount to use the samtools installation we completed above. Samtools is a dependency of bam-readcount. This tool uses Cmake to create its makefile, so compiling from source has an extra step here. Instead of using an official release from github we are cloning the latest code from the master branch. In general this practice should be avoided and you should use an official release instead.
```bash
# install bam-readcount
cd ~/workspace/bin
git clone https://github.com/genome/bam-readcount.git
mv bam-readcount bam-readcount-latest
cd bam-readcount-latest
export SAMTOOLS_ROOT=/home/ubuntu//workspace/bin/samtools-1.7
cmake -Wno-dev /home/ubuntu/workspace/bin/bam-readcount-latest
make

# test installation
~/workspace/bin/bam-readcount/bin/bam-readcount

```


