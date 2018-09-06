---
feature_text: |
  ## Precision Medicine
title: RNAseq Expression Estimation
categories:
    - Module 5
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-02-01
---

#### **Downloading RNA BAMs**
__________________________
First, if not done so already, make a separate folder named `/data/RNA_seq` and a subfolder called `fastqs_RNA` and download the RNA_seq data from `pmbio.org` to your instance.\
In order to prevent storage space from running out, you may want to unzip the files sequentially and delete the original zipped file once the unzipped file have been obtained.

You will need to have the following software installed, including HISAT, Sambamba, StringTie, Gffcompare, R. If you are missing any of the following software, or you run into problems with running commands using your currently installed version, please refer to the installation instructions below.\
For this tutorial, the software packages were installed in the directory `/data/RNA_seq/software`, please adjust the commands accordingly with your path to the packages.

#### **Software Installation**
__________________________
* HISAT 
  * `wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-Linux_x86_64.zip`
  * `unzip hisat2-2.0.4-Linux_x86_64.zip`
  * `ln -s /data/RNA_seq/software/hisat2-2.0.4/hisat2 /data/RNA_seq/software/bin/hisat2`
  * `cd /data/RNA_seq/software/bin and run ./hisat2 for testing`
* Sambamba 
  * `curl -L -k -o sambamba_v0.6.4_linux.tar.bz2 https://github.com/lomereiter/sambamba/releases/download/v0.6.4/sambamba_v0.6.4_linux.tar.bz2`
  * `tar --bzip2 -xvf sambamba_v0.6.4_linux.tar.bz2`
  * `ln -s /data/RNA_seq/software/sambamba_v0.6.4 /data/RNA_seq/software/sambamba`
  * `./sambamba`
* StringTie
  * `wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.0.Linux_x86_64.tar.gz`
  * `tar -xzvf stringtie-1.3.0.Linux_x86_64.tar.gz`
  * `ln -s /data/RNA_seq/software/stringtie-1.3.0.Linux_x86_64/stringtie /data/RNA_seq/software/bin/stringtie`
* Gffcompare 
  * `wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.9.8.Linux_x86_64.tar.gz`
  * `tar -xzvf gffcompare-0.9.8.Linux_x86_64.tar.gz`
  * `ln -s  /data/RNA_seq/software/gffcompare-0.9.8.Linux_x86_64/gffcompare  /data/RNA_seq/software/bin/gffcompare`
* Additional Step before R:
  * `sudo apt-get install gfortran`
  * `sudo apt-get install libreadline-dev`
  * `sudo apt install libpcre3-dev`
  * `sudo apt-get install libcurl4-openssl-dev`
* R
  * `export R_LIBS=`
  * `wget https://cran.r-project.org/src/base/R-3/R-3.4.0.tar.gz`
  * `tar -zxvf R-3.4.0.tar.gz`
  * `cd R-3.4.0`
  * `./configure --prefix=/data/RNA_seq/software --with-x=no`
  * `make`
  * `make install`
  * `sudo apt install r-base-core`
  * `Rscript --version`
  * `sudo apt-get install libssl-dev`
  * `sudo apt-get install libxml2-dev`
  * `/data/RNA_seq/software/bin/R`
  * In the interactive section:\
    ```
    install.packages("devtools",repos="http://cran.us.r-project.org")
    source("http://bioconductor.org/biocLite.R")
    biocLite(c("alyssafrazee/RSkittleBrewer","dplyr","genefilter","ballgown"))
    quit()
    ```



