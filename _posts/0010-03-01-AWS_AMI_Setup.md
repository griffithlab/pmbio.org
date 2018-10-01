---
feature_text: |
  ## Precision Medicine
title: AWS AMI Setup
categories:
    - Module-10-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0010-03-01
---

This module is primarily for the course developers to document how the AWS AMI was developed for the course. Students will start from an AMI where some system setup will be done already, but students will still learn to install all necesary bioinformatics files.

***

### Initial AWS setup for development and testing purposes

For development purposes we started with a very large instance (overkill). Future experimentation is needed to determine the appropriate size for actual student instances.

- Launch an EC2 instance:
- Select Ubuntu Server 16.04 LTS (HVM), SSD Volume Type - ami-2581aa40
- Choose r4.16xlarge (64 vCPUs, 488 GiB Memory, 25 Gigabit Network Performance)
- Add storage: 10,000 GiB (~10TB) EBS volume, not encrypted
- Configure security: Allow SSH access
- Login with key the usual way (e.g., ssh -i PMB.pem ubuntu@18.217.114.211)

### Formatting and mounting storage volumnes
With the tools we will be using now installed we need to mount and format the storage volume we allocated when we initialized the instance. Students will use this volume to install their own copies of tools used as well as input data and results. See [http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ebs-using-volumes.html](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ebs-using-volumes.html) for guidance on setting up fstab records for AWS.
```bash
# start sudo shell
sudo bash

# get structure
lsblk

# make a workspace directory in root
cd /
mkdir workspace

# format the mount
mkfs /dev/xvdb

# mount the drive with the allocated space
mount /dev/xvdb /workspace
chown -R ubuntu:ubuntu /workspace

# Make ephemeral storage mounts persistent
echo -e "LABEL=cloudimg-rootfs / ext4 defaults,discard 0 0\n/dev/xvdb /workspace ext4 defaults,nofail 0 2" | tee /etc/fstab

# make symlink for convenience
cd ~
ln -s /workspace workspace
```

### Software Dependencies
Many of the tools used also have underlying dependencies, in many linux distributions these packages will already be installed and available. In this AMI setup however we start from a very basic Ubuntu distrubtion and we will have to install these dependencies. Ubuntu is based on the Debian operating system and so we can use the Debian based package manager `apt-get` for installation. Below the stand-alone dependencies required for each bioinformatic tool used in the course is supplied. In this we will use the normal linux convention where our own compiled binaries and executables are installed in `/usr/local/bin`.

#### Pre-Installation
Describes the general system wide dependencies required for downloading and decompressing source and binary files related to the tools to be installed.
```bash
# start sudo shell
sudo bash

# general tools for installation
cd /usr/local/bin
apt-get update -y && apt-get install -y \
  wget \
  bzip2 \
  unzip \
  git \
  curl
```

#### Samtools 1.7
Describes dependencies for samtools 1.7, used in this course for general bam file manipulation.
```bash
# start sudo shell
sudo bash

# samtools dependencies
cd /usr/local/bin
apt-get update -y && apt-get install -y \
  build-essential \
  libncurses5-dev \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev

# samtools installation
wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2
tar --bzip2 -xvf samtools-1.7.tar.bz2
cd samtools-1.7
./configure --prefix=/usr/local/
make
make install
```

#### PICARD 2.18.14
Describes dependencies for PICARD 2.18.14, used in this course for general bam file manipulation and QC.
```bash
# start sudo shell
sudo bash

# picard dependencies
cd /usr/local/bin
apt-get update -y && apt-get install -y \
  openjdk-8-jdk

# picard installation
wget https://github.com/broadinstitute/picard/releases/download/2.18.14/picard.jar
export PICARD='/usr/local/bin/picard.jar'
```

#### BWA 0.7.17
Describes dependencies for BWA 0.7.17, used in this course for DNA alignment.
```bash
# start sudo shell
sudo bash

# bwa dependencies
cd /usr/local/bin
apt-get update -y && apt-get install -y \
  build-essential \
  libz-dev

# bwa installation
wget https://cytranet.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
tar --bzip2 -xvf bwa-0.7.17.tar.bz2
cd  bwa-0.7.17
make
ln -s /usr/local/bin/bwa-0.7.17/bwa /usr/local/bin/bwa
```

#### GATK 4.0.2.1
Describes dependencies for GATK 4.0.2.1, used in this course for .....
```bash
# start sudo shell
sudo bash

# install miniconda dependency
cd /usr/local/bin
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh # choose /usr/local/bin/miniconda as install location
source ~/.bashrc

# install additional dependencies
apt-get update -y && apt-get install -y \
  openjdk-8-jdk

# install GATK
wget https://github.com/broadinstitute/gatk/releases/download/4.0.2.1/gatk-4.0.2.1.zip
unzip gatk-4.0.2.1.zip
cd /usr/local/bin/gatk-4.0.2.1
conda env create -n gatk -f gatkcondaenv.yml
# to use gatk: source activate gatk
# for full functionality R and the libraries gsalib, ggplot2, reshape, gplots should be installed
```

#### VEP 93.4
Describes dependencies for VEP 93.4, used in this course for variant annotation.
- Note: VEP natively supports gnomad allele frequencies but it is unclear if this works for all variants or only for dbSNP subset of variants.
See: http://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html#assembly

- To access the full gnomAD data set, it is possible to use VEP's custom annotation feature to retrieve the frequency data directly from the gnomAD VCF files. See: http://useast.ensembl.org/info/docs/tools/vep/script/vep_example.html#gnomad

```bash
# start sudo shell
sudo bash

# Install VEP dependencies
cd /usr/local/bin
apt-get update -y && apt-get install -y \
  libdbi-perl \
  libdbd-mysql-perl \
  build-essential \
  zlib1g-dev \
  libmodule-build-perl \
  cpanminus
cpanm -i Bio::Root::Version

# install vep with the various plugins
mkdir -p /workspace/instructor/data/vep_cache
git clone https://github.com/Ensembl/ensembl-vep.git
wget https://github.com/Ensembl/ensembl-vep/archive/release/93.5.zip
unzip 93.5.zip
cd ensembl-vep-release-93.5/
perl INSTALL.pl --CACHEDIR /workspace/instructor/data/vep_cache # install cache, hg38:refseq,vep,merged
#TODO need to figure out which plugins we are using and download data for them
```

#### VarScan 2.4.2
Describes dependencies for VarScan 2.4.2, used in this course for variant calling
```bash
# start sudo shell
sudo bash

# varscan dependencies
cd /usr/local/bin
apt-get update -y && apt-get install -y \
  openjdk-8-jdk

# install varscan
curl -L -k -o VarScan.v2.4.2.jar https://github.com/dkoboldt/varscan/releases/download/2.4.2/VarScan.v2.4.2.jar
```

#### BCFtools 1.3.1
Describes dependencies for BCFtools 1.3.1, used in this course for manipulating VCF files.
```bash
# start sudo shell
sudo bash

# BCFtools dependencies
cd /usr/local/bin
apt-get update -y && apt-get install -y \
  build-essential \
  libz-dev

# install BCFtools
curl -L -k -o bcftools-1.3.1.tar.bz2                https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2
tar --bzip2 -xvf bcftools-1.3.1.tar.bz2
cd bcftools-1.3.1
make -j
make prefix=/usr/local/ install
```

#### Strelka 2.7.1
Describes dependencies for strelka, used in this course for variant calling.
```bash
# start sudo shell
sudo bash

# strelka dependencies
cd /usr/local/bin
apt-get update -y && apt-get install -y \
  python-dev

curl -L -k -o strelka-2.7.1.centos5_x86_64.tar.bz2 https://github.com/Illumina/strelka/releases/download/v2.7.1/strelka-2.7.1.centos5_x86_64.tar.bz2
tar --bzip2 -xvf strelka-2.7.1.centos5_x86_64.tar.bz2 # note uses python2
```

#### Sambamba 0.6.4
Describes dependencies for Sambamba 0.6.4, used in this course for ....
```bash
# start sudo shell
sudo bash

# install sambamba
cd /usr/local/bin
curl -L -k -o sambamba_v0.6.4_linux.tar.bz2 https://github.com/lomereiter/sambamba/releases/download/v0.6.4/sambamba_v0.6.4_linux.tar.bz2
tar --bzip2 -xvf sambamba_v0.6.4_linux.tar.bz2
ln -s /usr/local/bin/sambamba_v0.6.4 /usr/local/bin/sambamba
```

#### HISAT 2.0.4
Describes dependencies for HISAT 2.0.4, used in this course for RNA alignment.
```bash
# start sudo shell
sudo bash

# install hisat2
cd /usr/local/bin
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-Linux_x86_64.zip
unzip hisat2-2.0.4-Linux_x86_64.zip
ln -s /usr/local/bin/hisat2-2.0.4/hisat2 /usr/local/bin/hisat2
```

#### StringTie 1.3.0
Describes dependencies for StringTie 1.3.0, used in this course for transcript abundance estimates.
```bash
# start sudo shell
sudo bash

# install stringtie
cd /usr/local/bin
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.0.Linux_x86_64.tar.gz
tar -xzvf stringtie-1.3.0.Linux_x86_64.tar.gz
ln -s /usr/local/bin/stringtie-1.3.0.Linux_x86_64/stringtie /usr/local/bin/stringtie
```

#### Gffcompare 0.9.8
Describes dependencies for Gffcompare 0.9.8, used in this course for ....
```bash
# start sudo shell
sudo bash

# intall Gff compare
cd /usr/local/bin
wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.9.8.Linux_x86_64.tar.gz
tar -xzvf gffcompare-0.9.8.Linux_x86_64.tar.gz
ln -s /usr/local/bin/gffcompare-0.9.8.Linux_x86_64/gffcompare /usr/local/bin/gffcompare
```

#### R 3.5.1
Describes dependencies for R 3.5.1, used in this course for general file manipulation/analysis.
```bash
# start sudo shell
sudo bash

# Install R Dependencies
cd /usr/local/bin
apt-get update -y && apt-get install -y \
  gfortran \
  libreadline-dev \
  libpcre3-dev \
  libcurl4-openssl-dev \
  build-essential \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  openjdk-8-jdk

# Install R
wget https://cran.r-project.org/src/base/R-3/R-3.5.1.tar.gz
tar -zxvf R-3.5.1.tar.gz
cd R-3.5.1
./configure --prefix=/usr/local/ --with-x=no
make
make install

# devtools and BiocManager dependencies
apt-get update -y && apt-get install -y \
  libssl-dev \
  libxml2-dev

R --vanilla -e 'install.packages(c("devtools", "BiocManager"), repos="http://cran.us.r-project.org")'
```

#### copyCat 1.6.12
Describes dependencies for copyCat 1.6.12, used in this course for WGS copy number calling.
```bash
# start sudo shell
sudo bash

# Install copyCat dependencies (see R section for installing R)
cd /usr/local/bin
R --vanilla -e 'BiocManager::install(c("IRanges", "DNAcopy"))'

# Install copyCat
R --vanilla -e 'devtools::install_github("chrisamiller/copycat")'
```

#### CNVnator
Describes dependencies for CNVnator, used in this course for exome copy number calling.
```bash
# start sudo shell
sudo bash

# Install CNVnator dependency root
cd /usr/local/bin
apt-get update -y && apt-get install -y \
  build-essential \
  libncurses5-dev \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  build-essential \
  libxpm4
wget https://root.cern.ch/download/root_v6.14.04.Linux-ubuntu16-x86_64-gcc5.4.tar.gz
tar -xzvf root_v6.14.04.Linux-ubuntu16-x86_64-gcc5.4.tar.gz
export ROOTSYS=/usr/local/bin/root
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH

# Install CNVnator dependency yeppp
wget http://bitbucket.org/MDukhan/yeppp/downloads/yeppp-1.0.0.tar.bz2
tar -xvjf yeppp-1.0.0.tar.bz2
export YEPPPLIBDIR=/usr/local/bin/yeppp-1.0.0/binaries/linux/x86_64
export YEPPPINCLUDEDIR=/usr/local/bin/yeppp-1.0.0/library/headers
export LD_LIBRARY_PATH=$YEPPPLIBDIR:$LD_LIBRARY_PATH

# install CNVnator
wget https://github.com/abyzovlab/CNVnator/releases/download/v0.3.3/CNVnator_v0.3.3.zip
unzip CNVnator_v0.3.3.zip
cd CNVnator_v0.3.3/src/samtools
make
cd ../
make
```

#### cnvkit
Describes dependencies for cnvkit, used in this course for ...
```bash
# start sudo shell
sudo bash

# install cnvkit
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

### apache web serve setup
Set up apache web server for convenient access to files. First, edit config to allow files to be served from /data/.

```bash
# start sudo shell
sudo bash

vim /etc/apache2/apache2.conf
```

Add the following content to apache2.conf
```bash
<Directory /data/>
       Options Indexes FollowSymLinks
       AllowOverride None
       Require all granted
</Directory>
```

Edit vhost file

```bash
vim /etc/apache2/sites-available/000-default.conf
```

Change document root in 000-default.conf
```bash
DocumentRoot /workspace
```

Restart apache
```bash
service apache2 restart
```
