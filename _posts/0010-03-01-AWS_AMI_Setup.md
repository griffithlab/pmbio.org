---
feature_text: |
  ## Precision Medicine
title: AWS AMI Setup
categories:
    - Module 10. Appendix
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

### Software Dependencies
Many of the tools used also have underlying dependencies, in many linux distributions these packages will already be installed and available. In this AMI setup however we start from a very basic Ubuntu distrubtion and we will have to install these dependencies. Ubuntu is based on the Debian operating system and so we can use the Debian based package manager `apt-get` for installation. Below the stand-alone dependencies required for each bioinformatic tool used in the course is supplied. In this we will use the normal linux convention where our own compiled binaries and executables are installed in `/usr/local/bin`.

#### Pre-Installation
Describes the general system wide dependencies required for downloading and decompressing source and binary files related to the tools to be installed.
```bash
# general tools for installation
sudo apt-get update -y && sudo apt-get install -y \
     wget \
     bzip2 \
     unzip \
     git \
     curl
```

#### Samtools 1.7
Describes dependencies for samtools 1.7, used in this course for general bam file manipulation.
```bash
# Samtools
sudo apt-get update -y && sudo apt-get install -y \
     build-essential \
     libncurses5-dev \
     zlib1g-dev \
     libbz2-dev \
     liblzma-dev
```

#### PICARD 2.18.14
Describes dependencies for PICARD 2.18.14, used in this course for general bam file manipulation and QC.
```bash
# PICARD
sudo apt-get update -y && sudo apt-get install -y \
     openjdk-8-jdk
```

#### BWA 0.7.17
Describes dependencies for BWA 0.7.17, used in this course for DNA alignment.
```bash
# BWA
sudo apt-get update -y && sudo apt-get install -y \
    build-essential \
    libz-dev
```

#### GATK 4.0.2.1
Describes dependencies for GATK 4.0.2.1, used in this course for .....
```bash
# GATK 4
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda env create -n gatk -f gatkcondaenv.yml
source activate gatk
sudo apt-get update -y && sudo apt-get install -y \
     openjdk-8-jdk
# for full functionality R and the libraries gsalib, ggplot2, reshape, gplots should be installed
```

#### VEP 93.4
Describes dependencies for VEP 93.4, used in this course for variant annotation.
```bash
# VEP
sudo apt-get update -y && sudo apt-get install -y \
     libdbi-perl \
     libdbd-mysql-perl \
     build-essential \
     zlib1g-dev \
     libmodule-build-perl \
     cpanminus

cpanm -i Bio::Root::Version
```

#### VarScan 2.4.2
Describes dependencies for VarScan 2.4.2, used in this course for variant calling
```bash
# VarScan
sudo apt-get update -y && sudo apt-get install -y \
     openjdk-8-jdk
```

#### BCFtools 1.3.1
Describes dependencies for BCFtools 1.3.1, used in this course for manipulating VCF files.
```bash
# BCFtools
sudo apt-get update -y && sudo apt-get install -y \
     build-essential \
     libz-dev
```

#### Strelka 2.7.1
Describes dependencies for strelka, used in this course for variant calling.
```bash
# strelka
sudo apt-get update -y && sudo apt-get install -y \
     python-dev
```

#### Sambamba 0.6.4
Describes dependencies for Sambamba 0.6.4, used in this course for ....
```bash
# sambamba
```

#### HISAT 2.0.4
Describes dependencies for HISAT 2.0.4, used in this course for RNA alignment.
```bash
# hisat2
```

#### StringTie 1.3.0
Describes dependencies for StringTie 1.3.0, used in this course for transcript abundance estimates.
```bash
# stringtie
```

#### Gffcompare 0.9.8
Describes dependencies for Gffcompare 0.9.8, used in this course for ....
```bash
# Gffcompare
```

#### R 3.5.1
Describes dependencies for R 3.5.1, used in this course for general file manipulation/analysis.
```bash
# R
sudo apt-get update -y && sudo apt-get install -y \
     gfortran \
     libreadline-dev \
     libpcre3-dev \
     libcurl4-openssl-dev \
     build-essential \
     zlib1g-dev \
     libbz2-dev \
     liblzma-dev \
     openjdk-8-jdk

# devtools
sudo apt-get update -y && sudo apt-get install -y \
     libssl-dev \
     libxml2-dev
```

#### copyCat 1.6.12
Describes dependencies for copyCat 1.6.12, used in this course for WGS copy number calling.
```bash
# copyCat
R-3.5.1/bin/R --vanilla -e 'BiocManager::install(c("GenomicRanges", "Organism.dplyr"))'
```

#### CNVnator
Describes dependencies for CNVnator, used in this course for exome copy number calling.
```bash
# CNVnator
sudo apt-get update -y && sudo apt-get install -y \
     libxpm4
wget https://root.cern.ch/download/root_v6.14.04.Linux-ubuntu16-x86_64-gcc5.4.tar.gz
tar -xzvf root_v6.14.04.Linux-ubuntu16-x86_64-gcc5.4.tar.gz
export ROOTSYS=/bin/root
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH

wget http://bitbucket.org/MDukhan/yeppp/downloads/yeppp-1.0.0.tar.bz2
tar -xvjf yeppp-1.0.0.tar.bz2
export YEPPPLIBDIR=/bin/yeppp-1.0.0/binaries/linux/x86_64
export YEPPPINCLUDEDIR=/bin/yeppp-1.0.0/library/headers
export LD_LIBRARY_PATH=$YEPPPLIBDIR:$LD_LIBRARY_PATH
```

#### cnvkit
Describes dependencies for cnvkit, used in this course for ...
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

Notes:
- For performance reasons it may be desirable to create an instance with larger root volume and/or a separate tmp volume

### Perform basic linux configuration

Install mysql-server - set a root password (e.g., pmbiotest)

```bash
sudo apt-get install mysql-server libmysqlclient-dev
```

Format and mount an extra data volume. Create symlink from homedir for convenience. Create a tmp dir on the larger volume for tools (e.g., picard) that need more temp space than the default system temp dir.

```bash
lsblk
sudo mkfs -t ext4 /dev/xvdb
sudo mkdir /data
sudo mount /dev/xvdb /data
sudo chown -R ubuntu:ubuntu /data

cd ~
ln -s /data data

mkdir /data/tmp
```

Make ephemeral storage mounts persistent. See [http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ebs-using-volumes.html](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ebs-using-volumes.html) for guidance on setting up fstab records for AWS.

```bash
echo -e "LABEL=cloudimg-rootfs / ext4 defaults,discard 0 0\n/dev/xvdb /data ext4 defaults,nofail 0 2" | sudo tee /etc/fstab
```

Set up apache web server for convenient access to files. First, edit config to allow files to be served from /data/.

```bash
sudo vim /etc/apache2/apache2.conf
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
sudo vim /etc/apache2/sites-available/000-default.conf
```

Change document root in 000-default.conf
```bash
DocumentRoot /data
```

Restart apache
```bash
sudo service apache2 restart
```

### Formatting and mounting storage volumnes
With the tools we will be using now installed we need to mount and format the storage volume we allocated when we initialized the instance. Students will use this volume to install their own copies of tools used as well as input data and results.
```bash
# make a workspace directory in root
cd /
sudo mkdir workspace

# format the mount
sudo mkfs /dev/xvdb

# mount the drive with the allocated space
sudo mount /dev/xvdb /workspace

# make the mount persistent
echo -e "LABEL=cloudimg-rootfs / ext4 defaults,discard 0 0\n/dev/xvdb /workspace ext4 defaults,nofail 0 2" | sudo tee /etc/fstab

```
