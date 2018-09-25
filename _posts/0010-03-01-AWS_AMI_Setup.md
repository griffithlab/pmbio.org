---
feature_text: |
  ## Precision Medicine
title: AWS AMI Setup
categories:
    - Module 10. Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0010-03-01
---

This module is primarily for the course developers to document how the AWS AMI was developed for the course. Students will start from an AMI where some system setup will be done already, but students will still learn to install all necesary bioinformatics files. If students are interested

***

### Initial AWS setup for development and testing purposes

For development purposes we started with a very large instance (overkill). Future experimentation is needed to determine the appropriate size for actual student instances.

- Launch an EC2 instance:
- Select Ubuntu Server 16.04 LTS (HVM), SSD Volume Type - ami-2581aa40
- Choose r4.16xlarge (64 vCPUs, 488 GiB Memory, 25 Gigabit Network Performance)
- Add storage: 10,000 GiB (~10TB) EBS volume, not encrypted
- Configure security: Allow SSH access
- Login with key the usual way (e.g., ssh -i PMB.pem ubuntu@18.217.114.211)

### Perform basic linux configuration
Many of the tools used also have underlying dependencies, in many linux distributions these packages will already be installed and available. In this AMI setup however we start from a very basic Ubuntu distrubtion and we will have to install these dependencies. Ubuntu is based on the Debian operating system and so we can use the Debian based package manager `apt-get` for installation.

```bash
# general tools for installation
sudo apt-get update -y && sudo apt-get install -y \
     wget \
     bzip2
```

```bash
# Samtools
sudo apt-get update -y && sudo apt-get install -y \
     build-essential \
     libncurses5-dev \
     zlib1g-dev \
     libbz2-dev \
     liblzma-dev
```

```bash
# PICARD
sudo apt-get update -y && sudo apt-get install -y \
     openjdk-8-jdk
```

```bash
# BWA
sudo apt-get update -y && sudo apt-get install -y \
```

Notes:
- For performance reasons it may be desirable to create an instance with larger root volume and/or a separate tmp volume

```bash
sudo apt-get update
sudo apt-get upgrade
sudo apt-get -y install make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev default-jdk apache2 unzip tabix python-dev python-setuptools libffi-dev python-pip cpanminus
sudo pip install gsutil cmake
```

Install mysql-server - set a root password (e.g., pmbiotest)

```bash
sudo apt-get install mysql-server libmysqlclient-dev
```

Install perl dependencies (e.g., for VEP)

```bash
sudo cpanm DBI
sudo cpanm DBD::mysql
sudo cpanm Bio::Root::Version
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
