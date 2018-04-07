---
feature_text: |
  ## Precision Medicine
title: AWS AMI Setup
categories:
    - Module 1
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-02-01
---

This module is primarily for the course developers to document how the AWS AMI was developed for the course.

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

These steps will update ubuntu packages and install dependencies for software needed in the course. 

Notes:
- picard requires at least java 1.8.x (Didn't check default ubuntu java version/install status)
- samtools (and likely others) require: make, gcc, libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev
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


See [Installation](/module 1/0001/02/01/Software_Installation/) for software tools needed for next parts of this set up.

