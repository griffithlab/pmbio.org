---
feature_text: |
  ## Precision Medicine
title: AWS Developers Setup
categories:
    - Module 8
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-02-01
---

This module is primarily for the course developers to document how the AWS AMI and datasets were developed for the course.

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

### Prepare data

The plan is to provide students with raw fastq files as starting point. Download (using wget) the WGS/WES data for HCC1395/BL data from: [https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data](https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data).

```bash
cd ~/data
mkdir unaligned_bams
cd unaligned_bams
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_D1VCPACXX_6.bam
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_D1VCPACXX_6.bam
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_D1VCPACXX_7.bam
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_D1VCPACXX_8.bam
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_D1VCPACXX_1.bam
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_D1VCPACXX_2.bam
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_D1VCPACXX_3.bam
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_D1VCPACXX_4.bam
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_D1VCPACXX_5.bam
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_C1TD1ACXX_7_CGATGT.bam
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_C1TD1ACXX_7_ATCACG.bam
```

For newer RNAseq data (See [https://confluence.gsc.wustl.edu/pages/viewpage.action?spaceKey=CI&title=Cancer+Informatics+Test+Data](https://confluence.gsc.wustl.edu/pages/viewpage.action?spaceKey=CI&title=Cancer+Informatics+Test+Data)) download from inside MGI filesystem.

```bash
scp -i ~/PMB.pem /gscmnt/gc2764/cad/HCC1395/rna_seq/GRCh38/alignments/H_NJ-HCC1395-HCC1395_RNA.bam ubuntu@18.217.114.211:data/unaligned_bams/
scp -i ~/PMB.pem /gscmnt/gc2764/cad/HCC1395/rna_seq/GRCh38/alignments/H_NJ-HCC1395-HCC1395_BL_RNA.bam ubuntu@18.217.114.211:data/unaligned_bams/
```

Rename bam files for easier reference.

```bash
cd /home/ubuntu/data/unaligned_bams
mv gerald_D1VCPACXX_6.bam WGS_Norm_Lane1_gerald_D1VCPACXX_6.bam
mv gerald_D1VCPACXX_7.bam WGS_Norm_Lane2_gerald_D1VCPACXX_7.bam
mv gerald_D1VCPACXX_8.bam WGS_Norm_Lane3_gerald_D1VCPACXX_8.bam
mv gerald_D1VCPACXX_1.bam WGS_Tumor_Lane1_gerald_D1VCPACXX_1.bam
mv gerald_D1VCPACXX_2.bam WGS_Tumor_Lane2_gerald_D1VCPACXX_2.bam
mv gerald_D1VCPACXX_3.bam WGS_Tumor_Lane3_gerald_D1VCPACXX_3.bam
mv gerald_D1VCPACXX_4.bam WGS_Tumor_Lane4_gerald_D1VCPACXX_4.bam
mv gerald_D1VCPACXX_5.bam WGS_Tumor_Lane5_gerald_D1VCPACXX_5.bam
mv gerald_C1TD1ACXX_7_CGATGT.bam Exome_Norm_gerald_C1TD1ACXX_7_CGATGT.bam
mv gerald_C1TD1ACXX_7_ATCACG.bam Exome_Tumor_gerald_C1TD1ACXX_7_ATCACG.bam
mv H_NJ-HCC1395-HCC1395_BL_RNA.bam RNAseq_Norm_H_NJ-HCC1395-HCC1395_BL_RNA.bam 
mv H_NJ-HCC1395-HCC1395_RNA.bam RNAseq_Tumor_H_NJ-HCC1395-HCC1395_RNA.bam
```

Save read group information from each bam file (already separated by read group) for later use

```bash
cd /home/ubuntu/data/unaligned_bams
ls -1 | perl -ne 'chomp; print "samtools view -H $_ | grep -H --label=$_ \@RG\n"' | bash > readgroup_info.txt
cat readgroup_info.txt | perl -ne 'my ($bam, $id, $pl, $pu, $lb, $sm); if ($_=~/(\S+\.bam)\:/){$bam=$1} if ($_=~/(ID\:\d+)/){$id=$1} if ($_=~/(PL\:\w+)/){$pl=$1} if ($_=~/(PU\:\S+)/){$pu=$1} if($_=~/LB\:\"(.+)\"/){$lb=$1} if ($_=~/(SM\:\S+)/){$sm=$1} print "$bam\t$id\t$pl\t$pu\t$lb\t$sm\n";' > readgroup_info.clean.txt
```


Revert bams before conversion to fastq (best practice with MGI bams)
- Run times: ~41-318min
- Memory: 2.4-6.3GB

```bash
cd ~/data
mkdir reverted_bams
cd reverted_bams
mkdir Exome_Norm Exome_Tumor WGS_Norm WGS_Tumor RNAseq_Norm RNAseq_Tumor
java -Xmx16g -jar $PICARD RevertSam I=/home/ubuntu/data/unaligned_bams/Exome_Norm_gerald_C1TD1ACXX_7_CGATGT.bam OUTPUT_BY_READGROUP=true O=/home/ubuntu/data/reverted_bams/Exome_Norm/
java -Xmx16g -jar $PICARD RevertSam I=/home/ubuntu/data/unaligned_bams/Exome_Tumor_gerald_C1TD1ACXX_7_ATCACG.bam OUTPUT_BY_READGROUP=true O=/home/ubuntu/data/reverted_bams/Exome_Tumor/
java -Xmx16g -jar $PICARD RevertSam I=/home/ubuntu/data/unaligned_bams/WGS_Norm_Lane1_gerald_D1VCPACXX_6.bam OUTPUT_BY_READGROUP=true O=/home/ubuntu/data/reverted_bams/WGS_Norm/
java -Xmx16g -jar $PICARD RevertSam I=/home/ubuntu/data/unaligned_bams/WGS_Norm_Lane2_gerald_D1VCPACXX_7.bam OUTPUT_BY_READGROUP=true O=/home/ubuntu/data/reverted_bams/WGS_Norm/
java -Xmx16g -jar $PICARD RevertSam I=/home/ubuntu/data/unaligned_bams/WGS_Norm_Lane3_gerald_D1VCPACXX_8.bam OUTPUT_BY_READGROUP=true O=/home/ubuntu/data/reverted_bams/WGS_Norm/
java -Xmx16g -jar $PICARD RevertSam I=/home/ubuntu/data/unaligned_bams/WGS_Tumor_Lane1_gerald_D1VCPACXX_1.bam OUTPUT_BY_READGROUP=true O=/home/ubuntu/data/reverted_bams/WGS_Tumor/
java -Xmx16g -jar $PICARD RevertSam I=/home/ubuntu/data/unaligned_bams/WGS_Tumor_Lane2_gerald_D1VCPACXX_2.bam OUTPUT_BY_READGROUP=true O=/home/ubuntu/data/reverted_bams/WGS_Tumor/
java -Xmx16g -jar $PICARD RevertSam I=/home/ubuntu/data/unaligned_bams/WGS_Tumor_Lane3_gerald_D1VCPACXX_3.bam OUTPUT_BY_READGROUP=true O=/home/ubuntu/data/reverted_bams/WGS_Tumor/
java -Xmx16g -jar $PICARD RevertSam I=/home/ubuntu/data/unaligned_bams/WGS_Tumor_Lane4_gerald_D1VCPACXX_4.bam OUTPUT_BY_READGROUP=true O=/home/ubuntu/data/reverted_bams/WGS_Tumor/
java -Xmx16g -jar $PICARD RevertSam I=/home/ubuntu/data/unaligned_bams/WGS_Tumor_Lane5_gerald_D1VCPACXX_5.bam OUTPUT_BY_READGROUP=true O=/home/ubuntu/data/reverted_bams/WGS_Tumor/
java -Xmx16g -jar $PICARD RevertSam I=/home/ubuntu/data/unaligned_bams/RNAseq_Norm_H_NJ-HCC1395-HCC1395_BL_RNA.bam OUTPUT_BY_READGROUP=true O=/home/ubuntu/data/reverted_bams/RNAseq_Norm/
java -Xmx16g -jar $PICARD RevertSam I=/home/ubuntu/data/unaligned_bams/RNAseq_Tumor_H_NJ-HCC1395-HCC1395_RNA.bam OUTPUT_BY_READGROUP=true O=/home/ubuntu/data/reverted_bams/RNAseq_Tumor/
```

Convert bams to fastq files

- Run times: 18-58min

```bash
java -Xmx16g -jar $PICARD SamToFastq I=/data/reverted_bams/Exome_Norm/2891351068.bam F=/data/fastqs/Exome_Norm/2891351068_1.fastq F2=/data/fastqs/Exome_Norm/2891351068_2.fastq
java -Xmx16g -jar $PICARD SamToFastq I=/data/reverted_bams/Exome_Tumor/2891351066.bam F=/data/fastqs/Exome_Tumor/2891351066_1.fastq F2=/data/fastqs/Exome_Tumor/2891351066_2.fastq
java -Xmx16g -jar $PICARD SamToFastq I=/data/reverted_bams/WGS_Norm/2891323123.bam F=/data/fastqs/WGS_Norm/2891323123_1.fastq F2=/data/fastqs/WGS_Norm/2891323123_2.fastq
java -Xmx16g -jar $PICARD SamToFastq I=/data/reverted_bams/WGS_Norm/2891323124.bam F=/data/fastqs/WGS_Norm/2891323124_1.fastq F2=/data/fastqs/WGS_Norm/2891323124_2.fastq
java -Xmx16g -jar $PICARD SamToFastq I=/data/reverted_bams/WGS_Norm/2891323125.bam F=/data/fastqs/WGS_Norm/2891323125_1.fastq F2=/data/fastqs/WGS_Norm/2891323125_2.fastq
java -Xmx16g -jar $PICARD SamToFastq I=/data/reverted_bams/WGS_Tumor/2891322951.bam F=/data/fastqs/WGS_Tumor/2891322951_1.fastq F2=/data/fastqs/WGS_Tumor/2891322951_2.fastq
java -Xmx16g -jar $PICARD SamToFastq I=/data/reverted_bams/WGS_Tumor/2891323147.bam F=/data/fastqs/WGS_Tumor/2891323147_1.fastq F2=/data/fastqs/WGS_Tumor/2891323147_2.fastq
java -Xmx16g -jar $PICARD SamToFastq I=/data/reverted_bams/WGS_Tumor/2891323150.bam F=/data/fastqs/WGS_Tumor/2891323150_1.fastq F2=/data/fastqs/WGS_Tumor/2891323150_2.fastq
java -Xmx16g -jar $PICARD SamToFastq I=/data/reverted_bams/WGS_Tumor/2891323174.bam F=/data/fastqs/WGS_Tumor/2891323174_1.fastq F2=/data/fastqs/WGS_Tumor/2891323174_2.fastq
java -Xmx16g -jar $PICARD SamToFastq I=/data/reverted_bams/WGS_Tumor/2891323175.bam F=/data/fastqs/WGS_Tumor/2891323175_1.fastq F2=/data/fastqs/WGS_Tumor/2891323175_2.fastq
java -Xmx16g -jar $PICARD SamToFastq I=/data/reverted_bams/RNAseq_Norm/2895625992.bam F=/data/fastqs/RNAseq_Norm/2895625992_1.fastq F2=/data/fastqs/RNAseq_Norm/2895625992_2.fastq
java -Xmx16g -jar $PICARD SamToFastq I=/data/reverted_bams/RNAseq_Norm/2895626097.bam F=/data/fastqs/RNAseq_Norm/2895626097_1.fastq F2=/data/fastqs/RNAseq_Norm/2895626097_2.fastq
java -Xmx16g -jar $PICARD SamToFastq I=/data/reverted_bams/RNAseq_Tumor/2895626107.bam F=/data/fastqs/RNAseq_Tumor/2895626107_1.fastq F2=/data/fastqs/RNAseq_Tumor/2895626107_2.fastq
java -Xmx16g -jar $PICARD SamToFastq I=/data/reverted_bams/RNAseq_Tumor/2895626112.bam F=/data/fastqs/RNAseq_Tumor/2895626112_1.fastq F2=/data/fastqs/RNAseq_Tumor/2895626112_2.fastq
```

gzip all fastq files

```bash
gzip -v /data/fastqs/*/*.fastq
```

