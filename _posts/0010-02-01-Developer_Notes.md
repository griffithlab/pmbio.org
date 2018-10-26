---
feature_text: |
  ## Precision Medicine
title: Developer Notes
categories:
    - Module-10-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0010-02-01
---

This module is used to document background details that are generally considered too obscure for use in the main workshop but are helpful for the course developers to keep track of certain details.

### Set up reference genome files and store on genomedata.org
```

# set up genome references dir
cd /workspace/
mkdir -p references/genome/
cd references/genome/

# download original references file 
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.dict
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/20150713_location_of_centromeres_and_other_regions.txt
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla-extra.fa
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/README.20150309.GRCh38_full_analysis_set_plus_decoy_hla
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
mkdir all
mkdir chr6
mkdir chr17
mkdir chr6_and_chr17

# store original full genome and rename files for consistency
mv 20150713_location_of_centromeres_and_other_regions.txt all/location_of_centromeres_and_other_regions.txt
mv GRCh38_full_analysis_set_plus_decoy_hla-extra.fa all/ref_genome-extra.fa
mv GRCh38_full_analysis_set_plus_decoy_hla.dict all/ref_genome.dict
mv GRCh38_full_analysis_set_plus_decoy_hla.fa all/ref_genome.fa
mv GRCh38_full_analysis_set_plus_decoy_hla.fa.fai all/ref_genome.fa.fai
mv README.20150309.GRCh38_full_analysis_set_plus_decoy_hla all/README.txt

# split ref genome into pieces
mkdir split/
faSplit byname all/ref_genome.fa split/
mv split/chr6.fa chr6/ref_genome.fa
mv split/chr17.fa chr17/ref_genome.fa
cat chr6/ref_genome.fa chr17/ref_genome.fa > chr6_and_chr17/ref_genome.fa
rm -fr split

# create .fai files for each version of the reference
samtools faidx chr6/ref_genome.fa
samtools faidx chr17/ref_genome.fa
samtools faidx chr6_and_chr17/ref_genome.fa

# create .dict files for each version of the reference
java -jar /usr/local/bin/picard.jar CreateSequenceDictionary R=chr6/ref_genome.fa O=chr6/ref_genome.dict
java -jar /usr/local/bin/picard.jar CreateSequenceDictionary R=chr17/ref_genome.fa O=chr17/ref_genome.dict
java -jar /usr/local/bin/picard.jar CreateSequenceDictionary R=chr6_and_chr17/ref_genome.fa O=chr6_and_chr17/ref_genome.dict

# compress the reference files for storage on the file server
gzip all/ref_genome.fa
gzip chr6/ref_genome.fa
gzip chr17/ref_genome.fa
gzip chr6_and_chr17/ref_genome.fa

# create tarballs for convenient downloading
cd /workspace/references/genome/all
tar -cf ref_genome.tar *
rm -f README.txt location_of_centromeres_and_other_regions.txt ref_genome-extra.fa ref_genome.dict ref_genome.fa.gz ref_genome.fa.fai

cd /workspace/references/genome/chr6
tar -cf ref_genome.tar *
rm -f ref_genome.dict ref_genome.fa.gz ref_genome.fa.fai

cd /workspace/references/genome/chr17
tar -cf ref_genome.tar *
rm -f ref_genome.dict ref_genome.fa.gz ref_genome.fa.fai

cd /workspace/references/genome/chr6_and_chr17
tar -cf ref_genome.tar *
rm -f ref_genome.dict ref_genome.fa.gz ref_genome.fa.fai

```

### Create reference transcriptome files and store on genomedata.org for use in the course
Download transcriptome annotations (GTF files) from Ensembl. Make sure the version used matches the version of VEP installed on the AMI! Be sure to fix the chromosome names in the GTF to be compatible with our reference genome. Finally create sub-setted GTF files for the regions of interest used in the course.

```bash
# create dir for these annotation files
cd /home/ubuntu/workspace/inputs/references/
mkdir transcriptome
cd transcriptome

# download the GTF file for the desired version of Ensembl
wget ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz
gunzip Homo_sapiens.GRCh38.93.gtf.gz

# create a version of the GTF with the chr names fixed for our reference genome (i.e. with chr names, etc.)
# in order to do this we need a .dict file for the whole reference genome
wget http://genomedata.org/pmbio-workshop/references/genome/all/ref_genome.dict
convertEnsemblGTF.pl ref_genome.dict /opt/vep_cache/homo_sapiens/93_GRCh38/chr_synonyms.txt Homo_sapiens.GRCh38.93.gtf > Homo_sapiens.GRCh38.93.namefixed.gtf 
rm -f ref_genome.dict

# produce GTF files of various subsets
mkdir all chr6 chr17 chr6_and_chr17
mv Homo_sapiens.GRCh38.93.gtf all/ref_transcriptome_nochrs.gtf
mv Homo_sapiens.GRCh38.93.namefixed.gtf all/ref_transcriptome.gtf
cat all/ref_transcriptome.gtf | grep --color=never -w "^chr6" > chr6/ref_transcriptome.gtf
cat all/ref_transcriptome.gtf | grep --color=never -w "^chr17" > chr17/ref_transcriptome.gtf
cat chr6/ref_transcriptome.gtf chr17/ref_transcriptome.gtf > chr6_and_chr17/ref_transcriptome.gtf

# produce transcriptome (cDNA) fasta files of various subsets using our GTF and genome Fasta files


```




### Prepare original starting data

The plan is to provide students with raw down-sampled fastq files as starting point. These notes document our original source of data files and any transformations needed.

Download (using wget) the WGS/WES data for HCC1395/BL data from: [https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data](https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data).

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

Create tarball of all fastq files to host for later use
```bash
cd /data/
tar -cvf fastqs.tar fastqs/
```

### Create downsampled data sets to allow faster analysis a teaching setting
Starting with aligned data, the following steps were used to create down-sampled data files, example shown below for chr6 exome data:
```bash
#Starting off with slices of bam file with chr6 as the region specified. This contains all mapped reads that fall/overlap(?) with this region. We want to ensure that if read 1 is in this region that when reconstructing the fastq, we would also like to include read 2, thus we have the following approach.
cd <directory containing Exome alignment>

samtools view chr6_Exome_Norm_sorted_mrkdup_bqsr.bam | cut -d$'\t' -f1 > chr6_Exome_Norm_readnames.txt
samtools view chr6_Exome_Tumor_sorted_mrkdup_bqsr.bam | cut -d$'\t' -f1 > chr6_Exome_Tumor_readnames.txt

java -jar /usr/local/bin/picard.jar FilterSamReads I=Exome_Tumor_sorted_mrkdup_bqsr.bam O=chr6_Exome_Tumor_all_read_pairs.bam READ_LIST_FILE=chr6_Exome_Tumor_readnames.txt FILTER=includeReadList
java -jar /usr/local/bin/picard.jar FilterSamReads I=../final/Exome_Norm_sorted_mrkdup_bqsr.bam O=chr6_Exome_Norm_all_read_pairs.bam READ_LIST_FILE=chr6_Exome_Norm_readnames.txt FILTER=includeReadList

ls -1 | grep chr6_Exome_Tumor_all_read_pairs.bam | perl -ne 'chomp; print "samtools view -H $_ | grep -H --label=$_ \@RG\n"' | bash > chr6_Exome_Tumor_readgroup_info.txt
cat chr6_Exome_Tumor_readgroup_info.txt | perl -ne 'my ($bam, $id, $pl, $pu, $lb, $sm); if ($_=~/(\S+\.bam)\:/){$bam=$1} if ($_=~/(ID\:\d+)/){$id=$1} if ($_=~/(PL\:\w+)/){$pl=$1} if ($_=~/(PU\:\S+)/){$pu
=$1} if($_=~/LB\:\"(.+)\"/){$lb=$1} if ($_=~/(SM\:\S+)/){$sm=$1} print "$bam\t$id\t$pl\t$pu\t$lb\t$sm\n";' > chr6_Exome_Tumor_readgroup_info.clean.txt

ls -1 | grep chr6_Exome_Norm_all_read_pairs.bam | perl -ne 'chomp; print "samtools view -H $_ | grep -H --label=$_ \@RG\n"' | bash > chr6_Exome_Norm_readgroup_info.txt
cat chr6_Exome_Norm_readgroup_info.txt | perl -ne 'my ($bam, $id, $pl, $pu, $lb, $sm); if ($_=~/(\S+\.bam)\:/){$bam=$1} if ($_=~/(ID\:\d+)/){$id=$1} if ($_=~/(PL\:\w+)/){$pl=$1} if ($_=~/(PU\:\S+)/){$pu
=$1} if($_=~/LB\:\"(.+)\"/){$lb=$1} if ($_=~/(SM\:\S+)/){$sm=$1} print "$bam\t$id\t$pl\t$pu\t$lb\t$sm\n";' > chr6_Exome_Norm_readgroup_info.clean.txt

#reverting bams
mkdir -p reverted_bams
mkdir -p reverted_bams/Exome_Norm reverted_bams/Exome_Tumor

java -Xmx8g -jar /usr/local/bin/picard.jar RevertSam I=chr6_Exome_Tumor_all_read_pairs.bam OUTPUT_BY_READGROUP=true O=reverted_bams/Exome_Tumor/
java -Xmx8g -jar /usr/local/bin/picard.jar RevertSam I=chr6_Exome_Norm_all_read_pairs.bam OUTPUT_BY_READGROUP=true O=reverted_bams/Exome_Norm/

#Bam to Fastq:
mkdir -p fastqs
mkdir -p fastqs/Exome_Norm fastq/Exome_Tumor

java -Xmx8g -jar /usr/local/bin/picard.jar SamToFastq I=reverted_bams/Exome_Norm/2891351068.bam F=fastqs/Exome_Norm/2891351068_1.fastq F2=fastqs/Exome_Norm/2891351068_2.fastq
- java -Xmx8g -jar /data/bin/picard.jar SamToFastq I=reverted_bams/Exome_Tumor/2891351066.bam F=fastqs/Exome_Tumor/2891351066_1.fastq F2=fastqs/Exome_Tumor/2891351066_2.fastq

# You may want to move all the above data generated into a designated space for that specific chromosome e.g. mkdir chr6_subset_data

#After moving the data you would need to gzip files and tar ball them for uploading
```
