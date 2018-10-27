---
feature_text: |
  ## Precision Medicine
title: Annotation
categories:
    - Module-02-Inputs
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-03-01
---

### Obtain Additional GATK resource files needed
Use Google’s gsutil to download various annotation files that will be used by GATK and other resources. gsutil will be used to download these file from Google cloud storage. They could also be downloaded using wget from the course file server from here: http://genomedata.org/pmbio-workshop/references/gatk/.

```bash
cd /workspace/inputs/references/
mkdir -p gatk
cd gatk

# SNP calibration call sets - dbsnp, hapmap, omni, and 1000G
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf . 
gzip Homo_sapiens_assembly38.dbsnp138.vcf
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz .
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz .
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz .

# Indel calibration call sets - dbsnp, Mills
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz .
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz .

# Interval lists that can be used to parallelize certain GATK tasks
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list .
gsutil cp -r gs://genomics-public-data/resources/broad/hg38/v0/scattered_calling_intervals/ .

# list the files we just downloaded
ls -lh

```

### Download coordinates describing the Exome Reagent used to generate our exome data (SeqCapEZ_Exome_v3.0)
The reagent used to produce the exome sequenced data for this course was the SeqCapEZ_Exome_v3.0 from Roche Nimblegen. Nimblegen provides a set of files describing this reagent which can be downloaded from their site
[here](https://sequencing.roche.com/en/products-solutions/by-category/target-enrichment/hybridization/seqcap-ez-exome-v3-kit.html). We are most interested in the file named `SeqCap_EZ_Exome_v3_hg19_primary_targets.bed` which contains the chromosome, start, stop, and gene annotation for each probe used in the reagent. Let's go ahead and downlod these files for later use.

```bash
# change directories
mkdir -p /workspace/inputs/references/exome
cd /workspace/inputs/references/exome

# download the files
wget -c https://sequencing.roche.com/content/dam/rochesequence/worldwide/resources/SeqCapEZ_Exome_v3.0_Design_Annotation_files.zip
unzip SeqCapEZ_Exome_v3.0_Design_Annotation_files.zip

# remove the zip
rm -f SeqCapEZ_Exome_v3.0_Design_Annotation_files.zip
```

### Convert SeqCapEZ_Exome_v3.0

You might have noticed that these annotation files from Nimblegen are all from the hg19 genome assembly. Obviously this presents a problem as the analysis we're performing is using the newer hg38 assembly. This issue is actually not an uncommon situation and a variety of tools exist that are designed to make converting from one assembly to another easier. In the next section we will be using [UCSC liftover](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to perform this task.

To start we first need to download a chain file specific to the assembly conversion we want to perform (in our case hg19 -> hg38). These files provide a mapping between the two assemblies and can be downloaded from the [UCSC site](http://hgdownload.cse.ucsc.edu/downloads.html#liftover). Once we have our chain file we can run `liftOver` which will take following positional arguments.

1. Path to bed file to convert
2. Path to chain file for the desired conversion (downloaded from UCSC)
3. Path to output file that will contain the new coordinates
4. Path to output file containing any coordinates that do not convert successfully

It will also be usefull to have a version of this file without the gene annotations lets go ahead and create that here as well.

```bash
# change to the appropriate directory
cd /workspace/inputs/references/exome

# download the chain file
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# run liftover
liftOver SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_primary_targets.bed  hg19ToHg38.over.chain.gz SeqCap_EZ_Exome_v3_hg38_primary_targets.bed unMapped.bed

# create a version in standard bed format (chr, start, stop)
cut -f 1-3 SeqCap_EZ_Exome_v3_hg38_primary_targets.bed > SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed

# take a quick look at the format of this file
head SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed

# a very common question in exome (or any capture based approach) is: how big is my capture space?
# use bedtools to determine the size of the capture space represented by this bed file

# first sort the bedfile and store the sorted version
bedtools sort -i SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed > SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.bed

# now merge the bed file to collapse any overlapping regions so they are not double counted.
bedtools merge -i SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.bed > SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.bed

# finally use a Perl one liner to determine the size of each non-overlapping region and determine the cumulative sum
perl -ne 'chomp; @l=split("\t",$_); $size += $l[2]-$l[1]; print "$size\n"' SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.bed

# note that the size of the space targeted by the exome reagent is ~63 Mbp. Does that sound reasonable?

# now create a subset of this bed file
grep -w -P "^chr6|^chr17" SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.bed > exome_regions.bed

# clean up intermediate files
rm -fr SeqCapEZ_Exome_v3.0_Design_Annotation_files/ SeqCap_EZ_Exome_v3_hg38_primary_targets.bed SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.bed unMapped.bed

```

### Obtaining an interval list for the exome bed file

NOTE: The following interval list is being created for an annotation file the corresponds to the whole genome.  The `ref_genome.dict` we downloaded in the previous section only contains chr6 and chr17. We need to either trim the exome ROI file down to just these chrs as well. Or also download/create a full version of the dict file.

```bash
# first for the complete exome bed file
cd /workspace/inputs/references/
mkdir temp
cd temp
wget http://genomedata.org/pmbio-workshop/references/genome/all/ref_genome.dict
cd /workspace/inputs/references/exome
java -jar /usr/local/bin/picard.jar BedToIntervalList I=SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.bed O=SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.interval_list SD=/workspace/inputs/references/temp/ref_genome.dict
rm -fr /workspace/inputs/references/temp/

# next for our subset exome regions file
cd /workspace/inputs/references/exome
java -jar /usr/local/bin/picard.jar BedToIntervalList I=exome_regions.bed O=exome_regions.bed.interval_list SD=/workspace/inputs/references/genome/ref_genome.dict

```

### Obtaining transcriptome reference files
In this section we will download some reference transcriptome annotations in the mostly widely used formats (.fa and .gtf). As with the other annotation files, we have created a subsetted version of them to make downstream analysis more practical for a live tutorial demonstration.  To create these annotation we followed these basic steps

* Download complete GTF files from Ensembl represent all gene/transcript annotations (e.g. `Homo_sapiens.GRCh38.94.gtf.gz`) from Ensembl's [FTP site](http://www.ensembl.org/info/data/ftp/index.html/).
* Fix the chromosome names in this GTF. Remember that Ensembl uses names like `1`, `2`, etc. but our reference genome uses names like `chr1`, `chr2`, etc. We perform this conversion using chromosome synonym mappings from Ensembl and a simple script `convertEnsemblGTF.pl`. 
* Next we use a simple grep to pull out the subset of chromosomes we care about
* Finally we use the resulting subsetted GTF to create FASTA sequences for each transcript using the `gtf_to_fasta` tool.  This tools takes exon coordinates in a GTF file and a reference genome sequence and creates the reference transcript sequences.

Complete details of how these files were created can be found in the [Developer Notes](/module-10-appendix/0010/02/01/Developer_Notes/).

Now lets actually download these files and examine them...
```bash

...


```



