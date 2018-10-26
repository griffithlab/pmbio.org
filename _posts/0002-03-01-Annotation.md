---
feature_text: |
  ## Precision Medicine
title: Annotation
categories:
    - Module-02-Inputs
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-03-01
---

### Obtain Additional GATK resources needed
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

gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list .
gsutil cp -r gs://genomics-public-data/resources/broad/hg38/v0/scattered_calling_intervals/ .

# list the files we just downloaded
ls -lh

```

### SeqCapEZ_Exome_v3.0
The reagent used for exome sequencing of the data used in this course was the SeqCapEZ_Exome_v3.0 from Roche Nimblegen. Nimblegen provides a set of files describing this reagent which can be downloaded from their site
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

You might have noticed that these annotation files from Nimblegen are all from the hg19 genome assembly. Obviously this presents a problem as the analysis we're performing is using the newer hg38 assembly. This issue is actually not an uncommon situation and a variety of tools exist designed to make converting from one assembly to another easier. In the next section we will be using [UCSC liftover](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to perform this task.

To start we first need to download a chain file specific to the assembly conversion we want to perform (in our case hg19 -> hg38). These files provide a mapping between the two assemblies and can be downloaded from the [UCSC site](http://hgdownload.cse.ucsc.edu/downloads.html#liftover). Once we have our chain file we can run `liftOver` which will take following positional arguments.

1. Path to bed file to convert
2. Path to chain file downloaded from UCSC
3. Path to output file containing the new coordinates
4. Path to output file containing any coordinates that did not convert successfully

It will also be usefull to have a version of this file without the gene annotations lets go ahead and do that here has well.

```bash
# change to the appropriate directory
cd /workspace/inputs/references/exome

# download the chain file
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# run liftover
liftOver SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_primary_targets.bed  hg19ToHg38.over.chain.gz SeqCap_EZ_Exome_v3_hg38_primary_targets.bed unMapped.bed

# create a version in standard bed format (chr, start, stop)
cut -f 1-3 SeqCap_EZ_Exome_v3_hg38_primary_targets.bed > SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed
```

### Creating an interval list from bed file

NOTE: The following is broken because our ref_genome.dict now only contains chr6 and chr17. We need to either trim the exome ROI file down to just these chrs as well. Or also download/create a full version of the dict file.

```bash
cd /workspace/inputs/references/exome
java -jar $PICARD BedToIntervalList I=SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed O=SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.interval_list SD=/workspace/inputs/references/genome/ref_genome.dict
```
