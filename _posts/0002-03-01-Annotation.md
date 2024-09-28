---
feature_text: |
  ## Precision Medicine Bioinformatics
  Introduction to bioinformatics for DNA and RNA sequence analysis
title: Annotation
categories:
    - Module-02-Inputs
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-03-01
---

### Obtain Additional GATK resource files needed
Use Google’s gsutil to download various annotation files that will be used by GATK and other resources. gsutil will be used to download these file from Google cloud storage. 

```bash
cd /workspace/inputs/references/
mkdir -p gatk
cd gatk

# SNP calibration call sets - dbsnp, hapmap, omni, and 1000G
# Runtime: < 2min
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf .
# Runtime: ~ 2min
bgzip --threads 8 Homo_sapiens_assembly38.dbsnp138.vcf
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
[here](https://sequencing.roche.com/en/products-solutions/by-category/target-enrichment/hybridization/seqcap-ez-exome-v3-kit.html). We are most interested in the file named `SeqCap_EZ_Exome_v3_hg19_capture_targets.bed` which contains the chromosome, start, stop, and gene annotation for each probe used in the reagent. Let's go ahead and downlod these files for later use.

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

Here we are processing two of the downloaded files, `SeqCap_EZ_Exome_v3_hg19_primary_targets.bed` and `SeqCap_EZ_Exome_v3_hg19_capture_targets.bed`. The former corresponding to the exome regions that we are interested in and the latter corresponding to the regions that probes were able to be designed on.

```bash
# change to the appropriate directory
cd /workspace/inputs/references/exome

# download the chain file
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# run liftover
liftOver SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_primary_targets.bed  hg19ToHg38.over.chain.gz SeqCap_EZ_Exome_v3_hg38_primary_targets.bed unMapped.bed

liftOver SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed  hg19ToHg38.over.chain.gz SeqCap_EZ_Exome_v3_hg38_capture_targets.bed unMapped1.bed

# create a version in standard bed format (chr, start, stop)
cut -f 1-3 SeqCap_EZ_Exome_v3_hg38_primary_targets.bed > SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed

cut -f 1-3 SeqCap_EZ_Exome_v3_hg38_capture_targets.bed > SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.bed

# take a quick look at the format of these files
head SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed
head SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.bed

# a very common question in exome (or any capture based approach) is: how big is my capture space?
# use bedtools to determine the size of the capture space represented by this bed file

# first sort the bed files and store the sorted versions
bedtools sort -i SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed > SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.bed

bedtools sort -i SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.bed > SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.sort.bed

# now merge the bed files to collapse any overlapping regions so they are not double counted.
bedtools merge -i SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.bed > SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.bed

bedtools merge -i SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.sort.bed > SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.sort.merge.bed

# finally use a Perl one liner to determine the size of each non-overlapping region and determine the cumulative sum
perl -ne 'chomp; @l=split("\t",$_); $size += $l[2]-$l[1]; print "$size\n"' SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.bed
perl -ne 'chomp; @l=split("\t",$_); $size += $l[2]-$l[1]; print "$size\n"' SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.sort.merge.bed

# note that the size of the space targeted by the exome reagent is ~63 Mbp. Does that sound reasonable?

# now create a subset of these bed files
grep -w -P "^chr6|^chr17" SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.bed > exome_regions.bed

grep -w -P "^chr6|^chr17" SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.sort.merge.bed > probe_regions.bed

# clean up intermediate files
#rm -fr SeqCapEZ_Exome_v3.0_Design_Annotation_files/ SeqCap_EZ_Exome_v3_hg38_primary_targets.bed SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.bed unMapped.bed

```

### Obtaining an interval list for the exome bed files

NOTE: The following interval list is being created for an annotation file the corresponds to the whole genome.  The `ref_genome.dict` we downloaded in the previous section only contains chr6 and chr17. We need to either trim the exome ROI file down to just these chrs as well. Or also download/create a full version of the dict file.

```bash
# first for the complete exome and probe bed file
cd /workspace/inputs/references/
mkdir temp
cd temp
wget http://genomedata.org/pmbio-workshop/references/genome/all/ref_genome.dict
cd /workspace/inputs/references/exome
java -jar $PICARD BedToIntervalList I=SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.bed O=SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.interval_list SD=/workspace/inputs/references/temp/ref_genome.dict
java -jar $PICARD BedToIntervalList I=SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.sort.merge.bed O=SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.sort.merge.interval_list SD=/workspace/inputs/references/temp/ref_genome.dict
rm -fr /workspace/inputs/references/temp/

# next for our subset exome and probe regions file
cd /workspace/inputs/references/exome
java -jar /usr/local/bin/picard.jar BedToIntervalList I=exome_regions.bed O=exome_regions.bed.interval_list SD=/workspace/inputs/references/genome/ref_genome.dict
java -jar /usr/local/bin/picard.jar BedToIntervalList I=probe_regions.bed O=probe_regions.bed.interval_list SD=/workspace/inputs/references/genome/ref_genome.dict

```

### Obtaining transcriptome reference files
In this section we will download some reference transcriptome annotations in the mostly widely used formats (.fa and .gtf). As with the other annotation files, we have created a subsetted version of them to make downstream analysis more practical for a live tutorial demonstration.  To create these annotation files we followed these basic steps:

* Download complete GTF files from Ensembl represent all gene/transcript annotations (e.g. `Homo_sapiens.GRCh38.94.gtf.gz`) from Ensembl's [FTP site](http://www.ensembl.org/info/data/ftp/index.html/).
* Fix the chromosome names in this GTF. Remember that Ensembl uses names like `1`, `2`, etc. but our reference genome uses names like `chr1`, `chr2`, etc. We perform this conversion using chromosome synonym mappings from Ensembl and a simple script `convertEnsemblGTF.pl`.
* Next we use a simple grep to pull out the subset of chromosomes we care about
* Finally we use the resulting subsetted GTF to create FASTA sequences for each transcript using the `gtf_to_fasta` tool.  This tools takes exon coordinates in a GTF file and a reference genome sequence and creates the reference transcript sequences.

Complete details of how these files were created can be found in the [Developer Notes](/module-10-appendix/0010/02/01/Developer_Notes/).

Now lets actually download these files and examine them...
```bash
# make sure CHRS environment variable is set.  If this command doesn't give a value, please return to the Environment section of the course
echo $CHRS

# create a directory for transcriptome input files
mkdir -p /workspace/inputs/references/transcriptome
cd /workspace/inputs/references/transcriptome

# download the files
wget http://genomedata.org/pmbio-workshop/references/transcriptome/$CHRS/ref_transcriptome.gtf
wget http://genomedata.org/pmbio-workshop/references/transcriptome/$CHRS/ref_transcriptome.fa

# take a look at the contents of the gtf file. Press 'q' to exit the 'less' display.
less -p start_codon -S ref_transcriptome.gtf

# How many unique gene IDs are in the .gtf file?
# We can use a perl command-line command to find out:
perl -ne 'if ($_ =~ /(gene_id\s\"ENSG\w+\")/){print "$1\n"}' ref_transcriptome.gtf | sort | uniq | wc -l

```

* Using perl -ne '' will execute the code between single quotes, on the .gtf file, line-by-line.
* The $_ variable holds the contents of each line.
* The 'if ($_ =~//)' is a pattern-matching command which will look for the pattern "gene_id" followed by a space followed by "ENSG" and one or more word characters (indicated by \w+) surrounded by double quotes.
* The pattern to be matched is enclosed in parentheses. This allows us to print it out from the special variable $1.
* The output of this perl command will be a long list of ENSG Ids.
* By piping to sort, then uniq, then word count we can count the unique number of genes in the file.

```bash
# what are all the feature types listed in the third column of the GTF?
# how does the following command (3 commands piped together) answer that question?
cut -f 3 ref_transcriptome.gtf | sort | uniq -c

```

#### Review key definitions for this section

* Reference genome - The nucleotide sequence of the chromosomes of a species. Genes are the functional units of a reference genome and gene annotations describe the structure of transcripts expressed from those gene loci.
* Gene annotations - Descriptions of gene/transcript models for a genome. A transcript model consists of the coordinates of the exons of a transcript on a reference genome. Additional information such as the strand the transcript is generated from, gene name, coding portion of the transcript, alternate transcript start sites, and other information may be provided.
* GTF (.gtf) file - A common file format referred to as Gene Transfer Format used to store gene and transcript annotation information. You can learn more about this format here: http://genome.ucsc.edu/FAQ/FAQformat#format4
