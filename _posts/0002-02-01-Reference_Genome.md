---
feature_text: |
  ## Precision Medicine
title: Reference Genome
categories:
    - Module-02-Inputs
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-02-01
---

First let’s go over what a reference assembly actually is. In essence, a reference assembly is an attempt at a complete representation of the nucleotide sequence of an individual genome. Individual reads are *assembled* together to form contigs, minimizing gaps, for each chromosome of the species of interest. This reference assembly allows for a shortcut when sequencing future samples/individuals as they can be mapped to the reference, instead of building a new assembly. This has a number of benefits, the most obvious of which is that it is far more effecient than attempting to build a genome from scratch. However, there is no perfect reference assembly for an individual due to polymorphism (i.e., snps, hla-type, etc.). Further, due to the presence of repetitive elements and structural elements such as duplications, inverted repeats, tandem repeats, etc. a given assembly is almost always incomplete, and can always be improved upon. This leads to the publication of new assembly versions every so often such as GRCh37 (Feb. 2009) and GRCh38 (Dec. 2013) for the human reference genome. It is also good to be aware that different organizations can publish different reference assemblies, for example GRCh37 (NCBI) and hg19 (UCSC) are identical save for a few minor differences such as in the mitochondria sequence and naming of chromosomes (1 vs chr1). For a nice summary of genome versions and their release names refer to the [Assembly Releases and Versions FAQ](http://genome.ucsc.edu/FAQ/FAQreleases.html).

### Obtain a reference genome

We will use the 1000 genomes version of the human GRCh38 build. This reference includes extra decoy and HLA sequences in addition to the alternate haplotypes provided from the GRC consortium. The 1000 genomes project is one of several places that people routinely obtain human reference genome files. Some additional sources including those that host many non-human reference genomes are described later in this section.

We obtained the original reference genome files from the 1000 genomes FTP site here:
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/

We have created a copy of these files on our course file server.  Furthermore, we have created a smaller version of the reference to allow us to complete this analysis more quickly.  Using the whole reference genome would take too long for a workshop setting.  For example, aligning reads from a single lange of whole genome data to the whole reference genome can take several hours.

For this course we have selected two chromosomes: chr6 and chr17.  We chose these two chromosomes to illustrate fundamentals of bioinformatics analysis efficiently but also because of the significance of these two chromosomes to cancer biology.  Why are chr6 and chr17 particularly relevant to cancer?

Download the genome reference files for this course using the following commands. Note use of an environment variable `CHRS` to specify the custom reference genome we are using here.

```bash
# make sure CHRS environment variable is set.  If this command doesn't give a value, please return to the Environment section of the course
echo $CHRS

# create a directory for reference genome files and enter this dir
mkdir -p /workspace/inputs/references/genome
cd /workspace/inputs/references/genome

# dowload human reference genome files from the course data server
wget http://genomedata.org/pmbio-workshop/references/genome/$CHRS/ref_genome.tar

# unpack the archive using `tar -xvf` (`x` for extract, `v` for verbose, `f` for file) 
tar -xvf ref_genome.tar

# view contents
tree

# remove the archive
rm -f ref_genome.tar

# uncompress the reference genome FASTA file
gunzip ref_genome.fa.gz

# view contents
tree

```

### Explore the contents of the reference genome file
```bash
cd /workspace/inputs/references/genome

# View the first 10 lines of this file. Note the header line starting with `>`. Why does the sequence look like this?
head ref_genome.fa

# Pull out only the header lines
grep ">" ref_genome.fa

# How many lines and characters are in this file? 
wc ref_genome.fa

# How long are to two chromosomes combined (in bases and Mbp)? Use grep to skip the header lines for each chromosome.
grep -v ">" ref_genome.fa | wc 

# How long does that command take to run?
time grep -v ">" ref_genome.fa | wc

# View 10 lines from approximately the middle of this file
head -n 2500000 ref_genome.fa | tail

# What is the count of each base in the entire reference genome file (skipping the header lines for each sequence)?
cat ref_genome.fa | grep -v ">" | perl -ne 'chomp $_; $bases{$_}++ for split //; if (eof){print "$_ $bases{$_}\n" for sort keys %bases}'

# What does each of these bases refer to? What are the "unexpected bases"?

```

### EXERCISE
Use a commandline scripting approach of your choice to further examine our reference genome file and answer the following question. How many occurences of the EcoRI restriction site are present in the sequence?

{% include question.html question="Solution" answer='EcoRI site (GAATTC) count = 71525' %}




### Learn how to create our own Fasta Index (.fai) files and Dictionary (.dict) files
Index and dictionary files are widely used by other tools to access information in fasta files more efficiently (i.e. faster). These files were included with our reference files (sometimes the case) but it is useful to know how to generate these yourself. You may work with a custom reference in the future where you are required to create such "helper" files.

```bash
# first remove the .fai and .dict files that were downloaded. Do not remove the .fa file though!
cd /workspace/inputs/references/genome
rm -f ref_genome.fa.fai ref_genome.dict

# use samtools to create a fasta index file
samtools faidx ref_genome.fa

# view the contents of the index file
head ref_genome.fa.fai

# use picard to create a dictionary file
java -jar $PICARD CreateSequenceDictionary R=ref_genome.fa O=ref_genome.dict

# view the content of the dictionary file
cat ref_genome.dict

```

### EXERCISE
Figure out what the contents of the fasta index file refer to ...

### Obtain Additional GATK resources needed

```bash
mkdir -p /workspace/inputs/references/gatk
cd /workspace/inputs/references/gatk

#SNP calibration call sets - dbsnp, hapmap, omni, and 1000G
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf /workspace/inputs/references/gatk/
bgzip /workspace/inputs/references/gatk/Homo_sapiens_assembly38.dbsnp138.vcf
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz /workspace/inputs/references/gatk/
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz /workspace/inputs/references/gatk/
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz /workspace/inputs/references/gatk/

#Indel calibration call sets - dbsnp, Mills
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz /workspace/inputs/references/gatk/
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz /workspace/inputs/references/gatk/

#Interval lists
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list /workspace/inputs/references/gatk/
gsutil cp -r gs://genomics-public-data/resources/broad/hg38/v0/scattered_calling_intervals/ /workspace/inputs/references/gatk/

```

### Reference Genome Options

To do: Create a table documenting key reference genome options (builds/sources) and pros/cons
- e.g., 1000G, ensembl, UCSC, GDC

### EXERCISE ANSWERS
How many occurences of the EcoRI restriction site are present in the chromosome 22 sequence? The EcoRI restriction enzyme recognition sequence is 5'-GAATTC-'3. Since this is a palendrome, the reverse complement is the same and we only have to search for one sequence in our string. After accounting for end of line breaks and case sensitivity we find 71525 occurences of this sequence.

```bash
# example code
cd /workspace/inputs/references/genome/
cat ref_genome.fa | grep -v ">" | perl -ne 'chomp $_; $s = uc($_); print $_;' | perl -ne '$c += $_ =~ s/GAATTC/XXXXXX/g; if (eof){print "\nEcoRI site (GAATTC) count = $c\n\n";}'
```

