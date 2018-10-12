---
feature_text: |
  ## Precision Medicine
title: Somatic CNV Calling
categories:
    - Module-04-Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0004-04-01
---

[Copy number alterations (CNA)](https://en.wikipedia.org/wiki/Copy-number_variation) occur when sections of a genome are duplicated or deleted. This phenomenom can actually be quite usefull from an evolutionary standpoint, an example would be the duplication of opsin genes allowing some vertebrate species to see more colors. These types of events however can have a significant impact in the context of disease with perhaps the most famous being an amplification of chromosome 21 resulting in down sydrome. In this section we will go over identifying these types of alterations with [copycat](https://github.com/chrisamiller/copyCat), [CNVnator](https://github.com/abyzovlab/CNVnator), and [cnvkit](https://github.com/etal/cnvkit). However first let's examine exactly what we mean when the segment of a genome is duplicated or deleted and how these types of events can be identified in sequencing data.

{% include figure.html image="/assets/module_4/CNA_illustration.png" position="left" width="450" %}

To begin, the concept of a CNA is fairly straight forward, in the figure to the left we show a standard pair of chromosomes divided into 8 segments. In a copy number amplification a region of the genome is duplicated, in our figure region J has 2 extra segments right after each other. When looking at this region in the reference genome we would expect to see a sharp and dramatic increase in coverage at this region. Further with paired end data the insert size for reads spanning the amplification would be larger. In contrast a deletion is how it sounds and is simply a segment of the chromosome which is gone. In our figure we show a single copy deletion of segment C which would result in a sharp drop in coverage at that region in the sequencing data and a shorter average insert size at the breakpoints of the event.

### Somatic CNA for WGS
[copyCat](https://github.com/chrisamiller/copyCat) is an R package used for detecting somatic (experiment/control) copy number aberrations. It works by measuring the depth of coverage from a sequencing experiment. For example in a diploid organism such as human, a single copy number deletion should result in aproximately half the depth (number of reads) compared to the control. There are two obvious biases when using depth to determine copy nubmer. The first GC-content is well known to have an effect on coverage in illumina sequencing data. This goes all the way down to the PCR amplification of the library. The second thing which can introduce bias is the overall mapability of a region of the genome. A low complexity region in the genome for example will have coverage differences just because it is hard for alignment algorithms to align there. These biases are fairly straight forward to address in WGS data however in exome/targeted data it is much harder.

To begin with [copyCat](https://github.com/chrisamiller/copyCat) we will first need to obtain GC and mapability scores for our data to allow copyCat to adujust for these biases. Fortunately the author of copyCat has a script to create these annotation files. As our first step let's go ahead and download the script into the `workspace/bin` directory and extract the contents.

```bash
cd /workspace/bin
wget https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/copyCat/createCustomAnnotations.v1.zip
unzip createCustomAnnotations.v1.tar.gz
```

The script expects our fasta file to be split by chromosome, we can achieve this with the [faSplit](https://bioconda.github.io/recipes/ucsc-fasplit/README.html) utility. Go ahead and make a new directory called `/workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla_split` to store the result of the split. We then run [faSplit](https://bioconda.github.io/recipes/ucsc-fasplit/README.html) and give it the following positional parameters:

1. byname: tells the program to split the fasta by each record name (i.e. chromosome)
2. GRCh38_full_analysis_set_plus_decoy_hla.fa: location of the multi-record fasta
3. GRCh38_full_analysis_set_plus_decoy_hla_split/: directory to output the results

We also need to make sure we create an index for our new fasta with `samtools faidx`.
```bash
# make new directory and change directories
mkdir -p /workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla_split
cd /workspace/data/raw_data/references/

# split the long fasta by chromosome
faSplit byname GRCh38_full_analysis_set_plus_decoy_hla.fa GRCh38_full_analysis_set_plus_decoy_hla_split/

# index the fasta
samtools faidx /workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla_split/chr6.fa
```

Now that we've got everything set up we can run the script `runEachChr.sh` to create these mapability and gc content annotations, but first lets talk about what the script is actually doing. To create the mapability annotations it first takes the a single record fasta (i.e. for a single chromosome) and generates all possible combination of reads for a given read length (in our case the read length is 100 bp). It then takes these reads and aligns them to the multi record fasta (i.e. the whole genome) to determine which of these reads map uniquely to where they belong. If a read were to be mapped anywhere else from where it originally came from in the single record fasta that would be indicitive of a lower mapability at the region where that read was originally derived. The gc content annotation is the proportion of bases which are either a guanine or cytosine for a given region. The instructions for running this script can be found in the README however to sumarize the script will take the following positional arguments:

1. single record fasta (i.e. chr6)
2. multi record fasta (i.e. full genome)
3. average length of reads (for us this is 100)
4. entrypoint file (chromosome boundaries)
5. output directory

The script will create a directory called `copyCat_annotation` with annotations formatted and structured in the appropriate way. One thing to note is that while the script copied our entrypoints file over we only want to run copyCat on chromosome 6 and so we will need to edit our entrypoint file to reflect this. This script takes a couple hours to run so we provide the results which you can download from [genomedata.org](http://genomedata.org/). We also will need to add a `gaps.bed` file specifying the coordinates for regions to ignore (i.e. telomeres, etc.). A version of this is also made available by the author of copyCat and so will download this file and modify it slightly to correspond to only chromosme 6 and to say chr6 instead of 6. This file should be place at the top level of the copycat annotation directory.

```bash
# run the script to create the annotation dir
#bash runEachChr.sh /workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla_split/chr6.fa /workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa 100 hg38entrypoints.male /workspace/data/results/somatic/

# download the gaps.bed file
cd /workspace/data/results/somatic/copyCat_annotation
wget https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/copyCat/GRCh38/gaps.bed
cat gaps.bed | grep "^6" | awk '{print "chr"$0}' > tmp && mv tmp gaps.bed

# edit entrypoint file to contain only chromosome 6
grep "chr6" entrypoints.male > tmp && mv tmp entrypoints.male
```

The final step of the puzzle before we run [copyCat](https://github.com/chrisamiller/copyCat) is to obtain depth calculations corresponding to a specific window size for our sequencing data. We can use [mosdepth](https://academic.oup.com/bioinformatics/article/34/5/867/4583630) for this task. The paramters we give to mosdepth are described below, we also need to decompress the output from mosdepth and manipulate the output such that it has three columns with header names "Chr", "Start", "Counts.$average_read_size" (i.e. Counts.100).

1. --no-per-base: omit the per base depth calculation and only calculate depth per region (makes mosdepth run much faster)
2. -t: number of threads to use
3. -b: window size for calculating depth
4. output directory
5. bam file

```bash
# run mosdepth for tumor/normal
mosdepth --no-per-base -t 10 -b 1000 /workspace/data/results/somatic/WGS_Norm.mosdepth /workspace/data/results/align/WGS_Norm_merged_sorted_mrkdup.bam
bgzip -d WGS_Norm.mosdepth.regions.bed.gz
cat /workspace/data/results/somatic/WGS_Norm.mosdepth.regions.bed | cut -f 1,2,4 | awk 'BEGIN{print "Chr\tStart\tCounts.100"}1' > /workspace/data/results/somatic/WGS_Norm.mosdepth.regions.2.bed

mosdepth --no-per-base -t 10 -b 1000 /workspace/data/results/somatic/WGS_Tumor.mosdepth /workspace/data/results/align/WGS_Tumor_merged_sorted_mrkdup.bam
bgzip -d WGS_Tumor.mosdepth.regions.bed.gz
cat /workspace/data/results/somatic/WGS_Tumor.mosdepth.regions.bed | cut -f 1,2,4 | awk 'BEGIN{print "Chr\tStart\tCounts.100"}1' > /workspace/data/results/somatic/WGS_Tumor.mosdepth.regions.2.bed
```
With Everything now set up we can start `R` and load the copyCat library. From there we can run the `runPairedSampleAnalysis()` function to perform the analysis. Most of the parameters in the function are the defaults and are only provided for the sake of completeness, the parameters changed are as follows:

1. annotationDirectory: the path to the mapability and gc annotations we created from the bash script
2. outputDir: whre to write our output
3. normal: path to the normal mosdepth based depth calculations
4. tumor: path to the tumor mosdepth based depth calculations
5. maxCores: number of cores to use

```R
# Start R
R

# load the copyCat library
library(copyCat)

# run copyCat in paired mode
runPairedSampleAnalysis(annotationDirectory="/workspace/data/results/somatic/copyCat_annotation",
                        outputDirectory="/workspace/data/results/somatic/",
                        normal="/workspace/data/results/somatic/WGS_Norm.mosdepth.regions.2.bed",
                        tumor="/workspace/data/results/somatic/WGS_Tumor.mosdepth.regions.2.bed",
                        inputType="bins",
                        maxCores=4,
                        binSize=0,
                        perLibrary=1,
                        perReadLength=1,
                        verbose=TRUE,
                        minWidth=3,
                        minMapability=0.6,
                        dumpBins=TRUE,
                        doGcCorrection=TRUE,
                        samtoolsFileFormat="unknown",
                        purity=1,
                        normalSamtoolsFile=NULL,
                        tumorSamtoolsFile=NULL)
```

### cnvnator germline
### cnvkit exome
