---
feature_text: |
  ## Precision Medicine
title: Somatic CNV Calling
categories:
    - Module-05-Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-04-01
---

[Copy number alterations (CNA)](https://en.wikipedia.org/wiki/Copy-number_variation) occur when sections of a genome are duplicated or deleted. This phenomenom can actually be quite usefull from an evolutionary standpoint, an example would be the duplication of opsin genes allowing some vertebrate species to see more colors. These types of events however can have a significant impact in the context of disease with perhaps the most famous being an amplification of chromosome 21 resulting in down sydrome. In this section we will go over identifying these types of alterations with [copycat](https://github.com/chrisamiller/copyCat), [CNVnator](https://github.com/abyzovlab/CNVnator), and [cnvkit](https://github.com/etal/cnvkit). However first let's examine exactly what we mean when the segment of a genome is duplicated or deleted and how these types of events can be identified in sequencing data.

{% include figure.html image="/assets/module_4/CNA_illustration.png" position="left" width="450" %}

To begin, the concept of a CNA is fairly straight forward, in the figure to the left we show a standard pair of chromosomes divided into 8 segments. In a copy number amplification a region of the genome is duplicated, in our figure region J has 2 extra segments right after each other. When looking at this region in the reference genome we would expect to see a sharp and dramatic increase in coverage at this region. Further with paired end data the insert size for reads spanning the amplification would be larger. In contrast a deletion is how it sounds and is simply a segment of the chromosome which is gone. In our figure we show a single copy deletion of segment C which would result in a sharp drop in coverage at that region in the sequencing data and a shorter average insert size at the breakpoints of the event.

### Somatic CNA for WGS

[copyCat](https://github.com/chrisamiller/copyCat) is an R package used for detecting somatic (experiment/control) copy number aberrations. It works by measuring the depth of coverage from a sequencing experiment. For example in a diploid organism such as human, a single copy number deletion should result in aproximately half the depth (number of reads) compared to the control. Before we run [copyCat](https://github.com/chrisamiller/copyCat) we need to obtain depth calculations corresponding to a specific window size for our sequencing data. We can use [mosdepth](https://academic.oup.com/bioinformatics/article/34/5/867/4583630) for this task. The paramters we give to mosdepth are described below, to run with [copyCat](https://github.com/chrisamiller/copyCat) we also need to decompress the output from mosdepth and manipulate the output such that it has three columns with header names "Chr", "Start", "Counts.$average_read_size" (i.e. Counts.100).

**mosdepth paramters**:
1. --no-per-base: omit the per base depth calculation and only calculate depth per region (makes mosdepth run much faster)
2. -t: number of threads to use
3. -b: window size for calculating depth
4. output directory
5. bam file

```bash
# run mosdepth for tumor/normal
mosdepth --no-per-base -t 4 -b 10000 /workspace/data/results/somatic/WGS_Norm.mosdepth /workspace/data/results/align/WGS_Norm_merged_sorted_mrkdup.bam
bgzip -d WGS_Norm.mosdepth.regions.bed.gz
cat /workspace/data/results/somatic/WGS_Norm.mosdepth.regions.bed | cut -f 1,2,4 | awk 'BEGIN{print "Chr\tStart\tCounts.100"}1' > /workspace/data/results/somatic/WGS_Norm.mosdepth.regions.2.bed

mosdepth --no-per-base -t 4 -b 10000 /workspace/data/results/somatic/WGS_Tumor.mosdepth /workspace/data/results/align/WGS_Tumor_merged_sorted_mrkdup.bam
bgzip -d WGS_Tumor.mosdepth.regions.bed.gz
cat /workspace/data/results/somatic/WGS_Tumor.mosdepth.regions.bed | cut -f 1,2,4 | awk 'BEGIN{print "Chr\tStart\tCounts.100"}1' > /workspace/data/results/somatic/WGS_Tumor.mosdepth.regions.2.bed
```

Having run [mosdepth](https://academic.oup.com/bioinformatics/article/34/5/867/4583630) we have the core data we need to determine copy number. Let's go ahead and quickly examine a chromosome to get an idea of what we're looking at. The below code will make a crude adjustment to normalize the data based on the number of reads present in the samples and plot the result for chromosome 6.

```R
# start R
R

# set working directory
setwd("/workspace/data/results/somatic/")

# load libraries
library(ggplot2)
library(data.table)
library(viridis)

# read in the depth depth data
tumor_depth <- fread("WGS_Tumor.mosdepth.regions.2.bed")
normal_depth <- fread("WGS_Norm.mosdepth.regions.2.bed")

# merge the two
all_depth <- merge(tumor_depth, normal_depth, by=c("Chr", "Start"), suffixes=c("Tumor", "Normal"))

# normalize the data based on number of reads
total_tumor_reads <- sum(all_depth$Counts.100Tumor)
total_normal_reads <- sum(all_depth$Counts.100Normal)
all_depth$Counts.100Normal <- all_depth$Counts.100Normal * (total_tumor_reads/total_normal_reads)

# calculate the difference in depth between tumor and normal
all_depth$diff <- all_depth$Counts.100Tumor - all_depth$Counts.100Normal

# plot the result for chromosome 6
pdf(file="tumor_depth_cn.chr1.pdf", height=5, width=10)
ggplot(all_depth[all_depth$Chr == "chr6",], aes(x=Start, y=diff, color=diff)) + geom_point() + scale_color_viridis("Depth", option = "plasma") + theme_bw() +
    ylab("Relative Depth Difference") + xlab("Position")
dev.off()
```

you should see something like the plot below where there is a fairly clear indication of a copy number amplification on the p arm of the chromosome, and a copy deletion twoards the center of the chromosome. Keep in mind that what's plotted is the Tumor depth relative to the normal and not the actual copy number.

{% include figure.html image="/assets//module_4/mosdepth_CNA.png" %}


Nice job, that was pretty easy! But not so fast, there are two obvious biases when using depth to determine copy number. The first, GC-content is well known to have an effect on coverage in illumina sequencing data. This goes all the way back to the PCR amplification of the library affecting the number of reads generated. The second thing which can introduce bias is the overall mapability of a region of the genome. A low complexity region in the genome for example will have coverage differences just because it is hard for alignment algorithms to map reads there. These biases are fairly straight forward to address in WGS data however in exome/targeted data it is much more complicated.

The [copyCat](https://github.com/chrisamiller/copyCat) package is able to adjust for these biases in WGS data and is what we'll start with. To start with [copyCat](https://github.com/chrisamiller/copyCat) we need to obtain GC content and mapability scores for our data. Fortunately the author of copyCat has a script to create these annotation files. As our first step let's go ahead and download the script into the `workspace/bin` directory and extract the contents.

```bash
cd /workspace/bin
wget http://genomedata.org/pmbio-workshop/misc/createCustomAnnotations.v1.zip
unzip createCustomAnnotations.v1.tar.gz
```

The script expects our genome fasta file to be split by chromosome, we can achieve this with the [faSplit](https://bioconda.github.io/recipes/ucsc-fasplit/README.html) utility. Go ahead and make a new directory called `/workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla_split` to store the result of the split. We then run [faSplit](https://bioconda.github.io/recipes/ucsc-fasplit/README.html) and give it the following positional parameters:

- byname: tells the program to split the fasta by each record name (i.e. chromosome)
- GRCh38_full_analysis_set_plus_decoy_hla.fa: location of the multi-record fasta
- GRCh38_full_analysis_set_plus_decoy_hla_split/: directory to output the results

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

The script will create a directory called `copyCat_annotation` with annotations formatted and structured in the appropriate way. This script takes a couple hours to run so we provide the results which you can download from genomedata.org. We also will need to add a `gaps.bed` file specifying the coordinates for regions to ignore (i.e. telomeres, etc.). A version of this is also made available by the author of copyCat and so we will download this file as well and modify it slightly such that our chromosome names begin with **chr**. This file should be place at the top level of the copycat annotation directory.

```bash
## run the script to create the annotation dir
# cd /workspace/bin/createCustomAnnotations.v1
# bash runEachChr.sh /workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla_split/chr6.fa /workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa 100 hg38entrypoints.female /workspace/data/results/somatic/

# as mentioned the above takes some time so we'll just download the result for HG38
cd /workspace/data/results/somatic/
wget http://genomedata.org/pmbio-workshop/misc/copyCat_annotation.zip
unzip copyCat_annotation.zip

## download the gaps.bed file
# cd /workspace/data/results/somatic/copyCat_annotation
# https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/copyCat/GRCh38/gaps.bed
# cat gaps.bed | awk '{print "chr"$0}' > tmp && mv tmp gaps.bed
```

At the end your directory structure should look something like this but with more chromosomes:

```bash
├── entrypoints.female
├── gaps.bed
└── readlength.100
    ├── gcWinds
    │   ├── chr10.gc.gz
    │   ├── chr11.gc.gz
    └── mapability
        ├── chr10.dat.gz
        ├── chr11.dat.gz
```

With Everything now set up we can start `R` and load the copyCat library. From there we can run the `runPairedSampleAnalysis()` function to perform the analysis. Most of the parameters in the function are the defaults and are only provided for the sake of completeness, the parameters changed are as follows:

1. annotationDirectory: the path to the mapability and gc annotations we created from the bash script
2. outputDir: whre to write our output
3. normal: path to the normal mosdepth based depth calculations
4. tumor: path to the tumor mosdepth based depth calculations
5. maxCores: number of cores to use, 0 means all available cores

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
                        maxCores=0,
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

The analysis will take a few minutes to complete however once it's done there are a few files that we care about. First you'll notice there is now a plots directory at `/workspace/data/results/somatic/plots`, inside we can view the graphs [copyCat](https://github.com/chrisamiller/copyCat) created to visualize the gc bias correction, they should look something like this:

{% include figure.html image="/assets/module_4/normal.gccontent.lib1.readLength100.png" %}

As we can see there was quite an extreme bias between the number of reads mapped and the GC content of the reads particulary when the GC content falls below 30% or above 50% however the LOESS correction dealt with this quite nicely.

Next let's look at the actual copy number calls from [copyCat](https://github.com/chrisamiller/copyCat). There are three files output files that we care about here:

- **rd.bins.dat**: contains the chromsome, position and the inferred tumor CN where 2 represents a normal diploid state
- **segs.paired.dat**: contains the result of the segmentation algorithm with columns "chromosome", "start", "stop", "windows", "inferred tumor CN"
- **alts.paired.dat**: is identical to segs.paired.dat however those segments which are copy neutral are removed.

As a final step lets load R and plot the results for chromosome 6.

```R
# start R
R

# set working directory
setwd("/workspace/data/results/somatic/")

# load libraries
library(ggplot2)
library(data.table)
library(viridis)
library(scales)

# read in the data
copy_segment_alteration <- fread("alts.paired.dat")
colnames(copy_segment_alteration) <- c("Chr", "Start", "Stop", "Windows", "Tumor_CN")
cna_bin <- fread("rd.bins.dat")
colnames(cna_bin) <- c("Chr", "Pos", "CNA")

# create the plot
pdf(file="copycat_final.chr2.pdf", height=5, width=10)
ggplot() + geom_point(data=cna_bin[cna_bin$Chr == "chr2",], aes(x=Pos, y=CNA, color=CNA)) +
    geom_segment(data=copy_segment_alteration[copy_segment_alteration$Chr == "chr2",], aes(x=Start, xend=Stop, y=Tumor_CN, yend=Tumor_CN), color="black", size=1) +
    scale_y_continuous(limits=c(0, 15), oob=squish) + scale_color_viridis(limits=c(0, 4), option="plasma", oob=squish) + theme_bw() +
    geom_hline(yintercept = c(1, 2, 3), linetype="longdash")
dev.off()
```

{% include figure.html image="/assets/module_4/copyCat_final.png" %}

### cnvkit wgs
```bash
mkdir -p /workspace/data/results/somatic/cnvkit_wgs
cd /workspace/data/results/somatic/cnvkit_wgs

source activate cnvkit

cnvkit.py batch /workspace/data/results/align/WGS_Tumor_merged_sorted_mrkdup.bam --normal /workspace/data/results/align/WGS_Norm_merged_sorted_mrkdup.bam --fasta /workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa --access /workspace/data/raw_data/references/access-excludes.hg38.bed --output-reference /workspace/data/raw_data/references/my_reference.cnn --output-dir /workspace/data/results/somatic/cnvkit_wgs/ --method wgs -p 20 --diagram --scatter

source deactivate

convert WGS_Tumor_merged_sorted_mrkdup-scatter.pdf WGS_Tumor_merged_sorted_mrkdup-scatter.png
convert WGS_Tumor_merged_sorted_mrkdup-scatter.pdf WGS_Tumor_merged_sorted_mrkdup-scatter.jpg

convert WGS_Tumor_merged_sorted_mrkdup-diagram.pdf WGS_Tumor_merged_sorted_mrkdup-diagram.png
convert WGS_Tumor_merged_sorted_mrkdup-diagram.pdf WGS_Tumor_merged_sorted_mrkdup-diagram.jpg
```
```bash
mkdir -p /workspace/data/results/somatic/cnvkit_exome
cd /workspace/data/results/somatic/cnvkit_exome

source activate cnvkit

wget http://18.223.213.22/refseq/hglft_genome_304d_b78af0.bed

cnvkit.py access /workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa -x /workspace/data/results/somatic/copyCat_annotation/gaps.bed -o /workspace/data/raw_data/references/access-excludes.hg38.bed

cnvkit.py batch /workspace/data/results/align/Exome_Tumor_sorted_mrkdup.bam --normal /workspace/data/results/align/Exome_Norm_sorted_mrkdup.bam --targets hglft_genome_304d_b78af0.bed --fasta /workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa --access /workspace/data/raw_data/references/access-excludes.hg38.bed --output-reference /workspace/data/raw_data/references/my_reference.cnn --output-dir /workspace/data/results/somatic/cnvkit_exome/ --method hybrid -p 20 --diagram --scatter

source deactivate

convert Exome_Tumor_sorted_mrkdup-scatter.pdf Exome_Tumor_sorted_mrkdup-scatter.png
convert Exome_Tumor_sorted_mrkdup-scatter.pdf Exome_Tumor_sorted_mrkdup-scatter.jpg

convert Exome_Tumor_sorted_mrkdup-diagram.pdf Exome_Tumor_sorted_mrkdup-diagram.png
convert Exome_Tumor_sorted_mrkdup-diagram.pdf Exome_Tumor_sorted_mrkdup-diagram.jpg
```
