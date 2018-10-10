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
wget https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/createCustomAnnotations.v1.tar.gz
tar -xzvf createCustomAnnotations.v1.tar.gz
bash runEachChr.sh chr6 100 /workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla_split /workspace/data/results/somatic/
```

Next the script expects our fasta file to be split by chromosome, we can achieve this with the [faSplit](https://bioconda.github.io/recipes/ucsc-fasplit/README.html) utility. Go ahead and make a new directory called `/workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla_split` to store the result. We then run [faSplit](https://bioconda.github.io/recipes/ucsc-fasplit/README.html) and give it the following positional parameters:

1. byname: tells the program to split the fasta by each record name
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



```bash
# run mosdepth for tumor/normal
mosdepth -t 10 -b 100 /workspace/data/results/somatic/WGS_Norm.mosdepth /workspace/data/results/align/WGS_Norm_merged_sorted_mrkdup.bam
mosdepth -t 10 -b 100 /workspace/data/results/somatic/WGS_Tumor.mosdepth /workspace/data/results/align/WGS_Tumor_merged_sorted_mrkdup.bam
```

### cnvnator germline
### cnvkit exome
