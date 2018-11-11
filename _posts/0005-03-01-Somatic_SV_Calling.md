---
feature_text: |
  ## Precision Medicine Bioinformatics
  Introduction to bioinformatics for DNA and RNA sequence analysis
title: Somatic SV Calling
categories:
    - Module-05-Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-03-01
---

[Manta](https://github.com/Illumina/manta) is a structural variant caller maintained by Illumina and optimized for calling somatic variation in tumor/normal pairs. In this section we will use [Manta](https://github.com/Illumina/manta) to call structural variants in our WGS data but first let's go over what a structural variant actually is. Structural variants are rearrangements in DNA involving a breakpoint(s). Generally speaking structural variants can fall into four categories:

1. Insertions: a region is inserted into the DNA
2. Deletions: a region is deleted in the DNA
3. Inversions: a section of DNA is reversed
4. Translocations: a section of DNA is remved and re-inserted in a new region

Running manta is a two step process. The first step is configuring manta to run with `configManta.py`, it is at this step where we specify the inputs and any additional ouptions such as the output directory. The next step `runWorkflow.py` executes the manta workflow and will output the results based on the configuration in the first step. We should note here that manta requires python 2 in order to run, we've configured a conda environment for manta which will use python 2 so we need to activate that python environment first. From there we run the `configManta.py` script with the following options:

1. --normalBam: path to normal bam
2. --tumorBam: path to tumor bam
3. --referenceFasta: path to indexed reference fasta
4. --runDir: path to directory where results will be stored

```bash
# make directory for results
mkdir -p ~/workspace/somatic/manta_wgs
cd ~/workspace/somatic/manta_wgs

# source the manta config Environment
source activate manta

# run the manta config script
python /usr/local/bin/manta-1.4.0.centos6_x86_64/bin/configManta.py --normalBam=/workspace/align/WGS_Norm_merged_sorted_mrkdup_bqsr.bam --tumorBam=/workspace/align/WGS_Tumor_merged_sorted_mrkdup_bqsr.bam --referenceFasta=/workspace/inputs/references/genome/ref_genome.fa --runDir=/workspace/somatic/manta_wgs/
```

Next with everything configured we can go ahead and run manta. All we have to do here is run the `runWorfklow.py` script the configuration step created in the previous step and specify if we are running this locally or on a cluster with the `-m` The specific options we give and what they mean are provided below.

1. -m local: specifies we are running manta on a local computer
2. -j 8: number of jobs to launch at once
3. -g 60: run with 60 gigabytes

```bash
# run manta
python /workspace/somatic/manta_wgs/runWorkflow.py -m local -j 8 -g 60

# deactivate environment
source deactivate
```

The output of manta is quite verbose, however in gerneral the files we care about are in the `/workspace/somatic/manta_wgs/results/variants` folder. Specifically as we are interested in somatic events occuring in the tumor sample we care about the `somaticSV.vcf.gz` file.
