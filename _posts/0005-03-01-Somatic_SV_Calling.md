---
feature_text: |
  ## Precision Medicine
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

Running manta is a two step process. The first step is configuring manta to run with `configManta.py`, it is at this step where we specify the inputs and any additional ouptions such as the output directory. The next step `runWorkflow.py` executes the manta workflow and will output the results based on the configuration in the first step. We should note here that manta requires python 2 in order to run, we've configured a conda environment for manta which will use python2 so we need to activate that python environment first. 

```bash
# source the manta config Environment
source activate manta

# run the manta config script
python2 /usr/local/bin/manta-1.4.0.centos6_x86_64/bin/configManta.py --normalBam=/workspace/align/WGS_Norm_merged_sorted_mrkdup.bam --tumorBam=/workspace/align/WGS_Tumor_merged_sorted_mrkdup.bam --referenceFasta=/workspace/inputs/references/genome/ref_genome.fa --runDir=/workspace/somatic/manta/

# run manta
python2 /workspace/somatic/manta/runWorkflow.py

# deactivate environment
source deactivate
```
