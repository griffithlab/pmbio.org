---
feature_text: |
  ## Precision Medicine
title: Somatic SV Filtering/Annotation/Review
categories:
    - Module-05-Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-03-02
---

Now that we have our manta results, lets do some basic filtering annotation and visualization. First off Let's go ahead and look at how many SVs were called and if they passed or failed filters. We can accomplish this by cutting out the 7th column in the vcf file.

```bash
zcat somaticSV.vcf.gz | grep -v "^#" |cut -f 7 | sort | uniq -c
```

You should see something like this:

```bash
1 MaxMQ0Frac
1 MaxMQ0Frac;MinSomaticScore
29 MinSomaticScore
129 PASS
```

We can see that the majority of the SVs passed mantas default filters however 31 failed by either not meeting the min somatic score cutoff under mantas somatic model or the fraction of reads with a MAPQ0 (i.e. aligned to multiple locations) around the breakpoints of the SV exceeded the threshold. For a full explanation of why manta may have failed a somatic variant look at the manta documentation available [here](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#vcf-filter-fields).


Manta will detect the following SV types:

1. Insertions
2. Deletions
3. Inversions
4. Tandem Duplications
5. Interchromosomal Translocations

Lets do some command line magic to see the breakdown of SV types in our list with a bit of command line magic

```bash
zcat somaticSV.vcf.gz | grep -v "^#" | grep -iw "PASS" | cut -f 3 | tr ":" "\t" | cut -f 1 | sort | uniq -c
```

```bash
4 MantaBND
64 MantaDEL
39 MantaDUP
2 MantaINS
20 MantaINV
```

we can see from the output that we have 4 interchromosomal translocations, 64 deletions, 39 duplications, 2 insertions, and 20 inversions passing filters.

Let's look at an example of each of these in IGV to evaluate them. As an aside IGV has some nice documentation on interpreting SVs by insert size and pair orientation [here](https://software.broadinstitute.org/software/igv/interpreting_insert_size) and [here](https://software.broadinstitute.org/software/igv/interpreting_pair_orientations).

#### Interchromosomal translocation

chr6    695390

{% include figure.html image="/assets/module_5/translocation_igv.png" position="left" %}

#### Deletion

chr6    90088092

{% include figure.html image="/assets/module_5/deletion_igv.png" position="left" %}

#### Duplication

chr6    156876266

{% include figure.html image="/assets/module_5/duplication_igv.png" position="left" %}

#### Inversion

chr6    86719474

{% include figure.html image="/assets/module_5/inversion_igv.png" position="left" %}

#### Visualizing with svviz

Use a docker image to run svviz2
```bash

# create a directory for the svviz2 results
cd /workspace/somatic/manta_wgs/results/variants/
mkdir svviz2 

# use svviz2 within a docker image to produce visualiatons for our manta SV results
docker pull sridnona/svviz2:v2

# test the docker installation
docker run sridnona/svviz2:v2 /usr/local/bin/svviz2 --help

# start an interactive session with this docker image. make sure the /workspace volume is accessible inside the docker containers and mounted at the same path
docker run -i -v /workspace:/workspace -t sridnona/svviz2:v2 /bin/bash

# make sure we can run docker inside the container
svviz2 --help

# run svviz on the manta SV results, supplying the VCF, reference genome, normal WGS BAM, and tumor WGS BAM
svviz2 --ref /workspace/inputs/references/genome/ref_genome.fa --variants /workspace/somatic/manta_wgs/results/variants/somaticSV.vcf.gz /workspace/align/WGS_Norm_merged_sorted_mrkdup_bqsr.bam /workspace/align/WGS_Tumor_merged_sorted_mrkdup_bqsr.bam --outdir /workspace/somatic/manta_wgs/results/variants/svviz2

# leave the docker interactive session
exit
```

