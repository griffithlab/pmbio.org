---
feature_text: |
  ## Precision Medicine
title: Somatic CNV Calling
categories:
    - Module 04. Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0004-04-01
---

[Copy number alterations (CNA)](https://en.wikipedia.org/wiki/Copy-number_variation) occur when sections of a genome are duplicated or deleted. This phenomenom can actually be quite usefull from an evolutionary standpoint, an example would be the duplication of opsin genes allowing some vertebrate species to see more colors. These types of events however can have a significant impact in the context of disease with perhaps the most famous being an amplification of chromosome 21 resulting in down sydrome. In this section we will go over identifying these types of alterations with [VarScan](http://varscan.sourceforge.net/using-varscan.html) and [Delly](https://github.com/dellytools/delly). However first let's examine exactly what we mean when the segment of a genome is duplicated or deleted and how these types of events can be identified in sequencing data.

{% include figure.html image="/assets/module_3/CNA_illustration.png" position="left" width="450" %}

To begin, the concept of a CNA is fairly straight forward, in the figure to the left we show a standard pair of chromosomes divided into 8 segments. In a copy number amplification a region of the genome is duplicated, in our figure region G has 2 extra segments right after each other. When looking at this region in the reference genome we would expect to see a sharp and dramatic increase in coverage at this region. Further with paired end data the insert size for reads spanning the amplification would be larger. In contrast a deletion is how it sounds and is simply a segment of the chromosome which is gone. In our figure we show a single copy deletion of segment B which would result in a sharp drop in coverage at that region in the sequencing data and a shorter average insert size at the breakpoints of the event.

### VarScan

talk about how algorithm works

### Delly

talk about how algorithm works
