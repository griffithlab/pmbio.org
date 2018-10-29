---
feature_text: |
  ## Precision Medicine
title: Post Alignment Visualization
categories:
    - Module-03-Align
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-04-01
---

### View some of our alignments in IGV
Let's take a look at some of the aligments we just produced using IGV. Since we have configured our AWS instances to serve all of the files produced by the exercises we can load our BAM files by URL.  Since the BAM files are indexed, only the information we request to view will be transferred from the AWS instance to our browser for viewing.

To load the exome BAMs your URLs should look like these (don't forget to substitute your number for `#`):

* http://s#.pmbio.org/align/Exome_Norm_merged_sorted_mrkdup_bqsr.bam
* http://s#.pmbio.org/align/Exome_Tumor_merged_sorted_mrkdup_bqsr.bam

Let's go through some simple exercises for exploring the Exome BAMs in IGV.

#### Load the exome BAMs

* If IGV is already loaded, start a new session with `File` -> `New Session`
* Load the bam files (using the URLs above) with `File` -> `Load from URL`
* Rename the tracks (e.g., Exome Norm Coverage, Exome Norm BAM, Exome Tumor Coverage, Exome Tumor BAM) by right-clicking on each track and selecting `Rename Track`
* Save session with `File` -> `Save Session`

#### View alignments for an example gene

For example, let's take a look at TP53 on chr17. Simply type `TP53` into the search box and hit `Go` or enter.

{% include figure.html image="/assets/module_2/igv_exome_tumor_norm_tp53.png" width="1000" %}
 
Notice the coverage peaks roughly centered around each protein-coding exon. Notice that there are several variants. 

### EXERCISE

Explore the TP53 sequence and answer the following questions. How many potential variants can you find? Are they SNVs or indels? Are they coding or non-coding? Are the homozygous or heterozygous? Are they germline or somatic? How can you tell?

{% include question.html question="Hint" answer='Expand the Gene track to see where TP53 starts. In the collapsed view it appears to be contiguous with WRAP53 due to overlapping transcripts.' %}
{% include question.html question="Hint" answer='Zoom into the ~7-8kb region with coding exons (thicker than UTR). Look for colored bars in the coverage track for SNVs and stacks of black gaps (deletions) or purple bars (insertions).' %}
{% include question.html question="Answer" answer='There appear to be at least 5 SNVs, 1 insertion and 3 deletions. Two of the SNVs are in coding exons while the rest of the variants are intronic. Note: The indels all look like potential artifacts due to their proximity to homopolymer stretches or repetitive sequences. All of the SNVs appear homozygous (or hemizygous). One SNV appears to be somatic. In all cases the VAFs of the SNVs are at or near 100% and one of them is only observed in the tumor.' %}

This is one the variants you should have found in TP53. Which one?

{% include figure.html image="/assets/module_2/igv_exome_tumor_norm_tp53_somatic_example.png" width="1000" %}

Now, try viewing the reads as read pairs. How big do your fragments look?

{% include question.html question="Answer" answer='The fragment sizes vary. But, estimating by eye, and clicking on a few read pairs to get the insert size, it appears that fragments are in the 200-500bp range' %}

Here is a representative region showing reads in paired mode.

{% include figure.html image="/assets/module_2/igv_exome_tumor_norm_tp53_pairs_example.png" width="1000" %}

### MORE EXERCISES

* Looking at the exons of an example gene, what does the average coverage level look like?
* Try loading the .bed file for the exome reagent. How does the coverage pattern compare to the coordinates of targeted regions?

To load the WGS BAMs, your URLs should look like these:

* http://s#.pmbio.org/align/WGS_Norm_merged_sorted_mrkdup_bqsr.bam
* http://s#.pmbio.org/align/WGS_Tumor_merged_sorted_mrkdup_bqsr.bam

Some simple exercises for exploring the WGS BAMs in IGV

* After loading the exome BAMs and making any adjustments to track labeling, save your IGV session for convenience in case you want to load these BAMs again later.
* What does the average WGS coverage look like for Normal and Tumor?
* Color the reads according to read groups
* ...

