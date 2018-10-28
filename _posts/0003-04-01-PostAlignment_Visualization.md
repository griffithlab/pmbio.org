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

Some simple exercises for exploring the Exome BAMs in IGV

* After loading the exome BAMs and making any adjustments to track labeling, save your IGV session for convenience in case you want to load these BAMs again later.
* Find an example gene (e.g. BRCA1 on chr17)
* Try viewing the reads as read pairs. How big do your fragments look?
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

