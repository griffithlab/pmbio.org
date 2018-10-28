---
feature_text: |
  ## Precision Medicine
title: PreAlignment QC
categories:
    - Module-02-Inputs
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-06-01
---

Assessment of data quality can be an important stage of any project involving NGS data. It is common practice to perform a pilot experiment on a relatively small number of samples and assess data quantities and quality before proceeding with a larger more expensive experiment. Issues to consider during the quality assessment / quality control stage:

* Are there any fundamental library construction issues that may affect sensitivity and specificity of the results?
* Are the DNA/RNA fragments of expected/sufficient length? What does the fragment size distribution look like?
* Does the individual read length make sense for the fragment size?
* Are there a lot of adapter sequences in the reads?
* What is the complexity of the library?  Are there a lot of duplicate reads?
* What does the sequence base error profile look like? Are there problematic positions in the sequence run that affect many reads? Are there a lot of N's? Is the overall base error rate unusually high? Does the base quality at the ends of each read become problematic?
* Do you have enough data to achieve experimental objectives? All of the above points (and more) may influence the amount of effective data obtained. How this is assessed really depends on the objectives. For example, if you goal it to identify germline variants you might determine average coverage and aim to achieve >20-30x.  The rationale for such targets is not well documented in publications that examine the effect of coverage on variant calling sensitivity.  In tumor/normal sequencing for somatic variants you might have higher targets depending on tumor purity and heterogeneity, and the target might differ for tumor and normal (e.g. 30x normal, 50x tumor).  In RNA-seq you might assess detection of genes critical to your biological system or area of research interest.

### Use fastqc to produce base quality metrics for each FastQ file

### Use multiqc to produce a combined report of these QC files

### Explore the QC results in your browser
- Go to this URL and explore


