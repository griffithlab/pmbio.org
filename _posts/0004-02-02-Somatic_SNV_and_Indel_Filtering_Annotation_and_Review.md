---
feature_text: |
  ## Precision Medicine
title: Somatic SNV and Indel Filtering Annotation and Review
categories:
    - Module 04. Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0004-02-02
---

Following automated variant calling, heuristic filtering and manual review of aligned read sequences is required to identify a high-quality list of somatic variants.


#### **Somatic Variant Filtering**
__________________________  
Somatic variant filtering requires setting heuristic cutoffs for various sequencing data to eliminate called variants that might be false positives. This requires first defining the thresholds and then employing these thresholds on variant files. Typically we set thresholds for variant allele fraction (VAF), total coverage, and variant coverage.

**Variant allele fraction (VAF):** The VAF represents the relative frequency of an allele at a particular locus, which serves as a proxy for number of cells containing the variant (i.e., 5% VAF means that 10% of the cells in the sample contain the variant). Typically, at least 5% VAF is required to make confident calls (given 20X coverage); however, this threshold is experiment-specific. For experiments with higher average coverage, the minimum VAF threshold can be reduced accordingly.

**Total coverage:** The total coverage indicates the number of sequencing reads that align to the locus of interest.

**Variant coverage:** The variant coverage indicates the number of sequencing reads aligned to the locus of interest that contain the variant.


The following code allows you to filter on the VAF, total coverage, and variant coverage:

TO DO: INSERT CODE FOR FILTERING



#### **Somatic Variant Refinement**
__________________________  



### **Manual Review**



### **Automated Refinement**




#### **Somatic Variant Annotation**
__________________________  
