---
feature_text: |
  ## Precision Medicine
title: Germline SNV and Indel Filtering, Annotation, and Review
categories:
    - Module 05. Germline
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-02-02
---

Filter variants using GATK VQSR approach

For inspiration see: https://software.broadinstitute.org/gatk/documentation/article?id=2805

Note, for running VQSR on exome data. There will be a smaller number of variants per sample compared to WGS. These are typically insufficient to build a robust recalibration model. If running on only a few samples, GATK recommends that you analyze samples jointly in cohorts of at least 30 samples. If necessary, add exomes from 1000G Project or comparable. These should be processed with similar technical generation (technology, capture, read length, depth). 

See: https://drive.google.com/drive/folders/0BzI1CyccGsZiWlU5SXNvNnFGbVE (GATK Workshops -> 1709 -> GATKwr19-05-Variant_filtering.pdf)

Another possibility is hard filtering: https://software.broadinstitute.org/gatk/documentation/article?id=2806



