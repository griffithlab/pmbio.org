---
feature_text: |
  ## Precision Medicine
title: Annotation
categories:
    - Module-02-Inputs
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-03-01
---

### SeqCapEZ_Exome_v3.0
The reagent used for the exome sequencing of the data used in this course was the SeqCapEZ_Exome_v3.0 from Roche Nimblegen. 
[here](https://sequencing.roche.com/en/products-solutions/by-category/target-enrichment/hybridization/seqcap-ez-exome-v3-kit.html)
```bash
cd /workspace/data/results/inputs
wget -c https://sequencing.roche.com/content/dam/rochesequence/worldwide/resources/SeqCapEZ_Exome_v3.0_Design_Annotation_files.zip
```

### liftOver chain file
```bash
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

### run Liftover
```bash
cd /workspace/data/results/inputs
liftOver SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_primary_targets.bed  hg19ToHg38.over.chain.gz SeqCap_EZ_Exome_v3_hg38_primary_targets.bed unMapped.bed
cut -f 1-3 SeqCap_EZ_Exome_v3_hg38_primary_targets.bed > SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed
```
