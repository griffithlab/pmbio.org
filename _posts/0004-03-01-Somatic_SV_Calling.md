---
feature_text: |
  ## Precision Medicine
title: Somatic SV Calling
categories:
    - Module-04-Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0004-03-01
---

### Manta
```bash
python2 /usr/local/bin/manta-1.4.0.centos6_x86_64/bin/configManta.py --normalBam=/workspace/data/results/align/WGS_Norm_merged_sorted_mrkdup.bam --tumorBam=/workspace/data/results/align/WGS_Tumor_merged_sorted_mrkdup.bam --referenceFasta=/workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa --runDir=/workspace/data/results/somatic/manta/

python2 /workspace/data/results/somatic/manta/runWorkflow.py
```
