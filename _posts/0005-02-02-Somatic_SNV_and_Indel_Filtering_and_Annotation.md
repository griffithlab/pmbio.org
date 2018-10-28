---
feature_text: |
  ## Precision Medicine
title: Somatic SNV and Indel Filtering and Annotation
categories:
    - Module-05-Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-02-02
---

### **Basic Filtering on Somatic Variants**
First, let's do a basic filtering for `PASS` only variants on our merged and normalized vcf file:
```bash
cd /workspace/somatic
java -Xmx24g -jar /usr/local/bin/GenomeAnalysisTK.jar -T SelectVariants -R ~/workspace/inputs/references/genome/ref_genome.fa --excludeFiltered --variant ~/workspace/somatic/exome.merged.leftalignandtrim.decomposed.vcf -o ~/workspace/somatic/exome.merged.norm.pass_only.vcf
```

#### **VEP Annotation**
- `cd /data/bin/ensembl-vep` (previously installed if following tutorial)
- `unset PERL5LIB`

```bash
cd /workspace/somatic
vep -i ~/workspace/somatic/exome.merged.norm.pass_only.vcf --cache --dir /opt/vep_cache/ --format vcf --vcf --plugin Downstream --plugin Wildtype --symbol --terms SO --flag_pick --transcript_version -o ~/workspace/somatic/exome.merged.norm.annotated.vcf
```

### **Adding Bam-readcounts to VCF file**
We have added a python helper script that will take your vcf and DNA bam files and generates two bam-readcount output files, one for snv and one for indel.

```bash
cd /workspace/somatic
mkdir -p /workspace/somatic/bam_readcounts
# First activate the bam-readcount conda enironment
source activate bam-readcount
# Running bam-readcount on annotated vcf file, once using tumor bam and another time for normal bam file
python -u /usr/local/bin/bam_readcount_helper.py ~/workspace/somatic/exome.merged.norm.annotated.vcf 'TUMOR' ~/workspace/inputs/references/genome/ref_genome.fa /workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.bam ~/workspace/somatic/bam_readcounts/
python -u /usr/local/bin/bam_readcount_helper.py ~/workspace/somatic/exome.merged.norm.annotated.vcf 'NORMAL' ~/workspace/inputs/references/genome/ref_genome.fa /workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam ~/workspace/somatic/bam_readcounts/
```
Once we concatenate the snv readcount file as well as the indel readcount file, we can then run the vcf-annotation-tool to add these bam-readcounts to our vcf output.
```bash
cd ~/workspace/somatic/bam_readcounts/
cat TUMOR_bam_readcount_snv.tsv TUMOR_bam_readcount_indel.tsv > TUMOR_bam_readcount_combined.tsv
cat NORMAL_bam_readcount_snv.tsv NORMAL_bam_readcount_indel.tsv > NORMAL_bam_readcount_combined.tsv

mkdir -p ~/workspace/somatic/final/
cd ~/workspace/somatic/final/
vcf-readcount-annotator ~/workspace/somatic/exome.merged.norm.annotated.vcf ~/workspace/somatic/bam_readcounts/TUMOR_bam_readcount_combined.tsv DNA -s TUMOR -o ~/workspace/somatic/final/exome.final.vcf

```
### **Generating Table from VCF file**


### **Additional Filters **

**Please continue to the next section for instructions on how to perform manual review on these somatic variant results*
