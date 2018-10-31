---
feature_text: |
  ## Precision Medicine
title: Somatic SNV and Indel Filtering and Annotation
categories:
    - Module-05-Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-02-02
---

### Basic Filtering on Somatic Variants

First, let's do a basic filtering for `PASS` only variants on our merged and normalized vcf file:
```bash
cd /workspace/somatic
java -Xmx24g -jar /usr/local/bin/GenomeAnalysisTK.jar -T SelectVariants -R ~/workspace/inputs/references/genome/ref_genome.fa --excludeFiltered --variant ~/workspace/somatic/exome.merged.leftalignandtrim.decomposed.vcf -o ~/workspace/somatic/exome.merged.norm.pass_only.vcf
```

#### VEP Annotation
[VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html) stands for Variant Effect Predictor. We will use it to annotate our variants to determine the effect of the variants (e.g. SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions. (McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F.
The Ensembl Variant Effect Predictor.
Genome Biology Jun 6;17(1):122. (2016)
doi:10.1186/s13059-016-0974-4)

```bash
cd /workspace/somatic
vep -i ~/workspace/somatic/exome.merged.norm.pass_only.vcf --cache --dir /opt/vep_cache/ --format vcf --vcf --plugin Downstream --plugin Wildtype --symbol --terms SO --flag_pick --transcript_version -o ~/workspace/somatic/exome.merged.norm.annotated.vcf
```

### Adding Bam-readcounts to VCF file

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
vcf-readcount-annotator ~/workspace/somatic/exome.merged.norm.annotated.vcf ~/workspace/somatic/bam_readcounts/TUMOR_bam_readcount_combined.tsv DNA -s TUMOR -o ~/workspace/somatic/final/exome.tumordna_annotated.vcf.gz
vcf-readcount-annotator ~/workspace/somatic/final/exome.tumordna_annotated.vcf.gz ~/workspace/somatic/bam_readcounts/NORMAL_bam_readcount_combined.tsv DNA -s NORMAL -o ~/workspace/somatic/final/exome.annotated.vcf.gz
# Remove the intermediate file
rm exome.tumordna_annotated.vcf.gz
tabix -p vcf exome.annotated.vcf.gz
```
### Generating Table from VCF file

```bash
# Adjust the output fields accordingly
cd ~/workspace/somatic/final/
gunzip exome.annotated.vcf.gz
java -Xmx4g -jar /usr/local/bin/GenomeAnalysisTK.jar -T VariantsToTable -R ~/workspace/inputs/references/genome/ref_genome.fa --variant ~/workspace/somatic/final/exome.annotated.vcf -F CHROM -F POS -F ID -F REF -F ALT -F set -F AC -F AF -o variants.tsv
```
### Additional Filters
Filter vcf allele frequency

**Please continue to the next section for instructions on how to perform manual review on these somatic variant results*
