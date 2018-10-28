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
We have added a python helper script that will take your vcf and DNA bam files and generates twp bam-readcount output files, one for snv and one for indel.

`python -u /usr/bin/bam_readcount_helper.py <vcf_file> 'TUMOR' <reference fasta> <bam file> <output_dir>`

Once we concatenate the snv readcount file as well as the indel readcount file, we can then run the vcf-annotation-tool to add these bam-readcounts to our vcf output.

`head -1 <snv_readcount_file> > <combined_readcount_file>; tail -n +2 -q <indel_readcount_file> >> <combined_readcount_file>`

`vcf-readcount-annotator <vcf_file> <readcount_file_combined> DNA -s TUMOR -o <output_dir>`

### **Generating Table from VCF file**


### **Additional Filters **

**Please continue to the next section for instructions on how to perform manual review on these somatic variant results*
