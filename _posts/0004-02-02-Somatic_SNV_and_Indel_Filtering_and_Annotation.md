---
feature_text: |
  ## Precision Medicine
title: Somatic SNV and Indel Filtering and Annotation
categories:
    - Module-04-Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0004-02-02
---

#### **Basic Filtering on Somatic Variants**
Here, we are first doing a basic filtering for `PASS` only variants on our merged vcf file:

`java -Xmx4g -jar /data/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T SelectVariants -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa --excludeFiltered --variant /data/exome_chr6.merged.vcf.gz -o /data/exome_chr6.merged.pass.vcf.gz`


##### **VEP Annotation**
- `cd /data/bin/ensembl-vep` (previously installed if following tutorial)
- `unset PERL5LIB`
- `wget -O /data/vep_cache/Plugins/Wildtype.pm https://raw.githubusercontent.com/griffithlab/pVAC-Seq/master/pvacseq/VEP_plugins/Wildtype.pm --no-check-certificate`
- Need to make sure cache files installed: all build38 related (41, 43, 45)
- Fasta -> option 28 Homo Saipiens
- Plugins -> downstream option 14

`/usr/bin/perl /data/bin/ensembl-vep/vep.pl -i /data/exome_chr6.merged.pass.vcf.gz --cache --dir /data/vep_cache/ --format vcf --vcf --plugin Downstream --plugin Wildtype --symbol --terms SO --flag_pick -o /data/exome_chr6.merged.pass.annotated.vcf.gz`

#### **Adding Bam-readcounts to VCF file**
We have added a python helper script that will take your vcf and DNA bam files and generates twp bam-readcount output files, one for snv and one for indel.

`python -u /usr/bin/bam_readcount_helper.py <vcf_file> 'TUMOR' <reference fasta> <bam file> <output_dir>`

Once we concatenate the snv readcount file as well as the indel readcount file, we can then run the vcf-annotation-tool to add these bam-readcounts to our vcf output.

`head -1 <snv_readcount_file> > <combined_readcount_file>; tail -n +2 -q <indel_readcount_file> >> <combined_readcount_file>`

`vcf-readcount-annotator <vcf_file> <readcount_file_combined> DNA -s TUMOR -o <output_dir>`

#### **Generating Table from VCF file**


#### **Additional Filters **

**Please continue to the next section for instructions on how to perform manual review on these somatic variant results*
