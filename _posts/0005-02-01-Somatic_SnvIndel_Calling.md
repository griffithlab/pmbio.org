---
feature_text: |
  ## Precision Medicine
title: Somatic SNV/InDel Calling
categories:
    - Module-05-Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-02-01
---

#### Downloading Reference Files
__________________________  
In order to run the following variant detection algorithms, you will need to download a few reference files from `genome_data.org` to your directory for reference file storage:

```bash
mkdir -p ~/workspace/references
cd ~/workspace/references
wget genomedata.org/pmbio-workshop/references/NimbleGenExome_v3.interval_list
```

#### Running VARSCAN
__________________________  

The first variant caller that we will use here is [VARSCAN](http://varscan.sourceforge.net/), VarScan is a platform-independent mutation caller for targeted, exome, and whole-genome resequencing data and employs a robust heuristic/statistic approach to call variants that meet desired thresholds for read depth, base quality, variant allele frequency, and statistical significance:
```bash
mkdir -p ~/workspace/results/somatic/varscan

java -Xmx24g -jar /usr/local/bin/VarScan.v2.4.2.jar somatic <(samtools mpileup -l ~/workspace/results/inputs/SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed --no-BAQ -f ~/workspace/references/genome/GRCh38_full_analysis_set_plus_decoy_hla.fa ~/workspace/results/align/final/Exome_Norm_sorted_mrkdup_bqsr.bam ~/workspace/results/align/final/Exome_Tumor_sorted_mrkdup_bqsr.bam) ~/workspace/results/somatic/varscan/exome --mpileup 1 --output-vcf

cd ~/workspace/results/somatic/varscan/

java -Xmx24g -jar /usr/local/bin/VarScan.v2.4.2.jar processSomatic exome.snp.vcf exome.snp
java -Xmx24g -jar /usr/local/bin/VarScan.v2.4.2.jar processSomatic exome.indel.vcf exome.indel
find ~/workspace/data/results/somatic/varscan -name '*.vcf' -exec bgzip -f {} \;
find ~/workspace/data/results/somatic/varscan -name '*.vcf.gz' -exec tabix -f {} \;

/usr/local/bin/gatk VariantFiltration -R ~/workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa -V exome.snp.Somatic.vcf.gz --mask exome.snp.Somatic.hc.vcf.gz --mask-name "processSomatic" --filter-not-in-mask -O exome.snp.Somatic.hc.filter.vcf.gz
/usr/local/bin/gatk VariantFiltration -R ~/workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa -V exome.indel.Somatic.vcf.gz --mask exome.indel.Somatic.hc.vcf.gz --mask-name "processSomatic" --filter-not-in-mask -O exome.indel.Somatic.hc.filter.vcf.gz

bcftools concat -a -o exome.vcf.gz -O z exome.snp.Somatic.hc.filter.vcf.gz exome.indel.Somatic.hc.filter.vcf.gz
tabix -f ~/workspace/data/results/somatic/varscan/exome.vcf.gz
```

#### **Running STRELKA**
__________________________  
The second variant caller that we will use is [STRELKA](https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md). Strelka calls germline and somatic small variants from mapped sequencing reads and is optimized for rapid clinical analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs. Both germline and somatic callers include a final empirical variant rescoring step using a random forest model to reflect numerous features indicative of call reliability which may not be represented in the core variant calling probability model.

```bash
mkdir -p ~/workspace/data/results/somatic/strelka
cd ~

python2 /usr/local/bin/strelka-2.7.1.centos5_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam=/workspace/data/DNA_alignments/chr6+chr17/final/Exome_Norm_sorted_mrkdup_bqsr.bam --tumorBam=/workspace/data/DNA_alignments/chr6+chr17/final/Exome_Tumor_sorted_mrkdup_bqsr.bam --referenceFasta=/workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa --exome --runDir=/workspace/data/results/somatic/strelka

/data/strelka/exome/runWorkflow.py -m local -j 4 (Please specify according to the number of cpus avaliable. In this case I have 4 for my instance)
zcat exome/results/variants/somatic.snvs.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > exome/results/variants/somatic.snvs.gt.vcf
zcat exome/results/variants/somatic.indels.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > exome/results/variants/somatic.indels.gt.vcf
find exome/results/variants -name *.vcf -exec bgzip -f {} \;
find exome/results/variants -name *.vcf.gz -exec tabix -f {} \;

bcftools concat -a -o exome.vcf.gz -O z exome/results/variants/somatic.snvs.gt.vcf.gz exome/results/variants/somatic.indels.gt.vcf.gz
tabix exome.vcf.gz
```

#### **Running MuTect2**
__________________________
Before preceeding, you will need to download COSMIC reference file (e.g. in folder `/data/refseq/`) after registering on their website:
```bash
* sftp <your registered email address>@sftp-cancer.sanger.ac.uk
* get /files/grch38/cosmic/v79/VCF/CosmicCodingMuts.vcf.gz
* get /files/grch38/cosmic/v79/VCF/CosmicNonCodingVariants.vcf.gz
* Quit
* zgrep "^#" CosmicCodingMuts.vcf.gz > VCF_Header
* zgrep -v "^#" CosmicCodingMuts.vcf.gz | awk '{print "chr"$0}' | sed 's/^chrMT/chrM/' > CosmicCodingMuts.clean
* zgrep -v "^#" CosmicNonCodingVariants.vcf.gz | awk '{print "chr"$0}' | sed 's/^chrMT/chrM/' > CosmicNonCodingVariants.clean
* cat CosmicCodingMuts.clean CosmicNonCodingVariants.clean | sort -gk 2,2 > Cosmic_v79
* cat VCF_Header Cosmic_v79 > Cosmic_v79.vcf
* rm VCF_Header CosmicCodingMuts.clean CosmicNonCodingVariants.clean Cosmic_v79
```

With picard tools installed:
```bash
* `java -jar picard.jar SortVcf I=/data/refseq/Cosmic_v79.vcf O=/data/refseq/Cosmic_v79.dictsorted.vcf SEQUENCE_DICTIONARY=/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa`
* `java -Xmx4g -jar GenomeAnalysisTK.jar -T MuTect2 --disable_auto_index_creation_and_locking_when_reading_rods -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I:tumor /data/alignment/final/Exome_Tumor_sorted_mrkdup_bqsr.bam -I:Normal /data/alignment/final/Exome_Norm_sorted_mrkdup_bqsr.bam --dbsnp /data/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz --cosmic /data/refseq/Cosmic_v79.dictsorted.vcf.gz -o /data/mutect/exome.vcf.gz -L /data/refseq/NimbleGenExome_v3.interval_list`
* `echo /data/mutect/exome.vcf.gz > /data/mutect/exome_vcf.fof`
* `bcftools concat --allow-overlaps --remove-duplicates --file-list /data/mutect/exome_vcf.fof --output-type z --output /data/mutect/exome.vcf.gz`
* `tabix /data/mutect/exome.vcf.gz`
```
#### **Merge Variants**
__________________________
With outputs from all three algorithms, we can now merge the variants to generate a comprehensive list of detected variants:
```bash
* java -Xmx4g -jar GenomeAnalysisTK.jar -T CombineVariants -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -genotypeMergeOptions UNIQUIFY --variant:varscan /data/varscan/exome.vcf.gz --variant:strelka /data/strelka/exome.vcf.gz --variant:mutect /data/mutect/exome.vcf.gz -o /data/exome.unique.vcf.gz
* java -Xmx4g -jar GenomeAnalysisTK.jar -T CombineVariants -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -genotypeMergeOptions PRIORITIZE --rod_priority_list mutect,varscan,strelka --variant:varscan /data/varscan/exome.vcf.gz --variant:strelka strelka/exome.vcf.gz --variant:mutect /data/mutect/exome.vcf.gz -o /data/exome.merged.vcf
```

### **Left Align and Trim**
__________________________
Reference for explaining left align and trim:
https://genome.sph.umich.edu/wiki/Variant_Normalization#Left_alignment
```bash
java -Xmx4g -jar /data/bin/gatk-4.0.10.1/gatk-package-4.0.10.1-local.jar LeftAlignAndTrimVariants -V /data/exome_chr6.merged.vcf -O exome_chr6.merged.leftalignandtrim.vcf -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

Note that when running on chromosome 6 merged variants file, this gave 0 variants aligned.

### **Splitting Multi-allelic Variant**
__________________________
```bash
/data/bin/vt/vt decompose -s exome_chr6.merged.leftalignandtrim.vcf -o exome_chr6.merged.leftalignandtrim.decomposed.vcf
```

**Please continue to the next section for instructions on how to filter, annotation and review variants**
