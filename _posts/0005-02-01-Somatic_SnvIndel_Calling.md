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

python2 /usr/local/bin/strelka-2.7.1.centos5_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam=/workspace/data/DNA_alignments/chr6+chr17/final/Exome_Norm_sorted_mrkdup_bqsr.bam --tumorBam=/workspace/data/DNA_alignments/chr6+chr17/final/Exome_Tumor_sorted_mrkdup_bqsr.bam --referenceFasta=/workspace/data/raw_data/references/ref_genome.fa --exome --runDir=/workspace/data/results/somatic/strelka

#Please specify according to the number of cpus available or how many you would like to allocate to this job. In this case, four were given.
python2 /workspace/data/results/somatic/strelka/runWorkflow.py -m local -j 4
cd ~/workspace/data/results/somatic/strelka/results/variants
zcat somatic.snvs.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > somatic.snvs.gt.vcf
zcat somatic.indels.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > somatic.indels.gt.vcf
find ~/workspace/data/results/somatic/strelka/results/variants/ -name "*.vcf" -exec bgzip -f {} \;
find ~/workspace/data/results/somatic/strelka/results/variants/ -name "*.vcf.gz" -exec tabix -f {} \;

bcftools concat -a -o exome.vcf.gz -O z somatic.snvs.gt.vcf.gz somatic.indels.gt.vcf.gz

tabix exome.vcf.gz
```

#### **Running MuTect2**
__________________________
The final variant caller that we will also use results from is [MuTect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php). MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original MuTect (Cibulskis et al., 2013) with the assembly-based machinery of HaplotypeCaller.

```bash
#Obtaining germline resource from GATK
cd ~/workspace/data/raw_data/references
gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz ./

mkdir -p ~/workspace/data/results/somatic/mutect
cd ~/workspace/data/results/somatic/mutect

#Creating a panel of normals
/usr/local/bin/gatk Mutect2 -R ~/workspace/data/raw_data/references/ref_genome.fa -I ~/workspace/data/DNA_alignments/chr6+chr17/final/Exome_Norm_sorted_mrkdup_bqsr.bam -tumor-sample HCC1395BL_DNA -O Exome_Norm_PON.vcf.gz

#Running Mutect2
/usr/local/bin/gatk Mutect2 -R ~/workspace/data/raw_data/references/ref_genome.fa -I ~/workspace/data/DNA_alignments/chr6+chr17/final/Exome_Tumor_sorted_mrkdup_bqsr.bam -tumor HCC1395_DNA -I ~/workspace/data/DNA_alignments/chr6+chr17/final/Exome_Norm_sorted_mrkdup_bqsr.bam -normal HCC1395BL_DNA --germline-resource ~/workspace/data/raw_data/references/af-only-gnomad.hg38.vcf.gz --af-of-alleles-not-in-resource 0.00003125 --panel-of-normals ~/workspace/data/results/somatic/mutect/Exome_Norm_PON.vcf.gz -O ~/workspace/data/results/somatic/mutect/exome.vcf.gz -L chr6 -L chr17

echo ~/workspace/data/results/somatic/mutect/exome.vcf.gz > ~/workspace/data/results/somatic/mutect/exome_vcf.fof
bcftools concat --allow-overlaps --remove-duplicates --file-list ~/workspace/data/results/somatic/mutect/exome_vcf.fof --output-type z --output ~/workspace/data/results/somatic/mutect/mutect_exome.vcf.gz
mv mutect_exome.vcf.gz exome.vcf.gz
tabix ~/workspace/data/results/somatic/mutect/exome.vcf.gz
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
