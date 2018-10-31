---
feature_text: |
  ## Precision Medicine
title: Somatic SNV/InDel Calling
categories:
    - Module-05-Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0005-02-01
---

ADD INTRO

#### Running VARSCAN
__________________________  

The first variant caller that we will use here is [VARSCAN](http://varscan.sourceforge.net/), VarScan is a platform-independent mutation caller for targeted, exome, and whole-genome resequencing data and employs a robust heuristic/statistic approach to call variants that meet desired thresholds for read depth, base quality, variant allele frequency, and statistical significance:
```bash
mkdir -p ~/workspace/somatic/varscan
cd ~/workspace/somatic/varscan
# Runtime: ~5min
java -Xmx24g -jar /usr/local/bin/VarScan.v2.4.2.jar somatic <(samtools mpileup -l /workspace/inputs/references/exome/exome_regions.bed --no-BAQ -f /workspace/inputs/references/genome/ref_genome.fa /workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam /workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.bam) /workspace/somatic/varscan/exome --mpileup 1 --output-vcf

java -Xmx24g -jar /usr/local/bin/VarScan.v2.4.2.jar processSomatic exome.snp.vcf exome.snp
java -Xmx24g -jar /usr/local/bin/VarScan.v2.4.2.jar processSomatic exome.indel.vcf exome.indel
find ~/workspace/somatic/varscan -name '*.vcf' -exec bgzip -f {} \;
find ~/workspace/somatic/varscan -name '*.vcf.gz' -exec tabix -f {} \;

gatk VariantFiltration -R /workspace/inputs/references/genome/ref_genome.fa -V exome.snp.Somatic.vcf.gz --mask exome.snp.Somatic.hc.vcf.gz --mask-name "processSomatic" --filter-not-in-mask -O exome.snp.Somatic.hc.filter.vcf.gz
gatk VariantFiltration -R /workspace/inputs/references/genome/ref_genome.fa -V exome.indel.Somatic.vcf.gz --mask exome.indel.Somatic.hc.vcf.gz --mask-name "processSomatic" --filter-not-in-mask -O exome.indel.Somatic.hc.filter.vcf.gz

bcftools concat -a -o exome.vcf.gz -O z exome.snp.Somatic.hc.filter.vcf.gz exome.indel.Somatic.hc.filter.vcf.gz
tabix -f /workspace/somatic/varscan/exome.vcf.gz
```

#### **Running STRELKA**
__________________________  

The second variant caller that we will use is [STRELKA](https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md). Strelka calls germline and somatic small variants from mapped sequencing reads and is optimized for rapid clinical analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs. Both germline and somatic callers include a final empirical variant rescoring step using a random forest model to reflect numerous features indicative of call reliability which may not be represented in the core variant calling probability model.

```bash
mkdir -p ~/workspace/somatic/strelka
cd ~

source activate strelka
/usr/local/bin/strelka-2.7.1.centos5_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam=/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam --tumorBam=/workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.bam --referenceFasta=/workspace/inputs/references/genome/ref_genome.fa --exome --runDir=/workspace/somatic/strelka
source deactivate
#Please specify according to the number of cpus available or how many you would like to allocate to this job. In this case, four were given.
# Runtime: ~ 3min
python2 /workspace/somatic/strelka/runWorkflow.py -m local -j 8

cd ~/workspace/somatic/strelka/results/variants
zcat somatic.snvs.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > somatic.snvs.gt.vcf
zcat somatic.indels.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > somatic.indels.gt.vcf
find ~/workspace/somatic/strelka/results/variants/ -name "*.vcf" -exec bgzip -f {} \;
find ~/workspace/somatic/strelka/results/variants/ -name "*.vcf.gz" -exec tabix -f {} \;

bcftools concat -a -o exome.vcf.gz -O z somatic.snvs.gt.vcf.gz somatic.indels.gt.vcf.gz

tabix exome.vcf.gz
```

#### **Running MuTect2**
__________________________

The final variant caller that we will also use results from is [MuTect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php). MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original MuTect (Cibulskis et al., 2013) with the assembly-based machinery of HaplotypeCaller.

```bash
#Obtaining germline resource from GATK
cd ~/workspace/inputs/references
gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz .
gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi .

mkdir -p ~/workspace/somatic/mutect
cd ~/workspace/somatic/mutect

#Creating a panel of normals
# Runtime: ~ 17min
gatk --java-options "-Xmx24G" Mutect2 -R ~/workspace/inputs/references/genome/ref_genome.fa -I ~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam -tumor-sample HCC1395BL_DNA -O Exome_Norm_PON.vcf.gz

#Running Mutect2 Using latest version of GATK
# Runtime: ~20m
gatk -nt 8 --java-options "-Xmx24G" Mutect2 -R ~/workspace/inputs/references/genome/ref_genome.fa -I ~/workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.bam -tumor HCC1395_DNA -I ~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam -normal HCC1395BL_DNA --germline-resource ~/workspace/inputs/references/af-only-gnomad.hg38.vcf.gz --af-of-alleles-not-in-resource 0.00003125 --panel-of-normals ~/workspace/somatic/mutect/Exome_Norm_PON.vcf.gz -O ~/workspace/somatic/mutect/exome.vcf.gz -L chr6 -L chr17
#Need to change header sample names in order to combine variants with those from other algorithms
sed -i 's/HCC1395BL_DNA/NORMAL/' exome.vcf.gz
sed -i 's/HCC1395_DNA/TUMOR/' exome.vcf.gz
#Running mutect2 using gatk version 3.6
#java -Xmx12g -jar /usr/local/bin/GenomeAnalysisTK.jar -T MuTect2 --disable_auto_index_creation_and_locking_when_reading_rods -R ~/workspace/data/raw_data/references/ref_genome.fa -I:tumor ~/workspace/data/DNA_alignments/chr6+chr17/final/Exome_Tumor_sorted_mrkdup_bqsr.bam -I:Normal ~/workspace/data/DNA_alignments/chr6+chr17/final/Exome_Norm_sorted_mrkdup_bqsr.bam --dbsnp ~/workspace/data/raw_data/references/Homo_sapiens_assembly38.dbsnp138.vcf.gz --cosmic ~/workspace/data/raw_data/references/Cosmic_v79.dictsorted.vcf.gz -o ~/workspace/data/results/somatic/mutect/exome.vcf.gz -L ~/workspace/data/results/inputs/SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.interval_list

echo ~/workspace/somatic/mutect/exome.vcf.gz > ~/workspace/somatic/mutect/exome_vcf.fof
bcftools concat --allow-overlaps --remove-duplicates --file-list ~/workspace/somatic/mutect/exome_vcf.fof --output-type z --output ~/workspace/somatic/mutect/mutect_exome.vcf.gz
mv mutect_exome.vcf.gz exome.vcf.gz
tabix ~/workspace/somatic/mutect/exome.vcf.gz

```

#### **Merge Variants**
__________________________

With outputs from all three algorithms, we can now merge the variants to generate a comprehensive list of detected variants:
Installation of GATK version 3.6 is need for the further processing of our variants.

```bash
sudo bash
# Installing gatk 3.6
cd /usr/local/bin/
wget http://genomedata.org/pmbio-workshop/references/gatk/GenomeAnalysisTK-3.6-0-g89b7209.tar.bz2
tar --bzip2 -xvf GenomeAnalysisTK-3.6-0-g89b7209.tar.bz2

# Testing the installation
java -jar GenomeAnalysisTK.jar -h

exit
```

```bash
# Unzip the vcf.gz files before combining Variants
cd ~/workspace/somatic
gunzip ~/workspace/somatic/varscan/exome.vcf.gz
gunzip ~/workspace/somatic/strelka/results/variants/exome.vcf.gz
gunzip ~/workspace/somatic/mutect/exome.vcf.gz

# (UNIQUIFY command) java -Xmx4g -jar /usr/local/bin/GenomeAnalysisTK.jar -T CombineVariants -R ~/workspace/data/raw_data/references/ref_genome.fa -genotypeMergeOptions UNIQUIFY --variant:varscan ~/workspace/data/results/somatic/varscan/exome.vcf --variant:strelka ~/workspace/data/results/somatic/strelka/results/variants/exome.vcf --variant:mutect ~/workspace/data/results/somatic/mutect/new_gatk_files/exome.vcf -o ~/workspace/data/results/somatic/exome.unique.vcf.gz
java -Xmx24g -jar /usr/local/bin/GenomeAnalysisTK.jar -T CombineVariants -R ~/workspace/inputs/references/genome/ref_genome.fa -genotypeMergeOptions PRIORITIZE --rod_priority_list mutect,varscan,strelka --variant:varscan ~/workspace/somatic/varscan/exome.vcf --variant:strelka ~/workspace/somatic/strelka/results/variants/exome.vcf --variant:mutect ~/workspace/somatic/mutect/exome.vcf -o ~/workspace/somatic/exome.merged.vcf.gz
```

### **Left Align and Trim**
__________________________
Reference for explaining left align and trim:
https://genome.sph.umich.edu/wiki/Variant_Normalization#Left_alignment
```bash
cd ~/workspace/somatic/

gatk --java-options "-Xmx24G" LeftAlignAndTrimVariants -V ~/workspace/somatic/exome.merged.vcf.gz -O exome.merged.leftalignandtrim.vcf -R ~/workspace/inputs/references/genome/ref_genome.fa
```

Note that when running on chromosome 6 and 17 merged variants file, this gave 0 variants aligned.

### **Splitting Multi-allelic Variant**
__________________________
```bash
cd ~/workspace/somatic/

vt decompose -s /workspace/somatic/exome.merged.leftalignandtrim.vcf -o /workspace/somatic/exome.merged.leftalignandtrim.decomposed.vcf
```

**Please continue to the next section for instructions on how to filter, annotation and review variants**
