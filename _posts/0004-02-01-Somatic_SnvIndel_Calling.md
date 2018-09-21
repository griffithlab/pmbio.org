---
feature_text: |
  ## Precision Medicine
title: Somatic SNV/InDel Calling
categories:
    - Module 04. Somatic
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0004-02-01
---
TODO: Change directory names as needed for consistency across modules, also we will need a version for subset samples.
#### **Downloading Reference Files**
__________________________  
In order to run the following variant detection algorithms, you will need to download a few reference files from `genome_data.org` to your directory for reference file storage (e.g. `/data/refseq/`):
* `hglft_genome_304d_b78af0.bed`
* `NimbleGenExome_v3.interval_list`

#### **Running VARSCAN**
__________________________  

Given that you have VARSCAN properly installed, here are the commands for running VARSCAN in order, you may need to adjust according to how you have named and placed your directories:

* `java -Xmx4g -jar VarScan.v2.4.2.jar somatic <(samtools mpileup -l /data/refseq/hglft_genome_304d_b78af0.bed --no-BAQ -f /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa /data/alignment/final/Exome_Norm_sorted_mrkdup_bqsr.bam /data/alignment/final/Exome_Tumor_sorted_mrkdup_bqsr.bam) /data/varscan/exome --mpileup 1 --output-vcf`
* `cd /data/varscan/`
* `java -Xmx4g -jar VarScan.v2.4.2.jar processSomatic exome.snp.vcf exome.snp`
* `java -Xmx4g -jar VarScan.v2.4.2.jar processSomatic exome.indel.vcf exome.indel`
* `find /data/varscan -name '*.vcf' -exec bgzip -f {} \;`
* `find /data/varscan -name '*.vcf.gz' -exec tabix -f {} \;`

In order to continue to the next step, you may need to redownload GATK if you run into errors with your current installation:
1. Manually download GATK after accepting the license : `https://www.broadinstitute.org/gatk/download/auth?package=GATK`
2. Copy the download to your instance:
`scp -i <your pem file> Downloads/GenomeAnalysisTK-3.6.tar.bz2 ubuntu@<IP address of your instance>:/data/bin`.
3. Inside the instance, unzip the archive : `tar --bzip2 -xvf GenomeAnalysisTK-3.6.tar.bz2`
4. Test with JAVA8 : `java -jar GenomeAnalysisTK.jar -h`

Next steps:
* `java -Xmx4g -jar GenomeAnalysisTK.jar -T VariantFiltration -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa --variant exome.snp.Somatic.vcf.gz --mask exome.snp.Somatic.hc.vcf.gz --maskName "processSomatic" --filterNotInMask -o exome.snp.Somatic.hc.filter.vcf.gz`
* `java -Xmx4g -jar GenomeAnalysisTK.jar -T VariantFiltration -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa --variant exome.indel.Somatic.vcf.gz --mask exome.indel.Somatic.hc.vcf.gz --maskName "processSomatic" --filterNotInMask -o exome.indel.Somatic.hc.filter.vcf.gz`

To continue, you will need to have bcftools installed, please refer to the installation section above if needed.

* `bcftools concat -a -o exome.vcf.gz -O z exome.snp.Somatic.hc.filter.vcf.gz exome.indel.Somatic.hc.filter.vcf.gz`
* `tabix -f /data/varscan/exome.vcf.gz`

#### **Running STRELKA**
__________________________  
Given that you have STRELKA properly installed, here are the commands for running STRELKA in order, you may need to adjust according to how you have named and placed your directories:
* `/data/strelka/bin/strelka-2.7.1.centos5_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam=/data/alignment/final/Exome_Norm_sorted_mrkdup_bqsr.bam --tumorBam=/data/alignment/final/Exome_Tumor_sorted_mrkdup_bqsr.bam --referenceFasta=../reference/GRCh38_full_analysis_set_plus_decoy_hla.fa --exome --runDir=/data/strelka/exome`
* `/data/strelka/exome/runWorkflow.py -m local -j 4` (Please specify according to the number of cpus avaliable. In this case I have 4 for my instance)
*	`zcat exome/results/variants/somatic.snvs.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > exome/results/variants/somatic.snvs.gt.vcf`
* `zcat exome/results/variants/somatic.indels.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > exome/results/variants/somatic.indels.gt.vcf`
* `find exome/results/variants -name *.vcf -exec bgzip -f {} \;`
* `find exome/results/variants -name *.vcf.gz -exec tabix -f {}`

After installing BCFtools:
* `bcftools concat -a -o exome.vcf.gz -O z exome/results/variants/somatic.snvs.gt.vcf.gz exome/results/variants/somatic.indels.gt.vcf.gz`
* `tabix exome.vcf.gz`

#### **Running MuTect2**
__________________________
Before preceeding, you will need to download COSMIC reference file (e.g. in folder `/data/refseq/`) after registering on their website:
* `sftp <your registered email address>@sftp-cancer.sanger.ac.uk`
* `get /files/grch38/cosmic/v79/VCF/CosmicCodingMuts.vcf.gz`
* `get /files/grch38/cosmic/v79/VCF/CosmicNonCodingVariants.vcf.gz`
* `Quit`
* `zgrep "^#" CosmicCodingMuts.vcf.gz > VCF_Header`
* `zgrep -v "^#" CosmicCodingMuts.vcf.gz | awk '{print "chr"$0}' | sed 's/^chrMT/chrM/' > CosmicCodingMuts.clean`
* `zgrep -v "^#" CosmicNonCodingVariants.vcf.gz | awk '{print "chr"$0}' | sed 's/^chrMT/chrM/' > CosmicNonCodingVariants.clean`
* `cat CosmicCodingMuts.clean CosmicNonCodingVariants.clean | sort -gk 2,2 > Cosmic_v79`
* `cat VCF_Header Cosmic_v79 > Cosmic_v79.vcf`
* `rm VCF_Header CosmicCodingMuts.clean CosmicNonCodingVariants.clean Cosmic_v79`

With picard tools installed:
* `java -jar picard.jar SortVcf I=/data/refseq/Cosmic_v79.vcf O=/data/refseq/Cosmic_v79.dictsorted.vcf SEQUENCE_DICTIONARY=/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa`
* `java -Xmx4g -jar GenomeAnalysisTK.jar -T MuTect2 --disable_auto_index_creation_and_locking_when_reading_rods -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I:tumor /data/alignment/final/Exome_Tumor_sorted_mrkdup_bqsr.bam -I:Normal /data/alignment/final/Exome_Norm_sorted_mrkdup_bqsr.bam --dbsnp /data/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz --cosmic /data/refseq/Cosmic_v79.dictsorted.vcf.gz -o /data/mutect/exome.vcf.gz -L /data/refseq/NimbleGenExome_v3.interval_list`
* `echo /data/mutect/exome.vcf.gz > /data/mutect/exome_vcf.fof`
* `bcftools concat --allow-overlaps --remove-duplicates --file-list /data/mutect/exome_vcf.fof --output-type z --output /data/mutect/exome.vcf.gz`
* `tabix /data/mutect/exome.vcf.gz`

#### **Merge Variants**
__________________________
With outputs from all three algorithms, we can now merge the variants to generate a comprehensive list of detected variants:
* `java -Xmx4g -jar GenomeAnalysisTK.jar -T CombineVariants -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -genotypeMergeOptions UNIQUIFY --variant:varscan /data/varscan/exome.vcf.gz --variant:strelka /data/strelka/exome.vcf.gz --variant:mutect /data/mutect/exome.vcf.gz -o /data/exome.unique.vcf.gz`
* `java -Xmx4g -jar GenomeAnalysisTK.jar -T CombineVariants -R /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -genotypeMergeOptions PRIORITIZE --rod_priority_list mutect,varscan,strelka --variant:varscan /data/varscan/exome.vcf.gz --variant:strelka strelka/exome.vcf.gz --variant:mutect /data/mutect/exome.vcf.gz -o /data/exome.merged.vcf`

**Please continue to the next section for instructions on how to filter, annotation and review variants**
