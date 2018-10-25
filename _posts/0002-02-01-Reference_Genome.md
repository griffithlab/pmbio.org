---
feature_text: |
  ## Precision Medicine
title: Reference Genome
categories:
    - Module-02-Inputs
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-02-01
---

### Obtain a reference genome

We will use the 1000 genomes version of the human GRCh38 build. This reference includes extra decoy and HLA sequences in addition to the alternate haplotypes provided from the GRC consortium. The 1000 genomes project is one of several places that people routinely obtain human reference genome files. Some additional sources including those that host many non-human reference genomes are described later in this section.

We obtained the original reference genome files from the 1000 genomes FTP site here:
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/

We have created a copy of these files on our course file server.  Furthermore, we have created a smaller version of the reference to allow us to complete this analysis more quickly.  Using the whole reference genome would take too long for a workshop setting.  For example, aligning reads from a single lange of whole genome data to the whole reference genome can take several hours.

For this course we have selected two chromosomes: chr6 and chr17.  We chose these two chromosomes to illustrate fundamentals of bioinformatics analysis efficiently but also because of the significance of these two chromosomes to cancer biology.  Why are chr6 and chr17 particularly relevant to cancer?

Download the genome reference files for this course using the following commands. Note use of an environment variable `CHRS` to specify the custom reference genome we are using here.

```bash
echo `CHRS` #If this doesn't give a value, please return to the Environment section of the course

mkdir -p /workspace/references/genome
cd /workspace/references/genome

# dowload human reference genome files from the course data server
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.dict
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/20150713_location_of_centromeres_and_other_regions.txt
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla-extra.fa
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/README.20150309.GRCh38_full_analysis_set_plus_decoy_hla
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
```

### Obtain Additional GATK resources needed

```bash
cd ~/data/reference/

#SNP calibration call sets - dbsnp, hapmap, omni, and 1000G
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf ~/workspace/data/raw_data/references
bgzip ~/workspace/data/raw_data/references/Homo_sapiens_assembly38.dbsnp138.vcf
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz ~/workspace/data/raw_data/references
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz ~/workspace/data/raw_data/references
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz ~/workspace/data/raw_data/references

#Indel calibration call sets - dbsnp, Mills
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz ~/workspace/data/raw_data/references
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ~/workspace/data/raw_data/references

gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list ~/workspace/data/raw_data/references
gsutil cp -r gs://genomics-public-data/resources/broad/hg38/v0/scattered_calling_intervals/ ~/workspace/data/raw_data/references

```

### Reference Genome Options

To do: Create a table documenting key reference genome options (builds/sources) and pros/cons
- e.g., 1000G, ensembl, UCSC, GDC
