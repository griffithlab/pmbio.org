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

We will use the 1000 genomes version of the human GRCh38 build. This reference includes extra decoy and HLA sequences in addition to the alternate haplotypes provided from the GRC consortium. We should also download the dictionary (.dict) file, location of centromeres file, extra sequence fasta and README file for later reference.

- GRCh38_full_analysis_set_plus_decoy_hla.fa
- GRCh38_full_analysis_set_plus_decoy_hla.dict
- GRCh38_full_analysis_set_plus_decoy_hla-extra.fa
- 20150713_location_of_centromeres_and_other_regions.txt
- README.20150309.GRCh38_full_analysis_set_plus_decoy_hla

```bash
cd /workspace/data/raw_data/references
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
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf /home/ubuntu/data/reference
bgzip /home/ubuntu/data/reference/Homo_sapiens_assembly38.dbsnp138.vcf
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz /home/ubuntu/data/reference
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz /home/ubuntu/data/reference
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz /home/ubuntu/data/reference

#Indel calibration call sets - dbsnp, Mills
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz /home/ubuntu/data/reference
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz /home/ubuntu/data/reference

gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list /home/ubuntu/data/reference
gsutil cp -r gs://genomics-public-data/resources/broad/hg38/v0/scattered_calling_intervals/ /home/ubuntu/data/reference

```

### Reference Genome Options

To do: Create a table documenting key reference genome options (builds/sources) and pros/cons
- e.g., 1000G, ensembl, UCSC, GDC
