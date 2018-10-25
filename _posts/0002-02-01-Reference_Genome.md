---
feature_text: |
  ## Precision Medicine
title: Reference Genome
categories:
    - Module-02-Inputs
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-02-01
---

First let’s go over what a reference assembly actually is, in essence it’s just a representation of the nucleotide sequence from a cohort. These assemblies allow for a shortcut when mapping reads as they can be mapped to the assembly, rather than each other, to piece the genome of an individual together. This has a number of benefits, the most obvious of which is that it is far more effecient than attempting to build a genome from scratch. By its very nature however using this approach means there is no perfect reference assembly for an individual due to polymorphism (i.e. snps, hla-type, etc.). Further due to the presence of repetitive structural elements such as duplications, inverted repeats, tandem repeats, etc. a given assembly is almost always incomplete, and is constantly being improved upon. This leads to the publication of new assembly versions every so often such as grch37 (Feb. 2009) and grch38 (Dec. 2013). It is also good to be aware that different organization can publish different reference assemblies, for example grch37 (NCBI) and hg19 (UCSC) are identical save for a few minor differences such as in the mitochondria sequence and annotation of chromosomes (1 vs chr1). For a nice summary of genome versions and their release names refer to the [Assembly Releases and Versions FAQ](http://genome.ucsc.edu/FAQ/FAQreleases.html).

### Obtain a reference genome

We will use the 1000 genomes version of the human GRCh38 build. This reference includes extra decoy and HLA sequences in addition to the alternate haplotypes provided from the GRC consortium. We should also download the dictionary (.dict) file, location of centromeres file, extra sequence fasta and README file for later reference.

- GRCh38_full_analysis_set_plus_decoy_hla.fa
- GRCh38_full_analysis_set_plus_decoy_hla.dict
- GRCh38_full_analysis_set_plus_decoy_hla-extra.fa
- 20150713_location_of_centromeres_and_other_regions.txt
- README.20150309.GRCh38_full_analysis_set_plus_decoy_hla

```bash
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
