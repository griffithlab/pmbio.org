---
feature_text: |
  ## Precision Medicine
title: Alignment
categories:
    - Module 2
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-02-01
---

WGS and Exome fastq data will be aligned with BWA MEM using the following options:

- t: number of threads
- Y: use soft clipping for supplementary alignments
- R: read group header line info. See table below for sample details.

| Data | Readgroup ID | Platform | FC/BC/Lane | Library | Sample Name |
|-----|--------------|----------|------------|---------|-------------|
| Exome_Norm | 2891351068 | Illumina | C1TD1ACXX-CGATGT.7 | lib1 | HCC1395BL_DNA |
| Exome_Tumor | 2891351066 | Illumina | C1TD1ACXX-ATCACG.7 | lib1 | HCC1395_DNA |
| RNAseq_Norm_Lane1 | 2895625992 | Illumina | H3MYFBBXX-CTTGTA.4 | lib1 | HCC1395BL_RNA |
| RNAseq_Norm_Lane2 | 2895626097 | Illumina | H3MYFBBXX-CTTGTA.5 | lib1 | HCC1395BL_RNA |
| RNAseq_Tumor_Lane1 | 2895626107 | Illumina | H3MYFBBXX-GCCAAT.4 | lib1 | HCC1395_RNA |
| RNAseq_Tumor_Lane2 | 2895626112 | Illumina | H3MYFBBXX-GCCAAT.5 | lib1 | HCC1395_RNA |
| WGS_Norm_Lane1 | 2891323123 | Illumina | D1VCPACXX.6 | lib1 | HCC1395BL_DNA |
| WGS_Norm_Lane2 | 2891323124 | Illumina | D1VCPACXX.7 | lib2 | HCC1395BL_DNA |
| WGS_Norm_Lane3 | 2891323125 | Illumina | D1VCPACXX.8 | lib3 | HCC1395BL_DNA |
| WGS_Tumor_Lane1 | 2891322951 | Illumina | D1VCPACXX.1 | lib1 | HCC1395_DNA |
| WGS_Tumor_Lane2 | 2891323174 | Illumina | D1VCPACXX.2 | lib1 | HCC1395_DNA |
| WGS_Tumor_Lane3 | 2891323175 | Illumina | D1VCPACXX.3 | lib2 | HCC1395_DNA |
| WGS_Tumor_Lane4 | 2891323150 | Illumina | D1VCPACXX.4 | lib2 | HCC1395_DNA |
| WGS_Tumor_Lane5 | 2891323147 | Illumina | D1VCPACXX.5 | lib3 | HCC1395_DNA |
|-----|--------------|----------|------------|---------|-------------|


### Run bwa mem using the above information

Runtimes: Exome, 16min; 

```bash
cd ~/data
mkdir alignment
bwa mem -t 32 -Y -R "@RG\tID:2891351068\tPL:ILLUMINA\tPU:C1TD1ACXX-CGATGT.7\tLB:lib1\tSM:HCC1395BL_DNA" -o /home/ubuntu/data/alignment/Exome_Norm_2891351068.sam /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/data/fastqs/Exome_Norm/2891351068_1.fastq.gz /home/ubuntu/data/fastqs/Exome_Norm/2891351068_2.fastq.gz
```

### Convert sam to cram format

Runtimes: Exome, 13min;

```bash
samtools view -h -C -T /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -o Exome_Norm_2891351068.cram Exome_Norm_2891351068.sam
```

### Sort cram file




### Mark duplicates



