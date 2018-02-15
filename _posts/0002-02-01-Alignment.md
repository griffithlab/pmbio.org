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

| Bam | Readgroup ID | Platform | FC/BC/Lane | Library | Sample Name |
|-----|--------------|----------|------------|---------|-------------|
| Exome_Norm_gerald_C1TD1ACXX_7_CGATGT.bam | ID:2891351068 | PL:Illumina | PU:C1TD1ACXX-CGATGT.7 | libgroup-2891242742 | SM:H_NJ-HCC1395-HCC1395_BL |
| Exome_Tumor_gerald_C1TD1ACXX_7_ATCACG.bam | ID:2891351066 | PL:Illumina | PU:C1TD1ACXX-ATCACG.7 | libgroup-2891242741 | SM:H_NJ-HCC1395-HCC1395 |
| H_NJ-HCC1395-HCC1395_BL_RNA.bam | ID:2895625992 | PL:Illumina | PU:H3MYFBBXX-CTTGTA.4 | Pooled_RNA_2891006726-cD1-lib1 | SM:H_NJ-HCC1395-HCC1395_BL_RNA |
| H_NJ-HCC1395-HCC1395_BL_RNA.bam | ID:2895626097 | PL:Illumina | PU:H3MYFBBXX-CTTGTA.5 | Pooled_RNA_2891006726-cD1-lib1 | SM:H_NJ-HCC1395-HCC1395_BL_RNA |
| H_NJ-HCC1395-HCC1395_RNA.bam | ID:2895626112 | PL:Illumina | PU:H3MYFBBXX-GCCAAT.5 | Pooled_RNA_2891007020-cD1-lib1 | SM:H_NJ-HCC1395-HCC1395_RNA |
| H_NJ-HCC1395-HCC1395_RNA.bam | ID:2895626107 | PL:Illumina | PU:H3MYFBBXX-GCCAAT.4 | Pooled_RNA_2891007020-cD1-lib1 | SM:H_NJ-HCC1395-HCC1395_RNA |
| WGS_Norm_Lane1_gerald_D1VCPACXX_6.bam | ID:2891323123 | PL:Illumina | PU:D1VCPACXX.6 | H_NJ-HCC1395-HCC1395_BL-lig2-lib1 | SM:H_NJ-HCC1395-HCC1395_BL |
| WGS_Norm_Lane2_gerald_D1VCPACXX_7.bam | ID:2891323124 | PL:Illumina | PU:D1VCPACXX.7 | H_NJ-HCC1395-HCC1395_BL-lig2-lib2 | SM:H_NJ-HCC1395-HCC1395_BL |
| WGS_Norm_Lane3_gerald_D1VCPACXX_8.bam | ID:2891323125 | PL:Illumina | PU:D1VCPACXX.8 | H_NJ-HCC1395-HCC1395_BL-lig2-lib3 | SM:H_NJ-HCC1395-HCC1395_BL |
| WGS_Tumor_Lane1_gerald_D1VCPACXX_1.bam | ID:2891322951 | PL:Illumina | PU:D1VCPACXX.1 | H_NJ-HCC1395-HCC1395-lig2-lib1 | SM:H_NJ-HCC1395-HCC1395 |
| WGS_Tumor_Lane2_gerald_D1VCPACXX_2.bam | ID:2891323174 | PL:Illumina | PU:D1VCPACXX.2 | H_NJ-HCC1395-HCC1395-lig2-lib1 | SM:H_NJ-HCC1395-HCC1395 |
| WGS_Tumor_Lane3_gerald_D1VCPACXX_3.bam | ID:2891323175 | PL:Illumina | PU:D1VCPACXX.3 | H_NJ-HCC1395-HCC1395-lig2-lib2 | SM:H_NJ-HCC1395-HCC1395 |
| WGS_Tumor_Lane4_gerald_D1VCPACXX_4.bam | ID:2891323150 | PL:Illumina | PU:D1VCPACXX.4 | H_NJ-HCC1395-HCC1395-lig2-lib2 | SM:H_NJ-HCC1395-HCC1395 |
| WGS_Tumor_Lane5_gerald_D1VCPACXX_5.bam | ID:2891323147 | PL:Illumina | PU:D1VCPACXX.5 | H_NJ-HCC1395-HCC1395-lig2-lib3 | SM:H_NJ-HCC1395-HCC1395 |
|-----|--------------|----------|------------|---------|-------------|


### Run bwa mem using the above information

```bash
cd ~/data
mkdir alignment
bwa mem -t 32 -Y -R "@RG\tID:2891351068\tPL:ILLUMINA\tPU:C1TD1ACXX-CGATGT.7\tLB:libgroup-2891242742\tSM:H_NJ-HCC1395-HCC1395_BL" -o /home/ubuntu/data/alignment/2891351068.bam /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/data/fastqs/Exome_Norm/2891351068_1.fastq.gz /home/ubuntu/data/fastqs/Exome_Norm/2891351068_2.fastq.gz 


```



