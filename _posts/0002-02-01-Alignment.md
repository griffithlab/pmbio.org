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

| Data | Readgroup ID | Platform | FC[-BC].Lane | Library | Sample Name |
|-----|--------------|----------|------------|---------|-------------|
| Exome_Norm | 2891351068 | Illumina | C1TD1ACXX-CGATGT.7 | exome_norm_lib1 | HCC1395BL_DNA |
| Exome_Tumor | 2891351066 | Illumina | C1TD1ACXX-ATCACG.7 | exome_tumor_lib1 | HCC1395_DNA |
| RNAseq_Norm_Lane1 | 2895625992 | Illumina | H3MYFBBXX-CTTGTA.4 | rna_norm_lib1 | HCC1395BL_RNA |
| RNAseq_Norm_Lane2 | 2895626097 | Illumina | H3MYFBBXX-CTTGTA.5 | rna_norm_lib1 | HCC1395BL_RNA |
| RNAseq_Tumor_Lane1 | 2895626107 | Illumina | H3MYFBBXX-GCCAAT.4 | rna_tumor_lib1 | HCC1395_RNA |
| RNAseq_Tumor_Lane2 | 2895626112 | Illumina | H3MYFBBXX-GCCAAT.5 | rna_tumor_lib1 | HCC1395_RNA |
| WGS_Norm_Lane1 | 2891323123 | Illumina | D1VCPACXX.6 | wgs_norm_lib1 | HCC1395BL_DNA |
| WGS_Norm_Lane2 | 2891323124 | Illumina | D1VCPACXX.7 | wgs_norm_lib2 | HCC1395BL_DNA |
| WGS_Norm_Lane3 | 2891323125 | Illumina | D1VCPACXX.8 | wgs_norm_lib3 | HCC1395BL_DNA |
| WGS_Tumor_Lane1 | 2891322951 | Illumina | D1VCPACXX.1 | wgs_tumor_lib1 | HCC1395_DNA |
| WGS_Tumor_Lane2 | 2891323174 | Illumina | D1VCPACXX.2 | wgs_tumor_lib1 | HCC1395_DNA |
| WGS_Tumor_Lane3 | 2891323175 | Illumina | D1VCPACXX.3 | wgs_tumor_lib2 | HCC1395_DNA |
| WGS_Tumor_Lane4 | 2891323150 | Illumina | D1VCPACXX.4 | wgs_tumor_lib2 | HCC1395_DNA |
| WGS_Tumor_Lane5 | 2891323147 | Illumina | D1VCPACXX.5 | wgs_tumor_lib3 | HCC1395_DNA |
|-----|--------------|----------|------------|---------|-------------|

### Run bwa mem using the above information

Runtimes: Exome, 18 - 23min; WGS, 89-95min 

```bash
cd ~/data
mkdir alignment
cd alignment
bwa mem -t 32 -Y -R "@RG\tID:2891351068\tPL:ILLUMINA\tPU:C1TD1ACXX-CGATGT.7\tLB:exome_norm_lib1\tSM:HCC1395BL_DNA" -o /home/ubuntu/data/alignment/Exome_Norm_2891351068.sam /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/data/fastqs/Exome_Norm/2891351068_1.fastq.gz /home/ubuntu/data/fastqs/Exome_Norm/2891351068_2.fastq.gz
bwa mem -t 32 -Y -R "@RG\tID:2891351066\tPL:ILLUMINA\tPU:C1TD1ACXX-ATCACG.7\tLB:exome_tumor_lib1\tSM:HCC1395_DNA" -o /home/ubuntu/data/alignment/Exome_Tumor_2891351066.sam /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/data/fastqs/Exome_Tumor/2891351066_1.fastq.gz /home/ubuntu/data/fastqs/Exome_Tumor/2891351066_2.fastq.gz
bwa mem -t 32 -Y -R "@RG\tID:2891323123\tPL:ILLUMINA\tPU:D1VCPACXX.6\tLB:wgs_norm_lib1\tSM:HCC1395BL_DNA" -o /home/ubuntu/data/alignment/WGS_Norm_Lane1_2891323123.sam /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/data/fastqs/WGS_Norm/2891323123_1.fastq.gz /home/ubuntu/data/fastqs/WGS_Norm/2891323123_2.fastq.gz
bwa mem -t 32 -Y -R "@RG\tID:2891323124\tPL:ILLUMINA\tPU:D1VCPACXX.7\tLB:wgs_norm_lib2\tSM:HCC1395BL_DNA" -o /home/ubuntu/data/alignment/WGS_Norm_Lane2_2891323124.sam /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/data/fastqs/WGS_Norm/2891323124_1.fastq.gz /home/ubuntu/data/fastqs/WGS_Norm/2891323124_2.fastq.gz
bwa mem -t 32 -Y -R "@RG\tID:2891323125\tPL:ILLUMINA\tPU:D1VCPACXX.8\tLB:wgs_norm_lib3\tSM:HCC1395BL_DNA" -o /home/ubuntu/data/alignment/WGS_Norm_Lane3_2891323125.sam /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/data/fastqs/WGS_Norm/2891323125_1.fastq.gz /home/ubuntu/data/fastqs/WGS_Norm/2891323125_2.fastq.gz
bwa mem -t 32 -Y -R "@RG\tID:2891322951\tPL:ILLUMINA\tPU:D1VCPACXX.1\tLB:wgs_tumor_lib1\tSM:HCC1395_DNA" -o /home/ubuntu/data/alignment/WGS_Tumor_Lane1_2891322951.sam /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/data/fastqs/WGS_Tumor/2891322951_1.fastq.gz /home/ubuntu/data/fastqs/WGS_Tumor/2891322951_2.fastq.gz
bwa mem -t 32 -Y -R "@RG\tID:2891323174\tPL:ILLUMINA\tPU:D1VCPACXX.2\tLB:wgs_tumor_lib1\tSM:HCC1395_DNA" -o /home/ubuntu/data/alignment/WGS_Tumor_Lane2_2891323174.sam /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/data/fastqs/WGS_Tumor/2891323174_1.fastq.gz /home/ubuntu/data/fastqs/WGS_Tumor/2891323174_2.fastq.gz
bwa mem -t 32 -Y -R "@RG\tID:2891323175\tPL:ILLUMINA\tPU:D1VCPACXX.3\tLB:wgs_tumor_lib2\tSM:HCC1395_DNA" -o /home/ubuntu/data/alignment/WGS_Tumor_Lane3_2891323175.sam /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/data/fastqs/WGS_Tumor/2891323175_1.fastq.gz /home/ubuntu/data/fastqs/WGS_Tumor/2891323175_2.fastq.gz
bwa mem -t 32 -Y -R "@RG\tID:2891323150\tPL:ILLUMINA\tPU:D1VCPACXX.4\tLB:wgs_tumor_lib2\tSM:HCC1395_DNA" -o /home/ubuntu/data/alignment/WGS_Tumor_Lane4_2891323150.sam /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/data/fastqs/WGS_Tumor/2891323150_1.fastq.gz /home/ubuntu/data/fastqs/WGS_Tumor/2891323150_2.fastq.gz
bwa mem -t 32 -Y -R "@RG\tID:2891323147\tPL:ILLUMINA\tPU:D1VCPACXX.5\tLB:wgs_tumor_lib3\tSM:HCC1395_DNA" -o /home/ubuntu/data/alignment/WGS_Tumor_Lane5_2891323147.sam /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/data/fastqs/WGS_Tumor/2891323147_1.fastq.gz /home/ubuntu/data/fastqs/WGS_Tumor/2891323147_2.fastq.gz
```

### Convert sam to cram format

Runtimes: Exome, 13min;

```bash
cd ~/data/alignment
samtools view -h -C -T /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -o Exome_Norm_2891351068.cram Exome_Norm_2891351068.sam
```

### Sort cram file

```bash
cd ~/data/alignment
samtools sort -o Exome_Norm_2891351068.sorted.cram Exome_Norm_2891351068.cram
```

### Merge cram files for WGS data


### Index cram file

```bash
cd ~/data/alignment
samtools index Exome_Norm_2891351068.sorted.cram
```

### Mark duplicates



