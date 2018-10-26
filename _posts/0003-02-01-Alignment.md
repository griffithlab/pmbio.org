---
feature_text: |
  ## Precision Medicine
title: Alignment
categories:
    - Module-03-Align
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-02-01
---

WGS and Exome raw sequence (fastq) data for our tumor and normal samples will be aligned with BWA MEM using the following options:

- t: number of threads
- Y: use soft clipping for supplementary alignments
- R: read group header line info. See table in [Data Module](/module-02-inputs/0002/05/01/Data/) for sample details.

### Run bwa mem using the above information

```bash
cd /workspace/align
bwa mem -t 8 -Y -R "@RG\tID:Exome_Norm\tPL:ILLUMINA\tPU:C1TD1ACXX-CGATGT.7\tLB:exome_norm_lib1\tSM:HCC1395BL_DNA" -o /workspace/align/Exome_Norm.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/Exome_Norm/Exome_Norm_R1.fastq.gz /workspace/inputs/data/fastq/Exome_Norm/Exome_Norm_R2.fastq.gz
bwa mem -t 8 -Y -R "@RG\tID:Exome_Tumor\tPL:ILLUMINA\tPU:C1TD1ACXX-ATCACG.7\tLB:exome_tumor_lib1\tSM:HCC1395_DNA" -o /workspace/align/Exome_Tumor.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/Exome_Tumor/Exome_Tumor_R1.fastq.gz /workspace/inputs/data/fastq/Exome_Tumor/Exome_Tumor_R2.fastq.gz

bwa mem -t 8 -Y -R "@RG\tID:WGS_Norm_Lane1\tPL:ILLUMINA\tPU:D1VCPACXX.6\tLB:wgs_norm_lib1\tSM:HCC1395BL_DNA" -o /workspace/align/WGS_Norm_Lane1.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/WGS_Norm/WGS_Norm_Lane1_R1.fastq.gz /workspace/inputs/data/fastq/WGS_Norm/WGS_Norm_Lane1_R2.fastq.gz
bwa mem -t 8 -Y -R "@RG\tID:WGS_Norm_Lane2\tPL:ILLUMINA\tPU:D1VCPACXX.7\tLB:wgs_norm_lib2\tSM:HCC1395BL_DNA" -o /workspace/align/WGS_Norm_Lane2.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/WGS_Norm/WGS_Norm_Lane2_R1.fastq.gz /workspace/inputs/data/fastq/WGS_Norm/WGS_Norm_Lane2_R2.fastq.gz
bwa mem -t 8 -Y -R "@RG\tID:WGS_Norm_Lane3\tPL:ILLUMINA\tPU:D1VCPACXX.8\tLB:wgs_norm_lib3\tSM:HCC1395BL_DNA" -o /workspace/align/WGS_Norm_Lane3.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/WGS_Norm/WGS_Norm_Lane3_R1.fastq.gz /workspace/inputs/data/fastq/WGS_Norm/WGS_Norm_Lane3_R2.fastq.gz

bwa mem -t 8 -Y -R "@RG\tID:WGS_Tumor_Lane1\tPL:ILLUMINA\tPU:D1VCPACXX.1\tLB:wgs_tumor_lib1\tSM:HCC1395_DNA" -o /workspace/align/WGS_Tumor_Lane1.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane1_R1.fastq.gz /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane1_R2.fastq.gz
bwa mem -t 8 -Y -R "@RG\tID:WGS_Tumor_Lane2\tPL:ILLUMINA\tPU:D1VCPACXX.2\tLB:wgs_tumor_lib1\tSM:HCC1395_DNA" -o /workspace/align/WGS_Tumor_Lane2.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane2_R1.fastq.gz /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane2_R2.fastq.gz
bwa mem -t 8 -Y -R "@RG\tID:WGS_Tumor_Lane3\tPL:ILLUMINA\tPU:D1VCPACXX.3\tLB:wgs_tumor_lib2\tSM:HCC1395_DNA" -o /workspace/align/WGS_Tumor_Lane3.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane3_R1.fastq.gz /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane3_R2.fastq.gz
bwa mem -t 8 -Y -R "@RG\tID:WGS_Tumor_Lane4\tPL:ILLUMINA\tPU:D1VCPACXX.4\tLB:wgs_tumor_lib2\tSM:HCC1395_DNA" -o /workspace/align/WGS_Tumor_Lane4.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane4_R1.fastq.gz /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane4_R2.fastq.gz
bwa mem -t 8 -Y -R "@RG\tID:WGS_Tumor_Lane5\tPL:ILLUMINA\tPU:D1VCPACXX.5\tLB:wgs_tumor_lib3\tSM:HCC1395_DNA" -o /workspace/align/WGS_Tumor_Lane5.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane5_R1.fastq.gz /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane5_R2.fastq.gz
```

### Convert sam to bam format

```bash
cd /workspace/align
samtools view -@ 8 -h -b -o Exome_Norm.bam Exome_Norm.sam
samtools view -@ 8 -h -b -o Exome_Tumor.bam Exome_Tumor.sam
samtools view -@ 8 -h -b -o WGS_Norm_Lane1.bam WGS_Norm_Lane1.sam
samtools view -@ 8 -h -b -o WGS_Norm_Lane2.bam WGS_Norm_Lane2.sam
samtools view -@ 8 -h -b -o WGS_Norm_Lane3.bam WGS_Norm_Lane3.sam
samtools view -@ 8 -h -b -o WGS_Tumor_Lane1.bam WGS_Tumor_Lane1.sam
samtools view -@ 8 -h -b -o WGS_Tumor_Lane2.bam WGS_Tumor_Lane2.sam
samtools view -@ 8 -h -b -o WGS_Tumor_Lane3.bam WGS_Tumor_Lane3.sam
samtools view -@ 8 -h -b -o WGS_Tumor_Lane4.bam WGS_Tumor_Lane4.sam
samtools view -@ 8 -h -b -o WGS_Tumor_Lane5.bam WGS_Tumor_Lane5.sam
```

### Merge bam files

For the WGS data, we have multiple separate bams for each lane of data. Let's take this opportunity to merge them into a single bam for tumor and normal respectively.

```bash
cd /workspace/align
samtools merge -@ 8 WGS_Norm_merged.bam WGS_Norm_Lane1.bam WGS_Norm_Lane2.bam WGS_Norm_Lane3.bam
samtools merge -@ 8 WGS_Tumor_merged.bam WGS_Tumor_Lane1.bam WGS_Tumor_Lane2.bam WGS_Tumor_Lane3.bam WGS_Tumor_Lane4.bam WGS_Tumor_Lane5.bam
```


### Query name sort bam files

```bash
cd /workspace/align
java -Xmx60g -jar $PICARD SortSam I=Exome_Norm.bam O=Exome_Norm_namesorted.bam SO=queryname
java -Xmx60g -jar $PICARD SortSam I=Exome_Tumor.bam O=Exome_Tumor_namesorted.bam SO=queryname
java -Xmx60g -jar $PICARD SortSam I=WGS_Norm_merged.bam O=WGS_Norm_merged_namesorted.bam SO=queryname
java -Xmx60g -jar $PICARD SortSam I=WGS_Tumor_merged.bam O=WGS_Tumor_merged_namesorted.bam SO=queryname
```


### Mark duplicates

```bash
cd /workspace/align
java -Xmx64g -jar $PICARD MarkDuplicates I=Exome_Norm_namesorted.bam O=Exome_Norm_namesorted_mrkdup.bam ASSUME_SORT_ORDER=queryname METRICS_FILE=Exome_Norm_mrkdup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT
java -Xmx64g -jar $PICARD MarkDuplicates I=Exome_Tumor_namesorted.bam O=Exome_Tumor_namesorted_mrkdup.bam ASSUME_SORT_ORDER=queryname METRICS_FILE=Exome_Tumor_mrkdup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT
java -Xmx64g -jar $PICARD MarkDuplicates I=WGS_Norm_merged_namesorted.bam O=WGS_Norm_merged_namesorted_mrkdup.bam ASSUME_SORT_ORDER=queryname METRICS_FILE=WGS_Norm_mrkdup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT
java -Xmx64g -jar $PICARD MarkDuplicates I=WGS_Tumor_merged_namesorted.bam O=WGS_Tumor_merged_namesorted_mrkdup.bam ASSUME_SORT_ORDER=queryname METRICS_FILE=WGS_Tumor_mrkdup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT
```


### Position sort bam file

```bash
cd /workspace/align
java -Xmx64g -jar $PICARD SortSam I=Exome_Norm_namesorted_mrkdup.bam O=Exome_Norm_sorted_mrkdup.bam SO=coordinate
java -Xmx64g -jar $PICARD SortSam I=Exome_Tumor_namesorted_mrkdup.bam O=Exome_Tumor_sorted_mrkdup.bam SO=coordinate
java -Xmx64g -jar $PICARD SortSam I=WGS_Norm_merged_namesorted_mrkdup.bam O=WGS_Norm_merged_sorted_mrkdup.bam SO=coordinate
java -Xmx64g -jar $PICARD SortSam I=WGS_Tumor_merged_namesorted_mrkdup.bam O=WGS_Tumor_merged_sorted_mrkdup.bam SO=coordinate
```


### Create bam index for use with GATK, IGV, etc
Runtimes: Exome XXX; WGS 43-78min;

```bash
cd /workspace/align
java -Xmx64g -jar $PICARD BuildBamIndex I=Exome_Norm_sorted_mrkdup.bam
java -Xmx64g -jar $PICARD BuildBamIndex I=Exome_Tumor_sorted_mrkdup.bam
java -Xmx64g -jar $PICARD BuildBamIndex I=WGS_Norm_merged_sorted_mrkdup.bam
java -Xmx64g -jar $PICARD BuildBamIndex I=WGS_Tumor_merged_sorted_mrkdup.bam
```


### Perform Indel Realignment

If desired, add this step. See docs [here](https://software.broadinstitute.org/gatk/documentation/article?id=7156). But, note that as announced in the GATK v3.6 highlights, variant calling workflows that use HaplotypeCaller or MuTect2 now omit indel realignment. HaplotypeCaller includes a local read assembly that mostly deprecates/replaces the need for a separate indel realignment step. See the following [blog](https://software.broadinstitute.org/gatk/blog?id=7847) for a detailed discussion of this issue.

See [here](https://drive.google.com/drive/folders/1U6Zm_tYn_3yeEgrD1bdxye4SXf5OseIt) for latest versions of all GATK tutorials:


### Perform Base Quality Score Recalibration

Questions about GATK step.
Why that the BQSR commands below limit the modeling step to chr1-22. This is where the majority of known variants are located and the autosomes are expected to have more even coverage than sex chromosomes. However, once the model is built, we apply to all bases on all contigs.






#### Calculate BQSR Table


```bash
cd /workspace/align
gatk --java-options '-Xmx64g' BaseRecalibrator -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/Exome_Norm_sorted_mrkdup.bam -O /home/ubuntu/data/alignment/Exome_Norm_sorted_mrkdup_bqsr.table --known-sites /home/ubuntu/data/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites /home/ubuntu/data/reference/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites /home/ubuntu/data/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --preserve-qscores-less-than 6 --disable-bam-index-caching  -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22
gatk --java-options '-Xmx64g' BaseRecalibrator -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/Exome_Tumor_sorted_mrkdup.bam -O /home/ubuntu/data/alignment/Exome_Tumor_sorted_mrkdup_bqsr.table --known-sites /home/ubuntu/data/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites /home/ubuntu/data/reference/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites /home/ubuntu/data/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --preserve-qscores-less-than 6 --disable-bam-index-caching  -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22
gatk --java-options '-Xmx64g' BaseRecalibrator -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/WGS_Norm_merged_sorted_mrkdup.bam -O /home/ubuntu/data/alignment/WGS_Norm_merged_sorted_mrkdup_bqsr.table --known-sites /home/ubuntu/data/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites /home/ubuntu/data/reference/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites /home/ubuntu/data/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --preserve-qscores-less-than 6 --disable-bam-index-caching  -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22
gatk --java-options '-Xmx64g' BaseRecalibrator -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/WGS_Tumor_merged_sorted_mrkdup.bam -O /home/ubuntu/data/alignment/WGS_Tumor_merged_sorted_mrkdup_bqsr.table --known-sites /home/ubuntu/data/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites /home/ubuntu/data/reference/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites /home/ubuntu/data/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --preserve-qscores-less-than 6 --disable-bam-index-caching  -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22
```

#### Apply BQSR

Runtimes: Exome 39-48min; WGS 264-508min

```bash
cd /workspace/align
gatk --java-options '-Xmx64g' ApplyBQSR -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/Exome_Norm_sorted_mrkdup.bam -O /home/ubuntu/data/alignment/Exome_Norm_sorted_mrkdup_bqsr.bam --bqsr-recal-file /home/ubuntu/data/alignment/Exome_Norm_sorted_mrkdup_bqsr.table --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30
gatk --java-options '-Xmx64g' ApplyBQSR -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/Exome_Tumor_sorted_mrkdup.bam -O /home/ubuntu/data/alignment/Exome_Tumor_sorted_mrkdup_bqsr.bam --bqsr-recal-file /home/ubuntu/data/alignment/Exome_Tumor_sorted_mrkdup_bqsr.table --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30
gatk --java-options '-Xmx64g' ApplyBQSR -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/WGS_Norm_merged_sorted_mrkdup.bam -O /home/ubuntu/data/alignment/WGS_Norm_merged_sorted_mrkdup_bqsr.bam --bqsr-recal-file /home/ubuntu/data/alignment/WGS_Norm_merged_sorted_mrkdup_bqsr.table --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30
gatk --java-options '-Xmx64g' ApplyBQSR -R /home/ubuntu/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -I /home/ubuntu/data/alignment/WGS_Tumor_merged_sorted_mrkdup.bam -O /home/ubuntu/data/alignment/WGS_Tumor_merged_sorted_mrkdup_bqsr.bam --bqsr-recal-file /home/ubuntu/data/alignment/WGS_Tumor_merged_sorted_mrkdup_bqsr.table --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30
```


### Clean up un-needed sam/bam files

Keep final sorted, duplicated marked, bqsr bam/bai/table files and mrkdup.txt files. Delete everything else.

```bash
cd /workspace/align
mkdir final
mv *_sorted_mrkdup_bqsr.* final/
mv *.txt final/

rm *.sam
rm *.bam
rm *.bai

```

### Alignment Instructions for running on Chr6+Chr17

All of these may be missing WGS data processing
#### Run bwa mem
```bash
mkdir -p ~/workspace/data/raw_data/fastqs/
wget http://genomedata.org/pmbio-workshop/fastqs/chr6+chr17/Exome_Norm.tar
wget http://genomedata.org/pmbio-workshop/fastqs/chr6+chr17/Exome_Tumor.tar

mkdir -p ~/workspace/data/DNA_alignments/chr6+chr17
cd ~/workspace/data/DNA_alignments/chr6+chr17

bwa mem -t 4 -Y -R "@RG\tID:2891351068\tPL:ILLUMINA\tPU:C1TD1ACXX-CGATGT.7\tLB:exome_norm_lib1\tSM:HCC1395BL_DNA" -o ~/workspace/data/DNA_alignments/chr6+chr17/Exome_Norm.sam ~/workspace/data/raw_data/references/chr6+chr17.fa ~/workspace/data/raw_data/fastqs/chr6+chr17/Exome_Norm/Exome_Norm_R1.fastq.gz ~/workspace/data/raw_data/fastqs/chr6+chr17/Exome_Norm/Exome_Norm_R2.fastq.gz
bwa mem -t 4 -Y -R "@RG\tID:2891351066\tPL:ILLUMINA\tPU:C1TD1ACXX-ATCACG.7\tLB:exome_tumor_lib1\tSM:HCC1395_DNA" -o ~/workspace/data/DNA_alignments/chr6+chr17/Exome_Tumor.sam ~/workspace/data/raw_data/references/chr6+chr17.fa ~/workspace/data/raw_data/fastqs/chr6+chr17/Exome_Tumor/Exome_Tumor_R1.fastq.gz ~/workspace/data/raw_data/fastqs/chr6+chr17/Exome_Tumor/Exome_Tumor_R2.fastq.gz
```

#### Convert SAM to BAM format

```bash
cd ~/workspace/data/DNA_alignments/chr6+chr17
samtools view -h -b -o Exome_Norm.bam Exome_Norm.sam
samtools view -h -b -o Exome_Tumor.bam Exome_Tumor.sam
```
#### Merge Bam Files for WGS
```bash
cd ~/workspace/data/DNA_alignments/chr6+chr17
samtools merge -@ 4 WGS_Norm_merged.bam WGS_Norm_Lane1.bam WGS_Norm_Lane2.bam WGS_Norm_Lane3.bam
```
#### Query name sort bam files

```bash
cd ~/workspace/data/DNA_alignments/chr6+chr17
java -Xmx12g -jar $PICARD SortSam I=Exome_Norm.bam O=Exome_Norm_namesorted.bam SO=queryname
java -Xmx12g -jar $PICARD SortSam I=Exome_Tumor.bam O=Exome_Tumor_namesorted.bam SO=queryname
```
#### Mark Duplicates

```bash
cd ~/workspace/data/DNA_alignments/chr6+chr17
java -Xmx12g -jar $PICARD MarkDuplicates I=Exome_Norm_namesorted.bam O=Exome_Norm_namesorted_mrkdup.bam ASSUME_SORT_ORDER=queryname METRICS_FILE=Exome_Norm_mrkdup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT
java -Xmx12g -jar $PICARD MarkDuplicates I=Exome_Tumor_namesorted.bam O=Exome_Tumor_namesorted_mrkdup.bam ASSUME_SORT_ORDER=queryname METRICS_FILE=Exome_Tumor_mrkdup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT
```
#### Position sort bam file

```bash
cd ~/workspace/data/DNA_alignments/chr6+chr17
java -Xmx12g -jar $PICARD SortSam I=Exome_Norm_namesorted_mrkdup.bam O=Exome_Norm_sorted_mrkdup.bam SO=coordinate
java -Xmx12g -jar $PICARD SortSam I=Exome_Tumor_namesorted_mrkdup.bam O=Exome_Tumor_sorted_mrkdup.bam SO=coordinate
```
#### Create bam index for use with GATK, IGV, etc.

```bash
cd ~/workspace/data/DNA_alignments/chr6+chr17
java -Xmx12g -jar $PICARD BuildBamIndex I=Exome_Norm_sorted_mrkdup.bam
java -Xmx12g -jar $PICARD BuildBamIndex I=Exome_Tumor_sorted_mrkdup.bam
```
#### Calculate BQSR Table

```bash
cd ~/workspace/data/DNA_alignments/chr6+chr17
gatk --java-options '-Xmx24g' BaseRecalibrator -R ~/workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa -I ~/workspace/data/DNA_alignments/chr6+chr17/Exome_Norm_sorted_mrkdup.bam -O ~/workspace/data/DNA_alignments/chr6+chr17/Exome_Norm_sorted_mrkdup_bqsr.table --known-sites ~/workspace/data/raw_data/references/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites ~/workspace/data/raw_data/references/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites ~/workspace/data/raw_data/references/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --preserve-qscores-less-than 6 --disable-bam-index-caching -L chr6 -L chr17
gatk --java-options '-Xmx24g' BaseRecalibrator -R ~/workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa -I ~/workspace/data/DNA_alignments/chr6+chr17/Exome_Tumor_sorted_mrkdup.bam -O ~/workspace/data/DNA_alignments/chr6+chr17/Exome_Tumor_sorted_mrkdup_bqsr.table --known-sites ~/workspace/data/raw_data/references/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites ~/workspace/data/raw_data/references/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites
~/workspace/data/raw_data/references/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --preserve-qscores-less-than 6 --disable-bam-index-caching -L chr6 -L chr17

```

#### Apply BQSR

```bash
cd ~/workspace/data/DNA_alignments/chr6+chr17

gatk --java-options '-Xmx24g' ApplyBQSR -R ~/workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa -I ~/workspace/data/DNA_alignments/chr6+chr17/Exome_Norm_sorted_mrkdup.bam -O ~/workspace/data/DNA_alignments/chr6+chr17/Exome_Norm_sorted_mrkdup_bqsr.bam --bqsr-recal-file ~/workspace/data/DNA_alignments/chr6+chr17/Exome_Norm_sorted_mrkdup_bqsr.table --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30
gatk --java-options '-Xmx24g' ApplyBQSR -R ~/workspace/data/raw_data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa -I ~/workspace/data/DNA_alignments/chr6+chr17/Exome_Tumor_sorted_mrkdup.bam -O ~/workspace/data/DNA_alignments/chr6+chr17/Exome_Tumor_sorted_mrkdup_bqsr.bam --bqsr-recal-file ~/workspace/data/DNA_alignments/chr6+chr17/Exome_Tumor_sorted_mrkdup_bqsr.table --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30
```

#### Clean up un-needed sam/bam files

```bash
cd ~/workspace/data/DNA_alignments/chr6+chr17
mkdir final
mv *_sorted_mrkdup_bqsr.* final/
mv *.txt final/

rm *.sam
rm *.bam
rm *.bai
```
### Assess alignment quality

TO DO: Add Quality Control sections:

See docs here: https://github.com/genome/cancer-genomics-workflow/wiki/Alignment
