# change wd
cd /workspace/align

# run per lane alignments
bwa mem -t 8 -Y -R "@RG\tID:WGS_Tumor_Lane1\tPL:ILLUMINA\tPU:D1VCPACXX.1\tLB:wgs_tumor_lib1\tSM:HCC1395_DNA" -o /workspace/align/WGS_Tumor_Lane1.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane1_R1.fastq.gz /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane1_R2.fastq.gz &
bwa mem -t 8 -Y -R "@RG\tID:WGS_Tumor_Lane2\tPL:ILLUMINA\tPU:D1VCPACXX.2\tLB:wgs_tumor_lib1\tSM:HCC1395_DNA" -o /workspace/align/WGS_Tumor_Lane2.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane2_R1.fastq.gz /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane2_R2.fastq.gz &
bwa mem -t 8 -Y -R "@RG\tID:WGS_Tumor_Lane3\tPL:ILLUMINA\tPU:D1VCPACXX.3\tLB:wgs_tumor_lib2\tSM:HCC1395_DNA" -o /workspace/align/WGS_Tumor_Lane3.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane3_R1.fastq.gz /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane3_R2.fastq.gz &
bwa mem -t 8 -Y -R "@RG\tID:WGS_Tumor_Lane4\tPL:ILLUMINA\tPU:D1VCPACXX.4\tLB:wgs_tumor_lib2\tSM:HCC1395_DNA" -o /workspace/align/WGS_Tumor_Lane4.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane4_R1.fastq.gz /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane4_R2.fastq.gz &
bwa mem -t 8 -Y -R "@RG\tID:WGS_Tumor_Lane5\tPL:ILLUMINA\tPU:D1VCPACXX.5\tLB:wgs_tumor_lib3\tSM:HCC1395_DNA" -o /workspace/align/WGS_Tumor_Lane5.sam /workspace/inputs/references/genome/ref_genome.fa /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane5_R1.fastq.gz /workspace/inputs/data/fastq/WGS_Tumor/WGS_Tumor_Lane5_R2.fastq.gz &
wait

# convert sam to bam
samtools view -@ 8 -h -b -o WGS_Tumor_Lane1.bam WGS_Tumor_Lane1.sam
samtools view -@ 8 -h -b -o WGS_Tumor_Lane2.bam WGS_Tumor_Lane2.sam
samtools view -@ 8 -h -b -o WGS_Tumor_Lane3.bam WGS_Tumor_Lane3.sam
samtools view -@ 8 -h -b -o WGS_Tumor_Lane4.bam WGS_Tumor_Lane4.sam
samtools view -@ 8 -h -b -o WGS_Tumor_Lane5.bam WGS_Tumor_Lane5.sam

# merge bam files
samtools merge -@ 8 WGS_Tumor_merged.bam WGS_Tumor_Lane1.bam WGS_Tumor_Lane2.bam WGS_Tumor_Lane3.bam WGS_Tumor_Lane4.bam WGS_Tumor_Lane5.bam

# name sort bam files
java -Xmx30g -jar $PICARD SortSam I=WGS_Tumor_merged.bam O=WGS_Tumor_merged_namesorted.bam SO=queryname

# mark duplicates
java -Xmx30g -jar $PICARD MarkDuplicates I=WGS_Tumor_merged_namesorted.bam O=WGS_Tumor_merged_namesorted_mrkdup.bam ASSUME_SORT_ORDER=queryname METRICS_FILE=WGS_Tumor_mrkdup_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT

# position sort bam files
java -Xmx30g -jar $PICARD SortSam I=WGS_Tumor_merged_namesorted_mrkdup.bam O=WGS_Tumor_merged_sorted_mrkdup.bam SO=coordinate

# index the bam
java -Xmx30g -jar $PICARD BuildBamIndex I=WGS_Tumor_merged_sorted_mrkdup.bam

# calculate bqsr table
gatk --java-options '-Xmx30g' BaseRecalibrator -R /workspace/inputs/references/genome/ref_genome.fa -I /workspace/align/WGS_Tumor_merged_sorted_mrkdup.bam -O /workspace/align/WGS_Tumor_merged_sorted_mrkdup_bqsr.table --known-sites /workspace/inputs/references/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites /workspace/inputs/references/gatk/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites /workspace/inputs/references/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --preserve-qscores-less-than 6 --disable-bam-index-caching $GATK_REGIONS

# apply bqsr
gatk --java-options '-Xmx30g' ApplyBQSR -R /workspace/inputs/references/genome/ref_genome.fa -I /workspace/align/WGS_Tumor_merged_sorted_mrkdup.bam -O /workspace/align/WGS_Tumor_merged_sorted_mrkdup_bqsr.bam --bqsr-recal-file /workspace/align/WGS_Tumor_merged_sorted_mrkdup_bqsr.table --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30

# create index for final bam
java -Xmx30g -jar $PICARD BuildBamIndex I=WGS_Tumor_merged_sorted_mrkdup_bqsr.bam
